import traceback
from queue import Empty, Queue
from threading import Thread

import bpy
import numpy as np
import PyScaffolder


def read_verts(mesh):
    vn = len(mesh.vertices)
    v = np.zeros((vn * 3), dtype=np.float64)
    mesh.vertices.foreach_get("co", v)
    return np.reshape(v, (vn, 3))


def read_faces(mesh):
    fn = len(mesh.loop_triangles)
    f = np.zeros((fn * 3), dtype=np.int64)
    mesh.loop_triangles.foreach_get("vertices", f)
    return np.reshape(f, (fn, 3))


def generate_mesh_worker(result_queue, vertices, faces, parameter_values):
    try:
        params = PyScaffolder.Parameter()
        for name, value in parameter_values.items():
            setattr(params, name, value)

        result = PyScaffolder.generate_scaffold(vertices, faces, params)
        result_queue.put((
            "ok",
            result.v.tolist(),
            result.f.tolist(),
            result.porosity,
            result.surface_area,
            result.surface_area_ratio,
        ))
    except Exception:
        result_queue.put(("error", traceback.format_exc()))


class PYSCAFFOLDER_OP_generate_mesh(bpy.types.Operator):
    bl_idname = "object.pyscaffolder_generate_mesh"
    bl_label = "Generate mesh"
    bl_description = "Create a mesh from the target"
    bl_options = {"REGISTER", "UNDO"}

    _timer = None
    _max_steps = 100
    _step = 0
    _thread = None
    _result_queue = None
    _params = None

    @classmethod
    def poll(cls, context):
        props = context.scene.scaffolder_settings
        if props.target1 is not None:
            target_type = getattr(props.target1, "type", "")
            if target_type == "MESH":
                return True
        return False

    def _clean(self, context):
        self._step = 0
        if self._timer is not None:
            context.window_manager.event_timer_remove(self._timer)
            self._timer = None
        context.window_manager.progress_end()
        context.window.cursor_modal_restore()
        self._result_queue = None

    def cancel(self, context):
        self._clean(context)

    def execute(self, context):
        props = context.scene.scaffolder_settings
        if not props:
            return {"CANCELLED"}

        params = PyScaffolder.Parameter()
        params.coff = props.coff
        params.is_build_inverse = props.is_build_inverse
        params.grid_offset = props.grid_offset
        params.smooth_step = props.smooth_step
        params.k_slice = props.k_slice
        params.k_polygon = props.k_polygon
        params.grid_size = props.grid_size
        params.isolevel = props.isolevel
        params.qsim_percent = props.qsim_percent
        params.minimum_diameter = props.minimum_diameter
        params.surface_name = props.lua_file or props.surface_name
        params.is_intersect = props.is_intersect
        self._params = params

        self.report({"INFO"}, f"Generate {params.surface_name}")

        try:
            source_mesh = props.target1.to_mesh()
            try:
                source_mesh.calc_loop_triangles()
                vertices = read_verts(source_mesh)
                faces = read_faces(source_mesh)
            finally:
                props.target1.to_mesh_clear()

            parameter_values = {
                "coff": params.coff,
                "is_build_inverse": params.is_build_inverse,
                "grid_offset": params.grid_offset,
                "smooth_step": params.smooth_step,
                "k_slice": params.k_slice,
                "k_polygon": params.k_polygon,
                "grid_size": params.grid_size,
                "isolevel": params.isolevel,
                "qsim_percent": params.qsim_percent,
                "minimum_diameter": params.minimum_diameter,
                "surface_name": params.surface_name,
                "is_intersect": params.is_intersect,
            }
            self._result_queue = Queue()
            self._thread = Thread(
                target=generate_mesh_worker,
                args=(self._result_queue, vertices, faces, parameter_values),
                daemon=True,
            )
            self._thread.start()
        except Exception:
            self.report({"ERROR"}, traceback.format_exc())
            return {"CANCELLED"}

        return {"FINISHED"}

    def invoke(self, context, event):
        self._step = 0
        context.window_manager.progress_begin(0, self._max_steps)
        context.window.cursor_modal_set("WAIT")
        result = self.execute(context)
        if result != {"FINISHED"}:
            context.window_manager.progress_end()
            context.window.cursor_modal_restore()
            return result

        context.window_manager.modal_handler_add(self)
        self._timer = context.window_manager.event_timer_add(
            0.05, window=context.window
        )
        return {"RUNNING_MODAL"}

    def modal(self, context, event):
        if event.type in {"RIGHTMOUSE", "ESC"}:
            self.cancel(context)
            return {"CANCELLED"}

        if event.type == "TIMER":
            self._step = (self._step + 1) % (self._max_steps + 1)
            context.window_manager.progress_update(self._step)

            try:
                result = self._result_queue.get_nowait()
            except Empty:
                return {"RUNNING_MODAL"}

            if result[0] == "error":
                self.report({"ERROR"}, result[1])
                self._clean(context)
                return {"CANCELLED"}

            _, vertices, faces, porosity, surface_area, surface_area_ratio = result
            props = context.scene.scaffolder_settings
            props.result1 = (
                f"Porosity: {porosity:.2f}\n"
                f"Surface Area: {surface_area:.2f}\n"
                f"Surface Area Ratio: {surface_area_ratio:.2f}"
            )

            object_name = (
                f"{self._params.surface_name}_{self._params.coff:.3f}_"
                f"{self._params.isolevel}_{self._params.grid_size}"
            )
            mesh_name = f"mesh_{object_name}"
            new_mesh = bpy.data.meshes.new(mesh_name)
            new_mesh.from_pydata(vertices, [], faces)
            new_mesh.update()
            new_object = bpy.data.objects.new(object_name, new_mesh)

            new_collection = bpy.data.collections.get("Scaffolder")
            if new_collection is None:
                new_collection = bpy.data.collections.new("Scaffolder")
                context.scene.collection.children.link(new_collection)
            new_collection.objects.link(new_object)

            self._clean(context)
            return {"FINISHED"}

        return {"RUNNING_MODAL"}
