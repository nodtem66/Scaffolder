import bpy
import traceback
from . utils import read_verts, read_faces

class SCAFFOLDER_OP_generate_mesh(bpy.types.Operator):
    bl_idname = "object.scaffolder_generate_mesh"
    bl_label = "Generate mesh"
    bl_options = {"REGISTER", "UNDO"}
    
    _timer = None
    th = None
    prog = 0
    stop_early = False
    mesh = None
    
    @classmethod
    def poll(cls, context):
        props = context.scene.scaffolder_settings
        if props.target1 is not None:
            target_type = getattr(props.target1, 'type', '')
            if target_type == 'MESH':
                return True
        return False
    
    def modal(self, context, event):
        if event.type in {'RIGHTMOUSE', 'ESC'}:
            self.cancel(context)

            self.stop_early = True
            self.th.join()
            print('DONE EARLY')

            return {'CANCELLED'}

        if event.type == 'TIMER':
            #context.scene.ProgressWidget_progress = self.prog

            if not self.th.isAlive():
                self.th.join()
                print('DONE')
                return {'FINISHED'}

        return {'PASS_THROUGH'}
    
    def cancel(self, context):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
    
    def execute(self, context):
        import threading
        try:
            import PyScaffolder
        except:
            from . import PyScaffolder    
        props = context.scene.scaffolder_settings
        if props:
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
            params.surface_name = props.lua_file if props.lua_file else props.surface_name
            params.is_intersect = props.is_intersect
            
            self.report({'INFO'}, 'Generate %s' % (params.surface_name))
            def callback(v):
                props.progress1 = v
            try:
                mesh = props.target1.to_mesh()
                mesh.calc_loop_triangles()
                vn = len(mesh.vertices)
                v = read_verts(mesh)
                f = read_faces(mesh)
                
                def gen(self):
                    if ('__version__' in PyScaffolder.__dict__):
                        mesh = PyScaffolder.generate_scaffold(v, f, params, callback)
                    else:
                        mesh = PyScaffolder.generate_mesh(v, f, params, callback)
                    template_str = "Porosity: %.3f\nSurface Area: %.3f\nSurface Area Ratio: %.3f"
                    props.result1 = template_str % (mesh.porosity, mesh.surface_area, mesh.surface_area_ratio)
                    
                    mesh_name = 'mesh_%s_%.3f_%d_%d' % (params.surface_name, params.coff, params.isolevel, params.grid_size)
                    object_name = '%s_%.3f_%d_%d' % (params.surface_name, params.coff, params.isolevel, params.grid_size)
                    new_mesh = bpy.data.meshes.new(mesh_name)
                    new_mesh.from_pydata(mesh.v.tolist(), [], mesh.f.tolist())
                    new_mesh.update()
                    new_object = bpy.data.objects.new(object_name, new_mesh)
                    found = False
                    for c in bpy.data.collections:
                        if c.name == 'Scaffolder':
                            found = True
                            new_collection = c
                    if not found:
                        new_collection = bpy.data.collections.new('Scaffolder')
                        bpy.context.scene.collection.children.link(new_collection)
                    new_collection.objects.link(new_object)
                self.th = threading.Thread(target=gen, args=(self,))
                
                callback(1)
                self.th.start()
                
                wm = context.window_manager
                self._timer = wm.event_timer_add(0.3, context.window)
                wm.modal_handler_add(self)
                
            except:
                self.report({'ERROR'}, traceback.format_exc())
        
        return {'RUNNING_MODAL'}