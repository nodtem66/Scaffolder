import bpy
import traceback
from . utils import read_verts, read_faces
import numpy as np

class SCAFFOLDER_OP_slice_test(bpy.types.Operator):
    bl_idname = "object.scaffolder_slice_test"
    bl_label = "Slice Test"
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        props = context.scene.scaffolder_settings
        if props.target2 is not None:
            target_type = getattr(props.target2, 'type', '')
            if target_type == 'MESH':
                return True
        return False
    
    def execute(self, context):
        props = context.scene.scaffolder_settings
        if props:
            self.report({'INFO'}, 'Slice %s' % (props.direction))
            def callback(v):
                props.progress2 = v
            try:
                import PyScaffolder
            except:
                from . import PyScaffolder  
            try:
                mesh = props.target2.to_mesh()
                mesh.calc_loop_triangles()
                vn = len(mesh.vertices)
                v = read_verts(mesh)
                f = read_faces(mesh)
                direction = 3 if props.direction == 'All' else ord(props.direction.upper()) - ord('X')
                poresize = PyScaffolder.slice_test(v, f, props.k_slice, props.k_polygon, direction, callback)
                if len(poresize.minFeret) > 0 and len(poresize.maxFeret) > 0:
                    template_str = "Min diameter:\n   " + \
                        "Mean %.3f      Median %.3f\nMax diameter:\n   " + \
                        "Mean %.3f      Median %.3f"
                    props.result2 = template_str % (
                        np.mean(poresize.minFeret),
                        np.percentile(poresize.minFeret, 0.5),
                        np.mean(poresize.maxFeret),
                        np.percentile(poresize.maxFeret, 0.5),
                        )
                else:
                    props.result2 = "No pores"
            except:
                self.report({'ERROR'}, traceback.format_exc())
        return {'FINISHED'}