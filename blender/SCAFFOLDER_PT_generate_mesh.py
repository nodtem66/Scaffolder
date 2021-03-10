import bpy
from bpy.types import Panel
from bpy.props import *

class SCAFFOLDER_PT_generate_mesh(Panel):
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_label = "Generate Mesh"
    bl_category = "Scaffolder"
    
    def draw(self, context):
        
        layout = self.layout
        settings = context.scene.scaffolder_settings
        
        row = layout.row()
        row.label(icon="TRACKER")
        row.prop_search(settings, "target1", context.scene, "objects", text="Target")

        if settings.target1 is not None:
            row = layout.row()
            d = settings.target1.dimensions
            bbox = "[%.3f %.3f %.3f]" % (d.x, d.y, d.z)
            row.label(text="Bounding box:")
            row = layout.row()
            row.label(icon="BLANK1", text=bbox)
        
        row = layout.row()
        row.label(icon="SPHERE")
        col = row.column()
        col.prop(settings, "surface_name", text="")
        col.prop(settings, "lua_file", text="")
        
        row = layout.row()
        row.label(icon="EVENT_N")
        row.prop(settings, "unit_cell", text="#Unit cell (N)")

        row = layout.row()
        row.label(icon="EVENT_W")
        row.prop(settings, "coff", text="Angular freq. (w)")
        
        row = layout.row()
        row.label(icon="EVENT_T")
        row.prop(settings, "isolevel", text="Isolevel (t)")

        row = layout.row()
        row.label(icon="MESH_GRID")
        col = row.column()
        col.prop(settings, "grid_size", text="Grid size")
        col.prop(settings, "grid_offset", text="Grid offset")
        
        row = layout.row()
        row.label(icon="RNDCURVE")
        row.prop(settings, "smooth_step", text="Smooth step")
        
        row = layout.row()
        row.label(icon="FULLSCREEN_EXIT")
        col = row.column()
        col.prop(settings, "qsim_percent", text="Reduce size")
        col.prop(settings, "minimum_diameter", text="Cleaning diameter")
        
        row = layout.row()
        row.label(icon="MOD_MIRROR")
        row.prop(settings, "is_build_inverse", text=" Build Inverse Mesh")

        row = layout.row()
        row.label(icon="SELECT_INTERSECT")
        row.prop(settings, "is_intersect", text=" Intersect with target")
        
        if settings.result1:
            texts = settings.result1.split('\n')
            for t in texts:
                row = layout.row()
                row.label(text=t)
                
        row = layout.row()
        if settings.progress1 > 0 and settings.progress1 < 100:
            row.prop(settings, "readonly_progress1", text="Progress", slider=True)
        else:
            row.operator("object.scaffolder_generate_mesh")