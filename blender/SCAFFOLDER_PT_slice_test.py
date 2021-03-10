import bpy
from bpy.types import Panel
from bpy.props import *

class SCAFFOLDER_PT_slice_test(Panel):
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_label = "Slice Test"
    bl_category = "Scaffolder"
    
    def draw(self, context):
        
        layout = self.layout
        settings = context.scene.scaffolder_settings
        
        row = layout.row()
        row.label(icon="TRACKER")
        row.prop_search(settings, "target2", context.scene, "objects", text="Target")
        
        row = layout.row()
        row.label(icon="MOD_MULTIRES")
        col = row.column()
        col.prop(settings, "k_slice", text="K slice")
        
        row = layout.row()
        row.label(icon="STICKY_UVS_DISABLE")
        row.prop(settings, "k_polygon", text="K polygon")
        
        row = layout.row()
        row.label(icon="OBJECT_ORIGIN")
        col = row.column()
        col.prop(settings, "direction", text="")
        
        if settings.result2:
            texts = settings.result2.split('\n')
            for t in texts:
                row = layout.row()
                row.label(text=t)
                    
        row = layout.row()
        if settings.progress2 > 0 and settings.progress2 < 100:
            row.prop(settings, "readonly_progress2", text="Progress", slider=True)
        else:
            row.operator("object.scaffolder_slice_test")