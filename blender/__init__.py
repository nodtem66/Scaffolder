# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

bl_info = {
    "name" : "Scaffolder",
    "author" : "Jirawat I.",
    "description" : "Wrapper for Scaffolder program",
    "blender" : (2, 80, 0),
    "version" : (1, 5, 1),
    "location" : "View3D",
    "category" : "Object",
    "tracker_url": "https://github.com/nodtem66/Scaffolder"
}

# Blender imports
import bpy
from bpy.types import Panel
from bpy.props import *
import math
import time

class SCAFFOLDER_OP_generate_mesh(bpy.types.Operator):
    bl_idname = "object.scaffolder_generate_mesh"
    bl_label = "Generate mesh"
    bl_options = {"REGISTER", "UNDO"}
    
    def execute(self, context):
        props = context.scene.scaffolder_settings
        if props:
            self.report(
                {'ERROR'}, 'GS: %d' %
                (props.grid_size)
            )
            def callback(v):
                props.progress1 = v
            callback(50)
        return {'FINISHED'}

class SCAFFOLDER_OP_slice_test(bpy.types.Operator):
    bl_idname = "object.scaffolder_slice_test"
    bl_label = "Slice Test"
    bl_options = {"REGISTER", "UNDO"}
        
    def execute(self, context):
        props = context.scene.scaffolder_settings
        if props:
            self.report(
                {'ERROR'}, 'GS: %d' %
                (props.grid_size)
            )
            props.progress2 = 50
        return {'FINISHED'}

default_surface_names = [
    ("BCC", "bcc", "", 1),
    ("schwarzp", "schwarzp", "", 2),
    ("schwarzd", "schwarzd", "", 3),
    ("gyroid", "gyroid", "", 4),
    ("double-p", "double-p", "", 5),
    ("double-d", "double-d", "", 6),
    ("double-gyroiod", "double-gyroiod", "", 7),
    ("lidinoid", "lidinoid", "", 8),
    ("schoen_iwp", "schoen_iwp", "", 9),
    ("neovius", "neovius", "", 10),
    ("tubular_g_ab", "tubular_g_ab", "", 11),
    ("tubular_g_c", "tubular_g_c", "", 12)
]

default_direction = [
    ("X", "X", "", 0),
    ("Y", "Y", "", 1),
    ("Z", "Z", "", 2),
    ("All", "All", "", 3)
]

class SCAFFOLDER_settings(bpy.types.PropertyGroup):
    target1: PointerProperty(type=bpy.types.Object)
    target2: PointerProperty(type=bpy.types.Object)
    is_build_inverse: BoolProperty(name="Build inverse", default=False)
    
    grid_size: IntProperty(name="Grid size", default=100, min=10)
    grid_offset: IntProperty(name="Grid offset", default=5, min=0)
    smooth_step: IntProperty(name="smooth step", default=5, min=0)
    k_slice: IntProperty(name="k slice", default=100, min=1)
    k_polygon: IntProperty(name="k polygon", default=4, min=1)
    
    coff: FloatProperty(name="coff", default=math.pi)
    isolevel: FloatProperty(name="isolevel", default=0.0)
    qsim_percent: FloatProperty(name="qsim", default=0.0, min=0, max=1.0, subtype="FACTOR")
    minimum_diameter: FloatProperty(name="minimum diameter", default=0.25, min=0, max=1, subtype="FACTOR")
    
    lua_file: StringProperty(name="Lua file", default="", subtype="FILE_PATH")
    surface_name: EnumProperty(name="Surface name", items=default_surface_names)
    direction: EnumProperty(name="slice direction", items=default_direction)

    progress1: IntProperty(name="progress1", default=0, min=0, max=100)
    progress2: IntProperty(name="progress2", default=0, min=0, max=100)
    
    def get_progress1(self):
        return self.progress1
    
    def get_progress2(self):
        return self.progress2
    
    def null_set(self, value):
        pass
    
    readonly_progress1: IntProperty(min=0, max=100, subtype="PERCENTAGE", get=get_progress1, set=null_set)
    readonly_progress2: IntProperty(min=0, max=100, subtype="PERCENTAGE", get=get_progress2, set=null_set)

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
        
        row = layout.row()
        row.label(icon="SPHERE")
        col = row.column()
        col.prop(settings, "surface_name", text="")
        col.prop(settings, "lua_file", text="")
        
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
        if settings.progress1 > 0:
            row.prop(settings, "readonly_progress1", text="Progress", slider=True)
        else:
            row.operator("object.scaffolder_generate_mesh")

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
        
        row = layout.row()
        if settings.progress2 > 0:
            row.prop(settings, "readonly_progress2", text="Progress", slider=True)
        else:
            row.operator("object.scaffolder_slice_test")
     
classes = (
    SCAFFOLDER_settings,
    SCAFFOLDER_OP_generate_mesh,
    SCAFFOLDER_OP_slice_test,
    SCAFFOLDER_PT_generate_mesh,
    SCAFFOLDER_PT_slice_test
)

def register():
    for c in classes:
        bpy.utils.register_class(c)

def unregister():
    for c in classes:
        bpy.utils.unregister_class(c)


if __name__ == "__main__":
    register()
    bpy.types.Scene.scaffolder_settings = PointerProperty(type=SCAFFOLDER_settings)
