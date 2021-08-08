import bpy
from bpy.types import Panel
from bpy.props import *
import math

default_surface_names = [
    ("bcc", "bcc", "", 1),
    ("schwarzp", "schwarzp", "", 2),
    ("schwarzd", "schwarzd", "", 3),
    ("gyroid", "gyroid", "", 4),
    ("double-p", "double-p", "", 5),
    ("double-d", "double-d", "", 6),
    ("double-gyroid", "double-gyroid", "", 7),
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
    def set_coff(self, context):
        if self.target1 is not None:
            if self.target1.dimensions is not None:
                L = min(self.target1.dimensions)
                N = self.unit_cell
                self.coff = 2*math.pi*N/L
    target1: PointerProperty(name="target MESH", type=bpy.types.Object, update=set_coff)
    target2: PointerProperty(name="target MESH", type=bpy.types.Object)
    is_intersect: BoolProperty(name="Intersection", default=True)
    is_build_inverse: BoolProperty(name="Build inverse", default=False)
    
    grid_size: IntProperty(name="Grid size", default=100, min=10)
    grid_offset: IntProperty(name="Grid offset", default=3, min=0)
    smooth_step: IntProperty(name="smooth step", default=0, min=0)
    k_slice: IntProperty(name="k slice", default=100, min=1)
    k_polygon: IntProperty(name="k polygon", default=4, min=1)
    
    coff: FloatProperty(name="coff", default=math.pi)
    unit_cell: FloatProperty(name="unit_cell", default=1, min=0, update=set_coff)
    isolevel: FloatProperty(name="isolevel", default=0.0)
    qsim_percent: FloatProperty(name="qsim", default=0.0, min=0, max=1.0, subtype="FACTOR")
    minimum_diameter: FloatProperty(name="minimum diameter", default=0.25, min=0, max=1, subtype="FACTOR")
    
    lua_file: StringProperty(name="Lua file", default="", subtype="FILE_PATH")
    surface_name: EnumProperty(name="Surface name", items=default_surface_names)
    direction: EnumProperty(name="slice direction", items=default_direction)

    progress1: IntProperty(name="progress1", default=0, min=0, max=100)
    progress2: IntProperty(name="progress2", default=0, min=0, max=100)
    
    result1: StringProperty(name="result1")
    result2: StringProperty(name="result2")
    
    def get_progress1(self):
        return self.progress1
    
    def get_progress2(self):
        return self.progress2
    
    def null_set(self, value):
        pass
    
    readonly_progress1: IntProperty(min=0, max=100, subtype="PERCENTAGE", get=get_progress1, set=null_set)
    readonly_progress2: IntProperty(min=0, max=100, subtype="PERCENTAGE", get=get_progress2, set=null_set)