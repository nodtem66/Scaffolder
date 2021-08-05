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
    "location" : "View3D > Sidebar > Scaffolder",
    "category" : "Object",
    "tracker_url": "https://github.com/nodtem66/Scaffolder"
}

# Blender imports
import bpy
from bpy.props import PointerProperty

from . utils import is_module_available, ensure_site_packages, ShowMessageBox
from . SCAFFOLDER_OP_generate_mesh import SCAFFOLDER_OP_generate_mesh
from . SCAFFOLDER_OP_slice_test import SCAFFOLDER_OP_slice_test
from . SCAFFOLDER_PT_generate_mesh import SCAFFOLDER_PT_generate_mesh
from . SCAFFOLDER_PT_slice_test import SCAFFOLDER_PT_slice_test
from . SCAFFOLDER_settings import SCAFFOLDER_settings

ensure_site_packages([
    ("PyScaffolder", "PyScaffolder"),
])

try:
    import PyScaffolder
except:
    try:
        from . import PyScaffolder
    except:
        print("PyScaffolder was not found. Please install with pip install PyScaffolder or unpack wheel from PyPI to addon directory")

    
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
    bpy.types.Scene.scaffolder_settings = PointerProperty(type=SCAFFOLDER_settings)

def unregister():
    for c in classes:
        bpy.utils.unregister_class(c)


if __name__ == "__main__":
    register()