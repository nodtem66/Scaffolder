import numpy as np
import sys
import bpy

def ShowMessageBox(message = "", title = "Message Box", icon = 'INFO'):

    def draw(self, context):
        self.layout.label(text=message)

    bpy.context.window_manager.popup_menu(draw, title = title, icon = icon)

def is_module_available(module_name):
    if sys.version_info < (3, 0):
        # python 2
        import importlib
        torch_loader = importlib.find_loader(module_name)
    elif sys.version_info <= (3, 3):
        # python 3.0 to 3.3
        import pkgutil
        torch_loader = pkgutil.find_loader(module_name)
    elif sys.version_info >= (3, 4):
        # python 3.4 and above
        import importlib
        torch_loader = importlib.util.find_spec(module_name)

    return torch_loader is not None

def read_verts(mesh):
    vn = len(mesh.vertices)
    v = np.zeros((vn*3), dtype=np.float)
    mesh.vertices.foreach_get("co", v)
    return np.reshape(v, (vn, 3))

def read_faces(mesh):
    fn = len(mesh.loop_triangles)
    f = np.zeros((fn*3), dtype=np.int)
    mesh.loop_triangles.foreach_get("vertices", f)
    return np.reshape(f, (fn, 3))
