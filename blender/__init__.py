_needs_reload = "bpy" in locals()

import bpy
from bpy.props import PointerProperty

from . import (
    PYSCAFFOLDER_OP_generate_mesh,
    PYSCAFFOLDER_PT_generate_mesh,
    PYSCAFFOLDER_settings,
)

if _needs_reload:
    import importlib

    PYSCAFFOLDER_OP_generate_mesh = importlib.reload(PYSCAFFOLDER_OP_generate_mesh)


_classes = (
    PYSCAFFOLDER_settings.PYSCAFFOLDER_settings,
    PYSCAFFOLDER_OP_generate_mesh.PYSCAFFOLDER_OP_generate_mesh,
    PYSCAFFOLDER_PT_generate_mesh.PYSCAFFOLDER_PT_generate_mesh,
)


def register():
    for c in _classes:
        bpy.utils.register_class(c)
    bpy.types.Scene.scaffolder_settings = PointerProperty(
        type=PYSCAFFOLDER_settings.PYSCAFFOLDER_settings
    )


def unregister():
    for c in _classes:
        bpy.utils.unregister_class(c)
    del bpy.types.Scene.scaffolder_settings


if __name__ == "__main__":
    register()
