# Installation

Scaffolder supports the command-line usage by [standalone program](#standalone-program), [python scripting](#python-supports) via [PyPI](https://pypi.org/project/PyScaffolder/), and blender by [python-wrapped plugin](#blender-addon).

## Standalone program

The binary programs for various platforms are automatically built and deposited on the public repository: [Github](#github) and [Anaconda](#anaconda). 

### Github
1. Download `bin64.zip` from [Github Releases page](https://github.com/nodtem66/Scaffolder/releases)
2. Extract `bin64.zip` which contains the following files: 
```
├───macos
│       PyScaffolder.cpython-39-darwin.so
│       Scaffolder
│       Scaffolder.SliceTest
│
├───ubuntu
│       PyScaffolder.cpython-36m-x86_64-linux-gnu.so
│       Scaffolder
│       Scaffolder.SliceTest
│
└───windows
    └───Release
            PyScaffolder.cp37-win_amd64.pyd
            PyScaffolder.exp
            PyScaffolder.lib
            Scaffolder.exe
            Scaffolder.SliceTest.exe
```
3. Open command line at the directory according to your platform, such as `windows/Release` for window.
4. Test the program with `--help`
```bash
Scaffolder --help
Scaffolder.SliceTest --help
``` 
!!! note
    `PyScaffolder.cpython-<python_version>-<platform>` is a compiled python bytecode that can import to a python script as a module (see [Python supports](#python-supports)). 
### Anaconda
1. Install [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Conda can manage the environments and revisions, then:
    * If you want to install `Scaffolder` into current environment, use:
  ```bash
  conda install -c nodtem66 scaffolder
  ```
    * If you want to install `Scaffolder` into [a new environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) due to the conflicted from some dependencies, use:
  ```bash
  conda install -n new_env -c nodtem66 scaffolder
  conda activate new_env
  ```
3. Test the program with `--help`
```bash
Scaffolder --help
Scaffolder.SliceTest --help
```

## Python supports
The current version of `Scaffolder` supports python `2.7`, `3.6`, `3.7`, `3.8`, and `3.9` with only `x86_64` platform on:

* Window 10 (didn't test with the older version)
* Linux
* OSX 10.9

!!! note
    see the list of available version at [PyPI][]

Then install via `pip`
```bash
pip install PyScaffolder
```

The following script is an example usage of `PyScaffolder`.
```python
import PyScaffolder

# Load vertices and faces of a 20mm cube
v = [
	[0.0, 0.0, 0.0],
	[20.0, 0.0, 0.0],
	[0.0, -20.0, 0.0],
	[20.0, -20.0, 0.0],
	[0.0, 0.0, 20.0],
	[20.0, 0.0, 20.0],
	[0.0, -20.0, 20.0],
	[20.0, -20.0, 20.0]
]

f = [
	[1, 2, 0],
	[2, 1, 3],
	[7, 2, 3],
	[2, 7, 6],
	[1, 7, 3],
	[7, 1, 5],
	[7, 4, 6],
	[4, 7, 5],
	[4, 2, 6],
	[2, 4, 0],
	[4, 1, 0],
	[1, 4, 5]
]
print("Vertices: ", len(v))
print("Faces: ", len(f))

# Test Slice test with all direction
a = PyScaffolder.slice_test(v, f, direction=3)
print(len(a.minFeret))

# Set default parameter
params = PyScaffolder.Parameter()
params.coff = 1.0

# Test the mesh generation from v,f
a = PyScaffolder.[generate_mesh]()(v, f, params)
print(a.porosity)
print(a.surface_area_ratio)
print(len(a.v), len(a.f))
```

See [Python API](python.md) for details.

## Blender addon
1. Download `Scaffolder-blender.zip` from [Github](https://github.com/nodtem66/Scaffolder/releases/tag/v1.5.1), which contains following files (see also [/blender](https://github.com/nodtem66/Scaffolder/tree/master/blender)):
```
└───Scaffolder
    SCAFFOLDER_OP_generate_mesh.py
    SCAFFOLDER_OP_slice_test.py
    SCAFFOLDER_PT_generate_mesh.py
    SCAFFOLDER_PT_slice_test.py
    SCAFFOLDER_settings.py
    utils.py
    __init__.py
```
2. Open Blender > Edit > Preferences > Add-ons > Install, select `Scaffolder-blender.zip`, and check the box to enable it.
3. Press <kbd class="key-n">N</kbd> or View > Sidebar. Then Scaffolder will be found in one tab of a sidebar.

See [Blender plugin](blender.md) for video and more information.

<!--
!!! note "Find blender python version and scripting directory"
    1. Open Blender > Scripting tab
    2. In console window will show version of python: `PYTHON INTERACTIVE VERSION 3.x.x`.
    3. Then enter the following code to find scripting directory:
    ```python
    bpy.utils.user_resource('SCRIPTS', "addons")
    ```
-->
[PyPI]: https://pypi.org/project/PyScaffolder/#files