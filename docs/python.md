# Python API

`PyScaffolder` is a wrapper for `Scaffolder`, which is written in C++. 
It used [PyBind11](https://github.com/pybind/pybind11) to interface the C++ core function.
The main core functions are `generate_scaffold` and `slice_test` which resemble the standalone program: `Scaffolder` and `Scaffolder.Slice`.

## generate_scaffold
```python
def generate_scaffold(vertices, faces, params=Parameter(), callback=None)
```
:   Generate 3D mesh from input vertices, faces, and parameters
    ### Parameters
    * `vertices`: 2D array or numpy.array representing vertices in list of [x,y,z]   
    * `faces`: 2D array or numpy.array representing faces in list of [v1, v2, v3, v4]
    * `params`: *Optional* [Parameter](#parameter) see the detail in [Command line](cmd.md)
    * `callback`: *Optional* function with one `integer` parameter indicating the progression
    ### Return
    [MeshInfo](#meshinfo)
    ### Example
    ```python
    from PyScaffolder import generate_scaffold, Parameter

    # vertices and faces of a 20mm cube
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

    # the function required two parameters
    mesh_info = generate_scaffold(v, f)

    # use custom parameters
    parameter = Parameter()
    parameter.coff = 1.0
    
    mesh_info = generate_scaffold(v, f, parameter)

    # use callback to show % progression
    def progression(n):
      print(n)

    mesh_info = generate_scaffold(v, f, callback=progression)
    ```

## slice_test
```python
def slice_test(vertices, faces, k_slice=100, k_polygon=4, direction=0, callback=None)
```
:   Slice input mesh into pore sizes
    ### Parameters
    * `vertices`: 2D array or numpy.array representing vertices in list of [x,y,z]   
    * `faces`: 2D array or numpy.array representing faces in list of [v1, v2, v3, v4]
    * `k_slice`: *Optional* integer $k_{slice}$ the number of sliced layers
    * `k_polygon`: *Optional* integer $k_{polygon}$ the number of nearest-neighbor polygon to measure Feret diameter
    * `direction` *Optional* integer representing slicing axis (0=X, 1=Y, 2=Z, 3=All)
    * `callback`: *Optional* function with one `integer` parameter indicating the progression
    ###Return
    [PoreSize](#poresize)

    ###Example
    ```python
    import numpy as np
    import trimesh
    import PyScaffolder
    
    # load input.stl into mesh
    mesh = trimesh.load('input.stl')

    # slice all direction with default parameter
    pore_size = PyScaffolder.slice_test(mesh.vertices, mesh.faces, direction=3)
    print(len(pore_size.minFeret), len(pore_size.maxFeret))
    ```

## marching_cubes
```python
def marching_cubes(f, grid_size=(100,100,100), v_min=(0,0,0), delta=0.01, clean=False, callback=None)
```
:   Slice input mesh into pore sizes
    ### Parameters
    * `f`: numpy.array representing a discrete isosurface where F(x,y,z) = 0 is boundary   
    * `grid_size`: *Optional* array, list, tuple, or int representing the number of voxels 
    * `delta`: *Optional* array, list, tuple or double representing the dimension of a voxel
    * `v_min`: *Optional* array, list, tuple the coordinate of the corner of grid
    * `clean` *Optional* boolean that enable mesh cleaning after marching cubes
    * `callback`: *Optional* function with one `integer` parameter indicating the progression
    ###Return
    `#!python (v, f)` representing vertices and faces in `np.array`

    ###Example
    ```python
    import PyScaffolder
    Fxyz = [
        1,1,1,1,
        1,-1,-1,1,
        1,-1,-1,1,
        1,1,1,1
    ]*4
    (v, f) = PyScaffolder.marching_cubes(
        Fxyz,
        grid_size=4,
        delta=0.25,
        v_min=(-.5, -.5, -.5),
        callback = lambda x: print(x)
    )
    ```

## Parameter
```python
class Parameter:
    is_build_inverse = False
    is_intersect = True
    verbose = False
    shell = 0
    grid_offset = 5
    smooth_step = 5
    k_slice = 100
    k_polygon = 4
    fix_self_intersect = 0
    grid_size = 100
    isolevel = 0.0
    qsim_percent = 0
    coff = 3.141592653589793238462643383279502884
    minimum_diameter = 0.25
    surface_name = "bcc"
```
:   Parameters implement from standalone program (see [Command line](cmd.md))
    ### Example
    ```python
    from PyScaffolder import Parameter
    params = Parameter()
    params.coff = 1.0
    ```

## PoreSize
```python
class PoreSize:
    minFeret = np.array([])
	maxFeret = np.array([])
```
:   Collection of Minimum and Maximum Feret diameters, normally returned from [slice_test](#slice_test)

## MeshInfo
```python
class MeshInfo:
    v = np.array([[]])
	f = np.array([[]])
	porosity = 0.0
	surface_area = 0.0
	surface_area_ratio = 0.0
```
:   Collection of triangular mesh resulted from [generate_scaffold](#generate_scaffold)