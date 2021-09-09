# TPMS with Python through SDF and PyScaffolder

Author: [Jirawat Iamsamang](https://github.com/nodtem66)  
Colab: <https://github.com/nodtem66/Scaffolder/blob/master/docs/jupyter/test_1.ipynb>

## Abstract

SDF provides a class for discretizing and visualizing any implicit surfaces. The basic topologies (e.g. sphere, box) are already defined.

This notebook shows how to utilize this library to generate gyroid surface.



## Installation

* Currently, SDF is not in PyPI. So the [github of SDF](https://github.com/fogleman/sdf) needs to clone into local computer. See [Installation](https://github.com/fogleman/sdf#installation)

* By contrast, PyScaffolder was installed by `pip install PyScaffolder`



## Gyroid

The gyroid function is defined as shown in a following cell.  

The wrapper `@sdf3` will provide gyroid function with 3D points (`p`)).

Then these (x,y, z) points will multiply by `w` and calculate the iso-level of gyroid by vectorized numpy function. 


```python
%load_ext autoreload

%autoreload 2



import numpy as np

import PyScaffolder

from sdf import *



@sdf3

def gyroid(w = 3.14159, t=0):

    def f(p):

        q = w*p

        x, y, z = (q[:, i] for i in range(3))

        return np.cos(x)*np.sin(y) + np.cos(y)*np.sin(z) + np.cos(z)*np.sin(x) - t

    return f
```
    
## Generate with SKimage

SDF used `marching_cubes` from `skimage.measure` with a `ThreadPool`, so it's super fast to construct the 3D mesh.

Let's create a constructing function that intersect a gyroid and a unit box.


```python
f = box(1) & gyroid(w=12)
```


```python
# Generate with skimage.measure.marching_cubes

points = f.generate(step=0.01, verbose=True)

write_binary_stl('out_1.stl', points)
```

    min -0.565721, -0.565721, -0.565721
    max 0.565722, 0.565722, 0.565722
    step 0.01, 0.01, 0.01
    1601613 samples in 64 batches with 16 workers
    
    7 skipped, 0 empty, 57 nonempty
    233958 triangles in 0.659682 seconds
    

## Generate with PyScaffolder



However, this method occasionally results in incomplete mesh.  

Then let's try `Pyscaffolder.marching_cubes` which implements `dual marching cubes` from [@dominikwodniok/dualmc](https://github.com/dominikwodniok/dualmc).


```python
# Generate with PyScaffolder.marching_cubes

def marching_cubes(f, step=0.01, bounds=None, verbose=True, clean=True):

    from sdf.mesh import _estimate_bounds, _cartesian_product

    import time



    if not bounds:

        bounds = _estimate_bounds(f)

    (x0, y0, z0), (x1, y1, z1) = bounds

    

    try:

        dx, dy, dz = step

    except TypeError:

        dx = dy = dz = step

    

    if verbose:

        print('min %g, %g, %g' % (x0, y0, z0))

        print('max %g, %g, %g' % (x1, y1, z1))

        print('step %g, %g, %g' % (dx, dy, dz))



    X = np.arange(x0, x1, dx)

    Y = np.arange(y0, y1, dy)

    Z = np.arange(z0, z1, dz)



    P = _cartesian_product(X, Y, Z)

    try:

        # Since the PyScaffolder marching_cubes aceept FREP: F(x,y,z) > 0

        # Then the negative of implicit function is used

        Fxyz = (-f(P))

        # Reshape to Fortran array (column-based) due to implementation of dualmc starting from z axis to x

        Fxyz = Fxyz.reshape((len(X), len(Y), len(Z))).reshape(-1, order='F')

        start = time.time()

        (v, f) = PyScaffolder.marching_cubes(Fxyz, grid_size=[len(X), len(Y), len(Z)], v_min=bounds[0], delta=step, clean=clean)

        if verbose:

            seconds = time.time() - start

            print('\n%d triangles in %g seconds' % (len(points) // 3, seconds))

        # merge vertices and faces into points

        return v[f].reshape((-1, 3))

    except Exception as e:

        print(e)

        return np.array([])



points = marching_cubes(f, step=0.01, verbose=True, clean=True)

write_binary_stl('out_2.stl', points)
```
