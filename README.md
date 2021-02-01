# Scaffolder 
![Build Status](https://github.com/nodtem66/Scaffolder/workflows/Build/badge.svg) [![Build Status](https://dev.azure.com/n66/Public%20CI/_apis/build/status/nodtem66.Scaffolder?branchName=master)](https://dev.azure.com/n66/Public%20CI/_build/latest?definitionId=1&branchName=master)

![Scaffolder Logo](https://github.com/nodtem66/Scaffolder/raw/master/images/scaffolder_logo.jpg)

Generate scaffold from STL/PLY/OFF/OBJ file with implicit function (e.g., Schwarz P/ Gyroid).

```
Scaffolder - generate 3D scaffold from STL file
Usage:
  Scaffolder [OPTION...] INPUT OUTPUT PARAMETERS

  -h, --help                    Print help
  -i, --input INPUT             Input file (STL)
  -o, --output OUTPUT           Output filename with extension
                                stl,ply,obj,off [default: out]
      --params PARAMETERS       Combined parameters list:
                                surface[,coff,isolevel,grid_size,k_slice,k_polygon]
  -q, --quiet                   Disable verbose output [default: false]
  -c, --coff DOUBLE             Angular frequency (pore size adjustment)
                                default:PI
  -t, --isolevel DOUBLE         isolevel (porosity adjustment) [default: 0]
  -n, --surface NAME            implicit surface: rectlinear, schwarzp,
                                schwarzd, gyroid, double-p, double-d,
                                double-gyroiod, lidinoid, schoen_iwp, neovius, bcc,
                                tubular_g_ab, tubular_g_c [default: schwarzp]
  -g, --grid_size INT (0..60000)
                                Grid size [default: 100]
  -s, --shell INT (0..60000)    Outer thickness (layers) [default:0]
      --grid_offset INT (0..60000)
                                [default:3]
  -m, --microstructure          Analysis microstructure with Slice contour
                                technique ( [default: false]
      --export_microstructure   Analysis microstructure and export the 2D
                                contours (for debugging) [default: false]
      --k_slice INT (0..60000)  K_slice: the number of slicing layers in each
                                direction (used in microstructure analysis)
                                [default: 100]
      --k_polygon INT (>0)      K_polygon: the number of closest outer
                                contour (used in microstructure analysis) [default: 4]
  -z, --size_optimize DOUBLE (0..1)
                                Experimental Quadric simplification [default: 0]
      --smooth_step INT (0..60000)
                                Smooth with laplacian (default: 5)
      --dirty                   Disable autoclean [default false]
      --minimum_diameter DOUBLE (0..1)
                                used for removing small orphaned (between 0-1) [default: 0.25]
      --format FORMAT (default, csv)
                                Format of logging output [default: default]
      --output_inverse          additional output inverse scaffold [default: false]
      --fix_self_intersect      Experimental fix self-intersect faces
                                [default: false]

Example:

  Scaffolder input.stl output.stl bcc,3.14159,0,100
    Generated BCC scaffold with w=3.14159 (PI), t=0, and grid size=100

  Scaffolder input.stl output.stl custom.lua,3.14159,0,100,100,4 -m
    Generated and evaluated scaffold with custom.lua, w=3.14159 (PI), t=0,
    grid size=100, k_slice=100, k_polygon=4

  Scaffolder input.stl output.stl bcc,3.14159,0 -m -q --format csv
    Generated and evaluated BCC scaffold (w=3.14159, t=0) and reported in CSV


Lua file:

  Define the "surface" function that return the implicit function
  -----------------------------------------------------------------
  function surface (x, y, z)
    return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x) - t
  end
  -----------------------------------------------------------------

Special symbols can be used in lua file:

  w: value from -c, or --coff
  t: value from -t, or --thickness
  winding(x,y,z): function returning the winding number of position x,y,z
  signed_distance(x,y,z): function returning signed distance of position x,y,z
  and all functions from math module
```

## Screenshots

- **The figure of patterns implemented in this program**

![TPMS Patterns](https://github.com/nodtem66/Scaffolder/raw/master/images/patterns.jpg)

- **The examples of generated porous scaffold**

![Examples porous scaffold](https://github.com/nodtem66/Scaffolder/raw/master/images/examples.jpg)

## Dependencies
- [libigl](https://libigl.github.io/)
- [vcglib](https://github.com/cnr-isti-vclab/vcglib)
- [sol2](https://github.com/ThePhD/sol2)
- [tbb](https://github.com/oneapi-src/oneTBB)

## How it works
- Read STL file and finding the boundary box
- Generate the grid and calculate the winding number with STL mesh
- Use winding number to determine the condition for [implicit isosurface function](https://wewanttolearn.wordpress.com/2019/02/03/triply-periodic-minimal-surfaces/)
- Generate the isosurface field in the same-size grid
- Perform [Dual marching cube](https://github.com/dominikwodniok/dualmc) to construct the manifold
- Clean up the duplicated vertices or faces, and abandon the group of connected faces having the diameter below the setting
- Export to the target 3D format

## Coff and Thickness study
[Angular frequency and iso-level](https://colab.research.google.com/github/nodtem66/Scaffolder/blob/master/data/data_visualization.ipynb)

## References
- [libigl](https://github.com/libigl/libigl)
- [vcglib](https://github.com/cnr-isti-vclab/vcglib)
- [sol2](https://github.com/ThePhD/sol2)
- [tbb](https://github.com/oneapi-src/oneTBB) 
- [dualmc](https://github.com/dominikwodniok/dualmc)
- [cxxopts](https://github.com/jarro2783/cxxopts)
- [ProgressBar](https://github.com/prakhar1989/progress-cpp)
- [Minimal surface Blog](https://minimalsurfaces.blog/)

## Citation
[![DOI](https://zenodo.org/badge/227950058.svg)](https://zenodo.org/badge/latestdoi/227950058)

