# Scaffolder 
[![Build Status](https://travis-ci.org/nodtem66/Scaffolder.svg?branch=master)](https://travis-ci.org/nodtem66/Scaffolder) [![Build Status](https://dev.azure.com/n66/Public%20CI/_apis/build/status/nodtem66.Scaffolder?branchName=master)](https://dev.azure.com/n66/Public%20CI/_build/latest?definitionId=1&branchName=master)

Generate scaffold from STL file with implicit function (Schwarz P/ Gyroid).

```
Scaffolder - generate 3D scaffold from STL file
Usage:
  Scaffolder [OPTION...] [option args]

  -h, --help                    Print help
  -q, --quiet                   Disable verbose output [default: false]
  -m, --microstructure          Analysis microstructure ( [default: false]
      --m1                      Export and analysis microstructure 1 (Image
                                processing technique) [default: false]
      --m2                      Export and analysis microstructure 2 (Slice
                                coutour technique) [default: false]
  -f, --format arg              Output format (OFF,PLY,STL,OBJ) [default:
                                ply]
  -i, --input FILE              Input file (STL)
  -o, --output FILENAME         Output filename without extension [default:
                                out]
  -c, --coff DOUBLE             default:4*PI
  -s, --shell INT               [default:0]
  -n, --surface NAME            rectlinear, schwarzp, schwarzd, gyroid,
                                double-p, double-d, double-gyroiod, lidinoid,
                                schoen_iwp, neovius, bcc, tubular_g_ab, tubular_g_c
                                [default: schwarzp]
  -t, --thickness DOUBLE        Thickness [default: 0]
  -g, --grid_size INT           Grid size [default: 100]
      --grid_offset INT         [default:3]
      --smooth_step INT         Smooth with laplacian (default: 5)
      --method 0,1              Method of microstructure analysis: 0 (Image
                                processing technique) or 1 (Slice contour
                                technique) [default: 0]
      --output_inverse          additional output inverse scaffold
      --inverse                 Enable build inverse 3D scaffold (for pore
                                connectivity analysis)
      --dirty                   Disable autoclean
      --fix_self_intersect	    Enable experimental fixing for self-intersect problems (default: false)
      --qsim DOUBLE (0-1)       % Quadric simplification for reducing the resolution of 3D mesh (default: false)
      --minimum_diameter DOUBLE
                                used for removing small orphaned (between
                                0-1) [default: 0.25]
```

## Screenshots

** The figure of patterns implemented in this program

![TPMS Patterns](https://github.com/nodtem66/Scaffolder/raw/master/images/patterns.jpg)

** The examples of generated porous scaffold

![Examples porous scaffold](https://github.com/nodtem66/Scaffolder/raw/master/images/examples.jpg)

## Dependencies
- [libigl](https://libigl.github.io/)
- [vcglib](https://github.com/cnr-isti-vclab/vcglib)
- [diplib](https://github.com/DIPlib/diplib)

## How it works
- Read STL file and finding the boundary box
- Generate the grid and calculate the winding number with STL mesh
- Use winding number to determine the condition for [implicit isosurface function](https://wewanttolearn.wordpress.com/2019/02/03/triply-periodic-minimal-surfaces/)
- Generate the isosurface field in the same-size grid
- Perform [Dual marching cube](https://github.com/dominikwodniok/dualmc) to construct the manifold
- Clean up the duplicated vertices or faces, and abandon the group of connected faces having the diameter below the setting
- Export to the target 3D format

## coff and thickness study
[Angular frequency and iso-level](https://colab.research.google.com/github/nodtem66/Scaffold/data/data_visualization.ipynb)

## Reference 
- [dualmc](https://github.com/dominikwodniok/dualmc)
- [cxxopts](https://github.com/jarro2783/cxxopts)
- [ProgressBar](https://github.com/prakhar1989/progress-cpp)
- [Minimal surface Blog](https://minimalsurfaces.blog/)

## Citation
[![DOI](https://zenodo.org/badge/227950058.svg)](https://zenodo.org/badge/latestdoi/227950058)

