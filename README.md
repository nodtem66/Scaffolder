# Scaffolder 
![Build Status](https://github.com/nodtem66/Scaffolder/workflows/Build/badge.svg) [![Build Status](https://dev.azure.com/n66/PublicCI/_apis/build/status/nodtem66.Scaffolder%20(Anaconda)?branchName=master)](https://dev.azure.com/n66/PublicCI/_build/latest?definitionId=6&branchName=master)

[![Anaconda](https://anaconda.org/nodtem66/scaffolder/badges/version.svg)](https://anaconda.org/nodtem66/scaffolder) ![Last update](https://anaconda.org/nodtem66/scaffolder/badges/latest_release_date.svg) ![Platform](https://anaconda.org/nodtem66/scaffolder/badges/platforms.svg) ![Install](https://anaconda.org/nodtem66/scaffolder/badges/installer/conda.svg)

![Scaffolder Logo](https://github.com/nodtem66/Scaffolder/raw/master/docs/images/scaffolder_logo.jpg)

Generate scaffold from STL/PLY/OFF/OBJ file with implicit function (e.g., Schwarz P/ Gyroid).

[Documentation](https://nodtem66.github.io/Scaffolder) 

## Binary installation
* Download from [Github Releases page](https://github.com/nodtem66/Scaffolder/releases) or
* Install from Anaconda
```bash
conda install -c nodtem66 scaffolder
```

## Python supports
```bash
pip install PyScaffolder
```

## Blender addon
- Install the `Scaffolder-blender.zip` downloaded from [Release](https://github.com/nodtem66/Scaffolder/releases/tag/v1.5.1)
- The plugin will appear at `View > Sidebar` or `Press N`

## Screenshots

- **Blender plugin with PyScaffolder**

![Blender plugin](https://github.com/nodtem66/Scaffolder/raw/master/docs/images/blender-plugin.gif)

- **The figure of patterns implemented in this program**

![TPMS Patterns](https://github.com/nodtem66/Scaffolder/raw/master/docs/images/patterns.jpg)


- **The examples of generated porous scaffold**

![Examples porous scaffold](https://github.com/nodtem66/Scaffolder/raw/master/docs/images/examples.jpg)

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
[Computational method and program for generating a porous scaffold based on implicit surfaces](https://doi.org/10.1016/j.cmpb.2021.106088)
```bibtex
@article{IAMSAMANG2021106088,
title = {Computational method and program for generating a porous scaffold based on implicit surfaces},
journal = {Computer Methods and Programs in Biomedicine},
volume = {205},
pages = {106088},
year = {2021},
issn = {0169-2607},
doi = {https://doi.org/10.1016/j.cmpb.2021.106088},
url = {https://www.sciencedirect.com/science/article/pii/S0169260721001632},
author = {Jirawat Iamsamang and Phornphop Naiyanetr},
keywords = {Triply periodic minimal surface (TPMS), Implicit surface, Porous scaffold, Pore size, Porosity}
}
```
