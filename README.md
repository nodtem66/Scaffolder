# Scaffolder
![GitHub last commit](https://img.shields.io/github/last-commit/nodtem66/Scaffolder)

![Build Status](https://github.com/nodtem66/Scaffolder/actions/workflows/binary.yml/badge.svg) [![Build Status](https://dev.azure.com/n66/PublicCI/_apis/build/status%2Fnodtem66.Scaffolder%20(Anaconda%20release)?branchName=master)](https://dev.azure.com/n66/PublicCI/_build/latest?definitionId=8&branchName=master)

[![Anaconda](https://anaconda.org/nodtem66/scaffolder/badges/version.svg)](https://anaconda.org/nodtem66/scaffolder) ![Last update](https://anaconda.org/nodtem66/scaffolder/badges/latest_release_date.svg) ![Platform](https://anaconda.org/nodtem66/scaffolder/badges/platforms.svg)

[![Socket Badge](https://socket.dev/api/badge/pypi/package/PyScaffolder/1.5.2.post2?artifact_id=tar-gz)](https://socket.dev/pypi/package/PyScaffolder/overview/1.5.2.post2/tar-gz)

![Scaffolder Logo](https://github.com/nodtem66/Scaffolder/raw/master/docs/images/scaffolder_logo.jpg)


Transform a 3D model from STL/PLY/OFF/OBJ to a porous model with implicit function (e.g., Schwarz P/ Gyroid).

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
- [libigl](https://libigl.github.io/) - The computational geometry library
- [vcglib](https://github.com/cnr-isti-vclab/vcglib) - The mesh utility library
- [sol2](https://github.com/ThePhD/sol2) - Lua script integration
- [tbb](https://github.com/oneapi-src/oneTBB) - Threading library

## How it works
- Read STL file and finding the boundary box
- Generate the grid and calculate the winding number with STL mesh
- Use winding number to determine the condition for [implicit isosurface function](https://wewanttolearn.wordpress.com/2019/02/03/triply-periodic-minimal-surfaces/)
- Generate the isosurface field in the same-size grid
- Perform [Dual marching cube](https://github.com/dominikwodniok/dualmc) to construct the manifold
- Clean up the duplicated vertices or faces, and abandon the group of connected faces having the diameter below the setting
- Export to the target 3D format

## FAQ

### How can I find the dataset from a study of coffient and isolevel?
The raw dataset is available at [Mendeley Data](https://data.mendeley.com/datasets/sbxr7xxvnd/2).
The program that used to generate that data was released at [Github repository](https://github.com/nodtem66/Scaffolder). You can also find the interactive visualization at [Google Colab](https://colab.research.google.com/github/nodtem66/Scaffolder/blob/master/data/data_visualization.ipynb)

### Where is the implicit functions were defined in the C++ sourcecode?
https://github.com/nodtem66/Scaffolder/blob/master/include/implicit_function.h

### Can you suggest alternative softwares like this program?
- Rhino (Grasshopper)
- nTopology
- Hyperganic

## References
- [Minimal surface Blog](https://minimalsurfaces.blog/)
- Dual marching cube - [dualmc](https://github.com/dominikwodniok/dualmc)
- Command line parser - [cxxopts](https://github.com/jarro2783/cxxopts)
- Progress bar - [ProgressBar](https://github.com/prakhar1989/progress-cpp)

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
