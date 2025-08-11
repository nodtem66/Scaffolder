# Scaffolder
![GitHub Release](https://img.shields.io/github/v/release/nodtem66/scaffolder) ![GitHub last commit](https://img.shields.io/github/last-commit/nodtem66/Scaffolder) 
[![PyPI - Version](https://img.shields.io/pypi/v/PyScaffolder)](https://pypi.org/project/PyScaffolder/) [![Socket Badge](https://socket.dev/api/badge/pypi/package/PyScaffolder/1.5.3?artifact_id=tar-gz)](https://socket.dev/pypi/package/PyScaffolder/overview/1.5.3/tar-gz)

![Scaffolder Logo](https://github.com/nodtem66/Scaffolder/raw/master/docs/images/scaffolder_logo.jpg)


Transform a 3D model from STL/PLY/OFF/OBJ to a porous model with implicit function (e.g., Schwarz P/ Gyroid).

| [Documentation](https://nodtem66.github.io/Scaffolder) | [Command line options](https://github.com/nodtem66/Scaffolder/releases/tag/v1.5.2) | [Python API](https://nodtem66.github.io/Scaffolder/python/) | [Example: TPMS with Python](https://nodtem66.github.io/Scaffolder/tutorial_2/) |
| --- | ---- | ---- | ---- |

## Binary installation
* Download from [Github Releases page](https://github.com/nodtem66/Scaffolder/releases)
* Binary bundle consists of
  * Main program: `Scaffolder`
  * Pore analysis program: `Scaffold.SliceTest`
  * Python library: `PyScaffold*.pyd`

> [!note]
> 1. If you need only python library, you can install it from `pip` (see below)
> 2. From v1.5.3, I decided to end the support for binary executables as no one uses it. Only python library will be released.  

## Python supports
```bash
pip install PyScaffolder
```
>  [!note]
>  v1.5.3 support only Windows and Linux and Python [3.8 - 3.13](https://pypi.org/project/PyScaffolder/#files)

## Binary building
To build the binary executables, make sure you have installed the following softwares:
- Visual Studio >=2022 (Windows)
- GCC and Cmake (linux)
  ```bash
  sudo apt install cmake
  sudo apt install build-essential checkinstall zlib1g-dev libssl-dev -y
  ```
- XCode (MacOS)

### Steps
- Download or clone the source code from github
- Change current directory to the source code folder
- Create `build` folder using CMAKE
```bash
cmake -E make_directory ./build
cd ./build
```
- Configure Cmake variables
```bash
# Linux
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=11
# Window
cmake .. -DCMAKE_CXX_COMPILER=cl -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=11
# MacOS
cmake .. -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=11
```
- Start Building
```bash
cmake --build . --config Release
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
- Read STL file and find the boundary box of STL mesh
- Generate the grid and calculate the winding number with STL mesh
- Use winding number to determine the condition for [implicit isosurface function](https://wewanttolearn.wordpress.com/2019/02/03/triply-periodic-minimal-surfaces/)
- Generate the isosurface field in the grid
- Perform [Dual marching cube](https://github.com/dominikwodniok/dualmc) to construct the manifold
- Clean up the duplicated vertices or faces, and abandon the group of connected faces having the diameter below the setting
- Export to the target 3D format

## FAQ

### How can I find the dataset from the paper?
The raw dataset is available at [Mendeley Data](https://data.mendeley.com/datasets/sbxr7xxvnd/2).
The program that used to generate that data was released at [Github repository](https://github.com/nodtem66/Scaffolder). You can also find the interactive visualization at [Google Colab](https://colab.research.google.com/github/nodtem66/Scaffolder/blob/master/data/data_visualization.ipynb)

### Where are the implicit functions defined in the C++ source code?
https://github.com/nodtem66/Scaffolder/blob/master/include/implicit_function.h

### Can you suggest alternative softwares like this program?
- [Cesogen](https://git.sr.ht/~paulapatience/cesogen) (There you can find a curated list of the other software)
- [TPMS2STEP](https://github.com/Qiang-Zou/TPMS2STEP)
- [LattGen](https://www.sciencedirect.com/science/article/pii/S2665963824000538)
- Rhino (Grasshopper) (Commerical software)
- nTopology (Commerical software)
- Hyperganic (Commerical software)

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
