# Scaffolder
Generate scaffold from STL file with implicit function (Schwarz P/ Gyroid).
```
Usage: Scaffolder_2 [options]
Options:
    -q, --quiet            Disable verbose output
    -f, --format           Output format (OFF,PLY,STL,OBJ) [default: ply]
    -i, --input            Input file (STL) (Required)
    -o, --output           Output filename without extension [default: out]
    -g, --grid             Grid size [default: 100]
    --thickness            Thickness [default: 1.0]
    --border_offset        default:2
    --coff                 default:4*PI
    --minimum_diameter     used for removing small orphaned (between 0-1) [default: 0.25]
```

## Dependencies
- [libigl](https://libigl.github.io/)
- [vcglib](https://github.com/cnr-isti-vclab/vcglib)

## How it works
- Read STL file and finding the boundary box
- Generate the grid and calculate the winding number with STL mesh
- Use winding number to determine the condition for [implicit isosurface function](https://wewanttolearn.wordpress.com/2019/02/03/triply-periodic-minimal-surfaces/)
- Generate the isosurface field in the same-size grid
- Perform [Dual marching cube](https://github.com/dominikwodniok/dualmc) to construct the manifold
- Clean up the duplicated vertices or faces, and abandon the group of connected faces having the diameter below the setting
- Export to the target 3D format

## Reference 
- [dualmc](https://github.com/dominikwodniok/dualmc)
- [argparse](https://github.com/jamolnng/argparse)
- [Minimal surface Blog](https://minimalsurfaces.blog/)
