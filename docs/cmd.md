# Command line
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
```

## Examples

- Generated BCC scaffold with $w=\pi$, $t=0$, and $grid\_size=100$
```bash
  Scaffolder input.stl output.stl bcc,3.14159,0,100
```

- Generated [custom implicit function](#custom-implicit-function) (declaired in `custom.lua`) with $w=\pi$, $t=0$,
$grid\_size=100$, $k_{slice}=100$, $k_{polygon}=4$ and evaluated scaffold properties
```bash
  Scaffolder input.stl output.stl custom.lua,3.14159,0,100,100,4 -m
```

- Generated and evaluated BCC scaffold ($w=3.14159$, $t=0$) and reported in CSV
```bash
Scaffolder input.stl output.stl bcc,3.14159,0 -m -q --format csv
```

## Custom implicit function

If the built-in functions (`rectlinear, schwarzp, schwarzd, gyroid, double-p, double-d, double-gyroiod, lidinoid,
schoen_iwp, neovius, bcc, tubular_g_ab, tubular_g_c`) was not satisfied, the custom implicit function can be used
by creating a lua file and define the **"surface" function that return the [FRep](https://en.wikipedia.org/wiki/Function_representation)**
  
```lua
-- custom.lua 
-- define gyroid with
function surface (x, y, z)
    return -(sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x) - params.isolevel)
end
```

The special symbols can be used in lua file:
```lua
  params = { coff, isolevel, k_splice, k_polygon }
  bbox = { min, max, length, grid_density }
  winding(x,y,z) -- function returning the winding number of position x,y,z
  signed_distance(x,y,z) -- function returning signed distance of position x,y,z
```
and also all functions in [math module](https://www.lua.org/manual/5.1/manual.html#5.6)