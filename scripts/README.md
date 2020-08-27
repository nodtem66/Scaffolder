# This is experiment folder of scaffold generation and analysis
- Grid size: 300 best
- coff: depend on surface name and the length of pore
- thickness:  depend on the surface name and porosity
- Mechanical stress/stain: depend on the surface name, coff, thickness

## tetgen usage
from [Tetgen Manual](http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html)

- gBNEF export only `medit mesh` file
- p/0.01 decrease the dihedral angle to 0.01
- q1.414/0 set the radius-ratio
- R enable to mesh optimization

**Example:**
```
$ tetgen -gBNEFp/0.01q1.4R file.stl
// output will be file.1.mesh
```

## ElmerGrid usage
after .mesh file was generated by `tetgen` use this command to generate elmerSolver mesh
```
$ ElmerGrid 12 2 file -autoclean ...
// note that `file` not `file.mesh` 
```