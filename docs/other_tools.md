# Other tools

MSLattice, FLatt Pack, and Scaffolder were used to generate the 2 units of Gyroid (20mm×20mm×20mm).

The porosity was set to 70% (30% volume fraction).

MSLattice and FLatt Pack defined mesh density or discretization per unit cell, whereas Scaffolder defined grid size for entire structure.

Hence, the mesh density points and discretization were set to 100, while Scaffolder’s grid size was set to 200, 

the following relation between the number of unit cells (N), mesh length (L), unit cell size (l), and angular frequency (ω) is:
$$$
w=2\piNL=2\pil
$$$

Parameters used to generate Gyroid meshes
Program
Parameter
MSLattice
Uniform TPMS Lattices with solid network.
Architecture = Gyroid, Relative density = 30, Unit cell size = 10
Sample Length = 20, Sample Width = 20, Sample Height = 20
Mesh Density Points = 100
FLatt Pack
Geometry = Cuboid, Dimensions (x,y,z) = (20,20,20), Cells (x,y,z) = (2,2,2),
Phase=Network (Gyroid), Discretisation = 100,
Volume fraction =  Uniform (0.3), No skin
Scaffolder
ω=π/5, t=0.61, grid_size = 200
(Scaffolder cube20mm.stl out.stl gyroid,0.6283,0.61,200
--fix_self_intersect --intersect=false)



