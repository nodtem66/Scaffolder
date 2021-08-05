# Gyroid Sphere

In this tutorial, we will generate gyroid scaffold from [`sphere.stl`](stl/sphere.stl) whose bounding box is
$2 \times 2 \times 2$. In order to generate porous mesh, the implicit function (`gyroid` for this case), 
angular frequency ($w$), and isolevel ($t$) are required (See [Command line](cmd.md)).

## Angular frequency

According to [prior study](https://github.com/nodtem66/Scaffolder/blob/master/data/data_visualization.ipynb),
**$w$ is inversely proportational to the pore size** (see [Tips](tips.md)). That is when $w$ increases,
the pore size will be decreased. For gyroid, $w \approx \pi/\text{pore size}$. If we want pore size of 0.25
(unit as same as `sphere.stl`), we try to initially set $w = \pi/0.25 \approx 12$ 

## Isolevel

**$t$ is directly proportional to the porosity**. In the other word, if we want to increase porosity, we have to increase $t$!.
However, some values of $t$ are probably lead to the problematic 3D mesh. Thus we initially leave this value to the default ($t=0$).

## Grid size

Grid size is a resolution of the output 3D mesh. Unless we generate a 3D mesh for FEA, this parameter is leave as a default value of 100.
!!! note "Parameters summary"
    input file = `sphere.stl`, surface = `gyroid`, coff ($w$) = `12`, isolevel ($t$) = `0`, grid size = `100`

## Command line
```bash
Scaffolder sphere.stl out.stl gyroid,12,0,100 -m
```
We add `-m` or `--microstructure` to enable the pore size evaluation. Then the console output is reported as the following:
```linenums="1" hl_lines="47 48 49 50 59"
[Scaffolder v1.5.1]
-- Input file: .\sphere.stl
-- Output file: out.stl
-- Surface (-n): gyroid
-- Coff (-c): 12
-- Isolevel (-t): 0 
-- Grid size (-g): 100
--   Grid offset: 3
--   Shell: 0
-- Autoclean: True
--   Minimum diameter: 25%
--   Smooth step: 5
--   Fix self-intersect: False
--   Quadric Simplification (-z): 0
-- Analysis microstructure (-m): True
--   Slice grid (k_slice): 100
--   Nearest outer contours (k_polygon): 4
--   Export microstructure: False
--   Mean curvature: False
-- Export JPEG: No
--   Axis: X
-- Build: Yes(Intersect)
-- Bounding Box: [-0.9995, -0.9995, -0.9995] [0.9995, 0.9995, 0.9995]
-- Length: [1.999, 1.999, 1.999]
-- Grid delta: 0.01999
[Generating grid] OK
-- Grid size: 1191016 [106, 106, 106]
[Calculating Winding number] OK
-- Sign Distance: [ -0.834758, 0.998785]  Wind: [ -0.00238368, 1.00151]
[Generating isosurface Fxyz] OK
[Marching Cube]
[========================================] 100% 0.194s
-- Info: 118708 vertices 237636 faces
[libVCG Cleaning] OK
-- is_manifold: 1
[Topology Measurement]
-- Mesh is composed by 1 connected component(s)
-- border edge: 0
-- non-manifold edge: 0
-- non-manifold vertex: 0
-- self-intersect faces: 1831
-- Mesh is two-manifold
-- Mesh has 0 holes
-- Genus is 107
[=======================================>] 98% 2.1s7s
[Microstructure]
-- Avg Min Feret: 0.0893954
-- Avg Max Feret: 0.606401
-- Min Feret: [0.00073122 0.033939 0.105832 0.131947 0.615222]
-- Max Feret: [0.00132485 0.243246 0.630859 0.723238 1.91447]
-- Square Similarity: [0 0 0 0 0]
-- Circle Similarity: [0 0 0 0 0]
-- Triangle Similarity: [0 0 0 0 0]
-- Ellipse Similarity: [0 0 0 0 0]
-- Elongation Similarity: [0 0 0 0 0]
[Scaffold properties]
-- Volume: 2.04856
-- Surface Area: 29.5085
-- Porosity: 0.508139
-- Surface Area ratio: 2.35532
[Writing file] OK
[Finished]
```

## Adjust isolevel for 60% porosity

Suppose that the target porosity is 60%. Then we increase $t$ to 1 and generate 3D mesh again

```bash
Scaffolder sphere.stl out.stl gyroid,12,1,100 -m -q
```

This time, we add `-q` or `--quite` to prevent the verbose output. The program will save the report in a file (`<output_name>_<surface><timestamp>.txt`) instead.
We use binary search to adjust the new value of $t$ to find the optimal value with the 60% porosity. 

| $t$ | 0 | 1 | 0.5 | 0.25 | 0.375 | 0.3125 | 0.28 | 0.265
| -- | -- | -- | -- | -- | -- | -- | -- | -- |
| Porosity | 0.508 | 0.85 | 0.67 | 0.59 | 0.634 | 0.613 | 0.602 | 0.597

We can also use the classic numerical method called [Secant method](https://en.wikipedia.org/wiki/Secant_method) to find the optimal $t$ for the target porosity.

$$
F(t) = Porosity(t) - 0.6 = 0
$$

$$
t_{next} = t_{current} - F(t_{current})\frac{t_{current} - t_{old}}{F(t_{current}) - F(t_{old})}
$$

$$
t_{next} = t_{current} - (\text{Porosity}_{current} - 0.6)\frac{t_{current} - t_{old}}{\text{Porosity}_{current} - \text{Porosity}_{old}}
$$

| $t$ | 0 | 1 | 0.27 
| -- | -- | -- | -- | 
| Porosity | 0.508 | 0.85 | 0.5989