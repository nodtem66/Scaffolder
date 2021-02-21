import PyScaffolder

# Load STL mesh with igl
v = [
	[0.0, 0.0, 0.0],
	[20.0, 0.0, 0.0],
	[0.0, -20.0, 0.0],
	[20.0, -20.0, 0.0],
	[0.0, 0.0, 20.0],
	[20.0, 0.0, 20.0],
	[0.0, -20.0, 20.0],
	[20.0, -20.0, 20.0]
]

f = [
	[1, 2, 0],
	[2, 1, 3],
	[7, 2, 3],
	[2, 7, 6],
	[1, 7, 3],
	[7, 1, 5],
	[7, 4, 6],
	[4, 7, 5],
	[4, 2, 6],
	[2, 4, 0],
	[4, 1, 0],
	[1, 4, 5]
]

# Test Slice test with all direction
a = PyScaffolder.slice_test(v, f, direction=3)
print("Vertices: ", len(v))
print("Faces: ", len(f))
print(len(a.minFeret))

# Test generate surface with default parameter
params = PyScaffolder.Parameter()
params.coff = 1.0
a = PyScaffolder.generate_mesh(v, f, params)
print(a.porosity)
print(a.surface_area_ratio)
print(len(a.v), len(a.f))