import PyScaffolder
import unittest

from math import sqrt

class TestPyScaffolder(unittest.TestCase):

	v = [
		[0.0, 0.0, 0.0],
		[1.0, 0.0, 0.0],
		[0.0, -1.0, 0.0],
		[1.0, -1.0, 0.0],
		[0.0, 0.0, 1.0],
		[1.0, 0.0, 1.0],
		[0.0, -1.0, 1.0],
		[1.0, -1.0, 1.0]
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

	counter = 0

	def _callback(self, p):
		self.counter = p

	def test_version(self):
		self.assertIsNotNone(PyScaffolder.__version__)

	def test_slicer(self):
		a = PyScaffolder.slice_test(self.v, self.f, direction=3)
		self.assertEqual(len(self.v), 8)
		self.assertEqual(len(self.f), 12)
		self.assertEqual(len(a.minFeret), 0)
		self.assertEqual(len(a.maxFeret), 0)

	def test_slicer_with_callback(self):
		a = PyScaffolder.slice_test(self.v, self.f, direction=3, callback=self._callback)
		self.assertEqual(self.counter, 100)

	def test_scaffolder(self):			
		# Test generate surface with default parameter
		params = PyScaffolder.Parameter()
		params.coff = 12.0
		params.verbose = False
		a = PyScaffolder.generate_scaffold(self.v, self.f, params)
		self.assertGreater(a.porosity, 0)
		self.assertGreater(a.surface_area_ratio, 0)
		self.assertGreater(len(a.v), 0)
		self.assertEqual(len(a.v[0]), 3)
		self.assertGreater(len(a.f), 0)
		self.assertEqual(len(a.f[0]), 3)
	
	def test_scaffolder_with_callback(self):
		params = PyScaffolder.Parameter()
		params.coff = 12.0
		PyScaffolder.generate_scaffold(self.v, self.f, params, callback=self._callback)
		self.assertEqual(self.counter, 100)


	def test_marching_cubes(self):
		Fxyz = []
		N = 100
		for i in range(N):
			for j in range(N):
				for k in range(N):
					Fxyz.append(sqrt((i-50)**2+(j-50)**2+(k-50)**2) - 2**2)
		self.assertEqual(len(Fxyz), N**3)
		(v, f) = PyScaffolder.marching_cubes(Fxyz, grid_size=N, delta=0.02, v_min=(-1, -1, -1), clean=False)
		self.assertGreater(len(v), 0)
		self.assertGreater(len(f), 0)
		self.assertEqual(len(v[0]), 3)
		self.assertEqual(len(f[0]), 3)

	def test_marching_cubes_with_callback(self):
		Fxyz = [
			1,1,1,1,
			1,-1,-1,1,
			1,-1,-1,1,
			1,1,1,1
		]*4
		self.counter = 0
		(v, f) = PyScaffolder.marching_cubes(Fxyz, grid_size=4, delta=0.25, v_min=(-.5, -.5, -.5), clean=True, callback=self._callback)
		self.assertGreater(len(v), 0)
		self.assertGreater(len(f), 0)
		self.assertEqual(self.counter, 100)
		

if __name__ == '__main__':
    unittest.main()