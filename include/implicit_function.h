// Scaffolder_2.h : Include file for standard system include files,
// or project specific include files.
#ifndef IMPLICIT_FUNCTION_H
#define IMPLICIT_FUNCTION_H

#include <iostream>
#include <cmath>
#include <string>
#include <locale>
#include <algorithm>
#include <functional>
#include <sol/sol.hpp>

typedef double FT;
const FT pi = 3.141592653589793238462643383279502884L;
const FT eps = 1e-6;
const FT eps2 = 1e-2;
//typedef FT(Function_3)(FT, FT, FT);

class Function {
public:
	bool isLuaFunction = false;
	Function() {}
	virtual FT operator()(FT, FT, FT) = 0;
};

class LuaFunction : public Function {
private:
	sol::function* fn;
public:
	LuaFunction(sol::function& f) : fn(&f) { isLuaFunction = true; }
	FT operator()(FT x, FT y, FT z) {
		return (FT) (*fn)(x, y, z);
	}
};

class Fixed : public Function {
public:
	const FT val;
	Fixed(FT val) : val(val) {}
	FT operator()(FT x, FT y, FT z) {
		return val;
	}
};
/* Function return rectlinear style from slicer
 * where coff is required to adjust the diameter filament
 */
class Rectlinear : public Function {
private:
	const FT smooth_coff;
public:
	Rectlinear(FT smooth_coff = 0.05) : smooth_coff(smooth_coff) {}
	FT operator()(FT x, FT y, FT z) {
		return ((cos(x) + cos(y)) - 1.2) * ((cos(y + pi) + cos(z + pi)) - 1.2) - smooth_coff;
	}
};

/*
 * TPMS from http://www.msri.org/publications/sgp/jim/papers/morphbysymmetry/table/index.html
 */
class Schwarzp : public Function {
public:
	const FT t;
	Schwarzp(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return cos(x) + cos(y) + cos(z) + t;
	}
};

class DoubleP : public Function {
public:
	const FT t;
	DoubleP(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return cos(x) * cos(y) + cos(y) * cos(z) + cos(x) * cos(z) + 0.35 * (cos(2 * x) + cos(2 * y) + cos(2 * z)) + t;
	}
};

class Schwarzd : public Function {
public:
	const FT t;
	Schwarzd(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return sin(x) * sin(y) * sin(z) + sin(x) * cos(y) * cos(z) + cos(x) * sin(y) * cos(z) + cos(x) * cos(y) * sin(z) + t;
	}
};

class DoubleD : public Function {
public:
	const FT t;
	DoubleD(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return -1 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(x) * cos(z)) - 1 * (sin(x) * sin(y) * sin(z)) + t;
	}
};

class Gyroid : public Function {
public:
	const FT t = 0;
	Gyroid(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x) + t;
	}
};

class DoubleGyroid : public Function {
public:
	const FT t = 0;
	DoubleGyroid(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return -1 * (
			2.75 * (sin(2 * x) * sin(z) * cos(y) + sin(2 * y) * sin(x) * cos(z) + sin(2 * z) * sin(y) * cos(x)) -
			(cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * x) * cos(2 * z))
		) + t;
	}
};

class Lidinoid : public Function {
public:
	const FT t = 0;
	Lidinoid(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return (
			(sin(2 * x) * cos(y) * sin(z) + sin(2 * y) * cos(z) * sin(x) + sin(2 * z) * cos(x) * sin(y)) +
			(cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x))
		) + t;
	}
};

class Neovius : public Function {
public:
	const FT t;
	Neovius(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z) + t;
	}
};

class Schoen_iwp : public Function {
public:
	const FT t;
	Schoen_iwp(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x) + t;
	}
};

class BCC : public Function {
public:
	const FT t;
	BCC(FT t = 0) : t(t) {}
	FT operator ()(FT x, FT y, FT z) {
		return cos(x) + cos(y) + cos(z) - 2 * (cos(0.5 * x) * cos(0.5* y) + cos(0.5 * y) * cos(0.5 * z) + cos(0.5 * x) * cos(0.5 * z)) + t;
	}
};

class TGab : public Function {
public:
	TGab() {}
	FT operator ()(FT x, FT y, FT z) {
		return 20 * (cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x)) - 0.5 * (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x)) - 4;
	}
};

class TGc : public Function {
public:
	TGc() {}
	FT operator ()(FT x, FT y, FT z) {
		return -(10 * (cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x)) - 2 * (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x)) - 12);
	}
};

template<typename Function>
inline FT DFx(Function& F, FT x, FT y, FT z) {
	return (F(x + eps, y, z) - F(x, y, z)) / eps;
}
template<typename Function>
inline FT DFy(Function& F, FT x, FT y, FT z) {
	return (F(x, y + eps, z) - F(x, y, z)) / eps;
}
template<typename Function>
inline FT DFz(Function& F, FT x, FT y, FT z) {
	return (F(x, y, z + eps) - F(x, y, z)) / eps;
}
inline FT DF(FT dx, FT dy, FT dz) {
	return sqrt(dx * dx + dy * dy + dz * dz);
}
inline FT IsoThicken(Function& F, FT x, FT y, FT z, FT thickness = 1.0) {
	FT dx = DFx(F, x, y, z), dy = DFy(F, x, y, z), dz = DFz(F, x, y, z);
	FT const df = DF(dx, dy, dz);
	dx *= 0.5 * thickness / df;
	dy *= 0.5 * thickness / df;
	dz *= 0.5 * thickness / df;
	FT const iso1 = F(x + dx, y + dy, z + dz);
	FT const iso2 = F(x - dx, y - dy, z - dz);
	return iso1 * iso2;
}
inline FT max_3(FT a, FT b, FT c) {
	return (std::max)((std::max)(a, b), c);
}

/* Function return the isocuboid boundary
 * where minimum point origins at (origin_x, origin_y, origin_z) and 
 * maximum point is (origin_x+w, origin_y+h, origin_z+t)
 * return 0 if the point on the surface
 * return -1 if the point inside or 1 if outside
 */
class Iso_cuboid_condition : public Function {
private:
	const FT origin_x;
	const FT origin_y;
	const FT origin_z;
	const FT w;
	const FT h;
	const FT t;
public:
	Iso_cuboid_condition(FT x = 0, FT y = 0, FT z = 0, FT w = 1, FT h = 1, FT t = 1) :
		origin_x(x), origin_y(y), origin_z(z), w(w), h(h), t(t) {}
	FT operator()(FT x, FT y, FT z) {
		const FT dx = 0.5 * w, dy = 0.5 * h, dz = 0.5 * t;
		if (max_3(
			abs(x - origin_x - dx) / dx,
			abs(y - origin_y - dy) / dy,
			abs(z - origin_z - dz) / dz
		) > 1)
			return 1.001;
		return -.001;
	}
};
/* Function return Implicit isosurface function with following parameters:
 * double thickness: adding the boundary by +-(thickness/2) to the isosurface function
 */
class Implicit_function : public Function {
private:
	Function& isosurface;
	const double coff;
	const double thickness;
	const bool isLuaFunction;
public:
	Implicit_function(Function* isosurface, const double coff, const double thickness = 0):
		isosurface(*isosurface), coff(coff), thickness(thickness), isLuaFunction(isosurface->isLuaFunction) {}
	FT operator ()(FT x, FT y, FT z) {
		// Since we use FREP where F(x,y,z) >= 0 defined as a solid
		// then we inversed the implicit function to match FREP
		// For example, f(x,y,z) = x^2+y^2+z^2 - 1, We want the solid region of f(x, y, z) <= 0
		// so that F(x,y,z) = -f(x,y,z) >= 0
		if (thickness <= eps) {
			if (isLuaFunction) {
				return (isosurface)(x, y, z);
			}
			return -(isosurface)(x * coff, y * coff, z * coff);
		}
		return IsoThicken(isosurface, x * coff, y * coff, z * coff, thickness);
	}
};

Function* isosurface(std::string name, FT t) {
	if (name == "empty") {
		return new Fixed(1);
	}
	else if (name == "rectlinear") {
		return new Rectlinear();
	}
	else if (name == "schwarzp") {
		return new Schwarzp(t);
	}
	else if (name == "double-p") {
		return new DoubleP(t);
	}
	else if (name == "schwarzd") {
		return new Schwarzd(t);
	}
	else if (name == "double-d") {
		return new DoubleD(t);
	}
	else if (name == "gyroid") {
		return new Gyroid(t);
	}
	else if (name == "double-gyroid") {
		return new DoubleGyroid(t);
	}
	else if (name == "lidinoid") {
		return new Lidinoid(t);
	}
	else if (name == "neovius") {
		return new Neovius(t);
	}
	else if (name == "schoen_iwp") {
		return new Schoen_iwp(t);
	}
	else if (name == "bcc") {
		return new BCC(t);
	}
	else if (name == "tubular_g_ab") {
		return new TGab();
	}
	else if (name == "tubular_g_c") {
		return new TGc();
	}
	throw std::runtime_error("Implicit surface called `" + name + "` doesn't exist");
}
#endif