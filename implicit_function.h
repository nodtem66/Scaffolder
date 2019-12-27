// Scaffolder_2.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#ifndef IMPLICIT_FUNCTION_H
#define IMPLICIT_FUNCTION_H

#include <iostream>
#include <cmath>
#include <string>
#include <locale>
#include <algorithm>
#include <functional>

typedef double FT;
const FT pi = 3.141592653589793238462643383279502884L;
const FT eps = 1e-6;
const FT eps2 = 1e-2;
//typedef FT(Function_3)(FT, FT, FT);

class Function {
public:
	virtual FT operator()(FT, FT, FT) = 0;
};

/* Function return rectlinear style from slicer
 * where coff is required to adjust the diameter filament
 */
class Rectlinear : public Function {
private:
	const FT coff;
	const FT smooth_coff;
public:
	Rectlinear(FT coff = 1.0, FT smooth_coff = 0.05) : coff(coff), smooth_coff(smooth_coff) {}
	FT operator()(FT x, FT y, FT z) {
		return ((cos(x) + cos(y)) - 0.25 * coff) * ((cos(y + coff) + cos(z + coff)) - 0.25 * coff) - smooth_coff;
	}
};

class Schwarzp : public Function {
public:
	Schwarzp() {}
	FT operator ()(FT x, FT y, FT z) {
		return cos(x) + cos(y) + cos(z);
	}
};

class Schwarzd : public Function {
public:
	Schwarzd() {}
	FT operator ()(FT x, FT y, FT z) {
		return sin(x) * sin(y) * sin(z) + sin(x) * cos(y) * cos(z) + cos(x) * sin(y) * cos(z) + cos(x) * cos(y) * sin(z);
	}
};

class Gyroid : public Function {
public:
	Gyroid() {}
	FT operator ()(FT x, FT y, FT z) {
		return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x);
	}
};

class Lidinoid : public Function {
public:
	Lidinoid() {}
	FT operator ()(FT x, FT y, FT z) {
		return 0.5 * (sin(2 * x) * cos(y) * sin(z) + sin(2 * y) * cos(z) * sin(x) + sin(2 * z) * cos(x) * sin(y)) - 0.5 * (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x)) + 0.15;
	}
};

class Neovius : public Function {
public:
	Neovius() {}
	FT operator ()(FT x, FT y, FT z) {
		return 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z);
	}
};

class Schoen_iwp : public Function {
public:
	Schoen_iwp() {}
	FT operator ()(FT x, FT y, FT z) {
		return cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x) - cos(x) * cos(y) * cos(z);
	}
};
class PWHybrid : public Function {
public:
	PWHybrid() {}
	FT operator ()(FT x, FT y, FT z) {
		return 4 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x)) - 3 * cos(x) * cos(y) * cos(z);
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
public:
	Implicit_function(Function* isosurface, const double coff, const double thickness):
		isosurface(*isosurface), coff(coff), thickness(thickness) {}
	FT operator ()(FT x, FT y, FT z) {
		if (thickness <= eps)
			return (isosurface)(x * coff, y * coff, z * coff);
		return IsoThicken(isosurface, x * coff, y * coff, z * coff, thickness);
	}
};

inline void to_lower(std::string& s) {
	std::locale loc;
	for (std::string::size_type i = 0; i < s.length(); ++i)
		s[i] = std::tolower(s[i], loc);
}

Function* isosurface(std::string name, FT coff) {
	if (name == "rectlinear") {
		return new Rectlinear(coff);
	}
	else if (name == "schwarzd") {
		return new Schwarzd();
	}
	else if (name == "lidinoid") {
		return new Lidinoid();
	}
	else if (name == "gyroid") {
		return new Gyroid();
	}
	else if (name == "schwarzp") {
		return new Schwarzp();
	}
	else if (name == "pwhybrid") {
		return new PWHybrid();
	}
	else if (name == "neovius") {
		return new Neovius();
	}
	else if (name == "schoen_iwp") {
		return new Schoen_iwp();
	}
	throw std::runtime_error(name + " doesn't exist");
}
#endif