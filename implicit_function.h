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

typedef double FT;
const FT pi = 3.141592653589793238462643383279502884L;
const FT esp = 1e-6;
typedef FT(Function_3)(FT, FT, FT);

inline FT schwarzp(FT x, FT y, FT z) {
	return cos(x) + cos(y) + cos(z);
}
inline FT sphere(FT x, FT y, FT z, FT r) {
	return x * x + y * y + z * z - r * r;
}
inline FT DFx(Function_3 F, FT x, FT y, FT z) {
	return (F(x + esp, y, z) - F(x, y, z)) / esp;
}
inline FT DFy(Function_3 F, FT x, FT y, FT z) {
	return (F(x, y + esp, z) - F(x, y, z)) / esp;
}
inline FT DFz(Function_3 F, FT x, FT y, FT z) {
	return (F(x, y, z + esp) - F(x, y, z)) / esp;
}
inline FT DF(FT dx, FT dy, FT dz) {
	return sqrt(dx * dx + dy * dy + dz * dz);
}
inline FT IsoThicken(Function_3 F, FT x, FT y, FT z, FT thickness = 1.0) {
	FT dx = DFx(schwarzp, x, y, z), dy = DFy(schwarzp, x, y, z), dz = DFz(schwarzp, x, y, z);
	FT const df = DF(dx, dy, dz);
	dx *= 0.5 * thickness / df;
	dy *= 0.5 * thickness / df;
	dz *= 0.5 * thickness / df;
	FT const iso1 = schwarzp(x + dx, y + dy, z + dz);
	FT const iso2 = schwarzp(x - dx, y - dy, z - dz);
	return iso1 * iso2;
}
inline FT max_3(FT a, FT b, FT c) {
	return (std::max)((std::max)(a, b), c);
}

class Function {
public:
	virtual FT operator()(FT, FT, FT) = 0;
};
/* Function reture the isocuboid boundary
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
/* Function return schwarzp isosurface function with following parameters:
 * double thickness: adding the boundary by +-(thickness/2) to the isosurface function
 */
class Implicit_function : public Function {
private:
	const Function_3& isosurface;
	const double coff;
	const double thickness;
	//Function& condition;
public:
	Implicit_function(const Function_3& isosurface, const double coff, const double thickness) :
		isosurface(isosurface), coff(coff), thickness(thickness) {}
	FT operator()(FT x, FT y, FT z) {
		/*
		if (condition(x, y, z) <= 0) {
			return IsoThicken(isosurface, x * coff, y * coff, z * coff, thickness);
		}
		return 1.0;
		*/
		if (thickness <= esp) {
			return isosurface(x * coff, y * coff, z * coff);
		}
		return IsoThicken(isosurface, x * coff, y * coff, z * coff, thickness);
	}
};

void to_lower(std::string& s) {
	std::locale loc;
	for (std::string::size_type i = 0; i < s.length(); ++i)
		s[i] = std::tolower(s[i], loc);
}
#endif