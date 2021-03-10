#pragma once

#ifndef PYSCAFFOLDER_INCLUDED_H
#define PYSCAFFOLDER_INCLUDED_H

#include <iostream>
#include <vector>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"
#include "pybind11/eigen.h"
#include "pybind11/functional.h"

namespace PyScaffolder {

	struct PoreSize {
		PoreSize() {}
		Eigen::RowVectorXd minFeret;
		Eigen::RowVectorXd maxFeret;
	};

	struct MeshInfo {
		Eigen::MatrixXd v;
		Eigen::MatrixXi f;
		double porosity;
		double surface_area;
		double surface_area_ratio;
	};

	struct Parameter {
		bool is_build_inverse = false;
		bool is_intersect = true;
		uint16_t shell = 0;
		uint16_t grid_offset = 5;
		uint16_t smooth_step = 5;
		uint16_t k_slice = 100;
		uint16_t k_polygon = 4;
		uint16_t fix_self_intersect = 0;
		size_t grid_size = 100;
		double isolevel = 0.0;
		double qsim_percent = 0;
		double coff = 3.141592653589793238462643383279502884L;
		double minimum_diameter = 0.25;
		std::string surface_name = "bcc";
	};

	PoreSize slice_test(
		Eigen::MatrixXd v, 
		Eigen::MatrixXi f, 
		size_t k_slice = 100, 
		size_t k_polygon = 4, 
		int direction = 0, 
		const std::function<void(int)>& callback = NULL
	);

	MeshInfo generate_mesh(
		Eigen::MatrixXd v,
		Eigen::MatrixXi f,
		Parameter params,
		const std::function<void(int)>& callback = NULL
	);
}
#endif