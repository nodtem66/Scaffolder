#ifndef SCAFFOLDER_INCLUDED
#define SCAFFOLDER_INCLUDED

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <numeric>
#include <functional>

#include <igl/fast_winding_number.h>
#include <igl/signed_distance.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <Eigen/Core>

#include "Mesh.h"
#include "sol/sol.hpp"
#include "ProgressBar.hpp"

#define VERSION "v1.5.2"
#define PROGRESS_BAR_COLUMN 40

#define SCAFFOLDER_FORMAT_DEFAULT 0
#define SCAFFOLDER_FORMAT_CSV     1


typedef struct index_type {
    size_t x; size_t y; size_t z;
} index_type;
typedef std::map<size_t, bool> Queue_t;

inline size_t indexFromIJK(size_t i, size_t j, size_t k, Eigen::RowVector3i grid_size);
inline void indexToIJK(size_t index, Eigen::RowVector3i grid_size, index_type& r);
inline bool MarkAndSweepNeighbor(
    Eigen::VectorXd& W, 
    index_type& index, 
    Queue_t& queue, 
    Eigen::RowVector3i grid_size, 
    double value = 0.0, 
    bool findAbove = false
);

void marching_cube(
    TMesh& mesh, 
    Eigen::MatrixXd& Fxyz, 
    Eigen::RowVector3i grid_size,
    Eigen::RowVector3d& Vmin, 
    double delta, 
    bool verbose = true, 
    bool dirty = false
);
bool null_callback(int pos, const char* str);

bool qsim_callback(int pos, const char* str);
extern ProgressBar qsim_progress;

#endif