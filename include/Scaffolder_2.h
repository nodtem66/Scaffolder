#pragma once
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <numeric>

#include <igl/fast_winding_number.h>
#include <igl/signed_distance.h>
#include <Eigen/Core>

#include "cxxopts.hpp"
#include "dualmc/dualmc.h"

#include "ProgressBar.hpp"
#include "OptimalSlice.hpp"
#include "implicit_function.h"
#include "MeshOperation.h"
#include "utils.h"
#include "QuadricSimplification.h"
#include "sol/sol.hpp"

#define VERSION "v1.5.0-alpha"
#define PROGRESS_BAR_COLUMN 40

#define SCAFFOLDER_FORMAT_DEFAULT 0
#define SCAFFOLDER_FORMAT_CSV     1

#define LOG if (log_format == SCAFFOLDER_FORMAT_DEFAULT) log
#define CSV if (log_format == SCAFFOLDER_FORMAT_CSV) log

typedef struct index_type {
    size_t x; size_t y; size_t z;
} index_type;
typedef std::map<size_t, bool> Queue_t;

// Flatten between 1D and 3D
// https://stackoverflow.com/questions/7367770/how-to-flatten-or-index-3d-array-in-1d-array
inline size_t indexFromIJK(size_t i, size_t j, size_t k, Eigen::RowVector3i grid_size) {
    return i + grid_size(0) * (j + grid_size(1) * k);
}

inline void indexToIJK(size_t index, Eigen::RowVector3i grid_size, index_type& r) {
    r.z = index / grid_size(0) / grid_size(1);
    index -= r.z * grid_size(0) * grid_size(1);
    r.y = index / grid_size(0);
    r.x = index % grid_size(0);
}

inline bool MarkAndSweepNeighbor(Eigen::VectorXd& W, index_type& index, Queue_t& queue, Eigen::RowVector3i grid_size, double value = 0.0, bool findAbove = false) {
    bool isBorder = false;
    for (int8_t di = -1; di <= 1; di++) {
        for (int8_t dj = -1; dj <= 1; dj++) {
            for (int8_t dk = -1; dk <= 1; dk++) {
                if (di == 0 && dj == 0 && dk == 0) continue;
                const size_t id = indexFromIJK(index.x + di, index.y + dj, index.z + dk, grid_size);
                //std::cout << value << " " << W(id) << std::endl;
                if ((findAbove && W(id) >= value - eps2) || (!findAbove && W(id) <= value + eps2)) {
                    isBorder = true;
                    break;
                }
            }
            if (isBorder) break;
        }
        if (isBorder) break;
    }
    if (isBorder) {
        for (int8_t di = -1; di <= 1; di++) {
            for (int8_t dj = -1; dj <= 1; dj++) {
                for (int8_t dk = -1; dk <= 1; dk++) {
                    if (di == 0 && dj == 0 && dk == 0) continue;
                    const size_t id = indexFromIJK(index.x + di, index.y + dj, index.z + dk, grid_size);
                    if (W(id) >= 0.5 && W(id) < 1.1) {
                        queue.insert({ id, true });
                    }
                }
            }
        }
    }
    return isBorder;
}

inline void marching_cube(TMesh &mesh, Eigen::MatrixXd &Fxyz, Eigen::RowVector3i grid_size, Eigen::RowVector3d &Vmin, double delta, bool verbose = true, bool dirty = false) {
    {
        if (verbose) std::cout << "[Marching Cube] " << std::endl;
        dualmc::DualMC<double> builder;
        std::vector<dualmc::Vertex> mc_vertices;
        std::vector<dualmc::Quad> mc_quads;
        builder.build((double const*)Fxyz.data(), grid_size(0), grid_size(1), grid_size(2), 0, true, false, mc_vertices, mc_quads);
        TMesh::VertexIterator vi = vcg::tri::Allocator<TMesh>::AddVertices(mesh, mc_vertices.size());
        TMesh::FaceIterator fi = vcg::tri::Allocator<TMesh>::AddFaces(mesh, mc_quads.size() * 2);
        std::vector<TMesh::VertexPointer> vp(mc_vertices.size());
        for (size_t i = 0, len = mc_vertices.size(); i < len; i++, ++vi) {
            vp[i] = &(*vi);
            vi->P() = TMesh::CoordType(
                Vmin(0) + mc_vertices[i].x * delta,
                Vmin(1) + mc_vertices[i].y * delta,
                Vmin(2) + mc_vertices[i].z * delta
            );
        }
        for (size_t i = 0, len = mc_quads.size(); i < len; i++, ++fi) {
            fi->V(0) = vp[mc_quads[i].i0];
            fi->V(1) = vp[mc_quads[i].i1];
            fi->V(2) = vp[mc_quads[i].i2];
            ++fi;
            fi->V(0) = vp[mc_quads[i].i2];
            fi->V(1) = vp[mc_quads[i].i3];
            fi->V(2) = vp[mc_quads[i].i0];
        }
        if (!dirty) {
            vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
            vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
        }
    }
}

inline void mesh_to_eigen_vector(TMesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    V.resize(mesh.VN(), 3);
    size_t i = 0;
    std::vector<size_t> vertexId(mesh.vert.size());
    for (TMesh::VertexIterator it = mesh.vert.begin(); it != mesh.vert.end(); ++it) if (!it->IsD()) {
        vertexId[it - mesh.vert.begin()] = i;
        vcg::Point3d point = it->P();
        V(i, 0) = point[0];
        V(i, 1) = point[1];
        V(i, 2) = point[2];
        i++;
    }
    // Faces to Eigen matrixXi F1
    i = 0;
    F.resize(mesh.FN(), mesh.face.begin()->VN());
    for (TMesh::FaceIterator it = mesh.face.begin(); it != mesh.face.end(); ++it) if (!it->IsD()) {
        for (int k = 0; k < it->VN(); k++) {
            F(i, k) = vertexId[vcg::tri::Index(mesh, it->V(k))];
        }
        i++;
    }
}

ProgressBar qsim_progress(100, 40);
bool qsim_callback(int pos, const char* str) {
    if (pos >= 0 && pos <= 100) {
        qsim_progress.update(pos);
        qsim_progress.display();
    }
    if (pos >= 100)
        qsim_progress.done();
    return true;
}

void set_shorten_function(sol::state& lua) {
    lua.script("abs, acos, asin, atan, atan2 = math.abs, math.acos, math.atan, math.atan2");
    lua.script("ceil, cos, deg, exp, floor = math.ceil, math.cos, math.deg, math.exp, math.floor");
    lua.script("log, log10, max, min, mod = math.log, math.log10, math.max, math.min, math.mod");
    lua.script("pow, rad, sin, sqrt, tan = math.pow, math.rad, math.sin, math.sqrt, math.tan");
    lua.script("frexp, ldexp, random, randomseed = math.frexp, math.ldexp, math.random, math.randomseed");
    lua.script("local pi = math.pi");
}