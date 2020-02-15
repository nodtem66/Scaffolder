#include <igl/writePLY.h>
#include <igl/readSTL.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/fast_winding_number.h>
#include <Eigen/Core>
#include <string>
#include <sstream>

#include "cxxopts.hpp"
#include "diplib.h"
#include "diplib/file_io.h"
#include "diplib/regions.h"
#include "diplib/measurement.h"
#include "dualmc/dualmc.h"
//#include "H5Easy.hpp"
#include "ProgressBar.hpp"

#include "implicit_function.h"
#include "Mesh.h"

#define VERSION "v1.2"
#define PROGRESS_BAR_COLUMN 40

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
    r.z = index / (grid_size(0) * grid_size(1));
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
        if (verbose) std::cout << "[Marching Cube] ";
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
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
        }

        if (verbose) std::cout << "OK" << std::endl;
        if (verbose) std::cout << "-- Info: " << mc_vertices.size() << " vertices " << mc_quads.size() << " faces" << std::endl;
    }
}

inline void clean_mesh(TMesh& mesh, double minimum_diameter, uint16_t smooth_step, bool verbose = true) {
    if (verbose) std::cout << "[libVCG Cleaning] ";
    vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::UpdateBounding<TMesh>::Box(mesh);
    vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    if (mesh.fn > 0) {
        vcg::tri::UpdateNormal<TMesh>::PerFaceNormalized(mesh);
        vcg::tri::UpdateNormal<TMesh>::PerVertexAngleWeighted(mesh);
    }
    if (verbose) {
        std::cout
            << "OK" << std::endl
            << "-- IsWaterTight: " << vcg::tri::Clean<TMesh>::IsWaterTight(mesh) << std::endl
            << "-- Minimum_diameter: " << minimum_diameter * mesh.bbox.Diag() << std::endl;
    }
    if (smooth_step > 0) {
        if (verbose) std::cout << "[Laplacian smoothing] ";
        vcg::tri::Smooth<TMesh>::VertexCoordLaplacian(mesh, smooth_step, false, true);
        if (verbose) std::cout << "OK" << std::endl;
    }
}

int make_dir(std::string& str) {
    if (str.empty())
        return 0;
#ifdef _WIN32
    wchar_t* wc = new wchar_t[str.size() + 1];
    mbstowcs(wc, str.c_str(), str.size());
    return _wmkdir(wc);
#else
    #include <stdio.h>
    #include <io.h>
    return mkdir(str.c_str(), 0733);
#endif
}