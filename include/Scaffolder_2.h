#pragma once
#include <string>
#include <sstream>
#include <algorithm>
#include <ctime>

#include <igl/writePLY.h>
#include <igl/readSTL.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/fast_winding_number.h>
#include <Eigen/Core>

#include "cxxopts.hpp"
#include "diplib.h"
#include "diplib/file_io.h"
#include "diplib/regions.h"
#include "diplib/measurement.h"
#include "dualmc/dualmc.h"
#include "ProgressBar.hpp"
#include "OptimalSlice.hpp"

#include "implicit_function.h"
#include "Mesh.h"
#include "utils.h"

#define VERSION "v1.3"
#define PROGRESS_BAR_COLUMN 40

#define METHOD_IMAGE_PROCESS 0
#define METHOD_SLICE_CONTOUR 1

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
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
        }

        if (verbose) std::cout << "-- Info: " << mc_vertices.size() << " vertices " << mc_quads.size() << " faces" << std::endl;
    }
}

inline void clean_mesh(TMesh& mesh, double minimum_diameter, uint16_t smooth_step, bool verbose = true) {
    if (verbose) std::cout << "[libVCG Cleaning] ";
    vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
    vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::UpdateBounding<TMesh>::Box(mesh);
    vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::UpdateBounding<TMesh>::Box(mesh);
    if (verbose) std::cout << "OK" << std::endl;
    if (smooth_step > 0) {
        if (verbose) std::cout << "[Laplacian smoothing] ";
        vcg::tri::Smooth<TMesh>::VertexCoordLaplacian(mesh, smooth_step, false, true);
        if (verbose) std::cout << "OK" << std::endl;
    }
    vcg::tri::Clean<TMesh>::MergeCloseVertex(mesh, SLICE_PRECISION*1000);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
}

inline void fix_self_intersect_mesh(TMesh& mesh, double minimum_diameter, uint16_t max_iteration = 10, bool verbose = false) {
    
    std::vector<TMesh::FaceType*> faces;
    vcg::tri::Clean<TMesh>::SelfIntersections(mesh, faces);
    
    uint16_t iteration = 0;
    int maxSize = mesh.bbox.SquaredDiag();
    size_t nf = 0;
    while (iteration < max_iteration && (faces.size() > 0 || nf > 0)) {
        
        if (verbose) std::cout << "[Fix Self-intersect face]" << std::endl << "  -- Iteration " << iteration + 1 << std::endl;
        
        if (faces.size() > 0) {
            for (size_t i = 0; i < faces.size(); i++) {
                if (!faces[i]->IsD())
                    vcg::tri::Allocator<TMesh>::DeleteFace(mesh, *(faces[i]));
            }
            if (verbose) std::cout << "  -- self-intersect faces: " << faces.size() << std::endl;
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            vcg::tri::UpdateBounding<TMesh>::Box(mesh);
            if (verbose) std::cout << "  -- Remove faces [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            if (verbose) std::cout << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            if (verbose) std::cout << "  -- Remove non-manifold edges [OK]" << std::endl;
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                std::cout << "[Warning]: Fixed Self-intersecting failed: Mesh has some not 2-manifold edges" << std::endl;
                return;
            }
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            if (verbose) std::cout << "  -- Close holes [OK]" << std::endl;
        }
        
        vcg::tri::UpdateFlags<TMesh>::FaceBorderFromNone(mesh);
        vcg::tri::UpdateFlags<TMesh>::VertexBorderFromFaceBorder(mesh);
        vcg::tri::UpdateSelection<TMesh>::FaceFromBorderFlag(mesh);
        vcg::tri::UpdateSelection<TMesh>::VertexFromBorderFlag(mesh);
        nf = vcg::tri::UpdateSelection<TMesh>::FaceCount(mesh);
        if (nf > 0) {
            vcg::tri::UpdateSelection<TMesh>::VertexFromFaceLoose(mesh);
            vcg::tri::UpdateSelection<TMesh>::FaceFromVertexLoose(mesh);
            vcg::tri::UpdateSelection<TMesh>::VertexClear(mesh);
            vcg::tri::UpdateSelection<TMesh>::VertexFromFaceStrict(mesh);
            for (TMesh::FaceIterator it = mesh.face.begin(); it != mesh.face.end(); ++it) {
                if (!it->IsD() && it->IsS())
                    vcg::tri::Allocator<TMesh>::DeleteFace(mesh, *it);
            }
            for (TMesh::VertexIterator it = mesh.vert.begin(); it != mesh.vert.end(); ++it) {
                if (!it->IsD() && it->IsS())
                    vcg::tri::Allocator<TMesh>::DeleteVertex(mesh, *it);
            }
            vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            if (verbose) std::cout << "  -- Remove Border faces [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateBounding<TMesh>::Box(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            if (verbose) std::cout << "  -- Remove non-manifold edges [OK]" << std::endl;
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            if (verbose) std::cout << "  -- Close holes [OK]" << std::endl;

            vcg::tri::UpdateFlags<TMesh>::FaceBorderFromNone(mesh);
            vcg::tri::UpdateFlags<TMesh>::VertexBorderFromFaceBorder(mesh);
            vcg::tri::UpdateSelection<TMesh>::FaceFromBorderFlag(mesh);
            vcg::tri::UpdateSelection<TMesh>::VertexFromBorderFlag(mesh);
            nf = vcg::tri::UpdateSelection<TMesh>::FaceCount(mesh);
        }
        
        vcg::tri::Clean<TMesh>::SelfIntersections(mesh, faces);
        iteration++;
    }
    vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
    
    while (iteration < max_iteration && nf > 0) {
        if (verbose) std::cout << "[Fix Border Edges and holes]" << std::endl << "  -- Iteration " << iteration + 1 << std::endl;
        

        
        iteration++;
    }
    vcg::tri::UpdateBounding<TMesh>::Box(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
}

inline void report_mesh(TMesh& mesh) {
    int connectedComponentsNum = vcg::tri::Clean<TMesh>::CountConnectedComponents(mesh);
    std::cout
        << "[Topology Measurement] " << std::endl
        << "-- Mesh is composed by " << connectedComponentsNum << " connected component(s)" << std::endl;

    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
    vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);
    int vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(mesh, true);

    if (edgeNonManifoldNum == 0 && vertManifNum == 0) {
        int holeNum = vcg::tri::Clean<TMesh>::CountHoles(mesh);
        int genus = vcg::tri::Clean<TMesh>::MeshGenus(mesh.vn, edgeNum, mesh.fn, holeNum, connectedComponentsNum);

        std::cout
            << "-- Mesh is two-manifold " << std::endl
            << "-- Mesh has " << holeNum << " holes" << std::endl
            << "-- Genus is " << genus << std::endl;
    }
}

inline bool is_mesh_manifold(TMesh &mesh) {
    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
    vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);
    int vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(mesh, true);
    return edgeNonManifoldNum == 0 && vertManifNum == 0;
}