#pragma once
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/halfedge_indexed.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/hole.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#ifndef SLICE_PRECISION
#define SLICE_PRECISION 1e-8
#endif

class TFace;
class TVertex;

struct TUsedTypes : public vcg::UsedTypes< vcg::Use<TVertex>::AsVertexType, vcg::Use<TFace>::AsFaceType > {};

class TVertex : public vcg::Vertex< TUsedTypes,
    vcg::vertex::BitFlags,
    vcg::vertex::Coord3d,
    vcg::vertex::Normal3d,
    vcg::vertex::Mark > {};

class TFace : public vcg::Face<TUsedTypes,
    vcg::face::VertexRef,	// three pointers to vertices
    vcg::face::Normal3d,		// normal
    vcg::face::BitFlags,		// flags
    vcg::face::Mark,
    vcg::face::FFAdj			// three pointers to adjacent faces
> {};

/* the mesh is a container of vertices and a container of faces */
class TMesh : public vcg::tri::TriMesh< vector<TVertex>, vector<TFace> > {};

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
    vcg::tri::Clean<TMesh>::MergeCloseVertex(mesh, SLICE_PRECISION * 1000);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
}

inline bool fix_non_manifold(TMesh& mesh, double minimum_diameter, uint16_t max_iteration = 5, bool verbose = false) {
    size_t nf = 1;
    int maxSize = mesh.bbox.SquaredDiag();

    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
    vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);

    for (uint16_t iteration = 0; (edgeBorderNum > 0 || edgeNonManifoldNum > 0) && iteration < max_iteration; iteration++) {
        if (verbose) std::cout << "[Fix non-manifold]" << std::endl << "  -- Iteration " << iteration + 1 << std::endl;
        if (verbose) std::cout << "  -- non-manifold edges: " << edgeNonManifoldNum << std::endl;
        if (edgeNonManifoldNum > 0) {
            // fix floating faces
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            if (verbose) std::cout << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            // fix non-manifold edges
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            if (verbose) std::cout << "  -- Remove non-manifold edges [OK]" << std::endl;
            // fix holes
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                std::cout << "[Warning]: Fixed Self-intersecting failed: Mesh has some not 2-manifold edges" << std::endl;
                return false;
            }
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            if (verbose) std::cout << "  -- Close holes [OK]" << std::endl;
        }
        if (verbose) std::cout << "  -- border edges: " << edgeBorderNum << std::endl;
        if (edgeBorderNum > 0) {
            // select border vertices and faces
            vcg::tri::UpdateFlags<TMesh>::FaceBorderFromNone(mesh);
            vcg::tri::UpdateFlags<TMesh>::VertexBorderFromFaceBorder(mesh);
            vcg::tri::UpdateSelection<TMesh>::FaceFromBorderFlag(mesh);
            vcg::tri::UpdateSelection<TMesh>::VertexFromBorderFlag(mesh);
            // Dilate selection
            vcg::tri::UpdateSelection<TMesh>::VertexFromFaceLoose(mesh);
            vcg::tri::UpdateSelection<TMesh>::FaceFromVertexLoose(mesh);
            vcg::tri::UpdateSelection<TMesh>::VertexClear(mesh);
            vcg::tri::UpdateSelection<TMesh>::VertexFromFaceStrict(mesh);
            // delete all selected
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
            // fix floating faces
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            if (verbose) std::cout << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            // fix non-manifold edges
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            if (verbose) std::cout << "  -- Remove non-manifold edges [OK]" << std::endl;
            // fix holes
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                std::cout << "[Warning]: Fixed Self-intersecting failed: Mesh has some not 2-manifold edges" << std::endl;
                return false;
            }
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            if (verbose) std::cout << "  -- Close holes [OK]" << std::endl;
        }

        vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);
    }

    if (edgeBorderNum > 0 || edgeNonManifoldNum > 0) return false;
    return true;
}

inline bool fix_self_intersect_mesh(TMesh& mesh, double minimum_diameter, uint16_t max_iteration = 10, bool verbose = false) {

    std::vector<TMesh::FaceType*> faces;
    std::vector<TMesh::FaceType*>::iterator fc;
    vcg::tri::Clean<TMesh>::SelfIntersections(mesh, faces);

    uint16_t iteration = 0;
    int maxSize = mesh.bbox.SquaredDiag();
    size_t nf = 0, nv = 0;
    while (iteration < max_iteration && (faces.size() > 0 || nf > 0)) {

        if (verbose) std::cout << "[Fix Self-intersect face]" << std::endl << "  -- Iteration " << iteration + 1 << std::endl;
        vcg::tri::UpdateSelection<TMesh>::FaceClear(mesh);
        if (faces.size() > 0) {
            // Select self-intersect faces
            for (fc = faces.begin(); fc != faces.end(); fc++) {
                (*fc)->SetS();
            }
            // Dilate the faces and vertices
            for (uint16_t dilate_step = 0; dilate_step < iteration + 1; dilate_step++) {
                vcg::tri::UpdateSelection<TMesh>::VertexFromFaceLoose(mesh);
                vcg::tri::UpdateSelection<TMesh>::FaceFromVertexLoose(mesh);
            }
            // Select vertices from current faces and remove all selected
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
            if (verbose) std::cout << "  -- self-intersect faces: " << faces.size() << std::endl;
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            if (verbose) std::cout << "  -- Remove faces [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            if (verbose) std::cout << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            if (verbose) std::cout << "  -- Remove non-manifold edges [OK]" << std::endl;
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                std::cout << "[Warning]: Fixed Self-intersecting failed: Mesh has some not 2-manifold edges" << std::endl;
                return false;
            }
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            if (verbose) std::cout << "  -- Close holes [OK]" << std::endl;
            vcg::tri::UpdateSelection<TMesh>::FaceClear(mesh);
            vcg::tri::UpdateSelection<TMesh>::VertexClear(mesh);
        }

        vcg::tri::UpdateSelection<TMesh>::FaceClear(mesh);
        vcg::tri::UpdateSelection<TMesh>::VertexClear(mesh);
        vcg::tri::UpdateFlags<TMesh>::FaceBorderFromNone(mesh);
        vcg::tri::UpdateFlags<TMesh>::VertexBorderFromFaceBorder(mesh);
        vcg::tri::UpdateSelection<TMesh>::FaceFromBorderFlag(mesh);
        vcg::tri::UpdateSelection<TMesh>::VertexFromBorderFlag(mesh);
        nf = vcg::tri::UpdateSelection<TMesh>::FaceCount(mesh);
        nv = vcg::tri::UpdateSelection<TMesh>::VertexCount(mesh);
        if (verbose) std::cout << "  -- Count border edge v:" << nv << " f:" << nf << std::endl;
        if (nf > 0) {
            // Dilate the faces and vertices
            //for (uint16_t dilate_step = 0; dilate_step < iteration + 1; dilate_step++) {
                vcg::tri::UpdateSelection<TMesh>::VertexFromFaceLoose(mesh);
                vcg::tri::UpdateSelection<TMesh>::FaceFromVertexLoose(mesh);
            //}
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
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            if (verbose) std::cout << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            if (verbose) std::cout << "  -- Remove non-manifold edges [OK]" << std::endl;
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                std::cout << "[Warning]: Fixed Self-intersecting failed: Mesh has some not 2-manifold edges" << std::endl;
                return false;
            }
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            if (verbose) std::cout << "  -- Close holes [OK]" << std::endl;

            vcg::tri::UpdateSelection<TMesh>::FaceClear(mesh);
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
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    
    if (faces.size() > 0 && nf > 0) return false;
    return true;
}

inline void report_mesh(TMesh& mesh) {
    int connectedComponentsNum = vcg::tri::Clean<TMesh>::CountConnectedComponents(mesh);
    std::cout
        << "[Topology Measurement] " << std::endl
        << "-- Mesh is composed by " << connectedComponentsNum << " connected component(s)" << std::endl;

    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
    vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);
    int vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(mesh, false);

    std::cout
        << "-- border edge: " << edgeBorderNum << std::endl
        << "-- non-manifold edge: " << edgeNonManifoldNum << std::endl
        << "-- non-manifold vertex: " << vertManifNum << std::endl;

    if (edgeNonManifoldNum == 0 && vertManifNum == 0) {
        int holeNum = vcg::tri::Clean<TMesh>::CountHoles(mesh);
        int genus = vcg::tri::Clean<TMesh>::MeshGenus(mesh.vn, edgeNum, mesh.fn, holeNum, connectedComponentsNum);

        std::cout
            << "-- Mesh is two-manifold " << std::endl
            << "-- Mesh has " << holeNum << " holes" << std::endl
            << "-- Genus is " << genus << std::endl;
    }
}

inline bool is_mesh_manifold(TMesh& mesh) {
    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
    vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);
    int vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(mesh, false);
    return (edgeNonManifoldNum == 0 && vertManifNum == 0);
}