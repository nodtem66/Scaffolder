#pragma once
#include "Mesh.h"
#include "utils.h"

inline void clean_mesh(TMesh& mesh, double minimum_diameter, uint16_t smooth_step, std::ostream& log) {
    log << "[libVCG Cleaning] ";
    vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
    vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::UpdateBounding<TMesh>::Box(mesh);
    vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
    vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::UpdateBounding<TMesh>::Box(mesh);
    log << "OK" << std::endl;
    if (smooth_step > 0) {
        log << "[Laplacian smoothing] ";
        vcg::tri::Smooth<TMesh>::VertexCoordLaplacian(mesh, smooth_step, false, true);
        log << "OK" << std::endl;
    }
    vcg::tri::Clean<TMesh>::MergeCloseVertex(mesh, SLICE_PRECISION * 1000);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
}

inline bool fix_non_manifold_edges(TMesh& mesh, double minimum_diameter, uint16_t max_iteration = 5, std::ostream& log = util::null_stream) {
    size_t nf = 1;
    int maxSize = mesh.bbox.SquaredDiag();

    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
    vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);

    for (uint16_t iteration = 0; (edgeBorderNum > 0 || edgeNonManifoldNum > 0) && iteration < max_iteration; iteration++) {
        log << "[Fix non-manifold edges]" << std::endl << "  -- Iteration " << iteration + 1 << std::endl;
        log << "  -- non-manifold edges: " << edgeNonManifoldNum << std::endl;
        if (edgeNonManifoldNum > 0) {
            // fix floating faces
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            log << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            // fix non-manifold edges
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            log << "  -- Remove non-manifold edges [OK]" << std::endl;
            // fix holes
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                log << "[Warning]: closing holes failed: Mesh has some not 2-manifold edges" << std::endl;
            }
            else {
                vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
                log << "  -- Close holes [OK]" << std::endl;
            }
        }
        log << "  -- border edges: " << edgeBorderNum << std::endl;
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

            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            log << "  -- Remove Border faces [OK]" << std::endl;
            // fix floating faces
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            log << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            // fix non-manifold edges
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            log << "  -- Remove non-manifold edges [OK]" << std::endl;
            // fix holes
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                log << "[Warning]: closing holes failed: Mesh has some not 2-manifold edges" << std::endl;
                break;
            }
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            log << "  -- Close holes [OK]" << std::endl;
        }

        vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);
    }

    if (edgeBorderNum > 0 || edgeNonManifoldNum > 0) return false;
    return true;
}

inline bool fix_non_manifold_vertices(TMesh& mesh, double minimum_diameter, uint16_t max_iteration = 5, std::ostream& log = util::null_stream) {
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    
    int maxSize = mesh.bbox.SquaredDiag();
    int vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(mesh, true);

    for (uint16_t iteration = 0; (vertManifNum > 0) && iteration < max_iteration; iteration++) {
        log << "[Fix non-manifold vertices]" << std::endl << "  -- Iteration " << iteration + 1 << std::endl;
        log << "  -- non-manifold vertices: " << vertManifNum << std::endl;
        vcg::tri::Clean<TMesh>::SplitNonManifoldVertex(mesh, 0.5*(iteration+1));
        vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
        log << "  -- Remove Border faces [OK]" << std::endl;
        // fix floating faces
        vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
        vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
        log << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
        // fix non-manifold edges
        vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
        vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
        log << "  -- Remove non-manifold edges [OK]" << std::endl;
        // fix holes
        if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
            log << "[Warning]: closing holes failed: Mesh has some not 2-manifold edges" << std::endl;
            break;
        }
        vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
        log << "  -- Close holes [OK]" << std::endl;
        vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(mesh, true);
    }
    if (vertManifNum > 0) return false;
    return true;
}

inline bool fix_self_intersect_mesh(TMesh& mesh, double minimum_diameter, uint16_t max_iteration = 10, std::ostream& log = util::null_stream) {

    mesh.face.EnableVFAdjacency();

    std::vector<TMesh::FaceType*> faces;
    std::vector<TMesh::FaceType*>::iterator fc;
    vcg::tri::Clean<TMesh>::SelfIntersections(mesh, faces);

    uint16_t iteration = 0;
    int maxSize = mesh.bbox.SquaredDiag();
    size_t nf = 0, nv = 0;
    while (iteration < max_iteration && (faces.size() > 0 || nf > 0)) {

        log << "[Fix Self-intersect face]" << std::endl << "  -- Iteration " << iteration + 1 << std::endl;
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
            log << "  -- self-intersect faces: " << faces.size() << std::endl;
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            log << "  -- Remove faces [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            log << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            log << "  -- Remove non-manifold edges [OK]" << std::endl;
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                log << "[Warning]: closing holes failed: Mesh has some not 2-manifold edges" << std::endl;
                return false;
            }
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            log << "  -- Close holes [OK]" << std::endl;
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
        log << "  -- Count border edge v:" << nv << " f:" << nf << std::endl;
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
            log << "  -- Remove Border faces [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            log << "  -- Remove small components " << minimum_diameter * mesh.bbox.Diag() << " [OK]" << std::endl;
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            log << "  -- Remove non-manifold edges [OK]" << std::endl;
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                log << "[Warning]: closing holes failed: Mesh has some not 2-manifold edges" << std::endl;
                return false;
            }
            vcg::tri::Hole<TMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<TMesh>>(mesh, maxSize, false);
            log << "  -- Close holes [OK]" << std::endl;

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

    mesh.face.DisableVFAdjacency();

    if (faces.size() > 0 && nf > 0) return false;
    return true;
}

inline void report_mesh(std::ostream& log, TMesh& mesh) {
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    int connectedComponentsNum = vcg::tri::Clean<TMesh>::CountConnectedComponents(mesh);
    log << "[Topology Measurement] " << std::endl
        << "-- Mesh is composed by " << connectedComponentsNum << " connected component(s)" << std::endl;

    int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
    vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);

    int vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(mesh, false);

    std::vector<TMesh::FaceType*> faces;
    vcg::tri::Clean<TMesh>::SelfIntersections(mesh, faces);

    log << "-- border edge: " << edgeBorderNum << std::endl
        << "-- non-manifold edge: " << edgeNonManifoldNum << std::endl
        << "-- non-manifold vertex: " << vertManifNum << std::endl
        << "-- self-intersect faces: " << faces.size() << std::endl;

    if (edgeNonManifoldNum == 0 && vertManifNum == 0) {
        int holeNum = vcg::tri::Clean<TMesh>::CountHoles(mesh);
        int genus = vcg::tri::Clean<TMesh>::MeshGenus(mesh.vn, edgeNum, mesh.fn, holeNum, connectedComponentsNum);

        log << "-- Mesh is two-manifold " << std::endl
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