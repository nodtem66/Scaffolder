#ifndef SCAFFOLD_TMESH_H
#define SCAFFOLD_TMESH_H
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/halfedge_indexed.h>
#include <vcg/complex/algorithms/clean.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <Eigen/Core>

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
#endif