#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/halfedge_indexed.h>
#include <vcg/complex/algorithms/clean.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

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