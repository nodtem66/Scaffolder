#pragma once
#ifndef MESH_INCLUDED
#define MESH_INCLUDED
#include <openctm.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/halfedge_indexed.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/hole.h>
#include <vcg/math/histogram.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ctm.h>
#include <wrap/io_trimesh/export.h>

#ifndef SLICE_PRECISION
#define SLICE_PRECISION 1e-8
#endif

class TFace;
class TVertex;

struct TUsedTypes : public vcg::UsedTypes< vcg::Use<TVertex>::AsVertexType, vcg::Use<TFace>::AsFaceType > {};

// The main vertex class
// Some attributes are optional (ocf) and must enabled before use. See more at vcg::vecter::vector_ocf
// Each vertex either needs 32 bytes (ILP32) or 36 byte (LP64)
// see more: https://en.cppreference.com/w/cpp/language/types
class TVertex : public vcg::Vertex< TUsedTypes,
    vcg::vertex::InfoOcf,  // 4 byte
    vcg::vertex::Coord3f,  // 12 byte
    vcg::vertex::BitFlags, // 4 byte
    vcg::vertex::Normal3f, // 12 byte
    vcg::vertex::CurvaturefOcf, // 0 byte
    vcg::vertex::QualityfOcf, // 0 byte
    vcg::vertex::VFAdjOcf,  // 0 byte
    vcg::vertex::Mark  // 0 byte
> {};

// Each face needs 32 bytes (ILP32) or 36 byte (LP64)
class TFace : public vcg::Face<TUsedTypes,
    vcg::face::InfoOcf,     // 4 byte
    vcg::face::VertexRef,	// 12 byte, three pointers to vertices
    vcg::face::Normal3f,	// 12 byte normal
    vcg::face::BitFlags,    // 4 byte flags
    vcg::face::FFAdj,
    vcg::face::Mark,   
    vcg::face::QualityfOcf, // 0 byte
    vcg::face::VFAdjOcf     // 0 byte
> {};

/* the mesh is a container of vertices and a container of faces */
class TMesh : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<TVertex>, vcg::face::vector_ocf<TFace> > {};
#endif