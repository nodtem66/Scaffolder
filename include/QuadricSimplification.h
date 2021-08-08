/****************************************************************************
 * MeshLab                                                           o o     *
 * A versatile mesh processing toolbox                             o     o   *
 *                                                                _   O  _   *
 * Copyright(C) 2005                                                \/)\/    *
 * Visual Computing Lab                                            /\/|      *
 * ISTI - Italian National Research Council                           |      *
 *                                                                    \      *
 * All rights reserved.																											 *
 * This program is free software; you can redistribute it and/or modify      *
 * it under the terms of the GNU General Public License as published by      *
 * the Free Software Foundation; either version 2 of the License, or         *
 * (at your option) any later version.                                       *
 *                                                                           *
 * This program is distributed in the hope that it will be useful,           *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 * GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
 * for more details.                                                         *
 *                                                                           *
 ****************************************************************************/
#ifndef QUADSIM_INCLUDED
#define QUADSIM_INCLUDED

#include <vcg/container/simple_temporary_data.h>
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#include "Mesh.h"

namespace vcg {
namespace tri {

    typedef	SimpleTempData<TMesh::VertContainer, vcg::math::Quadric<double> > QuadricTemp;
    
    class QHelper
    {
    public:
        QHelper() {}
        static void Init() {}
        static math::Quadric<double>& Qd(TVertex& v) { return TD()[v]; }
        static math::Quadric<double>& Qd(TVertex* v) { return TD()[*v]; }
        static TVertex::ScalarType W(TVertex*) { return 1.0; }
        static TVertex::ScalarType W(TVertex&) { return 1.0; }
        static void Merge(TVertex&, TVertex const&) {}
        static QuadricTemp*& TDp() { static QuadricTemp* td; return td; }
        static QuadricTemp& TD() { return *TDp(); }
    };

    typedef BasicVertexPair<TVertex> VertexPair;

    class MyTriEdgeCollapse : public vcg::tri::TriEdgeCollapseQuadric< TMesh, VertexPair, MyTriEdgeCollapse, QHelper > {
    public:
        typedef  vcg::tri::TriEdgeCollapseQuadric< TMesh, VertexPair, MyTriEdgeCollapse, QHelper> TECQ;
        inline MyTriEdgeCollapse(const VertexPair& p, int i, BaseParameterClass* pp) :TECQ(p, i, pp) {}
    };
    
    void QuadricSimplification(TMesh& m, int  TargetFaceNum, vcg::tri::TriEdgeCollapseQuadricParameter& pp, vcg::CallBackPos* cb) {
        vcg::math::Quadric<double> QZero;
        QZero.SetZero();
        vcg::tri::QuadricTemp TD(m.vert, QZero);
        vcg::tri::QHelper::TDp() = &TD;

        if (pp.NormalCheck) pp.NormalThrRad = M_PI / 4.0;

        vcg::LocalOptimization<TMesh> DeciSession(m, &pp);
        cb(1, "Initializing simplification");
        DeciSession.Init<vcg::tri::MyTriEdgeCollapse >();

        DeciSession.SetTargetSimplices(TargetFaceNum);
        DeciSession.SetTimeBudget(0.1f); // this allows updating the progress bar 10 time for sec...
        //  if(TargetError< numeric_limits<double>::max() ) DeciSession.SetTargetMetric(TargetError);
        //int startFn=m.fn;
        int faceToDel = m.fn - TargetFaceNum;
        cb(2, "Are you ready... 3 2 1");
        while (DeciSession.DoOptimization() && m.fn > TargetFaceNum)
        {
            cb(100 - (98 * (m.fn - TargetFaceNum) / (faceToDel)), "Simplifying...");
        };

        DeciSession.Finalize<vcg::tri::MyTriEdgeCollapse >();

        vcg::tri::QHelper::TDp() = nullptr;
    }

    void mesh_quad_simplification(TMesh& mesh, double percent, vcg::CallBackPos* cb) {
        if (percent >= mesh.FN()) return;
        mesh.vert.EnableVFAdjacency();
        mesh.face.EnableQuality();
        mesh.face.EnableVFAdjacency();
        //vcg::tri::UpdateNormal<TMesh>::PerFace(mesh);
        vcg::tri::UpdateTopology<TMesh>::VertexFace(mesh);
        vcg::tri::UpdateFlags<TMesh>::FaceBorderFromVF(mesh);
        size_t TargetFaceNum = (percent > 0 && percent < 1) ? mesh.FN() * percent : percent;
        vcg::tri::TriEdgeCollapseQuadricParameter pp;
        pp.QualityThr = 0.3;
        pp.FastPreserveBoundary = true;
        pp.PreserveBoundary = false;
        pp.BoundaryQuadricWeight = 1.25;
        pp.PreserveTopology = true;
        pp.NormalCheck = false;
        pp.OptimalPlacement = true;
        pp.QualityWeight = false;
        pp.QualityQuadric = true;
        pp.QualityQuadricWeight = 0.001f;

        vcg::tri::QuadricSimplification(mesh, TargetFaceNum, pp, cb);

        if (true)
        {
            vcg::tri::Clean<TMesh>::RemoveFaceOutOfRangeArea(mesh, 0);
            vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            vcg::tri::Allocator<TMesh>::CompactVertexVector(mesh);
            vcg::tri::Allocator<TMesh>::CompactFaceVector(mesh);
        }
        mesh.vert.DisableVFAdjacency();
        mesh.face.DisableQuality();
        mesh.face.DisableVFAdjacency();
    }
} // end namespace tri
} // end namespace vcg
#endif