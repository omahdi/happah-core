// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/geometries/TriangleMesh.hpp"
#include "happah/geometries/TriangleMeshUtils.hpp"

namespace happah {

class RectangularTorusChain {
public:
     RectangularTorusChain(hpuint nHandles, hpreal side, hpreal spacing, hpreal hole);

/**********************************************************************************
 * Ordering of the vertices of a torus in the chain.
 *
 * 1 ------ 3
 * | 5 -- 7 |
 * | |    | |
 * | 4 -- 6 |
 * 0 ------ 2
 **********************************************************************************/
     template<class Vertex = VertexP3, class VertexFactory = VertexFactory<Vertex> >
     TriangleMesh<Vertex> toTriangleMesh(VertexFactory&& factory = VertexFactory()) const {
          //TODO: m_space < epsilon or m_hole < epsilon
          std::vector<Vertex> vertices;
          vertices.reserve(m_nHandles << 4);
          std::vector<hpuint> indices;
          hpuint nTriangles = m_nHandles * 36 - 4;
          indices.reserve(nTriangles * 3);

          hpreal o0 = 0.0;
          hpreal o1 = m_side;
          hpreal h0 = m_side * (1.0 - m_hole) / 2.0;
          hpreal h1 = m_side * (1.0 + m_hole) / 2.0;
          hpreal xo = 0.0;
          hpreal xh = h0;
          hpreal delta0 = m_side * m_hole;
          hpreal delta1 = h0 * 2.0 + m_spacing;
          for(hpuint i = 0; i < m_nHandles; ) {
               // outer vertices
               vertices.push_back(factory(Point3D(o0, xo, 0.0)));
               vertices.push_back(factory(Point3D(o1, xo, 0.0)));
               xo += m_side;
               vertices.push_back(factory(Point3D(o0, xo, 0.0)));
               vertices.push_back(factory(Point3D(o1, xo, 0.0)));
               xo += m_spacing;

               // inner vertices
               vertices.push_back(factory(Point3D(h0, xh, 0.0)));
               vertices.push_back(factory(Point3D(h1, xh, 0.0)));
               xh += delta0;
               vertices.push_back(factory(Point3D(h0, xh, 0.0)));
               vertices.push_back(factory(Point3D(h1, xh, 0.0)));
               xh += delta1;

               hpuint i0 = (i<<3), i1 = i0+1, i2 = i1+1, i3 = i2+1, i4 = i3+1, i5 = i4+1, i6 = i5+1, i7 = i6+1, i8 = i7+1, i9 = i8+1;
               hpuint is[] = { i0, i1, i5, i0, i5, i4, i0, i4, i6, i0, i6, i2, i3, i2, i6, i3, i6, i7, i3, i7, i5, i3, i5, i1, i3, i9, i2, i8, i2, i9 };
               if(++i < m_nHandles) indices.insert(indices.end(), is, is+30);
               else indices.insert(indices.end(), is, is+24);
          }
          
          TriangleMeshUtils::extrude(vertices, indices.begin(), indices.size(), std::back_inserter(indices), m_side);

          // sides
          hpuint deltav = m_nHandles << 3;
          {
               hpuint i0 = 0, i1 = i0+1, j0 = i0+deltav, j1 = i1+deltav;
               hpuint is[] = { i0, j0, i1, i1, j0, j1 };
               indices.insert(indices.end(), is, is+6);
          }
          for(hpuint i = 0; i < m_nHandles; ) {
               hpuint i0 = (i<<3), i1 = i0+1, i2 = i1+1, i3 = i2+1, i4 = i3+1, i5 = i4+1, i6 = i5+1, i7 = i6+1, i8 = i7+1, i9 = i8+1;
               hpuint j0 = i0+deltav, j1 = i1+deltav, j2 = i2+deltav, j3 = i3+deltav, j4 = i4+deltav, j5 = i5+deltav, j6 = i6+deltav, j7 = i7+deltav, j8 = i8+deltav, j9 = i9+deltav;
               hpuint is[] = { i0, i2, j0, i2, j2, j0, i1, j3, i3, i1, j1, j3, i4, j6, i6, i4, j4, j6, i5, j4, i4, i5, j5, j4, i7, j5, i5, i7, j7, j5, i6, j7, i7, i6, j6, j7, i3, j9, i9, i3, j3, j9, i2, i8, j2, i8, j8, j2 };
               if(++i < m_nHandles) indices.insert(indices.end(), is, is+48);
               else indices.insert(indices.end(), is, is+36);
          }
          {
               hpuint i0 = (m_nHandles << 3) - 6, i1 = i0+1, j0 = i0+deltav, j1 = i1+deltav;
               hpuint is[] = { i1, j1, i0, i0, j1, j0 };
               indices.insert(indices.end(), is, is+6);
          }

          return TriangleMesh<Vertex>(std::move(vertices), std::move(indices));
     } 

private:
     hpreal m_hole;
     hpuint m_nHandles;
     hpreal m_side;
     hpreal m_spacing;

};//RectangularTorusChain

}//namespace happah

