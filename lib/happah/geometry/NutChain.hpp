// Copyright 2017
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.08 - Hedwig Amberg    - Implemented nutchain triangle mesh.

#pragma once

#include <boost/range.hpp>

#include "happah/Happah.hpp"
#include "happah/geometry/TriangleMesh.hpp"
#include "happah/util/VertexFactory.hpp"

namespace happah {

//DEFINITIONS

class NutChain;

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const NutChain& chain, VertexFactory&& build = VertexFactory());

inline hpuint size(const NutChain& chain);

//DECLARATIONS

//TODO: make padding boost::optional<hpreal> and update make_triangle_mesh
//TODO: make_nut_chain with parameter sanity checking
class NutChain {
public:
     NutChain(hpuint nNuts, hpreal innerRadius, hpreal outerRadius, hpreal thickness, hpreal padding)
          : m_innerRadius(innerRadius), m_nNuts(nNuts), m_outerRadius(outerRadius), m_padding(padding), m_thickness(thickness) {}

     auto& getInnerRadius() const { return m_innerRadius; }

     auto& getNumberOfNuts() const { return m_nNuts; }

     auto& getOuterRadius() const { return m_outerRadius; }

     auto& getPadding() const { return m_padding; }

     auto& getThickness() const { return m_thickness; }

private:
     hpreal m_innerRadius;
     hpuint m_nNuts;
     hpreal m_outerRadius;
     hpreal m_padding;
     hpreal m_thickness;

};//NutChain

//Assume number of nuts is greater than zero.
template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const NutChain& chain, VertexFactory&& build) {
/*
 * T = top; M = middle; B = bottom
 * B14;T12 -------------------- M30 -------------------- B15;T13 ---- M33 ---- B48;T46 ---   ---
 *    |                       B16;T11                       |                     |             |
 *    |             B17;T09 --- M29 --- B18;T10             |                     |             |
 *   M00   B19;T07    M26                 M31    B20;T08    |       B34;T32       |      ...   M01
 *    |             B21;T05 --- M28 --- B22;T06             |                     |             |
 *    |                       B23;T04                       |                     |             |
 * B24;T02 -------------------- M27 -------------------- B25;T03 ---- M35 ---- B58;T36 ---   ---
 */

     auto indices = Indices();
     auto vertices = std::vector<Vertex>();
     auto nNuts = chain.getNumberOfNuts();
     auto outerLength = 0.707106781 * chain.getOuterRadius();
     auto innerLength = 0.707106781 * chain.getInnerRadius();
     auto padding = chain.getPadding();
     auto thickness = chain.getThickness();
     auto width = (outerLength - innerLength) / 2.0;

     vertices.reserve(2 + 30 * nNuts + 4 * (nNuts - 1));
     indices.reserve(3 * (8 + 56 * nNuts + 16 * (nNuts - 1)));

     vertices.assign({
          //midpoints at beginning and end of chain
          build(Point3D(0, outerLength / 2.0, thickness / 2.0)),
          build(Point3D(outerLength * nNuts + padding * (nNuts - 1), outerLength / 2.0, thickness / 2.0)),

          //top points
          build(Point3D(0,                         0,                         thickness)),
          build(Point3D(outerLength,               0,                         thickness)),
          build(Point3D(outerLength / 2.0,         width / 2.0,               thickness)),
          build(Point3D(width,                     width,                     thickness)),
          build(Point3D(width + innerLength,       width,                     thickness)),
          build(Point3D(width / 2.0,               outerLength / 2.0,         thickness)),
          build(Point3D(outerLength - width / 2.0, outerLength / 2.0,         thickness)),
          build(Point3D(width,                     width + innerLength,       thickness)),
          build(Point3D(width + innerLength,       width + innerLength,       thickness)),
          build(Point3D(outerLength / 2.0,         outerLength - width / 2.0, thickness)),
          build(Point3D(0,                         outerLength,               thickness)),
          build(Point3D(outerLength,               outerLength,               thickness)),

          //bottom points
          build(Point3D(0,                         outerLength,               0)),
          build(Point3D(outerLength,               outerLength,               0)),
          build(Point3D(outerLength / 2.0,         outerLength - width / 2.0, 0)),
          build(Point3D(width,                     width + innerLength,       0)),
          build(Point3D(width + innerLength,       width + innerLength,       0)),
          build(Point3D(width / 2.0,               outerLength / 2.0,         0)),
          build(Point3D(outerLength - width / 2.0, outerLength / 2.0,         0)),
          build(Point3D(width,                     width,                     0)),
          build(Point3D(width + innerLength,       width,                     0)),
          build(Point3D(outerLength / 2.0,         width / 2.0,               0)),
          build(Point3D(0,                         0,                         0)),
          build(Point3D(outerLength,               0,                         0)),

          //middle midpoints
          build(Point3D(width,               outerLength / 2.0,   thickness / 2.0)),
          build(Point3D(outerLength / 2.0,   0,                   thickness / 2.0)),
          build(Point3D(outerLength / 2.0,   width,               thickness / 2.0)),
          build(Point3D(outerLength / 2.0,   width + innerLength, thickness / 2.0)),
          build(Point3D(outerLength / 2.0,   outerLength,         thickness / 2.0)),
          build(Point3D(width + innerLength, outerLength / 2.0,   thickness / 2.0))
     });

     indices.assign({
          //triangles at beginning and end of chain
          hpuint(0), hpuint(2),  hpuint(12),
          hpuint(0), hpuint(12), hpuint(14),
          hpuint(0), hpuint(14), hpuint(24),
          hpuint(0), hpuint(24), hpuint(2),
          hpuint(1), hpuint(13 + 34 * (nNuts - 1)), hpuint( 3 + 34 * (nNuts - 1)),
          hpuint(1), hpuint( 3 + 34 * (nNuts - 1)), hpuint(25 + 34 * (nNuts - 1)),
          hpuint(1), hpuint(25 + 34 * (nNuts - 1)), hpuint(15 + 34 * (nNuts - 1)),
          hpuint(1), hpuint(15 + 34 * (nNuts - 1)), hpuint(13 + 34 * (nNuts - 1)),

          //top triangles
          hpuint(2),  hpuint(3),  hpuint(4),
          hpuint(2),  hpuint(4),  hpuint(5),
          hpuint(2),  hpuint(5),  hpuint(7),
          hpuint(2),  hpuint(7),  hpuint(12),
          hpuint(3),  hpuint(6),  hpuint(4),
          hpuint(3),  hpuint(8),  hpuint(6),
          hpuint(3),  hpuint(13), hpuint(8),
          hpuint(13), hpuint(10), hpuint(8),
          hpuint(13), hpuint(11), hpuint(10),
          hpuint(13), hpuint(12), hpuint(11),
          hpuint(12), hpuint(9),  hpuint(11),
          hpuint(12), hpuint(7),  hpuint(9),
          hpuint(4),  hpuint(6),  hpuint(5),
          hpuint(8),  hpuint(10), hpuint(6),
          hpuint(11), hpuint(9),  hpuint(10),
          hpuint(7),  hpuint(5),  hpuint(9),

          //bottom triangles
          hpuint(14), hpuint(15), hpuint(16),
          hpuint(14), hpuint(16), hpuint(17),
          hpuint(14), hpuint(17), hpuint(19),
          hpuint(14), hpuint(19), hpuint(24),
          hpuint(15), hpuint(18), hpuint(16),
          hpuint(15), hpuint(20), hpuint(18),
          hpuint(15), hpuint(25), hpuint(20),
          hpuint(25), hpuint(22), hpuint(20),
          hpuint(25), hpuint(23), hpuint(22),
          hpuint(25), hpuint(24), hpuint(23),
          hpuint(24), hpuint(21), hpuint(23),
          hpuint(24), hpuint(19), hpuint(21),
          hpuint(16), hpuint(18), hpuint(17),
          hpuint(20), hpuint(22), hpuint(18),
          hpuint(23), hpuint(21), hpuint(22),
          hpuint(19), hpuint(17), hpuint(21),

          //middle triangles
          hpuint(26), hpuint(9),  hpuint(5),
          hpuint(26), hpuint(5),  hpuint(21),
          hpuint(26), hpuint(21), hpuint(17),
          hpuint(26), hpuint(17), hpuint(9),
          hpuint(27), hpuint(3),  hpuint(2),
          hpuint(27), hpuint(2),  hpuint(24),
          hpuint(27), hpuint(24), hpuint(25),
          hpuint(27), hpuint(25), hpuint(3),
          hpuint(28), hpuint(5),  hpuint(6),
          hpuint(28), hpuint(6),  hpuint(22),
          hpuint(28), hpuint(22), hpuint(21),
          hpuint(28), hpuint(21), hpuint(5),
          hpuint(29), hpuint(10), hpuint(9),
          hpuint(29), hpuint(9),  hpuint(17),
          hpuint(29), hpuint(17), hpuint(18),
          hpuint(29), hpuint(18), hpuint(10),
          hpuint(30), hpuint(12), hpuint(13),
          hpuint(30), hpuint(13), hpuint(15),
          hpuint(30), hpuint(15), hpuint(14),
          hpuint(30), hpuint(14), hpuint(12),
          hpuint(31), hpuint(6),  hpuint(10),
          hpuint(31), hpuint(10), hpuint(18),
          hpuint(31), hpuint(18), hpuint(22),
          hpuint(31), hpuint(22), hpuint(6)
     });

     if(nNuts > 1) {
          auto temp = {
               //padding triangles
               hpuint(32), hpuint(13), hpuint(3),
               hpuint(32), hpuint(3),  hpuint(36),
               hpuint(32), hpuint(36), hpuint(46),
               hpuint(32), hpuint(46), hpuint(13),
               hpuint(33), hpuint(13), hpuint(46),
               hpuint(33), hpuint(46), hpuint(48),
               hpuint(33), hpuint(48), hpuint(15),
               hpuint(33), hpuint(15), hpuint(13),
               hpuint(34), hpuint(15), hpuint(48),
               hpuint(34), hpuint(48), hpuint(58),
               hpuint(34), hpuint(58), hpuint(25),
               hpuint(34), hpuint(25), hpuint(15),
               hpuint(35), hpuint(36), hpuint(3),
               hpuint(35), hpuint(3),  hpuint(25),
               hpuint(35), hpuint(25), hpuint(58),
               hpuint(35), hpuint(58), hpuint(36)
          };

          //padding points
          vertices.push_back(build(Point3D(outerLength + padding / 2.0, outerLength / 2.0, thickness)));
          vertices.push_back(build(Point3D(outerLength + padding / 2.0, outerLength,       thickness / 2.0)));
          vertices.push_back(build(Point3D(outerLength + padding / 2.0, outerLength / 2.0, 0)));
          vertices.push_back(build(Point3D(outerLength + padding / 2.0, 0,                 thickness / 2.0)));
          indices.insert(std::end(indices), std::begin(temp), std::end(temp));

          auto v0 = std::begin(vertices) + 2;
          auto v1 = std::end(vertices) - 4;
          auto i0 = std::begin(indices) + 3 * 8;
          auto i1 = std::end(indices) - 3 * 16;

          vertices.insert(std::end(vertices), v0, v1);
          indices.insert(std::end(indices), i0, i1);
          for(auto& vertex : boost::make_iterator_range(v1 + 4, std::end(vertices))) vertex.position.x += padding + outerLength;
          for(auto& i : boost::make_iterator_range(i1 + 3 * 16, std::end(indices))) i += 34;

          for(auto n = hpuint(2); n < nNuts; ++n) {
               v0 = v1;
               v1 = std::end(vertices);
               i0 = i1;
               i1 = std::end(indices);
               vertices.insert(v1, v0, v1);
               indices.insert(i1, i0, i1);
               for(auto& vertex : boost::make_iterator_range(v1, std::end(vertices))) vertex.position.x += padding + outerLength;
               for(auto& i : boost::make_iterator_range(i1, std::end(indices))) i += 34;
          }
     }

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const NutChain& chain) { return chain.getNumberOfNuts(); }

}//namespace happah

