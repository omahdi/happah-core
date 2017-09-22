// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.08 - Hedwig Amberg    - Implemented nutchain triangle mesh.

#pragma once

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

class NutChain {
     
public:
     NutChain(hpuint nNuts, hpreal outerLength, hpreal innerLength, hpreal thickness, hpreal padding)
          : m_nNuts(nNuts), m_innerLength(innerLength), m_outerLength(outerLength), m_padding(padding), m_thickness(thickness) {}

     hpreal getInnerLength() const { return m_innerLength; }

     hpuint getNumberOfNuts() const { return m_nNuts; }

     hpreal getOuterLength() const { return m_outerLength; }

     hpreal getPadding() const { return m_padding; }

     hpreal getThickness() const { return m_thickness; }

private:
     hpreal m_innerLength;
     hpuint m_nNuts;
     hpreal m_outerLength;
     hpreal m_padding;
     hpreal m_thickness;

};//NutChain

//Assume number of nuts is greater than zero.
template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const NutChain& chain, VertexFactory&& build) {
     auto indices = Indices();
     auto vertices = std::vector<Vertex>();
     auto nNuts = chain.getNumberOfNuts();
     auto outerLength = chain.getOuterLength();
     auto innerLength = chain.getInnerLength();
     auto padding = chain.getPadding();
     auto thickness = chain.getThickness();
     auto width = (outerLength - innerLength) / 2.0;
     auto a = 30;
     auto b = 3 * 56;

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
          hpuint(1), hpuint(13 + a * (nNuts - 1)), hpuint( 3 + a * (nNuts - 1)),
          hpuint(1), hpuint( 3 + a * (nNuts - 1)), hpuint(25 + a * (nNuts - 1)),
          hpuint(1), hpuint(25 + a * (nNuts - 1)), hpuint(15 + a * (nNuts - 1)),
          hpuint(1), hpuint(15 + a * (nNuts - 1)), hpuint(13 + a * (nNuts - 1)),

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

     for(auto n = hpuint(1); n < nNuts; ++n) {
          vertices.insert(std::end(vertices), std::end(vertices) - a, std::end(vertices));
          indices.insert(std::end(indices), std::end(indices) - b, std::end(indices));
          for(auto i = std::end(vertices) - a, end = std::end(vertices); i != end; ++i) (*i).position.x += padding + outerLength;
          for(auto i = std::end(indices) - b, end = std::end(indices); i != end; ++i) *i += a;
     }

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const NutChain& chain) { return chain.getNumberOfNuts(); }

}//namespace happah

