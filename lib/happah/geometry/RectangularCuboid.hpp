// Copyright 2017
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.08 - Hedwig Amberg    - Introduce RectangularCuboid.

#pragma once

#include "QuadMesh.hpp"
#include "happah/util/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

class RectangularCuboid;

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
QuadMesh<Vertex> make_quad_mesh(const RectangularCuboid& recc, VertexFactory&& build = VertexFactory());

//DEFINITIONS

class RectangularCuboid {
public:
     RectangularCuboid(hpreal width, hpreal height, hpreal depth)
          : m_width(width), m_height(height), m_depth(depth) {}

     auto& getWidth() const { return m_width; }
     
     auto& getHeight() const { return m_height; }
     
     auto& getDepth() const { return m_depth; }

private:
     hpreal m_width;
     hpreal m_height;
     hpreal m_depth;

};//RectangularCuboid

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
QuadMesh<Vertex> make_quad_mesh(const RectangularCuboid& recc, VertexFactory&& build){
     auto w = recc.getWidth();
     auto h = recc.getHeight();
     auto d = recc.getDepth();
     auto vertices = std::vector<Vertex>();
     vertices.assign({
          build(Point3D(0, 0, 0)),
          build(Point3D(w, 0, 0)),
          build(Point3D(w, h, 0)),
          build(Point3D(0, h, 0)),
          build(Point3D(0, 0, d)),
          build(Point3D(w, 0, d)),
          build(Point3D(w, h, d)),
          build(Point3D(0, h, d))
     });
     auto indices = Indices();
     indices.assign({
          0, 1, 2, 3,
          1, 5, 6, 2,
          0, 4, 5, 1,
          3, 7, 4, 0,
          2, 6, 7, 3,
          5, 4, 7, 6
     });
     return make_quad_mesh(std::move(vertices), std::move(indices));
}

}//namespace happah
