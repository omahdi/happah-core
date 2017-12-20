// Copyright 2017
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.08 - Hedwig Amberg    - Introduce RectangularCuboid.

#pragma once

#include "happah/geometry/QuadMesh.hpp"
#include "happah/geometry/TriangleMesh.hpp"
#include "happah/util/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

class RectangularCuboid;

namespace rcd {

struct minimal {};
struct symmetric {};

}//namespace rcd

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
QuadMesh<Vertex> make_quad_mesh(const RectangularCuboid& cuboid, VertexFactory&& build = VertexFactory());

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(rcd::minimal, const RectangularCuboid& cuboid, VertexFactory&& build = VertexFactory());

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(rcd::symmetric, const RectangularCuboid& cuboid, VertexFactory&& build = VertexFactory());

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const RectangularCuboid& cuboid, VertexFactory&& build = VertexFactory());

//DEFINITIONS

class RectangularCuboid {
public:
     RectangularCuboid(hpreal width, hpreal height, hpreal depth)
          : m_depth(depth), m_height(height), m_width(width) {}
     
     auto& getDepth() const { return m_depth; }
     
     auto& getHeight() const { return m_height; }
     
     auto& getWidth() const { return m_width; }

private:
     hpreal m_depth;
     hpreal m_height;
     hpreal m_width;

};//RectangularCuboid

template<class Vertex, class VertexFactory>
QuadMesh<Vertex> make_quad_mesh(const RectangularCuboid& cuboid, VertexFactory&& build){
     auto d = cuboid.getDepth();
     auto h = cuboid.getHeight();
     auto w = cuboid.getWidth();
     auto vertices = std::vector<Vertex>({
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
          3, 2, 1, 0,
          1, 2, 6, 5,
          0, 1, 5, 4,
          3, 0, 4, 7,
          2, 3, 7, 6,
          4, 5, 6, 7
     });

     return make_quad_mesh(std::move(vertices), std::move(indices));
}

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(rcd::minimal, const RectangularCuboid& cuboid, VertexFactory&& build) {
     auto d = cuboid.getDepth();
     auto h = cuboid.getHeight();
     auto w = cuboid.getWidth();
     auto vertices = std::vector<Vertex>({
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
          0, 2, 1,
          0, 3, 2,
          5, 2, 6,
          5, 1, 2,
          4, 1, 5,
          4, 0, 1,
          7, 0, 4,
          7, 3, 0,
          6, 3, 7,
          6, 2, 3,
          7, 5, 6,
          7, 4, 5
     });

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(rcd::symmetric, const RectangularCuboid& cuboid, VertexFactory&& build) {
     auto d = cuboid.getDepth();
     auto h = cuboid.getHeight();
     auto w = cuboid.getWidth();
     auto vertices = std::vector<Vertex>({
          build(Point3D(0, 0, 0)),
          build(Point3D(w, 0, 0)),
          build(Point3D(w, h, 0)),
          build(Point3D(0, h, 0)),
          build(Point3D(0, 0, d)),
          build(Point3D(w, 0, d)),
          build(Point3D(w, h, d)),
          build(Point3D(0, h, d)),
          build(Point3D(w / 2, 0, d / 2)),
          build(Point3D(w, h / 2, d / 2)),
          build(Point3D(w / 2, h, d / 2)),
          build(Point3D(0, h / 2, d / 2)),
          build(Point3D(w / 2, h / 2, 0)),
          build(Point3D(w / 2, h / 2, d))
     });
     auto indices = Indices();
     indices.assign({
          12, 0, 3,
          12, 3, 2,
          12, 2, 1,
          12, 1, 0,
          13, 7, 4,
          13, 4, 5,
          13, 5, 6,
          13, 6, 7,
          9, 5, 1,
          9, 1, 2,
          9, 2, 6,
          9, 6, 5,
          11, 7, 3,
          11, 3, 0,
          11, 0, 4,
          11, 4, 7,
          10, 6, 2,
          10, 2, 3,
          10, 3, 7,
          10, 7, 6,
          8, 4, 0,
          8, 0, 1,
          8, 1, 5,
          8, 5, 4
     });

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const RectangularCuboid& cuboid, VertexFactory&& build) { return make_triangle_mesh<Vertex>(rcd::symmetric{}, cuboid, std::forward<VertexFactory>(build)); }

}//namespace happah

