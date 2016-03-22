// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <type_traits>
#include <utility>

#include "happah/geometries/Mesh.h"
#include "happah/geometries/Geometry.h"

template<class Vertex>
class QuadMesh : public Geometry2D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;

public:
     using Indices = typename Mesh<Vertex>::Indices;
     using Vertices = typename Mesh<Vertex>::Vertices;

     QuadMesh(Vertices vertices, Indices indices)
          : Geometry2D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)) {}

     virtual ~QuadMesh() {}

};
using QuadMesh2D = QuadMesh<VertexP2>;
using QuadMesh3D = QuadMesh<VertexP3N>;

template<class Q, class Space = typename Q::SPACE, class Vertex = typename Q::VERTEX>
struct is_quad_mesh : std::integral_constant<bool, std::is_base_of<QuadMesh<Vertex>, Q>::value && std::is_base_of<typename Q::SPACE, Space>::value> {};

