// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <memory>
#include <vector>

#include "happah/geometries/Geometry.hpp"
#include "happah/geometries/Model.hpp"

template<class Vertex>
class VertexCloud : public Geometry0D<typename Vertex::SPACE>, public Model<Vertex> {
     static_assert(is_vertex<Vertex>::value, "A vertex cloud can only be parameterized by a vertex.");

public:
     using Vertices = typename Model<Vertex>::Vertices;

     VertexCloud(Vertices vertices) 
          : Model<Vertex>(std::move(vertices)) {}

     virtual ~VertexCloud() {}

};

typedef VertexCloud<VertexP2> PointCloud2D;
typedef std::shared_ptr<PointCloud2D> PointCloud2D_ptr;
typedef VertexCloud<VertexP3> PointCloud3D;
typedef std::shared_ptr<PointCloud3D> PointCloud3D_ptr;

template<class V, class Space = typename V::SPACE, class Vertex = typename V::VERTEX>
struct is_vertex_cloud : public std::integral_constant<bool, std::is_base_of<VertexCloud<Vertex>, V>::value && std::is_base_of<typename V::SPACE, Space>::value> {};

template<class VertexCloud, class Space = typename VertexCloud::SPACE, class Vertex = typename VertexCloud::VERTEX>
struct is_absolute_vertex_cloud : public std::integral_constant<bool, is_vertex_cloud<VertexCloud, Space, Vertex>::value && is_absolute_vertex<Vertex>::value> {};

template<class VertexCloud, class Space = typename VertexCloud::SPACE, class Vertex = typename VertexCloud::VERTEX>
struct is_relative_vertex_cloud : public std::integral_constant<bool, is_vertex_cloud<VertexCloud, Space, Vertex>::value && is_relative_vertex<Vertex>::value> {};

template<class VertexCloud, class Space = typename VertexCloud::SPACE, class Vertex = typename VertexCloud::VERTEX>
struct enable_if_vertex_cloud : public std::enable_if<is_vertex_cloud<VertexCloud, Space, Vertex>::value> {};

template<class VertexCloud, class Space = typename VertexCloud::SPACE, class Vertex = typename VertexCloud::VERTEX>
struct enable_if_absolute_vertex_cloud : public std::enable_if<is_absolute_vertex_cloud<VertexCloud, Space, Vertex>::value> {};

template<class VertexCloud, class Space = typename VertexCloud::SPACE, class Vertex = typename VertexCloud::VERTEX>
struct enable_if_relative_vertex_cloud : public std::enable_if<is_relative_vertex_cloud<VertexCloud, Space, Vertex>::value> {};

