// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/Happah.hpp"
#include "happah/geometries/Vertex.hpp"

template<class Vertex>
class Model {
     static_assert(is_vertex<Vertex>::value, "A model can only be parameterized by a vertex.");

public:
     using VERTEX = Vertex;
     using Vertices = std::vector<Vertex>;

     virtual ~Model() {}

     hpuint getNumberOfVertices() const { return m_vertices.size(); }

     auto& getVertex(hpuint index) const { return m_vertices[index]; }

     auto& getVertex(hpuint index) { return m_vertices[index]; }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

     auto getVertices(const happah::Indices& indices) const { return getVertices(indices.begin(), indices.end()); }

     template<class Iterator>
     auto getVertices(Iterator begin, Iterator end) const {
          Vertices vertices;
          vertices.reserve(std::distance(begin, end));
          while(begin != end) {
               vertices.push_back(m_vertices[*begin]);
               ++begin;
          }
          return vertices;
     }

protected:
     Vertices m_vertices;

     Model() {}

     Model(Vertices vertices)
          : m_vertices(std::move(vertices)) {}

};

template<class M, class Space = typename M::SPACE, class Vertex = typename M::VERTEX>
struct is_model : std::integral_constant<bool, std::is_base_of<Model<Vertex>, M>::value && std::is_base_of<typename Vertex::SPACE, Space>::value> {};

template<class Model, class Space = typename Model::SPACE, class Vertex = typename Model::VERTEX>
struct is_absolute_model : std::integral_constant<bool, is_model<Model, Space, Vertex>::value && is_absolute_vertex<Vertex>::value> {};

template<class Model, class Space = typename Model::SPACE, class Vertex = typename Model::VERTEX>
struct is_relative_model : std::integral_constant<bool, is_model<Model, Space, Vertex>::value && is_relative_vertex<Vertex>::value> {};

template<class Model, class Geometry>
struct is_relativizable<Model, Geometry, typename std::enable_if<is_relative_model<Model>::value && is_relativizable<typename Model::VERTEX, Geometry>::value>::type> : std::true_type {};

template<class Geometry, class Model, class Space = typename Geometry::SPACE, class Vertex = typename Model::VERTEX>
struct is_modelable : std::integral_constant<bool, is_geometry<Geometry, Space>::value && is_model<Model, Space, Vertex>::value && Model::DIMENSION <= Geometry::DIMENSION && (is_absolute_model<Model>::value || is_relativizable<Model, Geometry>::value)> {};

template<class Model, class Space = typename Model::SPACE, class Vertex = typename Model::VERTEX> 
struct enable_if_absolute_model : std::enable_if<is_absolute_model<Model, Space, Vertex>::value> {};

template<class Model, class Space = typename Model::SPACE, class Vertex = typename Model::VERTEX> 
struct enable_if_relative_model : std::enable_if<is_relative_model<Model, Space, Vertex>::value> {};

template<class Model, class Space = typename Model::SPACE, class Vertex = typename Model::VERTEX> 
struct enable_if_model : std::enable_if<is_model<Model, Space, Vertex>::value> {};

template<class Geometry, class Model, class Space = typename Geometry::SPACE, class Vertex = typename Model::VERTEX>
struct enable_if_modelable : std::enable_if<is_modelable<Geometry, Model, Space, Vertex>::value> {};

