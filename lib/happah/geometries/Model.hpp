// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/Happah.hpp"
#include "happah/geometries/Vertex.hpp"

namespace happah {

template<class Vertex>
class Model {
public:
     virtual ~Model() {}

     hpuint getNumberOfVertices() const { return m_vertices.size(); }

     auto& getVertex(hpuint index) const { return m_vertices[index]; }

     auto& getVertex(hpuint index) { return m_vertices[index]; }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

     auto getVertices(const happah::Indices& indices) const { return getVertices(indices.begin(), indices.end()); }

     template<class Iterator>
     auto getVertices(Iterator begin, Iterator end) const {
          std::vector<Vertex> vertices;
          vertices.reserve(std::distance(begin, end));
          while(begin != end) {
               vertices.push_back(m_vertices[*begin]);
               ++begin;
          }
          return vertices;
     }

protected:
     std::vector<Vertex> m_vertices;

     Model() {}

     Model(std::vector<Vertex> vertices)
          : m_vertices(std::move(vertices)) {}

};

}//namespace happah

