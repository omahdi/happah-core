// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/range/iterator_range.hpp>
#include <memory>
#include <vector>

namespace happah {

template<class Vertex>
class VertexCloud {
public:
     VertexCloud(std::vector<Vertex> vertices) 
          : m_vertices(std::move(vertices)) {}

     hpuint getNumberOfVertices() const { return m_vertices.size(); }

     auto& getVertex(hpindex v) const { return m_vertices[v]; }

     auto& getVertex(hpindex v) { return m_vertices[v]; }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

     auto getVertices(const Indices& indices) const { return getVertices(std::begin(indices), std::end(indices)); }

     template<class Iterator>
     auto getVertices(Iterator begin, Iterator end) const {
          auto vertices = std::vector<Vertex>();
          vertices.reserve(std::distance(begin, end));
          for(auto v : boost::make_iterator_range(begin, end)) vertices.push_back(m_vertices[v]);
          return vertices;
     }

private:
     std::vector<Vertex> m_vertices;

};//VertexCloud

}//namespace happah

