// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/geometry/TriangleMesh.hpp"

namespace happah {

template<class Vertex>
class TriangleArray {
public:
     TriangleArray(std::vector<Vertex> vertices)
          : m_vertices(std::move(vertices)) {}

     hpuint getNumberOfTriangles() const { return m_vertices.size() / 3; }

     auto& getVertices() const { return m_vertices; }

private:
     std::vector<Vertex> m_vertices;

};//TriangleArray

template<class Vertex>
TriangleArray<Vertex> make_triangle_array(std::vector<Vertex> vertices) { return { std::move(vertices) }; }

template<class Vertex>
TriangleArray<Vertex> make_triangle_array(const TriangleMesh<Vertex>& mesh) {
     auto vertices = std::vector<Vertex>();

     vertices.reserve(3 * size(mesh));
     for(auto& vertex : deindex(mesh.getVertices(), mesh.getIndices())) vertices.push_back(vertex);

     return make_triangle_array(std::move(vertices));
}

template<class Vertex>
hpuint size(const TriangleArray<Vertex>& array) { return array.getNumberOfTriangles(); }

}//namespace happah

