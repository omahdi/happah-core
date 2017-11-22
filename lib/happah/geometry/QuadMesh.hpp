// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.11 - Hedwig Amberg    - Introduce QuadMesh.

#pragma once

#include "happah/geometry/TriangleMesh.hpp"

namespace happah {

//DECLARATIONS

template<class Vertex>
class QuadMesh;

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(std::vector<Vertex> vertices, Indices indices);

//Convert a string representation in HPH format.
template<class Vertex = VertexP3>
QuadMesh<Vertex> make_quad_mesh(const std::string& mesh);

//Import data stored in the given file in HPH format.
template<class Vertex = VertexP3>
QuadMesh<Vertex> make_quad_mesh(const std::experimental::filesystem::path& mesh);

template<class Vertex = VertexP3>
TriangleMesh<Vertex> make_triangle_mesh(const QuadMesh<Vertex>&  mesh);

template<class Vertex>
hpuint size(const QuadMesh<Vertex>& mesh);

//DEFINITIONS

template<class Vertex>
class QuadMesh {
public:
     QuadMesh() {}

     //NOTE: Indices all have to be arranged counterclockwise.
     QuadMesh(std::vector<Vertex> vertices, Indices indices)
          : m_indices(std::move(indices)), m_vertices(std::move(vertices)) {}

     const Indices& getIndices() const { return m_indices; }

     Indices& getIndices() { return m_indices; }

     hpuint getNumberOfQuads() const { return (m_indices.size() >> 2); }

     hpuint getNumberOfVertices() const { return m_vertices.size(); }//TODO: number of vertices on mesh may be less than the number of vertices in vector

     std::tuple<const Vertex&, const Vertex&, const Vertex&> getQuad(hpuint t) const { return std::tie(getVertex(t, 0), getVertex(t, 1), getVertex(t, 2), getVertex(t, 3)); }

     auto& getVertex(hpindex v) const { return m_vertices[v]; }

     auto& getVertex(hpindex v) { return m_vertices[v]; }

     auto& getVertex(hpindex t, hpindex i) const { return m_vertices[m_indices[(t << 2) + i]]; }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

private:
     Indices m_indices;
     std::vector<Vertex> m_vertices;

};//QuadMesh

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(std::vector<Vertex> vertices, Indices indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(const std::string& mesh) { return format::hph::read<QuadMesh<Vertex> >(mesh); }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(const std::experimental::filesystem::path& mesh) { return format::hph::read<QuadMesh<Vertex> >(mesh); }

template<class Vertex = VertexP3>
TriangleMesh<Vertex> make_triangle_mesh(const QuadMesh<Vertex>&  mesh){
     auto q_indices = mesh.getIndices();
     auto indices = Indices();
     auto vertices = mesh.getVertices();
     visit_quartets(q_indices, [&](hpuint v0, hpuint v1, hpuint v2, hpuint v3) {
          indices.push_back(v1);
          indices.push_back(v3);
          indices.push_back(v0);
          indices.push_back(v3);
          indices.push_back(v1);
          indices.push_back(v2);
     });
     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

template<class Vertex>
hpuint size(const QuadMesh<Vertex>& mesh) { return mesh.getNumberOfQuads(); }

}//namespace happah
