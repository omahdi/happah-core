// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.11 - Hedwig Amberg    - Introduce QuadMesh.
// 2017.12 - Hedwig Amberg    - add make_neighbors.
// 2017.12 - Hedwig Amberg    - add SpokesWalker

#pragma once

#include "happah/geometry/TriangleMesh.hpp"

namespace happah {

//DECLARATIONS

template<class Vertex>
class QuadMesh;

namespace qum {
     
class SpokesWalker;

class SpokesEnumerator;

}//namespace qum

template<class Vertex>
Indices make_neighbors(const QuadMesh<Vertex>& mesh);

Indices make_neighbors(quad, const Indices& indices);

//Return the index of the ith neighbor of the qth quad.
inline quid make_neighbor_index(const Indices& neighbors, quid q, quat i);

//Assuming the qth and rth quads are neighbors, return the offset of the rth quad among the 4 neighbors of the qth quad.
quat make_neighbor_offset(const Indices& neighbors, quid q, quid r);

inline quid make_quad_index(const Indices& indices, hpindex v);

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(std::vector<Vertex> vertices, Indices indices);

//Convert a string representation in HPH format.
template<class Vertex = VertexP3>
QuadMesh<Vertex> make_quad_mesh(const std::string& mesh);

//Import data stored in the given file in HPH format.
template<class Vertex = VertexP3>
QuadMesh<Vertex> make_quad_mesh(const std::experimental::filesystem::path& mesh);

quat make_vertex_offset(const Indices& indices, quid q, hpindex v);

inline qum::SpokesEnumerator make_spokes_enumerator(const Indices& neighbors, quid q, quat i);

template<class Vertex>
qum::SpokesEnumerator make_spokes_enumerator(const QuadMesh<Vertex>& mesh, const Indices& neighbors, hpindex v);

inline qum::SpokesWalker make_spokes_walker(const Indices& neighbors, quid q, quat i);

template<class Vertex = VertexP3>
TriangleMesh<Vertex> make_triangle_mesh(const QuadMesh<Vertex>& mesh);

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

     auto& getVertex(quid q, quat i) const { return m_vertices[m_indices[(q << 2) + i]]; }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

private:
     Indices m_indices;
     std::vector<Vertex> m_vertices;

};//QuadMesh

namespace qum {

class SpokesWalker {
public:
     SpokesWalker(const Indices& neighbors, quid q, quat i) : m_i(i), m_neighbors(neighbors), m_q(q) {}

     auto operator==(const SpokesWalker& walker) const { return m_i == walker.m_i && m_q == walker.m_q; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return std::make_tuple(m_q, m_i); }

     auto& operator++() {
          static constexpr hpuint o[4] = { 3, 0, 1, 2 };

          auto q = quid(m_neighbors[(m_q << 2) + o[m_i]]);

          m_i = make_neighbor_offset(m_neighbors, q, m_q);
          m_q = q;
          return *this;
     }

     auto& operator--() {
          static constexpr hpuint o[4] = { 1, 2, 3, 0 };

          auto q = quid(m_neighbors[(m_q << 2) + m_i]);

          m_i = o[make_neighbor_offset(m_neighbors, q, m_q)];
          m_q = q;
          return *this;
     }

     auto& operator+=(hpuint n) {
          while(n--) ++(*this);
          return *this;
     }

     auto& operator-=(hpuint n) {
          while(n--) --(*this);
          return *this;
     }

     auto operator+(hpuint n) const {
          auto copy = *this;
          return copy += n;
     }

     auto operator-(hpuint n) const {
          auto copy = *this;
          return copy -= n;
     }

private:
     quat m_i;
     const Indices& m_neighbors;
     quid m_q;

};//SpokesWalker

class SpokesEnumerator {
public:
     SpokesEnumerator(SpokesWalker i) : m_begin(i), m_i(std::move(i)) {}

     explicit operator bool() const { return m_i != m_begin; }

     auto operator*() const { return *m_i; }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     SpokesWalker m_begin;
     SpokesWalker m_i;

};//SpokesEnumerator

}//namespace qum

template<class Vertex>
Indices make_neighbors(const QuadMesh<Vertex>& mesh) { return make_neighbors(quad{}, mesh.getIndices()); }

inline quid make_neighbor_index(const Indices& neighbors, quid q, quat i) { return quid(neighbors[(q << 2) + i]); }

inline quid make_quad_index(const Indices& indices, hpindex v) { return quid(std::distance(std::begin(indices), std::find(std::begin(indices), std::end(indices), v)) >> 2); };

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(std::vector<Vertex> vertices, Indices indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(const std::string& mesh) { return format::hph::read<QuadMesh<Vertex> >(mesh); }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(const std::experimental::filesystem::path& mesh) { return format::hph::read<QuadMesh<Vertex> >(mesh); }

inline qum::SpokesEnumerator make_spokes_enumerator(const Indices& neighbors, quid q, quat i) { return { { neighbors, q, i } }; }

template<class Vertex>
qum::SpokesEnumerator make_spokes_enumerator(const QuadMesh<Vertex>& mesh, const Indices& neighbors, hpindex v) {
     auto& indices = mesh.getIndices();
     auto q = make_quad_index(indices, v);
     auto i = make_vertex_offset(indices, q, v);

     return { { neighbors, q, i } };
}

inline qum::SpokesWalker make_spokes_walker(const Indices& neighbors, quid q, quat i) { return { neighbors, q, i }; }

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(const QuadMesh<Vertex>& mesh) {
     auto indices = Indices();

     visit_quartets(mesh.getIndices(), [&](auto i0, auto i1, auto i2, auto i3) {
          indices.push_back(i1);
          indices.push_back(i3);
          indices.push_back(i0);
          indices.push_back(i3);
          indices.push_back(i1);
          indices.push_back(i2);
     });

     return make_triangle_mesh(mesh.getVertices(), std::move(indices));
}

template<class Vertex>
hpuint size(const QuadMesh<Vertex>& mesh) { return mesh.getNumberOfQuads(); }

}//namespace happah

