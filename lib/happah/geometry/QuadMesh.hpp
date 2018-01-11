// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.11 - Hedwig Amberg    - Introduce QuadMesh.
// 2017.12 - Hedwig Amberg    - add make_neighbors.
// 2017.12 - Hedwig Amberg    - add SpokesWalker

#pragma once

#include "happah/Happah.hpp"
#include "happah/geometry/TriangleMesh.hpp"

namespace happah {

//DECLARATIONS

template<class Vertex>
class QuadMesh;

namespace qum {
     
class SpokesWalker;

class SpokesEnumerator;

}//namespace qum

Quadruples<hpindex> make_neighbors(const Quadruples<hpindex>& indices);

template<class Vertex>
Quadruples<hpindex> make_neighbors(const QuadMesh<Vertex>& mesh);

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(std::vector<Vertex> vertices, Quadruples<hpindex> indices);

//Convert a string representation in HPH format.
template<class Vertex = VertexP3>
QuadMesh<Vertex> make_quad_mesh(const std::string& mesh);

//Import data stored in the given file in HPH format.
template<class Vertex = VertexP3>
QuadMesh<Vertex> make_quad_mesh(const std::experimental::filesystem::path& mesh);

inline qum::SpokesEnumerator make_spokes_enumerator(const Quadruples<hpindex>& neighbors, hpindex q, quat i);

template<class Vertex>
qum::SpokesEnumerator make_spokes_enumerator(const QuadMesh<Vertex>& mesh, const Quadruples<hpindex>& neighbors, hpindex v);

inline qum::SpokesWalker make_spokes_walker(const Quadruples<hpindex>& neighbors, hpindex q, quat i);

template<class Vertex>
qum::SpokesWalker make_spokes_walker(const QuadMesh<Vertex>& mesh, const Quadruples<hpindex>& neighbors, hpindex v);

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
     QuadMesh(std::vector<Vertex> vertices, Quadruples<hpindex> indices)
          : m_indices(std::move(indices)), m_vertices(std::move(vertices)) {}

     auto& getIndices() const { return m_indices; }

     auto& getIndices() { return m_indices; }

     hpuint getNumberOfQuads() const { return (m_indices.size() >> 2); }

     hpuint getNumberOfVertices() const { return m_vertices.size(); }//TODO: number of vertices on mesh may be less than the number of vertices in vector

     auto getQuad(hpindex q) const { return std::tie(getVertex(q, 0), getVertex(q, 1), getVertex(q, 2), getVertex(q, 3)); }

     auto& getVertex(hpindex v) const { return m_vertices[v]; }

     auto& getVertex(hpindex v) { return m_vertices[v]; }

     auto& getVertex(hpindex q, quat i) const { return m_vertices[m_indices(q, i)]; }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

private:
     Quadruples<hpindex> m_indices;
     std::vector<Vertex> m_vertices;

};//QuadMesh

namespace qum {

class SpokesWalker {
public:
     SpokesWalker(const Quadruples<hpindex>& neighbors, hpindex q, quat i) : m_i(i), m_neighbors(neighbors), m_q(q) {}

     auto operator==(const SpokesWalker& walker) const { return m_i == walker.m_i && m_q == walker.m_q; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return std::make_tuple(m_q, m_i); }

     auto& operator++() {
          static const quat o[4] = { QUAT3, QUAT0, QUAT1, QUAT2 };

          auto q = m_neighbors(m_q, o[m_i]);

          m_i = make_offset(m_neighbors, q, m_q);
          m_q = q;
          return *this;
     }

     auto& operator--() {
          static const quat o[4] = { QUAT1, QUAT2, QUAT3, QUAT0 };

          auto q = m_neighbors(m_q, m_i);

          m_i = o[make_offset(m_neighbors, q, m_q)];
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
     const Quadruples<hpindex>& m_neighbors;
     hpindex m_q;

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
Quadruples<hpindex> make_neighbors(const QuadMesh<Vertex>& mesh) { return make_neighbors(mesh.getIndices()); }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(std::vector<Vertex> vertices, Quadruples<hpindex> indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(const std::string& mesh) { return format::hph::read<QuadMesh<Vertex> >(mesh); }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(const std::experimental::filesystem::path& mesh) { return format::hph::read<QuadMesh<Vertex> >(mesh); }

inline qum::SpokesEnumerator make_spokes_enumerator(const Quadruples<hpindex>& neighbors, hpindex q, quat i) { return { make_spokes_walker(neighbors, q, i) }; }

template<class Vertex>
qum::SpokesEnumerator make_spokes_enumerator(const QuadMesh<Vertex>& mesh, const Quadruples<hpindex>& neighbors, hpindex v) { return { make_spokes_walker(mesh, neighbors, v) }; }

inline qum::SpokesWalker make_spokes_walker(const Quadruples<hpindex>& neighbors, hpindex q, quat i) { return { neighbors, q, i }; }

template<class Vertex>
qum::SpokesWalker make_spokes_walker(const QuadMesh<Vertex>& mesh, const Quadruples<hpindex>& neighbors, hpindex v) {
     auto& indices = mesh.getIndices();
     auto q = make_index(indices, v);
     auto i = make_offset(indices, q, v);

     return { neighbors, q, i };
}

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(const QuadMesh<Vertex>& mesh) {
     auto indices = Triples<hpindex>();

     visit(mesh.getIndices(), [&](auto i0, auto i1, auto i2, auto i3) {
          indices.insert(std::end(indices), {
               i1, i3, i0,
               i3, i1, i2
          });
     });

     return make_triangle_mesh(mesh.getVertices(), std::move(indices));
}

template<class Vertex>
hpuint size(const QuadMesh<Vertex>& mesh) { return mesh.getNumberOfQuads(); }

}//namespace happah

