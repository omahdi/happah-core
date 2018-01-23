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
#include "happah/geometry/TriangleArray.hpp"
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
auto deindex(const QuadMesh<Vertex>& mesh);

template<class Vertex, class Point = typename Vertex::SPACE::POINT>
std::tuple<Point, Point> make_axis_aligned_bounding_box(const QuadMesh<Vertex>& mesh);

Quadruples<quax> make_neighbors(const Quadruples<hpindex>& indices);

template<class Vertex>
Quadruples<quax> make_neighbors(const QuadMesh<Vertex>& mesh);

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(std::vector<Vertex> vertices, Quadruples<hpindex> indices);

//Convert a string representation in HPH format.
template<class Vertex = VertexP3>
QuadMesh<Vertex> make_quad_mesh(const std::string& mesh);

//Import data stored in the given file in HPH format.
template<class Vertex = VertexP3>
QuadMesh<Vertex> make_quad_mesh(const std::experimental::filesystem::path& mesh);

inline qum::SpokesEnumerator make_spokes_enumerator(const Quadruples<quax>& neighbors, quax x);

template<class Vertex>
qum::SpokesEnumerator make_spokes_enumerator(const QuadMesh<Vertex>& mesh, const Quadruples<quax>& neighbors, hpindex v);

inline qum::SpokesWalker make_spokes_walker(const Quadruples<quax>& neighbors, quax x);

template<class Vertex>
qum::SpokesWalker make_spokes_walker(const QuadMesh<Vertex>& mesh, const Quadruples<quax>& neighbors, hpindex v);

template<class Vertex>
TriangleArray<Vertex> make_triangle_array(const QuadMesh<Vertex>& mesh);

template<class Vertex>
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

     auto& getVertex(quax x) const { return m_vertices[m_indices[x]]; }

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
     SpokesWalker(const Quadruples<quax>& neighbors, quax x)
          : m_neighbors(neighbors), m_x(x) {}

     auto operator==(const SpokesWalker& walker) const { return m_x == walker.m_x; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return m_x; }

     auto& operator++() {
          m_x = m_neighbors[m_x.getPrevious()];
          return *this;
     }

     auto& operator--() {
          m_x = m_neighbors[m_x].getNext();
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
     const Quadruples<quax>& m_neighbors;
     quax m_x;

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
auto deindex(const QuadMesh<Vertex>& mesh) { return deindex(mesh.getVertices(), mesh.getIndices()); }

template<class Vertex, class Point>
std::tuple<Point, Point> make_axis_aligned_bounding_box(const QuadMesh<Vertex>& mesh) { return make_axis_aligned_bounding_box(mesh.getVertices()); }

template<class Vertex>
Quadruples<quax> make_neighbors(const QuadMesh<Vertex>& mesh) { return make_neighbors(mesh.getIndices()); }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(std::vector<Vertex> vertices, Quadruples<hpindex> indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(const std::string& mesh) { return format::hph::read<QuadMesh<Vertex> >(mesh); }

template<class Vertex>
QuadMesh<Vertex> make_quad_mesh(const std::experimental::filesystem::path& mesh) { return format::hph::read<QuadMesh<Vertex> >(mesh); }

inline qum::SpokesEnumerator make_spokes_enumerator(const Quadruples<quax>& neighbors, quax x) { return { make_spokes_walker(neighbors, x) }; }

template<class Vertex>
qum::SpokesEnumerator make_spokes_enumerator(const QuadMesh<Vertex>& mesh, const Quadruples<quax>& neighbors, hpindex v) { return { make_spokes_walker(mesh, neighbors, v) }; }

inline qum::SpokesWalker make_spokes_walker(const Quadruples<quax>& neighbors, quax x) { return { neighbors, x }; }

template<class Vertex>
qum::SpokesWalker make_spokes_walker(const QuadMesh<Vertex>& mesh, const Quadruples<quax>& neighbors, hpindex v) { return { neighbors, find(mesh.getIndices(), v) }; }

template<class Vertex>
TriangleArray<Vertex> make_triangle_array(const QuadMesh<Vertex>& mesh) {
     auto vertices = std::vector<Vertex>();

     vertices.reserve(6 * size(mesh));
     visit(deindex(mesh), [&](auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
          vertices.insert(std::end(vertices), {
               vertex1, vertex3, vertex0,
               vertex3, vertex1, vertex2
          });
     });

     return make_triangle_array(std::move(vertices));
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

