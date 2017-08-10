// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <stack>
#include <string>

#include "happah/format/hph.hpp"
#include "happah/util/ProxyArray.hpp"

namespace happah {

//DECLARATIONS

template<class Vertex>
class TriangleMesh;

namespace trm {

class SpokesWalker;

class FanEnumerator;

class RingEnumerator;

class SpokesEnumerator;

class VerticesEnumerator;

}//namespace trm

bool is_neighbor(const Indices& neighbors, hpuint t, hpuint u);

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const std::vector<Vertex>& vertices);

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const TriangleMesh<Vertex>& mesh);
     
inline hpindex make_edge_offset(hpindex e);

Indices make_fan(trm::FanEnumerator e);

template<class Transformer>
auto make_fan(EnumeratorTransformer<trm::FanEnumerator, Transformer> e);

inline trm::FanEnumerator make_fan_enumerator(const Indices& neighbors, hpuint t, hpuint i);

template<class Vertex>
trm::FanEnumerator make_fan_enumerator(const TriangleMesh<Vertex>& mesh, const Indices& neighbors, hpuint v);

//Return the index of the ith neighbor of the tth triangle.
hpindex make_neighbor_index(const Indices& neighbors, hpuint t, hpuint i);

//Assuming the tth and uth triangles are neighbors, return the offset of the uth triangle among the three neighbors of the tth triangle.
hpuint make_neighbor_offset(const Indices& neighbors, hpuint t, hpuint u);

Indices make_neighbors(const Indices& indices);

template<class Vertex>
Indices make_neighbors(const TriangleMesh<Vertex>& mesh);

template<class Transformer>
auto make_ring(EnumeratorTransformer<trm::RingEnumerator, Transformer> e);

inline trm::RingEnumerator make_ring_enumerator(const Indices& neighbors, hpuint t, hpuint i);

template<class Transformer>
EnumeratorTransformer<trm::RingEnumerator, Transformer> make_ring_enumerator(const Indices& neighbors, hpuint t, hpuint i, Transformer&& transform);

template<class Vertex>
auto make_ring_enumerator(const TriangleMesh<Vertex>& mesh, const Indices& neighbors, hpuint v);

inline trm::SpokesEnumerator make_spokes_enumerator(const Indices& neighbors, hpuint t, hpuint i);

template<class Vertex>
trm::SpokesEnumerator make_spokes_enumerator(const TriangleMesh<Vertex>& mesh, const Indices& neighbors, hpuint v);

inline trm::SpokesWalker make_spokes_walker(const Indices& neighbors, hpindex t, hpindex i);

inline hpindex make_triangle_index(hpindex e);

inline hpindex make_triangle_index(const Indices& indices, hpindex v);

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(std::vector<Vertex> vertices, Indices indices);

//Convert a string representation in HPH format.
template<class Vertex = VertexP3>
TriangleMesh<Vertex> make_triangle_mesh(const std::string& mesh);

//Import data stored in the given file in HPH format.
template<class Vertex = VertexP3>
TriangleMesh<Vertex> make_triangle_mesh(const std::experimental::filesystem::path& mesh);

//NOTE: Border has to be sorted.
template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const Indices& neighbors, const Indices& border, hpindex t, const Vertex& vertex0, const Vertex& vertex1, const Vertex& vertex2, VertexFactory&& build);

hpuint make_valence(trm::FanEnumerator e);

hpuint make_valence(trm::RingEnumerator e);

hpuint make_valence(trm::SpokesEnumerator e);

template<class Vertex>
Indices make_valences(const TriangleMesh<Vertex>& mesh);

hpindex make_vertex_offset(const Indices& indices, hpindex t, hpindex v);

inline trm::VerticesEnumerator make_vertices_enumerator(const Indices& neighbors);

Indices seal(Indices neighbors);

template<class Vertex>
hpuint size(const TriangleMesh<Vertex>& mesh);

//NOTE: Border has to be sorted.
template<class Visitor>
void visit(trm::SpokesWalker e, const Indices& border, Visitor&& visit);

template<class Visitor>
void visit_edges(const Indices& neighbors, Visitor&& visit);

template<class Visitor>
void visit_fan(trm::FanEnumerator e, Visitor&& visit);

template<class Visitor>
void visit_ring(trm::RingEnumerator e, Visitor&& visit);

template<class Visitor>
void visit_spokes(trm::SpokesEnumerator e, Visitor&& visit);

template<class Visitor>
void visit_vertices(const Indices& neighbors, Visitor&& visit);

//DEFINITIONS

template<class Vertex>
class TriangleMesh {
public:
     TriangleMesh() {}

     //NOTE: Indices all have to be arranged counterclockwise.
     TriangleMesh(std::vector<Vertex> vertices, Indices indices)
          : m_indices(std::move(indices)), m_vertices(std::move(vertices)) {}

     const Indices& getIndices() const { return m_indices; }

     Indices& getIndices() { return m_indices; }

     hpuint getNumberOfTriangles() const { return m_indices.size() / 3; }

     hpuint getNumberOfVertices() const { return m_vertices.size(); }//TODO: number of vertices on mesh may be less than the number of vertices in vector

     std::tuple<const Vertex&, const Vertex&, const Vertex&> getTriangle(hpuint t) const { return std::tie(getVertex(t, 0), getVertex(t, 1), getVertex(t, 2)); }

     auto& getVertex(hpindex v) const { return m_vertices[v]; }

     auto& getVertex(hpindex v) { return m_vertices[v]; }

     auto& getVertex(hpindex t, hpindex i) const { return m_vertices[m_indices[3 * t + i]]; }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

private:
     Indices m_indices;
     std::vector<Vertex> m_vertices;

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const TriangleMesh<Vertex>& mesh) {
          using happah::format::hph::operator<<;

          stream << mesh.m_vertices << '\n';
          stream << mesh.m_indices;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, TriangleMesh<Vertex>& mesh) {
          using happah::format::hph::operator>>;

          stream >> mesh.m_vertices;
          stream >> mesh.m_indices;
          return stream;
     }

};//TriangleMesh

namespace trm {

class SpokesWalker {
public:
     SpokesWalker(const Indices& neighbors, hpindex t, hpindex i)
          : m_i(i), m_neighbors(neighbors), m_t(t) {}

     auto operator==(const SpokesWalker& walker) const { return m_i == walker.m_i && m_t == walker.m_t; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return std::make_tuple(m_t, m_i); }

     auto& operator++() {
          static constexpr hpindex o[3] = { 2, 0, 1 };
          auto t = m_neighbors[3 * m_t + o[m_i]];
          m_i = make_neighbor_offset(m_neighbors, t, m_t);
          m_t = t;
          return *this;
     }

     auto& operator--() {
          static constexpr hpindex o[3] = { 1, 2, 0 };
          auto t = m_neighbors[3 * m_t + m_i];
          m_i = o[make_neighbor_offset(m_neighbors, t, m_t)];
          m_t = t;
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
     hpindex m_i;
     const Indices& m_neighbors;
     hpindex m_t;

};//SpokesWalker

/*
 *   auto e = make_..._enumerator(...);
 *   do {
 *        ...
 *   } while(++e);
 */
class FanEnumerator {
public:
     FanEnumerator(SpokesWalker i)
          : m_begin(i), m_i(std::move(i)) {}

     explicit operator bool() const { return m_i != m_begin; }

     auto operator*() const { return std::get<0>(*m_i); }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     SpokesWalker m_begin;
     SpokesWalker m_i;

};//FanEnumerator

class RingEnumerator {
public:
     RingEnumerator(SpokesWalker i)
          : m_begin(i), m_i(std::move(i)) {}

     explicit operator bool() const { return m_i != m_begin; }

     auto operator*() const {
          static constexpr hpindex o[3] = { 1, 2, 0 };
          auto t = 0u, i = 0u;
          std::tie(t, i) = *m_i;
          return std::make_tuple(t, o[i]);
     }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     SpokesWalker m_begin;
     SpokesWalker m_i;

};//RingEnumerator

class SpokesEnumerator {
public:
     SpokesEnumerator(SpokesWalker i)
          : m_begin(i), m_i(std::move(i)) {}

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

class VerticesEnumerator {
public:
     VerticesEnumerator(const Indices& neighbors)
          : m_i(0), m_neighbors(neighbors), m_visited(neighbors.size(), false) {}

     explicit operator bool() const { return m_i < m_neighbors.size(); }

     std::tuple<hpuint, hpuint> operator*() const {
          auto t = m_i / 3;
          auto i = m_i - 3 * t;
          return std::make_tuple(t, i);
     }

     auto& operator++() {
          auto t = m_i / 3;
          auto i = m_i - 3 * t;
          visit_spokes(make_spokes_enumerator(m_neighbors, t, i), [&](auto t, auto i) { m_visited[3 * t + i] = true; });
          if(i == 0) {
               if(!m_visited[++m_i]) return *this;
               if(!m_visited[++m_i]) return *this;
          } else if(i == 1) {
               if(!m_visited[++m_i]) return *this;
          }
          while(++m_i < m_neighbors.size()) {
               if(!m_visited[m_i]) return *this;
               if(!m_visited[++m_i]) return *this;
               if(!m_visited[++m_i]) return *this;
          }
          return *this;
     }

private:
     hpuint m_i;
     const Indices& m_neighbors;
     boost::dynamic_bitset<> m_visited;

};//VerticesEnumerator

}//namespace trm

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const std::vector<Vertex>& vertices) {
     auto min = Point3D(std::numeric_limits<hpreal>::min());
     auto max = Point3D(std::numeric_limits<hpreal>::min());

     for(auto& vertex : vertices) {
          if(vertex.position.x < min.x) min.x = vertex.position.x;
          if(vertex.position.y < min.y) min.y = vertex.position.y;
          if(vertex.position.z < min.z) min.z = vertex.position.z;
          
          if(vertex.position.x > max.x) max.x = vertex.position.x;
          if(vertex.position.y > max.y) max.y = vertex.position.y;
          if(vertex.position.z > max.z) max.z = vertex.position.z;
     }

     return std::make_tuple(min, max);
}

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const TriangleMesh<Vertex>& mesh) { return make_axis_aligned_bounding_box(mesh.getVertices()); }

inline hpindex make_edge_offset(hpindex e) { return e - 3 * make_triangle_index(e); }
     
template<class Transformer>
auto make_fan(EnumeratorTransformer<trm::FanEnumerator, Transformer> e) {
     using T = decltype(*e);

     auto fan = std::vector<T>();
     do fan.push_back(*e); while(++e);
     return fan;
}

inline hpindex make_neighbor_index(const Indices& neighbors, hpuint t, hpuint i) { return neighbors[3 * t + i]; }

template<class Vertex>
Indices make_neighbors(const TriangleMesh<Vertex>& mesh) { return make_neighbors(mesh.getIndices()); }

template<class Transformer>
auto make_ring(EnumeratorTransformer<trm::RingEnumerator, Transformer> e) {
     using T = decltype(*e);

     auto ring = std::vector<T>();
     do ring.push_back(*e); while(++e);
     return ring;
}

inline trm::RingEnumerator make_ring_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

inline trm::FanEnumerator make_fan_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

template<class Vertex>
trm::FanEnumerator make_fan_enumerator(const TriangleMesh<Vertex>& mesh, const Indices& neighbors, hpuint v) {
     auto& indices = mesh.getIndices();
     auto t = make_triangle_index(indices, v);
     auto i = make_vertex_offset(indices, t, v);
     return { { neighbors, t, i } };
}

template<class Transformer>
EnumeratorTransformer<trm::RingEnumerator, Transformer> make_ring_enumerator(const Indices& neighbors, hpuint t, hpuint i, Transformer&& transform) { return { make_ring_enumerator(neighbors, t, i), std::forward<Transformer>(transform) }; }

template<class Vertex>
auto make_ring_enumerator(const TriangleMesh<Vertex>& mesh, const Indices& neighbors, hpuint v) {
     auto& indices = mesh.getIndices();
     auto t = make_triangle_index(indices, v);
     auto i = make_vertex_offset(indices, t, v);
     return make_ring_enumerator(neighbors, t, i, [&](auto u, auto j) { return mesh.getVertex(u, j); });
}

inline trm::SpokesEnumerator make_spokes_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

template<class Vertex>
trm::SpokesEnumerator make_spokes_enumerator(const TriangleMesh<Vertex>& mesh, const Indices& neighbors, hpuint v) { return { { neighbors, make_triangle_index(mesh.getIndices(), v), v } }; }

inline trm::SpokesWalker make_spokes_walker(const Indices& neighbors, hpindex t, hpindex i) { return { neighbors, t, i }; }

inline hpindex make_triangle_index(hpindex e) { return e / 3; }

inline hpindex make_triangle_index(const Indices& indices, hpindex v) { return std::distance(std::begin(indices), std::find(std::begin(indices), std::end(indices), v)) / 3; }

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(std::vector<Vertex> vertices, Indices indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(const std::string& mesh) { return format::hph::read<TriangleMesh<Vertex> >(mesh); }

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(const std::experimental::filesystem::path& mesh) { return format::hph::read<TriangleMesh<Vertex> >(mesh); }

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const Indices& neighbors, const Indices& border, hpindex t, const Vertex& vertex0, const Vertex& vertex1, const Vertex& vertex2, VertexFactory&& build) {
     auto vertices = std::vector<Vertex>();
     auto indices = Indices(neighbors.size(), std::numeric_limits<hpindex>::max());
     auto todo = std::stack<hpindex>();

     auto push = [&](auto vertex, auto t, auto i) {
          auto n = vertices.size();
          vertices.push_back(vertex);
          visit(make_spokes_walker(neighbors, t, i), border, [&](auto t, auto i) { indices[3 * t + i] = n; });
     };

     assert(*std::max_element(std::begin(neighbors), std::end(neighbors)) < std::numeric_limits<hpuint>::max());//NOTE: Implementation assumes a closed topology.

     push(vertex0, t, 0);
     push(vertex1, t, 1);
     push(vertex2, t, 2);
     todo.emplace(3 * t + 0);
     todo.emplace(3 * t + 1);
     todo.emplace(3 * t + 2);

     while(!todo.empty()) {
          static constexpr hpuint o0[3] = { 0, 1, 2 };
          static constexpr hpuint o1[3] = { 2, 0, 1 };
          static constexpr hpuint o2[3] = { 1, 2, 0 };

          auto e = todo.top();
          todo.pop();
          if(std::binary_search(std::begin(border), std::end(border), e)) continue;
          auto u = make_triangle_index(e);
          auto j = make_edge_offset(e);
          auto v = make_neighbor_index(neighbors, u, j);
          auto k = make_neighbor_offset(neighbors, v, u);
          if(indices[3 * v + o1[k]] != std::numeric_limits<hpindex>::max()) continue;
          auto temp = std::begin(indices) + 3 * u;
          push(build(u, j, vertices[temp[o0[j]]], vertices[temp[o1[j]]], vertices[temp[o2[j]]]), v, o1[k]);
          todo.emplace(3 * v + o1[k]);
          todo.emplace(3 * v + o2[k]);
     }

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

template<class Vertex>
Indices make_valences(const TriangleMesh<Vertex>& mesh) {
     auto valences = Indices(mesh.getVertices().size(), 0);
     auto i = std::begin(mesh.getIndices());
     auto middle = i + 3 * size(mesh);
     auto end = std::end(mesh.getIndices());

     while(i != middle) {
          ++valences[i[0]];
          ++i;
     }
     while(i != end) {
          ++valences[i[0]];
          i += 3;
     }

     return valences;
}

template<class Vertex>
hpuint size(const TriangleMesh<Vertex>& mesh) { return mesh.getNumberOfTriangles(); }

inline trm::VerticesEnumerator make_vertices_enumerator(const Indices& neighbors) { return { neighbors }; }

template<class Visitor>
void visit(trm::SpokesWalker e, const Indices& border, Visitor&& visit) {
     auto begin = e;
     do apply(visit, *e); while((++e) != begin && !std::binary_search(std::begin(border), std::end(border), 3 * std::get<0>(*e) + std::get<1>(*e)));
     if(e == begin) return;
     while(!std::binary_search(std::begin(border), std::end(border), 3 * std::get<0>(*begin) + std::get<1>(*begin))) {
          --begin;
          apply(visit, *begin);
     }
}

template<class Visitor>
void visit_edges(const Indices& neighbors, Visitor&& visit) {
     boost::dynamic_bitset<> visited(neighbors.size(), false);

     auto do_visit_edges = [&](auto t, auto i) {
          visit(t, i);
          auto u = neighbors[3 * t + i];
          if(u == UNULL) return;
          auto j = make_neighbor_offset(neighbors, u, t);
          visited[3 * u + j] = true;
     };

     for(auto t : boost::irange(0lu, neighbors.size() / 3)) {
          if(!visited[3 * t]) do_visit_edges(t, 0);
          if(!visited[3 * t + 1]) do_visit_edges(t, 1);
          if(!visited[3 * t + 2]) do_visit_edges(t, 2);
     }
}

template<class Visitor>
void visit_fan(trm::FanEnumerator e, Visitor&& visit) { do apply(visit, *e); while(++e); }

template<class Visitor>
void visit_ring(trm::RingEnumerator e, Visitor&& visit) { do apply(visit, *e); while(++e); }

template<class Visitor>
void visit_spokes(trm::SpokesEnumerator e, Visitor&& visit) { do apply(visit, *e); while(++e); }

template<class Visitor>
void visit_vertices(const Indices& neighbors, Visitor&& visit) {
     for(auto e = make_vertices_enumerator(neighbors); e; ++e) {
          auto t = 0u, i = 0u;
          std::tie(t, i) = *e;
          visit(t, i);
     }
}

}//namespace happah

