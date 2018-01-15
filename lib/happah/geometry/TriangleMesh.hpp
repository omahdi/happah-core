// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <stack>
#include <string>

#include "happah/Happah.hpp"
#include "happah/format/hph.hpp"
#include "happah/util/ProxyArray.hpp"

namespace happah {

//DECLARATIONS

template<class Vertex>
class TriangleMesh;

namespace trm {

template<class Iterator>
class EdgesEnumerator;

class SpokesWalker;

class SpokesEnumerator;

class VerticesEnumerator;

template<class Iterator>
EdgesEnumerator<Iterator> make_edges_enumerator(Iterator begin, Iterator end);

}//namespace trm

template<class Vertex>
auto deindex(const TriangleMesh<Vertex>& mesh);

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const std::vector<Vertex>& vertices);

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const TriangleMesh<Vertex>& mesh);

template<class Vertex>
auto make_center(const TriangleMesh<Vertex>& mesh);
     
inline auto make_diamonds_enumerator(const Triples<hpindex>& neighbors);

template<class Vertex>
auto make_diamonds_enumerator(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors);

inline auto make_edges_enumerator(const Triples<hpindex>& neighbors);

inline auto make_fan_enumerator(trm::SpokesEnumerator e);

inline auto make_fan_enumerator(const Triples<hpindex>& neighbors, hpindex t, trit i);

template<class Vertex>
auto make_fan_enumerator(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, hpindex v);

Triples<hpindex> make_neighbors(const Triples<hpindex>& indices);

template<class Vertex>
Triples<hpindex> make_neighbors(const TriangleMesh<Vertex>& mesh);

inline auto make_ring_enumerator(trm::SpokesEnumerator e);

inline auto make_ring_enumerator(const Triples<hpindex>& neighbors, hpindex t, trit i);

template<class Vertex>
auto make_ring_enumerator(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, hpindex v);

inline trm::SpokesEnumerator make_spokes_enumerator(trm::SpokesWalker walker);

inline trm::SpokesEnumerator make_spokes_enumerator(const Triples<hpindex>& neighbors, hpindex t, trit i);

template<class Vertex>
trm::SpokesEnumerator make_spokes_enumerator(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, hpindex v);

inline trm::SpokesWalker make_spokes_walker(const Triples<hpindex>& neighbors, hpindex t, trit i);

template<class Vertex>
trm::SpokesWalker make_spokes_walker(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, hpindex v);

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(std::vector<Vertex> vertices, Triples<hpindex> indices);

//Convert a string representation in HPH format.
template<class Vertex = VertexP3>
TriangleMesh<Vertex> make_triangle_mesh(const std::string& mesh);

//Import data stored in the given file in HPH format.
template<class Vertex = VertexP3>
TriangleMesh<Vertex> make_triangle_mesh(const std::experimental::filesystem::path& mesh);

//NOTE: Border has to be sorted.
template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const Triples<hpindex>& neighbors, const Indices& border, hpindex t, const Vertex& vertex0, const Vertex& vertex1, const Vertex& vertex2, VertexFactory&& build);

template<class Vertex>
Indices make_valences(const TriangleMesh<Vertex>& mesh);

inline trm::VerticesEnumerator make_vertices_enumerator(const Triples<hpindex>& neighbors);

Triples<hpindex> seal(Triples<hpindex> neighbors);

hpuint size(trm::SpokesEnumerator e);

template<class Vertex>
hpuint size(const TriangleMesh<Vertex>& mesh);

template<class Vertex, class VertexRule, class EdgeRule>
TriangleMesh<Vertex> subdivide(alg::loop, const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, VertexRule&& vertexRule, EdgeRule&& edgeRule);

template<class Vertex>
TriangleMesh<Vertex> subdivide(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors);

//NOTE: Border has to be sorted.
template<class Visitor>
void visit(trm::SpokesWalker e, const Indices& border, Visitor&& visit);

//DEFINITIONS

template<class Vertex>
class TriangleMesh {
public:
     TriangleMesh() {}

     //NOTE: Indices all have to be arranged counterclockwise.
     TriangleMesh(std::vector<Vertex> vertices, Triples<hpindex> indices)
          : m_indices(std::move(indices)), m_vertices(std::move(vertices)) {}

     auto& getIndices() const { return m_indices; }

     auto& getIndices() { return m_indices; }

     hpuint getNumberOfTriangles() const { return m_indices.size() / 3; }

     hpuint getNumberOfVertices() const { return m_vertices.size(); }//TODO: number of vertices on mesh may be less than the number of vertices in vector

     auto getTriangle(hpindex t) const { return std::tie(getVertex(t, 0), getVertex(t, 1), getVertex(t, 2)); }

     auto& getVertex(hpindex v) const { return m_vertices[v]; }

     auto& getVertex(hpindex v) { return m_vertices[v]; }

     auto& getVertex(trix x) const { return m_vertices[m_indices[x]]; }

     auto& getVertex(hpindex t, trit i) const { return m_vertices[m_indices[3 * t + i]]; }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

private:
     Triples<hpindex> m_indices;
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

template<class Iterator>
class EdgesEnumerator {
public:
     EdgesEnumerator(Iterator begin, Iterator end)
          : m_i(begin), m_e(0), m_end(end) {}

     explicit operator bool() const { return m_i != m_end; }

     auto operator*() const {
          auto t = hpindex(m_e / 3);
          auto i = trit(m_e - 3 * t);

          return std::make_tuple(t, i);
     }

     auto& operator++() {
          do ++m_e; while(++m_i != m_end && (*m_i == std::numeric_limits<hpindex>::max() || 3 * *m_i < m_e));
          return *this;
     }

private:
     Iterator m_i;
     hpindex m_e;
     Iterator m_end;

};//EdgesEnumerator

class SpokesWalker {
public:
     SpokesWalker(const Triples<hpindex>& neighbors, hpindex t, trit i)
          : m_i(i), m_neighbors(neighbors), m_t(t) {}

     auto operator==(const SpokesWalker& walker) const { return m_i == walker.m_i && m_t == walker.m_t; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return std::make_tuple(m_t, m_i); }

     auto& operator++() {
          static const trit o[3] = { TRIT2, TRIT0, TRIT1 };

          auto t = m_neighbors(m_t, o[m_i]);

          m_i = make_offset(m_neighbors, t, m_t);
          m_t = t;
          return *this;
     }

     auto& operator--() {
          static const trit o[3] = { TRIT1, TRIT2, TRIT0 };

          auto t = m_neighbors(m_t, m_i);

          m_i = o[make_offset(m_neighbors, t, m_t)];
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
     trit m_i;
     const Triples<hpindex>& m_neighbors;
     hpindex m_t;

};//SpokesWalker

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
     VerticesEnumerator(const Triples<hpindex>& neighbors)
          : m_n(0), m_neighbors(neighbors), m_visited(neighbors.size(), false) {}

     explicit operator bool() const { return m_n < m_neighbors.size(); }

     auto operator*() const {
          auto t = hpindex(m_n / 3);
          auto i = trit(m_n - 3 * t);

          return std::make_tuple(t, i);
     }

     auto& operator++() {
          auto t = hpindex(m_n / 3);
          auto i = trit(m_n - 3 * t);

          visit(make_spokes_enumerator(m_neighbors, t, i), [&](auto t, auto i) { m_visited[3 * t + i] = true; });
          while(++m_n < m_neighbors.size() && m_visited[m_n]);
          return *this;
     }

private:
     hpindex m_n;
     const Triples<hpindex>& m_neighbors;
     boost::dynamic_bitset<> m_visited;

};//VerticesEnumerator

template<class Iterator>
EdgesEnumerator<Iterator> make_edges_enumerator(Iterator begin, Iterator end) { return { begin, end }; }

}//namespace trm

template<class Vertex>
auto deindex(const TriangleMesh<Vertex>& mesh) { return deindex(mesh.getVertices(), mesh.getIndices()); }

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

template<class Vertex>
auto make_center(const TriangleMesh<Vertex>& mesh) {
     using Point = typename Vertex::SPACE::POINT;

     auto center = Point(0);
     for(auto& vertex : mesh.getVertices()) center += vertex.position;
     center /= mesh.getNumberOfVertices();

     return center;
}

inline auto make_edges_enumerator(const Triples<hpindex>& neighbors) { return trm::make_edges_enumerator(std::begin(neighbors), std::end(neighbors)); }

inline auto make_diamonds_enumerator(const Triples<hpindex>& neighbors) { return transform(make_edges_enumerator(neighbors), [&](auto t, auto i) {
     static const trit o0[3] = { TRIT2, TRIT0, TRIT1 };
     static const trit o1[3] = { TRIT1, TRIT2, TRIT0 };

     auto u = neighbors(t, i);
     auto j = make_offset(neighbors, u, t);

     return std::make_tuple(t, i, t, i, u, o0[j], t, o1[i], t, o0[i]);
}); }

template<class Vertex>
auto make_diamonds_enumerator(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors) { return transform(make_diamonds_enumerator(neighbors), [&](auto t, auto i, auto t0, auto i0, auto t1, auto i1, auto t2, auto i2, auto t3, auto i3) { return std::make_tuple(t, i, std::ref(mesh.getVertex(t0, i0)), std::ref(mesh.getVertex(t1, i1)), std::ref(mesh.getVertex(t2, i2)), std::ref(mesh.getVertex(t3, i3))); }); }

template<class Vertex>
Triples<hpindex> make_neighbors(const TriangleMesh<Vertex>& mesh) { return make_neighbors(mesh.getIndices()); }

inline auto make_fan_enumerator(trm::SpokesEnumerator e) { return transform(std::move(e), [&](auto t, auto i) { return t; }); }

inline auto make_fan_enumerator(const Triples<hpindex>& neighbors, hpindex t, trit i) { return make_fan_enumerator(make_spokes_enumerator(neighbors, t, i)); }

template<class Vertex>
auto make_fan_enumerator(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, hpindex v) { return make_fan_enumerator(make_spokes_enumerator(mesh, neighbors, v)); }

inline auto make_ring_enumerator(trm::SpokesEnumerator e) { return transform(std::move(e), [&](auto t, auto i) {
     static const trit o[3] = { TRIT1, TRIT2, TRIT0 };

     return std::make_tuple(t, o[i]);
}); }

inline auto make_ring_enumerator(const Triples<hpindex>& neighbors, hpindex t, trit i) { return make_ring_enumerator(make_spokes_enumerator(neighbors, t, i)); }

template<class Vertex>
auto make_ring_enumerator(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, hpindex v) { return transform(make_ring_enumerator(make_spokes_enumerator(mesh, neighbors, v)), [&](auto t, auto i) { return mesh.getVertex(t, i); }); }

inline trm::SpokesEnumerator make_spokes_enumerator(trm::SpokesWalker walker) { return { std::move(walker) }; }

inline trm::SpokesEnumerator make_spokes_enumerator(const Triples<hpindex>& neighbors, hpindex t, trit i) { return { make_spokes_walker(neighbors, t, i) }; }

template<class Vertex>
trm::SpokesEnumerator make_spokes_enumerator(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, hpindex v) { return { make_spokes_walker(mesh, neighbors, v) }; }

inline trm::SpokesWalker make_spokes_walker(const Triples<hpindex>& neighbors, hpindex t, trit i) { return { neighbors, t, i }; }

template<class Vertex>
trm::SpokesWalker make_spokes_walker(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, hpindex v) {
     auto& indices = mesh.getIndices();
     auto x = find(indices, v);

     return { neighbors, x.getTriple(), x.getOffset() };
}

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(std::vector<Vertex> vertices, Triples<hpindex> indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(const std::string& mesh) { return format::hph::read<TriangleMesh<Vertex> >(mesh); }

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(const std::experimental::filesystem::path& mesh) { return format::hph::read<TriangleMesh<Vertex> >(mesh); }

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const Triples<hpindex>& neighbors, const Indices& border, hpindex t, const Vertex& vertex0, const Vertex& vertex1, const Vertex& vertex2, VertexFactory&& build) {
     auto vertices = std::vector<Vertex>();
     auto indices = Triples<hpindex>(neighbors.size(), std::numeric_limits<hpindex>::max());
     auto todo = std::stack<trix>();
     auto visited = boost::dynamic_bitset<>(3 * neighbors.size(), false);

     auto push = [&](auto vertex, auto t, auto i) {
          auto n = vertices.size();

          vertices.push_back(vertex);
          visit(make_spokes_walker(neighbors, t, i), border, [&](auto t, auto i) { indices[3 * t + i] = n; });
     };

     assert(*std::max_element(std::begin(neighbors), std::end(neighbors)) < std::numeric_limits<hpuint>::max());//NOTE: Implementation assumes a closed topology.

     push(vertex0, t, TRIT0);
     push(vertex1, t, TRIT1);
     push(vertex2, t, TRIT2);
     todo.emplace(t, TRIT0);
     todo.emplace(t, TRIT1);
     todo.emplace(t, TRIT2);

     while(!todo.empty()) {
          static const trit o1[3] = { TRIT1, TRIT2, TRIT0 };
          static const trit o2[3] = { TRIT2, TRIT0, TRIT1 };

          auto e = todo.top();
          todo.pop();
          if(visited[e]) continue;
          visited[e] = true;
          if(std::binary_search(std::begin(border), std::end(border), e)) continue;
          auto u = e.getTriple();
          auto j = e.getOffset();
          auto v = neighbors(u, j);
          auto k = make_offset(neighbors, v, u);
          if(indices[3 * v + o2[k]] == std::numeric_limits<hpindex>::max()) {
               auto temp = std::begin(indices) + 3 * u;
               push(build(u, j, vertices[temp[o1[j]]], vertices[temp[j]], vertices[temp[o2[j]]]), v, trit(o2[k]));
          }
          todo.emplace(v, o1[k]);
          todo.emplace(v, o2[k]);
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

inline trm::VerticesEnumerator make_vertices_enumerator(const Triples<hpindex>& neighbors) { return { neighbors }; }

template<class Vertex>
hpuint size(const TriangleMesh<Vertex>& mesh) { return mesh.getNumberOfTriangles(); }

template<class Vertex, class VertexRule, class EdgeRule>
TriangleMesh<Vertex> subdivide(alg::loop, const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors, VertexRule&& vertexRule, EdgeRule&& edgeRule) {
     auto vertices = std::vector<Vertex>();
     auto indices = Triples<hpindex>();
     auto vs = Triples<hpindex>(3 * size(mesh), std::numeric_limits<hpindex>::max());

     vertices.reserve(mesh.getNumberOfVertices() + ((3 * size(mesh)) >> 1));
     indices.reserve((size(mesh) << 2) * 3);

     auto v = 0u;
     for(auto& vertex : mesh.getVertices()) {
          auto ring = make(make_ring_enumerator(mesh, neighbors, v++));

          vertices.emplace_back(vertexRule(vertex, std::begin(ring), std::end(ring)));
     }

     visit(make_diamonds_enumerator(mesh, neighbors), [&](auto t, auto i, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
          auto u = neighbors(t, i);
          auto j = make_offset(neighbors, u, t);

          vs(t, i) = vertices.size();
          vs(u, j) = vertices.size();
          vertices.emplace_back(edgeRule(vertex0, vertex2, vertex1, vertex3));
     });

     auto i = std::begin(vs) - 1;
     visit(mesh.getIndices(), [&](auto v0, auto v2, auto v5) {
          auto v1 = *(++i);
          auto v4 = *(++i);
          auto v3 = *(++i);

          indices.insert(std::end(indices), {
               v0, v1, v3,
               v1, v2, v4,
               v1, v4, v3,
               v4, v5, v3
          });
     });

     return { std::move(vertices), std::move(indices) };
}

template<class Vertex>
TriangleMesh<Vertex> subdivide(const TriangleMesh<Vertex>& mesh, const Triples<hpindex>& neighbors) {
     return subdivide(alg::loop{}, mesh, neighbors, [](auto& center, auto begin, auto end) {
          auto valence = std::distance(begin, end);

          auto mean = Vertex();
          while(begin != end) {
               mean.position += (*begin).position;
               ++begin;
          }
          mean.position *= 1.f / valence;

          auto alpha = 3.f / 8.f + 2.f / 8.f * (hpreal)glm::cos(2 * glm::pi<hpreal>() / valence);
          alpha = 5.f / 8.f - alpha * alpha;

          return alpha * mean.position + (1 - alpha) * center.position;
     }, [](auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) { return (3.f * (vertex0.position + vertex1.position) + (vertex2.position + vertex3.position)) / 8.f; });
}

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

}//namespace happah

