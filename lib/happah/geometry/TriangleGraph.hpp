// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>
#include <boost/range/irange.hpp>
#include <random>
#include <string>
#include <stack>
#include <glm/gtc/constants.hpp>

#include <happah/Eigen.hpp>

#include "happah/format/hph.hpp"
#include "happah/geometry/TriangleMesh.hpp"
#include "happah/util/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

struct Edge;

template<class Vertex>
class TriangleGraph;

namespace trg {

template<class Iterator>
class EdgesEnumerator;

class SpokesWalker;

class SpokesEnumerator;

template<class Iterator>
EdgesEnumerator<Iterator> make_edges_enumerator(Iterator begin, Iterator end, hpuint nTriangles);

}//namespace trg

//Returns information about valences of branch nodes in cut.  A cut consists of multiple sequences of valence two nodes followed by a branch node and ending with a sequence of valence two nodes.  The output stores the length of the sequence of valence two nodes (which can be zero) followed by the valence of the branch node.  The last element in the output stores the length of the sequence of valence two nodes at the end of the cut (which can be zero).
template<class Vertex>
std::tuple<Indices, Indices, Indices> analyze(const TriangleGraph<Vertex>& graph, const Indices& cut);

//Returns path that cuts the mesh into a disk.
template<class Picker>
Indices cut(const Triples<Edge>& edges, hpindex t, Picker&& pick);

Indices cut(const Triples<Edge>& edges);

template<class Vertex>
Indices cut(const TriangleGraph<Vertex>& graph);

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const TriangleGraph<Vertex>& graph);
     
inline auto make_diamonds_enumerator(const Triples<Edge>& edges, hpuint nTriangles);

template<class Vertex>
auto make_diamonds_enumerator(const TriangleGraph<Vertex>& graph);

Triples<Edge> make_edges(const Triples<hpindex>& indices);

inline auto make_edges_enumerator(const Triples<Edge>& edges, hpuint nTriangles);

inline auto make_fan_enumerator(trg::SpokesEnumerator e);

inline auto make_fan_enumerator(const Triples<Edge>& edges, trix e);

template<class Vertex>
auto make_fan_enumerator(const TriangleGraph<Vertex>& graph, hpindex v);

template<class Vertex>
Triples<hpindex> make_indices(const TriangleGraph<Vertex>& graph);

Triples<trix> make_neighbors(const Triples<Edge>& edges, hpuint nTriangles);

template<class Vertex>
Triples<trix> make_neighbors(const TriangleGraph<Vertex>& graph);

inline auto make_ring_enumerator(trg::SpokesEnumerator e);

inline auto make_ring_enumerator(const Triples<Edge>& edges, trix e);

template<class Vertex>
auto make_ring_enumerator(const TriangleGraph<Vertex>& graph, hpindex v);

inline trg::SpokesEnumerator make_spokes_enumerator(trg::SpokesWalker walker);

inline trg::SpokesEnumerator make_spokes_enumerator(const Triples<Edge>& edges, trix e);

template<class Vertex>
auto make_spokes_enumerator(const TriangleGraph<Vertex>& graph, hpindex v);

inline trg::SpokesWalker make_spokes_walker(const Triples<Edge>& edges, trix e);

template<class Vertex>
trg::SpokesWalker make_spokes_walker(const TriangleGraph<Vertex>& graph, hpindex v);

template<class Vertex>
TriangleGraph<Vertex> make_triangle_graph(std::vector<Vertex> vertices, const Triples<hpindex>& indices);

template<class Vertex>
TriangleGraph<Vertex> make_triangle_graph(const TriangleMesh<Vertex>& mesh);

template<class Vertex0, class Vertex1 = VertexP2, class VertexFactory = VertexFactory<Vertex1> >
TriangleGraph<Vertex1> make_triangle_graph(const TriangleGraph<Vertex0>& graph, const Indices& cut, const std::vector<Point2D>& polygon, const std::vector<Point2D>& interior, VertexFactory&& build = VertexFactory());

std::vector<Point2D> parametrize(const Indices& lengths, const std::vector<Point3D>& polyline);

//Assume polygon is convex.  Cut is an array of indices of edges ordered by their position in a linear traversal of the cut that transform the given graph into a disk.
template<class Vertex>
std::vector<Point2D> parametrize(const TriangleGraph<Vertex>& graph, const Indices& cut, const std::vector<Point2D>& polygon);

hpuint size(trg::SpokesEnumerator e);

template<class Vertex>
hpuint size(const TriangleGraph<Vertex>& graph);

/*
 * NOTE: Assume closed mesh.
 * The ith triangle in the input graph is replaced by the (4i)th, (4i+1)th, (4i+2)th, and (4i+3)th triangles in the output mesh.  The order of the output triangles is given by the diagram below.  The order of the corresponding vertices is { { 0, 1, 3 }, { 1, 2, 4 }, { 1, 4, 3 }, { 4, 5, 3 } } and is the same ordering as in the BINARY_UNIFORM triangle refinement scheme.
 *
 *   INPUT            
 *                   /\
 *                  /  \
 *                 /    \
 *                /      \
 *               /        \
 *              /          \
 *              ------------
 *
 *   OUTPUT          10
 *                   / \
 *                  /   \
 *                 11----9
 *                 8-----7
 *               2  \   /  5
 *              / \  \ /  / \
 *             /   \  6  /   \
 *            0-----1   3-----4
 */
template<class Vertex, class VertexRule, class EdgeRule>
TriangleGraph<Vertex> subdivide(alg::loop, const TriangleGraph<Vertex>& graph, VertexRule&& vertexRule, EdgeRule&& edgeRule);

template<class Vertex>
TriangleGraph<Vertex> subdivide(const TriangleGraph<Vertex>& graph);

//Remove any useless branches in a path.
template<class Vertex>
Indices trim(const TriangleGraph<Vertex>& graph, Indices path);

template<class Vertex>
Indices undegenerate(const TriangleGraph<Vertex>& graph, const Indices& cut);

template<class Vertex>
hpuint validate_cut(const TriangleGraph<Vertex>& graph, const Indices& path);

//Make sure the next edge is really in the next ring.
template<bool loop, class Vertex>
bool validate_path(const TriangleGraph<Vertex>& graph, const Indices& path);

//DEFINITIONS

struct Edge {
     Edge(trix id, hpindex vertex)
          : m_id(id), m_vertex(vertex) {}

     Edge(trix id, trix opposite, hpindex vertex)
          : m_id(id), m_opposite(opposite), m_vertex(vertex) {}

     auto& getId() const { return m_id; }

     auto getNext() const { return m_id.getNext(); }

     auto getOffset() const { return m_id.getOffset(); }

     auto& getOpposite() const { return m_opposite; }

     auto getPrevious() const { return m_id.getPrevious(); }

     auto getTriangle() const { return m_id.getTriple(); }

     auto getVertex() const { return m_vertex; }

     void setOpposite(trix x) { m_opposite = x; }

     void setVertex(hpindex v) { m_vertex = v; }

private:
     trix m_id;
     trix m_opposite;
     hpindex m_vertex;

};//Edge

template<class Vertex>
class TriangleGraph {
public:
     TriangleGraph() {}

     TriangleGraph(std::vector<Vertex> vertices, Triples<Edge> edges, hpuint nTriangles)
          : m_edges(std::move(edges)), m_nTriangles(nTriangles), m_outgoing(vertices.size()), m_vertices(std::move(vertices)) { std::for_each(std::begin(m_edges), std::begin(m_edges) + 3 * m_nTriangles, [&](auto& edge) { m_outgoing[edge.getVertex()] = edge.getOpposite(); }); }

     auto& getEdge(hpindex e) const { return m_edges[e]; }

     auto& getEdge(hpindex t, trit i) const { return getEdge(3 * t + i); }

     auto& getEdges() const { return m_edges; }

     hpuint getNumberOfEdges() const { return m_edges.size(); }

     hpuint getNumberOfTriangles() const { return m_nTriangles; }

     hpuint getNumberOfVertices() const { return m_vertices.size(); }//TODO: number of vertices in vector may be greater than number of vertices in graph

     auto& getOutgoing() const { return m_outgoing; }

     auto getOutgoing(hpindex v) const { return m_outgoing[v]; }

     auto getOutgoing(trix x) const { return getOutgoing(x.getTriple(), x.getOffset()); }

     auto getOutgoing(hpindex t, trit i) const {
          static const trit o[3] = { TRIT2, TRIT0, TRIT1 };

          return m_outgoing[getEdge(t, o[i]).getVertex()];
     }

     auto getTriangle(hpindex t) const { return std::tie(getVertex(t, TRIT0), getVertex(t, TRIT1), getVertex(t, TRIT2)); }

     auto& getVertex(hpindex v) const { return m_vertices[v]; }

     auto& getVertex(hpindex v) { return m_vertices[v]; }

     auto& getVertex(trix x) const { return getVertex(x.getTriple(), x.getOffset()); }

     auto& getVertex(trix x) { return getVertex(x.getTriple(), x.getOffset()); }

     auto& getVertex(hpindex t, trit i) const {
          static const trit o[3] = { TRIT2, TRIT0, TRIT1 };

          return getVertex(getEdge(t, o[i]).getVertex());
     }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

private:
     Triples<Edge> m_edges;
     hpuint m_nTriangles;
     std::vector<trix> m_outgoing;
     std::vector<Vertex> m_vertices;

};//TriangleGraph

namespace trg {

template<class Iterator>
class EdgesEnumerator {
public:
     EdgesEnumerator(Iterator begin, Iterator end, hpuint nTriangles)
          : m_i(begin), m_e(0), m_end(end), m_nTriangles(nTriangles) {}

     explicit operator bool() const { return m_i != m_end; }

     auto& operator*() const { return *m_i; }

     auto& operator++() {
          do ++m_e; while(++m_i != m_end && (*m_i).getOpposite() < m_e && (*m_i).getOpposite() < 3 * m_nTriangles);
          return *this;
     }

private:
     Iterator m_i;
     hpindex m_e;
     Iterator m_end;
     hpuint m_nTriangles;

};//EdgesEnumerator

class SpokesWalker {
public:
     SpokesWalker(const Triples<Edge>& edges, trix e)
          : m_e(e), m_edges(edges) {}

     auto operator==(const SpokesWalker& walker) const { return m_e == walker.m_e; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return m_e; }

     auto& operator++() {
          m_e = m_edges[m_e.getPrevious()].getOpposite();
          return *this;
     }

     auto& operator--() {
          m_e = m_edges[m_e].getOpposite().getNext();
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
     trix m_e;
     const Triples<Edge>& m_edges;

};//SpokesWalker

/*
 *   auto e = make_..._enumerator(...);
 *   do {
 *        ...
 *   } while(++e);
 */

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

template<class Iterator>
EdgesEnumerator<Iterator> make_edges_enumerator(Iterator begin, Iterator end, hpuint nTriangles) { return { begin, end, nTriangles }; }

}//namespace trg

template<class Vertex>
std::tuple<Indices, Indices, Indices> analyze(const TriangleGraph<Vertex>& graph, const Indices& cut) {
     auto indices = Indices();
     auto valences = Indices();
     auto pairings = Indices();
     auto& edges = graph.getEdges();
     auto cache = std::vector<hpuint>(graph.getNumberOfVertices(), 0);
     auto n = hpindex(0);

     for(auto e : cut) ++cache[edges[e].getVertex()];
     for(auto e : cut) {
          auto valence = cache[edges[e].getVertex()];
          if(valence > hpuint(2)) {
               valences.push_back(valence);
               indices.push_back(n);
          }
          ++n;
     }

     cache.assign(graph.getNumberOfVertices(), std::numeric_limits<hpuint>::max());
     for(auto i = std::begin(indices), end = std::end(indices) - 1; i != end; ++i) {
          auto& edge = edges[cut[*(i + 1)]];
          auto& begin = cache[edge.getVertex()];
          auto k = begin, j = hpindex(0);

          while(k != std::numeric_limits<hpuint>::max() && edge.getOpposite() != cut[indices[k] + 1]) {
               j = k;
               k = pairings[k];
          }
          if(k == std::numeric_limits<hpuint>::max()) {
               auto& temp = cache[edges[cut[*i]].getVertex()];
               pairings.push_back(temp);
               temp = pairings.size() - 1;
          } else {
               if(k == begin) begin = pairings[k];
               else pairings[j] = pairings[k];
               pairings[k] = pairings.size();
               pairings.push_back(k);
          }
     }

     auto k = cache[edges[cut[indices.front()]].getVertex()];
     assert(k != std::numeric_limits<hpuint>::max());
     pairings[k] = pairings.size();
     pairings.push_back(k);

     return std::make_tuple(std::move(valences), std::move(indices), std::move(pairings));
}

template<class Picker>
Indices cut(const Triples<Edge>& edges, hpindex t, Picker&& pick) {
     auto neighbors = Indices(edges.size() << 1, std::numeric_limits<hpindex>::max());

     neighbors[6 * t + 0] = 6 * t + 2;
     neighbors[6 * t + 1] = 6 * t + 1;
     neighbors[6 * t + 2] = 6 * t + 0;
     neighbors[6 * t + 3] = 6 * t + 2;
     neighbors[6 * t + 4] = 6 * t + 1;
     neighbors[6 * t + 5] = 6 * t + 0;

     for(auto e = pick(neighbors); e != std::numeric_limits<hpindex>::max(); e = pick(neighbors)) {
          auto& edge0 = edges[e];
          auto& edge1 = edges[edge0.getOpposite()];
          auto p = neighbors[e << 1];
          auto n = neighbors[(e << 1) + 1];

          assert(neighbors[e << 1] != std::numeric_limits<hpindex>::max());
          assert(neighbors[edge0.getOpposite() << 1] == std::numeric_limits<hpindex>::max());
          assert(neighbors[edge1.getNext() << 1] == std::numeric_limits<hpindex>::max());
          assert(neighbors[edge1.getPrevious() << 1] == std::numeric_limits<hpindex>::max());

          neighbors[(p << 1) + 1] = edge1.getNext();
          neighbors[edge1.getNext() << 1] = p;
          neighbors[(edge1.getNext() << 1) + 1] = edge1.getPrevious();
          neighbors[edge1.getPrevious() << 1] = edge1.getNext();
          neighbors[(edge1.getPrevious() << 1) + 1] = n;
          neighbors[n << 1] = edge1.getPrevious();
          neighbors[e << 1] = std::numeric_limits<hpindex>::max();
     }

     auto temp = std::begin(neighbors);
     while(*temp == std::numeric_limits<hpindex>::max()) temp += 2;
     auto begin = temp[1];
     auto e = begin;
     auto i = std::begin(neighbors);
     do {
          *i = e;
          i += 2;
     } while((e = neighbors[(e << 1) + 1]) != begin);
     for(auto j = std::begin(neighbors); j != i; j += 2) j[1] = std::numeric_limits<hpindex>::max();
     neighbors.resize(std::distance(std::begin(neighbors), defrag(std::begin(neighbors), i)));
     neighbors.shrink_to_fit();

     return neighbors;
}

template<class Vertex>
Indices cut(const TriangleGraph<Vertex>& graph) { return cut(graph.getEdges()); }

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const TriangleGraph<Vertex>& graph) { return make_axis_aligned_bounding_box(graph.getVertices()); }

inline auto make_edges_enumerator(const Triples<Edge>& edges, hpuint nTriangles) { return trg::make_edges_enumerator(std::begin(edges), std::end(edges), nTriangles); }

inline auto make_diamonds_enumerator(const Triples<Edge>& edges, hpuint nTriangles) { return transform(make_edges_enumerator(edges, nTriangles), [&](auto& edge) { return std::make_tuple(std::ref(edge), edge.getVertex(), edges[edge.getNext()].getVertex(), edges[edge.getPrevious()].getVertex(), edges[edge.getOpposite().getNext()].getVertex()); }); }

template<class Vertex>
auto make_diamonds_enumerator(const TriangleGraph<Vertex>& graph) { return transform(make_diamonds_enumerator(graph.getEdges(), size(graph)), [&](auto& edge, auto v0, auto v1, auto v2, auto v3) { return std::make_tuple(std::ref(edge), std::ref(graph.getVertex(v0)), std::ref(graph.getVertex(v1)), std::ref(graph.getVertex(v2)), std::ref(graph.getVertex(v3))); }); }

inline auto make_fan_enumerator(trg::SpokesEnumerator e) { return transform(std::move(e), [&](auto e) { return e.getTriple(); }); }

inline auto make_fan_enumerator(const Triples<Edge>& edges, trix e) { return make_fan_enumerator(make_spokes_enumerator(edges, e)); }

template<class Vertex>
auto make_fan_enumerator(const TriangleGraph<Vertex>& graph, hpindex v) { return make_fan_enumerator(make_spokes_enumerator(make_spokes_walker(graph, v))); }

template<class Vertex>
Triples<hpindex> make_indices(const TriangleGraph<Vertex>& graph) {
     auto indices = Triples<hpindex>();
     auto nTriangles = size(graph);

     indices.reserve(3 * nTriangles);
     visit(graph.getEdges(), nTriangles, [&](auto& edge0, auto& edge1, auto& edge2) {
          indices.emplace_back(edge2.getVertex());
          indices.emplace_back(edge0.getVertex());
          indices.emplace_back(edge1.getVertex());
     });

     return indices;
}

template<class Vertex>
Triples<trix> make_neighbors(const TriangleGraph<Vertex>& graph) { return make_neighbors(graph.getEdges(), size(graph)); }

inline auto make_ring_enumerator(trg::SpokesEnumerator e) { return transform(std::move(e), [&](auto e) {
     static const trit o[3] = { TRIT1, TRIT2, TRIT0 };

     return trix(e.getTriple(), o[e.getOffset()]);
}); }

inline auto make_ring_enumerator(const Triples<Edge>& edges, trix e) { return make_ring_enumerator(make_spokes_enumerator(edges, e)); }

template<class Vertex>
auto make_ring_enumerator(const TriangleGraph<Vertex>& graph, hpindex v) { return transform(make_ring_enumerator(make_spokes_enumerator(make_spokes_walker(graph, v))), [&](auto x) -> auto& { return graph.getVertex(x); }); }

inline trg::SpokesEnumerator make_spokes_enumerator(trg::SpokesWalker walker) { return { std::move(walker) }; }

inline trg::SpokesEnumerator make_spokes_enumerator(const Triples<Edge>& edges, trix e) { return { make_spokes_walker(edges, e) }; }

template<class Vertex>
auto make_spokes_enumerator(const TriangleGraph<Vertex>& graph, hpindex v) { return transform(make_spokes_enumerator(make_spokes_walker(graph, v)), [&](auto e) { return graph.getEdge(e); }); }

inline trg::SpokesWalker make_spokes_walker(const Triples<Edge>& edges, trix e) { return { edges, e }; }

template<class Vertex>
trg::SpokesWalker make_spokes_walker(const TriangleGraph<Vertex>& graph, hpindex v) { return { graph.getEdges(), graph.getOutgoing(v) }; }

template<class Vertex>
TriangleGraph<Vertex> make_triangle_graph(std::vector<Vertex> vertices, const Triples<hpindex>& indices) { return { std::move(vertices), make_edges(indices), hpuint(indices.size() / 3) }; }

template<class Vertex>
TriangleGraph<Vertex> make_triangle_graph(const TriangleMesh<Vertex>& mesh) { return make_triangle_graph(mesh.getVertices(), mesh.getIndices()); }

template<class Vertex0, class Vertex1, class VertexFactory>
TriangleGraph<Vertex1> make_triangle_graph(const TriangleGraph<Vertex0>& graph, const Indices& cut, const std::vector<Point2D>& polygon, const std::vector<Point2D>& interior, VertexFactory&& build) {
     auto edges = graph.getEdges();
     auto vertices = std::vector<Vertex1>();
     auto p = Indices(graph.getNumberOfVertices(), hpindex(0));
     auto n = hpindex(0);

     auto lambda = [&](auto e, auto f) {
          auto& edge = edges[e];
          auto walker = make_spokes_walker(edges, edge.getNext());

          edge.setOpposite(trix());
          while(*walker != f) {
               edges[edges[*walker].getPrevious()].setVertex(n);
               --walker;
          }
          edges[edges[*walker].getPrevious()].setVertex(n++);
     };

     for(auto& e : cut) p[edges[e].getVertex()] = std::numeric_limits<hpindex>::max();
     for(auto& v : p) if(v == hpindex(0)) v = n++;

     vertices.reserve(interior.size() + polygon.size());
     for(auto& point : interior) vertices.push_back(build(point));
     for(auto& point : polygon) vertices.push_back(build(point));
     for(auto& edge : edges) edge.setVertex(p[edge.getVertex()]);
     for(auto i = std::begin(cut), end = std::end(cut) - 1; i != end; ++i) lambda(i[0], i[1]);
     lambda(cut.back(), cut.front());

     return { std::move(vertices), std::move(edges), size(graph) };
}

template<class Vertex>
std::vector<Point2D> parametrize(const TriangleGraph<Vertex>& graph, const Indices& cut, const std::vector<Point2D>& polygon) {
     using Vector = Eigen::Matrix<hpreal, Eigen::Dynamic, 1>;

     auto& edges = graph.getEdges();
     auto sums = std::vector<hpreal>(graph.getNumberOfVertices(), hpreal(0));
     auto lambdas = std::vector<hpreal>(edges.size(), hpreal(0));
     auto l = std::begin(lambdas) - 1;
     auto m = std::begin(lambdas) - 1;
     auto a = std::vector<Eigen::Triplet<hpreal> >();
     auto p = Indices(graph.getNumberOfVertices(), hpuint(0));
     auto n = hpuint(0);

     for(auto& e : cut) p[edges[e].getVertex()] = std::numeric_limits<hpuint>::max();
     for(auto& v : p) if(v == hpuint(0)) {
          a.emplace_back(n, n, hpreal(1));
          v = n++;
     }

     visit(edges, [&](auto& edge0, auto& edge1, auto& edge2) {
          auto v1 = edge0.getVertex();
          auto v2 = edge1.getVertex();
          auto v0 = edge2.getVertex();
          auto& point0 = graph.getVertex(v0).position;
          auto& point1 = graph.getVertex(v1).position;
          auto& point2 = graph.getVertex(v2).position;
          auto l0 = glm::length2(point1 - point0);
          auto l1 = glm::length2(point2 - point1);
          auto l2 = glm::length2(point0 - point2);
          auto m0 = glm::sqrt(l0);
          auto m1 = glm::sqrt(l1);
          auto m2 = glm::sqrt(l2);
          auto w0 = std::tan(std::acos((l0 + l2 - l1) / (hpreal(2) * m0 * m2)) / hpreal(2));
          auto w1 = std::tan(std::acos((l1 + l0 - l2) / (hpreal(2) * m1 * m0)) / hpreal(2));
          auto w2 = std::tan(std::acos((l2 + l1 - l0) / (hpreal(2) * m2 * m1)) / hpreal(2));
          auto x0 = w0 / m0;
          auto x1 = w1 / m1;
          auto x2 = w2 / m2;
          auto y0 = w0 / m2;
          auto y1 = w1 / m0;
          auto y2 = w2 / m1;

          *(++l) += x0;
          sums[v0] += x0 + y0;
          *(++l) += x1;
          sums[v1] += x1 + y1;
          *(++l) += x2;
          sums[v2] += x2 + y2;
          lambdas[edge0.getOpposite()] += y1;
          lambdas[edge1.getOpposite()] += y2;
          lambdas[edge2.getOpposite()] += y0;
     });

     auto bx = Vector(Vector::Zero(n));
     auto by = Vector(Vector::Zero(n));

     for(auto& edge : edges) {
          auto v = edges[edge.getOpposite()].getVertex();
          auto i = p[v];
          auto j = p[edge.getVertex()];
          auto lambda = *(++m) / sums[v];
          if(i == std::numeric_limits<hpuint>::max()) continue;
          if(j == std::numeric_limits<hpuint>::max()) {
               auto walker = make_spokes_walker(edges, edge.getOpposite());
               while(std::find(std::begin(cut), std::end(cut), *(++walker)) == std::end(cut));
               auto k = std::distance(std::begin(cut), std::find(std::begin(cut), std::end(cut), edges[*walker].getOpposite()));
               auto& point = polygon[k];
               bx[i] += lambda * point.x;
               by[i] += lambda * point.y;
          } else a.emplace_back(i, j, -lambda);
     }

     auto A = make_sparse_matrix(n, n, a);
     
     Eigen::SparseLU<Eigen::SparseMatrix<hpreal> > solver;
     solver.analyzePattern(A);
     solver.factorize(A);

     auto x = Vector(solver.solve(bx));
     auto y = Vector(solver.solve(by));
     auto points = std::vector<Point2D>();

     points.reserve(n);

     for(auto i = hpuint(0); i < n; ++i) points.emplace_back(x[i], y[i]);

     return points;
}

template<class Vertex>
hpuint size(const TriangleGraph<Vertex>& graph) { return graph.getNumberOfTriangles(); }

template<class Vertex, class VertexRule, class EdgeRule>
TriangleGraph<Vertex> subdivide(alg::loop, const TriangleGraph<Vertex>& graph, VertexRule&& vertexRule, EdgeRule&& edgeRule) {
     auto vertices = std::vector<Vertex>();
     auto edges = Triples<Edge>();
     auto vs = Indices(graph.getNumberOfEdges(), std::numeric_limits<hpindex>::max());

     edges.reserve((size(graph) << 2) * 3);
     vertices.reserve(graph.getNumberOfVertices() + (graph.getNumberOfEdges() >> 1));

     auto v = 0u;
     for(auto& vertex : graph.getVertices()) {
          auto ring = make(make_ring_enumerator(graph, v++));

          vertices.emplace_back(vertexRule(vertex, std::begin(ring), std::end(ring)));
     }

     visit(make_diamonds_enumerator(graph), [&](auto e, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
          auto& edge = graph.getEdge(e);

          vs[e] = vertices.size();
          vs[edge.getOpposite()] = vertices.size();
          vertices.emplace_back(edgeRule(vertex0, vertex2, vertex1, vertex3));
     });

     //TODO
     /*auto i = std::begin(vs) - 1;
     for(auto j = std::begin(edges), end = std::begin(edges) + 3 * size(graph); j != end; j += 3) {
          auto v0 = j[2].vertex;
          auto v2 = j[0].vertex;
          auto v5 = j[1].vertex;
          auto v1 = *(++i);
          auto v4 = *(++i);
          auto v3 = *(++i);

          indices.insert(std::end(indices), {
               v0, v1, v3,
               v1, v2, v4,
               v1, v4, v3,
               v4, v5, v3
          });
     }*/

     return { std::move(vertices), std::move(edges), size(graph) << 2 };
}

template<class Vertex>
TriangleGraph<Vertex> subdivide(const TriangleGraph<Vertex>& graph) {
     return subdivide(alg::loop{}, graph, [](auto& center, auto begin, auto end) {
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

template<class Vertex>
Indices trim(const TriangleGraph<Vertex>& graph, Indices path) {
     auto i = std::begin(path);
     auto j = next(i);

     auto test = [](auto e) { return e != std::numeric_limits<hpindex>::max(); };
     auto next = [&](auto i) {
          auto j = std::find_if(i + 1, std::end(path), test);
          if(j == std::end(path)) j = std::find_if(std::begin(path), i, test);
          return j;
     };
     auto previous = [&](auto i) {
          auto j = std::find_if(std::make_reverse_iterator(i), std::rend(path), test);
          if(j == std::rend(path)) return next(std::begin(path));
          else return j.base() - 1;
     };

     while(true) {
          if(*j == graph.getEdge(*i).getOpposite()) {
               *i = std::numeric_limits<hpindex>::max();
               *j = std::numeric_limits<hpindex>::max();
               i = previous(i);
               j = next(i);
               if(j == i) break;
          } else {
               auto temp = i;
               i = next(i);
               if(std::distance(temp, i) <= 0) break;
               j = next(i);
          }
     }

     path.resize(std::distance(std::begin(path), defrag(path)));
     path.shrink_to_fit();

     return path;
}

namespace detail {

//TODO: cut should be std::vector<trix>
template<class Iterator>
void undegenerate(const Triples<Edge>& edges, const Indices& cut, hpindex e, Iterator begin, Iterator end, Indices& result) {
     while(begin != end) {
          auto walker = make_spokes_walker(edges, edges[e].getNext());

          while(*walker != *begin) {
               auto f = *walker;
               auto v = edges[f].getVertex();

               --walker;
               if(edges[*end].getVertex() == v) {
                    auto temp = make_spokes_walker(edges, edges[f].getOpposite());
                    while(edges[*(++temp)].getOpposite() != *end) if(std::find(std::begin(cut), std::end(cut), *temp) != std::end(cut)) goto leave;
                    result.push_back(f);
                    return;
                    leave:;
               }
               for(auto j = end - 1; j != begin; --j) if(edges[*j].getVertex() == v) {
                    result.push_back(f);
                    undegenerate(edges, cut, f, j + 1, end, result);
                    return;
               }
          }
          e = *begin;
          result.push_back(e);
          ++begin;
     }
     result.push_back(*end);
}

}//namespace detail

template<class Vertex>
Indices undegenerate(const TriangleGraph<Vertex>& graph, const Indices& cut) {
     auto& edges = graph.getEdges();
     auto analysis = analyze(graph, cut);
     auto& indices = std::get<1>(analysis);
     auto& pairings = std::get<2>(analysis);
     auto result = Indices();
     auto b = hpindex(0);
     auto p = std::begin(pairings);
     auto lengths = Indices();
     auto f = cut[indices[0]];

     lengths.reserve(indices.size());
     lengths.push_back(-1);
     for(auto i = std::begin(indices) + 1, end = std::end(indices); i != end; ++i, ++b, ++p) {
          if(*p < b) for(auto j = lengths[*p + 1], end = lengths[*p]; j != end; --j) result.push_back(edges[result[j]].getOpposite());
          else {
               auto i0 = *(i - 1);
               auto i1 = *i;
               auto temp = Indices();
               auto n = result.size();
               auto m = pairings[(*p == hpindex(0)) ? pairings.size() - 1 : (*p - 1)];
               auto e = edges[(m < b) ? result[lengths[m] + 1] : cut[(indices[m] == cut.size() - 1) ? 0 : indices[m] + 1]].getOpposite();

               temp.reserve(i1 - i0);
               detail::undegenerate(edges, cut, f, std::begin(cut) + (i0 + 1), std::begin(cut) + i1, temp);
               for(auto& e : temp) e = edges[e].getOpposite();
               detail::undegenerate(edges, cut, e, std::rbegin(temp), std::rend(temp) - 1, result);
               std::reverse(std::begin(result) + n, std::end(result));
               for(auto& e : boost::make_iterator_range(std::begin(result) + n, std::end(result))) e = edges[e].getOpposite();
          }
          f = result.back();
          lengths.push_back(result.size() - 1);
     }
     for(auto i = lengths[*p + 1], end = lengths[*p]; i != end; --i) result.push_back(edges[result[i]].getOpposite());

     return result;
}

//TODO: path should be std::vector<trix>?
template<class Vertex>
hpuint validate_cut(const TriangleGraph<Vertex>& graph, const Indices& path) {
     auto cache = path;
     auto todo = std::stack<hpindex>();
     auto visited = boost::dynamic_bitset<>(graph.getNumberOfTriangles(), false);

     auto contains = [&](auto e) { return std::binary_search(std::begin(cache), std::end(cache), e); };
     auto push = [&](auto e) { if(!contains(e)) todo.push(graph.getEdge(e).getOpposite().getTriple()); };

     if(path.size() == 0) return 0;
     if(!validate_path<true>(graph, path)) return 1;
     std::sort(std::begin(cache), std::end(cache));

     //Make sure each edge is exactly once on the path.
     if(std::unique(std::begin(cache), std::end(cache)) != std::end(cache)) return 2;

     //Make sure opposite edge is always on path.
     for(auto e : path) if(!contains(graph.getEdge(e).getOpposite())) return 3;

     //Make sure the path does not cross itself.
     for(auto i = std::begin(path), end = std::end(path) - 1; i != end; ++i) {
          if(graph.getEdge(i[0]).getOpposite() == i[1]) continue;
          auto e = make_spokes_walker(graph.getEdges(), graph.getEdge(i[1]).getId());
          while(!contains(*(++e)));
          if(graph.getEdge(*e).getOpposite() != i[0]) return 4;
     }

     //Make sure the path encloses a single connected region containing all the triangles.
     todo.push(0);
     while(!todo.empty()) {
          auto t = todo.top();
          todo.pop();
          if(visited[t]) continue;
          visited[t] = true;
          push(3 * t + 0);
          push(3 * t + 1);
          push(3 * t + 2);
     }
     if(!visited.all()) return 5;

     return 0;
}

template<bool loop, class Vertex>
bool validate_path(const TriangleGraph<Vertex>& graph, const Indices& path) {
     auto exists = [&](auto e0, auto e1) {
          auto e = make_spokes_enumerator(graph.getEdges(), graph.getEdge(e1).getId());
          do if(*e == e1) return true; while(++e);
          return false;
     };

     if(path.size() == 0) return true;
     for(auto i = std::begin(path), end = std::end(path) - 1; i != end; ++i) if(!exists(i[0], i[1])) return false;
     if(loop && !exists(path.back(), path.front())) return false;

     return true;
}

}//namespace happah

