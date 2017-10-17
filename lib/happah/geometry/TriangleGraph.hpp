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
#include <unordered_set>
#include <glm/gtc/constants.hpp>

#include <happah/Eigen.hpp>

#include "happah/format/hph.hpp"
#include "happah/geometry/TriangleMesh.hpp"

namespace happah {

//DECLARATIONS

struct Edge;

template<class Vertex>
class TriangleGraph;

namespace trg {

class SpokesWalker;

class FanEnumerator;

class RingEnumerator;

class SpokesEnumerator;

}//namespace trg

//Returns path that cuts the mesh into a disk.
template<class Picker>
Indices basic_cut(const std::vector<Edge>& edges, hpindex t, Picker&& pick);

template<class Vertex, class RandomDev = std::random_device>
Indices cut(const TriangleGraph<Vertex>& graph, hpindex t = std::numeric_limits<hpindex>::max(), RandomDev&& random = RandomDev());

template<class RandomDev = std::random_device>
Indices cut(const std::vector<Edge>& edges, hpindex t = std::numeric_limits<hpindex>::max(), RandomDev&& random = RandomDev());

//Assume polygon is convex.  Cut is an array of indices of edges ordered by their position in a linear traversal of the cut that transform the given graph into a disk.
template<class Vertex>
std::vector<Point2D> embed(const TriangleGraph<Vertex>& graph, const Indices& cut, const std::vector<Point2D>& polygon);

/*
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
TriangleMesh<Vertex> loopivide(const TriangleGraph<Vertex>& graph, VertexRule&& vertexRule, EdgeRule&& edgeRule);

template<class Vertex>
TriangleMesh<Vertex> loopivide(const TriangleGraph<Vertex>& graph);

template<class Vertex>
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const TriangleGraph<Vertex>& graph);
     
//Return the index of this edge in the edges array.
inline hpuint make_edge_index(const Edge& edge);

template<class Vertex>
boost::optional<hpindex> make_edge_index(const TriangleGraph<Vertex>& graph, hpindex v0, hpindex v1);

//Return the offset of this edge among the three edges of its adjacent triangle.
inline hpuint make_edge_offset(const Edge& edge);

std::vector<Edge> make_edges(const Indices& indices);

Indices make_fan(trg::FanEnumerator e);

template<class Transformer>
auto make_fan(EnumeratorTransformer<trg::FanEnumerator, Transformer> e);

inline trg::FanEnumerator make_fan_enumerator(const std::vector<Edge>& edges, hpuint e);

template<class Vertex>
trg::FanEnumerator make_fan_enumerator(const TriangleGraph<Vertex>& graph, hpuint v);

template<class Vertex>
Indices make_indices(const TriangleGraph<Vertex>& graph);

//Return the index of the ith neighbor of the tth triangle.
inline hpindex make_neighbor_index(const std::vector<Edge>& edges, hpuint t, hpuint i);

//Return the index of the ith neighbor of the tth triangle.
template<class Vertex>
hpuint make_neighbor_index(const TriangleGraph<Vertex>& graph, hpuint t, hpuint i);

//Assuming the tth and uth triangles are neighbors, return the offset of the uth triangle among the three neighbors of the tth triangle.
hpuint make_neighbor_offset(const std::vector<Edge>& edges, hpuint t, hpuint u);

//Assuming the tth and uth triangles are neighbors, return the offset of the uth triangle among the three neighbors of the tth triangle.
template<class Vertex>
hpuint make_neighbor_offset(const TriangleGraph<Vertex>& graph, hpuint t, hpuint u);

Indices make_neighbors(const std::vector<Edge>& edges, hpuint nTriangles);

template<class Vertex>
Indices make_neighbors(const TriangleGraph<Vertex>& graph);

template<class Transformer>
auto make_ring(EnumeratorTransformer<trg::RingEnumerator, Transformer> e);

inline trg::RingEnumerator make_ring_enumerator(const std::vector<Edge>& edges, hpuint e);

template<class Transformer>
EnumeratorTransformer<trg::RingEnumerator, Transformer> make_ring_enumerator(const std::vector<Edge>& edges, hpuint e, Transformer&& transform);

template<class Vertex>
auto make_ring_enumerator(const TriangleGraph<Vertex>& graph, hpuint v);

inline trg::SpokesEnumerator make_spokes_enumerator(const std::vector<Edge>& edges, hpuint e);

template<class Transformer>
EnumeratorTransformer<trg::SpokesEnumerator, Transformer> make_spokes_enumerator(const std::vector<Edge>& edges, hpuint e, Transformer&& transform);

template<class Vertex>
auto make_spokes_enumerator(const TriangleGraph<Vertex>& graph, hpuint v);

inline trg::SpokesWalker make_spokes_walker(const std::vector<Edge>& edges, hpindex e);

inline hpindex make_triangle_index(const Edge& edge);

template<class Vertex>
TriangleGraph<Vertex> make_triangle_graph(std::vector<Vertex> vertices, const Indices& indices);

template<class Vertex>
TriangleGraph<Vertex> make_triangle_graph(const TriangleMesh<Vertex>& mesh);

hpuint make_valence(trg::FanEnumerator e);

hpuint make_valence(trg::RingEnumerator e);

hpuint make_valence(trg::SpokesEnumerator e);

template<class Vertex>
hpuint size(const TriangleGraph<Vertex>& graph);

//Remove any useless branches in a path.
template<class Vertex>
Indices trim(const TriangleGraph<Vertex>& graph, Indices path);

template<class Vertex>
hpuint validate_cut(const TriangleGraph<Vertex>& graph, const Indices& path);

//Make sure the next edge is really in the next ring.
template<bool loop, class Vertex>
bool validate_path(const TriangleGraph<Vertex>& graph, const Indices& path);

template<class Visitor>
void visit_diamonds(const std::vector<Edge>& edges, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_diamonds(const TriangleGraph<Vertex>& graph, Visitor&& visit);

template<class Visitor>
void visit_fan(trg::FanEnumerator e, Visitor&& visit);

template<class Visitor>
void visit_ring(trg::RingEnumerator e, Visitor&& visit);

template<class Visitor>
void visit_spokes(trg::SpokesEnumerator e, Visitor&& visit);

template<class Visitor>
void visit_vertices(const std::vector<Edge>& edges, Visitor&& visit);

//DEFINITIONS

struct Edge {
     hpuint next;
     hpuint opposite;
     hpuint previous;
     hpuint vertex;//vertex to which edge points

     Edge(hpuint vertex, hpuint next, hpuint opposite, hpuint previous)
          : next(next), opposite(opposite), previous(previous), vertex(vertex) {}

};//Edge

template<class Vertex>
class TriangleGraph {
public:
     TriangleGraph() {}

     TriangleGraph(TriangleGraph&&) = default;
     TriangleGraph(const TriangleGraph&) = default;
     TriangleGraph& operator=(const TriangleGraph&) = default;
     TriangleGraph& operator=(TriangleGraph&&) = default;

     TriangleGraph(std::vector<Vertex> vertices, std::vector<Edge> edges, hpuint nTriangles)
          : m_edges(std::move(edges)), m_nTriangles(nTriangles), m_outgoing(vertices.size(), std::numeric_limits<hpindex>::max()), m_vertices(std::move(vertices)) { std::for_each(std::begin(m_edges), std::begin(m_edges) + 3 * m_nTriangles, [&](auto& edge) { m_outgoing[edge.vertex] = edge.opposite; }); }

     template<class Iterator>
     void exsect(Iterator begin, Iterator end)  {
          --end;
          while(++begin != end) {
               auto i = m_outgoing[*begin];
               auto p = *(begin - 1);
               auto n = *(begin + 1);
               while(m_edges[i].vertex != p) i = m_edges[m_edges[m_edges[i].next].next].opposite;
               i = m_edges[m_edges[m_edges[i].next].next].opposite;
               while(m_edges[i].vertex != n) {
                    splitEdge(i);
                    i = m_edges[m_edges[m_edges[i].next].next].opposite;
               }
               i = m_edges[m_edges[m_edges[i].next].next].opposite;
               while(m_edges[i].vertex != p) {
                    splitEdge(i);
                    i = m_edges[m_edges[m_edges[i].next].next].opposite;
               }
          }
     }

     auto& getEdge(hpuint e) const { return m_edges[e]; }

     auto& getEdge(hpindex t, hpindex i) const { return getEdge(3 * t + i); }

     auto& getEdges() const { return m_edges; }

     auto getNumberOfEdges() const { return m_edges.size(); }

     auto getNumberOfTriangles() const { return m_nTriangles; }

     hpuint getNumberOfVertices() const { return m_vertices.size(); }//TODO: number of vertices in vector may be greater than number of vertices in graph

     auto& getOutgoing() const { return m_outgoing; }

     auto getOutgoing(hpuint v) const { return m_outgoing[v]; }

     std::tuple<const Vertex&, const Vertex&, const Vertex&> getTriangle(hpuint t) const { return std::tie(getVertex(t, 0), getVertex(t, 1), getVertex(t, 2)); }

     auto& getVertex(hpindex v) const { return m_vertices[v]; }

     auto& getVertex(hpindex v) { return m_vertices[v]; }

     auto& getVertex(hpindex t, hpindex i) const {
          static constexpr hpindex o[3] = { 2u, 0u, 1u };

          return getVertex(getEdge(t, o[i]).vertex);
     }

     auto& getVertices() const { return m_vertices; }

     auto& getVertices() { return m_vertices; }

/**********************************************************************************
 * split edge
 *
 *            BEFORE
 *       ___    v0
 *         //  /||   \\
 *        //  / ||    \\
 *       //     ||     \\
 *      //e2    ||    e5\\
 *     //       ||       \\ |
 *    //        ||        \\|
   v2         e0||e1         v3
 *   /\\        ||        //
 *  /  \\       ||       //
 *      \\e4    ||    e3//
 *       \\     ||     //
 *        \\  / || /  //
 *         \\/  ||/  //__
 *              v1
 *
 *
 *            AFTER
 *       ___    v0
 *         //  /||   \\
 *        //  / ||    \\
 *       //e2   ||n4   \\
 *      //      ||      \\
 *     //     n2|| /   e5\\ |
 *    //   n1   ||/  n5   \\|
   v2 ==========vn========== v3
 *   /\\   n0  /||   n3   //
 *  /  \\     / ||e1     //
 *      \\e4    ||      //
 *       \\   e0||   e3//
 *        \\  / || /  //
 *         \\/  ||/  //__
 *              v1
 *
 **********************************************************************************/
     //NOTE: This works only on absolute meshes.  For relative meshes need base.
     void splitEdge(hpuint edge, hpreal u = 0.5) {
          auto border = 3 * m_nTriangles;
          auto e0 = edge;
          auto& edge0 = m_edges[e0];
          assert(edge < border && edge0.opposite < border);//TODO: edge to split is border edge or opposite is border
          auto e1 = edge0.opposite;
          auto& edge1 = m_edges[e1];
          auto e2 = edge0.next;
          auto& edge2 = m_edges[e2];
          auto e3 = edge1.next;
          auto& edge3 = m_edges[e3];
          auto e4 = edge2.next;
          auto& edge4 = m_edges[e4];
          auto e5 = edge3.next;
          auto& edge5 = m_edges[e5];

          auto v0 = edge0.vertex;
          auto& vertex0 = getVertex(v0);
          auto v1 = edge1.vertex;
          auto& vertex1 = getVertex(v1);
          auto v2 = edge2.vertex;
          auto v3 = edge3.vertex;
          hpuint vn = m_vertices.size();
          Vertex vertexn(vertex0);//TODO: improve; possibilities: vertex as parameter, VertexUtils::mix(v1,v2)
          vertexn.position = vertex0.position * u + vertex1.position * (1.0f - u);
          m_vertices.push_back(vertexn);

          hpuint n0 = m_edges.size(), n1 = n0 + 1, n2 = n1 + 1, n3 = n2 + 1, n4 = n3 + 1, n5 = n4 + 1;

          edge0.next = n0;
          edge0.vertex = vn;
          edge1.previous = n3;
          edge2.next = n1;
          edge2.previous = n2;
          edge3.next = n3;
          edge4.previous = n0;
          edge5.previous = n5;
          edge5.next = n4;

          if(border < m_edges.size()) {
               for(auto i = m_edges.begin(), end = m_edges.end(); i != end; ++i) {
                    Edge& temp = *i;
                    if(temp.next >= border) temp.next += 6;
                    if(temp.opposite >= border) temp.opposite += 6;
                    if(temp.previous >= border) temp.previous += 6;
               }
          }

          Edge edges[] = { Edge(v2, e4, n1, e0), Edge(vn, n2, n0, e2), Edge(v0, e2, n4, n1), Edge(vn, e1, n5, e3), Edge(vn, n5, n2, e5), Edge(v3, e5, n3, n4) };
          m_edges.insert(m_edges.begin() + border, edges, edges + 6);

          m_outgoing[v0] = n4;
          m_outgoing.push_back(n2);
     }

/**********************************************************************************
 * split triangle
 *
 *          BEFORE
 *       ___  v2
 *         //    \\
 *        //      \\
 *       //        \\
 *      //e2      e1\\
 *     //            \\ |
 *    //      e0      \\|
 * v0 ================== v1
 *
 *
 *                    AFTER
 *               ___    v1
 *                 //  /||   \\
 *                //  / ||    \\
 *               //     ||     \\
 *              //    g2||f1    \\
 *             //       || /     \\ |
 *            //        ||/       \\|
 *           //         vn         \\
 *          //e2       //\\       e1\\
 *         //        //    \\        \\
 *        //       //        \\       \\
 *       //      //            \\      \\
 *      //   f2//                \\g1   \\
 *     //    //g0                f0\\    \\
 *    //   //                        \\   \\
 *   //  //                            \\  \\
 *  // //                                \\ \\
 * ////___              e0                 \\\\
 * v0 ====================================== v1
 *
 **********************************************************************************/
     void splitTriangle(hpuint triangle, hpreal u = 1.0/3.0, hpreal v = 1.0/3.0) {
          auto border = 3 * m_nTriangles;
          auto e0 = 3 * triangle;
          auto e1 = e0 + 1;
          auto e2 = e1 + 1;
          auto e = std::begin(m_edges) + e0;
          auto& edge0 = *e;
          auto& edge1 = *(++e);
          auto& edge2 = *(++e);
          auto v0 = edge2.vertex;
          auto v1 = edge0.vertex;
          auto v2 = edge1.vertex;
          auto& vertex0 = getVertex(v0);
          auto& vertex1 = getVertex(v1);
          auto& vertex2 = getVertex(v2);
          auto vn = hpuint(m_vertices.size());

          hpuint f0 = m_edges.size(), g0 = f0 + 1, f1 = g0 + 1, g1 = f1 + 1, f2 = g1 + 1, g2 = f2 + 1;

          Vertex vertexn(vertex0);//TODO: improve; possibilities: vertex as parameter, VertexUtils::mix(v1,v2), lambda function
          vertexn.position = vertex0.position * u + vertex1.position * v + vertex2.position * (1.0f - u - v);
          m_vertices.push_back(vertexn);
          m_outgoing.push_back(g0);

          if(border < m_edges.size()) {//TODO: refactor into method insertEdges
               for(auto& edge : m_edges) {
                    if(edge.next >= border) edge.next += 6;
                    if(edge.opposite >= border) edge.opposite += 6;
                    if(edge.previous >= border) edge.previous += 6;
               }
          }

          edge0.next = f0;
          edge0.previous = g0;
          edge1.next = f1;
          edge1.previous = g1;
          edge2.next = f2;
          edge2.previous = g2;

          Edge edges[] = { Edge(vn, g0, g1, e0), Edge(v0, e0, f2, f0), Edge(vn, g1, g2, e1), Edge(v1, e1, f0, f1), Edge(vn, g2, g0, e2), Edge(v2, e2, f1, f2) };
          m_edges.insert(m_edges.begin() + border, edges, edges + 6);
     }

private:
     std::vector<Edge> m_edges;
     hpuint m_nTriangles;
     Indices m_outgoing;
     std::vector<Vertex> m_vertices;

};//TriangleGraph

namespace trg {

class SpokesWalker {
public:
     SpokesWalker(const std::vector<Edge>& edges, hpindex e)
          : m_e(e), m_edges(edges) {}

     auto operator==(const SpokesWalker& walker) const { return m_e == walker.m_e; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return m_e; }

     auto& operator++() {
          m_e = m_edges[m_edges[m_e].previous].opposite;
          return *this;
     }

     auto& operator--() {
          m_e = m_edges[m_edges[m_e].opposite].next;
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
     hpindex m_e;
     const std::vector<Edge>& m_edges;

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

     auto operator*() const { return hpindex(*m_i / 3); }

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

          auto e = *m_i;
          auto t = make_triangle_index(e);
          auto i = make_edge_offset(e);
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

}//namespace trg

template<class Picker>
Indices basic_cut(const std::vector<Edge>& edges, hpindex t, Picker&& pick) {
     auto neighbors = Indices(edges.size() << 1, std::numeric_limits<hpindex>::max());

     neighbors[6 * t + 0] = 6 * t + 2;
     neighbors[6 * t + 1] = 6 * t + 1;
     neighbors[6 * t + 2] = 6 * t + 0;
     neighbors[6 * t + 3] = 6 * t + 2;
     neighbors[6 * t + 4] = 6 * t + 1;
     neighbors[6 * t + 5] = 6 * t + 0;

     for(auto e = pick(neighbors); e != std::numeric_limits<hpindex>::max(); e = pick(neighbors)) {
          auto& edge0 = edges[e];
          auto& edge1 = edges[edge0.opposite];
          auto p = neighbors[e << 1];
          auto n = neighbors[(e << 1) + 1];

          assert(neighbors[e << 1] != std::numeric_limits<hpindex>::max());
          assert(neighbors[edge0.opposite << 1] == std::numeric_limits<hpindex>::max());
          assert(neighbors[edge1.next << 1] == std::numeric_limits<hpindex>::max());
          assert(neighbors[edge1.previous << 1] == std::numeric_limits<hpindex>::max());

          neighbors[(p << 1) + 1] = edge1.next;
          neighbors[edge1.next << 1] = p;
          neighbors[(edge1.next << 1) + 1] = edge1.previous;
          neighbors[edge1.previous << 1] = edge1.next;
          neighbors[(edge1.previous << 1) + 1] = n;
          neighbors[n << 1] = edge1.previous;
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

template<class RandomDev>
Indices cut(const std::vector<Edge>& edges, hpindex t, RandomDev&& random) {
     auto cache = boost::dynamic_bitset<>(edges.size(), false);
     auto range = std::mt19937();

     range.seed(random());
     auto vdist = std::uniform_int_distribution<std::mt19937::result_type>(0, (size(edges)/3) - 1);
     //const auto xt = (t != std::numeric_limits<hpindex>::max()) ? t : vdist(range);
     hpindex xt = 0;
     cache[3*xt + 0] = true;
     cache[3*xt + 1] = true;
     cache[3*xt + 2] = true;

     return basic_cut(edges, xt, [&](auto& neighbors) {
          //for(auto e : boost::irange(0u, hpindex(mesh.getEdges().size())))
          //     if(neighbors[e << 1] != std::numeric_limits<hpindex>::max() && neighbors[mesh.getEdge(e).opposite << 1] == std::numeric_limits<hpindex>::max()) return e;
          if(cache.count() == 0u) return std::numeric_limits<hpindex>::max();
          auto distribution = std::uniform_int_distribution<std::mt19937::result_type>(1, cache.count());
          auto n = distribution(range);
          auto e = cache.find_first();
          while(--n) e = cache.find_next(e);
          auto& edge = edges[edges[e].opposite];
          auto e0 = edges[edge.previous].opposite;
          auto e1 = edges[edge.next].opposite;
          auto b0 = neighbors[e0 << 1] == std::numeric_limits<hpindex>::max();
          auto b1 = neighbors[e1 << 1] == std::numeric_limits<hpindex>::max();
          if(b0) cache[edge.previous] = true;
          else cache[e0] = false;
          if(b1) cache[edge.next] = true;
          else cache[e1] = false;
          cache[e] = false;
          return hpindex(e);
     });
}


template<class Vertex, class RandomDev>
Indices cut(const TriangleGraph<Vertex>& graph, hpindex t, RandomDev&& random) { return cut(graph.getEdges(), t, random); }

template<class Vertex>
std::vector<Point2D> embed(const TriangleGraph<Vertex>& graph, const Indices& cut, const std::vector<Point2D>& polygon) {

     auto vertices = graph.getVertices();
     auto edges = graph.getEdges();
     
     auto nVertices = size(vertices);
     auto nEdges = size(edges);
     auto nTriangles = graph.getNumberOfTriangles();

     //Reorder Vertices 
     auto v = Indices();
     v.assign(nVertices, 0);
     hpuint insert_front = 0;
     hpuint insert_back = nVertices - 1;
     for(hpuint I = 0; I < nVertices; ++I){
          bool in_cut = false;
          for(auto e = std::begin(cut); e != std::end(cut); ++e){
               if(I == edges[*e].vertex){
                    in_cut = true;
                    break;
               }
          }
          if(in_cut){
               v[insert_back--] = I;
          } else {
               v[insert_front++] = I;
          }
          
     }
     hpuint n = insert_front;

     //Calculate all angles
     auto angles = std::vector<hpreal>();
     angles.reserve(3 * nTriangles);
     for(auto t = 0; t < nTriangles; ++t){
          auto v0 = graph.getVertex(t, 0);
          auto v1 = graph.getVertex(t, 1);
          auto v2 = graph.getVertex(t, 2);
          
          auto calc_angle = [&](Vertex x_i, Vertex a, Vertex b){
               return acos(glm::dot(glm::normalize(a.position - x_i.position), glm::normalize(b.position - x_i.position)));
          };
          
          angles[3 * t + 0] = calc_angle(v0, v2, v1);
          angles[3 * t + 1] = calc_angle(v1, v0, v2);
          angles[3 * t + 2] = calc_angle(v2, v1, v0);
     }
     
     auto make_sparse_matrix_ = [&](hpuint n, hpuint m, const std::vector<Eigen::Triplet<hpreal> >& triplets){
          auto matrix = Eigen::SparseMatrix<hpreal>(n, m);
          matrix.setFromTriplets(std::begin(triplets), std::end(triplets));
          matrix.makeCompressed();
          return matrix;
     };
     
     auto make_edge_index_ = [&](hpuint i, hpuint j){
          for(hpuint e = 0; e < nEdges; ++e){
               if( (edges[e].vertex == v[j]) && (edges[edges[e].opposite].vertex == v[i]) ){
                    return e;
               }
          }
          return 0u;
     };
     
     auto alpha = [&](hpuint i, hpuint j){
          auto e = make_edge_index_(i, j);
          static constexpr hpuint o[3] = { 1, 2, 0 };
          auto opposite = edges[e].opposite;
          auto u = opposite / 3;
          return angles[3 * u + o[opposite - 3 * u]];
     };
     
     auto beta = [&](hpuint j, hpuint i){
          auto e = make_edge_index_(i, j);
          return angles[e];
     };
     
     auto N = [&](hpuint i){
          auto neighbors = Indices();
          for(auto e : edges){
               if(edges[e.opposite].vertex == v[i]){
                    hpuint J = e.vertex;
                    hpuint j;
                    for(j = 0; j < nVertices; ++j){
                         if (v[j] == J){
                              break;
                         }
                    }
                    bool contained = false;
                    for(hpuint k : neighbors){
                         if(k == j){
                              contained = true;
                              break;
                         }
                    }
                    if(!contained){
                         neighbors.push_back(j);
                    }
               }
          }
          return neighbors;
     };
     std::vector<std::vector<hpindex>> _N;
     _N.reserve(nVertices);
     for (unsigned i = 0; i < nVertices; i++)
          _N.push_back(N(i));
     
     //Calculate values delta
     auto w = std::vector<Eigen::Triplet<hpreal>>();
     auto sum_w = std::vector<hpreal>();
     sum_w.assign(nVertices, hpreal(0));
     for(hpuint i = 0; i < nVertices; ++i){
          //std::cout << "loop i=" << i << "\n";
          for(hpuint j : _N[i]){
               //std::cout << "- neighbor j=" << j << "\n";
               hpreal r_ij = glm::distance(vertices[v[i]].position, vertices[v[j]].position);
               hpreal val = (tan(alpha(i, j) / 2.0) + tan(beta(j, i) / 2.0)) / r_ij;
               w.push_back(Eigen::Triplet<hpreal>(i, j, val));
               sum_w[i] += val;
          }
     }
     std::cout << "pre: W = make_sparse_matrix\n";
     auto W = make_sparse_matrix_(nVertices, nVertices, w);
     
     std::cout << "computing deltas\n";
     auto delta = std::vector<Eigen::Triplet<hpreal>>();
     for(hpuint i = 0; i < nVertices; ++i){
          //std::cout << "deltas: i=" << i << "\n";
          for(hpuint j : _N[i]){
               //std::cout << "- j=" << j << "\n";
               delta.push_back(Eigen::Triplet<hpreal>(i, j, W.coeff(i, j) / sum_w[i]));
          }
     }
     std::cout << "pre: Delta = make_sparse_matrix\n";
     auto Delta = make_sparse_matrix_(nVertices, nVertices, delta);

     std::unordered_set<hpindex> cut_cache(size(cut));
     for (auto e : cut)
          cut_cache.emplace(e);
     auto in_cut = [&cut_cache] (hpuint edge) { return cut_cache.count(edge) > 0; };
     //auto in_cut = [&](hpuint edge){
     //     for(auto e = std::begin(cut); e != std::end(cut); ++e){
     //          if(*e == edge){
     //               return true;
     //          }
     //     }
     //     return false;
     //};
     
     auto fixed = [&](hpuint i, hpuint j){
          for(hpuint e = 0; e < size(cut); ++e){
               if(edges[cut[e]].vertex == v[j]){
                    auto f = edges[cut[e]].next;
                    bool found = false;
                    while(!in_cut(f)){
                         if(edges[f].vertex == v[i]){
                              found = true;
                              break;
                         }
                         f = edges[edges[f].opposite].next;
                    }
                    if(found){
                         return polygon[(e+1) % size(polygon)];
                    }
               }
          }
          return Point2D(0, 0);
     };
     
     //fill in solution vectors
     auto u_bar = Eigen::Matrix<hpreal, Eigen::Dynamic, 1>(n);
     u_bar.setZero();
     auto v_bar = Eigen::Matrix<hpreal, Eigen::Dynamic, 1>(n);
     v_bar.setZero();
     
     std::cout << "pre: loop setting u_bar, v_bar\n";
     for(hpuint i = 0; i < n; ++i){
          for(hpuint j : _N[i]){
               if(j >= n){
                    auto point = fixed(i, j);
                    u_bar[i] += Delta.coeff(i, j) * (point).x;
                    v_bar[i] += Delta.coeff(i, j) * (point).y;
               }
          }
     }
     
     //fill in matrix
     std::cout << "pre: loop filling in matrix\n";
     auto a = std::vector<Eigen::Triplet<hpreal>>();
     for(hpuint i = 0; i < n; ++i){
          //std::cout << "- i=" << i << "\n";
          a.push_back(Eigen::Triplet<hpreal>(i, i, hpreal(1)));
          for(hpuint k : _N[i])
               if (k < n)
                    a.push_back(Eigen::Triplet<hpreal>(i, k, - Delta.coeff(i, k)));
     }
     
     std::cout << "pre: A = make_sparse_matrix\n";
     auto A = make_sparse_matrix_(n, n, a);
     
     Eigen::SparseLU<Eigen::SparseMatrix<hpreal>> solver;
     std::cout << "pre: analyze pattern\n";
     solver.analyzePattern(A);
     std::cout << "pre: factorize\n";
     solver.factorize(A);
     
     std::cout << "pre: solve for U\n";
     auto U = solver.solve(u_bar).eval();
     std::cout << "pre: solve for V\n";
     auto V = solver.solve(v_bar).eval();
     
     auto points = std::vector<Point2D>(nVertices);
     
     //std::cout << "pre: setting (inner) points\n";
     for(hpuint i = 0; i < n; ++i){
          points[v[i]].x = U[i];
          points[v[i]].y = V[i];
     }
     //std::cout << "pre: setting outer points\n";
     for(hpuint i = n; i < nVertices; ++i){
          auto j = std::begin(_N[i]);
          while(*j >= n){ ++j; }
          //points[v[i]] = fixed(*j, i);
     }
     
     return points;
}

template<class Vertex, class VertexRule, class EdgeRule>
TriangleMesh<Vertex> loopivide(const TriangleGraph<Vertex>& graph, VertexRule&& vertexRule, EdgeRule&& edgeRule) {
     auto vertices = std::vector<Vertex>();
     auto indices = Indices();
     auto es = Indices(graph.getNumberOfEdges(), std::numeric_limits<hpindex>::max());

     vertices.reserve(graph.getNumberOfVertices() + (graph.getNumberOfEdges() >> 1));
     indices.reserve((graph.getNumberOfTriangles() << 2) * 3);

     auto v = 0u;
     for(auto& vertex : graph.getVertices()) {
          auto ring = make_ring(graph, v++);
          vertices.emplace_back(vertexRule(vertex, std::begin(ring), std::end(ring)));
     }

     visit_diamonds(graph, [&](auto e, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
          auto& edge = graph.getEdge(e);
          es[e] = vertices.size();
          es[edge.opposite] = vertices.size();
          vertices.emplace_back(edgeRule(vertex0, vertex2, vertex1, vertex3));
     });

     auto e = std::begin(es) - 1;
     visit_triplets(graph.getIndices(), [&](auto v0, auto v2, auto v5) {
          auto v1 = *(++e);
          auto v4 = *(++e);
          auto v3 = *(++e);

          indices.insert(std::end(indices), {
               v0, v1, v3,
               v1, v2, v4,
               v1, v4, v3,
               v4, v5, v3
          });
     });

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

template<class Vertex>
TriangleMesh<Vertex> loopivide(const TriangleGraph<Vertex>& graph) {
     return loopivide(graph, [](auto& center, auto begin, auto end) {
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
std::tuple<Point3D, Point3D> make_axis_aligned_bounding_box(const TriangleGraph<Vertex>& graph) { return make_axis_aligned_bounding_box(graph.getVertices()); }

inline hpindex make_edge_index(const Edge& edge) { return 3 * make_triangle_index(edge) + make_edge_offset(edge); }
     
template<class Vertex>
boost::optional<hpindex> make_edge_index(const TriangleGraph<Vertex>& graph, hpindex v0, hpindex v1) {
     auto e = make_spokes_enumerator(graph.getEdges(), graph.getOutgoing(v0));
     do if(graph.getEdge(*e).vertex == v1) return *e; while(++e);
     return boost::none;
}

inline hpindex make_edge_offset(const Edge& edge) { return 3 - edge.next - edge.previous + 6 * make_triangle_index(edge); }

template<class Transformer>
auto make_fan(EnumeratorTransformer<trg::FanEnumerator, Transformer> e) {
     using T = decltype(*e);

     auto fan = std::vector<T>();
     do fan.push_back(*e); while(++e);
     return fan;
}

template<class Vertex>
Indices make_fan(const TriangleGraph<Vertex>& graph, hpuint v) { return make_fan(make_fan_enumerator(graph, v)); }

inline trg::FanEnumerator make_fan_enumerator(const std::vector<Edge>& edges, hpuint e) { return { { edges, e } }; }
     
template<class Vertex>
trg::FanEnumerator make_fan_enumerator(const TriangleGraph<Vertex>& graph, hpuint v) { return { { graph.getEdges(), graph.getOutgoing(v) } }; }

template<class Vertex>
Indices make_indices(const TriangleGraph<Vertex>& graph) {
     auto indices = Indices();
     const auto nTriangles = size(graph);

     indices.reserve(3 * nTriangles);
     visit_triplets(std::begin(graph.getEdges()), nTriangles, 3, [&] (const auto& edge0, const auto& edge1, const auto& edge2) {
          indices.emplace_back(edge2.vertex);
          indices.emplace_back(edge0.vertex);
          indices.emplace_back(edge1.vertex);
     });

     return indices;
}

inline hpindex make_neighbor_index(const std::vector<Edge>& edges, hpuint t, hpuint i) { return make_triangle_index(edges[3 * t + i].opposite); }

template<class Vertex>
hpuint make_neighbor_index(const TriangleGraph<Vertex>& graph, hpuint t, hpuint i) { return make_neighbor_index(graph.getEdges(), t, i); }

template<class Vertex>
hpuint make_neighbor_offset(const TriangleGraph<Vertex>& graph, hpuint t, hpuint u) { return make_neighbor_offset(graph.getEdges(), t, u); }

template<class Vertex>
Indices make_neighbors(const TriangleGraph<Vertex>& graph) { return make_neighbors(graph.getEdges(), size(graph)); }

template<class Transformer>
auto make_ring(EnumeratorTransformer<trg::RingEnumerator, Transformer> e) {
     using T = decltype(*e);

     auto ring = std::vector<T>();
     do ring.push_back(*e); while(++e);
     return ring;
}

inline trg::RingEnumerator make_ring_enumerator(const std::vector<Edge>& edges, hpuint e) { return { { edges, e } }; }

template<class Transformer>
EnumeratorTransformer<trg::RingEnumerator, Transformer> make_ring_enumerator(const std::vector<Edge>& edges, hpuint e, Transformer&& transform) { return { make_ring_enumerator(edges, e), std::forward<Transformer>(transform) }; }

template<class Vertex>
auto make_ring_enumerator(const TriangleGraph<Vertex>& graph, hpuint v) { return make_ring_enumerator(graph.getEdges(), graph.getOutgoing(v), [&](auto t, auto i) { return graph.getVertex(t, i); }); }

inline trg::SpokesEnumerator make_spokes_enumerator(const std::vector<Edge>& edges, hpuint e) { return { { edges, e } }; }

template<class Transformer>
EnumeratorTransformer<trg::SpokesEnumerator, Transformer> make_spokes_enumerator(const std::vector<Edge>& edges, hpuint e, Transformer&& transform) { return { make_spokes_enumerator(edges, e), std::forward<Transformer>(transform) }; }

template<class Vertex>
auto make_spokes_enumerator(const TriangleGraph<Vertex>& graph, hpuint v) { return make_spokes_enumerator(graph.getEdges(), graph.getOutgoing(v), [&](auto e) { return graph.getEdge(e); }); }

inline trg::SpokesWalker make_spokes_walker(const std::vector<Edge>& edges, hpindex e) { return { edges, e }; }

inline hpindex make_triangle_index(const Edge& edge) { return make_triangle_index(edge.next); }

template<class Vertex>
TriangleGraph<Vertex> make_triangle_graph(std::vector<Vertex> vertices, const Indices& indices) { return { std::move(vertices), make_edges(indices), hpuint(indices.size() / 3) }; }

template<class Vertex>
TriangleGraph<Vertex> make_triangle_graph(const TriangleMesh<Vertex>& mesh) { return make_triangle_graph(mesh.getVertices(), mesh.getIndices()); }

template<class Vertex>
hpuint size(const TriangleGraph<Vertex>& graph) { return graph.getNumberOfTriangles(); }

template<class Vertex>
Indices trim(const TriangleGraph<Vertex>& graph, Indices path) {
     auto i = std::begin(path);
     auto j = next(i);

     auto test = [](auto e) { return e != std::numeric_limits<hpindex>::max(); };
     auto next = [&](auto i) {
          auto j = std::find_if(i + 1, std::end(path), test);
          if(j == std::end(path)) j = std::find_if(std::begin(path), i, test);
          return j;
     };;
     auto previous = [&](auto i) {
          auto j = std::find_if(std::make_reverse_iterator(i), std::rend(path), test);
          if(j == std::rend(path)) return next(std::begin(path));
          else return j.base() - 1;
     };;
   
     while(true) {
          if(*j == graph.getEdge(*i).opposite) {
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

template<class Vertex>
hpuint validate_cut(const TriangleGraph<Vertex>& graph, const Indices& path) {
     auto cache = path;
     auto todo = std::stack<hpindex>();
     auto visited = boost::dynamic_bitset<>(graph.getNumberOfTriangles(), false);

     auto contains = [&](auto e) { return std::binary_search(std::begin(cache), std::end(cache), e); };
     auto push = [&](auto e) { if(!contains(e)) todo.push(make_triangle_index(graph.getEdge(e).opposite)); };

     if(path.size() == 0) return 0;
     if(!validate_path<true>(graph, path)) return 1;
     std::sort(std::begin(cache), std::end(cache));

     //Make sure each edge is exactly once on the path.
     if(std::unique(std::begin(cache), std::end(cache)) != std::end(cache)) return 2;

     //Make sure opposite edge is always on path.
     for(auto e : path) if(!contains(graph.getEdge(e).opposite)) return 3;

     //Make sure the path does not cross itself.
     for(auto i = std::begin(path), end = std::end(path) - 1; i != end; ++i) {
          if(graph.getEdge(i[0]).opposite == i[1]) continue;
          auto e = make_spokes_walker(graph.getEdges(), i[1]);
          while(!contains(*(++e)));
          if(graph.getEdge(*e).opposite != i[0]) return 4;
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
          auto e = make_spokes_enumerator(graph.getEdges(), e1);
          do if(*e == e1) return true; while(++e);
          return false;
     };

     if(path.size() == 0) return true;
     for(auto i = std::begin(path), end = std::end(path) - 1; i != end; ++i) if(!exists(i[0], i[1])) return false;
     if(loop && !exists(path.back(), path.front())) return false;

     return true;
}

template<class Visitor>
void visit_diamonds(const std::vector<Edge>& edges, Visitor&& visit) {
     auto visited = boost::dynamic_bitset<>(edges.size(), false);

     auto e = hpindex(-1);
     for(auto& edge : edges) {
          if(visited[++e]) continue;
          assert(!visited[edge.opposite]);
          visit(e, edge.vertex, edges[edge.next].vertex, edges[edge.previous].vertex, edges[edges[edge.opposite].next].vertex);
          visited[edge.opposite] = true;
     }
}

template<class Vertex, class Visitor>
void visit_diamonds(const TriangleGraph<Vertex>& graph, Visitor&& visit) { visit_diamonds(graph.getEdges(), [&](auto e, auto v0, auto v1, auto v2, auto v3) { visit(e, graph.getVertex(v0), graph.getVertex(v1), graph.getVertex(v2), graph.getVertex(v3)); }); }

template<class Visitor>
void visit_fan(trg::FanEnumerator e, Visitor&& visit) { do apply(visit, *e); while(++e); }

template<class Visitor>
void visit_ring(trg::RingEnumerator e, Visitor&& visit) { do apply(visit, *e); while(++e); }

template<class Visitor>
void visit_spokes(trg::SpokesEnumerator e, Visitor&& visit) { do apply(visit, *e); while(++e); }

template<class Visitor>
void visit_vertices(const std::vector<Edge>& edges, Visitor&& visit) {
     boost::dynamic_bitset<> visited(edges.size(), false);

     for(auto e : boost::irange(0lu, edges.size())) {
          if(visited[e]) continue;
          visit(e);
          visit_spokes(make_spokes_enumerator(edges, e), [&](auto e) { visited[e] = true; });
     }
}

}//namespace happah

