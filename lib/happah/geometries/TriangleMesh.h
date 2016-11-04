// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <cilk/cilk.h>
#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_map>

#include "happah/geometries/Geometry.h"
#include "happah/geometries/Mesh.h"
#include "happah/math/Space.h"
#include "happah/utils/DeindexedArray.h"
#include "happah/utils/visitors.h"
#include "happah/geometries/TriangleMeshUtils.h"

namespace happah { 

namespace mode { namespace TriangleMesh { //TODO: maybe rename to output? 

using namespace mode::Mesh;
constexpr hpuint TRIANGLES = 10;

}//namespace Mesh
}//namespace mode

namespace view { namespace TriangleMesh {

using namespace view::Mesh;
constexpr hpuint TRIANGLES = 10;

}//namespace Mesh
}//namespace view

enum class Format { DIRECTED_EDGE, SIMPLE };

namespace {

namespace Mode = happah::mode::TriangleMesh;
namespace View = happah::view::TriangleMesh;

}//namespace

template<class Vertex, Format t_format = Format::SIMPLE>
class TriangleMesh;

template<class T, class Visitor>
void visit_triangles(const std::vector<T>& ts, const std::vector<hpuint>& indices, Visitor&& visit) { visit_triplets(deindex(ts, indices).begin(), indices.size() / 3, 3, std::forward<Visitor>(visit)); }

//TODO: shortestpathfinder and weigher now use trianglemeshutils mode and view although only mesh required as input
template<class Vertex>
class TriangleMesh<Vertex, Format::SIMPLE> : public Geometry2D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;

public:
     using Indices = typename Mesh<Vertex>::Indices;
     using Vertices = typename Mesh<Vertex>::Vertices;

     TriangleMesh(Vertices vertices, Indices indices)
          : Geometry2D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)) {}

     template<Format format>
     TriangleMesh(const TriangleMesh<Vertex, format>& mesh)
          : TriangleMesh(mesh.getVertices(), mesh.getIndices()) {}

     template<Format format>
     TriangleMesh(TriangleMesh<Vertex, format>&& mesh)
          : TriangleMesh(std::move(mesh.getVertices()), std::move(mesh.getIndices())) {}

     ~TriangleMesh() {}

};//TriangleMesh
using TriangleMesh2D = TriangleMesh<VertexP2>;
using TriangleMesh2D_ptr = std::shared_ptr<TriangleMesh2D>;
using TriangleMesh3D = TriangleMesh<VertexP3N>;
using TriangleMesh3D_ptr = std::shared_ptr<TriangleMesh3D>;

struct Edge {
     hpuint next;
     hpuint opposite;
     hpuint previous;
     hpuint vertex;//vertex to which edge points

     Edge(hpuint vertex, hpuint next, hpuint opposite, hpuint previous)
          : next(next), opposite(opposite), previous(previous), vertex(vertex) {}

};

std::vector<Edge> make_edges(const std::vector<hpuint>& indices);

template<class Test>
boost::optional<hpuint> find_if_in_spokes(const std::vector<Edge>& edges, hpuint begin, Test&& test) {
     auto e = begin;
     do {
          auto& edge = edges[e];
          if(test(edge)) return e;
          e = edges[edges[edge.next].next].opposite;
     } while(e != begin);
     return boost::none;
}

boost::optional<hpuint> find_in_ring(const std::vector<Edge>& edges, hpuint begin, hpuint v);

template<class Visitor, bool closed = false>
void visit_fan(const std::vector<hpuint>& neighbors, hpuint t, hpuint i, Visitor&& visit) {
     auto current = t;

     if(!closed) do {
          hpuint previous;
          visit_triplet(neighbors, current, [&](hpuint n0, hpuint n1, hpuint n2) { previous = (i == 0) ? n0 : (i == 1) ? n1 : n2; });
          if(previous == UNULL) break;
          visit_triplet(neighbors, previous, [&](hpuint n0, hpuint n1, hpuint n2) { i = (n0 == current) ? 1 : (n1 == current) ? 2 : 0; });
          current = previous;
     } while(current != t);
     t = current;

     do {
          hpuint next;
          visit(current);
          visit_triplet(neighbors, current, [&](hpuint n0, hpuint n1, hpuint n2) { next = (i == 0) ? n2 : (i == 1) ? n0 : n1; });
          if(!closed && next == UNULL) break;
          visit_triplet(neighbors, next, [&](hpuint n0, hpuint n1, hpuint n2) { i = (n0 == current) ? 0 : (n1 == current) ? 1 : 2; });
          current = next;
     } while(current != t);
}

template<class Visitor>
void visit_fans(const std::vector<hpuint>& neighbors, Visitor&& visit) {
     boost::dynamic_bitset<> visited(neighbors.size(), false);

     auto update_visited = [&](hpuint current, hpuint next, hpuint n0, hpuint n1, hpuint n2) {
          if(n0 == next) visited[3 * current + 1] = true;
          else if(n1 == next) visited[3 * current + 2] = true;
          else {
               assert(n2 == next);
               visited[3 * current] = true;
          }
     };

     auto do_visit_fans = [&](hpuint t, hpuint i) {
          std::vector<hpuint> fan;
          visit_fan(neighbors, t, i, [&](hpuint n) { fan.push_back(n); });
          visit(fan);
          if(fan.size() == 1) return;
          visit_pairs(fan.begin(), fan.size() - 1, 1, [&](hpuint current, hpuint next) {
               visit_triplet(neighbors, current, [&](hpuint n0, hpuint n1, hpuint n2) { update_visited(current, next, n0, n1, n2); });
          });
          auto current = fan.back();
          auto next = fan.front();
          visit_triplet(neighbors, current, [&](hpuint n0, hpuint n1, hpuint n2) { if(n0 == next || n1 == next || n2 == next) update_visited(current, next, n0, n1, n2); });
     };

     for(auto t = 0lu, end = neighbors.size() / 3; t != end; ++t) {
          if(!visited[3 * t]) do_visit_fans(t, 0);
          if(!visited[3 * t + 1]) do_visit_fans(t, 1);
          if(!visited[3 * t + 2]) do_visit_fans(t, 2);
     }
}

template<class Visitor>
void visit_spokes(const std::vector<Edge>& edges, hpuint begin, Visitor&& visit) {
     auto e = begin;
     do {
          auto& edge = edges[e];
          visit(edge);
          e = edges[edges[edge.next].next].opposite;
     } while(e != begin);
}

template<class Visitor>
void visit_fans(const std::vector<Edge>& edges, Visitor&& visit) {
     boost::dynamic_bitset<> visited(edges.size(), false);
     for(auto begin = 0lu, end = edges.size(); begin < end; ++begin) {
          if(visited[begin]) continue;
          std::vector<hpuint> triangles;
          visit_spokes(edges, begin, [&](const Edge& edge) {
               triangles.emplace_back(edge.next / 3);
               visited[edge.previous] = true;
          });
          visit(triangles.begin(), triangles.end());
     }
}

template<class Visitor>
void visit_fans(const std::vector<Edge>& edges, const std::vector<hpuint>& outgoing, Visitor&& visit) {
     for(auto begin : outgoing) {
          std::vector<hpuint> triangles;
          visit_spokes(edges, begin, [&](const Edge& edge) { triangles.emplace_back(edge.next / 3); });
          visit(triangles.begin(), triangles.end());
     }
}

template<class Visitor>
void visit_ring(const std::vector<Edge>& edges, hpuint begin, Visitor&& visit) { visit_spokes(edges, begin, [&](const Edge& edge) { visit(edge.vertex); }); }

template<class Visitor>
void visit_rings(const std::vector<Edge>& edges, Visitor&& visit) {
     boost::dynamic_bitset<> visited(edges.size(), false);
     for(auto begin = 0lu, end = edges.size(); begin < end; ++begin) {
          if(visited[begin]) continue;
          std::vector<hpuint> vertices;
          visit_spokes(edges, begin, [&](const Edge& edge) {
               vertices.emplace_back(edge.vertex);
               visited[edge.previous] = true;
          });
          visit(vertices.begin(), vertices.end());
     }
}

template<class Visitor>
void visit_rings(const std::vector<Edge>& edges, const std::vector<hpuint>& outgoing, Visitor&& visit) {
     for(auto begin : outgoing) {
          std::vector<hpuint> vertices;
          visit_spokes(edges, begin, [&](const Edge& edge) { vertices.emplace_back(edge.vertex); });
          visit(vertices.begin(), vertices.end());
     }
}

template<class Visitor>
void visit_thorns(const std::vector<Edge>& edges, hpuint t, Visitor&& visit) {
     auto e = edges.begin() + 3 * t;
     auto n0 = (*e).opposite / 3;
     auto n1 = (*(++e)).opposite / 3;
     auto n2 = (*(++e)).opposite / 3;
     visit(n0, n1, n2);
}

template<class Vertex>
class TriangleMesh<Vertex, Format::DIRECTED_EDGE> : public Geometry2D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;

public:
     using Indices = typename Mesh<Vertex>::Indices;
     using Vertices = typename Mesh<Vertex>::Vertices;

private:
     template<int t_dummy, hpuint t_view, hpuint... t_modes>//TODO: remove
     class Iterator {
     public:
          using difference_type = hpuint;
          using value_type = std::tuple<typename Iterator<t_dummy, t_view, t_modes>::value_type...>;

          Iterator(const TriangleMesh& triangles, hpuint triangle)
               : m_i(Iterator<t_dummy, t_view, t_modes>(triangles, triangle)...) {}

          difference_type operator-(const Iterator& iterator) { return std::get<0>(m_i) - std::get<0>(iterator.m_i); }

          Iterator& operator++() {
               preincrement();
               return *this;
          }

          Iterator& operator--() {
               predecrement();
               return *this;
          }

          Iterator operator++(int) {
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          Iterator operator--(int) {
               Iterator iterator(*this);
               --(*this);
               return iterator;
          }

          bool operator!=(const Iterator& iterator) const { return std::get<0>(iterator.m_i) != std::get<0>(m_i); }

          value_type operator*() const { return getValue(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

     private:
          std::tuple<Iterator<t_dummy, t_view, t_modes>...> m_i;

          template<unsigned long... Is>
          value_type getValue(std::index_sequence<Is...>) const { return std::make_tuple(*(std::get<Is>(m_i))...); }

          BUILD_TUPLE_HANDLER_METHODS(predecrement, doPredecrement)

          void predecrement() { predecrement(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<hpuint mode>
          void doPredecrement(Iterator<t_dummy, t_view, mode>& i) { --i; }

          BUILD_TUPLE_HANDLER_METHODS(preincrement, doPreincrement)

          void preincrement() { preincrement(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<hpuint mode>
          void doPreincrement(Iterator<t_dummy, t_view, mode>& i) { ++i; }

     };//Iterator

     //NOTE: Iterator iterates over vertex' outgoing edges counter-clockwise if directed edges are arranged counter-clockwise.
     template<int t_dummy>
     class Iterator<t_dummy, View::VERTEX, Mode::EDGES> {
     public:
          using value_type = std::pair<const Edge&, hpuint>;

          Iterator(const TriangleMesh& triangleMesh, hpuint vertex)
               : m_edges(triangleMesh.m_edges), m_edge(triangleMesh.m_outgoing[vertex]) {}
         
          Iterator& operator++() {
               m_edge = m_edges[m_edges[m_edge].previous].opposite;
               return *this;
          }

          Iterator& operator--() {
               m_edge = m_edges[m_edges[m_edge].opposite].next;
               return *this;
          }

          Iterator operator++(int) { 
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          Iterator operator--(int) {
               Iterator iterator(*this);
               --(*this);
               return iterator;
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_edge != m_edge; }

          bool operator==(const Iterator& iterator) const { return iterator.m_edge == m_edge; }

          value_type operator*() const { return std::make_pair(std::cref(m_edges[m_edge]), m_edge); }
          //TODO: don't like that need to return m_edge

          Iterator& operator=(const Iterator& iterator) {
               m_edge = iterator.m_edge;
               return *this;
          }

     private:
          const std::vector<Edge>& m_edges;
          hpuint m_edge;

     };//Iterator

     //NOTE: Iterator iterates over vertex' neighborhood counter-clockwise if directed edges are arranged counter-clockwise.
     template<int t_dummy>
     class Iterator<t_dummy, View::VERTEX, Mode::VERTICES> {
     public:
          using value_type = hpuint;

          Iterator(const TriangleMesh& triangleMesh, hpuint vertex)
               : m_edges(triangleMesh.m_edges), m_i(m_edges.cbegin() + triangleMesh.m_outgoing[vertex]) {}

          Iterator& operator++() {
               m_i = m_edges.cbegin() + m_edges[(*m_i).previous].opposite;
               return *this;
          }

          Iterator& operator--() { 
               m_i = m_edges.cbegin() + m_edges[(*m_i).opposite].next;
               return *this;
          }

          Iterator operator++(int) {
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          Iterator operator--(int) {
               Iterator iterator(*this);
               --(*this);
               return iterator;
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_i != m_i; }

          bool operator==(const Iterator& iterator) const { return iterator.m_i == m_i; }

          value_type operator*() const { return (*m_i).vertex; }

     private:
          const std::vector<Edge>& m_edges;
          typename std::vector<Edge>::const_iterator m_i;

     };//Iterator

public:
     //NOTE: Indices all have to be arranged counterclockwise.
     TriangleMesh(Vertices vertices, Indices indices)
          : Geometry2D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)), m_edges(make_edges(this->m_indices)), m_outgoing(this->m_vertices.size(), -1) { std::for_each(m_edges.begin(), m_edges.begin() + this->m_indices.size(), [&](const Edge& edge) { m_outgoing[edge.vertex] = edge.opposite; }); }

     template<Format format>
     TriangleMesh(const TriangleMesh<Vertex, format>& mesh)
          : TriangleMesh(mesh.getVertices(), mesh.getIndices()) {}

     template<Format format>
     TriangleMesh(TriangleMesh<Vertex, format>&& mesh)
          : TriangleMesh(std::move(mesh.getVertices()), std::move(mesh.getIndices())) {}

     virtual ~TriangleMesh() {}

     template<hpuint view, hpuint... modes>
     Iterator<0, view, modes...> cbegin(hpuint index) const { return Iterator<0, view, modes...>(*this, index); }//TODO: remove

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

     const Edge& getEdge(hpuint e) const { return m_edges[e]; }

     boost::optional<hpuint> getEdgeIndex(hpuint v0, hpuint v1) const { return find_in_ring(m_edges, m_outgoing[v0], v1); }

     const std::vector<Edge>& getEdges() const { return m_edges; }

     hpuint getNumberOfEdges() const { return m_edges.size(); }

     hpuint getNumberOfTriangles() const { return this->m_indices.size() / 3; }

     const std::vector<hpuint>& getOutgoing() const { return m_outgoing; }

     hpuint getOutgoing(hpuint v) const { return m_outgoing[v]; }

     std::tuple<const Vertex&, const Vertex&, const Vertex&> getTriangle(hpuint t) const { return std::tie(this->m_vertices[3 * t], this->m_vertices[3 * t + 1], this->m_vertices[3 * t + 2]); }

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
          auto border = this->m_indices.size();
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
          auto& vertex0 = this->getVertex(v0);
          auto v1 = edge1.vertex;
          auto& vertex1 = this->getVertex(v1);
          auto v2 = edge2.vertex;
          auto v3 = edge3.vertex;
          hpuint vn = this->m_vertices.size();
          Vertex vertexn(vertex0);//TODO: improve; possibilities: vertex as parameter, VertexUtils::mix(v1,v2)
          vertexn.position = vertex0.position * u + vertex1.position * (1.0f - u);
          this->m_vertices.push_back(vertexn);
          
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
          
          hpuint found = 0;
          for(auto i = this->m_indices.begin(), end = this->m_indices.end(); i != end; ++i) {//TODO: replace with simd
               auto j = i;
               hpuint u0 = *i;
               hpuint u1 = *(++i);
               hpuint u2 = *(++i);
               if((u0 == v0 && u1 == v2 && u2 == v1) || (u0 == v2 && u1 == v1 && u2 == v0) || (u0 == v1 && u1 == v0 && u2 == v2)) {
                    *j = vn;
                    *(++j) = v2;
                    *(++j) = v1;
                    if((++found) == 2) break;
               } else if((u0 == v0 && u1 == v1 && u2 == v3) || (u0 == v1 && u1 == v3 && u2 == v0) || (u0 == v3 && u1 == v0 && u2 == v1)) {
                    *j = vn;
                    *(++j) = v1;
                    *(++j) = v3;
                    if((++found) == 2) break;
               }
          }
          hpuint indices[] = { vn, v0, v2, vn, v3, v0 };
          this->m_indices.insert(this->m_indices.end(), indices, indices + 6);

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
          auto border = this->m_indices.size();
          auto offset = 3 * triangle;
          auto i = this->m_indices.cbegin() + offset;
          auto v0 = *i;
          auto& vertex0 = this->getVertex(v0);
          auto v1 = *(++i);
          auto& vertex1 = this->getVertex(v1);
          auto v2 = *(++i);
          auto& vertex2 = this->getVertex(v2);
          hpuint vn = this->m_vertices.size();

          auto e0 = *getEdgeIndex(v0, v1);
          auto& edge0 = m_edges[e0];
          auto e1 = edge0.next;
          auto& edge1 = m_edges[e1];
          auto e2 = edge1.next;
          auto& edge2 = m_edges[e2];
          hpuint f0 = m_edges.size(), g0 = f0 + 1, f1 = g0 + 1, g1 = f1 + 1, f2 = g1 + 1, g2 = f2 + 1;

          Vertex vertexn(vertex0);//TODO: improve; possibilities: vertex as parameter, VertexUtils::mix(v1,v2)
          vertexn.position = vertex0.position * u + vertex1.position * v + vertex2.position * (1.0f - u - v);
          this->m_vertices.push_back(vertexn);
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

          this->m_indices[offset + 2] = vn;
          hpuint indices[] = { vn, v1, v2, vn, v2, v0 };
          this->m_indices.insert(this->m_indices.end(), indices, indices + 6);

          //TODO: if triangle or edges are split and there are caches, caches must be updated
     }

private:
     std::vector<Edge> m_edges;
     std::vector<hpuint> m_outgoing;

};//TriangleMesh

template<class Mesh, class Space = typename Mesh::SPACE, class Vertex = typename Mesh::VERTEX, typename = void>
struct is_triangle_mesh : public std::false_type {};

template<class Mesh, class Space, class Vertex>
struct is_triangle_mesh<Mesh, Space, Vertex, typename std::enable_if<std::is_base_of<TriangleMesh<Vertex>, Mesh>::value && std::is_base_of<typename Mesh::SPACE, Space>::value>::type> : public std::true_type {};

template<class Vertex, Format t_format = Format::SIMPLE>
static TriangleMesh<Vertex, t_format> make_triangle_mesh(std::vector<Vertex> vertices, std::vector<hpuint> indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex, class Visitor>
void visit_spokes(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint begin, Visitor&& visit) { visit_spokes(mesh.getEdges(), begin, std::forward<Visitor>(visit)); }

template<class Vertex, class Visitor>
void visit_fans(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit) {
     for(auto begin : mesh.getOutgoing()) {
          std::vector<Vertex> triangles;
          visit_spokes(mesh.getEdges(), begin, [&](const Edge& edge) {
               auto triangle = mesh.getTriangle(edge.next / 3);
               triangles.emplace_back(std::get<0>(triangle));
               triangles.emplace_back(std::get<1>(triangle));
               triangles.emplace_back(std::get<2>(triangle));
          });
          visit(triangles.begin(), triangles.end());
     }
}

template<class Vertex, class Visitor>
void visit_ring(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v, Visitor&& visit) { visit_spokes(mesh.getEdges(), mesh.getOutgoing(v), [&](const Edge& edge) { visit(mesh.getVertex(edge.vertex)); }); }

template<class Vertex, class Visitor>
void visit_rings(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit) {
     for(auto begin : mesh.getOutgoing()) {
          std::vector<Vertex> vertices;
          visit_spokes(mesh.getEdges(), begin, [&](const Edge& edge) { vertices.emplace_back(mesh.getVertex(edge.vertex)); });
          visit(vertices.begin(), vertices.end());
     }
}

template<class Vertex>
hpuint get_degree(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v) {
     auto degree = 0u;
     visit_spokes(mesh.getEdges(), mesh.getOutgoing(v), [&](const Edge& edge) { ++degree; });
     return degree;
}

}//namespace happah

