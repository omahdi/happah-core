// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <string>

#include "happah/geometries/Geometry.h"
#include "happah/geometries/Mesh.h"
#include "happah/math/Space.h"
#include "happah/readers/ReaderHPH.h"
#include "happah/utils/DeindexedArray.h"
#include "happah/utils/visitors.h"
#include "happah/geometries/TriangleMeshUtils.h"

namespace happah {

//DECLARATIONS

enum class Format { DIRECTED_EDGE, SIMPLE };

template<class Vertex, Format t_format = Format::SIMPLE>
class TriangleMesh;

struct Edge;

template<Format t_format>
class FanEnumerator;

template<Format t_format>
class SpokesEnumerator;

namespace trm {

template<Format t_format>
class RingEnumerator;

}//namespace trm

template<Format t_format>
class VerticesEnumerator;

template<class Test>
boost::optional<std::tuple<hpuint, hpuint, FanEnumerator<Format::SIMPLE> > > find_fan(const Indices& neighbors, Test&& test);

boost::optional<hpuint> find_in_ring(const std::vector<Edge>& edges, hpuint e, hpuint v);//TODO: interface not clear

bool is_neighbor(const Indices& neighbors, hpuint t, hpuint u);

//Return the index of this edge in the edges array.
hpuint make_edge_index(const Edge& edge);

//Return the offset of this edge among the three edges of its adjacent triangle.
hpuint make_edge_offset(const Edge& edge);

std::vector<Edge> make_edges(const Indices& indices);

template<Format format>
Indices make_fan(FanEnumerator<format> e);

Indices make_fan(const Indices& neighbors, hpuint t, hpuint i);

Indices make_fan(const std::vector<Edge>& edges, hpuint nTriangles, hpuint t, hpuint i);

template<class Vertex>
Indices make_fan(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint i);

FanEnumerator<Format::SIMPLE> make_fan_enumerator(const Indices& neighbors, hpuint t, hpuint i);

//Return the index of the ith neighbor of the tth triangle.
hpuint make_neighbor_index(const Indices& neighbors, hpuint t, hpuint i);

//Return the index of the ith neighbor of the tth triangle.
hpuint make_neighbor_index(const std::vector<Edge>& edges, hpuint t, hpuint i);

//Return the index of the ith neighbor of the tth triangle.
template<class Vertex>
hpuint make_neighbor_index(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint i);

//Assuming the tth and uth triangles are neighbors, return the offset of the uth triangle among the three neighbors of the tth triangle.
hpuint make_neighbor_offset(const Indices& neighbors, hpuint t, hpuint u);

//Assuming the tth and uth triangles are neighbors, return the offset of the uth triangle among the three neighbors of the tth triangle.
hpuint make_neighbor_offset(const std::vector<Edge>& edges, hpuint t, hpuint u);

//Assuming the tth and uth triangles are neighbors, return the offset of the uth triangle among the three neighbors of the tth triangle.
template<class Vertex>
hpuint make_neighbor_offset(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint u);

std::vector<hpuint> make_neighbors(const Indices& indices);

std::vector<hpuint> make_neighbors(const std::vector<Edge>& edges, hpuint nTriangles);

trm::RingEnumerator<Format::SIMPLE> make_ring_enumerator(const Indices& neighbors, hpuint t, hpuint i);

hpuint make_triangle_index(const Edge& edge);

template<class Vertex, Format format = Format::SIMPLE>
TriangleMesh<Vertex, format> make_triangle_mesh(std::vector<Vertex> vertices, Indices indices);

template<class Vertex = VertexP3, Format format = Format::SIMPLE>
TriangleMesh<Vertex, format> make_triangle_mesh(const std::string& path);

hpuint make_valence(const Indices& neighbors, hpuint t, hpuint i);

template<class Vertex>
hpuint make_valence(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v);

VerticesEnumerator<Format::SIMPLE> make_vertices_enumerator(const Indices& neighbors);

template<class Vertex, Format format>
hpuint size(const TriangleMesh<Vertex, format>& mesh);

template<class Visitor>
void visit_diamonds(const std::vector<Edge>& edges, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_diamonds(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit);

template<class Visitor>
void visit_edges(const Indices& neighbors, Visitor&& visit);

template<Format format, class Visitor>
void visit_fan(FanEnumerator<format> e, Visitor&& visit);

//Visit the fan about the ith vertex of the tth triangle.
template<class Visitor>
void visit_fan(const Indices& neighbors, hpuint t, hpuint i, Visitor&& visit);

template<class Visitor>
void visit_fan(const std::vector<Edge>& edges, hpuint nTriangles, hpuint t, hpuint i, Visitor&& visit);

template<class Visitor>
void visit_fans(const std::vector<Edge>& edges, hpuint nTriangles, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_fans(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit);

template<class Visitor>
void visit_ring(const std::vector<Edge>& edges, hpuint nTriangles, hpuint e, Visitor&& visit);

template<class Visitor>
void visit_ring(const std::vector<Edge>& edges, hpuint nTriangles, hpuint t, hpuint i, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_ring(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint i, Visitor&& visit);

template<class Visitor>
void visit_rings(const Indices& neighbors, Visitor&& visit);

template<class Visitor>
void visit_rings(const std::vector<Edge>& edges, hpuint nTriangles, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_rings(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit);

template<class Visitor>
void visit_spokes(const std::vector<Edge>& edges, hpuint nTriangles, hpuint e, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_spokes(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint e, Visitor&& visit);

template<class Visitor, bool closed = false>
void visit_subfan(const Indices& neighbors, hpuint t, hpuint i, hpuint u, Visitor&& visit);

template<class Visitor>
void visit_thorns(const std::vector<Edge>& edges, hpuint t, Visitor&& visit);

template<class Visitor>
void visit_vertices(const Indices& neighbors, Visitor&& visit);

template<class Visitor>
void visit_vertices(const std::vector<Edge>& edges, hpuint nTriangles, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_vertices(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit);

//DEFINITIONS

template<class Vertex>
class TriangleMesh<Vertex, Format::SIMPLE> : public Geometry2D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;
     using Vertices = typename Mesh<Vertex>::Vertices;

public:
     TriangleMesh() {}

     TriangleMesh(Vertices vertices, Indices indices)
          : Geometry2D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)) {}

     template<Format format>
     TriangleMesh(const TriangleMesh<Vertex, format>& mesh)
          : TriangleMesh(mesh.getVertices(), mesh.getIndices()) {}

     template<Format format>
     TriangleMesh(TriangleMesh<Vertex, format>&& mesh)
          : TriangleMesh(std::move(mesh.getVertices()), std::move(mesh.getIndices())) {}

private:
     template<class Stream>
     friend Stream& operator<<(Stream& stream, const TriangleMesh<Vertex, Format::SIMPLE>& mesh) {
          stream << mesh.m_vertices << '\n';
          stream << mesh.m_indices;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, TriangleMesh<Vertex, Format::SIMPLE>& mesh) {
          stream >> mesh.m_vertices;
          stream >> mesh.m_indices;
          return stream;
     }

};//TriangleMesh
using TriangleMesh2D = TriangleMesh<VertexP2>;
using TriangleMesh3D = TriangleMesh<VertexP3N>;

struct Edge {
     hpuint next;
     hpuint opposite;
     hpuint previous;
     hpuint vertex;//vertex to which edge points

     Edge(hpuint vertex, hpuint next, hpuint opposite, hpuint previous)
          : next(next), opposite(opposite), previous(previous), vertex(vertex) {}

};

template<class Vertex>
class TriangleMesh<Vertex, Format::DIRECTED_EDGE> : public Geometry2D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;
     using Vertices = typename Mesh<Vertex>::Vertices;

public:
     TriangleMesh() {}

     //NOTE: Indices all have to be arranged counterclockwise.
     TriangleMesh(Vertices vertices, Indices indices)
          : Geometry2D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)), m_edges(make_edges(this->m_indices)), m_outgoing(this->m_vertices.size(), UNULL) { std::for_each(std::begin(m_edges), std::begin(m_edges) + this->m_indices.size(), [&](auto& edge) { m_outgoing[edge.vertex] = edge.opposite; }); }

     template<Format format>
     TriangleMesh(const TriangleMesh<Vertex, format>& mesh)
          : TriangleMesh(mesh.getVertices(), mesh.getIndices()) {}

     template<Format format>
     TriangleMesh(TriangleMesh<Vertex, format>&& mesh)
          : TriangleMesh(std::move(mesh.getVertices()), std::move(mesh.getIndices())) {}

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

     auto getEdgeIndex(hpuint v0, hpuint v1) const { return find_in_ring(m_edges, m_outgoing[v0], v1); }

     auto& getEdges() const { return m_edges; }

     auto getNumberOfEdges() const { return m_edges.size(); }

     auto getNumberOfTriangles() const { return this->m_indices.size() / 3; }

     auto& getOutgoing() const { return m_outgoing; }

     auto getOutgoing(hpuint v) const { return m_outgoing[v]; }

     std::tuple<const Vertex&, const Vertex&, const Vertex&> getTriangle(hpuint t) const { return std::tie(getVertex(t, 0), getVertex(t, 1), getVertex(t, 2)); }

     using Model<Vertex>::getVertex;

     auto& getVertex(hpuint t, hpuint i) const { return this->m_vertices[this->m_indices[3 * t + i]]; }

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

          Vertex vertexn(vertex0);//TODO: improve; possibilities: vertex as parameter, VertexUtils::mix(v1,v2), lambda function
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
     }

private:
     std::vector<Edge> m_edges;
     Indices m_outgoing;

};//TriangleMesh

template<>
class FanEnumerator<Format::SIMPLE> {
public:
     FanEnumerator(const Indices& neighbors, hpuint t, hpuint i)
          : m_current(t), m_flag(false), m_i(i), m_j(i), m_neighbors(neighbors), m_t(t) {
          //if(!closed) do {
          do {
               hpuint previous;
               visit_triplet(m_neighbors, m_current, [&](hpuint n0, hpuint n1, hpuint n2) { previous = (m_i == 0) ? n0 : (m_i == 1) ? n1 : n2; });
               if(previous == std::numeric_limits<hpuint>::max()) break;
               visit_triplet(m_neighbors, previous, [&](hpuint n0, hpuint n1, hpuint n2) { m_i = (n0 == m_current) ? 1 : (n1 == m_current) ? 2 : 0; });
               m_current = previous;
          } while(m_current != m_t);
          m_t = m_current;
          m_j = m_i;
     }

     std::tuple<hpuint, hpuint> front() const { return std::make_tuple(m_t, m_i); }

     explicit operator bool() const { return m_current != UNULL && !(m_current == m_t && m_flag); }

     std::tuple<hpuint, hpuint> operator*() const { return std::make_tuple(m_current, m_j); }

     auto& operator++() {
          m_flag = true;
          auto next = std::numeric_limits<hpuint>::max();
          if(m_current != std::numeric_limits<hpuint>::max()) visit_triplet(m_neighbors, m_current, [&](hpuint n0, hpuint n1, hpuint n2) { next = (m_j == 0) ? n2 : (m_j == 1) ? n0 : n1; });
          if(next != std::numeric_limits<hpuint>::max()) visit_triplet(m_neighbors, next, [&](hpuint n0, hpuint n1, hpuint n2) { m_j = (n0 == m_current) ? 0 : (n1 == m_current) ? 1 : 2; });
          m_current = next;
          return *this;
     }

     auto& operator+=(hpuint n) {
          auto& me = *this;
          repeat(n, [&]() { ++me; });
          return me;
     }

     auto operator+(hpuint n) const {
          auto copy = *this;
          return copy += n;
     }

private:
     hpuint m_current;
     bool m_flag;
     hpuint m_i;
     hpuint m_j;
     const Indices& m_neighbors;
     hpuint m_t;

};//FanEnumerator<Format::SIMPLE>

template<>
class SpokesEnumerator<Format::SIMPLE> {
public:
     SpokesEnumerator(const Indices& neighbors, hpuint p, hpuint i)
          : m_e(neighbors, p, i), m_valid(true) {}

     operator bool() const { return m_valid; }

     auto operator*() const {
          static constexpr hpindex o[3] = { 2, 0, 1 };
          if(m_e) return *m_e;
          else return std::make_tuple(m_p, o[m_i]);
     }

     auto& operator++() {
          std::tie(m_p, m_i) = *m_e;
          m_valid = bool(m_e);
          ++m_e;
          m_valid &= bool(m_e) || (m_valid && std::get<0>(*m_e) == UNULL);
          return *this;
     }

private:
     FanEnumerator<Format::SIMPLE> m_e;
     hpuint m_i;
     hpuint m_p;
     bool m_valid;

};//SpokesEnumerator<Format::SIMPLE>

namespace trm {

template<>
class RingEnumerator<Format::SIMPLE> {
public:
     RingEnumerator(const Indices& neighbors, hpuint t, hpuint i)
          : m_e(neighbors, t, i), m_flag(true) {}//TODO: use spokes enumerator

     explicit operator bool() const { return m_flag; }

     std::tuple<hpuint, hpuint> operator*() const {
          auto t = 0u, i = 0u;
          std::tie(t, i) = *m_e;
          if(m_e) return std::make_tuple(t, ((i + 1) == 3) ? 0 : i + 1);
          else return std::make_tuple(m_t, (m_i == 0) ? 2 : (m_i == 1) ? 0 : 1);
     }

     auto& operator++() {
          std::tie(m_t, m_i) = *m_e;
          m_flag = !!m_e;
          ++m_e;
          m_flag &= !!m_e || (m_flag && !m_e && std::get<0>(*m_e) == UNULL);
          return *this;
     }

private:
     FanEnumerator<Format::SIMPLE> m_e;
     bool m_flag;
     hpuint m_i;
     hpuint m_t;

};//RingEnumerator<Format::SIMPLE>

}//namespace trm

template<>
class VerticesEnumerator<Format::SIMPLE> {
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
          visit_fan(m_neighbors, t, i, [&](auto u, auto j) { m_visited[3 * u + j] = true; });
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

template<class Test>
boost::optional<std::tuple<hpuint, hpuint, FanEnumerator<Format::SIMPLE> > > find_fan(const Indices& neighbors, Test&& test) {
     auto e = make_vertices_enumerator(neighbors);
     while(e) {
          auto t = 0u, i = 0u;
          std::tie(t, i) = *e;
          auto fan = make_fan_enumerator(neighbors, t, i);
          if(test(t, i, fan)) return std::make_tuple(t, i, fan);
          ++e;
     }
     return boost::none;
}

template<Format format>
Indices make_fan(FanEnumerator<format> e) {
     auto fan = Indices();
     while(e) {
          auto t = 0u, i = 0u;
          std::tie(t, i) = *e;
          fan.push_back(t);
          fan.push_back(i);
          ++e;
     }
     return fan;
}

template<class Vertex>
Indices make_fan(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint i) { return make_fan(mesh.getEdges(), mesh.getNumberOfTriangles(), t, i); }

template<class Vertex>
hpuint make_neighbor_index(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint i) { return make_neighbor_index(mesh.getEdges(), t, i); }

template<class Vertex>
hpuint make_neighbor_offset(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint u) { return make_neighbor_offset(mesh.getEdges(), t, u); }

template<class Vertex, Format format>
TriangleMesh<Vertex, format> make_triangle_mesh(std::vector<Vertex> vertices, Indices indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex, Format format>
TriangleMesh<Vertex, format> make_triangle_mesh(const std::string& path) { return ReaderHPH::read<TriangleMesh<Vertex, format> >(path); }

template<class Vertex>
hpuint make_valence(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v) {
     auto valence = 0u;
     visit_spokes(mesh, mesh.getOutgoing(v), [&](auto& edge) { ++valence; });
     return valence;
}

template<class Vertex, Format format>
hpuint size(const TriangleMesh<Vertex, format>& mesh) { return size(mesh.getIndices()) / 3; }

template<class Visitor>
void visit_diamonds(const std::vector<Edge>& edges, Visitor&& visit) {
     boost::dynamic_bitset<> visited(edges.size(), false);

     for(auto& edge : edges) {
          auto t = make_triangle_index(edge);
          auto i = make_edge_offset(edge);
          if(visited[3 * t + i]) continue;
          visit(t, i, edge.vertex, edges[edge.next].vertex, edges[edge.previous].vertex, edges[edges[edge.opposite].next].vertex);
          visited[edge.opposite] = true;
     }
}

template<class Vertex, class Visitor>
void visit_diamonds(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit) { visit_diamonds(mesh.getEdges(), [&](auto t, auto i, auto v0, auto v1, auto v2, auto v3) { visit(t, i, mesh.getVertex(v0), mesh.getVertex(v1), mesh.getVertex(v2), mesh.getVertex(v3)); }); }

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

template<Format format, class Visitor>
void visit_fan(FanEnumerator<format> e, Visitor&& visit) {
     while(e) {
          auto u = 0u, j = 0u;
          std::tie(u, j) = *e;
          visit(u, j);
          ++e;
     }
}

template<class Visitor>
void visit_fan(const Indices& neighbors, hpuint t, hpuint i, Visitor&& visit) { visit_fan(make_fan_enumerator(neighbors, t, i), std::forward<Visitor>(visit)); }

template<class Visitor>
void visit_fan(const std::vector<Edge>& edges, hpuint nTriangles, hpuint t, hpuint i, Visitor&& visit) {
     visit_spokes(edges, nTriangles, 3 * t + i, [&](auto& edge) {
          auto u = make_triangle_index(edge);
          if(u >= nTriangles) return;
          auto j = make_edge_offset(edge);
          visit(u, j);
     });
}

template<class Visitor>
void visit_fans(const std::vector<Edge>& edges, hpuint nTriangles, Visitor&& visit) {
     visit_vertices(edges, nTriangles, [&](auto t, auto i) {
          auto fan = make_fan(edges, nTriangles, t, i);
          visit(t, i, std::begin(fan), std::end(fan));
     });
}

template<class Vertex, class Visitor>
void visit_fans(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit) { visit_fans(mesh.getEdges(), mesh.getNumberOfTriangles(), std::forward<Visitor>(visit)); }

template<class Visitor>
void visit_ring(const std::vector<Edge>& edges, hpuint nTriangles, hpuint e, Visitor&& visit) { visit_spokes(edges, nTriangles, e, [&](auto& edge) { visit(edge.vertex); }); }

template<class Visitor>
void visit_ring(const std::vector<Edge>& edges, hpuint nTriangles, hpuint t, hpuint i, Visitor&& visit) { visit_ring(edges, nTriangles, 3 * t + i, std::forward<Visitor>(visit)); }

template<class Vertex, class Visitor>
void visit_ring(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint i, Visitor&& visit) { visit_ring(mesh.getEdges(), mesh.getNumberOfTriangles(), t, i, [&](auto v) { visit(mesh.getVertex(v)); }); }

template<class Visitor>
void visit_rings(const Indices& neighbors, Visitor&& visit) { 
     for(auto e = make_vertices_enumerator(neighbors); e; ++e) {
          auto t = 0u, i = 0u;
          std::tie(t, i) = *e;
          visit(t, i, make_ring_enumerator(neighbors, t, i));
     }
}

template<class Visitor>
void visit_rings(const std::vector<Edge>& edges, hpuint nTriangles, Visitor&& visit) {
     visit_vertices(edges, nTriangles, [&](auto t, auto i) {
          auto ring = Indices();
          visit_ring(edges, nTriangles, t, i, [&](auto v) { ring.push_back(v); });
          visit(t, i, std::begin(ring), std::end(ring));
     });
}

template<class Vertex, class Visitor>
void visit_rings(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit) {
     visit_vertices(mesh, [&](auto t, auto i) {
          auto ring = std::vector<Vertex>();
          visit_ring(mesh, t, i, [&](auto& vertex) { ring.push_back(vertex); });
          visit(t, i, std::begin(ring), std::end(ring));
     });
}

template<class Visitor>
void visit_spokes(const std::vector<Edge>& edges, hpuint nTriangles, hpuint e, Visitor&& visit) {
     auto begin = e;
     do {
          auto temp = edges[e].opposite;
          if(temp >= 3 * nTriangles) break;
          e = edges[temp].next;
     } while(e != begin);
     begin = e;

     do {
          auto& edge = edges[e];
          visit(edge);
          if(e >= 3 * nTriangles) break;
          e = edges[edge.previous].opposite;
     } while(e != begin);
}

template<class Vertex, class Visitor>
void visit_spokes(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint e, Visitor&& visit) { visit_spokes(mesh.getEdges(), mesh.getNumberOfTriangles(), e, std::forward<Visitor>(visit)); }

template<class Visitor, bool closed>
void visit_subfan(const Indices& neighbors, hpuint t, hpuint i, hpuint u, Visitor&& visit) {
     while(t != u) {
          visit(t, i);
          hpuint next;
          visit_triplet(neighbors, t, [&](hpuint n0, hpuint n1, hpuint n2) { next = (i == 0) ? n2 : (i == 1) ? n0 : n1; });
          visit_triplet(neighbors, next, [&](hpuint n0, hpuint n1, hpuint n2) { i = (n0 == t) ? 0 : (n1 == t) ? 1 : 2; });
          t = next;
     }
     visit(u, i);
}

template<class Visitor>
void visit_thorns(const std::vector<Edge>& edges, hpuint t, Visitor&& visit) {
     auto e = std::begin(edges) + 3 * t;
     visit(make_triangle_index(e[0]), make_triangle_index(e[1]), make_triangle_index(e[2]));
}

template<class Visitor>
void visit_vertices(const Indices& neighbors, Visitor&& visit) {
     for(auto e = make_vertices_enumerator(neighbors); e; ++e) {
          auto t = 0u, i = 0u;
          std::tie(t, i) = *e;
          visit(t, i);
     }
}

template<class Visitor>
void visit_vertices(const std::vector<Edge>& edges, hpuint nTriangles, Visitor&& visit) {
     boost::dynamic_bitset<> visited(3 * nTriangles, false);

     auto do_visit_vertices = [&](auto t, auto i) {
          visit(t, i);
          visit_fan(edges, nTriangles, t, i, [&](auto u, auto j) { visited[3 * u + j] = true; });
     };

     for(auto t : boost::irange(0u, nTriangles)) {
          if(!visited[3 * t]) do_visit_vertices(t, 0);
          if(!visited[3 * t + 1]) do_visit_vertices(t, 1);
          if(!visited[3 * t + 2]) do_visit_vertices(t, 2);
     }
}

template<class Vertex, class Visitor>
void visit_vertices(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit) { visit_vertices(mesh.getEdges(), mesh.getNumberOfTriangles(), std::forward<Visitor>(visit)); }

//WORKSPACE

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

template<class Mesh, class Space = typename Mesh::SPACE, class Vertex = typename Mesh::VERTEX, typename = void>
struct is_triangle_mesh : public std::false_type {};

template<class Mesh, class Space, class Vertex>
struct is_triangle_mesh<Mesh, Space, Vertex, typename std::enable_if<std::is_base_of<TriangleMesh<Vertex>, Mesh>::value && std::is_base_of<typename Mesh::SPACE, Space>::value>::type> : public std::true_type {};

}//namespace happah

