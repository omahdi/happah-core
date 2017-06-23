// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <random>
#include <string>
#include <stack>

#include "happah/format/hph.h"
#include "happah/geometries/Geometry.h"
#include "happah/geometries/Mesh.h"
#include "happah/geometries/TriangleMeshUtils.h"
#include "happah/math/Space.h"
#include "happah/utils/DeindexedArray.h"
#include "happah/utils/visitors.h"

namespace happah {

//DECLARATIONS

enum class Format { DIRECTED_EDGE, SIMPLE };

template<class Vertex, Format t_format = Format::SIMPLE>
class TriangleMesh;

struct Edge;

namespace trm {

template<Format t_format>
class SpokesWalker;

template<Format t_format>
class FanEnumerator;

template<Format t_format>
class RingEnumerator;

template<Format t_format>
class SpokesEnumerator;

template<Format t_format>
class VerticesEnumerator;

}//namespace trm

//Returns path that cuts the mesh into a disk.
template<class Vertex>
Indices cut(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh);

template<class Vertex, class Picker>
Indices cut(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpindex t, Picker&& pick);

bool is_neighbor(const Indices& neighbors, hpuint t, hpuint u);

//Return the index of this edge in the edges array.
hpuint make_edge_index(const Edge& edge);

template<class Vertex>
boost::optional<hpindex> make_edge_index(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpindex v0, hpindex v1);

hpindex make_edge_offset(hpindex e);

//Return the offset of this edge among the three edges of its adjacent triangle.
hpuint make_edge_offset(const Edge& edge);

std::vector<Edge> make_edges(const Indices& indices);

template<Format format>
Indices make_fan(trm::FanEnumerator<format> e);

template<Format format, class Iterator, class T = typename std::iterator_traits<Iterator>::value_type>
std::vector<T> make_fan(trm::FanEnumerator<format> e, Iterator begin);

template<class Vertex>
Indices make_fan(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v);

template<class Vertex>
Indices make_fan(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v);

trm::FanEnumerator<Format::SIMPLE> make_fan_enumerator(const Indices& neighbors, hpuint t, hpuint i);

trm::FanEnumerator<Format::DIRECTED_EDGE> make_fan_enumerator(const std::vector<Edge>& edges, hpuint e);

template<class Vertex>
trm::FanEnumerator<Format::SIMPLE> make_fan_enumerator(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v);

template<class Vertex>
trm::FanEnumerator<Format::DIRECTED_EDGE> make_fan_enumerator(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v);

//Return the index of the ith neighbor of the tth triangle.
hpindex make_neighbor_index(const Indices& neighbors, hpuint t, hpuint i);

//Return the index of the ith neighbor of the tth triangle.
hpindex make_neighbor_index(const std::vector<Edge>& edges, hpuint t, hpuint i);

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

Indices make_neighbors(const Indices& indices);

Indices make_neighbors(const std::vector<Edge>& edges, hpuint nTriangles);

template<class Vertex>
Indices make_neighbors(const TriangleMesh<Vertex, Format::SIMPLE>& mesh);

template<class Vertex>
Indices make_neighbors(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh);

template<Format format>
Indices make_ring(trm::RingEnumerator<format> e);

template<Format format, class Iterator, class T = typename std::iterator_traits<Iterator>::value_type>
std::vector<T> make_ring(trm::RingEnumerator<format> e, Iterator begin);

template<Format format, class T>
std::vector<T> make_ring(trm::RingEnumerator<format> e, const std::vector<T>& ts);

Indices make_ring(const std::vector<Edge>& edges, hpindex e);

template<class Vertex>
std::vector<Vertex> make_ring(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v);

template<class Vertex>
std::vector<Vertex> make_ring(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v);

trm::RingEnumerator<Format::SIMPLE> make_ring_enumerator(const Indices& neighbors, hpuint t, hpuint i);

trm::RingEnumerator<Format::DIRECTED_EDGE> make_ring_enumerator(const std::vector<Edge>& edges, hpuint e);

template<class Vertex>
trm::RingEnumerator<Format::SIMPLE> make_ring_enumerator(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v);

template<class Vertex>
trm::RingEnumerator<Format::DIRECTED_EDGE> make_ring_enumerator(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v);

trm::SpokesEnumerator<Format::SIMPLE> make_spokes_enumerator(const Indices& neighbors, hpuint t, hpuint i);

trm::SpokesEnumerator<Format::DIRECTED_EDGE> make_spokes_enumerator(const std::vector<Edge>& edges, hpuint e);

template<class Vertex>
trm::SpokesEnumerator<Format::SIMPLE> make_spokes_enumerator(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v);

template<class Vertex>
trm::SpokesEnumerator<Format::DIRECTED_EDGE> make_spokes_enumerator(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v);

trm::SpokesWalker<Format::SIMPLE> make_spokes_walker(const Indices& neighbors, hpindex t, hpindex i);

trm::SpokesWalker<Format::DIRECTED_EDGE> make_spokes_walker(const std::vector<Edge>& edges, hpindex e);

hpindex make_triangle_index(hpindex e);

hpindex make_triangle_index(const Indices& indices, hpindex v);

hpindex make_triangle_index(const Edge& edge);

template<class Vertex, Format format = Format::SIMPLE>
TriangleMesh<Vertex, format> make_triangle_mesh(std::vector<Vertex> vertices, Indices indices);

template<class Vertex = VertexP3, Format format = Format::SIMPLE>
TriangleMesh<Vertex, format> make_triangle_mesh(const std::string& path);

template<Format format>
hpuint make_valence(trm::FanEnumerator<format> e);

template<Format format>
hpuint make_valence(trm::RingEnumerator<format> e);

template<Format format>
hpuint make_valence(trm::SpokesEnumerator<format> e);

hpuint make_valence(const Indices& neighbors, hpindex t, hpindex i);

template<class Vertex>
hpuint make_valence(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v);

template<class Vertex>
hpuint make_valence(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v);

template<class Vertex, Format format>
Indices make_valences(const TriangleMesh<Vertex, format>& mesh);

hpindex make_vertex_offset(const Indices& indices, hpindex t, hpindex v);

trm::VerticesEnumerator<Format::SIMPLE> make_vertices_enumerator(const Indices& neighbors);

template<class Vertex, Format format>
hpuint size(const TriangleMesh<Vertex, format>& mesh);

//Remove any useless branches in a path.
template<class Vertex>
Indices trim(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Indices path);

template<class Vertex>
hpuint validate_cut(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, const Indices& path);

//Make sure the next edge is really in the next ring.
template<bool loop, class Vertex>
bool validate_path(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, const Indices& path);

template<class Visitor>
void visit_diamonds(const std::vector<Edge>& edges, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_diamonds(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit);

template<class Visitor>
void visit_edges(const Indices& neighbors, Visitor&& visit);

template<Format format, class Visitor>
void visit_fan(trm::FanEnumerator<format> e, Visitor&& visit);

//Visit the fan about the ith vertex of the tth triangle.
template<class Vertex, class Visitor>
void visit_fan(const Indices& neighbors, hpuint t, hpuint i, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_fan(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_fan(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v, Visitor&& visit);

template<Format format, class Visitor>
void visit_ring(trm::RingEnumerator<format> e, Visitor&& visit);

template<class Visitor>
void visit_ring(const Indices& neighbors, hpuint t, hpuint i, Visitor&& visit);

template<class Visitor>
void visit_ring(const std::vector<Edge>& edges, hpuint e, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_ring(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v, Visitor&& visit);

//TODO: refactor visit_rings
template<class Visitor>
void visit_rings(const Indices& neighbors, Visitor&& visit);

template<class Visitor>
void visit_rings(const std::vector<Edge>& edges, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_rings(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit);

template<Format format, class Visitor>
void visit_spokes(trm::SpokesEnumerator<format> e, Visitor&& visit);

template<class Visitor>
void visit_spokes(const Indices& neighbors, hpuint t, hpuint i, Visitor&& visit);

template<class Visitor>
void visit_spokes(const std::vector<Edge>& edges, hpuint e, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_spokes(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v, Visitor&& visit);

template<class Vertex, class Visitor>
void visit_spokes(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v, Visitor&& visit);

template<class Visitor>
void visit_subfan(const Indices& neighbors, hpuint t, hpuint i, hpuint u, Visitor&& visit);

template<class Visitor>
void visit_vertices(const Indices& neighbors, Visitor&& visit);

template<class Visitor>
void visit_vertices(const std::vector<Edge>& edges, Visitor&& visit);

//DEFINITIONS

template<class Vertex>
class TriangleMesh<Vertex, Format::SIMPLE> {
public:
     TriangleMesh() {}

     TriangleMesh(std::vector<Vertex> vertices, Indices indices)
          : m_indices(std::move(indices)), m_vertices(std::move(vertices)) {}

     const Indices& getIndices() const { return m_indices; }

     Indices& getIndices() { return m_indices; }

     hpuint getNumberOfTriangles() const { return m_indices.size() / 3; }

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
     friend Stream& operator<<(Stream& stream, const TriangleMesh<Vertex, Format::SIMPLE>& mesh) {
          using happah::format::hph::operator<<;

          stream << mesh.m_vertices << '\n';
          stream << mesh.m_indices;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, TriangleMesh<Vertex, Format::SIMPLE>& mesh) {
          using happah::format::hph::operator>>;

          stream >> mesh.m_vertices;
          stream >> mesh.m_indices;
          return stream;
     }

};//TriangleMesh

struct Edge {
     hpuint next;
     hpuint opposite;
     hpuint previous;
     hpuint vertex;//vertex to which edge points

     Edge(hpuint vertex, hpuint next, hpuint opposite, hpuint previous)
          : next(next), opposite(opposite), previous(previous), vertex(vertex) {}

};

//TODO: remove indices in directed_edge format
template<class Vertex>
class TriangleMesh<Vertex, Format::DIRECTED_EDGE> {
public:
     TriangleMesh() {}

     //NOTE: Indices all have to be arranged counterclockwise.
     TriangleMesh(std::vector<Vertex> vertices, Indices indices)
          : m_edges(make_edges(indices)), m_indices(std::move(indices)), m_outgoing(vertices.size(), UNULL), m_vertices(std::move(vertices)) { std::for_each(std::begin(m_edges), std::begin(m_edges) + m_indices.size(), [&](auto& edge) { m_outgoing[edge.vertex] = edge.opposite; }); }

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

     auto& getEdges() const { return m_edges; }

     const Indices& getIndices() const { return m_indices; }

     Indices& getIndices() { return m_indices; }

     auto getNumberOfEdges() const { return m_edges.size(); }

     auto getNumberOfTriangles() const { return m_indices.size() / 3; }

     hpuint getNumberOfVertices() const { return m_vertices.size(); }

     auto& getOutgoing() const { return m_outgoing; }

     auto getOutgoing(hpuint v) const { return m_outgoing[v]; }

     std::tuple<const Vertex&, const Vertex&, const Vertex&> getTriangle(hpuint t) const { return std::tie(getVertex(t, 0), getVertex(t, 1), getVertex(t, 2)); }

     auto& getVertex(hpindex v) const { return m_vertices[v]; }

     auto& getVertex(hpindex v) { return m_vertices[v]; }

     auto& getVertex(hpindex t, hpindex i) const { return m_vertices[m_indices[3 * t + i]]; }

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
          auto border = m_indices.size();
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

          hpuint found = 0;
          for(auto i = m_indices.begin(), end = m_indices.end(); i != end; ++i) {//TODO: replace with simd
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
          m_indices.insert(m_indices.end(), indices, indices + 6);

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
          auto border = m_indices.size();
          auto offset = 3 * triangle;
          auto i = m_indices.cbegin() + offset;
          auto v0 = *i;
          auto& vertex0 = getVertex(v0);
          auto v1 = *(++i);
          auto& vertex1 = getVertex(v1);
          auto v2 = *(++i);
          auto& vertex2 = getVertex(v2);
          hpuint vn = m_vertices.size();

          auto e0 = *make_edge_index(*this, v0, v1);
          auto& edge0 = m_edges[e0];
          auto e1 = edge0.next;
          auto& edge1 = m_edges[e1];
          auto e2 = edge1.next;
          auto& edge2 = m_edges[e2];
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

          m_indices[offset + 2] = vn;
          hpuint indices[] = { vn, v1, v2, vn, v2, v0 };
          m_indices.insert(m_indices.end(), indices, indices + 6);
     }

private:
     std::vector<Edge> m_edges;
     Indices m_indices;
     Indices m_outgoing;
     std::vector<Vertex> m_vertices;

};//TriangleMesh

namespace trm {

template<>
class SpokesWalker<Format::SIMPLE> {
public:
     SpokesWalker(const Indices& neighbors, hpindex t, hpindex i)
          : m_i(i), m_neighbors(neighbors), m_t(t) {}

     auto operator==(const SpokesWalker& walker) const { return m_i == walker.m_i && m_t == walker.m_t; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return hpindex(3 * m_t + m_i); }

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

};//SpokesWalker<Format::SIMPLE>

template<>
class SpokesWalker<Format::DIRECTED_EDGE> {
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

};//SpokesWalker<Format::DIRECTED_EDGE>

/*
 *   auto e = make_..._enumerator(...);
 *   do {
 *        ...
 *   } while(++e);
 */
template<Format t_format>
class FanEnumerator {
public:
     FanEnumerator(SpokesWalker<t_format> i)
          : m_begin(i), m_i(std::move(i)) {}

     explicit operator bool() const { return m_i != m_begin; }

     auto operator*() const { return hpindex(*m_i / 3); }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     SpokesWalker<t_format> m_begin;
     SpokesWalker<t_format> m_i;

};//FanEnumerator

template<Format t_format>
class RingEnumerator {
public:
     RingEnumerator(SpokesWalker<t_format> i)
          : m_begin(i), m_i(std::move(i)) {}

     explicit operator bool() const { return m_i != m_begin; }

     auto operator*() const {
          static constexpr hpindex o[3] = { 1, 2, 0 };
          auto e = *m_i;
          auto t = hpindex(e / 3);
          auto i = e - 3 * t;
          return hpindex(3 * t + o[i]);
     }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     SpokesWalker<t_format> m_begin;
     SpokesWalker<t_format> m_i;

};//RingEnumerator

template<Format t_format>
class SpokesEnumerator {
public:
     SpokesEnumerator(SpokesWalker<t_format> i)
          : m_begin(i), m_i(std::move(i)) {}

     explicit operator bool() const { return m_i != m_begin; }

     auto operator*() const { return *m_i; }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     SpokesWalker<t_format> m_begin;
     SpokesWalker<t_format> m_i;

};//SpokesEnumerator

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
          visit_spokes(m_neighbors, t, i, [&](auto e) { m_visited[e] = true; });
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
Indices cut(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh) {
     auto cache = boost::dynamic_bitset<>(mesh.getNumberOfEdges(), false);
     auto range = std::mt19937();

     cache[0] = true;
     cache[1] = true;
     cache[2] = true;
     range.seed(std::random_device()());

     return cut(mesh, 0, [&](auto& neighbors) {
          //for(auto e : boost::irange(0u, hpindex(mesh.getEdges().size())))
          //     if(neighbors[e << 1] != std::numeric_limits<hpindex>::max() && neighbors[mesh.getEdge(e).opposite << 1] == std::numeric_limits<hpindex>::max()) return e;
          if(cache.count() == 0u) return std::numeric_limits<hpindex>::max();
          auto distribution = std::uniform_int_distribution<std::mt19937::result_type>(1, cache.count());
          auto n = distribution(range);
          auto e = cache.find_first();
          while(--n) e = cache.find_next(e);
          auto& edge = mesh.getEdge(mesh.getEdge(e).opposite);
          auto e0 = mesh.getEdge(edge.previous).opposite;
          auto e1 = mesh.getEdge(edge.next).opposite;
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

template<class Vertex, class Picker>
Indices cut(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpindex t, Picker&& pick) {
     auto neighbors = Indices(mesh.getNumberOfEdges() << 1, std::numeric_limits<hpindex>::max());

     neighbors[6 * t + 0] = 6 * t + 2;
     neighbors[6 * t + 1] = 6 * t + 1;
     neighbors[6 * t + 2] = 6 * t + 0;
     neighbors[6 * t + 3] = 6 * t + 2;
     neighbors[6 * t + 4] = 6 * t + 1;
     neighbors[6 * t + 5] = 6 * t + 0;

     for(auto e = pick(neighbors); e != std::numeric_limits<hpindex>::max(); e = pick(neighbors)) {
          auto& edge0 = mesh.getEdge(e);
          auto& edge1 = mesh.getEdge(edge0.opposite);
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

template<class Vertex>
boost::optional<hpindex> make_edge_index(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpindex v0, hpindex v1) {
     auto e = make_spokes_enumerator(mesh, v0);
     do if(mesh.getEdge(*e).vertex == v1) return *e; while(++e);
     return boost::none;
}

template<Format format>
Indices make_fan(trm::FanEnumerator<format> e) {
     auto fan = Indices();
     do fan.push_back(*e); while(++e);
     return fan;
}

template<Format format, class Iterator, class T = typename std::iterator_traits<Iterator>::value_type>
std::vector<T> make_fan(trm::FanEnumerator<format> e, Iterator begin) {
     auto fan = std::vector<T>();
     do fan.push_back(begin[*e]); while(++e);
     return fan;
}

template<class Vertex>
Indices make_fan(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v) { return make_fan(make_fan_enumerator(mesh, neighbors, v)); }

template<class Vertex>
Indices make_fan(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v) { return make_fan(make_fan_enumerator(mesh, v)); }

template<class Vertex>
trm::FanEnumerator<Format::SIMPLE> make_fan_enumerator(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v) { return { { neighbors, make_triangle_index(mesh.getIndices(), v), v } }; }

template<class Vertex>
trm::FanEnumerator<Format::DIRECTED_EDGE> make_fan_enumerator(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v) { return { { mesh.getEdges(), mesh.getOutgoing(v) } }; }

template<class Vertex>
hpuint make_neighbor_index(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint i) { return make_neighbor_index(mesh.getEdges(), t, i); }

template<class Vertex>
hpuint make_neighbor_offset(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint t, hpuint u) { return make_neighbor_offset(mesh.getEdges(), t, u); }

template<class Vertex>
Indices make_neighbors(const TriangleMesh<Vertex, Format::SIMPLE>& mesh) { return make_neighbors(mesh.getIndices()); }

template<class Vertex>
Indices make_neighbors(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh) { return make_neighbors(mesh.getEdges(), size(mesh)); }

template<Format format>
Indices make_ring(trm::RingEnumerator<format> e) {
     auto ring = Indices();
     do ring.push_back(*e); while(++e);
     return ring;
}

template<Format format, class Iterator, class T = typename std::iterator_traits<Iterator>::value_type>
std::vector<T> make_ring(trm::RingEnumerator<format> e, Iterator begin) {
     auto ring = std::vector<T>();
     do ring.push_back(begin[*e]); while(++e);
     return ring;
}

template<Format format, class T>
std::vector<T> make_ring(trm::RingEnumerator<format> e, const std::vector<T>& ts) { return make_ring(e, std::begin(ts)); }

template<class Vertex>
std::vector<Vertex> make_ring(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v) { return make_ring(make_ring_enumerator(mesh, neighbors, v), std::begin(deindex(mesh.getVertices(), mesh.getVertices()))); }

template<class Vertex>
std::vector<Vertex> make_ring(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v) { return make_ring(make_ring_enumerator(mesh, v), std::begin(deindex(mesh.getVertices(), mesh.getIndices()))); }//TODO: EnumeratorTransformer<Enumerator, Transformer>

template<class Vertex>
trm::RingEnumerator<Format::SIMPLE> make_ring_enumerator(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v) { return { { neighbors, make_triangle_index(mesh.getIndices(), v), v } }; }

template<class Vertex>
trm::RingEnumerator<Format::DIRECTED_EDGE> make_ring_enumerator(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v) { return { { mesh.getEdges(), mesh.getOutgoing(v) } }; }

template<class Vertex>
trm::SpokesEnumerator<Format::SIMPLE> make_spokes_enumerator(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v) { return { { neighbors, make_triangle_index(mesh.getIndices(), v), v } }; }

template<class Vertex>
trm::SpokesEnumerator<Format::DIRECTED_EDGE> make_spokes_enumerator(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v) { return { { mesh.getEdges(), mesh.getOutgoing(v) } }; }

template<class Vertex, Format format>
TriangleMesh<Vertex, format> make_triangle_mesh(std::vector<Vertex> vertices, Indices indices) { return { std::move(vertices), std::move(indices) }; }

template<class Vertex, Format format>
TriangleMesh<Vertex, format> make_triangle_mesh(const std::string& path) { return format::hph::read<TriangleMesh<Vertex, format> >(path); }

template<Format format>
hpuint make_valence(trm::FanEnumerator<format> e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

template<Format format>
hpuint make_valence(trm::RingEnumerator<format> e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

template<Format format>
hpuint make_valence(trm::SpokesEnumerator<format> e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

template<class Vertex>
hpuint make_valence(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v) { return make_valence(make_spokes_enumerator(mesh, neighbors, v)); }

template<class Vertex>
hpuint make_valence(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v) { return make_valence(make_spokes_enumerator(mesh, v)); }

template<class Vertex, Format format>
Indices make_valences(const TriangleMesh<Vertex, format>& mesh) {
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

template<class Vertex, Format format>
hpuint size(const TriangleMesh<Vertex, format>& mesh) { return mesh.getNumberOfTriangles(); }

template<class Vertex>
Indices trim(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Indices path) {
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
          if(*j == mesh.getEdge(*i).opposite) {
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
hpuint validate_cut(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, const Indices& path) {
     auto cache = path;
     auto todo = std::stack<hpindex>();
     auto visited = boost::dynamic_bitset<>(mesh.getNumberOfTriangles(), false);

     auto contains = [&](auto e) { return std::binary_search(std::begin(cache), std::end(cache), e); };
     auto push = [&](auto e) { if(!contains(e)) todo.push(make_triangle_index(mesh.getEdge(e).opposite)); };

     if(path.size() == 0) return 0;
     if(!validate_path<true>(mesh, path)) return 1;
     std::sort(std::begin(cache), std::end(cache));

     //Make sure each edge is exactly once on the path.
     if(std::unique(std::begin(cache), std::end(cache)) != std::end(cache)) return 2;

     //Make sure opposite edge is always on path.
     for(auto e : path) if(!contains(mesh.getEdge(e).opposite)) return 3;

     //Make sure the path does not cross itself.
     for(auto i = std::begin(path), end = std::end(path) - 1; i != end; ++i) {
          if(mesh.getEdge(i[0]).opposite == i[1]) continue;
          auto e = make_spokes_walker(mesh.getEdges(), i[1]);
          while(!contains(*(++e)));
          if(mesh.getEdge(*e).opposite != i[0]) return 4;
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
bool validate_path(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, const Indices& path) {
     auto exists = [&](auto e0, auto e1) {
          auto e = make_spokes_enumerator(mesh.getEdges(), e1);
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
void visit_diamonds(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit) { visit_diamonds(mesh.getEdges(), [&](auto e, auto v0, auto v1, auto v2, auto v3) { visit(e, mesh.getVertex(v0), mesh.getVertex(v1), mesh.getVertex(v2), mesh.getVertex(v3)); }); }

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
void visit_fan(trm::FanEnumerator<format> e, Visitor&& visit) { do visit(*e); while(++e); }

template<class Visitor>
void visit_fan(const Indices& neighbors, hpuint t, hpuint i, Visitor&& visit) { visit_fan(make_fan_enumerator(neighbors, t, i), std::forward<Visitor>(visit)); }

template<class Visitor>
void visit_fan(const std::vector<Edge>& edges, hpuint e, Visitor&& visit) { visit_fan(make_fan_enumerator(edges, e), std::forward<Visitor>(visit)); }

template<Format format, class Visitor>
void visit_ring(trm::RingEnumerator<format> e, Visitor&& visit) { do visit(*e); while(++e); }

template<class Visitor>
void visit_ring(const Indices& neighbors, hpuint t, hpuint i, Visitor&& visit) { visit_ring(make_ring_enumerator(neighbors, t, i), std::forward<Visitor>(visit)); }

template<class Visitor>
void visit_ring(const std::vector<Edge>& edges, hpuint e, Visitor&& visit) { visit_ring(make_ring_enumerator(edges, e), std::forward<Visitor>(visit)); }

template<class Vertex, class Visitor>
void visit_ring(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v, Visitor&& visit) { visit_ring(make_ring_enumerator(mesh, v), [&](auto v) { visit(mesh.getVertex(mesh.getIndices()[v])); }); }

template<class Visitor>
void visit_rings(const Indices& neighbors, Visitor&& visit) { 
     for(auto e = make_vertices_enumerator(neighbors); e; ++e) {
          auto t = 0u, i = 0u;
          std::tie(t, i) = *e;
          visit(t, i, make_ring_enumerator(neighbors, t, i));
     }
}

template<class Visitor>
void visit_rings(const std::vector<Edge>& edges, Visitor&& visit) { visit_vertices(edges, [&](auto v) { visit(v, make_ring_enumerator(edges, v)); }); }

template<class Vertex, class Visitor>
void visit_rings(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor&& visit) { for(auto v : boost::irange(0u, mesh.getNumberOfVertices())) visit(v, make_ring_enumerator(mesh, v)); }

//TODO: index v can be of two forms: 1. 3t+i or 2. index into vertices array; how to keep consistent to avoid confusion?
template<Format format, class Visitor>
void visit_spokes(trm::SpokesEnumerator<format> e, Visitor&& visit) { do visit(*e); while(++e); }

template<class Visitor>
void visit_spokes(const Indices& neighbors, hpuint t, hpuint i, Visitor&& visit) { visit_spokes(make_spokes_enumerator(neighbors, t, i), std::forward<Visitor>(visit)); }

template<class Visitor>
void visit_spokes(const std::vector<Edge>& edges, hpuint e, Visitor&& visit) { visit_spokes(make_spokes_enumerator(edges, e), std::forward<Visitor>(visit)); }

template<class Vertex, class Visitor>
void visit_spokes(const TriangleMesh<Vertex, Format::SIMPLE>& mesh, const Indices& neighbors, hpuint v, Visitor&& visit) {
     visit_spokes(make_spokes_enumerator(mesh, neighbors, v), [&](auto e) {
          auto t = make_triangle_index(e);
          auto i = make_edge_offset(e);
          visit(t, i);
     });
}

template<class Vertex, class Visitor>
void visit_spokes(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint v, Visitor&& visit) { visit_spokes(make_spokes_enumerator(mesh, v), [&](auto e) { visit(mesh.getEdge(e)); }); }

template<class Visitor>
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
void visit_vertices(const Indices& neighbors, Visitor&& visit) {
     for(auto e = make_vertices_enumerator(neighbors); e; ++e) {
          auto t = 0u, i = 0u;
          std::tie(t, i) = *e;
          visit(t, i);
     }
}

template<class Visitor>
void visit_vertices(const std::vector<Edge>& edges, Visitor&& visit) {
     boost::dynamic_bitset<> visited(edges.size(), false);

     for(auto e : boost::irange(0lu, edges.size())) {
          if(visited[e]) continue;
          visit(e);
          visit_spokes(edges, e, [&](auto e) { visited[e] = true; });
     }
}

}//namespace happah

