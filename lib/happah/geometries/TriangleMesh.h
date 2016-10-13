// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cilk/cilk.h>
#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_map>

#include "happah/geometries/Geometry.h"
#include "happah/geometries/Mesh.h"
#include "happah/math/Space.h"
#include "happah/utils/DeindexedArray.h"
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

enum class Cache { NEIGHBORS, RINGS };//TODO: figure out better name for triangle neighbors? because can be confused with vertex neighbors

enum class Format { DIRECTED_EDGE, SIMPLE };

namespace {

namespace Mode = happah::mode::TriangleMesh;
namespace View = happah::view::TriangleMesh;

}//namespace

template<class Vertex, Format t_format = Format::SIMPLE, Cache... t_caches>
class TriangleMesh;

//TODO: rename iterators in meshes consistently
//TODO: shortestpathfinder and weigher now use trianglemeshutils mode and view although only mesh required as input
//TODO: view "vertex" makes less sense here; maybe get rid of it
template<class Vertex>
class TriangleMesh<Vertex, Format::SIMPLE> : public Geometry2D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;

public:
     using Indices = typename Mesh<Vertex>::Indices;
     using Vertices = typename Mesh<Vertex>::Vertices;

private:
     template<int t_dummy = 0, hpuint... t_modes>//TODO: introduce views here as well
     class Iterator {
     public:
          using difference_type = hpuint;
          using value_type = std::tuple<typename Iterator<t_dummy, t_modes>::value_type...>;

          Iterator(const TriangleMesh& triangles, hpuint triangle)
               : m_i(Iterator<t_dummy, t_modes>(triangles, triangle)...) {}

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
          std::tuple<Iterator<t_dummy, t_modes>...> m_i;

          template<unsigned long... Is>
          value_type getValue(std::index_sequence<Is...>) const { return std::make_tuple(*(std::get<Is>(m_i))...); }

          BUILD_TUPLE_HANDLER_METHODS(predecrement, doPredecrement)

          void predecrement() { predecrement(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<hpuint mode>
          void doPredecrement(Iterator<t_dummy, mode>& i) { --i; }

          BUILD_TUPLE_HANDLER_METHODS(preincrement, doPreincrement)

          void preincrement() { preincrement(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<hpuint mode>
          void doPreincrement(Iterator<t_dummy, mode>& i) { ++i; }

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, Mode::TRIANGLES> {
     public:
          using difference_type = hpuint;
          using value_type = std::tuple<hpuint, hpuint, hpuint>;

          Iterator(const TriangleMesh& triangles, hpuint triangle) 
               : m_triangle(triangle), m_triangles(triangles) {}

          difference_type operator-(const Iterator& iterator) const { return m_triangle - iterator.m_triangle; }

          Iterator& operator++() {
               ++m_triangle;
               return *this;
          }

          Iterator& operator--() {
               --m_triangle;
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

          bool operator!=(const Iterator& iterator) const { return iterator.m_triangle != m_triangle; }

          value_type operator*() const { return m_triangles.getNeighbors(m_triangle); }

     private:
          hpuint m_triangle;
          const TriangleMesh& m_triangles;

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, Mode::VERTICES> {
          using ProxyIterator = typename DeindexedArray<std::vector<Vertex>, std::vector<hpuint> >::const_iterator;

     public:
          using difference_type = typename ProxyIterator::difference_type;
          using value_type = std::tuple<const Vertex&, const Vertex&, const Vertex&>;

          Iterator(const TriangleMesh& triangles, hpuint triangle) 
               : m_i(i(triangles, triangle)) {}

          difference_type operator-(const Iterator& iterator) const { return (m_i - iterator.m_i) / 3; }

          Iterator& operator++() { 
               m_i += 3;
               return *this;
          }

          Iterator& operator--() { 
               m_i -= 3;
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

          value_type operator*() const { return std::make_tuple(std::cref(*m_i), std::cref(*(m_i + 1)), std::cref(*(m_i + 2))); }

     private:
          ProxyIterator m_i;

          static ProxyIterator i(const TriangleMesh& triangles, hpuint triangle) { return deindex(triangles.m_vertices, triangles.m_indices).cbegin() + 3 * triangle; }

     };//Iterator

public:
     TriangleMesh(Vertices vertices, Indices indices)
          : Geometry2D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)) {}

     template<Format format, Cache... caches>
     TriangleMesh(const TriangleMesh<Vertex, format, caches...>& mesh)
          : TriangleMesh(mesh.getVertices(), mesh.getIndices()) {}

     template<Format format, Cache... caches>
     TriangleMesh(TriangleMesh<Vertex, format, caches...>&& mesh)
          : TriangleMesh(std::move(mesh.getVertices()), std::move(mesh.getIndices())) {}

     virtual ~TriangleMesh() {}

     std::tuple<hpuint, hpuint, hpuint> getNeighbors(hpuint triangle) const { return TriangleMeshUtils::getNeighbors(this->m_indices, triangle); }

     template<hpuint... modes>
     Iterator<0, modes...> cbegin() const { return Iterator<0, modes...>(*this, 0); };

     template<hpuint... modes>
     Iterator<0, modes...> cend() const { return Iterator<0, modes...>(*this, this->m_indices.size() / 3); };

     void computeFlatNormals() {//TODO: remove
          for(auto i = this->m_indices.cbegin(), end = this->m_indices.cend(); i != end; ++i) {
               Vertex& v0 = this->m_vertices[*i];
               Vertex& v1 = this->m_vertices[*(++i)];
               Vertex& v2 = this->m_vertices[*(++i)];
               v2.normal = glm::cross(v1.position - v0.position, v2.position - v0.position);
          }
     }

};//TriangleMesh
using TriangleMesh2D = TriangleMesh<VertexP2>;
using TriangleMesh2D_ptr = std::shared_ptr<TriangleMesh2D>;
using TriangleMesh3D = TriangleMesh<VertexP3N>;
using TriangleMesh3D_ptr = std::shared_ptr<TriangleMesh3D>;

template<class Vertex>
class TriangleMesh<Vertex, Format::DIRECTED_EDGE> : public Geometry2D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;

     struct Edge {
          hpuint next;
          hpuint opposite;
          hpuint previous;
          hpuint vertex;//vertex to which edge points

          Edge(hpuint vertex, hpuint next, hpuint previous)
               : next(next), opposite(-1), previous(previous), vertex(vertex) {}

          Edge(hpuint vertex, hpuint next, hpuint opposite, hpuint previous)
               : next(next), opposite(opposite), previous(previous), vertex(vertex) {}

     }; 

public:
     using Indices = typename Mesh<Vertex>::Indices;
     using Vertices = typename Mesh<Vertex>::Vertices;

private:
     template<int t_dummy, hpuint t_view, hpuint... t_modes>
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
     //NOTE: Indices all have to be arranged counterclockwise or clockwise.
     TriangleMesh(Vertices vertices, Indices indices)
          : Geometry2D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)), m_border(0), m_outgoing(this->m_vertices.size(), -1) {
          hpuint nVertices = this->m_vertices.size();
          hpuint nEdges = 6 * this->m_vertices.size();//NOTE: 2 * 3 * nVertices is an approximation of the upper bound based on Euler's formula.
          m_edges.reserve(nEdges);

          using Key = std::pair<hpuint, hpuint>;
          using Value = hpuint;

          auto getHash = [](const Key& k) -> uint64_t {
               int32_t d = k.first - k.second;
               int32_t min = k.second + (d & d >> 31);
               int32_t max = k.first - (d & d >> 31);
               return ((uint64_t)max << 32 | min);
          };

          auto isKeysEqual = [](const Key& k1, const Key& k2) { return (k1.first == k2.first && k1.second == k2.second) || (k1.first == k2.second && k1.second == k2.first); };

          using Map = std::unordered_map<Key, Value, decltype(getHash), decltype(isKeysEqual)>;
          Map map(nEdges, getHash, isKeysEqual);
          hpuint edge;

          std::vector<hpuint> degrees(nVertices, 0);

          auto addEdge = [&](hpuint va, hpuint vb, hpuint next, hpuint previous) {
               m_outgoing[va] = edge;
               Key key(va, vb);
               auto i = map.find(key);
               if(i == map.end()) {
                    ++(degrees[va]);
                    ++(degrees[vb]);
                    map[key] = edge;
                    m_edges.push_back({ vb, next, previous });
               } else {
                    hpuint opposite = (*i).second;
                    m_edges[opposite].opposite = edge;
                    m_edges.push_back({ vb, next, opposite, previous });
               }
          };

          edge = 0;
          for(auto i = this->m_indices.cbegin(), end = this->m_indices.cend(); i != end; ++i) {
               hpuint v1 = *i;
               hpuint v2 = *(++i);
               hpuint v3 = *(++i);

               hpuint h1 = edge;
               hpuint h2 = edge + 1;
               hpuint h3 = edge + 2;

               addEdge(v1, v2, h2, h3);
               ++edge;
               addEdge(v2, v3, h3, h1);
               ++edge;
               addEdge(v3, v1, h1, h2);
               ++edge;
          }

          m_maxDegree = *std::max_element(degrees.begin(), degrees.end());
          
          edge = 0;
          auto i = m_edges.begin();
          auto end = m_edges.end();
          while((*i).opposite != -1 && i != end) {
               ++i;
               ++edge;
          }

          m_border = m_edges.size();
          if(i != end) {
               hpuint first = edge;
               hpuint next = m_border;
               hpuint previous = next - 2;
               do {
                    (*i).opposite = next;
                    i = m_edges.begin() + (*i).previous;
                    m_edges.push_back(Edge((*i).vertex, ++next, edge, ++previous));
                    edge = m_edges[(*i).opposite].previous;
                    i = m_edges.begin() + edge;
               } while(edge != first);
               m_edges[m_border].previous = m_edges.size() - 1;
               m_edges.back().next = m_border;
          }
     }

     template<Format format, Cache... caches>
     TriangleMesh(const TriangleMesh<Vertex, format, caches...>& mesh)
          : TriangleMesh(mesh.getVertices(), mesh.getIndices()) {}

     template<Format format, Cache... caches>
     TriangleMesh(TriangleMesh<Vertex, format, caches...>&& mesh)
          : TriangleMesh(std::move(mesh.getVertices()), std::move(mesh.getIndices())) {}

     virtual ~TriangleMesh() {}

     template<hpuint view, hpuint... modes>
     Iterator<0, view, modes...> cbegin(hpuint index) const { return Iterator<0, view, modes...>(*this, index); }

     template<hpuint view, hpuint... modes>
     bool contains(hpuint me, const typename Iterator<0, view, modes...>::value_type& value) const {
          auto i = cbegin<view, modes...>(me);
          auto first = i;
          do if(*i == value) return true;
          while((++i) != first);
          return false;
     }

     template<class Iterator>
     void exsect(Iterator begin, Iterator end)  {
          --end;
          while(++begin != end) {
               auto i = cbegin<View::VERTEX, Mode::EDGES>(*begin);
               auto p = *(begin - 1);
               auto n = *(begin + 1);
               while((*i).first.vertex != p) ++i;
               while((*(++i)).first.vertex != n) splitEdge((*i).second);
               while((*(++i)).first.vertex != p) splitEdge((*i).second);
          }
     }

     hpuint getDegree(hpuint vertex) const {
          hpuint degree = 0;
          auto i = cbegin<View::VERTEX, Mode::VERTICES>(vertex);
          auto first = i;
          do ++degree;
          while(++i != first);
          return degree;
     }

     Edge& getEdge(hpuint e) { return m_edges[e]; }//TODO: should I really expose edge?

     const Edge& getEdge(hpuint e) const { return m_edges[e]; }

     boost::optional<hpuint> getEdgeIndex(hpuint v0, hpuint v1) const {//TODO: get edge id?
          auto i = cbegin<View::VERTEX, Mode::EDGES>(v0);
          auto first = i;
          do if((*i).first.vertex == v1) return (*i).second;
          while(++i != first);
          return boost::none;
     }

     hpuint getMaxDegree() const { return m_maxDegree; }

     std::tuple<hpuint, hpuint, hpuint> getNeighbors(hpuint triangle) const {
          hpuint n1 = happah::UNULL;
          hpuint n2 = happah::UNULL;
          hpuint n3 = happah::UNULL;
          hpuint temp = 3 * triangle;
          auto h = m_edges.cbegin() + temp;
          hpuint opposite;

          opposite = (*h).opposite;
          if(opposite < m_border) n1 = opposite / 3;
          opposite = (*(++h)).opposite;
          if(opposite < m_border) n2 = opposite / 3;
          opposite = (*(++h)).opposite;
          if(opposite < m_border) n3 = opposite / 3;
          
          return std::make_tuple(n1, n2, n3);
     }

     hpuint getNumberOfEdges() const { return m_edges.size(); }

     hpuint getNumberOfTriangles() const { return this->m_indices.size() / 3; }

     const std::vector<hpuint>& getOutgoing() const { return m_outgoing; }

     hpuint getOutgoing(hpuint v) const { return m_outgoing[v]; }

     std::vector<hpuint> getRing(hpuint vertex) const {
          std::vector<hpuint> neighbors;
          neighbors.reserve(m_maxDegree);
          auto i = cbegin<View::VERTEX, Mode::VERTICES>(vertex);
          auto first = i;
          do neighbors.push_back(*i);
          while(++i != first);
          return std::move(neighbors);
     }

     bool isRing(hpuint vertex, hpuint neighbor) const {
          auto i = cbegin<View::VERTEX, Mode::VERTICES>(vertex);
          auto first = i;
          do if(*i == neighbor) return true;
          while(++i != first);
          return false;
     }

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
          auto e0 = edge;
          auto& edge0 = m_edges[e0];
          assert(edge < m_border && edge0.opposite < m_border);//TODO: edge to split is border edge or opposite is border
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

          if(m_border < m_edges.size()) {
               for(auto i = m_edges.begin(), end = m_edges.end(); i != end; ++i) {
                    Edge& temp = *i;
                    if(temp.next >= m_border) temp.next += 6;
                    if(temp.opposite >= m_border) temp.opposite += 6;
                    if(temp.previous >= m_border) temp.previous += 6;
               }
          }

          Edge edges[] = { Edge(v2, e4, n1, e0), Edge(vn, n2, n0, e2), Edge(v0, e2, n4, n1), Edge(vn, e1, n5, e3), Edge(vn, n5, n2, e5), Edge(v3, e5, n3, n4) };
          m_edges.insert(m_edges.begin() + m_border, edges, edges + 6);
          m_border += 6;
          
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

          m_maxDegree = std::max(m_maxDegree, getDegree(v2));
          m_maxDegree = std::max(m_maxDegree, getDegree(v3));
          m_maxDegree = std::max(m_maxDegree, 4u);
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

          if(m_border < m_edges.size()) {//TODO: refactor into method insertEdges
               for(auto& edge : m_edges) {
                    if(edge.next >= m_border) edge.next += 6;
                    if(edge.opposite >= m_border) edge.opposite += 6;
                    if(edge.previous >= m_border) edge.previous += 6;
               }
          }

          edge0.next = f0;
          edge0.previous = g0;
          edge1.next = f1;
          edge1.previous = g1;
          edge2.next = f2;
          edge2.previous = g2;

          Edge edges[] = { Edge(vn, g0, g1, e0), Edge(v0, e0, f2, f0), Edge(vn, g1, g2, e1), Edge(v1, e1, f0, f1), Edge(vn, g2, g0, e2), Edge(v2, e2, f1, f2) };
          m_edges.insert(m_edges.begin() + m_border, edges, edges + 6);
          m_border += 6;

          this->m_indices[offset + 2] = vn;
          hpuint indices[] = { vn, v1, v2, vn, v2, v0 };
          this->m_indices.insert(this->m_indices.end(), indices, indices + 6);

          //TODO: if triangle or edges are split and there are caches, caches must be updated

          m_maxDegree = std::max(m_maxDegree, 3u);
     }

private:
     hpuint m_border;
     std::vector<Edge> m_edges;
     hpuint m_maxDegree;
     std::vector<hpuint> m_outgoing;

};//TriangleMesh

template<class Vertex, Format t_format>
class TriangleMesh<Vertex, t_format, Cache::NEIGHBORS> : public virtual TriangleMesh<Vertex, t_format> {
public:
     using Indices = typename Mesh<Vertex>::Indices;
     using Vertices = typename Mesh<Vertex>::Vertices;

     TriangleMesh(Vertices vertices, Indices indices)
          : TriangleMesh<Vertex, t_format>(std::move(vertices), std::move(indices)), m_neighbors(TriangleMeshUtils::getNeighbors(this->m_indices)) {}

     template<Format format, Cache... caches>
     TriangleMesh(const TriangleMesh<Vertex, format, caches...>& mesh)
          : TriangleMesh(mesh.getVertices(), mesh.getIndices()) {}

     template<Format format, Cache... caches>
     TriangleMesh(TriangleMesh<Vertex, format, caches...>&& mesh)
          : TriangleMesh(std::move(mesh.getVertices()), std::move(mesh.getIndices())) {}

     std::tuple<hpuint, hpuint, hpuint> getNeighbors(hpuint triangle) const {
          auto n = m_neighbors.cbegin() + 3 * triangle;
          hpuint n0 = *n;
          hpuint n1 = *(++n);
          hpuint n2 = *(++n);
          return std::make_tuple(n0, n1, n2);
     }

private:
     std::vector<hpuint> m_neighbors;

};//TriangleMesh

template<class Vertex, Format t_format>
class TriangleMesh<Vertex, t_format, Cache::RINGS> : public virtual TriangleMesh<Vertex, t_format> {
public:
     using Indices = typename Mesh<Vertex>::Indices;
     using Vertices = typename Mesh<Vertex>::Vertices;

     TriangleMesh(Vertices vertices, Indices indices)
          : TriangleMesh<Vertex, t_format>(std::move(vertices), std::move(indices)) {
          hpuint nVertices = this->m_vertices.size();
          m_rings.reserve(nVertices);
          for(hpuint i = 0; i < nVertices; ++i) m_rings.push_back(TriangleMesh<Vertex, t_format>::getRing(i));
     }

     template<Format format, Cache... caches>
     TriangleMesh(const TriangleMesh<Vertex, format, caches...>& mesh)
          : TriangleMesh(mesh.getVertices(), mesh.getIndices()) {}

     template<Format format, Cache... caches>
     TriangleMesh(TriangleMesh<Vertex, format, caches...>&& mesh)
          : TriangleMesh(std::move(mesh.getVertices()), std::move(mesh.getIndices())) {}

     const std::vector<hpuint>& getRing(hpuint vertex) const { return m_rings[vertex]; }
          
private:
     std::vector<std::vector<hpuint> > m_rings;

};//TriangleMesh

template<Cache... t_caches>
struct Caches {};

template<class Mesh, class Space = typename Mesh::SPACE, class Vertex = typename Mesh::VERTEX, typename = void>
struct is_triangle_mesh : public std::false_type {};

template<class Mesh, class Space, class Vertex>
struct is_triangle_mesh<Mesh, Space, Vertex, typename std::enable_if<std::is_base_of<TriangleMesh<Vertex>, Mesh>::value && std::is_base_of<typename Mesh::SPACE, Space>::value>::type> : public std::true_type {};

template<class Vertex, Cache... caches>
static TriangleMesh<Vertex, Format::SIMPLE, caches...> make_triangle_mesh(std::vector<Vertex> vertices, std::vector<hpuint> indices, Caches<caches...> = Caches<caches...>()) { return TriangleMesh<Vertex, Format::SIMPLE, caches...>(std::move(vertices), std::move(indices)); }

template<class Vertex, class Visitor>
void visit_rings_and_fans(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, Visitor visit) {
     for(auto begin : mesh.getOutgoing()) {
          std::vector<hpuint> vertices, triangles;
          auto e = begin;
          do {
               auto& edge = mesh.getEdge(e);
               vertices.emplace_back(edge.vertex);
               triangles.emplace_back(e / 3);
               e = mesh.getEdge(mesh.getEdge(edge.next).next).opposite;
          } while(e != begin);
          visit(vertices.cbegin(), triangles.cbegin(), vertices.size());
     }
}

}//namespace happah

