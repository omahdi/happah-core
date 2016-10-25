// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <vector>

#include "happah/Happah.h"
#include "happah/math/TriangleDecomposition.h"
#include "happah/utils/Arrays.h"
#include "happah/utils/ShortestPathFinder.h"
#include "happah/weighers/EdgeLengthWeigher.h"
#include "happah/weighers/TraversableEdgeLengthWeigher.h"

namespace happah {

template<class Mesh>
class HexagonDecomposition {
     using Indices = std::vector<hpuint>;
     using IndicesArrays = Arrays<hpuint>;

     class TriangleDecomposer {
          using Weigher = TraversableEdgeLengthWeigher<Mesh>;

     public:
          static TriangleDecomposition<Mesh> decompose(const HexagonDecomposition& decomposition) {
               TriangleDecomposer decomposer(decomposition);
               return decomposer.decompose();
          }

     private:
          IndicesArrays m_boundaries;
          const HexagonDecomposition& m_decomposition;
          Indices m_indices;//every three indices form a triangle border
          Mesh& m_mesh;
          Indices m_neighbors;
          boost::dynamic_bitset<> m_reverse;
          std::set<hpuint> m_wallEdges;
          std::set<hpuint> m_wallVertices;
          Weigher m_weigher;

          TriangleDecomposer(const HexagonDecomposition& decomposition)
               : m_boundaries(decomposition.m_boundaries), m_decomposition(decomposition), m_mesh(m_decomposition.m_mesh), m_weigher(m_mesh) {
               namespace Mode = mode::Mesh;
               using ShortestPathFinder = ShortestPathFinder<Weigher, Mesh>;
               namespace View = view::Mesh;

               const hpuint nHexagons = m_decomposition.getNumberOfHexagons();
               const hpuint nTriangles = 6 * nHexagons;
               m_indices.reserve(3 * nTriangles);
               m_neighbors.reserve(3 * nTriangles);//TODO: introduce iterators into hexagon decomposition
               m_reverse.resize(3 * nTriangles);

               for(hpuint hexagon = 0; hexagon < nHexagons; ++hexagon) {
                    auto boundary = m_decomposition.getBoundary(hexagon);
                    extendWall(boundary.cbegin(), boundary.cend(), true);
               }
               splitDangerousEdges(m_wallVertices.cbegin(), m_wallVertices.cend());

               ShortestPathFinder shortestPathFinder(m_mesh, m_weigher);

               auto i = m_decomposition.m_indices.cbegin();
               auto n = m_decomposition.m_neighbors.cbegin();
               for(hpuint hexagon = 0; hexagon < nHexagons; ++hexagon) {
                    auto center = m_decomposition.getCenter(hexagon);
                    while(get_degree(m_mesh, center) < 6) {
                         auto i = m_mesh.template cbegin<View::VERTEX, Mode::EDGES>(center);
                         auto first = i;
                         do {
                              auto e = (*i).first.next;
                              if(m_wallEdges.find(e) != m_wallEdges.end()) continue;
                              m_mesh.splitEdge(e);
                         } while(++i != first);
                    }
                    auto first = m_boundaries.size();
                    auto b5 = m_boundaries[*(i + 5)];
                    hpuint ep = (m_decomposition.m_reverse[6 * hexagon + 5]) ? *(b5.first + 1) : *(b5.second - 2);
                    auto tp = 6 * hexagon + 5;
                    for(hpuint triangle = 6 * hexagon, end = triangle + 6; triangle < end; ++triangle) {
                         auto b = m_boundaries[*i];
                         hpuint e, en;
                         if(m_decomposition.m_reverse[triangle]) {
                              e = *(b.second - 1);
                              en = *(b.second - 2);
                         } else {
                              e = *b.first;
                              en = *(b.first + 1);
                         }
                         auto j = m_mesh.template cbegin<View::VERTEX, Mode::EDGES>(e);//make sure path can only come in on one side
                         while((*j).first.vertex != ep) ++j;
                         while((*(++j)).first.vertex != en) {
                              m_weigher.removeEdge((*j).second);
                              m_weigher.removeEdge((*j).first.opposite);
                         }
                         while((*(++j)).first.vertex != ep) m_weigher.unremoveEdge((*j).first.opposite);
                         ep = (m_decomposition.m_reverse[triangle]) ? *(b.first + 1) : *(b.second - 2);
                         m_indices.push_back(m_boundaries.size());
                         m_neighbors.push_back(tp);
                         tp = triangle;
                         m_indices.push_back(*i);
                         m_neighbors.push_back(6 * (*n) + m_decomposition.getBoundaryOffset(*n, *i));
                         m_indices.push_back(m_boundaries.size() + 1);
                         m_neighbors.push_back(triangle + 1);
                         m_reverse[3 * triangle] = true;
                         m_reverse[3 * triangle + 1] = m_decomposition.m_reverse[triangle];
                         std::vector<hpuint> temp;//get approximation of path
                         if(!shortestPathFinder.getShortestPath(center, e, m_wallVertices.cbegin(), m_wallVertices.cend(), temp)) {
                              std::cerr << "Failed to find path from center.\n";
                              continue;
                         }
                         m_mesh.exsect(temp.cbegin(), temp.cend());
                         m_weigher.resize(m_mesh.getNumberOfEdges());//TODO: better m_weigher.handleMeshChangedEvent(m_mesh); or recreate weigher and shortest path finder
                         if(!shortestPathFinder.getShortestPath(center, e, m_wallVertices.cbegin(), m_wallVertices.cend(), IndicesArrays::ArrayAppender(m_boundaries))) {
                              std::cerr << "Failed to find path from center.\n";
                              continue;
                         }
                         auto last = *(--m_boundaries.end());
                         extendWall(last.first, last.second);
                         splitDangerousEdges(last.first, last.second);
                         ++i;
                         ++n;
                    }
                    m_indices.back() = first;
                    m_neighbors.back() = 6 * hexagon;
               }
          }

          template<class Iterator>
          void extendWall(Iterator begin, Iterator end, bool loop = false) {
               for(auto i0 = begin, i1 = i0 + 1; i1 != end; i0 = i1, ++i1) {
                    auto e = *m_mesh.getEdgeIndex(*i0, *i1);
                    m_wallEdges.insert(e);
                    m_wallEdges.insert(m_mesh.getEdge(e).opposite);
                    m_wallVertices.insert(*i0);
               }
               if(loop) {
                    auto e = *m_mesh.getEdgeIndex(*begin, *(end - 1));
                    m_wallEdges.insert(e);
                    m_wallEdges.insert(m_mesh.getEdge(e).opposite);
                    m_wallVertices.insert(*(end - 1));
               }
          }

          template<class Iterator>
          void splitDangerousEdges(Iterator begin, Iterator end) {
               do {
                    auto i = m_mesh.template cbegin<View::VERTEX, Mode::EDGES>(*begin);
                    auto first = i;
                    do {
                         auto temp = *i;
                         if(m_wallEdges.find(temp.second) != m_wallEdges.end()) continue;
                         if(m_wallVertices.find(temp.first.vertex) == m_wallVertices.end()) continue;
                         m_mesh.splitEdge(temp.second);
                    } while(++i != first);
               } while(++begin != end);
               m_weigher.resize(m_mesh.getNumberOfEdges());
          }

          TriangleDecomposition<Mesh> decompose() { return { m_decomposition.m_mesh, std::move(m_boundaries), std::move(m_indices), std::move(m_neighbors), std::move(m_reverse) }; }

     };//TriangleDecomposer

public:
     //NOTE: The contract is that any changes to the mesh leave the given boundaries intact.
     HexagonDecomposition(Mesh& mesh)
          : m_mesh(mesh) {}

     //NOTE: Every six indices point to six borders of hexagon in boundaries array.  Every six neighbors point to the six hexagons that are neighbors.  The ith border is the border with the ith neighbor.  Reverse is a flag per border segment that says whether the border should be traversed in reverse order to be counterclockwise.  Some borders are shared, which is why we need these flags.
     HexagonDecomposition(Mesh& mesh, IndicesArrays boundaries, Indices indices, Indices neighbors, boost::dynamic_bitset<> reverse)
          : m_boundaries(std::move(boundaries)), m_indices(std::move(indices)), m_mesh(mesh), m_neighbors(std::move(neighbors)), m_reverse(std::move(reverse)) {
          namespace Mode = mode::Mesh;
          namespace View = view::Mesh;

          auto check = [&](hpuint hexagon, hpuint boundaryIndex, hpuint neighborIndex) -> bool {
               auto n = m_neighbors.cbegin() + (6 * hexagon);
               auto b = m_indices.cbegin() + (6 * hexagon);
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               ++n; ++b;
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               ++n; ++b;
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               ++n; ++b;
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               ++n; ++b;
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               ++n; ++b;
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               return false;
          };

          auto contains = [&](hpuint e0, hpuint e1) -> bool { return (bool)find_in_ring(m_mesh.getEdges(), m_mesh.getOutgoing(e0), e1); };

          hpuint h = 0;
          hpuint offset = 0;
          auto i = m_indices.cbegin();
          for(auto n = m_neighbors.cbegin(), end = m_neighbors.cend(); n != end; ++n, ++i, ++h, ++offset) {
               assert(check(*n, *i, h));
               auto b = m_boundaries[*i];
               auto first = (m_reverse[offset]) ? *(b.second - 1) : *b.first;
               auto e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(check(*(++n), *(++i), h));
               b = m_boundaries[*i];
               auto e1 = (m_reverse[++offset]) ? *(b.second - 1) : *b.first;
               assert(contains(e0, e1));
               e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(check(*(++n), *(++i), h));
               b = m_boundaries[*i];
               e1 = (m_reverse[++offset]) ? *(b.second - 1) : *b.first;
               assert(contains(e0, e1));
               e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(check(*(++n), *(++i), h));
               b = m_boundaries[*i];
               e1 = (m_reverse[++offset]) ? *(b.second - 1) : *b.first;
               assert(contains(e0, e1));
               e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(check(*(++n), *(++i), h));
               b = m_boundaries[*i];
               e1 = (m_reverse[++offset]) ? *(b.second - 1) : *b.first;
               assert(contains(e0, e1));
               e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(check(*(++n), *(++i), h));
               b = m_boundaries[*i];
               e1 = (m_reverse[++offset]) ? *(b.second - 1) : *b.first;
               assert(contains(e0, e1));
               e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(contains(e0, first));
          }
     }

     //NOTE: Boundary is arranged counterclockwise.
     std::vector<hpuint> getBoundary(hpuint hexagon) const {
          using r = std::reverse_iterator<IndicesArrays::const_iterator::iterator>;

          std::vector<hpuint> boundary;
          hpuint offset = 6 * hexagon;
          auto i = m_indices.cbegin() + offset;
          auto b0 = m_boundaries[*i];
          auto b1 = m_boundaries[*(++i)];
          auto b2 = m_boundaries[*(++i)];
          auto b3 = m_boundaries[*(++i)];
          auto b4 = m_boundaries[*(++i)];
          auto b5 = m_boundaries[*(++i)];
          if(m_reverse[offset]) boundary.insert(boundary.end(), r(b0.second), r(b0.first + 1));
          else boundary.insert(boundary.end(), b0.first, b0.second - 1);
          if(m_reverse[++offset]) boundary.insert(boundary.end(), r(b1.second), r(b1.first + 1));
          else boundary.insert(boundary.end(), b1.first, b1.second - 1);
          if(m_reverse[++offset]) boundary.insert(boundary.end(), r(b2.second), r(b2.first + 1));
          else boundary.insert(boundary.end(), b2.first, b2.second - 1);
          if(m_reverse[++offset]) boundary.insert(boundary.end(), r(b3.second), r(b3.first + 1));
          else boundary.insert(boundary.end(), b3.first, b3.second - 1);
          if(m_reverse[++offset]) boundary.insert(boundary.end(), r(b4.second), r(b4.first + 1));
          else boundary.insert(boundary.end(), b4.first, b4.second - 1);
          if(m_reverse[++offset]) boundary.insert(boundary.end(), r(b5.second), r(b5.first + 1));
          else boundary.insert(boundary.end(), b5.first, b5.second - 1);

          return boundary;
     }

     std::tuple<hpuint, hpuint, hpuint, hpuint, hpuint, hpuint> getCorners(hpuint hexagon) const {
          hpuint offset = 6 * hexagon;
          auto i = m_indices.cbegin() + offset;
          auto b0 = m_boundaries[*i];
          auto b1 = m_boundaries[*(++i)];
          auto b2 = m_boundaries[*(++i)];
          auto b3 = m_boundaries[*(++i)];
          auto b4 = m_boundaries[*(++i)];
          auto b5 = m_boundaries[*(++i)];
          auto c0 = (m_reverse[offset]) ? *(b0.second - 1) : *b0.first;
          auto c1 = (m_reverse[++offset]) ? *(b1.second - 1) : *b1.first;
          auto c2 = (m_reverse[++offset]) ? *(b2.second - 1) : *b2.first;
          auto c3 = (m_reverse[++offset]) ? *(b3.second - 1) : *b3.first;
          auto c4 = (m_reverse[++offset]) ? *(b4.second - 1) : *b4.first;
          auto c5 = (m_reverse[++offset]) ? *(b5.second - 1) : *b5.first;
          return std::make_tuple(c0, c1, c2, c3, c4, c5);
     }

     hpuint getCenter(hpuint hexagon) const {
          using r = std::reverse_iterator<Indices::const_iterator>;
          using Weigher = EdgeLengthWeigher<Mesh>;
          using Weight = typename Weigher::Weight;

          auto boundary = getBoundary(hexagon);
          hpuint c0, c1, c2, c3, c4, c5;
          std::tie(c0, c1, c2, c3, c4, c5) = getCorners(hexagon);
          auto vertices = getVertices(hexagon);
          auto weigher = make_edge_length_weigher(m_mesh);
          auto shortestPathFinder = make_shortest_path_finder(m_mesh, weigher);
          auto minVariance = Weigher::MAX_WEIGHT;
          hpuint center;
          for(auto vertex : vertices) {//TODO: replace with faster algorithm
               Indices path;
               Weight l0, l1, l2, l3, l4, l5;

               shortestPathFinder.getShortestPath(vertex, c0, boundary.cbegin(), boundary.cend(), path);
               if(path.size() == 2 && path[0] == path[1]) break; 
               l0 = weigher.weigh(r(path.end()), r(path.begin()));
               path.clear();
               shortestPathFinder.getShortestPath(vertex, c1, boundary.cbegin(), boundary.cend(), path);
               if(path.size() == 2 && path[0] == path[1]) break; 
               l1 = weigher.weigh(r(path.end()), r(path.begin()));
               path.clear();
               shortestPathFinder.getShortestPath(vertex, c2, boundary.cbegin(), boundary.cend(), path);
               if(path.size() == 2 && path[0] == path[1]) break; 
               l2 = weigher.weigh(r(path.end()), r(path.begin()));
               path.clear();
               shortestPathFinder.getShortestPath(vertex, c3, boundary.cbegin(), boundary.cend(), path);
               if(path.size() == 2 && path[0] == path[1]) break; 
               l3 = weigher.weigh(r(path.end()), r(path.begin()));
               path.clear();
               shortestPathFinder.getShortestPath(vertex, c4, boundary.cbegin(), boundary.cend(), path);
               if(path.size() == 2 && path[0] == path[1]) break; 
               l4 = weigher.weigh(r(path.end()), r(path.begin()));
               path.clear();
               shortestPathFinder.getShortestPath(vertex, c5, boundary.cbegin(), boundary.cend(), path);
               if(path.size() == 2 && path[0] == path[1]) break; 
               l5 = weigher.weigh(r(path.end()), r(path.begin()));

               auto mean = (l0 + l1 + l2 + l3 + l4 + l5) / 6.0;
               auto v0 = mean - l0;
               v0 *= v0;
               auto v1 = mean - l1;
               v1 *= v1;
               auto v2 = mean - l2;
               v2 *= v2;
               auto v3 = mean - l3;
               v3 *= v3;
               auto v4 = mean - l4;
               v4 *= v4;
               auto v5 = mean - l5;
               v5 *= v5;
               auto variance = (v0 + v1 + v2 + v3 + v4 + v5) / 6.0;

               if(variance < minVariance) {
                    minVariance = variance;
                    center = vertex;
               }
          }
          return center;
     }

     std::vector<hpuint> getNeighbors(hpuint hexagon) const {
          auto n = m_neighbors.cbegin() + (6 * hexagon);
          return std::vector<hpuint>(n, n + 6);
     }

     std::vector<hpuint> getVertices(hpuint hexagon) const {
          namespace Mode = mode::Mesh;
          namespace View = view::Mesh;

          std::stack<hpuint> todo;
          std::vector<hpuint> boundary = getBoundary(hexagon);

          for(auto b0 = boundary.cbegin(), b1 = b0 + 1, end = boundary.cend(); b1 != end; b0 = b1, ++b1) {
               auto i = this->m_mesh.template cbegin<View::VERTEX, Mode::EDGES>(*b0);
               while((*i).first.vertex != *b1) ++i;
               todo.push((*(++i)).first.vertex);
          }
          {
               auto i = this->m_mesh.template cbegin<View::VERTEX, Mode::EDGES>(boundary.back());
               while((*i).first.vertex != boundary[0]) ++i;
               todo.push((*(++i)).first.vertex);
          }

          std::sort(boundary.begin(), boundary.end());

          std::vector<hpuint> vertices;
          while(!todo.empty()) {
               auto t = todo.top();
               todo.pop();
               if(std::binary_search(boundary.cbegin(), boundary.cend(), t)) continue;
               auto i = this->m_mesh.template cbegin<View::VERTEX, Mode::EDGES>(t);
               auto first = i;
               do {
                    auto v = (*i).first.vertex;
                    if(std::binary_search(boundary.cbegin(), boundary.cend(), v)) break;
                    else if(std::find(vertices.cbegin(), vertices.cend(), v) == vertices.end()) {
                         vertices.push_back(v);
                         todo.push(v);
                    }
                    ++i;
               } while(i != first);
               if(i == first) continue;
               i = first;
               while((--i) != first) {
                    auto v = (*i).first.vertex;
                    if(std::binary_search(boundary.cbegin(), boundary.cend(), v)) break;
                    else if(std::find(vertices.cbegin(), vertices.cend(), v) == vertices.end()) {
                         vertices.push_back(v);
                         todo.push(v);
                    }
               }
          }
          return vertices;
     }

     hpuint getNumberOfHexagons() const { return m_indices.size() / 6; }

     TriangleDecomposition<Mesh> toTriangleDecomposition() { return TriangleDecomposer::decompose(*this); }

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const HexagonDecomposition& decomposition) {
          stream << decomposition.m_boundaries;
          stream << decomposition.m_indices;
          stream << decomposition.m_neighbors;
          stream << decomposition.m_reverse;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, HexagonDecomposition& decomposition) {
          stream >> decomposition.m_boundaries;
          stream >> decomposition.m_indices;
          stream >> decomposition.m_neighbors;
          stream >> decomposition.m_reverse;
          return stream;
     }

private:
     IndicesArrays m_boundaries;
     Indices m_indices;//every six indices form a hexagon border
     Mesh& m_mesh;
     Indices m_neighbors;
     boost::dynamic_bitset<> m_reverse;//whether the border needs to be reversed to be counterclockwise

     hpuint getBoundaryOffset(hpuint h, hpuint b) const {
          auto i = m_indices.cbegin() + (6 * h);
          if(*i == b) return 0;
          else if(*(++i) == b) return 1;
          else if(*(++i) == b) return 2;
          else if(*(++i) == b) return 3;
          else if(*(++i) == b) return 4;
          else {
               assert(*(++i) == b);
               return 5;
          }
     }

};//HexagonDecomposition

}//namespace happah

