// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/Happah.h"
#include "happah/utils/Arrays.h"
#include "happah/utils/ShortestPathFinder.h"
#include "happah/weighers/EdgeLengthWeigher.h"

namespace happah {

template<class Mesh>
class TriangleDecomposition {
     using Indices = std::vector<hpuint>;
     using IndicesArrays = Arrays<hpuint>;
     using Point = typename Mesh::SPACE::POINT;

public:
     TriangleDecomposition(Mesh& mesh)
          : m_mesh(mesh) {}

     //NOTE: The contract is that any changes to the mesh leave the given boundaries intact.
     TriangleDecomposition(Mesh& mesh, IndicesArrays boundaries, Indices indices, Indices neighbors, boost::dynamic_bitset<> reverse)
          : m_boundaries(std::move(boundaries)), m_indices(std::move(indices)), m_mesh(mesh), m_neighbors(std::move(neighbors)), m_reverse(std::move(reverse)) {
          auto check = [&](hpuint t, hpuint boundaryIndex, hpuint neighborIndex) -> bool {
               auto n = m_neighbors.cbegin() + (3 * t);
               auto b = m_indices.cbegin() + (3 * t);
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               ++n; ++b;
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               ++n; ++b;
               if(*n == neighborIndex && *b == boundaryIndex) return true;
               ++n; ++b;
               return false;
          };

          auto contains = [&](hpuint e0, hpuint e1) -> bool {
               auto e = make_ring_enumerator(mesh, e0);
               do { if(*e == e1) return true;
               } while(++e);
               return false;
          };

          hpuint t = 0;
          hpuint offset = 0;
          auto i = m_indices.cbegin();
          for(auto n = m_neighbors.cbegin(), end = m_neighbors.cend(); n != end; ++n, ++i, ++t, ++offset) {
               assert(check(*n, *i, t));
               auto b = m_boundaries[*i];
               auto first = (m_reverse[offset]) ? *(b.second - 1) : *b.first;
               auto e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(check(*(++n), *(++i), t));
               b = m_boundaries[*i];
               auto e1 = (m_reverse[++offset]) ? *(b.second - 1) : *b.first;
               assert(contains(e0, e1));
               e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(check(*(++n), *(++i), t));
               b = m_boundaries[*i];
               e1 = (m_reverse[++offset]) ? *(b.second - 1) : *b.first;
               assert(contains(e0, e1));
               e0 = (m_reverse[offset]) ? *(b.first + 1) : *(b.second - 2);
               assert(contains(e0, first));
          }
     }

     //NOTE: Boundary is arranged counterclockwise.
     std::vector<hpuint> getBoundary(hpuint t) const {
          using r = std::reverse_iterator<IndicesArrays::const_iterator::iterator>;

          std::vector<hpuint> boundary;
          auto offset = 3 * t;
          auto i = m_indices.cbegin() + offset;
          auto b0 = m_boundaries[*i];
          auto b1 = m_boundaries[*(++i)];
          auto b2 = m_boundaries[*(++i)];
          if(m_reverse[offset]) boundary.insert(boundary.end(), r(b0.second), r(b0.first + 1));
          else boundary.insert(boundary.end(), b0.first, b0.second - 1);
          if(m_reverse[++offset]) boundary.insert(boundary.end(), r(b1.second), r(b1.first + 1));
          else boundary.insert(boundary.end(), b1.first, b1.second - 1);
          if(m_reverse[++offset]) boundary.insert(boundary.end(), r(b2.second), r(b2.first + 1));
          else boundary.insert(boundary.end(), b2.first, b2.second - 1);

          return boundary;
     }

     std::tuple<hpuint, hpuint, hpuint> getCorners(hpuint t) const {
          using r = std::reverse_iterator<IndicesArrays::const_iterator::iterator>;

          auto offset = 3 * t;
          auto i = m_indices.cbegin() + offset;
          auto b0 = m_boundaries[*i];
          auto b1 = m_boundaries[*(++i)];
          auto b2 = m_boundaries[*(++i)];
          auto c0 = (m_reverse[offset]) ? *r(b0.second) : *b0.first;
          auto c1 = (m_reverse[++offset]) ? *r(b1.second) : *b1.first;
          auto c2 = (m_reverse[++offset]) ? *r(b2.second) : *b2.first;

          return std::make_tuple(c0, c1, c2);
     }

     std::tuple<hpuint, hpuint, hpuint> getNeighbors(hpuint t) const {
          auto n = m_neighbors.cbegin() + (3 * t);
          return std::make_tuple(*n, *(n + 1), *(n + 2));
     }

     hpuint getNumberOfTriangles() const { return m_indices.size() / 3; }

     //NOTE: Points have to be sorted by increasing w and then increasing v.
     template<class Iterator>
     std::vector<Point> sample(hpuint t, Iterator pointsBegin, Iterator pointsEnd, hpreal epsilon = EPSILON) const {
          using r = std::reverse_iterator<IndicesArrays::const_iterator::iterator>;
          using Weigher = EdgeLengthWeigher<Mesh>;
          using Weight = typename Weigher::Weight;

          auto offset = 3 * t;
          auto i = m_indices.cbegin() + offset;
          auto b0 = m_boundaries[*i];
          auto b1 = m_boundaries[*(++i)];
          auto b2 = m_boundaries[*(++i)];
          auto b00 = b0.first;
          auto b01 = b0.second;
          auto b10 = b1.first;
          auto b11 = b1.second;
          auto b20 = b2.first;
          auto b21 = b2.second;
          auto r0 = m_reverse[offset];
          auto r1 = m_reverse[++offset];
          auto r2 = m_reverse[++offset];

          auto c1 = (r1) ? Indices(r(b11), r(b10)) : Indices(b10, b11);
          auto c2 = (r2) ? Indices(b20, b21) : Indices(r(b21), r(b20));

          Mesh mesh(m_mesh);//TODO: can we avoid copy?
          auto weigher = make_edge_length_weigher(mesh);
          auto shortestPathFinder = make_shortest_path_finder(mesh, weigher);

          std::vector<Point> samples;
          samples.reserve(std::distance(pointsBegin, pointsEnd));

          auto l1 = weigher.weigh(c1.cbegin(), c1.cend());
          auto l2 = weigher.weigh(c2.cbegin(), c2.cend());

          auto p = pointsBegin;
          auto q = p + 1;

          hpuint e0, e1; 
          auto handle = [&](const Indices& path, Weight l) {
               while(p != q && std::abs((*p).y) < epsilon) {
                    samples.push_back(mesh.getVertices()[e0].position);
                    ++p;
               }
               while(p != q && std::abs((*p).x) > epsilon) {
                    samples.push_back(weigher.getPosition(r(path.cend()), r(path.cbegin()), l * ((*p).y / ((*p).x + (*p).y))));
                    ++p;
               }
               while(p != q) {
                    samples.push_back(mesh.getVertices()[e1].position);
                    ++p;
               }
          };

          auto position = [&](std::vector<hpuint>& b, Weight distance) -> hpuint {
               for(auto i = b.begin(), j = i + 1, end = b.end(); j != end; i = j, ++j) {
                    auto temp = weigher.weigh(*i, *j);
                    if(temp > distance) {
                         auto e = *make_edge_index(mesh, *i, *j);
                         mesh.splitEdge(e, distance / temp);
                         auto v = mesh.getEdge(e).vertex;
                         b.insert(j, v);
                         return v;
                    } else distance -= temp;
               }
               if(std::abs(distance) < epsilon) return b.back();
               throw std::runtime_error("Distance must be less than or equal to the path length.");
          };

          if(std::abs((*p).z) < epsilon) {
               while(q != pointsEnd && std::abs((*p).z - (*q).z) < epsilon) ++q;
               auto c0 = (r0) ? Indices(b00, b01) : Indices(r(b01), r(b00));//NOTE: Copy is unnecessary except to keep code clean.  Path is reversed because the lambda function handle processes the path backwards.  The lambda function processes the path backwards because the shortest path calculated in subsequent commands returns the path in reverse order.
               e0 = c0.back();
               e1 = c0.front();
               auto l0 = weigher.weigh(r(c0.cend()), r(c0.cbegin()));
               handle(c0, l0);
          }

          auto wall = getBoundary(t);
          while(p != pointsEnd && std::abs(1.0f - (*p).z) > epsilon) {
               Indices path;

               while(q != pointsEnd && std::abs((*p).z - (*q).z) < epsilon) ++q;
               e0 = position(c2, l2 * (*p).z);
               e1 = position(c1, l1 * (*p).z);
               assert(e0 != e1);
               if(std::distance(p, q) > 2) {
                    assert(shortestPathFinder.getShortestPath(e0, e1, wall.cbegin(), wall.cend(), path));
                    mesh.exsect(r(path.end()), r(path.begin()));
                    path.clear();
                    assert(shortestPathFinder.getShortestPath(e0, e1, wall.cbegin(), wall.cend(), path));
                    wall.insert(wall.end(), r(path.cend()), r(path.cbegin()));
                    auto l = weigher.weigh(r(path.cend()), r(path.cbegin()));
                    handle(path, l);
               } else {
                    samples.push_back(mesh.getVertices()[e0].position);
                    samples.push_back(mesh.getVertices()[e1].position);
                    p += 2;
               }
          }

          while(p != pointsEnd) {
               samples.push_back(mesh.getVertices()[c2.back()].position);
               ++p;
          }

          return samples;
     }

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const TriangleDecomposition& decomposition) {
          stream << decomposition.m_boundaries;
          stream << decomposition.m_indices;
          stream << decomposition.m_neighbors;
          stream << decomposition.m_reverse;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, TriangleDecomposition& decomposition) {
          stream >> decomposition.m_boundaries;
          stream >> decomposition.m_indices;
          stream >> decomposition.m_neighbors;
          stream >> decomposition.m_reverse;
          return stream;
     }

private:
     IndicesArrays m_boundaries;
     Indices m_indices;//every three indices form a triangle border
     Mesh& m_mesh;
     Indices m_neighbors;
     boost::dynamic_bitset<> m_reverse;//whether the border needs to be reversed to be counterclockwise

};//TriangleDecomposition

}//namespace happah

