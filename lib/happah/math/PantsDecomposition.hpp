// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/math/HexagonDecomposition.hpp"
#include "happah/utils/Arrays.hpp"
#include "happah/utils/IteratorJoiner.hpp"
#include "happah/weighers/EdgeLengthWeigher.hpp"
#include "happah/weighers/HoleyWallWeigher.hpp"
#include "happah/weighers/HoleyWallsWeigher.hpp"

namespace happah {

template<class Mesh>
class PantsDecomposition {
     using Indices = std::vector<hpuint>;
     using IndicesArrays = Arrays<hpuint>;
     using Boundary = IndicesArrays::const_iterator::value_type;

     /**********************************************************************************
      * Waist is (0); legs are (1), (2).  Number in parenthesis is offset in boundaries
      * array.  Offset next to line is offset in cuts array.  Arrow points in the 
      * direction that the cut was calculated.  Be aware that the shortest path finder
      * returns the path in the opposite direction.
      *
      *       (2)
      *      /   <
      *   2 /     \ 0
      *    <       \
      *  (0)<------(1)
      *        1
      * 
      * The part of the pants reached by starting at the first edge of the first cut
      * rotating counterclockwise is the first hexagon in the neighbors array.
      *
      * The hexagon's neighbors are arranged in the following fashion:
      *  [ cut 0, leg 2, cut 2, waist 0, cut 1, leg 1 ] at the top and
      *  [ cut 0, leg 1, cut 1, waist 0, cut 2, leg 2 ] at the bottom.
      *
      * 2i, 2i+1 are the hexagons of the ith pants.
      **********************************************************************************/
     class HexagonDecomposer {
     public:
          static HexagonDecomposition<Mesh> decompose(const PantsDecomposition& decomposition) {
               HexagonDecomposer decomposer(decomposition);
               return decomposer.decompose();
          }

     private:
          IndicesArrays m_boundaries;
          const PantsDecomposition& m_decomposition;
          Indices m_indices;//every six indices form a hexagon border
          Mesh& m_mesh;
          Indices m_neighbors;
          boost::dynamic_bitset<> m_reverse;

          HexagonDecomposer(const PantsDecomposition& decomposition)
               : m_decomposition(decomposition), m_mesh(decomposition.m_mesh) {
               hpuint nPants = decomposition.getNumberOfPants();
               hpuint nCuts = 3 * nPants;
               m_boundaries.reserve((decomposition.m_boundaries.size() << 1) + nCuts, decomposition.m_boundaries.data().size());
               hpuint nHexagons = 2 * nPants;
               hpuint nNeighbors = 6 * nHexagons;
               m_neighbors.reserve(nNeighbors);
               m_indices.reserve(nNeighbors);
               m_reverse.resize(nNeighbors, false);
               for(hpuint p = 0, end = nPants; p != nPants; ++p) {//TODO: iterator that returns tuples with all infos
                    if(m_decomposition.isHandle(p)) cutHandlesPants(p);
                    else cutBasePatchsPants(p);

                    const hpuint h0 = p << 1;
                    const hpuint h1 = h0 + 1;
                    m_reverse[6 * h0] = true;//cuts
                    m_reverse[6 * h0 + 2] = true;
                    m_reverse[6 * h1 + 2] = true;
                    m_reverse[6 * h0 + 3] = m_decomposition.m_reverse[3 * p];//waist
                    m_reverse[6 * h1 + 3] = m_decomposition.m_reverse[3 * p];
                    m_reverse[6 * h0 + 5] = m_decomposition.m_reverse[3 * p + 1];//leg1
                    m_reverse[6 * h1 + 1] = m_decomposition.m_reverse[3 * p + 1];
                    m_reverse[6 * h0 + 1] = m_decomposition.m_reverse[3 * p + 2];//leg2
                    m_reverse[6 * h1 + 5] = m_decomposition.m_reverse[3 * p + 2];
               }
          }

          static std::pair<hpuint, hpuint> getEnds(IndicesArrays::const_iterator cuts, hpuint boundary) {
               if(boundary == 0) return getEnds0(cuts);
               else if(boundary == 1) return getEnds1(cuts);
               else {
                    assert(boundary == 2);
                    return getEnds2(cuts);
               }
          }

          static std::pair<hpuint, hpuint> getEnds0(IndicesArrays::const_iterator cuts) { return std::make_pair(*((cuts+1).begin()), *((cuts+2).begin())); }

          static std::pair<hpuint, hpuint> getEnds1(IndicesArrays::const_iterator cuts) { return std::make_pair(*(cuts.end()-1), *((cuts+1).end()-1)); }

          static std::pair<hpuint, hpuint> getEnds2(IndicesArrays::const_iterator cuts) { return std::make_pair(*((cuts+2).end()-1), *(cuts.begin())); }

          //@return: true if neighbor is right side up and false if it is flipped
          static bool getOrientation(IndicesArrays::const_iterator cuts, hpuint boundary, hpuint e0, hpuint e1) {
               hpuint f0, f1;
               std::tie(f0, f1) = getEnds(cuts, boundary);
               if(f0 == e0 && f1 == e1) return false;
               else {
                    assert(f0 == e1 && f1 == e0);
                    return true;
               }
          }

          std::tuple<Indices, Indices, Indices> cut(hpuint p0, hpuint p1, hpuint p2, hpuint e1[2], hpuint e2[2], Indices& wall) {
               Indices temp2 = doCut(p1, e1, e1 + 2, e2, e2 + 2, wall);
               wall.insert(wall.end(), temp2.begin(), temp2.end());
               Indices temp0, temp1;
               std::tie(temp0, temp1) = doCut(p0, *getOther(temp2.back(), e1), *getOther(temp2.front(), e2), wall);
               return std::make_tuple(temp0, temp1, temp2);
          }

          std::tuple<Indices, Indices, Indices> cut(hpuint p0, hpuint p1, hpuint p2, hpuint e2[2], Indices& wall) {
               auto boundary0 = m_decomposition.getBoundary(p0);
               auto boundary1 = m_decomposition.getBoundary(p1);
               auto weigher = make_holey_wall_weigher(p0, boundary0.first, boundary0.second);
               auto shortestPathFinder = make_shortest_path_finder(m_mesh, weigher);

               Indices temp0, temp1, temp2;

               if(!shortestPathFinder.getShortestPath(boundary0.first, boundary0.second, boundary1.first, boundary1.second, wall.begin(), wall.end(), temp0)) std::cerr << "Failed to cut from boundary 1 to boundary 0.\n";
               if(!shortestPathFinder.getShortestPath(boundary0.first, boundary0.second, e2, e2 + 2, wall.begin(), wall.end(), temp1)) std::cerr << "Failed to cut from boundary 1 to end 2.\n";

               {
                    auto weigher = make_edge_length_weigher(m_mesh);
                    auto l0 = weigher.weigh(temp0.begin(), temp0.end());
                    auto l1 = weigher.weigh(temp1.begin(), temp1.end());
                    if(l0 < l1) {
                         auto f0 = ShortestPathFinderUtils::getFarthestPoint(weigher, temp0.back(), boundary0.first, boundary0.second);
                         auto f1 = ShortestPathFinderUtils::getFarthestPoint(weigher, temp0.front(), boundary1.first, boundary1.second);
                         wall.insert(wall.end(), temp0.begin(), temp0.end());
                         temp1 = doCut(p0, f0, f0 + 1, e2, e2 + 2, wall);
                         wall.insert(wall.end(), temp1.begin(), temp1.end());
                         auto f2 = getOther(temp1.front(), e2);
                         temp2 = doCut(p1, f1, f1 + 1, f2, f2 + 1, wall);
                    } else {
                         auto f0 = ShortestPathFinderUtils::getFarthestPoint(weigher, temp1.back(), boundary0.first, boundary0.second);
                         wall.insert(wall.end(), temp1.begin(), temp1.end());
                         std::tie(temp0, temp2) = doCut(p1, *f0, *getOther(temp1.front(), e2), wall);
                         std::reverse(temp0.begin(), temp0.end());
                    }
               }

               return std::make_tuple(temp0, temp1, temp2);
          }

          void cutBasePatchsPants(hpuint p) {
               hpuint n0, n1, n2;
               std::tie(n0, n1, n2) = m_decomposition.getNeighbors(p);
               hpuint p0 = 3 * p;//TODO: rename to b0
               hpuint p1 = 3 * p + 1;
               hpuint p2 = 3 * p + 2;
               auto b0 = m_decomposition.m_indices[p0];//TODO: rename to i0
               auto b1 = m_decomposition.m_indices[p1];
               auto b2 = m_decomposition.m_indices[p2];
               auto boundary0 = m_decomposition.m_boundaries[b0];
               auto boundary1 = m_decomposition.m_boundaries[b1];
               auto boundary2 = m_decomposition.m_boundaries[b2];
               bool d0 = n0 < p;
               bool d1 = n1 < p;
               bool d2 = n2 < p;
               hpuint e0[] = { UNULL, UNULL };
               hpuint e1[] = { UNULL, UNULL };
               hpuint e2[] = { UNULL, UNULL };
               auto i0 = m_decomposition.getBoundaryOffset(n0, b0);//TODO: rename to o0
               auto i1 = m_decomposition.getBoundaryOffset(n1, b1);
               auto i2 = m_decomposition.getBoundaryOffset(n2, b2);
               if(d0) {
                    std::tie(e0[0], e0[1]) = getEnds(n0, i0);
                    assert(e0[0] != e0[1]);
               }
               if(d1) {
                    std::tie(e1[0], e1[1]) = getEnds(n1, i1);
                    assert(e1[0] != e1[1]);
               }
               if(d2) {
                    std::tie(e2[0], e2[1]) = getEnds(n2, i2);
                    assert(e2[0] != e2[1]);
               }

               Indices wall;
               wall.insert(wall.end(), boundary0.first, boundary0.second);
               wall.insert(wall.end(), boundary1.first, boundary1.second);
               wall.insert(wall.end(), boundary2.first, boundary2.second);

               Indices cut0, cut1, cut2;
               if(d0 && d1 && d2) {
                    auto weigher = make_holey_wall_weigher(p1, boundary1.first, boundary1.second);
                    auto shortestPathFinder = make_shortest_path_finder(m_mesh, weigher);
     
                    Indices temp0, temp1;

                    if(!shortestPathFinder.getShortestPath(e1, e1 + 2, e2, e2 + 2, wall.begin(), wall.end(), temp0)) std::cerr << "Failed to cut from boundary to end 0.\n";
                    if(!shortestPathFinder.getShortestPath(e1, e1 + 2, e0, e0 + 2, wall.begin(), wall.end(), temp1)) std::cerr << "Failed to cut from boundary to end 1.\n";

                    {
                         auto weigher = make_edge_length_weigher(m_mesh);
                         auto l0 = weigher.weigh(temp0.begin(), temp0.end());
                         auto l1 = weigher.weigh(temp1.begin(), temp1.end());
                         if(l0 < l1) {
                              auto f1 = getOther(temp0.back(), e1);
                              auto f2 = getOther(temp0.front(), e2);
                              wall.insert(wall.end(), temp0.begin(), temp0.end());
                              if(!shortestPathFinder.getShortestPath(f1, f1 + 1, e0, e0 + 2, wall.begin(), wall.end(), cut1)) std::cerr << "Failed to find cut 1.\n";
                              auto f0 = getOther(cut1.front(), e0);
                              wall.insert(wall.end(), cut1.begin(), cut1.end());
                              cut2 = doCut(p2, f2, f2 + 1, f0, f0 + 1, wall);
                              m_boundaries.push_back(temp0.begin(), temp0.end());
                              m_boundaries.push_back(cut1.begin(), cut1.end());
                              m_boundaries.push_back(cut2.begin(), cut2.end());
                         } else {
                              auto f0 = getOther(temp1.front(), e0);
                              auto f1 = getOther(temp1.back(), e1);
                              wall.insert(wall.end(), temp1.begin(), temp1.end());
                              if(!shortestPathFinder.getShortestPath(f1, f1 + 1, e2, e2 + 2, wall.begin(), wall.end(), cut0)) std::cerr << "Failed to find cut 1.\n";
                              auto f2 = getOther(cut0.front(), e2);
                              wall.insert(wall.end(), cut0.begin(), cut0.end());
                              cut2 = doCut(p2, f2, f2 + 1, f0, f0 + 1, wall);
                              m_boundaries.push_back(cut0.begin(), cut0.end());
                              m_boundaries.push_back(temp1.begin(), temp1.end());
                              m_boundaries.push_back(cut2.begin(), cut2.end());
                         }
                    }
               } else if(!d0 && d1 && d2) {
                    std::tie(cut1, cut2, cut0) = cut(p0, p1, p2, e1, e2, wall);
                    m_boundaries.push_back(cut0.begin(), cut0.end());
                    m_boundaries.push_back(cut1.rbegin(), cut1.rend());
                    m_boundaries.push_back(cut2.rbegin(), cut2.rend());
               } else if(d0 && !d1 && d2) {
                    std::tie(cut1, cut0, cut2) = cut(p1, p0, p2, e0, e2, wall);
                    m_boundaries.push_back(cut0.begin(), cut0.end());
                    m_boundaries.push_back(cut1.begin(), cut1.end());
                    m_boundaries.push_back(cut2.rbegin(), cut2.rend());
               } else if(d0 && d1 && !d2) {
                    std::tie(cut2, cut0, cut1) = cut(p2, p0, p1, e0, e1, wall);
                    m_boundaries.push_back(cut0.rbegin(), cut0.rend());
                    m_boundaries.push_back(cut1.rbegin(), cut1.rend());
                    m_boundaries.push_back(cut2.begin(), cut2.end());
               } else if(!d0 && !d1 && d2) {
                    std::tie(cut1, cut2, cut0) = cut(p0, p1, p2, e2, wall);
                    m_boundaries.push_back(cut0.begin(), cut0.end());
                    m_boundaries.push_back(cut1.rbegin(), cut1.rend());
                    m_boundaries.push_back(cut2.rbegin(), cut2.rend());
               } else if(!d0 && d1 && !d2) {
                    std::tie(cut2, cut1, cut0) = cut(p0, p2, p1, e1, wall);
                    m_boundaries.push_back(cut0.rbegin(), cut0.rend());
                    m_boundaries.push_back(cut1.rbegin(), cut1.rend());
                    m_boundaries.push_back(cut2.rbegin(), cut2.rend());
               } else if(d0 && !d1 && !d2) {
                    std::tie(cut0, cut1, cut2) = cut(p1, p2, p0, e0, wall);
                    m_boundaries.push_back(cut0.begin(), cut0.end());
                    m_boundaries.push_back(cut1.begin(), cut1.end());
                    m_boundaries.push_back(cut2.begin(), cut2.end());
               } else {
                    assert(!d0 && !d1 && !d2);
                    auto weigher = make_holey_wall_weigher(p1, boundary1.first, boundary1.second);
                    auto shortestPathFinder = make_shortest_path_finder(m_mesh, weigher);

                    if(!shortestPathFinder.getShortestPath(boundary1.first, boundary1.second, boundary0.first, boundary0.second, wall.begin(), wall.end(), cut1)) std::cerr << "Failed to cut from boundary 1 to boundary 0.\n";
                    if(!shortestPathFinder.getShortestPath(boundary1.first, boundary1.second, boundary2.first, boundary2.second, wall.begin(), wall.end(), cut0)) std::cerr << "Failed to cut from boundary 1 to boundary 2.\n";

                    {
                         auto weigher = make_edge_length_weigher(m_mesh);
                         auto l0 = weigher.weigh(cut0.begin(), cut0.end());
                         auto l1 = weigher.weigh(cut1.begin(), cut1.end());
                         if(l0 < l1) {
                              auto f1 = ShortestPathFinderUtils::getFarthestPoint(weigher, cut0.back(), boundary1.first, boundary1.second);
                              auto f2 = ShortestPathFinderUtils::getFarthestPoint(weigher, cut0.front(), boundary2.first, boundary2.second);
                              wall.insert(wall.end(), cut0.begin(), cut0.end());
                              std::tie(cut1, cut2) = doCut(p0, *f1, *f2, wall);
                              m_boundaries.push_back(cut0.begin(), cut0.end());
                              m_boundaries.push_back(cut1.rbegin(), cut1.rend());
                              m_boundaries.push_back(cut2.rbegin(), cut2.rend());
                         } else {
                              auto f1 = ShortestPathFinderUtils::getFarthestPoint(weigher, cut1.back(), boundary1.first, boundary1.second);
                              auto f0 = ShortestPathFinderUtils::getFarthestPoint(weigher, cut1.front(), boundary0.first, boundary0.second);
                              wall.insert(wall.end(), cut1.begin(), cut1.end());
                              std::tie(cut2, cut0) = doCut(p2, *f0, *f1, wall);
                              m_boundaries.push_back(cut0.rbegin(), cut0.rend());
                              m_boundaries.push_back(cut1.begin(), cut1.end());
                              m_boundaries.push_back(cut2.begin(), cut2.end());
                         }
                    }
               }

               const hpuint h0 = p << 1;//top hexagon
               const hpuint h1 = h0 + 1;//bottom hexagon
               hpuint neighbors[] = { h1, UNULL, h1, UNULL, h1, UNULL, h0, UNULL, h0, UNULL, h0, UNULL };
               hpuint nBoundaries = m_boundaries.size();
               auto cuts = m_boundaries.cend() - 3;
               hpuint c0 = nBoundaries - 3;
               hpuint c1 = nBoundaries - 2;
               hpuint c2 = nBoundaries - 1;
               hpuint t00, t10, t20, t01, t11, t21;
               hpuint e00, e10, e01, e11, e02, e12;
               std::tie(e00, e10) = getEnds0(cuts);
               std::tie(e01, e11) = getEnds1(cuts);
               std::tie(e02, e12) = getEnds2(cuts);
               if(d0) {
                    if(getOrientation(getCuts(n0), i0, e00, e10)) {
                         neighbors[3] = (n0 << 1);
                         neighbors[9] = (n0 << 1) + 1;
                         setNeighbors(n0, i0, h0, h1);
                         t00 = getTopBoundaryIndex(n0, i0);
                         t01 = getBottomBoundaryIndex(n0, i0);
                    } else {
                         neighbors[3] = (n0 << 1) + 1;
                         neighbors[9] = (n0 << 1);
                         setNeighbors(n0, i0, h1, h0);
                         t00 = getBottomBoundaryIndex(n0, i0);
                         t01 = getTopBoundaryIndex(n0, i0);
                    }
               } else {
                    if(m_decomposition.m_reverse[3 * p]) {
                         t01 = nBoundaries;
                         t00 = t01 + 1;
                    } else {
                         t00 = nBoundaries;
                         t01 = t00 + 1;
                    }
                    split(boundary0, e00, e10);
                    nBoundaries += 2;
               }
               if(d1) {
                    if(getOrientation(getCuts(n1), i1, e01, e11)) {
                         neighbors[5] = (n1 << 1);
                         neighbors[7] = (n1 << 1) + 1;
                         setNeighbors(n1, i1, h0, h1);
                         t10 = getTopBoundaryIndex(n1, i1);
                         t11 = getBottomBoundaryIndex(n1, i1);
                    } else {
                         neighbors[5] = (n1 << 1) + 1;
                         neighbors[7] = (n1 << 1);
                         setNeighbors(n1, i1, h1, h0);
                         t10 = getBottomBoundaryIndex(n1, i1);
                         t11 = getTopBoundaryIndex(n1, i1);
                    }
               } else {
                    if(m_decomposition.m_reverse[3 * p + 1]) {
                         t11 = nBoundaries;
                         t10 = t11 + 1;
                    } else {
                         t10 = nBoundaries;
                         t11 = t10 + 1;
                    }
                    split(boundary1, e01, e11);
                    nBoundaries += 2;
               }
               if(d2) {
                    if(getOrientation(getCuts(n2), i2, e02, e12)) {
                         neighbors[1] = (n2 << 1);
                         neighbors[11] = (n2 << 1) + 1;
                         setNeighbors(n2, i2, h0, h1);
                         t20 = getTopBoundaryIndex(n2, i2);
                         t21 = getBottomBoundaryIndex(n2, i2);
                    } else {
                         neighbors[1] = (n2 << 1) + 1;
                         neighbors[11] = (n2 << 1);
                         setNeighbors(n2, i2, h1, h0);
                         t20 = getBottomBoundaryIndex(n2, i2);
                         t21 = getTopBoundaryIndex(n2, i2);
                    }
               } else {
                    if(m_decomposition.m_reverse[3 * p + 2]) {
                         t21 = nBoundaries;
                         t20 = t21 + 1;
                    } else {
                         t20 = nBoundaries;
                         t21 = t20 + 1;
                    }
                    split(boundary2, e02, e12);
               }
               m_neighbors.insert(m_neighbors.end(), neighbors, neighbors + 12);
               hpuint indices[] = { c0, t20, c2, t00, c1, t10, c0, t11, c1, t01, c2, t21 };
               m_indices.insert(m_indices.end(), indices, indices + 12);
          }

          void cutHandlesPants(hpuint p) {
               using Iterator = IndicesArrays::const_iterator::iterator;

               auto boundary0 = m_decomposition.m_boundaries[m_decomposition.m_indices[3 * p]];
               auto boundary12 = m_decomposition.m_boundaries[m_decomposition.m_indices[3 * p + 1]];

               {
                    //NOTE: Finds inner cut between pants' legs.
                    auto weigher = happah::make_holey_wall_weigher(m_mesh, boundary12.first, boundary12.second);
                    weigher.removeVertices(boundary0.first, boundary0.second);
                    auto shortestPathFinder = make_shortest_path_finder(m_mesh, weigher);
                    auto cut0 = shortestPathFinder.template getShortestLoop<Iterator, false>(boundary12.first, boundary12.second);
                    assert(cut0[0] == cut0.back());//NOTE: The ends of the inner cut must meet.
                    m_boundaries.push_back(cut0.cbegin(), cut0.cend());
               }

               Iterator farthestPoint;
               {
                    auto weigher = make_edge_length_weigher(m_mesh);
                    farthestPoint = ShortestPathFinderUtils::getFarthestPoint(weigher, *((--m_boundaries.end()).begin()), boundary12.first, boundary12.second);
               }

               auto n0 = m_decomposition.m_neighbors[3 * p];
               auto b0 = m_decomposition.m_indices[3 * p];
               bool d0 = n0 < p;
               hpuint e0[2];
               hpuint i0 = m_decomposition.getBoundaryOffset(n0, b0);

               {
                    auto weigher = happah::make_holey_wall_weigher(m_mesh, boundary12.first, boundary12.second);
                    weigher.template puncture<false>(farthestPoint);
                    auto shortestPathFinder = make_shortest_path_finder(m_mesh, weigher);

                    if(d0) {
                         std::tie(e0[0], e0[1]) = getEnds(n0, i0);
                         assert(e0[0] != e0[1]);

                         Indices path0, path1;
                         if(!shortestPathFinder.getShortestPath(*farthestPoint, e0, e0 + 2, (--m_boundaries.end()).begin(), m_boundaries.end().begin(), path0)) std::cerr << "ERROR: Failed to find handle's pants' first leg to waist cut (0).\n";
                         weigher.plug(farthestPoint);
                         weigher.template puncture<true>(farthestPoint);
                         if(!shortestPathFinder.getShortestPath(*farthestPoint, path0[0], (--m_boundaries.end()).begin(), m_boundaries.end().begin(), path1)) std::cerr << "ERROR: Failed to find handle's pants' first leg to waist cut (1).\n";

                         auto l0 = weigher.weigh(path0.cbegin(), path0.cend());
                         auto l1 = weigher.weigh(path1.cbegin(), path1.cend());
                         auto target = (path0[0] == e0[0]) ? e0[1] : e0[0];
                         if(l1 > l0) {
                              m_boundaries.push_back(path0.cbegin(), path0.cend());
                              if(!shortestPathFinder.getShortestPath(*farthestPoint, target, (m_boundaries.end()-2).begin(), m_boundaries.end().begin(), IndicesArrays::ArrayAppender(m_boundaries))) std::cerr << "ERROR: Failed to find handle's pants' second leg to waist cut.\n";
                         } else {
                              weigher.plug(farthestPoint);
                              weigher.template puncture<false>(farthestPoint);
                              if(!shortestPathFinder.getShortestPath(*farthestPoint, target, (m_boundaries.end()-2).begin(), m_boundaries.end().begin(), IndicesArrays::ArrayAppender(m_boundaries))) std::cerr << "ERROR: Failed to find handle's pants' second leg to waist cut.\n";
                              m_boundaries.push_back(path1.cbegin(), path1.cend());
                         }
                    } else {
                         if(!shortestPathFinder.getShortestPath(*farthestPoint, boundary0.first, boundary0.second, (--m_boundaries.end()).begin(), m_boundaries.end().begin(), IndicesArrays::ArrayAppender(m_boundaries))) std::cerr << "ERROR: Failed to find handle's pants' first leg to waist cut.\n";
                         weigher.plug(farthestPoint);
                         weigher.template puncture<true>(farthestPoint);
                         if(!shortestPathFinder.getShortestPath(*farthestPoint, boundary0.first, boundary0.second, (m_boundaries.end()-2).begin(), m_boundaries.end().begin(), IndicesArrays::ArrayAppender(m_boundaries))) std::cerr << "ERROR: Failed to find handle's pants' second leg to waist cut.\n";
                    }
               }

               const hpuint h0 = p << 1;
               const hpuint h1 = h0 + 1;
               hpuint neighbors[] = { h1, h0, h1, UNULL, h1, h0, h0, h1, h0, UNULL, h0, h1 };
               hpuint nBoundaries = m_boundaries.size();
               auto cuts = m_boundaries.cend() - 3;
               hpuint c0 = nBoundaries - 3;
               hpuint c1 = nBoundaries - 2;
               hpuint c2 = nBoundaries - 1;
               hpuint t00, t01, t12;
               hpuint e00, e10, e01, e11;
               std::tie(e00, e10) = getEnds0(cuts);
               std::tie(e01, e11) = getEnds1(cuts);
               if(d0) {
                    if(getOrientation(getCuts(n0), i0, e00, e10)) {
                         neighbors[3] = (n0 << 1);
                         neighbors[9] = (n0 << 1) + 1;
                         setNeighbors(n0, i0, h0, h1);
                         t00 = getTopBoundaryIndex(n0, i0);
                         t01 = getBottomBoundaryIndex(n0, i0);
                    } else {
                         neighbors[3] = (n0 << 1) + 1;
                         neighbors[9] = (n0 << 1);
                         setNeighbors(n0, i0, h1, h0);
                         t00 = getBottomBoundaryIndex(n0, i0);
                         t01 = getTopBoundaryIndex(n0, i0);
                    }
               } else {
                    if(m_decomposition.m_reverse[3 * p]) {
                         t01 = nBoundaries;
                         t00 = t01 + 1;
                    } else {
                         t00 = nBoundaries;
                         t01 = t00 + 1;
                    }
                    split(boundary0, e00, e10);
                    nBoundaries += 2;
               }
               t12 = nBoundaries;
               split(boundary12, e01, e11);
               m_neighbors.insert(m_neighbors.end(), neighbors, neighbors + 12);
               hpuint indices[] = { c0, t12, c2, t00, c1, t12, c0, t12 + 1, c1, t01, c2, t12 + 1 };
               m_indices.insert(m_indices.end(), indices, indices + 12);
          }

          //TODO: need to think about where there are dangerous edges in the decomposition into a hexagon mesh if any
          HexagonDecomposition<Mesh> decompose() { return { m_mesh, std::move(m_boundaries), std::move(m_indices), std::move(m_neighbors), std::move(m_reverse) }; }

          template<class Iterator0, class Iterator1> 
          Indices doCut(hpuint p0, Iterator0 begin0, Iterator0 end0, Iterator1 begin1, Iterator1 end1, Indices& wall) {
               auto weigher = make_holey_wall_weigher(p0);
               auto shortestPathFinder = make_shortest_path_finder(m_mesh, weigher);
               Indices temp;
               if(!shortestPathFinder.getShortestPath(begin0, end0, begin1, end1, wall.begin(), wall.end(), temp)) std::cerr << "Failed to cut from end 0 to end 1.\n";
               return temp;
          }

          std::tuple<Indices, Indices> doCut(hpuint p0, hpuint e1, hpuint e2, Indices& wall) {
               auto boundary0 = m_decomposition.getBoundary(p0);
               auto weigher = make_holey_wall_weigher(p0, boundary0.first, boundary0.second);
               auto shortestPathFinder = make_shortest_path_finder(m_mesh, weigher);

               Indices temp0, temp1;

               if(!shortestPathFinder.getShortestPath(boundary0.first, boundary0.second, &e1, &e1 + 1, wall.begin(), wall.end(), temp0)) std::cerr << "Failed to cut from boundary to end 0.\n";
               if(!shortestPathFinder.getShortestPath(boundary0.first, boundary0.second, &e2, &e2 + 1, wall.begin(), wall.end(), temp1)) std::cerr << "Failed to cut from boundary to end 1.\n";

               {
                    auto weigher = make_edge_length_weigher(m_mesh);
                    auto l0 = weigher.weigh(temp0.begin(), temp0.end());
                    auto l1 = weigher.weigh(temp1.begin(), temp1.end());
                    if(l0 < l1) {
                         wall.insert(wall.end(), temp0.begin(), temp0.end());
                         temp1.clear();
                         auto f0 = ShortestPathFinderUtils::getFarthestPoint(weigher, temp0.back(), boundary0.first, boundary0.second);
                         if(!shortestPathFinder.getShortestPath(*f0, e2, wall.begin(), wall.end(), temp1)) std::cerr << "Failed to cut from boundary to end 2.\n";
                    } else {
                         wall.insert(wall.end(), temp1.begin(), temp1.end());
                         temp0.clear();
                         auto f0 = ShortestPathFinderUtils::getFarthestPoint(weigher, temp1.back(), boundary0.first, boundary0.second);
                         if(!shortestPathFinder.getShortestPath(*f0, e1, wall.begin(), wall.end(), temp0)) std::cerr << "Failed to cut from boundary to end 1.\n";
                    }
               }

               return std::make_tuple(temp0, temp1);
          }

          hpuint getBottomBoundaryIndex(hpuint p, hpuint i) {
               if(i == 0) return m_indices[12 * p + 9];
               else if(i == 1) return m_indices[12 * p + 7];
               else {
                    assert(i == 2);
                    return m_indices[12 * p + 11];
               }
          }

          IndicesArrays::const_iterator getCuts(hpuint p) const { return m_boundaries.cbegin() + m_indices[12 * p]; }

          std::pair<hpuint, hpuint> getEnds(hpuint p, hpuint boundary) { return getEnds(getCuts(p), boundary); }

          auto getOther(hpuint e, hpuint ei[2]) {
               if(e == ei[0]) return ei + 1;
               else {
                    assert(e == ei[1]);
                    return ei;
               }
          }

          hpuint getTopBoundaryIndex(hpuint p, hpuint i) {
               if(i == 0) return m_indices[12 * p + 3];
               else if(i == 1) return m_indices[12 * p + 5];
               else {
                    assert(i == 2);
                    return m_indices[12 * p + 1];
               }
          }

          auto make_holey_wall_weigher(hpuint b) {
               auto boundary = m_decomposition.getBoundary(b);
               return make_holey_wall_weigher(b, boundary.first, boundary.second);
          }

          template<class Iterator>
          auto make_holey_wall_weigher(hpuint b, Iterator begin, Iterator end) {
               auto weigher = happah::make_holey_wall_weigher(m_mesh, begin, end);
               if(m_decomposition.m_reverse[b]) while(begin != end) weigher.template puncture<true>(begin++);
               else while(begin != end) weigher.template puncture<false>(begin++);
               return weigher;
          }
          
          void setNeighbors(hpuint n, hpuint boundary, hpuint top, hpuint bottom) {
               hpuint i = 12 * n;
               if(boundary == 0) {
                    m_neighbors[i + 3] = top;
                    m_neighbors[i + 9] = bottom;
               } else if(boundary == 1) {
                    m_neighbors[i + 5] = top;
                    m_neighbors[i + 7] = bottom;
               } else {
                    assert(boundary == 2);
                    m_neighbors[i + 1] = top;
                    m_neighbors[i + 11] = bottom;
               }
          }

          void split(const Boundary& boundary, hpuint e0, hpuint e1) {
               auto begin = boundary.first;
               auto end = boundary.second;

               while(begin != end && *begin != e1) ++begin;
               assert(*begin == e1);

               {
                    IndicesArrays::ArrayAppender boundaries(m_boundaries);
                    while(*begin != e0) {
                         boundaries.push_back(*begin);
                         ++begin;
                         if(begin == end) begin = boundary.first + 1; 
                    }
                    boundaries.push_back(*begin);
               }

               assert(*begin == e0);

               {
                    IndicesArrays::ArrayAppender boundaries(m_boundaries);
                    while(*begin != e1) {
                         boundaries.push_back(*begin);
                         ++begin;
                         if(begin == end) begin = boundary.first + 1; 
                    }
                    boundaries.push_back(*begin);
               }

               assert(*begin == e1);
          }

     };//HexagonDecomposer

public:
     //NOTE: The contract is that any changes to the mesh leave the given boundaries intact.
     PantsDecomposition(Mesh& mesh, IndicesArrays boundaries, Indices indices, Indices neighbors, boost::dynamic_bitset<> reverse)
          : m_boundaries(std::move(boundaries)), m_indices(std::move(indices)), m_mesh(mesh), m_neighbors(std::move(neighbors)), m_reverse(std::move(reverse)) {
          for(auto b = m_boundaries.cbegin(), end = m_boundaries.cend(); b != end; ++b) {
               auto boundary = *b;
               assert(*boundary.first == *(boundary.second - 1));
          }
     }

     auto getBoundary(hpuint b) const { return m_boundaries[m_indices[b]]; }

     auto getBoundary(hpuint p, hpuint o) const { return m_boundaries[m_indices[3 * p + o]]; }

     std::tuple<hpuint, hpuint, hpuint> getNeighbors(hpuint p) const {
          auto n = m_neighbors.cbegin() + (3 * p);
          return std::make_tuple(*n, *(n + 1), *(n + 2));
     }

     hpuint getNumberOfPants() const { return m_indices.size() / 3; }

     bool isHandle(hpuint p) const {
          auto i = m_indices.cbegin() + (3 * p);
          return  *(i + 1) == *(i + 2) || *i == *(i + 1) || *i == *(i + 2);//TODO: current implementation assumes common boundary is at i+1 and i+2
     }

     bool isNeighbor(hpuint p, hpuint n) const {
          auto i = m_neighbors.cbegin() + 3 * p;
          if(*i == n) return true;
          else if(*(++i) == n) return true;
          else if(*(++i) == n) return true;
          else return false;
     }

     void move(hpuint p) {//simple move
          assert(isHandle(p));
          assert(false);//TODO
     }

     //association = true if leg 1 to leg 1; false if leg 1 of p0 to leg 2 0f p1
     void move(const hpuint p0, const hpuint i0, bool association) {//associativity move
          using Iterator = IndicesArrays::const_iterator::iterator;

          const auto p1 = m_neighbors[3 * p0 + i0];
          assert(p0 != p1);

          const auto i1 = (i0 == 2) ? 0 : i0 + 1;
          const auto i2 = (i1 == 2) ? 0 : i1 + 1;

          auto b0 = m_indices[3 * p0 + i0];
          const auto j0 = getBoundaryOffset(p1, b0);
          const auto j1 = (j0 == 2) ? 0 : j0 + 1;
          const auto j2 = (j1 == 2) ? 0 : j1 + 1;

          auto path0 = (association) ? connect(p0, i1, i2, p1, j1, j2) : connect(p0, i1, i2, p1, j2, j1);
          auto path1 = (association) ? connect(p0, i2, i1, p1, j2, j1) : connect(p0, i2, i1, p1, j1, j2);
          //TODO: split dangerous edges between path0 and path1

          auto boundary10 = m_boundaries[m_indices[3 * p0 + i1]];
          auto boundary20 = m_boundaries[m_indices[3 * p0 + i2]];
          auto boundary11 = m_boundaries[m_indices[3 * p1 + j1]];
          auto boundary21 = m_boundaries[m_indices[3 * p1 + j2]];

          IteratorJoiner<Iterator> wall;
          wall.push_back(boundary10.first, boundary10.second);
          wall.push_back(boundary20.first, boundary20.second);
          wall.push_back(boundary11.first, boundary11.second);
          wall.push_back(boundary21.first, boundary21.second);

          Indices path2;
          { 
               using Weigher = EdgeLengthWeigher<Mesh>;
               using ShortestPathFinder = ShortestPathFinder<Weigher, Mesh>;

               Weigher weigher(m_mesh);
               ShortestPathFinder shortestPathFinder(m_mesh, weigher);
               if(!shortestPathFinder.getShortestPath(path0.cbegin() + 1, path0.cend() - 1, path1.cbegin() + 1, path1.cend() - 1, wall.begin(), wall.end(), path2)) std::cerr << "ERROR: Failed to connect two paths.\n";
          }
          assert(path2.size() > 2);

          Indices waist;
          {
               using Weigher = HoleyWallWeigher<Mesh, Iterator>;
               using ShortestPathFinder = ShortestPathFinder<Weigher, Mesh>;

               wall.push_back(path0.cbegin(), path0.cend());
               wall.push_back(path1.cbegin(), path1.cend());
               Weigher weigher(m_mesh, path2.cbegin(), path2.cend(), false);
               weigher.removeVertices(wall.begin(), wall.end());
               ShortestPathFinder shortestPathFinder(m_mesh, weigher);
               waist = shortestPathFinder.getShortestLoop(path2.cbegin(), path2.cend());
          }
          assert(waist[0] == waist.back());//NOTE: The ends of the waist must meet.

          auto h0 = isHandle(p0);
          auto h1 = isHandle(p1);

          m_boundaries.erase(m_boundaries.begin() + b0);
          for(auto& i : m_indices) if(i > b0) --i;
          auto bn = m_boundaries.size();
          m_boundaries.push_back(waist.cbegin(), waist.cend());

          auto b2 = m_indices[3 * p0 + i2];
          auto n2 = m_neighbors[3 * p0 + i2];
          auto r2 = m_reverse[3 * p0 + i2];

          //NOTE: The new pants p0 contains path0.
          m_indices[3 * p0 + i0] = bn;
          m_indices[3 * p1 + j0] = bn;
          m_reverse[3 * p0 + i0] = true;
          m_reverse[3 * p1 + j0] = false;
          if(h0) m_neighbors[3 * p0 + i1] = p1;
          if(association) {
               m_indices[3 * p0 + i2] = m_indices[3 * p1 + j1];
               m_indices[3 * p1 + j1] = b2;
               if(h1) m_neighbors[3 * p1 + j2] = p0;
               m_neighbors[3 * p0 + i2] = m_neighbors[3 * p1 + j1];
               m_neighbors[3 * p1 + j1] = n2;
               m_reverse[3 * p0 + i2] = m_reverse[3 * p1 + j1];
               m_reverse[3 * p1 + j1] = r2;
          } else {
               m_indices[3 * p0 + i2] = m_indices[3 * p1 + j2];
               m_indices[3 * p1 + j2] = b2;
               if(h1) m_neighbors[3 * p1 + j1] = p0;
               m_neighbors[3 * p0 + i2] = m_neighbors[3 * p1 + j2];
               m_neighbors[3 * p1 + j2] = n2;
               m_reverse[3 * p0 + i2] = m_reverse[3 * p1 + j2];
               m_reverse[3 * p1 + j2] = r2;
          }
     }

     HexagonDecomposition<Mesh> toHexagonDecomposition() { return HexagonDecomposer::decompose(*this); }

private:
     IndicesArrays m_boundaries;
     Indices m_indices;
     Mesh& m_mesh;
     Indices m_neighbors;
     boost::dynamic_bitset<> m_reverse;//whether the border needs to be reversed to be counterclockwise

     //NOTE: Connect p0's i0 boundary to p1's j0 boundary.
     Indices connect(hpuint p0, hpuint i0, hpuint i1, hpuint p1, hpuint j0, hpuint j1) const {//TODO: maybe lamda function
          using Iterator = IndicesArrays::const_iterator::iterator;
          using Weigher = HoleyWallsWeigher<Mesh, Iterator>;
          using ShortestPathFinder = ShortestPathFinder<Weigher, Mesh>;

          auto i = m_indices.cbegin() + 3 * p0;
          auto boundary00 = m_boundaries[*(i + i0)];
          auto boundary10 = m_boundaries[*(i + i1)];

          auto j = m_indices.cbegin() + 3 * p1;
          auto boundary01 = m_boundaries[*(j + j0)];
          auto boundary11 = m_boundaries[*(j + j1)];

          IteratorJoiner<Iterator> wall;
          wall.push_back(boundary10.first, boundary10.second);
          wall.push_back(boundary11.first, boundary11.second);

          boost::dynamic_bitset<> loops(2, true);
          Weigher weigher(m_mesh, { boundary00.first, boundary00.second, boundary01.first, boundary01.second }, loops);
          if(m_reverse[3 * p0 + i0]) for(auto b = boundary00.first; b != boundary00.second; ++b) weigher.template puncture<true>(b);
          else for(auto b = boundary00.first; b != boundary00.second; ++b) weigher.template puncture<false>(b);
          if(m_reverse[3 * p1 + j0]) for(auto b = boundary01.first; b != boundary01.second; ++b) weigher.template puncture<false>(b);
          else for(auto b = boundary01.first; b != boundary01.second; ++b) weigher.template puncture<true>(b);
          ShortestPathFinder shortestPathFinder(m_mesh, weigher);

          Indices path;
          if(!shortestPathFinder.getShortestPath(boundary00.first, boundary00.second, boundary01.first, boundary01.second, wall.begin(), wall.end(), path)) std::cerr << "ERROR: Failed to connect legs.\n";
          return path;
     }

     hpuint getBoundaryOffset(hpuint p, hpuint b) const {
          auto i = m_indices.cbegin() + (3 * p);
          if(*i == b) return 0;
          else if(*(++i) == b) return 1;
          else {
               assert(*(++i) == b);
               return 2;
          }
     }

};//PantsDecomposition

}//namespace happah

