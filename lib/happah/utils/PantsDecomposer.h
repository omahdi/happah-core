// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/geometries/Mesh.h"
#include "happah/math/PantsDecomposition.h"
#include "happah/utils/Arrays.h"
#include "happah/utils/IteratorJoiner.h"
#include "happah/utils/ShortestPathFinder.h"
#include "happah/weighers/EdgeLengthWeigher.h"
#include "happah/weighers/HoleyWallWeigher.h"

//REFERENCE: [Li2008] Xin Li, Xianfeng Gu, and Hong Qin.  Surface Matching Using Consistent Pants Decomposition. ACM, 2008.

namespace happah {

template<class Mesh>
class PantsDecomposer {
     using Indices = std::vector<hpuint>;
     using IndicesArrays = Arrays<hpuint>;

public:
     static PantsDecomposition<Mesh> decompose(Mesh& mesh) {
          PantsDecomposer decomposer(mesh);
          return decomposer.decompose();
     }

private:
     IndicesArrays m_boundaries;
     Indices m_indices;
     Mesh& m_mesh;
     Indices m_neighbors;
     boost::dynamic_bitset<> m_reverse;//whether the border needs to be reversed to be counterclockwise

     PantsDecomposer(Mesh& mesh) 
          : m_mesh(mesh) {
          using Weigher = EdgeLengthWeigher<Mesh>;
          using ShortestPathFinder = ShortestPathFinder<Weigher, Mesh>;

          hpuint genus = *mesh.getGenus();
          assert(genus > 1);

          auto nBoundaries = 3 * (genus - 1);
          auto nPants = (genus - 1) << 1;

          m_boundaries.reserve(nBoundaries);
          m_indices.reserve(3 * nPants);
          m_neighbors.reserve(3 * nPants);
          m_reverse.resize(3 * nPants, false);

          Weigher weigher(mesh);
          ShortestPathFinder shortestPathFinder(mesh, weigher);

          auto& handleTunnelLoops = *(mesh.getHandleTunnelLoops());

          if(genus == 2) {
               std::vector<hpuint> path;

               auto h0 = handleTunnelLoops.cbegin();
               auto h1 = h0 + 1;
               auto h2 = h1 + 1;
               auto h3 = h2 + 1;

               splitDangerousEdges(h0.begin(), h1.end(), h2.begin(), h3.end());
               if(!shortestPathFinder.getShortestPath(h2.begin(), h3.end(), h0.begin(), h1.end(), path)) std::cerr << "ERROR: Failed to find shortest path connecting two handles.\n";

               auto waist = getShortestLoop(path.begin(), path.end(), h0.begin(), h3.end());
               insert(waist.begin(), waist.end(), 0);
               insert(0u, 0u, 1u);

               return;
          }

          IndicesArrays cap, holes;
          Indices holesIndices;
          holesIndices.reserve(genus);
          for(hpuint i = 1; i <= genus; ++i)
               holesIndices.push_back(i);

          std::cout << "sorting holes\n";
          {
               //NOTE: Makes copy of handles and tunnels with indices in each loop sorted.
               holes.reserve(genus - 1, handleTunnelLoops.data().size());
               auto i = handleTunnelLoops.begin();
               auto end = handleTunnelLoops.end();
               auto holeBegin = (*i).first;
               auto holeEnd = (*(++i)).second;
               cap.push_back(holeBegin, holeEnd);
               ++i;
               while(i != end) {
                    holeBegin = (*i).first;
                    holeEnd = (*(++i)).second;
                    holes.push_back(holeBegin, holeEnd);
                    auto hole = --holes.end();
                    std::sort(hole.begin(), hole.end());
                    ++i;
               }
          }

          std::cout << "creating cap\n";
          std::vector<hpuint> handles;
          handles.reserve(genus);
          handles.push_back(0);
          hpuint limit = genus - 1;
          while(limit > 0) {
               if(limit == 1) splitDangerousEdges(holes.begin().begin(), holes.end().begin(), cap.begin().begin(), cap.end().begin());
               if(!shortestPathFinder.getShortestPath(cap.begin().begin(), cap.end().begin(), holes.begin().begin(), holes.end().begin(), IndicesArrays::ArrayAppender(cap))) std::cerr << "ERROR: Failed to find shortest path connecting hole to cap.\n";
               hpuint target = *((--cap.end()).begin());
               auto holeIndex = holesIndices.begin();
               auto nHoles = holes.size();
               for(auto i = holes.begin(), end = holes.end(); i != end; ++i) {
                    auto hole = *i;
                    if(std::binary_search(hole.first, hole.second, target)) {
                         cap.reserve(1, cap.data().size() + std::distance(hole.first, hole.second));
                         cap.push_back(hole.first, hole.second);
                         holes.erase(i);
                         handles.push_back(*holeIndex);
                         holesIndices.erase(holeIndex);
                         break;
                    }
                    ++holeIndex;
               }
               assert(holes.size() == (nHoles - 1));
               --limit;
          }

          auto handleIndex = handles.end() - 1;

          std::cout << "finding first waist" << std::endl;
          try {
               //NOTE: Finds first waist.
               auto path = --(--cap.end());
               //NOTE: If path length is less than three, then path is a dangerous edge, which should theoretically not happen because of the above call to split dangerous edges.
               auto waist = getShortestLoop(path.begin(), path.end(), cap.begin().begin(), cap.end().begin());
               insert(waist.begin(), waist.end(), *handleIndex);//waist is counterclockwise
               --handleIndex;
               cap.erase(path, cap.end());
               cap.push_back(waist.begin(), waist.end());
               --genus;
          } catch(...) { std::cerr << "ERROR: Failed to find first waist.\n"; }

          std::cout << "finding middle waists" << std::endl;
          limit = genus;
          while(genus > 1) {
               //NOTE: Finds middle waists.
               std::vector<hpuint> hole;
               {
                    //NOTE: Removes next hole to process.
                    auto temp = cap.begin() + ((genus - 1) << 1);
                    hole.insert(hole.end(), temp.begin(), temp.end());
                    cap.erase(temp - 1, temp + 1);
               }

               {
                    //NOTE: Connects waists to rest of cap.
                    auto temp = cap.begin() + (((genus - 2) << 1) + 1);
                    if(!shortestPathFinder.getShortestPath(temp.begin(), cap.end().begin(), cap.begin().begin(), temp.begin(), hole.begin(), hole.end(), IndicesArrays::ArrayAppender(cap))) std::cerr << "ERROR: Failed to connect waists to rest of cap.\n";
               }

               //NOTE: Connects hole to cap and finds waist.
               splitDangerousEdges(hole.begin(), hole.end(), cap.begin().begin(), cap.end().begin());
               if(shortestPathFinder.getShortestPath(cap.begin().begin(), cap.end().begin(), hole.begin(), hole.end(), IndicesArrays::ArrayAppender(cap))) {
                    cap.push_back(hole.begin(), hole.end());
                    auto path = --(--cap.end());
                    //NOTE: If path length is less than three, then path is a dangerous edge, which should theoretically not happen because of the above call to split dangerous edges.
                    auto waist = getShortestLoop(path.begin(), path.end(), cap.begin().begin(), cap.end().begin());
                    insert(waist.begin(), waist.end(), *handleIndex);//waist is counterclockwise
                    --handleIndex;
                    cap.erase(--path, cap.end());
                   
                    //NOTE: Connects new waist to other waists.
                    if(!shortestPathFinder.getShortestPath(waist.begin(), waist.end(), (cap.begin() + (((genus - 2) << 1) + 1)).begin(), cap.end().begin(), IndicesArrays::ArrayAppender(cap))) std::cerr << "ERROR: Failed to connect waist to the other waists.\n";
                    cap.push_back(waist.begin(), waist.end());
               }

               --genus;
          }

          std::cout << "finding last waist\n";
          try {
               //NOTE: Finds last waist.
               auto temp = ++(cap.begin());
               splitDangerousEdges(cap.begin().begin(), temp.begin(), temp.begin(), cap.end().begin());
               if(!shortestPathFinder.getShortestPath(temp.begin(), cap.end().begin(), cap.begin().begin(), temp.begin(), IndicesArrays::ArrayAppender(cap))) std::cerr << "ERROR: Failed to connect waists to hole.\n";
               auto path = --cap.end();
               //NOTE: If path length is less than three, then path is a dangerous edge, which should theoretically not happen because of the above call to split dangerous edges.
               auto waist = getShortestLoop(path.begin(), path.end(), cap.begin().begin(), cap.end().begin());
               insert(waist.begin(), waist.end(), *handleIndex);
               --handleIndex;
          } catch(...) { std::cerr << "ERROR: Failed to find last waist.\n"; }

          hpuint nHandlesWaists = *mesh.getGenus();
          if(nHandlesWaists > 3) {
               cap.erase(--cap.end());
               cap.erase(cap.begin());

               //NOTE: Finds base patch's waists.
               while(nHandlesWaists > 3) {
                    cap.erase(cap.end() - 2);
                    auto waist1 = --cap.end();
                    hpuint p2 = m_indices.size() / 3 - 1;
                    auto waist2 = m_boundaries[m_indices[3 * p2]];
                    if(!shortestPathFinder.getShortestPath(waist1.begin(), waist1.end(), waist2.first, waist2.second, cap.begin().begin(), waist1.begin(), IndicesArrays::ArrayAppender(cap))) std::cerr << "ERROR: Failed to connect pair of waists together.\n";
                    cap.push_back(waist2.first, waist2.second);
                    auto pair = cap.end() - 3;
                    splitDangerousEdges(cap.begin().begin(), pair.begin(), pair.begin(), cap.end().begin());
                    if(!shortestPathFinder.getShortestPath(cap.begin().begin(), pair.begin(), pair.begin(), cap.end().begin(), IndicesArrays::ArrayAppender(cap))) std::cerr << "ERROR: Failed to connect pair of waists to cap.\n";
                    auto path = --cap.end();
                    //NOTE: If path length is less than three, then path is a dangerous edge, which should theoretically not happen because of the above call to split dangerous edges.
                    auto waist = getShortestLoop(path.begin(), path.end(), cap.begin().begin(), cap.end().begin());
                    insert(waist.begin(), waist.end(), nHandlesWaists - 2, 0u, p2, 0u);
                    cap.erase(cap.end() - 4, cap.end());
                    --nHandlesWaists;
               }
               insert(0u, 0u, 1u, 0u, unsigned(m_indices.size() / 3 - 1), 0u);
          } else if(nHandlesWaists == 3) insert(0u, 0u, 1u, 0u, 2u, 0u);
     }

     PantsDecomposition<Mesh> decompose() const { return { m_mesh, std::move(m_boundaries), std::move(m_indices), std::move(m_neighbors), std::move(m_reverse) }; }

     Indices getLegCut(hpuint w, hpuint handle) {
          auto waist = m_boundaries[w];

          //get handle loop
          auto h = (*m_mesh.getHandleTunnelLoops()).begin() + (handle << 1);
          std::vector<hpuint> handleLoop(h.begin(), h.end());

          std::vector<hpuint> temp;

          {//get green curve in [Li2008, Fig. 7b]
               using Weigher = HoleyWallWeigher<Mesh, std::vector<hpuint>::const_iterator>;
               using ShortestPathFinder = ShortestPathFinder<Weigher, Mesh>;

               Weigher weigher(m_mesh, handleLoop.cbegin(), handleLoop.cend(), true);
               for(auto i = handleLoop.cbegin(), end = handleLoop.cend(); i != end; ++i) weigher.template puncture<true>(i);
               ShortestPathFinder shortestPathFinder(m_mesh, weigher);
               if(!shortestPathFinder.getShortestPath(handleLoop.cbegin(), handleLoop.cend(), handleLoop.cbegin(), handleLoop.cend(), waist.first, waist.second, temp)) std::cerr << "ERROR: Failed to find handle's temporary tunnel loop.\n";
          }

          IteratorJoiner<IndicesArrays::iterator::iterator> wall;
          wall.push_back(waist.first, waist.second);

          if(temp.front() != temp.back()) {//if green curve does not form a loop, connect the ends with a shortest path
               using Weigher = EdgeLengthWeigher<Mesh>;
               using ShortestPathFinder = ShortestPathFinder<Weigher, Mesh>;

               wall.push_back(temp.begin(), temp.end());
               Weigher weigher(m_mesh);
               ShortestPathFinder shortestPathFinder(m_mesh, weigher);
               hpuint end1 = temp.back();
               temp.pop_back();
               if(!shortestPathFinder.getShortestPath(temp.front(), end1, wall.begin(), wall.end(), temp)) std::cerr << "ERROR: Failed to connect temporary tunnel loop ends.\n";
               wall.pop_back();
          }

          splitDangerousEdges(temp.begin(), temp.end(), wall.begin(), wall.end());
          wall.push_back(temp.begin(), temp.end());

          //get red curve in [Li2008, Fig. 7b]
          return getShortestLoop(temp.begin(), temp.end(), wall.begin(), wall.end());
     }

     template<class SourcesIterator, class WallIterator>
     std::vector<hpuint> getShortestLoop(SourcesIterator sourcesBegin, SourcesIterator sourcesEnd, WallIterator wallsBegin, WallIterator wallsEnd) {
          using Weigher = HoleyWallWeigher<Mesh, SourcesIterator>;
          using ShortestPathFinder = ShortestPathFinder<Weigher, Mesh>;

          Weigher weigher(m_mesh, sourcesBegin, sourcesEnd);
          weigher.removeVertices(wallsBegin, wallsEnd);
          ShortestPathFinder shortestPathFinder(m_mesh, weigher);
          return shortestPathFinder.getShortestLoop(sourcesBegin, sourcesEnd);
     }

     //NOTE: Inserts handles' pants.  The given boundary (begin, end) is the waist and the cut forming the legs must be computed.
     template<class Iterator>
     void insert(Iterator begin, Iterator end, hpuint handle) {
          auto nBoundaries = m_boundaries.size();
          m_boundaries.push_back(begin, end);
          auto leg = getLegCut(nBoundaries, handle);
          m_boundaries.push_back(leg.begin(), leg.end());

          hpuint p = m_indices.size() / 3;

          hpuint indices[] = { nBoundaries, nBoundaries + 1, nBoundaries + 1 };
          m_indices.insert(m_indices.end(), indices, indices + 3);
          hpuint neighbors[] = { UNULL, p, p };
          m_neighbors.insert(m_neighbors.end(), neighbors, neighbors + 3);

          m_reverse[3 * p + 2] = true;
     }

     void insert(hpuint p0, hpuint b0, hpuint handle) {
          hpuint nBoundaries = m_boundaries.size();
          hpuint p = m_indices.size() / 3;
          auto w = m_indices[3 * p0 + b0];

          auto leg = getLegCut(w, handle);
          m_boundaries.push_back(leg.begin(), leg.end());

          hpuint indices[] = { w, nBoundaries, nBoundaries };
          m_indices.insert(m_indices.end(), indices, indices + 3);
          hpuint neighbors[] = { p0, p, p };
          m_neighbors.insert(m_neighbors.end(), neighbors, neighbors + 3);

          m_neighbors[3 * p0 + b0] = p;

          m_reverse[3 * p] = true;
          m_reverse[3 * p + 2] = true;
     }

     //NOTE: Inserts pants whose boundaries are the b0th boundary of the p0th pants, the b1th boundary of the p1th pants, and the b2th boundary of the p2th pants.
     void insert(hpuint p0, hpuint b0, hpuint p1, hpuint b1, hpuint p2, hpuint b2) {
          hpuint p = m_indices.size() / 3;
          hpuint indices[] = { m_indices[3 * p0 + b0], m_indices[3 * p1 + b1], m_indices[3 * p2 + b2] };
          m_indices.insert(m_indices.end(), indices, indices + 3);
          hpuint neighbors[] = { p0, p1, p2 };
          m_neighbors.insert(m_neighbors.end(), neighbors, neighbors + 3);
          m_neighbors[3 * p0 + b0] = p;
          m_neighbors[3 * p1 + b1] = p;
          m_neighbors[3 * p2 + b2] = p;
          m_reverse[3 * p] = true;
          m_reverse[3 * p + 1] = true;
          m_reverse[3 * p + 2] = true;
     }

     //NOTE: Inserts pants in base patch whose boundaries are the given boundary (begin, end), the b1th boundary of the p1th pants, and the b2th boundary of the p2th pants.
     template<class Iterator>
     void insert(Iterator begin, Iterator end, hpuint p1, hpuint b1, hpuint p2, hpuint b2) {
          auto nBoundaries = m_boundaries.size();
          hpuint p = m_indices.size() / 3;

          m_boundaries.push_back(begin, end);

          hpuint indices[] = { nBoundaries, m_indices[3 * p1 + b1], m_indices[3 * p2 + b2] };
          m_indices.insert(m_indices.end(), indices, indices + 3);
          hpuint neighbors[] = { UNULL, p1, p2 };
          m_neighbors.insert(m_neighbors.end(), neighbors, neighbors + 3);
          m_neighbors[3 * p1 + b1] = p;
          m_neighbors[3 * p2 + b2] = p;
          m_reverse[3 * p + 1] = true;
          m_reverse[3 * p + 2] = true;
     }

     template<class SourcesIterator, class TargetsIterator>
     void splitDangerousEdges(SourcesIterator sourcesBegin, SourcesIterator sourcesEnd, TargetsIterator targetsBegin, TargetsIterator targetsEnd) {
          std::vector<hpuint> targets(targetsBegin, targetsEnd);
          std::sort(targets.begin(), targets.end());
          hpuint previous = -1;
          while(sourcesBegin != sourcesEnd) {
               hpuint current = *sourcesBegin;
               if(std::binary_search(targets.begin(), targets.end(), current)) {
                    previous = current;
                    ++sourcesBegin;
                    continue;
               }
               SourcesIterator temp = sourcesBegin+1;
               hpuint next = (temp == sourcesEnd) ? -1 : *temp;
               bool none = true;
               visit_spokes(m_mesh, current, [&](const Edge& edge) {
                    auto target = edge.vertex;
                    if(target != previous && target != next && std::binary_search(targets.begin(), targets.end(), target)) {
                         m_mesh.splitEdge(edge.opposite);
                         none = false;
                    }
               });
               if(none) {
                    previous = current;
                    ++sourcesBegin;
               }
          }
     }

};//PantsDecomposer

}//namespace happah

