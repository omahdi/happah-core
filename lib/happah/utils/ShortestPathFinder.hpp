// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <vector>

#include "happah/geometries/Mesh.hpp"

namespace happah {

struct ShortestPathFinderUtils {
     //NOTE: Input is a loop.  We assume source is on the loop.
     template<class Iterator, class Weigher>
     static Iterator getFarthestPoint(const Weigher& weigher, hpuint source, Iterator begin, Iterator end) {
          using Weight = typename Weigher::Weight;

          if(*begin == *(end - 1)) --end;

          std::vector<Weight> distances(std::distance(begin, end), Weigher::MAX_WEIGHT);
          auto s0 = begin;
          auto i0 = distances.begin();
          while(*s0 != source) { ++s0; ++i0; }
          *i0 = 0;
          Iterator target;
          Weight maxDistance = 0;

          if(s0 == begin) {
               {
                    Weight total = 0;
                    auto previous = s0;
                    auto s = s0 + 1;
                    auto i = i0 + 1;
                    while(s != end) {
                         total += weigher.weigh(*previous, *s);
                         *i = total;
                         ++previous; ++s; ++i;
                    }
               }
               {
                    auto s = end - 1;
                    auto i = distances.end() - 1;
                    auto total = weigher.weigh(*begin, *s);
                    auto previous = end;
                    while(s != begin) {
                         Weight distance = std::min(*i, total);
                         if(distance > maxDistance) {
                              target = s;
                              maxDistance = distance;
                         }
                         --previous; --s; --i;
                         total += weigher.weigh(*previous, *s);
                    }
               }
          } else {
               {
                    Weight total = 0;
                    {
                         Iterator previous = s0;
                         Iterator s = s0 + 1;
                         auto i = i0 + 1;
                         while(s != end) {
                              total += weigher.weigh(*previous, *s);
                              *i = total;
                              ++previous; ++s; ++i;
                         }
                    }
                    {
                         auto i = distances.begin();
                         Iterator previous = begin - 1;
                         Iterator s = begin;
                         total += weigher.weigh(*(end - 1), *begin);
                         while(s != s0) {
                              *i = total;
                              ++previous; ++s; ++i;
                              total += weigher.weigh(*previous, *s);
                         }
                    }
               }
               {
                    Weight total = 0;
                    {
                         Iterator previous = s0;
                         Iterator s = s0 - 1;
                         auto i = i0 - 1;
                         while(previous != begin) {
                              total += weigher.weigh(*previous, *s);
                              Weight distance = std::min(*i, total);
                              if(distance > maxDistance) {
                                   target = s;
                                   maxDistance = distance;
                              }
                              --previous; --s; --i;
                         }
                    }
                    {
                         Iterator s = end - 1;
                         auto i = distances.end() - 1;
                         total += weigher.weigh(*begin, *s);
                         Iterator previous = end;
                         while(s != s0) {
                              Weight distance = std::min(*i, total);
                              if(distance > maxDistance) {
                                   target = s;
                                   maxDistance = distance;
                              }
                              --previous, --s; --i;
                              total += weigher.weigh(*previous, *s);
                         }
                    }
               }
          }

          return target;
     }

     //NOTE: There is no check if the given vertex is on the path.
     template<class Iterator>
     static std::pair<hpuint, hpuint> getNeighborhood(Iterator begin, Iterator end, hpuint n) {
          if(*begin == n) {
               Iterator last = end-1;
               return (*begin == *last) ? std::make_pair(*(last-1),*(begin+1)) : std::make_pair(happah::UNULL, *(begin+1));
          }
          Iterator previous = begin;
          Iterator current = previous+1;
          Iterator next = current+1;
          while(*current != n && next != end) {
               previous = current;
               current = next;
               ++next;
          }
          return (next == end) ? std::make_pair(*previous, happah::UNULL) : std::make_pair(*previous, *next);
     }

     template<class Path>
     static void getPath(const std::vector<hpuint>& predecessors, hpuint target, Path& path) {
          hpuint current = target;
          hpuint predecessor = predecessors[current];
          do {
               path.push_back(current);
               current = predecessor;
               assert(current != UNULL);
               predecessor = predecessors[current];
          } while(predecessor != current && current != target);
          path.push_back(current);
     }

     template<class Weigher, class Iterator>
     static typename Weigher::Weight getPathLength(const Weigher& weigher, Iterator pathBegin, Iterator pathEnd) {//TODO: move to Weigher structures; Weigher.weigh(begin, end);
          typename Weigher::Weight length = 0;
          for(auto p = pathBegin, q = p + 1; q != pathEnd; p = q, ++q) length += weigher.weigh(*p, *q);
          return length;
     }

};

//NOTE: This implementation works on non-directed as well as directed meshes.
template<class Weigher, class Mesh>
class ShortestPathFinder {
     using Weight = typename Weigher::Weight;
     //TODO: adapt implementation to serial version if compiling for non intel processor
     //TODO: see also Priority Queues and Dijkstra’s Algorithm by Chen, Chowdhury, Roche, Tong, Ramachandran
     //TODO: sequence heap; see Fast priority queues for cached memory by Peter Sanders

public:
     ShortestPathFinder(const Mesh& mesh, Weigher& weigher)//TODO: fix this back to const if possible or make copy; getShortestLoop changes weigher
          : m_mesh(mesh), m_weigher(weigher) {}

     template<class SourcesIterator, bool direction = true>
     std::vector<hpuint> getShortestLoop(SourcesIterator sourcesBegin, SourcesIterator sourcesEnd) {//TODO: refactor
          using r = std::reverse_iterator<std::vector<hpuint>::const_iterator>;

          //TODO: HoleyWallWeigher; ++ move hole over one, -- move hole over one back
          Weight shortestLoopLength = Weigher::MAX_WEIGHT;
          std::vector<hpuint> loop;
          bool closed = *sourcesBegin == *(sourcesEnd - 1);
          for(auto i = sourcesBegin + 1, end = (closed) ? sourcesEnd : sourcesEnd - 1; i != end; ++i) {
               m_weigher.template puncture<direction>(i);
               std::vector<hpuint> temp;
               if(getShortestPath(*i, *i, temp)) {//TODO: stop if path length is longer than current shortest
                    Weight loopLength = ShortestPathFinderUtils::getPathLength(m_weigher, r(temp.cend()), r(temp.cbegin()));
                    if(loopLength < shortestLoopLength) {
                         loop = std::move(temp);
                         shortestLoopLength = loopLength;
                    }
               }
               m_weigher.plug(i);
          }
          if(loop.size() > 0) return loop;
          throw std::runtime_error("Failed to find shortest loop.");
     }

     template<class Paths>
     void getShortestPaths(hpuint source, Paths& paths) {
          hpuint* sourcesBegin = &source;
          hpuint* sourcesEnd = sourcesBegin + 1;
          doGetShortestPath<0>(sourcesBegin, sourcesEnd, (hpuint*)NULL, (hpuint*)NULL, (hpuint*)NULL, (hpuint*)NULL, paths);
     }

     template<class Iterator, class Paths>
     void getShortestPaths(Iterator sourcesBegin, Iterator sourcesEnd, Paths& paths) { doGetShortestPath<0>(sourcesBegin, sourcesEnd, (hpuint*)NULL, (hpuint*)NULL, (hpuint*)NULL, (hpuint*)NULL, paths); }

     template<class Path>
     bool getShortestPath(hpuint source, hpuint target, Path&& path) {
          hpuint* sourcesBegin = &source;
          hpuint* sourcesEnd = sourcesBegin + 1;
          return getShortestPath(sourcesBegin, sourcesEnd, target, path);
     }

     template<class Iterator, class Path>
     bool getShortestPath(hpuint source, Iterator targetsBegin, Iterator targetsEnd, Path&& path) {
          hpuint* sourcesBegin = &source;
          hpuint* sourcesEnd = sourcesBegin + 1;
          return getShortestPath(sourcesBegin, sourcesEnd, targetsBegin, targetsEnd, path);
     }

     template<class SourcesIterator, class TargetsIterator, class Path>
     bool getShortestPath(SourcesIterator sourcesBegin, SourcesIterator sourcesEnd, TargetsIterator targetsBegin, TargetsIterator targetsEnd, Path&& path) { return doGetShortestPath<2>(sourcesBegin, sourcesEnd, targetsBegin, targetsEnd, (hpuint*)NULL, (hpuint*)NULL, path); }

     template<class SourcesIterator, class TargetsIterator, class WallsIterator, class Path>
     bool getShortestPath(SourcesIterator sourcesBegin, SourcesIterator sourcesEnd, TargetsIterator targetsBegin, TargetsIterator targetsEnd, WallsIterator wallsBegin, WallsIterator wallsEnd, Path&& path) { return doGetShortestPath<2>(sourcesBegin, sourcesEnd, targetsBegin, targetsEnd, wallsBegin, wallsEnd, path); }

     template<class TargetsIterator, class WallsIterator, class Path>
     bool getShortestPath(hpuint source, TargetsIterator targetsBegin, TargetsIterator targetsEnd, WallsIterator wallsBegin, WallsIterator wallsEnd, Path&& path) { 
          hpuint* sourcesBegin = &source;
          hpuint* sourcesEnd = sourcesBegin + 1;
          return doGetShortestPath<2>(sourcesBegin, sourcesEnd, targetsBegin, targetsEnd, wallsBegin, wallsEnd, path); 
     }

     template<class WallsIterator, class Path>
     bool getShortestPath(hpuint source, hpuint target, WallsIterator wallsBegin, WallsIterator wallsEnd, Path&& path) { 
          hpuint* sourcesBegin = &source;
          hpuint* sourcesEnd = sourcesBegin + 1;
          hpuint* targetsBegin = &target;
          hpuint* targetsEnd = targetsBegin + 1;
          return doGetShortestPath<1>(sourcesBegin, sourcesEnd, targetsBegin, targetsEnd, wallsBegin, wallsEnd, path); 
     }

     template<class Iterator, class Path>
     bool getShortestPath(Iterator sourcesBegin, Iterator sourcesEnd, hpuint target, Path& path) { 
          hpuint* targetsBegin = &target;
          hpuint* targetsEnd = targetsBegin + 1;
          return doGetShortestPath<1>(sourcesBegin, sourcesEnd, targetsBegin, targetsEnd, (hpuint*)NULL, (hpuint*)NULL, path);
     }

private:
     const Mesh& m_mesh;
     Weigher& m_weigher;

     template<hpuint t_nTargets, class SourcesIterator, class TargetsIterator, class WallsIterator, class Path>
     bool doGetShortestPath(SourcesIterator sourcesBegin, SourcesIterator sourcesEnd, TargetsIterator targetsBegin, TargetsIterator targetsEnd, WallsIterator wallsBegin, WallsIterator wallsEnd, Path& path) {
          auto nTargets = std::distance(targetsBegin, targetsEnd);
          hpuint targets[nTargets];
          std::copy(targetsBegin, targetsEnd, targets);
          auto temp = targets + nTargets;
          std::sort(targets, temp);
          auto nVertices = m_mesh.getVertices().size();
          std::vector<hpuint> predecessors;
          if(t_nTargets == 0) path.resize(nVertices, UNULL);
          else predecessors.resize(nVertices, UNULL);
          std::vector<Weight> distances(nVertices, Weigher::MAX_WEIGHT);
          boost::dynamic_bitset<> todo(nVertices);
          todo.set();
          for(auto i = wallsBegin; i != wallsEnd; ++i) todo[*i] = false;
          do {
               auto source = *sourcesBegin;
               distances[source] = 0;
               todo[source] = true;//NOTE: Source may be on the wall.
               if(t_nTargets == 0) path[source] = source;
               else predecessors[source] = source;
          } while((++sourcesBegin) != sourcesEnd);
          for(auto i = 0lu; i < nVertices; ++i) {
               hpuint vertex;
               vertex = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));
               auto distance = distances[vertex];
               if(distance == Weigher::MAX_WEIGHT) break;
               switch(t_nTargets) {
               case 1:
                    if(*targetsBegin == vertex && predecessors[vertex] != vertex) {
                         ShortestPathFinderUtils::getPath(predecessors, vertex, path);
                         return true;
                    }
                    break;
               case 2:
                    if(std::binary_search(targets, temp, vertex) && predecessors[vertex] != vertex) {
                         ShortestPathFinderUtils::getPath(predecessors, vertex, path);
                         return true;
                    }
                    break;
               default:
                    break;    
               }
               visit_spokes(make_spokes_enumerator(m_mesh.getEdges(), m_mesh.getOutgoing(vertex)),
                    [&](auto ei) {
                         const auto neighbor = m_mesh.getEdge(ei).vertex;
                         if(todo[neighbor] || std::binary_search(targets, temp, neighbor)) {
                              auto delta = m_weigher.weigh(vertex, neighbor);
                              if((distance + delta) < distances[neighbor]) {
                                   if(t_nTargets == 0) path[neighbor] = vertex;
                                   else predecessors[neighbor] = vertex;
                                   distances[neighbor] = distance + delta;
                              }
                         }
                    });
               distances[vertex] = Weigher::MAX_WEIGHT;
               todo[vertex] = false;
          }
          return t_nTargets == 0;
     }

};//ShortestPathFinder

template<class Weigher, class Mesh>
ShortestPathFinder<Weigher, Mesh> make_shortest_path_finder(const Mesh& mesh, Weigher& weigher) { return { mesh, weigher }; }

}//namespace happah

