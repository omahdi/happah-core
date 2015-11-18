// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/utils/ShortestPathFinder.h"
#include "happah/weighers/EdgeLengthWeigher.h"

class MeshUtils {
public:
     template<class Mesh, class Iterator>
     static void sort(const Mesh& mesh, Iterator begin, Iterator end) {
          using r = std::reverse_iterator<std::vector<hpuint>::const_iterator>;
          using Weigher = EdgeLengthWeigher<Mesh>;
          using Weight = typename Weigher::Weight;
          
          //TODO: memory use can be improved using tree structure
          //TODO: this is actually a depth first search; just need to remember where to pick up when moving up
          std::vector<std::vector<hpuint> > paths;
          paths.push_back(std::vector<hpuint>(1, *begin));
          std::vector<std::set<hpuint> > rests;
          rests.push_back(std::set<hpuint>(begin+1, end));
          Weigher weigher(mesh);
          std::vector<hpuint> shortestPath;
          Weight shortestPathLength = Weigher::MAX_WEIGHT;
          while(true) {
               while(!rests.empty() && rests.back().size() > 0) {
                    auto neighbors = mesh.getRing(paths.back().back());
                    std::vector<hpuint> nexts;
                    for(hpuint neighbor : neighbors) {
                         auto next = rests.back().find(neighbor);
                         if(next != rests.back().end()) nexts.push_back(*next);
                    }
                    if(nexts.size() == 0) {
                         paths.pop_back();
                         rests.pop_back();
                    } else if(nexts.size() == 1) {
                         paths.back().push_back(nexts.front());
                         rests.back().erase(rests.back().find(nexts.front()));
                    } else {
                         std::vector<hpuint> path(paths.back());
                         std::set<hpuint> rest(rests.back());
                         paths.pop_back();
                         rests.pop_back();
                         for(hpuint next : nexts) {
                              std::vector<hpuint> cpath(path);
                              std::set<hpuint> crest(rest);
                              cpath.push_back(next);
                              crest.erase(crest.find(next));
                              paths.push_back(cpath);
                              rests.push_back(crest);
                         }
                    }
               }
               if(rests.empty()) break;
               Weight pathLength = ShortestPathFinderUtils::getPathLength(weigher, r(paths.back().cend()), r(paths.back().cbegin()));
               if(pathLength < shortestPathLength) {
                    shortestPath.clear();
                    shortestPath.insert(shortestPath.end(), paths.back().begin(), paths.back().end());
                    shortestPathLength = pathLength;
               }
               paths.pop_back();
               rests.pop_back();
          }
          std::copy(shortestPath.begin(), shortestPath.end(), begin);
     }

};

