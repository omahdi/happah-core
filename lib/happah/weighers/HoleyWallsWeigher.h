// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

#include "happah/Happah.h"
#include "happah/weighers/HoleyWallWeigher.h"
#include "happah/weighers/EdgeLengthWeigher.h"

namespace happah {

template<class Mesh, class Iterator>
class HoleyWallsWeigher {
public:
     using Weight = typename HoleyWallWeigher<Mesh, Iterator>::Weight;
     static const Weight MAX_WEIGHT;

     HoleyWallsWeigher(const Mesh& mesh, const std::vector<Iterator>& walls, const boost::dynamic_bitset<>& loops) 
          : m_mesh(mesh) {
          {
               auto w = 0u;
               for(auto w0 = walls.cbegin(), end = walls.cend(); w0 != end; w0 += 2, ++w) m_weighers.emplace_back(mesh, *w0, *(w0 + 1), loops[w]);
          }
          {
               auto w = m_weighers.begin();
               for(auto w0 = walls.cbegin(), end = walls.cend(); w0 != end; w0 += 2, ++w) {
                    for(auto i0 = *w0, i1 = *(w0 + 1); i0 != i1; ++i0) m_cache.emplace(*i0, std::ref(*w));
               }
          }
     }

     void plug(Iterator i) { m_cache.find(*i)->second.get().plug(i); }

     //NOTE: Make sure that when calling this method and the wall is not a loop, i is not the first or last element on the path.
     template<bool direction>
     void puncture(Iterator i) { m_cache.find(*i)->second.get().template puncture<direction>(i); }

     Weight weigh(hpuint v0, hpuint v1) const { 
          if(auto i = make_edge_index(m_mesh, v0, v1)) return weigh(v0, v1, *i);
          else return MAX_WEIGHT;
     }

     Weight weigh(hpuint v0, hpuint v1, hpuint edge) const {
          for(auto& weigher : m_weighers) if(!weigher.isTraversable(edge)) return MAX_WEIGHT;
          auto& vertices = m_mesh.getVertices();
          return glm::length(vertices[v0].position - vertices[v1].position);
     }
private:
     std::unordered_map<hpuint, std::reference_wrapper<HoleyWallWeigher<Mesh, Iterator> > > m_cache;
     const Mesh& m_mesh;
     std::vector<HoleyWallWeigher<Mesh, Iterator> > m_weighers;

};//HoleyWallsWeigher
template<class Mesh, class Iterator>
const typename HoleyWallsWeigher<Mesh, Iterator>::Weight HoleyWallsWeigher<Mesh, Iterator>::MAX_WEIGHT = HoleyWallWeigher<Mesh, Iterator>::MAX_WEIGHT;

}//namespace happah

