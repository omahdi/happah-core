// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

#include "happah/Happah.hpp"
#include "happah/geometries/Mesh.hpp"
#include "happah/weighers/TraversableEdgeLengthWeigher.hpp"

namespace happah {

template<class Mesh, class Iterator>
class HoleyWallWeigher : public TraversableEdgeLengthWeigher<Mesh> {
public:
     using Weight = hpreal;
     static const Weight MAX_WEIGHT;

     //NOTE: Wall cannot cross itself or contain common subsegments.
     HoleyWallWeigher(const Mesh& mesh, Iterator begin, Iterator end, bool loop = false)
          : TraversableEdgeLengthWeigher<Mesh>(mesh, begin, end), m_begin(begin), m_closed(*begin == *(end - 1)), m_end(end), m_loop(loop || m_closed), m_sources(begin, end) {
          if(m_closed) --end;
          while(begin != end) {
               m_cache[*begin] = begin;
               ++begin;
          }
          std::sort(m_sources.begin(), m_sources.end()); 
     }

     void plug(Iterator i) { TraversableEdgeLengthWeigher<Mesh>::removeVertex(*i); }

     //NOTE: Make sure that when calling this method and the wall is not a loop, i is not the first or last element on the path.
     template<bool direction>
     void puncture(Iterator i) {
          assert(m_loop || (i != m_begin && i != m_end - 1));
          auto previous = getPrevious(i);
          auto next = getNext(i);
          puncture<direction>(previous, i, next);
     }

private:
     Iterator m_begin;//wall
     bool m_closed;
     Iterator m_end;
     std::unordered_map<hpuint, Iterator> m_cache;
     bool m_loop;
     std::vector<hpuint> m_sources;

     Iterator getNext(Iterator i) const {
          if(m_loop) return (i == m_end - 1) ? ((m_closed) ? m_begin + 1 : m_begin) : i + 1;
          else return i + 1;
     }

     Iterator getPrevious(Iterator i) const {
          if(m_loop) return (i == m_begin) ? ((m_closed) ? m_end - 2 : m_end - 1) : i - 1;
          else return i - 1;
     }

     bool isWall(hpuint v) const { return std::binary_search(m_sources.cbegin(), m_sources.cend(), v); }

     /**********************************************************************************
      * A given path may not intersect itself but it may come so close to itself that
      * a neighbor of a given vertex is another vertex on the path.  There are two ways
      * how a path may pass along itself when moving from the origin to the destination;
      * the direction of the arrows pointing from the origin to the destination of the 
      * neighboring portions are either parallel or in opposite directions:
      *
      *           ^    ^                        ^    v
      *           |    |                        |    |
      *           ^    ^ vertexNext             ^    v vertexPrevious
      *           |    |                        |    |
      *   current ^--->^ vertex         current ^--->v vertex
      *           |    |                        |    |
      *                  vertexPrevious                vertexNext
      * 
      * In the left case, regardless of the direction of the flow through the path,
      * there is no conflict along the edge connecting the path to itself.  In the
      * right case, however, regardless of the direction of the flow through the path,
      * the flow cannot pass through the connecting edge because one way flow conflict
      * with each other.  Hence, the connecting edge in the right case must be removed.
      * To determine whether we are in the left or right case, we consider the previous
      * and next neighbors of vertex along the path and check the order of the neighbors
      * of vertex.  Starting at vertexNext and moving counterclockwise, if we see
      * vertexPrevious first, then we are in the right case.  If we see current first,
      * then we aee in the left case.
      **********************************************************************************/
     template<bool direction>
     void puncture(Iterator previous, Iterator current, Iterator next) {
          auto p = *previous;
          auto c = *current;
          auto n = *next;

          auto& edges = this->m_mesh.getEdges();

          hpuint v;
          auto i = this->m_mesh.getOutgoing(c);
          while(edges[i].vertex != n) i = edges[edges[i].previous].opposite;

          i = edges[edges[i].previous].opposite;
          while((v = edges[i].vertex) != p) {
               if(isWall(v) && getFlow(current, toIterator(v))) {
                    this->m_removed[i] = 1;
                    this->m_removed[edges[i].opposite] = 1;
               } else {
                    this->m_removed[i] = (direction) ? 1 : 0;
                    this->m_removed[edges[i].opposite] = (direction) ? 0 : 1;
               }
               i = edges[edges[i].previous].opposite;
          }
          this->m_removed[i] = 1;
          this->m_removed[edges[i].opposite] = 1;
          i = edges[edges[i].previous].opposite;
          while((v = edges[i].vertex) != n) {
               if(isWall(v) && getFlow(current, toIterator(v))) {
                    this->m_removed[i] = 1;
                    this->m_removed[edges[i].opposite] = 1;
               } else {
                    this->m_removed[i] = (direction) ? 0 : 1;
                    this->m_removed[edges[i].opposite] = (direction) ? 1 : 0;
               }
               i = edges[edges[i].previous].opposite;
          }
     }

     //NOTE: true if flow in opposite directions and false if in the same direction
     bool getFlow(Iterator current, Iterator vertex) const {
          auto c = *current;
          auto pv = *getPrevious(vertex);
          auto v = *vertex;
          auto nv = *getNext(vertex);

          auto& edges = this->m_mesh.getEdges();

          auto i = this->m_mesh.getOutgoing(v);
          while(edges[i].vertex != nv) i = edges[edges[i].previous].opposite;
          i = edges[edges[i].previous].opposite;
          while(edges[i].vertex != pv && edges[i].vertex != c) i = edges[edges[i].previous].opposite;
          return edges[i].vertex == pv;
     }

     Iterator toIterator(hpuint v) const { return (*m_cache.find(v)).second; }

};//HoleyWallWeigher
template<class Mesh, class Iterator>
const typename HoleyWallWeigher<Mesh, Iterator>::Weight HoleyWallWeigher<Mesh, Iterator>::MAX_WEIGHT = std::numeric_limits<typename HoleyWallWeigher<Mesh, Iterator>::Weight>::max();

template<class Mesh, class Iterator>
HoleyWallWeigher<Mesh, Iterator> make_holey_wall_weigher(const Mesh& mesh, Iterator begin, Iterator end, bool loop = false) { return { mesh, begin, end, loop }; }

}//namespace happah

