// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>

#include "happah/Happah.h"
#include "happah/geometries/Mesh.h"

namespace happah {

template<class Mesh>
class TraversableEdgeLengthWeigher {
public:
     using Weight = hpreal;
     static const Weight MAX_WEIGHT;

     TraversableEdgeLengthWeigher(const Mesh& mesh)
          : m_mesh(mesh), m_removed(mesh.getNumberOfEdges()) {}

     TraversableEdgeLengthWeigher(const Mesh& mesh, boost::dynamic_bitset<> removed)
          : m_mesh(mesh), m_removed(std::move(removed)) {}

     template<class Iterator>
     TraversableEdgeLengthWeigher(const Mesh& mesh, Iterator begin, Iterator end)
          : m_mesh(mesh), m_removed(mesh.getNumberOfEdges()) { removeVertices(begin, end); }

     bool isTraversable(hpuint e) const { return !m_removed[e]; }

     void removeEdge(hpuint e) { m_removed[e] = true; }

     void removeEdge(hpuint v0, hpuint v1) { if(auto e = make_edge_index(m_mesh, v0, v1)) removeEdge(*e); }

     void removeVertex(hpuint v) { this->template set<true>(v); }

     template<class Iterator>
     void removeVertices(Iterator begin, Iterator end) {
          while(begin != end) {
               removeVertex(*begin);
               ++begin;
          }
     }

     void resize(hpuint n) { m_removed.resize(n, false); }

     void unremoveEdge(hpuint e) { m_removed[e] = false; }

     void unremoveVertex(hpuint v) { this->template set<false>(v); }

     Weight weigh(hpuint v0, hpuint v1) const { 
          if(auto i = make_edge_index(m_mesh, v0, v1)) return weigh(v0, v1, *i);
          else return MAX_WEIGHT;
     }

     Weight weigh(hpuint v0, hpuint v1, hpuint edge) const {
          if(m_removed[edge]) return MAX_WEIGHT;
          auto& vertices = m_mesh.getVertices();
          return glm::length(vertices[v0].position - vertices[v1].position);
     }

     template<class Iterator>
     Weight weigh(Iterator begin, Iterator end) const {
          Weight weight = 0;
          for(auto p = begin, q = p + 1; q != end; p = q, ++q) weight += weigh(*p, *q);
          return weight;
     }

protected:
     const Mesh& m_mesh;
     boost::dynamic_bitset<> m_removed;

     template<bool value>
     void set(hpuint v) {
          visit_spokes(m_mesh, m_mesh.getOutgoing(v), [&](const Edge& edge) {
               removeEdge(edge.opposite);
               removeEdge(m_mesh.getEdge(edge.opposite).opposite);
          });
     }

};//TraversableEdgeLengthWeigher
template<class Mesh>
const typename TraversableEdgeLengthWeigher<Mesh>::Weight TraversableEdgeLengthWeigher<Mesh>::MAX_WEIGHT = std::numeric_limits<typename TraversableEdgeLengthWeigher<Mesh>::Weight>::max();

}//namespace happah

