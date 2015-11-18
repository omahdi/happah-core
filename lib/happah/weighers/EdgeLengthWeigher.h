// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.h"
#include "happah/geometries/Mesh.h"

namespace happah {

template<class Mesh>
class EdgeLengthWeigher {
     static_assert(is_absolute_mesh<Mesh>::value, "Weighing without a base can only be done on absolute meshes.");
     using Point = typename Mesh::SPACE::POINT;

public:
     using Weight = hpreal;
     static const Weight MAX_WEIGHT;

     EdgeLengthWeigher(const Mesh& mesh)
          : m_mesh(mesh) {}

     Point getPosition(hpuint v0, hpuint v1, Weight distance) const {
          auto& vertices = m_mesh.getVertices();
          return vertices[v0].position + distance * glm::normalize(vertices[v1].position - vertices[v0].position);
     }

     template<class Iterator>
     Point getPosition(Iterator begin, Iterator end, Weight distance) const {
          for(auto p = begin, q = p + 1; q != end; p = q, ++q) {
               auto temp = weigh(*p, *q);
               if(temp > distance) return getPosition(*p, *q, distance);
               else distance -= temp;
          }
          throw std::runtime_error("Distance must be less than or equal to the path length.");
     }

     Weight weigh(hpuint v0, hpuint v1) const {
          auto& vertices = m_mesh.getVertices();
          return glm::length(vertices[v0].position - vertices[v1].position); 
     }

     Weight weigh(hpuint v0, hpuint v1, hpuint edge) const { return weigh(v0, v1); }

     template<class Iterator>
     Weight weigh(Iterator begin, Iterator end) const {
          Weight weight = 0;
          for(auto p = begin, q = p + 1; q != end; p = q, ++q) weight += weigh(*p, *q);
          return weight;
     }

private:
     const Mesh& m_mesh;

};//EdgeLengthWeigher
template<class Mesh>
const typename EdgeLengthWeigher<Mesh>::Weight EdgeLengthWeigher<Mesh>::MAX_WEIGHT = std::numeric_limits<typename EdgeLengthWeigher<Mesh>::Weight>::max();

template<class Mesh>
EdgeLengthWeigher<Mesh> make_edge_length_weigher(const Mesh& mesh) { return { mesh }; }

}//namespace happah

