// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <tuple>
#include <vector>

#include "happah/Happah.hpp"
#include "happah/utils/VertexFactory.hpp"

namespace happah {

class TriangleMeshUtils {
public:
     template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
     static void extrude(std::vector<Vertex>& vertices, std::vector<hpuint>& indices, hpreal z, VertexFactory&& factory = VertexFactory()) {
          hpuint nIndices = indices.size();
          indices.reserve(nIndices << 1);
          extrude(vertices, indices.begin(), nIndices, std::back_inserter(indices), z, factory);
     }

     template<class Vertex, class IndicesInputIterator, class IndicesOutputIterator, class VertexFactory = VertexFactory<Vertex> >
     static void extrude(std::vector<Vertex>& vertices, IndicesInputIterator k, hpuint nIndices, IndicesOutputIterator l, hpreal z, VertexFactory factory = VertexFactory()) {
          hpuint nVertices = vertices.size();
          vertices.reserve(nVertices << 1);
          auto i = vertices.begin();
          auto j = std::back_inserter(vertices);
          hpuint n = nVertices;
          while(n > 0) {
               *j = factory(Point3D((*i).position.x, (*i).position.y, (*i).position.z + z));
               ++i, ++j;
               --n;
          }
          k += nIndices;
          while(nIndices > 0) {
               --k;
               *l = (*k) + nVertices;
               ++l;
               --nIndices;
          }
     }

};//TriangleMeshUtils

}//namespace happah

