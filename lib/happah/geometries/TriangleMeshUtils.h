// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <tuple>
#include <vector>

#include "happah/Happah.h"
#include "happah/utils/VertexFactory.h"

namespace happah {

class TriangleMeshUtils {
public:
     template<class Vertex>
     static void computeAveragedNormals(std::vector<Vertex>& vertices, const std::vector<hpuint>& indices) {
          static_assert(contains_normal<Vertex>::value, "The computation of normals makes sense only for vertices that contain normals.");
          //TODO
     }

     //TODO: even better deindexedarray iterator or something like that
     //NOTE: The last vertex in each triple is the provoking vertex.
     template<class Vertex>
     static void computeFlatNormals(std::vector<Vertex>& vertices, const std::vector<hpuint>& indices) {
          static_assert(contains_normal<Vertex>::value, "The computation of normals makes sense only for vertices that contain normals.");
          
          for(auto i = indices.cbegin(), end = indices.cend(); i != end; ++i) {
               Vertex& v0 = vertices[*i];
               Vertex& v1 = vertices[*(++i)];
               Vertex& v2 = vertices[*(++i)];
               v2.normal = glm::cross(v1.position - v0.position, v2.position - v0.position);
          }
     }

     template<class Vertex, class Base>
     static void computeFlatNormals(std::vector<Vertex>& vertices, const std::vector<hpuint>& indices, const Base& base) {
          static_assert(is_relative_vertex<Vertex>::value && contains_normal<Vertex>::value, "The computation of normals using a base makes sense only for vertices that contain normals.");
          using PointFactory = typename Base::Utils::PointFactory;

          PointFactory factory(base);
          for(auto i = indices.cbegin(), end = indices.cend(); i != end; ++i) {
               Vertex& v0 = vertices[*i];
               Vertex& v1 = vertices[*(++i)];
               Vertex& v2 = vertices[*(++i)];
               auto p0 = factory.build(v0.abscissa, v0.ordinate);
               v2.normal = glm::cross(factory.build(v1.abscissa, v1.ordinate) - p0, factory.build(v2.abscissa, v2.ordinate) - p0);
          }
     }

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

     /** 
      * Given a set of triangle indices and an edge specified by the indices of its endpoints, this finds the triangle containing the edge, if any.
      *
      * @note An edge is shared by at most two triangles.  This function returns the first adjacent triangle it finds.
      */
     static hpuint findTriangle(const std::vector<hpuint>& indices, hpuint v0, hpuint v1, hpuint& v2);

     static std::tuple<hpuint, hpuint, hpuint> getNeighbors(const std::vector<hpuint>& indices, hpuint triangle);

};//TriangleMeshUtils

}//namespace happah

