// Copyright 2016
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Tobias Ribizel - Karlsruhe Institute of Technology - upsj@upsj.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2016.10 - Pawel Herman     - Refactored for reusability and readability.
// 2016.09 - Tobias Ribizel   - Initial commit.

#pragma once

#include "happah/geometries/TriangleMesh.h"

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <tuple>

namespace happah {

class VertexRule {
public:
     template<class Iterator, class Vertex = typename Iterator::value_type>
     Vertex operator()(const Vertex& center, Iterator begin, Iterator end) {
          auto valence = end - begin;

          Vertex mean;
          while(begin != end) {
               mean.position += (*begin).position;
               ++begin;
          }
          mean.position *= 1.f / valence;

          auto alpha = 3.f / 8.f + 2.f / 8.f * (float)glm::cos(2 * glm::pi<double>() / valence);
          alpha = 5.f / 8.f - alpha * alpha;

          return {alpha * mean.position + (1 - alpha) * center.position};
     }     

};//VertexRule

class EdgeRule {
public:
     template<class Vertex>
     Vertex operator()(const Vertex& vertex0, const Vertex& vertex1, const Vertex& vertex2, const Vertex& vertex3) const { return (3.f * (vertex0.position + vertex1.position) + (vertex2.position + vertex3.position)) / 8.f; }

};//EdgeRule

template<class VertexRule, class EdgeRule>
class LoopSubdivider {
public:
     LoopSubdivider(VertexRule vertexRule, EdgeRule edgeRule)
          : m_edgeRule(std::move(edgeRule)), m_vertexRule(std::move(vertexRule)) {}

     /*
      * The ith triangle in the input mesh is replaced by the (4i)th, (4i+1)th, (4i+2)th, and (4i+3)th triangles in the output mesh.  The order of the output triangles is given by the diagram below.  The order of the corresponding vertices is { { 0, 1, 3 }, { 1, 2, 4 }, { 1, 4, 3 }, { 4, 5, 3 } }.
      *
      *   INPUT            w
      *                   / \
      *                  /   \
      *                 /     \
      *                /       \
      *               /         \
      *              /           \
      *             u-------------v
      *
      *   OUTPUT          10
      *                   / \
      *                  /   \
      *                 11----9
      *                 8-----7
      *               2  \   /  5
      *              / \  \ /  / \
      *             /   \  6  /   \
      *            0-----1   3-----4
      */
     template<class Vertex>
     TriangleMesh<Vertex> subdivide(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh) {
          auto& vertices0 = mesh.getVertices();
          auto& indices0 = mesh.getIndices();

          std::vector<Vertex> vertices1;
          std::vector<hpuint> indices1;

          auto nVertices = mesh.getNumberOfVertices();
          auto nTriangles = mesh.getNumberOfTriangles();
          auto nEdges = mesh.getNumberOfEdges();

          vertices1.reserve(nVertices + nEdges / 2);
          indices1.reserve(nTriangles * 3 * 4);

          std::vector<hpuint> es(nEdges, -1);//TODO: in case there are future stringent memory requirements, this vector can be avoided

          auto v = 0u;
          for(auto& vertex : vertices0) {
               std::vector<hpuint> temp;
               visit_ring(mesh.getEdges(), mesh.getNumberOfTriangles(), mesh.getOutgoing(v++), [&](hpuint w) { temp.push_back(w); });
               auto ring = deindex(vertices0, temp);
               vertices1.emplace_back(m_vertexRule(vertex, begin(ring), end(ring)));
          }

          for(auto e = 0u; e < nEdges; ++e) {
               if(es[e] != UNULL) continue;
               auto edge0 = mesh.getEdge(e);
               auto edge1 = mesh.getEdge(edge0.opposite);
               hpuint v0, v1, v2, v3;
               v0 = edge0.vertex;
               v1 = edge1.vertex;
               v2 = mesh.getEdge(edge0.next).vertex;
               v3 = mesh.getEdge(edge1.next).vertex;
               es[e] = vertices1.size();
               es[edge0.opposite] = vertices1.size();
               vertices1.emplace_back(m_edgeRule(vertices0[v0], vertices0[v1], vertices0[v2], vertices0[v3]));
          }

          auto t = 0u;
          for(auto i = indices0.begin(), end = indices0.end(); i != end; ++i, ++t) {
               auto u = *i;        // vertex 0
               auto v = *(++i);    // vertex 2
               auto w = *(++i);    // vertex 5
               auto uv = es[t];    // vertex 1
               auto vw = es[++t];  // vertex 4
               auto wu = es[++t];  // vertex 3

               indices1.insert(indices1.end(), {
                    // see TriangleRefinementScheme.cpp
                    u,  uv, wu, // 0 1 3
                    uv,  v, vw, // 1 2 4
                    uv, vw, wu, // 1 4 3
                    vw,  w, wu, // 4 5 3
               });
          }

          return {vertices1, indices1};
     }

private:
     EdgeRule m_edgeRule;
     VertexRule m_vertexRule;

};//LoopSubdivider

LoopSubdivider<VertexRule, EdgeRule> make_loop_subdivider() { return {VertexRule(), EdgeRule()}; }

template<class VertexRule, class EdgeRule>
LoopSubdivider<VertexRule, EdgeRule> make_loop_subdivider(VertexRule vertexRule, EdgeRule edgeRule) { return {std::forward(vertexRule), std::forward(edgeRule)}; }

}//namespace happah

