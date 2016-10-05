// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2016.10 - Pawel Herman     - Refactoring for reusability and readability.
// 2016.09 -                  - Initial commit.

#pragma once

#include "happah/geometries/TriangleMesh.h"

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <tuple>

namespace happah {

template<class Iterator>
class VertexRule {
     using Vertex = typename Iterator::value_type;

public:
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

};

template<class Vertex>
class LoopSubdivider {
private:
    using OutputMesh = TriangleMesh<Vertex, Format::SIMPLE>;
    using InputMesh = TriangleMesh<Vertex, Format::DIRECTED_EDGE>;
    using Edge = std::pair<hpuint, hpuint>;

    std::unordered_map<Edge, hpuint, boost::hash<Edge>> m_edge_index;//TODO: remove

    hpuint edge_index(hpuint v, hpuint w) {
        if (v >= w) std::swap(v, w);
        return m_edge_index[Edge(v, w)];
    }

public:
     //TODO: LoopSubdivider();
     //TODO: template<class VertexRule, class EdgeRule>
     //LoopSubdivider(VertexRule vertexRule, EdgeRule edgeRule);

     OutputMesh subdivide(const InputMesh& mesh) {
          using VertexRule = VertexRule<typename DeindexedArray<typename InputMesh::Vertices, typename InputMesh::Indices>::const_iterator>;

          auto& vertices0 = mesh.getVertices();
          auto& indices0 = mesh.getIndices();

          hpuint n = vertices0.size();
          hpuint f = mesh.getNumberOfTriangles();
          hpuint e = mesh.getNumberOfEdges();

          std::vector<Vertex> vertices1;
          std::vector<hpuint> indices1;

          vertices1.reserve(n + e / 2);
          indices1.reserve(f * 3 * 4);

          m_edge_index.clear();//TODO: remove
          m_edge_index.reserve(e / 2);

          // compute new vertex points
          auto v = 0u;
          VertexRule vertexRule;
          for(auto& vertex : vertices0) {
               auto temp = mesh.getRing(v++);
               auto ring = deindex(vertices0, temp);
               vertices1.emplace_back(vertexRule(vertex, begin(ring), end(ring)));
          }

          // compute new edge points
          EdgeRule edgeRule;
          hpuint edge_vertex = vertices1.size();
          for (hpuint iedge = 0; iedge < e; ++iedge) {
               auto edge = mesh.getEdge(iedge);
               auto edge_rev = mesh.getEdge(edge.opposite);
               hpuint v, w, x, y;
               v = edge.vertex;
               w = edge_rev.vertex;
               // only process each full edge once
               if (v < w) {
                    // find neighboring vertices
                    x = mesh.getEdge(edge.next).vertex;
                    y = mesh.getEdge(mesh.getEdge(edge.opposite).next).vertex;
                    m_edge_index.insert(std::make_pair(Edge(v, w), edge_vertex));
                    vertices1.emplace_back(edgeRule(vertices0[v], vertices0[w], vertices0[x], vertices0[y]));
                    ++edge_vertex;
               }
          }

          // build new index buffer
          for (hpuint iface = 0; iface < f; ++iface) {
               auto u = indices0[3 * iface];     // vertex 0
               auto v = indices0[3 * iface + 1]; // vertex 2
               auto w = indices0[3 * iface + 2]; // vertex 5
               auto uv = edge_index(u, v);        // vertex 1
               auto vw = edge_index(v, w);        // vertex 4
               auto wu = edge_index(w, u);        // vertex 3

               indices1.insert(indices1.end(), {
                    // see TriangleRefinementScheme.cpp
                    u,  uv, wu, // 0 1 3
                    uv,  v, vw, // 1 2 4
                    uv, vw, wu, // 1 4 3
                    vw,  w, wu, // 4 5 3
               });
          }

          return OutputMesh(vertices1, indices1);
     }

};//LoopSubdivider

}//namespace happah

