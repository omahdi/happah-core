// Copyright 2016
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Tobias Ribizel - Karlsruhe Institute of Technology - upsj@upsj.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2016.10 - Pawel Herman     - Refactored for reusability and readability.
// 2016.09 - Tobias Ribizel   - Initial commit.

#pragma once

#include "happah/geometries/TriangleMesh.h"

namespace happah {

/*
 * The ith triangle in the input mesh is replaced by the (4i)th, (4i+1)th, (4i+2)th, and (4i+3)th triangles in the output mesh.  The order of the output triangles is given by the diagram below.  The order of the corresponding vertices is { { 0, 1, 3 }, { 1, 2, 4 }, { 1, 4, 3 }, { 4, 5, 3 } } and is the same ordering as in the BINARY_UNIFORM triangle refinement scheme.
 *
 *   INPUT            
 *                   /\
 *                  /  \
 *                 /    \
 *                /      \
 *               /        \
 *              /          \
 *              ------------
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
template<class Vertex, class VertexRule, class EdgeRule>
TriangleMesh<Vertex> loopivide(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, VertexRule&& vertexRule, EdgeRule&& edgeRule) {
     auto vertices = std::vector<Vertex>();
     auto indices = Indices();
     auto es = Indices(mesh.getNumberOfEdges(), std::numeric_limits<hpindex>::max());

     vertices.reserve(mesh.getNumberOfVertices() + mesh.getNumberOfEdges() >> 1);
     indices.reserve((mesh.getNumberOfTriangles() << 2) * 3);

     auto v = 0u;
     for(auto& vertex : mesh.getVertices()) {
          auto ring = make_ring(mesh, v++);
          vertices.emplace_back(vertexRule(vertex, std::begin(ring), std::end(ring)));
     }

     visit_diamonds(mesh, [&](auto e, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
          auto& edge = mesh.getEdge(e);
          es[e] = vertices.size();
          es[edge.opposite] = vertices.size();
          vertices.emplace_back(edgeRule(vertex0, vertex2, vertex1, vertex3));
     });

     auto e = std::begin(es) - 1;
     visit_triplets(mesh.getIndices(), [&](auto v0, auto v2, auto v5) {
          auto v1 = *(++e);
          auto v4 = *(++e);
          auto v3 = *(++e);

          indices.insert(std::end(indices), {
               v0, v1, v3,
               v1, v2, v4,
               v1, v4, v3,
               v4, v5, v3
          });
     });

     return make_triangle_mesh(vertices, indices);
}

template<class Vertex>
TriangleMesh<Vertex> loopivide(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh) {
     return loopivide(mesh, [](auto& center, auto begin, auto end) {
          auto valence = std::distance(begin, end);

          auto mean = Vertex();
          while(begin != end) {
               mean.position += (*begin).position;
               ++begin;
          }
          mean.position *= 1.f / valence;

          auto alpha = 3.f / 8.f + 2.f / 8.f * (hpreal)glm::cos(2 * glm::pi<hpreal>() / valence);
          alpha = 5.f / 8.f - alpha * alpha;

          return alpha * mean.position + (1 - alpha) * center.position;
     }, [](auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) { return (3.f * (vertex0.position + vertex1.position) + (vertex2.position + vertex3.position)) / 8.f; });
}

}//namespace happah

