// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/range/iterator_range.hpp>
#include <cmath>
#include <glm/gtc/constants.hpp>

#include "happah/geometries/TriangleGraph.hpp"
#include "happah/geometries/SurfaceSplineBEZ.hpp"
#include "happah/utils/visitors.hpp"

namespace happah {

//Convert a closed(!) triangle graph into a quartic polynomial spline.
template<class Vertex>
auto make_spline_surface(const TriangleGraph<Vertex>& graph) {
     using Space = typename Vertex::SPACE;
     using Vector = typename Space::VECTOR;

     auto surface = QuarticSurfaceSplineBEZ<Space>(graph.getNumberOfTriangles());

     auto set_boundary_point = [&](auto t, auto i, auto k, auto&& point) {
          auto u = make_neighbor_index(graph, t, i);
          auto j = make_neighbor_offset(graph, u, t);
          surface.setBoundaryPoint(t, i, k, u, j, point);
     };

     visit_diamonds(graph, [&](auto e, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
          auto t = make_triangle_index(e);
          auto i = make_edge_offset(e);
          set_boundary_point(t, i, 1, (1.0f / 6.0f) * (2.0f * vertex0.position + vertex1.position + 2.0f * vertex2.position + vertex3.position));
     });

     for(auto v : boost::irange(0u, graph.getNumberOfVertices())) {
          auto ring = make_ring(make_ring_enumerator(graph, v));
          auto begin = std::begin(ring);
          auto end = std::end(ring);
          auto t = make_triangle_index(graph.getOutgoing(v));
          auto i = make_edge_offset(graph.getOutgoing(v));
          auto& center = graph.getVertex(v);
          auto fan = Indices();
          visit_spokes(make_spokes_enumerator(graph.getEdges(), graph.getOutgoing(v)), [&](auto e) {
               auto u = make_triangle_index(e);
               auto j = make_edge_offset(e);
               fan.push_back(u);
               fan.push_back(j);
          });
          auto valence = std::distance(begin, end);

          auto corner = center.position;
          if(valence == 6) {
               corner *= 6.0f;
               for(auto& vertex : boost::make_iterator_range(begin, end)) corner += vertex.position;
               corner /= 12.0f;

               auto make_boundary_point = [&](auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3, auto& vertex4) -> auto { return (1.0f / 24.0f) * (12.0f * center.position + vertex0.position + 3.0f * vertex1.position + 4.0f * vertex2.position + 3.0f * vertex3.position + vertex4.position); };

               set_boundary_point(fan[0], fan[1], 0, make_boundary_point(begin[4], begin[5], begin[0], begin[1], begin[2]));
               set_boundary_point(fan[2], fan[3], 0, make_boundary_point(begin[5], begin[0], begin[1], begin[2], begin[3]));
               set_boundary_point(fan[4], fan[5], 0, make_boundary_point(begin[0], begin[1], begin[2], begin[3], begin[4]));
               set_boundary_point(fan[6], fan[7], 0, make_boundary_point(begin[1], begin[2], begin[3], begin[4], begin[5]));
               set_boundary_point(fan[8], fan[9], 0, make_boundary_point(begin[2], begin[3], begin[4], begin[5], begin[0]));
               set_boundary_point(fan[10], fan[11], 0, make_boundary_point(begin[3], begin[4], begin[5], begin[0], begin[1]));
          } else {
               auto omega = 3.0f / 8.0f + std::cos(glm::two_pi<hpreal>() / valence) / 4.0f;
               omega = (5.0f / 8.0f) - omega * omega;
               omega = 3.0f * valence / (8.0f * omega);
               corner *= omega;
               for(auto& vertex : boost::make_iterator_range(begin, end)) corner += vertex.position;
               corner /= (valence + omega);

               auto f = std::begin(fan);
               auto delta = glm::two_pi<hpreal>() / valence;
               for(auto middle = begin; middle != end; ++middle, f += 2) {
                    auto theta = 0.0f;
                    auto tangent = Vector(0.0);
                    auto update_tangent = [&](auto& vertex) {
                         tangent += (2.0f / valence) * std::cos(theta) * vertex.position;
                         theta += delta;
                    };
                    std::for_each(middle, end, update_tangent);
                    std::for_each(begin, middle, update_tangent);//TODO: instead of recomputing the tagent, simply rotate the first one
                    tangent = glm::normalize(tangent);
                    auto vector = get_boundary_point(surface, t, i, 1) - corner;
                    auto r = std::fmin(glm::length2(vector) / std::abs(glm::dot(tangent, vector)), glm::length(vector)) / 2.0f;
                    set_boundary_point(f[0], f[1], 0, corner + r * tangent);
               }
          }
          surface.setCorner(t, i, corner);
          visit_pairs(fan, [&](auto u, auto j) { surface.setCorner(u, j, t, i); });

          auto set_interior_point = [&](auto t, auto i, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
               auto point = (1.0f / 24.0f) * (10.0f * center.position + vertex0.position + 6.0f * vertex1.position + 6.0f * vertex2.position + vertex3.position);
               surface.setInteriorPoint(t, i, point);
          };

          set_interior_point(fan[0], fan[1], begin[valence - 1], begin[0], begin[1], begin[2]);
          if(valence > 3) {
               auto middle = begin;
               visit_pairs(std::begin(fan) + 2, valence - 3, 2, [&](auto u, auto j) {
                    set_interior_point(u, j, middle[0], middle[1], middle[2], middle[3]);
                    ++middle;
               });
          }
          auto n = (valence << 1) - 4;
          set_interior_point(fan[n], fan[n + 1], begin[valence - 3], begin[valence - 2], begin[valence - 1], begin[0]);
          set_interior_point(fan[n + 2], fan[n + 3], begin[valence - 2], begin[valence - 1], begin[0], begin[1]);
     }

     return surface;
}

}//namespace happah
