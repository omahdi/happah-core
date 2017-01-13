// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/range/iterator_range.hpp>
#include <cmath>
#include <glm/gtc/constants.hpp>

#include "happah/geometries/TriangleMesh.h"
#include "happah/geometries/SurfaceSplineBEZ.h"

namespace happah {

//Convert a closed(!) triangle mesh into a quartic polynomial spline.
template<class Vertex>
auto make_quartic_surface_spline_bez(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpreal epsilon = EPSILON) {
     using Space = typename Vertex::SPACE;
     using Vector = typename Space::VECTOR;

     auto surface = QuarticSurfaceSplineBEZ<Space>(mesh.getNumberOfTriangles());

     auto set_boundary_point = [&](auto t, auto i, auto k, auto&& point) {
          surface.setBoundaryPoint(t, i, k, point);
          auto u = make_neighbor_index(mesh, t, i);
          auto j = make_neighbor_offset(mesh, u, t);
          surface.setBoundaryPoint(u, j, 2 - k, t, i);
     };

     visit_diamonds(mesh, [&](auto t, auto i, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) { set_boundary_point(t, i, 1, (1.0f / 6.0f) * (2.0f * vertex0.position + vertex1.position + 2.0f * vertex2.position + vertex3.position)); });

     visit_rings(mesh, [&](auto t, auto i, auto begin, auto end) {
          auto& center = mesh.getVertex(t, i);
          auto fan = make_fan(mesh, t, i);
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
                    auto vector = surface.getBoundaryPoint(t, i, 1) - corner;
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
     });

     return surface;
}

}//namespace happah

