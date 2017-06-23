// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>

#include "happah/geometries/TriangleMesh.h"

namespace happah {

class CirclePackingMetric {
public:
     CirclePackingMetric(std::vector<hpreal> radii, std::vector<hpreal> weights);

     const std::vector<hpreal>& getRadii() const;

     const std::vector<hpreal>& getWeights() const;

private:
     std::vector<hpreal> m_radii;
     std::vector<hpreal> m_weights;

};//CirclePackingMetric

CirclePackingMetric make_circle_packing_metric(std::vector<hpreal>& radii, std::vector<hpreal>& weights);

template<class Vertex>
CirclePackingMetric make_circle_packing_metric(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpreal threshold = 0.05, hpreal epsilon = 0.01) {
     auto lengths = std::vector<hpreal>(mesh.getNumberOfEdges());
     auto weights = std::vector<hpreal>(mesh.getNumberOfEdges());
     auto radii = std::vector<hpreal>(mesh.getNumberOfVertices(), 0);
     auto e = 0u, v = 0u;

     e = hpindex(-1);
     for(auto& length : lengths) {
          auto& edge0 = mesh.getEdge(++e);
          auto& edge1 = mesh.getEdge(edge0.opposite);
          auto& vertex0 = mesh.getVertex(edge0.vertex);
          auto& vertex1 = mesh.getVertex(edge1.vertex);
          length = glm::length(vertex1.position - vertex0.position);
     };

     v = hpindex(-1);
     for(auto& radius : radii) {
          auto valence = 0u;
          visit_spokes(mesh.getEdges(), mesh.getOutgoing(++v), [&](auto e) {
               auto& edge = mesh.getEdge(e);
               radius += (lengths[e] + lengths[edge.previous] - lengths[edge.next]) / 2.0;
               ++valence;
          });
          radius /= valence;
     }

     e = hpindex(-1);
     for(auto& weight : weights) {
          auto& edge0 = mesh.getEdge(++e);
          auto& edge1 = mesh.getEdge(edge0.opposite);
          auto r0 = radii[edge0.vertex];
          auto r1 = radii[edge1.vertex];
          auto l = lengths[e];
          weight = (r0 + r1 < l || r0 + l < r1 || r1 + l < r0) ? 0.0 : (std::cosh(r0) * std::cosh(r1) - std::cosh(l)) / (std::sinh(r0) * std::sinh(r1));
          if(weight < 0) weight = 0;
          assert(weight >= 0.0 && weight <= 1.0);
     }

     auto todo = false;
     auto counter = 100;
     do {
          todo = false;

          e = hpindex(-1);
          for(auto& length : lengths) {
               auto& edge0 = mesh.getEdge(++e);
               auto& edge1 = mesh.getEdge(edge0.opposite);
               auto r0 = radii[edge0.vertex];
               auto r1 = radii[edge1.vertex];
               length = std::acosh(std::cosh(r0) * std::cosh(r1) - std::sinh(r0) * std::sinh(r1) * weights[e]);
               assert(length > 0.0);
          }

          v = hpindex(-1);
          auto max = std::numeric_limits<hpreal>::min();
          for(auto& radius : radii) {
               auto curvature = glm::two_pi<hpreal>();
               visit_spokes(mesh.getEdges(), mesh.getOutgoing(++v), [&](auto e) {
                    auto& edge = mesh.getEdge(e);
                    auto l0 = lengths[e];
                    auto l1 = lengths[edge.previous];
                    auto l2 = lengths[edge.next];
                    assert(l0 + l1 > l2 && l0 + l2 > l1 && l1 + l2 > l0);
                    curvature -= std::acos((std::cosh(l0) * std::cosh(l1) - std::cosh(l2)) / (std::sinh(l0) * std::sinh(l1)));
               });
               radius = std::log(std::tanh(radius / 2.0));
               radius -= epsilon * curvature;
               radius = 2.0 * std::atanh(std::exp(radius));
               max = std::max(std::abs(curvature), max);
               todo |= std::abs(curvature) > threshold;
          }
          if(counter == 100) {
               counter = 1;
               std::cout << max << '\n';
          } else ++counter;
     } while(todo);

     return make_circle_packing_metric(radii, weights);
}

}//namespace happah

