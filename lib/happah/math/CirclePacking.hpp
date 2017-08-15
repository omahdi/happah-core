// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <glm/gtc/constants.hpp>

#include "happah/geometry/Circle.hpp"
#include "happah/geometry/TriangleGraph.hpp"
#include "happah/math/Space.hpp"
#include "happah/util/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

class CirclePacking;

hpreal angle_sum(const CirclePacking& packing, const Indices& neighbors, hpindex t, hpindex i);

inline hpreal length(const CirclePacking& packing, hpindex t, hpindex i);

inline CirclePacking make_circle_packing(std::vector<hpreal> radii, std::vector<hpreal> weights, Indices indices);

CirclePacking make_circle_packing(std::vector<hpreal> weights, Indices indices, const Indices& neighbors, hpreal epsilon = EPSILON);

inline Indices make_neighbors(const CirclePacking& packing);

template<class Vertex = VertexP2, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const CirclePacking& packing, const Indices& neighbors, const Indices& border, hpindex t, VertexFactory&& build = VertexFactory());

hpreal validate(const CirclePacking& packing, const Indices& neighbors);

//DEFINITIONS

class CirclePacking {
public:
     CirclePacking(std::vector<hpreal> radii, std::vector<hpreal> weights, Indices indices)
          : m_indices(std::move(indices)), m_radii(std::move(radii)), m_weights(std::move(weights)) {}

     const Indices& getIndices() const { return m_indices; }

     const std::vector<hpreal>& getRadii() const { return m_radii; }

     hpreal getRadius(hpindex t, hpindex i) const { return m_radii[m_indices[3 * t + i]]; }

     hpreal getWeight(hpindex e) const { return m_weights[e]; }

     const std::vector<hpreal>& getWeights() const { return m_weights; }

     void setRadius(hpindex t, hpindex i, hpreal r) { m_radii[m_indices[3 * t + i]] = r; }

private:
     Indices m_indices;
     std::vector<hpreal> m_radii;
     std::vector<hpreal> m_weights; 

};//CirclePacking

inline hpreal length(const CirclePacking& packing, hpindex t, hpindex i) {
     constexpr hpindex o[3] = { 1, 2, 0 };

     auto r0 = packing.getRadius(t, i);
     auto r1 = packing.getRadius(t, o[i]);
     return std::acosh(std::cosh(r0) * std::cosh(r1) - std::sinh(r0) * std::sinh(r1) * packing.getWeight(3 * t + i));
}

inline CirclePacking make_circle_packing(std::vector<hpreal> radii, std::vector<hpreal> weights, Indices indices) { return { std::move(radii), std::move(weights), std::move(indices) }; }

inline Indices make_neighbors(const CirclePacking& packing) { return make_neighbors(packing.getIndices()); }

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const CirclePacking& packing, const Indices& neighbors, const Indices& border, hpindex t, VertexFactory&& build) {
     auto l0 = length(packing, t, 0);
     auto l1 = length(packing, t, 1);
     auto l2 = length(packing, t, 2);
     auto temp = (std::cosh(l0) * std::cosh(l2) - std::cosh(l1)) / (std::sinh(l0) * std::sinh(l2));
     auto x1 = std::exp(l0);
     auto x2 = std::exp(l2);

     x1 = (x1 - 1) / (x1 + 1);
     x2 = (x2 - 1) / (x2 + 1);

     return make_triangle_mesh(neighbors, border, t, build(Point2D(0, 0)), build(Point2D(x1, 0)), build(x2 * Point2D(temp, std::sqrt(1.0 - temp * temp))), [&](auto t, auto i, auto& vertex0, auto& vertex1, auto& vertex2) {
          constexpr hpuint o1[3] = { 1, 2, 0 };
          constexpr hpuint o2[3] = { 2, 0, 1 };

          auto u = make_neighbor_index(neighbors, t, i);
          auto j = make_neighbor_offset(neighbors, u, t);
          auto circle0 = poincare_to_euclidean(make_circle(vertex0.position, length(packing, u, o2[j])));
          auto circle1 = poincare_to_euclidean(make_circle(vertex1.position, length(packing, u, o1[j])));
          assert(intersect(circle0, circle1) != boost::none);
          auto intersections = *intersect(circle0, circle1);
          if(auto intersection = boost::get<Point2D>(&intersections)) return build(*intersection);
          else return build(std::get<0>(boost::get<std::tuple<Point2D, Point2D> >(intersections)));
     });
}

//WORKSPACE

template<class Vertex>
CirclePacking make_circle_packing(const TriangleGraph<Vertex>& graph, hpreal threshold = 0.05, hpreal epsilon = 0.01) {
     auto lengths = std::vector<hpreal>(graph.getNumberOfEdges());
     auto weights = std::vector<hpreal>(graph.getNumberOfEdges());
     auto radii = std::vector<hpreal>(graph.getNumberOfVertices(), 0);
     auto e = 0u, v = 0u;

     e = hpindex(-1);
     for(auto& length : lengths) {
          auto& edge0 = graph.getEdge(++e);
          auto& edge1 = graph.getEdge(edge0.opposite);
          auto& vertex0 = graph.getVertex(edge0.vertex);
          auto& vertex1 = graph.getVertex(edge1.vertex);
          length = glm::length(vertex1.position - vertex0.position);
     };

     v = hpindex(-1);
     for(auto& radius : radii) {
          auto valence = 0u;
          visit_spokes(make_spokes_enumerator(graph.getEdges(), graph.getOutgoing(++v)), [&](auto e) {
               auto& edge = graph.getEdge(e);
               radius += (lengths[e] + lengths[edge.previous] - lengths[edge.next]) / 2.0;
               ++valence;
          });
          radius /= valence;
     }

     e = hpindex(-1);
     for(auto& weight : weights) {
          auto& edge0 = graph.getEdge(++e);
          auto& edge1 = graph.getEdge(edge0.opposite);
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
               auto& edge0 = graph.getEdge(++e);
               auto& edge1 = graph.getEdge(edge0.opposite);
               auto r0 = radii[edge0.vertex];
               auto r1 = radii[edge1.vertex];
               length = std::acosh(std::cosh(r0) * std::cosh(r1) - std::sinh(r0) * std::sinh(r1) * weights[e]);
               assert(length > 0.0);
          }

          v = hpindex(-1);
          auto max = std::numeric_limits<hpreal>::min();
          for(auto& radius : radii) {
               auto curvature = glm::two_pi<hpreal>();
               visit_spokes(make_spokes_enumerator(graph.getEdges(), graph.getOutgoing(++v)), [&](auto e) {
                    auto& edge = graph.getEdge(e);
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

     return make_circle_packing(radii, weights, make_indices(graph));
}

template<class Vertex>
CirclePacking make_circle_packing_euclidean(const TriangleGraph<Vertex>& graph, hpreal threshold = 0.05, hpreal epsilon = 0.01) {
     auto lengths = std::vector<hpreal>(graph.getNumberOfEdges());
     auto weights = std::vector<hpreal>(graph.getNumberOfEdges());
     auto radii = std::vector<hpreal>(graph.getNumberOfVertices(), 0);
     auto e = 0u, v = 0u;

     e = hpindex(-1);
     for(auto& length : lengths) {
          auto& edge0 = graph.getEdge(++e);
          auto& edge1 = graph.getEdge(edge0.opposite);
          auto& vertex0 = graph.getVertex(edge0.vertex);
          auto& vertex1 = graph.getVertex(edge1.vertex);
          length = glm::length(vertex1.position - vertex0.position);
     };

     v = hpindex(-1);
     for(auto& radius : radii) {
          auto valence = 0u;
          visit_spokes(graph.getEdges(), graph.getOutgoing(++v), [&](auto e) {
               auto& edge = graph.getEdge(e);
               radius += (lengths[e] + lengths[edge.previous] - lengths[edge.next]) / 2.0;
               ++valence;
          });
          radius /= valence;
     }

     e = hpindex(-1);
     for(auto& weight : weights) {
          auto& edge0 = graph.getEdge(++e);
          auto& edge1 = graph.getEdge(edge0.opposite);
          auto r0 = radii[edge0.vertex];
          auto r1 = radii[edge1.vertex];
          auto l = lengths[e];
          weight = (r0 + r1 < l || r0 + l < r1 || r1 + l < r0) ? 0.0 : (pow(r0,2)+pow(r1,2)-pow(l,2))/(2*r0*r1);
          if(weight < 0) weight = 0;
          assert(weight >= 0.0 && weight <= 1.0);
     }

     auto todo = false;
     auto counter = 100;
     do {
          todo = false;

          e = hpindex(-1);
          for(auto& length : lengths) {
               auto& edge0 = graph.getEdge(++e);
               auto& edge1 = graph.getEdge(edge0.opposite);
               auto r0 = radii[edge0.vertex];
               auto r1 = radii[edge1.vertex];
               length = sqrt(pow(r0,2)+pow(r1,2)-2*r0*r1*weights[e]);
               assert(length > 0.0);
          }

          v = hpindex(-1);
          auto max = std::numeric_limits<hpreal>::min();
          auto S = 0;
          for(auto& radius : radii) {
               auto curvature = glm::two_pi<hpreal>();
               visit_spokes(make_spokes_enumerator(graph.getEdges(), graph.getOutgoing(++v)), [&](auto e) {
                    auto& edge = graph.getEdge(e);
                    auto l0 = lengths[e];
                    auto l1 = lengths[edge.previous];
                    auto l2 = lengths[edge.next];
                    assert(l0 + l1 > l2 && l0 + l2 > l1 && l1 + l2 > l0);
                    curvature -= std::acos((std::cosh(l0) * std::cosh(l1) - std::cosh(l2)) / (std::sinh(l0) * std::sinh(l1)));
               });
               radius = std::log(radius);
               radius -= epsilon * curvature;
	       S += radius;
               radius = std::exp(radius);
               max = std::max(std::abs(curvature), max);
               todo |= std::abs(curvature) > threshold;
          }
          std::cout<<"S = "<<S<<std::endl;
          assert(S == 0);
          if(counter == 100) {
               counter = 1;
               std::cout << max << '\n';
          } else ++counter;
     } while(todo);

     return make_circle_packing(radii, weights, make_indices(graph));
}

}//namespace happah

