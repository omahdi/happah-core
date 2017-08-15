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
     static constexpr hpindex o[3] = { 1, 2, 0 };

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
          static constexpr hpuint o1[3] = { 1, 2, 0 };
          static constexpr hpuint o2[3] = { 2, 0, 1 };

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
void draw_fundamental_domain(const TriangleGraph<Vertex>& graph, const std::vector<Point2D> & positions) {
     auto & edges = graph.getEdges();
     auto & vertices = graph.getVertices();  
     auto indices = make_indices(graph);
     std::ofstream file("fundamental_domain.off", std::ios::out | std::ios::trunc);
     if (file) {
          file<<"OFF"<<std::endl;
	  file<<graph.getVertices().size()<<" "<<graph.getNumberOfTriangles()<<" 0"<<std::endl;
          for (int i = 0; i < vertices.size(); i++) { 
	       file<<positions[i].x;
	       file<<" ";
	       file<<positions[i].y;
	       file<<" ";
               file<<1<<std::endl;
	  }
          visit_triplets(indices, [&](auto i0, auto i1, auto i2) { file << "3 " << i0 << ' ' << i1 << ' ' << i2 << '\n'; });
	  file.close();
     } else {
          std::cout<<"problem when oppening the file."<<std::endl;
     }
}

template<class Vertex>
std::vector<Point2D> make_fundamental_domain(const TriangleGraph<Vertex>& graph, const CirclePacking & packing, const Indices cut, const hpindex edge0 = 0u) {
     auto positions = std::vector<Point2D>(graph.getNumberOfVertices()); 
     auto & edges = graph.getEdges();
     auto & vertices = graph.getVertices();
     auto lengths = make_lengths(packing, graph); 
     auto flattened = boost::dynamic_bitset<>(graph.getNumberOfEdges(), false);
     auto flattened_vertices = boost::dynamic_bitset<>(graph.getNumberOfVertices(),false);
     auto vertices_in_cut = boost::dynamic_bitset<>(graph.getNumberOfVertices(),false);
     
     auto in_cut = boost::dynamic_bitset<>(graph.getNumberOfEdges(), false);
     for(auto& e : cut) {
          in_cut[e] = true;
	  vertices_in_cut[edges[e].vertex] = true;
     }

     auto & radii = packing.getRadii();
     auto & weights = packing.getWeights(); 

     auto e0 = edges[edge0];
     auto edge1 = e0.next;
     auto e1 = edges[edge1];
     auto edge2 = e0.previous;
     auto e2 = edges[edge2];
     
     auto vertex0 = e0.vertex;
     auto vertex1 = e1.vertex;
     auto vertex2 = e2.vertex;
    
     auto l01 = lengths[edge1]; 
     auto l02 = lengths[edge0];
     auto l12 = lengths[edge2];
     assert(l01 + l12 > l02 && l01 + l02 > l12 && l12 + l02 > l01);
     auto theta_0 = std::acos(((std::cosh(l01)*std::cosh(l02))-std::cosh(l12))/(std::sinh(l01)*std::sinh(l02)));
     positions[vertex0].x = 0;
     positions[vertex0].y = 0;
     auto coeff1 = (std::exp(l01)-1)/(std::exp(l01)+1);
     positions[vertex1].x = coeff1;
     positions[vertex1].y = 0; 
     auto coeff2 = (std::exp(l02)-1)/(std::exp(l02)+1);
     positions[vertex2].x = coeff2*std::cos(theta_0); 
     positions[vertex2].y = coeff2*std::sin(theta_0);
     

     flattened[edge0] = true;
     flattened[edge1] = true;
     flattened[edge2] = true;
     flattened_vertices[vertex0] = true;
     flattened_vertices[vertex1] = true;
     flattened_vertices[vertex2] = true;
     std::stack<hpindex> todo; 
     if (!in_cut[edge0]) todo.push(e0.opposite);
     if (!in_cut[edge1]) todo.push(e1.opposite);
     if (!in_cut[edge2]) todo.push(e2.opposite);
     int nb = 0;
     int nb_edges = 0;
     for (int i = 0; i < edges.size(); i++) {
	  if (flattened[i]) nb_edges++;
     }
     std::cout<<"nb_edges = "<<nb_edges<<std::endl;
     int cptr = 0;
     std::cout<<"size todo = "<<todo.size()<<std::endl;
     while (!todo.empty()) {
          auto edge0 = todo.top();
          todo.pop();
          auto e0 = edges[edge0];
	  auto edge1 = e0.next;
          auto e1 = edges[edge1];
	  auto edge2 = e0.previous;
          auto e2 = edges[edge2];
          if (flattened[edge0] && flattened[edge1] && flattened[edge2]) continue;
          assert(flattened_vertices[e0.vertex]);
          assert(flattened_vertices[edges[e0.opposite].vertex]);
          auto v0 = e0.vertex;
          auto v1 = e1.vertex; 
          if (!vertices_in_cut[v1] && flattened_vertices[v1]) continue;
          auto v2 = e2.vertex;
          assert(flattened_vertices[v0] && flattened_vertices[v2]);
          auto c0 = positions[v0];
          auto r0 = lengths[edge1];
          auto c1 = positions[v2];
     	  auto r1 = lengths[edge2];
          auto circle0 = make_circle(c0,r0);
          auto circle1 = make_circle(c1,r1);
          auto circle0_eucl = poincare_to_euclidian(circle0);
          auto circle1_eucl = poincare_to_euclidian(circle1);
          auto intersections = intersect_circle(circle0_eucl, circle1_eucl);
          Point2D A;
          A.x = positions[v2].x - positions[v0].x;
          A.y = positions[v2].y - positions[v0].y;
          glm::vec3 Avec = glm::vec3(A.x, A.y,0);
          Point2D B;
          B.x = std::get<0>(intersections).x - positions[v0].x;
          B.y = std::get<0>(intersections).y - positions[v0].y;
          glm::vec3 Bvec = glm::vec3(B.x, B.y,0);
          Point2D C;
          C.x = std::get<1>(intersections).x - positions[v0].x;
          C.y = std::get<1>(intersections).y - positions[v0].y;
          glm::vec3 Cvec = glm::vec3(C.x, C.y,0);
          if (glm::cross(Avec,Bvec)[2] > 0) {
	       positions[v1] = std::get<0>(intersections);
          } else {
               if (glm::cross(Avec,Cvec)[2] < 0) {
		    nb++;
		    if (glm::cross(Avec,Cvec)[2] < glm::cross(Avec,Bvec)[2] ) {
			 positions[v1] = std::get<0>(intersections);
		    } else {
			 positions[v1] = std::get<1>(intersections);
		    }
               }	
               positions[v1] = std::get<1>(intersections);
          }
          flattened[edge0] = true;  
          flattened[edge1] = true; 
          flattened[edge2] = true; 
          flattened_vertices[v1] = true;
          if (!in_cut[edge1]) todo.push(e1.opposite);
          if (!in_cut[edge2]) todo.push(e2.opposite);
	  cptr++;
     }
     nb_edges = 0;
     for (int i = 0; i < edges.size(); i++) {
	  if (flattened[i]) nb_edges++;
     }
     int nb_vertices = 0;
     for (int i = 0; i < vertices.size(); i++) {
	  if (flattened_vertices[i]) nb_vertices++;
     }
     std::cout<<"nb_vertices flattened = "<<nb_vertices<<std::endl;
     std::cout<<"nb edges = "<<graph.getNumberOfEdges()<<std::endl;
     std::cout<<"nb_edges flattened = "<<nb_edges<<std::endl;
     std::cout<<"nb = "<<nb<<std::endl;
     std::cout<<"nb vertices = "<<graph.getNumberOfVertices()<<std::endl;
     std::cout<<"cptr = "<<cptr<<std::endl;
     return positions;
}

template<class Vertex>
std::vector<hpreal> make_lengths(const CirclePacking& packing, const TriangleGraph<Vertex>& graph) {
     auto lengths = std::vector<hpreal>(graph.getNumberOfEdges());
     auto& edges = graph.getEdges();
     auto& radii = packing.getRadii();
     auto& weights = packing.getWeights();
     for(auto& edge : edges) {
          auto r0 = radii[edge.vertex];
          auto r1 = radii[edges[edge.opposite].vertex];
          auto index = make_edge_index(edge);
          lengths[index] = std::acosh(std::cosh(r0) * std::cosh(r1) - std::sinh(r0) * std::sinh(r1) * weights[index]);
     }
     return lengths;
}

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
