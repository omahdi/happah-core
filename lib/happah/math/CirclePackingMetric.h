// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <glm/gtc/constants.hpp>
#include <happah/math/Space.h>
#include "happah/geometries/TriangleGraph.h"

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


std::tuple<Point2D, Point2D> intersect_circle(const Point2D & c1, const hpreal & r1, const Point2D & c2, const hpreal & r2) {
     auto x1 = c1.x;
     auto y1 = c1.y;
     auto x2 = c2.x;
     auto y2 = c2.y;
     auto d = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
     assert(d < r1+r2); // circles must intersect
     assert(d>std::abs(r1-r2)); // circles mustn't be contained within the other. 
     assert((x1 != x2) || (y1 != y2));
     Point2D result0;
     Point2D result1;
     if (x1 == x2) {
          auto Y = (pow(r1,2)-pow(r2,2)-pow(y1,2)+pow(y2,2))/(-2*(y1-y2));
          auto a = 1;
          auto b = -2*x1;
          auto c = pow(x1,2) +pow(Y,2) - 2*Y*y2 + pow(y2,2) - pow(r2,2);
          auto delta = pow(b,2) -4*a*c;
          assert(delta > 0);
          auto sol1_x = (-b -sqrt(delta))/(2*a) ;
          auto sol2_x = (-b +sqrt(delta))/(2*a) ;
          result0.x = sol1_x;
          result0.y = Y;
          result1.x = sol2_x;
          result1.y = Y;
     } else if (y1 == y2) {
          auto X = (pow(r1,2)-pow(r2,2)-pow(x1,2)+pow(x2,2))/(2*(x2-x1));
          auto a = 1;
          auto b = -2*y1;
          auto c = pow(X,2) - 2*X*x2 + pow(x2,2) + pow(y1,2) - pow(r2,2);
          auto delta = pow(b,2) -4*a*c;
          assert(delta > 0);
          auto sol1_y = (-b -sqrt(delta))/(2*a) ;
          auto sol2_y = (-b +sqrt(delta))/(2*a) ;
          result0.x = X;
          result0.y = sol1_y;
          result1.x = X;
          result1.y = sol2_y;
     } else {
          auto X = (pow(r1,2)-pow(r2,2)-pow(x1,2)+pow(x2,2)-pow(y1,2)+pow(y2,2))/(2*(x2-x1));
          auto Y = (y2-y1)/(x2-x1);
          auto a = 1 + pow(Y,2);
          auto b = 2*Y*x2 - 2*X*Y - 2*y2;
          auto c = pow(X,2) - 2*X*x2 + pow(x2,2) + pow(y2,2) - pow(r2,2);
          auto delta = pow(b,2) -4*a*c;
          assert(delta > 0);
          auto sol1_y = (-b -sqrt(delta))/(2*a) ;
          auto sol2_y = (-b +sqrt(delta))/(2*a) ;
          auto sol1_x = X - Y*sol1_y;
          auto sol2_x = X - Y*sol2_y;
          result0.x = sol1_x;
          result0.y = sol1_y;
          result1.x = sol2_x;
          result1.y = sol2_y;
     }
     return std::make_tuple(result0, result1);
}

template<class Vertex>
void draw_fundamental_domain(const TriangleGraph<Vertex>& graph, const std::vector<Point2D> & positions) {
     auto & edges = graph.getEdges();
     auto & vertices = graph.getVertices();  
     auto indices = make_indices(graph.getEdges());
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

void poincare_to_euclidian (Point2D & C, hpreal & r) {
     auto u = (std::exp(r)-1)/(std::exp(r)+1);
     auto dhyp = sqrt(pow(C.x,2)+pow(C.y,2));
     auto temp1 = (2-2*pow(u,2))/(1-pow(u,2)*pow(dhyp,2));
     C.x = temp1*C.x;
     C.y = temp1*C.y;
     auto deucl = sqrt(pow(C.x,2)+pow(C.y,2));
   //  assert(dhyp != deucl); I don't think it has to be different
     auto temp2 = pow(deucl,2) - (pow(dhyp,2)-pow(u,2))/(1-pow(u,2)*pow(dhyp,2));
     r = sqrt(temp2);
}

template<class Vertex>
std::vector<Point2D> make_fundamental_domain(const TriangleGraph<Vertex>& graph, const CirclePackingMetric & metric, const Indices cut, const hpindex edge0 = 0u) {
     auto positions = std::vector<Point2D>(graph.getNumberOfVertices()); 
     auto & edges = graph.getEdges();
     auto & vertices = graph.getVertices();
     auto lengths = make_length(graph, metric ); 
     auto flattened = boost::dynamic_bitset<>(graph.getNumberOfEdges(), false);
     auto flattened_vertices = boost::dynamic_bitset<>(graph.getNumberOfVertices(),false);
     auto vertices_in_cut = boost::dynamic_bitset<>(graph.getNumberOfVertices(),false);
     
     auto in_cut = boost::dynamic_bitset<>(graph.getNumberOfEdges(), false);
     for(auto& e : cut) {
          in_cut[e] = true;
	  vertices_in_cut[edges[e].vertex] = true;
     }

     auto & radii = metric.getRadii();
     auto & weights = metric.getWeights(); 

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
          poincare_to_euclidian(c0,r0);
          poincare_to_euclidian(c1,r1);
          auto intersections = intersect_circle(c0, r0, c1, r1);
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
std::vector<hpreal> make_length(const TriangleGraph<Vertex>& graph, const CirclePackingMetric & metric ) {
     auto lengths = std::vector<hpreal>(graph.getNumberOfEdges());
     auto & edges = graph.getEdges();
     auto & radii = metric.getRadii();
     auto & weights = metric.getWeights();
     for (auto& e : edges) {
          auto e_opp = edges[e.opposite];
          auto r0 = radii[e.vertex];
          auto r1 = radii[e_opp.vertex];
	  auto index = make_edge_index(e);
	  auto length = std::acosh(std::cosh(r0) * std::cosh(r1) - std::sinh(r0) * std::sinh(r1) * weights[index]);
          lengths[index] = length;
     }
     return lengths;
}

template<class Vertex>
CirclePackingMetric make_circle_packing_metric(const TriangleGraph<Vertex>& graph, hpreal threshold = 0.05, hpreal epsilon = 0.01) {
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

     return make_circle_packing_metric(radii, weights);
}

template<class Vertex>
CirclePackingMetric make_circle_packing_metric_euclidean(const TriangleGraph<Vertex>& graph, hpreal threshold = 0.05, hpreal epsilon = 0.01) {
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

     return make_circle_packing_metric(radii, weights);
}

}//namespace happah

