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

hpreal vector_product(const Point2D P1, const Point2D P2) {
     auto x1 = P1.x;
     auto y1 = P1.y;
     auto x2 = P2.x;
     auto y2 = P2.y;
     return x1*y2 - y1*x2;
}

std::vector<Point2D> intersect_circle(const Point2D & c1, const hpreal & r1, const Point2D & c2, const hpreal & r2) {
     std::vector<Point2D> result;
     result.reserve(2);
     auto x1 = c1.x;
     auto y1 = c1.y;
     auto x2 = c2.x;
     auto y2 = c2.y;
     auto d = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
     assert(d < r1+r2); // circles must intersect
     assert(d>std::abs(r1-r2)); // circles mustn't be contained within the other. 
     assert((x1 != x2) || (y1 != y2));
     if (x1 == x2) {
          auto Y = (pow(r1,2)-pow(r2,2)-pow(y1,2)+pow(y2,2))/(-2*(y1-y2));
          auto a = 1;
          auto b = -2*x1;
          auto c = pow(x1,2) +pow(Y,2) - 2*Y*y2 + pow(y2,2) - pow(r2,2);
          auto delta = pow(b,2) -4*a*c;
          assert(delta > 0);
          auto sol1_x = (-b -sqrt(delta))/(2*a) ;
          auto sol2_x = (-b +sqrt(delta))/(2*a) ;
          result[0].x = sol1_x;
          result[0].y = Y;
          result[1].x = sol2_x;
          result[1].y = Y;
     } else if (y1 == y2) {
          auto X = (pow(r1,2)-pow(r2,2)-pow(x1,2)+pow(x2,2))/(2*(x2-x1));
          auto a = 1;
          auto b = -2*y1;
          auto c = pow(X,2) - 2*X*x2 + pow(x2,2) + pow(y1,2) - pow(r2,2);
          auto delta = pow(b,2) -4*a*c;
          assert(delta > 0);
          auto sol1_y = (-b -sqrt(delta))/(2*a) ;
          auto sol2_y = (-b +sqrt(delta))/(2*a) ;
          result[0].x = X;
          result[0].y = sol1_y;
          result[1].x = X;
          result[1].y = sol2_y;
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
          result[0].x = sol1_x;
          result[0].y = sol1_y;
          result[1].x = sol2_x;
          result[1].y = sol2_y;
     }
     return result;
}

template<class Vertex>
void draw_fundamental_domain(const TriangleGraph<Vertex>& graph, const std::vector<Point3D> & positions) {
     auto & edges = graph.getEdges();
     auto & vertices = graph.getVertices();  
     auto& indices = graph.getIndices();
     std::ofstream file("fundamental_domain.off", std::ios::out | std::ios::trunc);
     if (file) {
          file<<"OFF"<<std::endl;
	  file<<graph.getVertices().size()<<" "<<graph.getNumberOfTriangles()<<" 0"<<std::endl;
          for (int i = 0; i < vertices.size(); i++) { 
	       file<<positions[i].x;
	       file<<" ";
	       file<<positions[i].y;
	       file<<" ";
               file<<positions[i].z<<std::endl;
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
     assert(dhyp != deucl);
     auto temp2 = pow(deucl,2) - (pow(dhyp,2)-pow(u,2))/(1-pow(u,2)*pow(dhyp,2));
     r = sqrt(temp2);
}

//v0 and v1 are already flattened, v2 is to flatten
template<class Vertex>
void flatten_third_point(const TriangleGraph<Vertex>& graph, const CirclePackingMetric & metric, std::vector<Point2D> & positions, const std::vector<hpreal> & lengths, hpindex e0, hpindex e1, hpindex e2, int & nb) { 
     auto & edges = graph.getEdges();    
     auto v0 = edges[e0].vertex;
     auto v1 = edges[e1].vertex; 
     auto v2 = edges[e2].vertex;
     auto c0 = positions[v0];
     auto r0 = lengths[e0];
     auto c1 = positions[v1];
     auto r1 = lengths[e2];
     poincare_to_euclidian(c0,r0);
     poincare_to_euclidian(c1,r1);
     auto intersections = intersect_circle(c0, r0, c1, r1);
     Point2D A;
     A.x = positions[v1].x - positions[v0].x;
     A.y = positions[v1].y - positions[v0].y;
     Point2D B;
     B.x = intersections[0].x - positions[v0].x;
     B.y = intersections[0].y - positions[v0].y;
     Point2D C;
     C.x = intersections[1].x - positions[v0].x;
     C.y = intersections[1].y - positions[v0].y;
     if (vector_product(A,B) > 0) {
	  positions[v2] = intersections[0];
     } else {
          //assert(vector_product(A,C) >0);
          if (vector_product(A,C) < 0) {	
               nb++;
	       if (vector_product(A,B) > vector_product(A,C)) {
		    positions[v2] = intersections[0]; 
	       } else {
		    positions[v2] = intersections[1];
	       }
	  } else {
               positions[v2] = intersections[1];
	  }
     }
}

template<class Vertex>
std::vector<Point3D> flatten_fundamental_domain(const TriangleGraph<Vertex>& graph, const CirclePackingMetric & metric) {
     auto positions = std::vector<Point2D>(); 
     auto & edges = graph.getEdges();
     auto & vertices = graph.getVertices();
     positions.reserve(graph.getNumberOfVertices());
     auto lengths = make_length(graph, metric ); 
     auto flattened = boost::dynamic_bitset<>(graph.getVertices().size(), false);
     hpindex edge0 = 0u;
     flatten_seed_face(graph,positions,edge0, metric, lengths);
     flattened[edges[edge0].vertex] = true;
     flattened[edges[edges[edge0].next].vertex] = true;
     flattened[edges[edges[edge0].previous].vertex] = true;
     std::stack<hpindex> todo; 
     todo.push(edges[edge0].opposite);
     todo.push(edges[edges[edge0].next].opposite);
     todo.push(edges[edges[edge0].previous].opposite);
     int nb = 0;
     while (!todo.empty()) {
          auto e0 = todo.top();
          todo.pop();
          auto edge0 = edges[e0];
	  auto e1 = edge0.next;
	  auto e2 = edge0.previous;
          auto v0 = edges[e0].vertex;
          auto v1 = edges[e1].vertex;
	  auto v2 = edges[e2].vertex;
          if (flattened[v0] && flattened[v1] && flattened[v2]) continue;
	  if (!flattened[v1]) {
	       assert(!flattened[v1] && flattened[v0] && flattened[v2]);
               flatten_third_point(graph, metric, positions, lengths, e2, e0, e1, nb);
               flattened[v1] = true;
               todo.push(edges[e2].opposite);
               todo.push(edges[e1].opposite);
          } else if (!flattened[v0]) {
	       assert(!flattened[v0] && flattened[v2] && flattened[v1]);
               flatten_third_point(graph, metric, positions, lengths, e1, e2, e0, nb);
               flattened[v0] = true;
               todo.push(edges[e1].opposite);
               todo.push(edges[e0].opposite);
          } else {
               assert(!flattened[v2] && flattened[v1] && flattened[v0]);
               flatten_third_point(graph, metric, positions, lengths, e0, e1, e2, nb);
               flattened[v2] = true;
               todo.push(edges[e0].opposite);
               todo.push(edges[e2].opposite);
          }
     }
     std::cout<<"nb vertices = "<<vertices.size()<<std::endl;
     std::cout<<"nb = "<<nb<<std::endl;
     auto positions3D = std::vector<Point3D>(); 
     positions3D.reserve(graph.getNumberOfVertices());
     for (int i = 0; i < vertices.size(); i++) {
	  positions3D[i].x = positions[i].x;
	  positions3D[i].y = positions[i].y;
	  positions3D[i].z = 1;
     }
     return positions3D;
}

template<class Vertex>
void flatten_seed_face(const TriangleGraph<Vertex>& graph, std::vector<Point2D> & positions, const hpindex edge0, const CirclePackingMetric & metric, const std::vector<hpreal>& lengths) {
     auto & vertices = graph.getVertices();  
     auto & edges = graph.getEdges();
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
               visit_spokes(graph.getEdges(), graph.getOutgoing(++v), [&](auto e) {
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

}//namespace happah

