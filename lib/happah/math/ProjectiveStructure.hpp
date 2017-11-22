// Copyright 2017
//   Obada Mahdi    -                                       - omahdi@gmail.com
//   Pawel Herman   - Karlsruhe Institute of Technology     - pherman@ira.uka.de
//   Louis Lutzweiler                                       - louis.lutzweiler@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.10 - Louis Lutzweiler - Added make_projective_structure from a genus zero triangle mesh.

#pragma once

#include <vector>

#include "happah/Happah.hpp"
#include "happah/geometry/TriangleGraph.hpp"
#include "happah/geometry/TriangleMesh.hpp"
#include "happah/geometry/Vertex.hpp"
#include "happah/util/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

class ProjectiveStructure;

/// Constructs a convex hyperbolic polygon with given interior angles
/// \p thetas and returns the Euclidean coordinates of its corner vertices in the
/// Poincare model.
///
/// Angles \p thetas must be in [0,pi).
///
/// See Theorem 7.16.2 in [1, p. 155f] for the mathematics.
///
/// \note [1]: Beardon, "The Geometry of Discrete Groups"
std::vector<Point2D> make_convex_polygon(const std::vector<hpreal>& angles, hpreal epsilon = EPSILON);

std::vector<Point2D> make_convex_polygon(const Indices& valences, hpreal epsilon = EPSILON);

inline ProjectiveStructure make_projective_structure(Indices neighbors, std::vector<hpreal> transitions);

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleMesh<Vertex>& mesh, const Point3D& center, const Indices& neighbors);

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph);

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph, const Indices& cut);

ProjectiveStructure make_projective_structure(const Indices& valences, const Indices& pairings);

std::vector<Point3D> make_sun(const Indices& valences, const Indices& pairings);

//NOTE: Border has to be sorted.
template<class Vertex = VertexP3, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const ProjectiveStructure& structure, const Indices& border, hpindex t, const Point3D& p0, const Point3D& p1, const Point3D& p2, VertexFactory&& factory = VertexFactory());

hpreal validate(const ProjectiveStructure& structure);

//DEFINITIONS

class ProjectiveStructure {
public:
     ProjectiveStructure(Indices neighbors, std::vector<hpreal> transitions)
          : m_neighbors(std::move(neighbors)), m_transitions(std::move(transitions)) {}

     const Indices& getNeighbors() const { return m_neighbors; }

     const std::vector<hpreal>& getTransitions() const { return m_transitions; }

private:
     Indices m_neighbors;
     std::vector<hpreal> m_transitions;

};//ProjectiveStructure

inline ProjectiveStructure make_projective_structure(Indices neighbors, std::vector<hpreal> transitions) { return { std::move(neighbors), std::move(transitions) }; }

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleMesh<Vertex>& mesh, const Point3D& center, const Indices& neighbors) {
     auto transitions = std::vector<hpreal>(9 * size(mesh));
     auto t = hpindex(0);

     visit_triplets(mesh.getIndices(), [&](auto i0, auto i1, auto i2) {
          auto p0 = mesh.getVertex(i0).position - center;
          auto p1 = mesh.getVertex(i1).position - center;
          auto p2 = mesh.getVertex(i2).position - center;
          
          auto make_transition = [&](auto i, auto& p0, auto& p1, auto& p2) {
               static constexpr hpuint o[3] = { 2, 0, 1 };

               auto u = make_neighbor_index(neighbors, t, i);
               auto j = make_neighbor_offset(neighbors, u, t);
               auto matrix = glm::inverse(hpmat3x3(p0, p1, p2));
               auto p3 = mesh.getVertex(u, o[j]).position - center;
               auto transition = matrix * p3;
               auto temp = std::begin(transitions) + (9 * u + 3 * j);

               temp[0] = transition.x;
               temp[1] = transition.y;
               temp[2] = transition.z;
          };
     
          make_transition(0, p0, p2, p1);
          make_transition(1, p1, p0, p2);
          make_transition(2, p2, p1, p0);

          ++t;
     });
     
     return ProjectiveStructure(neighbors, transitions);
}

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph) { return make_projective_structure(graph, undegenerate(graph, trim(graph, cut(graph)))); }
//ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph) { return make_projective_structure(graph, trim(graph, cut(graph))); }

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph, const Indices& cut, std::vector<Point3D> &test, std::vector<Point3D> &triangle, hpuint& hhh) {
     using happah::format::hph::operator<<;

     auto& edges = graph.getEdges();
     auto analysis = analyze(graph, cut);
     auto& valences = std::get<0>(analysis);
     auto& indices = std::get<1>(analysis);
     auto& pairings = std::get<2>(analysis);
     auto lengths = std::vector<hpuint>();
     auto transitions = std::vector<hpreal>();
     auto sun = make_sun(valences, pairings);

     transitions.reserve(9 * size(graph));
     lengths.reserve(valences.size());
     for(auto i = std::begin(indices), end = std::end(indices) - 1; i != end; ++i) lengths.push_back(*(i + 1) - *i - 1);
     lengths.push_back(cut.size() - indices.back() + indices.front() - 1);
     
     auto polyline = std::vector<Point3D>();
     auto j = std::begin(sun);
     auto* point0 = &sun[0];
     
     auto do_parametrize = [&](auto& point0, auto& point1, auto m) {
          auto delta = (hpreal(1.0) / hpreal(m + 1)) * (point1 - point0);
          auto point = point0;

          polyline.emplace_back(point.x, point.y, point.z);
          while(m--) {
               point += delta;
               polyline.emplace_back(point.x, point.y, point.z);
          }
     };

     for(auto i : boost::make_iterator_range(std::begin(lengths), std::end(lengths) - 1)) {
          auto& point1 = *(++j);

          do_parametrize(*point0, point1, i);
          point0 = &point1;
     }
     do_parametrize(sun[valences.size() - 1], sun[0], lengths.back());
     
     auto polygon = std::vector<Point2D>();
     polygon.reserve(polyline.size());
     for(auto& point : polyline) polygon.emplace_back(point.x / point.z, point.y / point.z);
     
     auto interior = parametrize(graph, cut, polygon);
     auto points = std::vector<Point3D>();
     auto p = Indices(graph.getNumberOfVertices(), hpuint(0));
     auto n = hpuint(0);

     for(auto p : polyline) std::cout << p << std::endl;
     std::cout << "----------------------" << std::endl;
     for(auto p : polygon) std::cout << p << std::endl;
     //std::cout << "before: ";
     //for(auto& point : polyline) std::cout << point << ' ';
     //std::cout << '\n';

     //std::cout << "cut" << cut.size() << "polygon" << polygon.size() << std::endl;
     assert(cut.size() == polygon.size());
     for(auto& e : cut) p[edges[e].vertex] = std::numeric_limits<hpuint>::max();
     for(auto& v : p) if(v == hpuint(0)) v = n++;
     points.reserve(n);

     auto angles = std::vector<hpreal>();
     angles.reserve(polyline.size());
     auto min_angle = std::atan2(polyline[0].y, polyline[0].x);
     auto t = hpuint(0);
     for(auto& point : polyline) {
          auto a = std::atan2(polyline[t].y, polyline[t].x);
          a = a < min_angle ? a + 2 * M_PI : a;
          angles.push_back(a);
          ++t;
     }
     
     auto cache = std::vector<hpreal>();
     cache.reserve(4 * polyline.size());
     auto center = Point3D(0, 0, 1);
     t = 0;
     for(auto& point : polyline) {
          auto q = polyline[(t + 1) % polyline.size()];
          auto n = glm::cross(q - center, polyline[t] - center);
          cache.push_back(glm::dot(n, center));
          cache.push_back(n.x);
          cache.push_back(n.y);
          cache.push_back(n.z);
          ++t;
     }
     
     //project interior to mesh
     for(auto& point : interior) {
          auto a = std::atan2(point.y, point.x);
          a = a < min_angle ? a + 2 * M_PI : a;
          auto t = std::lower_bound(angles.begin(), angles.end(), a) - angles.begin() - 1;
          auto lambda = cache[4 * t] / (cache[4 * t + 1] * point.x + cache[4 * t + 2] * point.y + cache[4 * t + 3]);
          points.push_back(lambda * Point3D(point, 1));
          //std::cout << "triangle #" << t << std::endl;
          //std::cout << "lambda=" << lambda << std::endl;
          //std::cout << "projecting " << point.x << "," << point.y << " to " << lambda * point.x << "," << lambda * point.y << "," << lambda << std::endl;
     }
     

     auto get_point = [&](auto& e) {
          auto w = make_spokes_walker(edges, e); 
          auto temp = std::end(cut);
          while(temp == std::end(cut)) temp = std::find(std::begin(cut), std::end(cut), edges[*(++w)].opposite);
          auto n = std::distance(std::begin(cut), temp);
          n = (n < indices.front()) ? size(cut) - indices.front() + n : n - indices.front();
          //n = (n - indices.front()) % size(indices);
          //n = (n + indices.front()) % size(polyline);
          return polyline[n];
     };

     assert(polyline.size() == cut.size());
     
     //std::cout << "after: ";
     //for(auto& point : polyline) std::cout << point << ' ';
     //std::cout << '\n';
     auto found_first_inner_triangle = false;
     for(auto& edge : edges) {
          auto& i0 = edge.next;
          auto& edge1 = edges[edge.opposite];
          auto& edge2 = edges[i0];
          auto& edge3 = edges[edge1.next];
          auto& p0 = edge.vertex;
          auto& p1 = edge1.vertex;
          auto& p2 = edge2.vertex;
          auto& p3 = edge3.vertex;
          auto c0 = p[p0] == std::numeric_limits<hpuint>::max();
          auto c1 = p[p1] == std::numeric_limits<hpuint>::max();
          auto c2 = p[p2] == std::numeric_limits<hpuint>::max();
          auto c3 = p[p3] == std::numeric_limits<hpuint>::max();
          auto point0 = (c0) ? get_point(i0) : points[p[p0]];
          auto point1 = (c1) ? get_point(edge1.opposite) : points[p[p1]];
          auto point2 = (c2) ? get_point(edge.previous) : points[p[p2]];
          auto point3 = Point3D();
          
          if(!c0 && !c1 && !c2 && !found_first_inner_triangle){
               found_first_inner_triangle = true;
               hhh = edge.next / 3;
               triangle.push_back(points[p[p1]]);
               triangle.push_back(points[p[p0]]);
               triangle.push_back(points[p[p2]]);
               std::cout << "found it" << std::endl;
          }



          auto transition = Point3D(0);
          
          if(std::find(std::begin(cut), std::end(cut), make_edge_index(edge)) != std::end(cut)) {//edge is on cut
               //std::cout << "on cut" << std::endl;
               /*
	          auto polyindex = get_polyline_index(p0);
	          auto a0 = polyline[polyindex + 1];
	          auto a1 = polyline[polyindex];
               auto a2 = sun[(sun.size() >> 1) + polyindex];
	          auto b0 = polyline[pairings[polyindex + 1]];
     	     auto b1 = polyline[pairings[polyindex]];
	          auto b2 = center;

               auto c = glm::inverse(hpmat3x3(b0, b1, b2)) * ((c3) ? get_point(edge3.next) : points[p[p3]]);
               point3 = hpmat3x3(a0, a1, a2) * c;
               */
          } else {//edge is in interior
               std::cout << "in interior" << std::endl;
               point3 = (c3) ? get_point(edge3.next) : points[p[p3]];
          }

          transition = glm::inverse(hpmat3x3(point0, point3, point1)) * point2;
          std::cout << "p0: " << point0 << ' ' << c0 << '\n';
          std::cout << "p1: " << point1 << ' ' << c1 << '\n';
          std::cout << "p2: " << point2 << ' ' << c2 << '\n';
          std::cout << "p3: " << point3 << ' ' << c3 << '\n';
          std::cout << "transition " << transition.x << ", " << transition.y << ", " << transition.z << std::endl;
          std::cout << "-------------" << std::endl;
          if(std::find(std::begin(cut), std::end(cut), make_edge_index(edge)) == std::end(cut)) assert(glm::length(transition) > EPSILON);
          
          transitions.push_back(transition.x);
          transitions.push_back(transition.y);
          transitions.push_back(transition.z);
     }

     std::cout << "number transitions " << transitions.size() << '\n';
     test.insert(std::end(test), std::begin(polyline), std::end(polyline));
     test.insert(std::end(test), std::begin(points), std::end(points));
     return make_projective_structure(make_neighbors(graph), std::move(transitions));
}

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const ProjectiveStructure& structure, const Indices& border, hpindex t, const Point3D& point0, const Point3D& point1, const Point3D& point2, VertexFactory&& build) {
     auto& neighbors = structure.getNeighbors();
     auto& transitions = structure.getTransitions();
     return make_triangle_mesh(neighbors, border, t, build(point0), build(point1), build(point2), [&](auto t, auto i, auto& vertex0, auto& vertex1, auto& vertex2) {
          auto u = make_neighbor_index(neighbors, t, i);
          auto j = make_neighbor_offset(neighbors, u, t);
          auto transition = std::begin(transitions) + 3 * (3 * u + j);
          return build(transition[0] * vertex1.position + transition[1] * vertex2.position + transition[2] * vertex0.position);
     });
}

}//namespace happah

