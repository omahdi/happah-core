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

//NOTE: Border has to be sorted.
template<class Vertex = VertexP3, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const ProjectiveStructure& structure, const Indices& border, hpindex t, const Point3D& p0, const Point3D& p1, const Point3D& p2, VertexFactory&& factory = VertexFactory());

hpreal validate(const ProjectiveStructure& structure);

namespace detail {

std::tuple<std::vector<Point2D>, hpreal> make_sun(const Indices& valences);

}//namespace detail

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
     auto transitions = std::vector<hpreal>();
     auto t = hpindex(0);
     auto n = std::begin(neighbors) - 1;

     auto make_transition = [&](auto u, auto& point0, auto& point2, auto& point3) {
          static constexpr hpuint o[3] = { 2, 0, 1 };

          auto j = make_neighbor_offset(neighbors, u, t);
          auto point1 = mesh.getVertex(u, o[j]).position - center;

          return glm::inverse(hpmat3x3(point0, point1, point2)) * point3;
     };

     transitions.reserve(9 * size(mesh));
     visit_triplets(mesh.getIndices(), [&](auto i0, auto i1, auto i2) {
          auto point0 = mesh.getVertex(i0).position - center;
          auto point1 = mesh.getVertex(i1).position - center;
          auto point2 = mesh.getVertex(i2).position - center;
          auto transition0 = make_transition(*(++n), point1, point0, point2);
          auto transition1 = make_transition(*(++n), point2, point1, point0);
          auto transition2 = make_transition(*(++n), point0, point2, point1);

          transitions.insert(std::end(transitions), {
               transition0.x, transition0.y, transition0.z,
               transition1.x, transition1.y, transition1.z,
               transition2.x, transition2.y, transition2.z
          });

          ++t;
     });
     
     return make_projective_structure(neighbors, std::move(transitions));
}

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph) { return make_projective_structure(graph, undegenerate(graph, trim(graph, cut(graph)))); }

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph, const Indices& cut) {
     auto& edges = graph.getEdges();
     auto analysis = analyze(graph, cut);
     auto& valences = std::get<0>(analysis);
     auto& indices = std::get<1>(analysis);
     auto& pairings = std::get<2>(analysis);
     auto lengths = std::vector<hpuint>();
     auto transitions = std::vector<hpreal>();
     auto w = hpreal(0);
     auto sun = std::vector<Point2D>();

     std::tie(sun, w) = detail::make_sun(valences);
     lengths.reserve(valences.size());
     for(auto i = std::begin(indices), end = std::end(indices) - 1; i != end; ++i) lengths.push_back(*(i + 1) - *i - 1);
     lengths.push_back(cut.size() - indices.back() + indices.front() - 1);

     auto polygon = std::vector<Point2D>(cut.size(), Point2D(0));//parametrize(lengths, polyline);TODO
     auto interior = parametrize(graph, cut, polygon);
     auto p = Indices(graph.getNumberOfVertices(), hpuint(0));
     auto n = hpuint(0);

     transitions.reserve(9 * size(graph));
     assert(cut.size() == polygon.size());
     for(auto& e : cut) p[edges[e].vertex] = std::numeric_limits<hpuint>::max();
     for(auto& v : p) if(v == hpuint(0)) v = n++;

     for(auto& edge : edges) {
          //continue if cut on edge (if both ends are not in interior)
          //make transition
     }

     for(auto e : cut) {
          //transform opposite vertex using pairing transition
          //make transition
     }

     return make_projective_structure(make_neighbors(graph), std::move(transitions));
}


//template<class Vertex>
//ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph) { return make_projective_structure(graph, undegenerate(graph, trim(graph, cut(graph)))); }
//ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph) { return make_projective_structure(graph, trim(graph, cut(graph))); }

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph, const Indices& cut, std::vector<Point2D> &test, std::vector<Point2D> &triangle, hpuint& hhh) {
     using happah::format::hph::operator<<;

     auto& edges = graph.getEdges();
     auto analysis = analyze(graph, cut);
     auto& valences = std::get<0>(analysis);
     auto& indices = std::get<1>(analysis);
     auto& pairings = std::get<2>(analysis);
     auto lengths = std::vector<hpuint>();
     auto transitions = std::vector<hpreal>();
     auto w = hpreal(0);
     auto sun = std::vector<Point2D>();

     std::tie(sun, w) = detail::make_sun(valences);
     transitions.reserve(9 * size(graph));
     lengths.reserve(valences.size());
     for(auto i = std::begin(indices), end = std::end(indices) - 1; i != end; ++i) lengths.push_back(*(i + 1) - *i - 1);
     lengths.push_back(cut.size() - indices.back() + indices.front() - 1);
     
     auto polyline = std::vector<Point2D>();
     auto j = std::begin(sun);
     auto* point0 = &sun[0];
     auto& sun0 = sun[valences.size() - 1];
     auto& sun2 = sun[0];
     auto sun1 = sun2 - (hpreal(indices.front()) / hpreal(lengths.back() + 1)) * (sun2 - sun0);
     
     auto do_parametrize = [&](auto& point0, auto& point1, auto m) {
          auto delta = (hpreal(1.0) / hpreal(m + 1)) * (point1 - point0);
          auto point = point0;

          polyline.emplace_back(point.x, point.y);
          while(m--) {
               point += delta;
               polyline.emplace_back(point.x, point.y);
          }
     };

     if(indices.front() > 0) do_parametrize(sun1, sun2, indices.front() - 1);
     for(auto i : boost::make_iterator_range(std::begin(lengths), std::end(lengths) - 1)) {
          auto& point1 = *(++j);

          do_parametrize(*point0, point1, i);
          point0 = &point1;
     }
     if(indices.front() > 0) do_parametrize(sun0, sun1, lengths.back() - indices.front());
     else do_parametrize(sun0, sun2, lengths.back());
     
     
     auto interior = parametrize(graph, cut, polyline);
     auto p = Indices(graph.getNumberOfVertices(), hpuint(0));
     auto n = hpuint(0);

     
     for(auto& e : cut) p[edges[e].vertex] = std::numeric_limits<hpuint>::max();
     for(auto& v : p) if(v == hpuint(0)) v = n++;
     
     auto get_point = [&](auto& e) {
          auto w = make_spokes_walker(edges, e); 
          auto temp = std::end(cut);
          while(temp == std::end(cut)) temp = std::find(std::begin(cut), std::end(cut), edges[*(++w)].opposite);
          return polyline[std::distance(std::begin(cut), temp)];
     };

     assert(polyline.size() == cut.size());
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
          auto point0 = Point3D((c0) ? get_point(i0) : interior[p[p0]], 1);
          auto point1 = Point3D((c1) ? get_point(edge1.opposite) : interior[p[p1]], 1);
          auto point2 = Point3D((c2) ? get_point(edge.previous) : interior[p[p2]], 1);
          auto point3 = Point3D();
          
          if(!c0 && !c1 && !c2 && !found_first_inner_triangle){
               found_first_inner_triangle = true;
               hhh = edge.next / 3;
               triangle.push_back(interior[p[p1]]);
               triangle.push_back(interior[p[p0]]);
               triangle.push_back(interior[p[p2]]);
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
               point3 = Point3D((c3) ? get_point(edge3.next) : interior[p[p3]], 1);
          }

          transition = glm::inverse(hpmat3x3(point0, point3, point1)) * point2;
          std::cout << "p0: " << point0 << ' ' << c0 << '\n';
          std::cout << "p1: " << point1 << ' ' << c1 << '\n';
          std::cout << "p2: " << point2 << ' ' << c2 << '\n';
          std::cout << "p3: " << point3 << ' ' << c3 << '\n';
          std::cout << "transition " << transition.x << ", " << transition.y << ", " << transition.z << std::endl;
          if(std::find(std::begin(cut), std::end(cut), make_edge_index(edge)) == std::end(cut)) assert(glm::length(transition) > EPSILON);
          std::cout << "-------------" << std::endl;
          
          transitions.push_back(transition.x);
          transitions.push_back(transition.y);
          transitions.push_back(transition.z);
     }

     std::cout << "number transitions " << transitions.size() << '\n';
     std::cout << "indices.front " << indices.front() << '\n';
     test.insert(std::end(test), std::begin(polyline), std::end(polyline));
     test.insert(std::end(test), std::begin(interior), std::end(interior));
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

