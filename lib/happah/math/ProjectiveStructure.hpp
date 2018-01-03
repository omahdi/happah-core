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
#include "happah/util/visitors.hpp"

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

inline ProjectiveStructure make_projective_structure(Triples<hpindex> neighbors, std::vector<hpreal> transitions);

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleMesh<Vertex>& mesh, const Point3D& center, const Triples<hpindex>& neighbors);

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
     ProjectiveStructure(Triples<hpindex> neighbors, std::vector<hpreal> transitions)
          : m_neighbors(std::move(neighbors)), m_transitions(std::move(transitions)) {}

     auto& getNeighbors() const { return m_neighbors; }

     auto& getTransitions() const { return m_transitions; }

private:
     Triples<hpindex> m_neighbors;
     std::vector<hpreal> m_transitions;

};//ProjectiveStructure

inline ProjectiveStructure make_projective_structure(Triples<hpindex> neighbors, std::vector<hpreal> transitions) { return { std::move(neighbors), std::move(transitions) }; }

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleMesh<Vertex>& mesh, const Point3D& center, const Triples<hpindex>& neighbors) {
     auto transitions = std::vector<hpreal>();
     auto t = hpindex(0);
     auto n = std::begin(neighbors) - 1;

     auto make_transition = [&](auto u, auto& point0, auto& point2, auto& point3) {
          static constexpr hpuint o[3] = { 2, 0, 1 };

          auto j = make_offset(neighbors, u, t);
          auto point1 = mesh.getVertex(u, trit(o[j])).position - center;

          return glm::inverse(hpmat3x3(point0, point1, point2)) * point3;
     };

     transitions.reserve(9 * size(mesh));
     visit(mesh.getIndices(), [&](auto i0, auto i1, auto i2) {
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
     auto center = Point3D(0, 0, 1);

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
     for(auto& edge : edges) {
          auto& edge1 = edges[edge.opposite];
          auto& edge2 = edges[edge.next];
          auto& edge3 = edges[edge1.next];
          auto& i0 = p[edge.vertex];
          auto& i1 = p[edge1.vertex];
          auto& i2 = p[edge2.vertex];
          auto& i3 = p[edge3.vertex];
          auto c0 = i0 == std::numeric_limits<hpuint>::max();
          auto c1 = i1 == std::numeric_limits<hpuint>::max();
          auto c2 = i2 == std::numeric_limits<hpuint>::max();
          auto c3 = i3 == std::numeric_limits<hpuint>::max();
          auto point0 = Point3D((c0) ? get_point(edge.next) : interior[i0], 1);
          auto point1 = Point3D((c1) ? get_point(edge1.opposite) : interior[i1], 1);
          auto point2 = Point3D((c2) ? get_point(edge.previous) : interior[i2], 1);
          auto point3 = Point3D();
          auto transition = Point3D();

          if(std::find(std::begin(cut), std::end(cut), edge1.opposite) == std::end(cut)) {
               point3 = Point3D((c3) ? get_point(edge3.next) : interior[i3], 1);
               transition = glm::inverse(hpmat3x3(point0, point3, point1)) * point2;
          }
          transitions.push_back(transition.x);
          transitions.push_back(transition.y);
          transitions.push_back(transition.z);
     }
     
     auto A = hpmat3x3();
     auto o = std::begin(sun) + valences.size();
     auto r = std::begin(pairings);
     auto s = std::begin(sun);

     auto do_transition = [&](auto e) {
          auto& edge = edges[e];
          auto& edge1 = edges[edge.opposite];
          auto& edge2 = edges[edge.next];
          auto& edge3 = edges[edge1.next];
          auto& i2 = p[edge2.vertex];
          auto& i3 = p[edge3.vertex];
          auto c2 = i2 == std::numeric_limits<hpuint>::max();
          auto c3 = i3 == std::numeric_limits<hpuint>::max();
          auto point0 = Point3D(get_point(edge.next), 1);
          auto point1 = Point3D(get_point(e), 1);
          auto point2 = Point3D((c2) ? get_point(edge.previous) : interior[i2], 1);
          auto point3 = Point3D((c3) ? get_point(edge3.next) : interior[i3], 1);
          auto transition = glm::inverse(hpmat3x3(point0, A * point3, point1)) * point2;

          transitions[3 * e + 0] = transition.x;
          transitions[3 * e + 1] = transition.y;
          transitions[3 * e + 2] = transition.z;
     };

     auto make_A = [&](auto point0, auto point1, auto point2, auto point3, auto point4) { return hpmat3x3(Point3D(point0, 1), Point3D(point1, 1), Point3D(point2 * w, w)) * glm::inverse(hpmat3x3(Point3D(point4, 1), Point3D(point3, 1), center)); };

     visit_pairs(std::begin(indices), indices.size() - 1, 1, [&](auto i0, auto i1) {
          A = make_A(s[0], s[1], o[0], sun[*r], sun[(*r == valences.size() - 1) ? 0 : *r + 1]);
          for(auto e : boost::make_iterator_range(std::begin(cut) + (i0 + 1), std::begin(cut) + (i1 + 1))) do_transition(e);
          ++o;
          ++r;
          ++s;
     });
     A = make_A(sun[valences.size() - 1], sun[0], sun.back(), sun[pairings.back()], sun[pairings.back() + 1]);
     for(auto e : boost::make_iterator_range(std::begin(cut) + indices.back() + 1, std::end(cut))) do_transition(e);
     for(auto e : boost::make_iterator_range(std::begin(cut), std::begin(cut) + indices.front() + 1)) do_transition(e);

     return make_projective_structure(make_neighbors(graph), std::move(transitions));
}

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const ProjectiveStructure& structure, const Indices& border, hpindex t, const Point3D& point0, const Point3D& point1, const Point3D& point2, VertexFactory&& build) {
     auto& neighbors = structure.getNeighbors();
     auto& transitions = structure.getTransitions();

     return make_triangle_mesh(neighbors, border, t, build(point0), build(point1), build(point2), [&](auto t, auto i, auto& vertex0, auto& vertex1, auto& vertex2) {
          auto u = neighbors(t, i);
          auto j = make_offset(neighbors, u, t);
          auto transition = std::begin(transitions) + 3 * (3 * u + j);

          return build(transition[0] * vertex1.position + transition[1] * vertex2.position + transition[2] * vertex0.position);
     });
}

}//namespace happah

