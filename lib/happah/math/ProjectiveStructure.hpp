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

std::tuple<std::vector<Point2D>, hpreal> make_sun(const Indices& valences);

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
ProjectiveStructure make_projective_structure(const TriangleGraph<Vertex>& graph) { return make_projective_structure(graph, trim(graph, cut(graph))); }

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

     std::tie(sun, w) = make_sun(valences);
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

