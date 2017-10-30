// Copyright 2017
//   Obada Mahdi    -                                       - omahdi@gmail.com
//   Pawel Herman   - Karlsruhe Institute of Technology     - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/Happah.hpp"
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

