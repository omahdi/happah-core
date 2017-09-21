// Copyright 2017
//   Obada Mahdi    -                                       - omahdi@gmail.com
//   Pawel Herman   - Karlsruhe Institute of Technology     - pherman@ira.uka.de
//   Louis Lutzweiler                                       - louis.lutzweiler@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.10 - Louis Lutzweiler - Added make_projective_structure.

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

     template<class Stream>
     friend Stream& operator<<(Stream& s, const ProjectiveStructure& ps) {
          using format::hph::operator<<;
          s << ps.m_neighbors;
          s << ps.m_transitions;
          return s;
     }
     template<class Stream>
     friend Stream& operator>>(Stream& s, ProjectiveStructure& ps) {
          using format::hph::operator>>;
          s >> ps.m_neighbors;
          s >> ps.m_transitions;
          return s;
     }

private:
     Indices m_neighbors;
     std::vector<hpreal> m_transitions;

};//ProjectiveStructure

inline ProjectiveStructure make_projective_structure(Indices neighbors, std::vector<hpreal> transitions) { return { std::move(neighbors), std::move(transitions) }; }

template<class Vertex>
ProjectiveStructure make_projective_structure(const TriangleMesh<Vertex>& mesh, const Point3D& center, const Indices& neighbors) {
     auto triangles = deindex(mesh.getVertices(), mesh.getIndices());
     auto num = mesh.getNumberOfTriangles();
     auto transitions = std::vector<hpreal>(9 * num);
     static constexpr hpuint o[3] = { 2, 0, 1 };
     auto t = hpindex(0);

     visit_triplets(mesh.getIndices(), [&](auto i0, auto i1, auto i2) {
          auto& vertex0 = mesh.getVertex(i0);
          auto& vertex1 = mesh.getVertex(i1);
          auto& vertex2 = mesh.getVertex(i2);
          
          auto mat0 = hpmat3x3(vertex0.position - center, vertex2.position - center, vertex1.position - center);
          auto mat1 = hpmat3x3(vertex1.position - center, vertex0.position - center, vertex2.position - center);
          auto mat2 = hpmat3x3(vertex2.position - center, vertex1.position - center, vertex0.position - center);
          mat0 = glm::inverse(mat0);
          mat1 = glm::inverse(mat1);
          mat2 = glm::inverse(mat2);
 
          auto make_transition = [&](int e, hpmat3x3 matrix) {
               auto u = make_neighbor_index(neighbors, t, e);
               auto j = make_neighbor_offset(neighbors, u, t);
               auto vertex3 = triangles[3 * u + o[j]];

               auto transition = matrix * (vertex3.position - center);

               transitions.at(3*(3*u+j)) = transition.x;
               transitions.at(3*(3*u+j) + 1) = transition.y;
               transitions.at(3*(3*u+j) + 2) = transition.z;
          };
     
          make_transition(0, mat0);
          make_transition(1, mat1);
          make_transition(2, mat2);

          ++t;
     });
     
     return ProjectiveStructure(neighbors, transitions);
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

