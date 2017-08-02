// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <stack>
#include <vector>

#include "happah/Happah.hpp"
#include "happah/geometries/TriangleMesh.hpp"
#include "happah/geometries/Vertex.hpp"
#include "happah/utils/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

class ProjectiveStructure;

ProjectiveStructure make_projective_structure(Indices neighbors, std::vector<hpreal> transitions);

//NOTE: Border has to be sorted.
template<class Vertex = VertexP3, class VertexFactory = happah::VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const ProjectiveStructure& structure, const Indices& border, hpreal t, const Point3D& p0, const Point3D& p1, const Point3D& p2, VertexFactory&& factory = VertexFactory());

bool validate(const ProjectiveStructure& structure, hpreal epsilon = EPSILON);

//DEFINITIONS

class ProjectiveStructure {
public:
     ProjectiveStructure(Indices neighbors, std::vector<hpreal> transitions);

     const Indices& getNeighbors() const;

     const std::vector<hpreal>& getTransitions() const;

private:
     Indices m_neighbors;
     std::vector<hpreal> m_transitions;

};//ProjectiveStructure

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const ProjectiveStructure& structure, const Indices& border, hpreal t, const Point3D& p0, const Point3D& p1, const Point3D& p2, VertexFactory&& factory) {
     auto& neighbors = structure.getNeighbors();
     auto& transitions = structure.getTransitions();
     auto vertices = std::vector<Vertex>();
     auto indices = Indices(neighbors.size(), std::numeric_limits<hpindex>::max());
     auto todo = std::stack<hpindex>();

     auto push = [&](auto position, auto t, auto i) {
          auto n = vertices.size();
          vertices.push_back(factory(position));
          visit(make_spokes_walker(neighbors, t, i), border, [&](auto t, auto i) { indices[3 * t + i] = n; });
     };

     assert(*std::max_element(std::begin(neighbors), std::end(neighbors)) < std::numeric_limits<hpuint>::max());//NOTE: Implementation assumes a closed topology.

     push(p0, t, 0);
     push(p1, t, 1);
     push(p2, t, 2);
     todo.emplace(3 * t + 0);
     todo.emplace(3 * t + 1);
     todo.emplace(3 * t + 2);

     while(!todo.empty()) {
          static constexpr hpuint o0[3] = { 0, 1, 2 };
          static constexpr hpuint o1[3] = { 2, 0, 1 };
          static constexpr hpuint o2[3] = { 1, 2, 0 };

          auto e = todo.top();
          todo.pop();
          if(std::binary_search(std::begin(border), std::end(border), e)) continue;
          auto u = make_triangle_index(e);
          auto j = make_edge_offset(e);
          auto v = make_neighbor_index(neighbors, u, j);
          auto k = make_neighbor_offset(neighbors, v, u);
          if(indices[3 * v + o1[k]] != std::numeric_limits<hpindex>::max()) continue;
          auto transition = std::begin(transitions) + 3 * (3 * v + k);
          auto temp = std::begin(indices) + 3 * u;
          auto& v0 = vertices[temp[o0[j]]];
          auto& v1 = vertices[temp[o1[j]]];
          auto& v2 = vertices[temp[o2[j]]];
          push(transition[0] * v0.position + transition[1] * v1.position + transition[2] * v2.position, v, o1[k]);
          todo.emplace(3 * v + o1[k]);
          todo.emplace(3 * v + o2[k]);
     }

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

}//namespace happah

