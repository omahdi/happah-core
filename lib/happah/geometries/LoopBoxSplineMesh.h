// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/LoopBoxSpline.h"
#include "happah/geometries/TriangleMesh.h"

namespace happah {

template<class Vertex>
class LoopBoxSplineMesh {
public:
     LoopBoxSplineMesh(std::vector<Vertex> controlPoints, Indices indices)
          : m_controlPoints(std::move(controlPoints)), m_indices(std::move(indices)) {}

     auto& getControlPoints() const { return m_controlPoints; }

     auto& getIndices() const { return m_indices; }

private:
     std::vector<Vertex> m_controlPoints;
     Indices m_indices;

};//LoopBoxSplineMesh

template<class Vertex>
LoopBoxSplineMesh<Vertex> make_loop_box_spline_mesh(const TriangleMesh<Vertex, Format::SIMPLE>& mesh) {
     auto neighbors = make_neighbors(mesh);
     auto indices = Indices();
     auto valences = make_valences(mesh);
     auto t = 0u;

     visit_triplets(mesh.getIndices(), [&](auto i0, auto i1, auto i2) {
          if(valences[i0] == 6 && valences[i1] == 6 && valences[i2] == 6) {//TODO: and not on border
               auto ring0 = make_ring(make_ring_enumerator(neighbors, t, 0), mesh.getIndices());
               auto ring1 = make_ring(make_ring_enumerator(neighbors, t, 1), mesh.getIndices());
               auto ring2 = make_ring(make_ring_enumerator(neighbors, t, 2), mesh.getIndices());
               assert(ring0.size() == 6 && ring1.size() == 6 && ring2.size() == 6);
               auto temp = make_loop_box_spline_control_points(i0, std::begin(ring0), i1, std::begin(ring1), i2, std::begin(ring2));
               indices.insert(std::end(indices), std::begin(temp), std::end(temp));
          }
          ++t;
     });

     return { mesh.getVertices(), std::move(indices) };
}

template<class Vertex>
hpuint size(const LoopBoxSplineMesh<Vertex>& mesh) { return mesh.getIndices().size() / 12; }

}//namespace happah

