// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>

#include "happah/geometries/TriangleMesh.h"
#include "happah/geometries/SurfaceSplineBEZ.h"

namespace happah {

template<class Vertex>
QuarticSurfaceSplineBEZ<typename Vertex::SPACE> make_quartic_surface_spline_bez(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh) {
     using Indices = std::vector<hpuint>;
     using Iterator = Indices::const_iterator;

     boost::dynamic_bitset<> visited(mesh.getNumberOfVertices(), false);

     Indices indices;
     std::vector<typename Vertex::SPACE::POINT> controlPoints;

     //TODO

     return { std::move(controlPoints), std::move(indices) };
}

}//namespace happah

