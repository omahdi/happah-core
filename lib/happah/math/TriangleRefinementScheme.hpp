// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/math/Space.hpp"

namespace happah {

struct TriangleRefinementScheme {
     static const TriangleRefinementScheme BINARY_UNIFORM;
     static const TriangleRefinementScheme POWELL_SABIN_12_SPLIT;
     static const TriangleRefinementScheme POWELL_SABIN_6_SPLIT;
     static const TriangleRefinementScheme QUATERNARY_UNIFORM;
     static const TriangleRefinementScheme TERNARY_UNIFORM;

     std::vector<hpuint> indices;//every three indices form a triangle in counterclockwise order
     std::vector<Point3D> points;//TODO: Point2D (w = 1 - u - v)

     TriangleRefinementScheme(std::vector<hpuint> indices, std::vector<Point3D> points)
          : indices(std::move(indices)), points(std::move(points)) {}

     hpuint getNumberOfTriangles() const { return indices.size() / 3; }

};//TriangleRefinementScheme

}//namespace happah

