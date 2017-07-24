// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/TriangleRefinementScheme.hpp"

namespace happah {

const TriangleRefinementScheme TriangleRefinementScheme::BINARY_UNIFORM = TriangleRefinementScheme({
     0, 1, 3, 
     1, 2, 4, 
     1, 4, 3, 
     4, 5, 3
}, std::vector<Point3D>({
     Point3D(1, 0, 0),
     Point3D(0.5, 0.5, 0),
     Point3D(0, 1, 0),
     Point3D(0.5, 0, 0.5),
     Point3D(0, 0.5, 0.5),
     Point3D(0, 0, 1)
}));

const TriangleRefinementScheme TriangleRefinementScheme::POWELL_SABIN_12_SPLIT = TriangleRefinementScheme({
     0, 1, 3,
     3, 1, 5,
     5, 1, 4,
     4, 1, 2,
     6, 0, 3,
     6, 3, 5,
     6, 5, 7,
     7, 5, 8,
     8, 5, 4,
     8, 4, 2,
     6, 7, 9,
     7, 8, 9
}, std::vector<Point3D>({
     Point3D(1, 0, 0),
     Point3D(0.5, 0.5, 0),
     Point3D(0, 1, 0),
     Point3D(0.5, 0.25, 0.25),
     Point3D(0.25, 0.5, 0.25),
     Point3D(1.0/3.0, 1.0/3.0, 1.0/3.0),
     Point3D(0.5, 0, 0.5),
     Point3D(0.25, 0.25, 0.5),
     Point3D(0, 0.5, 0.5),
     Point3D(0, 0, 1)
}));

const TriangleRefinementScheme TriangleRefinementScheme::POWELL_SABIN_6_SPLIT = TriangleRefinementScheme({
     0, 4, 3,
     0, 1, 4,
     1, 2, 4,
     4, 2, 5,
     3, 4, 6,
     6, 4, 5
}, std::vector<Point3D>({
     Point3D(1, 0, 0),
     Point3D(0.5, 0.5, 0),
     Point3D(0, 1, 0),
     Point3D(0.5, 0, 0.5),
     Point3D(1.0/3.0, 1.0/3.0, 1.0/3.0),
     Point3D(0, 0.5, 0.5),
     Point3D(0, 0, 1)
}));

const TriangleRefinementScheme TriangleRefinementScheme::QUATERNARY_UNIFORM = TriangleRefinementScheme({
     0, 1, 5,
     1, 2, 6,
     2, 3, 7,
     3, 4, 8,
     1, 6, 5,
     2, 7, 6,
     3, 8, 7,
     5, 6, 9,
     6, 7, 10,
     7, 8, 11,
     6, 10, 9,
     7, 11, 10,
     9, 10, 12,
     10, 11, 13,
     10, 13, 12,
     12, 13, 14
}, std::vector<Point3D>({
     Point3D(1, 0, 0),
     Point3D(0.75, 0.25, 0),
     Point3D(0.5, 0.5, 0),
     Point3D(0.25, 0.75, 0),
     Point3D(0, 1, 0),
     Point3D(0.75, 0, 0.25),
     Point3D(0.5, 0.25, 0.25),
     Point3D(0.25, 0.5, 0.25),
     Point3D(0, 0.75, 0.25),
     Point3D(0.5, 0, 0.5),
     Point3D(0.25, 0.25, 0.5),
     Point3D(0, 0.5, 0.5),
     Point3D(0.25, 0, 0.75),
     Point3D(0, 0.25, 0.75),
     Point3D(0, 0, 1)
}));

const TriangleRefinementScheme TriangleRefinementScheme::TERNARY_UNIFORM = TriangleRefinementScheme({
     0, 1, 4,
     1, 2, 5,
     2, 3, 6,
     1, 5, 4,
     2, 6, 5,
     4, 5, 7,
     5, 6, 8,
     5, 8, 7,
     7, 8, 9     
}, std::vector<Point3D>({
     Point3D(1, 0, 0),
     Point3D(2.0/3.0, 1.0/3.0, 0),
     Point3D(1.0/3.0, 2.0/3.0, 0),
     Point3D(0, 1, 0),
     Point3D(2.0/3.0, 0, 1.0/3.0),
     Point3D(1.0/3.0, 1.0/3.0, 1.0/3.0),
     Point3D(0, 2.0/3.0, 1.0/3.0),
     Point3D(1.0/3.0, 0, 2.0/3.0),
     Point3D(0, 1.0/3.0, 2.0/3.0),
     Point3D(0, 0, 1)
}));

TriangleRefinementScheme::TriangleRefinementScheme(std::vector<hpuint> indices, std::vector<Point3D> points)
     : indices(std::move(indices)), points(std::move(points)) {}

hpuint TriangleRefinementScheme::getNumberOfTriangles() const { return indices.size() / 3; }

}//namespace happah

