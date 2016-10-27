// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/Space.h"

hpreal length2(const Point3D& point) { return glm::length2(point); }

Point3D mix(const Point3D& point, hpreal lambda) { return point * lambda; }

Point2D mix(const Point2D& p0, hpreal u, const Point2D& p1, hpreal v, const Point2D& p2, hpreal w) { return p0 * u + p1 * v + p2 * w; }

