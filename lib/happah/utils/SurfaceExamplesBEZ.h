// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <memory>
#include <vector>

#include "happah/math/Space.h"

class SurfaceExamplesBEZ {
public:
     using ControlPoints = std::vector<Point3D>;

     static const ControlPoints CUBIC_CONTROL_POINTS;
     static const ControlPoints QUADRATIC_CONTROL_POINTS;
     static const ControlPoints QUARTIC_CONTROL_POINTS;

};

