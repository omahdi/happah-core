// Copyright 2015-2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "LoopBoxSpline.h"

LoopBoxSpline::LoopBoxSpline(std::array<glm::vec3, 12> controlPoints)
     : m_controlPoints(std::move(controlPoints)) {}

//LoopBoxSpline::LoopBoxSpline(std::vector<glm::vec3> controlPoints) { std::copy_n(std::make_move_iterator(controlPoints.begin()), 12, m_controlPoints.begin()); }

const std::array<glm::vec3, 12>& LoopBoxSpline::getControlPoints() const { return m_controlPoints; }

//    10  11
//   7   8   9
// 3   4   5   6
//   0   1   2
TriangleMesh LoopBoxSpline::getControlPolygon() const { return {{m_controlPoints.begin(), m_controlPoints.end()},{0,1,4,1,2,5,3,0,4,4,1,5,5,2,6,3,4,7,4,5,8,5,6,9,7,4,8,8,5,9,7,8,10,8,9,11,10,8,11}}; }

