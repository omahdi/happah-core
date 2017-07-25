// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/Circle.hpp"

namespace happah {

Circle::Circle(Point2D center, hpreal radius)
     : m_center(std::move(center)), m_radius(radius) {}

const Point2D& Circle::getCenter() const { return m_center; }

hpreal Circle::getRadius() const { return m_radius; }

Circle make_circle(Point2D center, hpreal radius) { return { std::move(center), radius }; }

}//namespace happah

