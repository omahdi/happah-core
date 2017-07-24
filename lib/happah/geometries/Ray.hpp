// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <utility>

#include "happah/geometries/Geometry.hpp"

namespace happah {

//TODO: show ray cloud and line mesh as cylinders using shader programs
template<class Space>
class Ray : public Geometry1D<Space> {
     using Point = typename Space::POINT;
     using Vector = typename Space::VECTOR;

public:
     static bool isInHalfspace(const Point& origin, const Vector& direction, const Point& point) {
          //TODO
          return false;
     }

     //Ray(Point p0, const Point& p1)
     //     : m_origin{std::move(p0)}, m_direction{p1-p0} {}

     Ray(Point origin, Vector direction)
          : m_origin(std::move(origin)), m_direction(std::move(direction)) {}

     const Vector& getDirection() const { return m_direction; }

     const Point& getOrigin() const { return m_origin; }

     Point getPoint(hpreal t) const { return m_origin + t * m_direction; }

     bool isInHalfspace(const Point& point) const { return isInHalfspace(m_origin, m_direction, point); }

     void normalize() { m_direction = glm::normalize(m_direction); }

     void setDirection(Vector direction) { m_direction = std::move(direction); }

     void setOrigin(Point origin) { m_origin = std::move(origin); }

private:
     Point m_origin;
     Vector m_direction;
     //NOTE: The order of the member variables was changed to match the order of the member variables in VertexPN allowing reinterpretation casts between Ray and VertexPN.

};
using Ray2D = Ray<Space2D>;
using Ray3D = Ray<Space3D>;

Ray3D make_ray(Vector3D direction);

}//namespace happah

