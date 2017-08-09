// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/math/Space.hpp"

namespace happah {

//DECLARATIONS

template<class Space>
class Ray;

using Ray2D = Ray<Space2D>;
using Ray3D = Ray<Space3D>;

template<class Space>
auto make_point(const Ray<Space>& ray, hpreal t);

inline Ray3D make_ray(Vector3D direction);

template<class Space>
Ray<Space> normalize(Ray<Space> ray);

//DEFINITIONS

template<class Space>
class Ray {
     using Point = typename Space::POINT;
     using Vector = typename Space::VECTOR;

public:
     Ray(Point origin, Vector direction)
          : m_origin(std::move(origin)), m_direction(std::move(direction)) {}

     auto& getDirection() const { return m_direction; }

     auto& getOrigin() const { return m_origin; }

     void setDirection(Vector direction) { m_direction = std::move(direction); }

     void setOrigin(Point origin) { m_origin = std::move(origin); }

private:
     Point m_origin;
     Vector m_direction;
     //NOTE: The order of the member variables was changed to match the order of the member variables in VertexPN allowing reinterpretation casts between Ray and VertexPN.

};//Ray

template<class Space>
auto make_point(const Ray<Space>& ray, hpreal t) { return ray.getOrigin() + t * ray.getDirection(); }

inline Ray3D make_ray(Vector3D direction) { return { { 0, 0, 0 }, std::move(direction) }; }

template<class Space>
Ray<Space> normalize(Ray<Space> ray) {
     ray.setDirection(glm::normalize(ray.getDirection()));
     return ray;
}

}//namespace happah

