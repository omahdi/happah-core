// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/optional.hpp>
#include <type_traits>
#include <vector>

#include "happah/Happah.hpp"
#include "happah/geometries/Ray.hpp"

namespace happah {

//DECLARATIONS

class Plane;

//NOTE: Ray's direction must be normalized.
boost::optional<hpreal> intersect(const Plane& plane, const Ray3D& ray, hpreal epsilon = EPSILON);

std::tuple<Vector3D, Vector3D> make_axes(const Plane& plane, hpreal epsilon = EPSILON);

template<class T>
std::vector<hpreal> make_factors(const std::vector<T>& ts, const Plane& plane, hpreal epsilon = EPSILON);

//NOTE: Returns the transformation that transforms the x-y plane to this plane.
std::tuple<hpmat3x3, hpvec3> make_transformation(const Plane& plane, hpreal epsilon = happah::EPSILON);

//NOTE: Point must be in local coordinate system of the plane.
Point2D project(const Plane& plane, const Point3D& point);

template<class T>
std::vector<Point2D> restrict(const std::vector<T>& ts, const Plane& plane, const std::vector<hpreal>& factors);

Ray3D unproject(const Plane& plane, const Point2D& point);

//DEFINITIONS

class Plane {
public:
     Plane(Point3D origin, Vector3D normal)
          : m_normal(std::move(normal)), m_origin(std::move(origin)) {}

     auto& getNormal() const { return m_normal; }

     auto& getOrigin() const { return m_origin; }

     void setNormal(Vector3D normal) { m_normal = std::move(normal); }

     void setOrigin(Point3D origin) { m_origin = std::move(origin); }

private:
     Vector3D m_normal;
     Point3D m_origin;

};//Plane

template<class T>
std::vector<hpreal> make_factors(const std::vector<T>& ts, const Plane& plane, hpreal epsilon) {
     auto factors = std::vector<hpreal>();

     factors.reserve(ts.size());

     for(auto& t : ts) {
          auto ray = make_ray(t);
          ray.normalize();
          if(auto lambda = intersect(plane, ray, epsilon)) {
               auto intersection = ray.getPoint(*lambda);
               factors.push_back(std::sqrt(length2(intersection) / length2(t)));
          } else throw std::runtime_error("Plane does not properly intersect all tetrahedra.");
     }

     return factors;
}

template<class T>
std::vector<Point2D> restrict(const std::vector<T>& ts, const Plane& plane, const std::vector<hpreal>& factors) {
     auto points = std::vector<Point2D>();
     auto f = factors.begin();

     points.reserve(ts.size());

     for(auto& t : ts) {
          points.push_back(project(plane, mix(t, *f)));
          ++f;
     }

     return points;
}

}//namespace happah

