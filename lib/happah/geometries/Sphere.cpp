// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/Sphere.hpp"

namespace happah {

boost::optional<boost::variant<Point3D, std::tuple<Point3D, Point3D> > > intersect(const Sphere& sphere, const Ray3D& ray, hpreal epsilon) {
     using T = boost::variant<Point3D, std::tuple<Point3D, Point3D> >;

     auto temp = sphere.getCenter() - ray.getOrigin();
     auto foot = glm::dot(temp, ray.getDirection());
     auto r2 = sphere.getRadius() * sphere.getRadius();
     auto v2 = glm::dot(temp, temp);
     if(foot < 0 && v2 > r2) return boost::none;
     auto d2 = v2 - foot * foot;
     if(d2 > r2) return boost::none;
     auto delta = glm::sqrt(r2 - d2);
     if(v2 > r2) return T(std::make_tuple(make_point(ray, foot - delta), make_point(ray, foot + delta)));
     else return T(make_point(ray, foot + delta));
}

}//namespace happah

