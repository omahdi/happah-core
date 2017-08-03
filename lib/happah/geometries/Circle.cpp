// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/Circle.hpp"

namespace happah {

boost::optional<boost::variant<Point2D, std::tuple<Point2D, Point2D> > > intersect(const Circle& circle0, const Circle& circle1, hpreal epsilon) {
     auto& center0 = circle0.getCenter();
     auto& center1 = circle1.getCenter();
     auto normal = center1 - center0;
     auto r0 = circle0.getRadius();
     auto r1 = circle1.getRadius();
     auto d = glm::length(normal);
     if(d < epsilon && std::abs(r0 - r1) < epsilon) return boost::none;//circles are the same
     if(d > r0 + r1) return boost::none;//circle0 and circle1 are disjoint
     if(d < std::abs(r0 - r1)) return boost::none; //either circle0 or circle1 is contained in the other circle
     normal /= d;
     if(std::abs(d - (r0 + r1)) < epsilon) return boost::variant<Point2D, std::tuple<Point2D, Point2D> >(center0 + r0 * normal);
     auto a = (r0 * r0 - r1 * r1 + d * d) / (d + d);
     auto h = std::sqrt(r0 * r0 - a * a);
     auto tangent = h * Vector2D(-normal.y, normal.x);
     auto temp = center0 + a * normal;
     return boost::variant<Point2D, std::tuple<Point2D, Point2D> >(std::make_tuple(temp + tangent, temp - tangent));
}

Circle poincare_to_euclidean(const Circle& circle) { 
     auto r = circle.getRadius();
     auto C = circle.getCenter();
     auto u = (std::exp(r)-1)/(std::exp(r)+1);
     auto dhyp = glm::length2(C);
     auto temp1 = (2-2*u*u)/(1-u*u*dhyp);
     Point2D eucl_center;
     eucl_center.x = temp1*C.x;
     eucl_center.y = temp1*C.y;
     auto deucl = glm::length2(eucl_center);
     auto temp2 = deucl - (dhyp-u*u)/(1-u*u*dhyp);
     auto eucl_r = sqrt(temp2);
     return make_circle(eucl_center, eucl_r);
}

}//namespace happah

