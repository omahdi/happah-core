// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/Circle.hpp"

namespace happah {

boost::optional<boost::variant<Point2D, std::tuple<Point2D, Point2D> > > intersect(const Circle& circle0, const Circle& circle1, hpreal epsilon) {
     using T = boost::variant<Point2D, std::tuple<Point2D, Point2D> >;

     auto& center0 = circle0.getCenter();
     auto& center1 = circle1.getCenter();
     auto normal = center1 - center0;
     auto r0 = circle0.getRadius();
     auto r1 = circle1.getRadius();
     auto d = glm::length(normal);
     if(d < epsilon && std::abs(r0 - r1) < epsilon) return boost::none;//circles are the same
     if(d > r0 + r1) return boost::none;//circle0 and circle1 are disjoint
     if(d < std::abs(r0 - r1)) return boost::none;//either circle0 or circle1 is contained in the other circle
     normal /= d;
     if(std::abs(d - (r0 + r1)) < epsilon) return T(center0 + r0 * normal);
     auto a = (r0 * r0 - r1 * r1 + d * d) / (d + d);
     auto h = std::sqrt(r0 * r0 - a * a);
     auto tangent = h * Vector2D(-normal.y, normal.x);
     auto temp = center0 + a * normal;
     return T(std::make_tuple(temp + tangent, temp - tangent));
}

Circle poincare_to_euclidean(const Circle& circle) {
     auto& center = circle.getCenter();
     auto c2 = glm::length2(center);
     auto t = std::tanh(circle.getRadius() / hpreal(2.0));
     auto t2 = t * t;
     auto d = hpreal(1.0) - t2 * c2;
     return make_circle(((hpreal(1.0) - t2) / d) * center, t * (hpreal(1.0) - c2) / d);
}

}//namespace happah

