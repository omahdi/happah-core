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

boost::optional<std::tuple<Point2D, Point2D> > intersect(const Circle& circle0, const Circle& circle1) {
     auto x1 = circle0.getCenter().x;
     auto y1 = circle0.getCenter().y;
     auto x2 = circle1.getCenter().x;
     auto y2 = circle1.getCenter().y;
     auto r1 = circle0.getRadius();
     auto r2 = circle1.getRadius();
     auto d = glm::length(circle0.getCenter() - circle1.getCenter());
     if (d > r1+r2) return boost::none; // circles must intersect
     if (d < std::abs(r1-r2)) return boost::none; // circles mustn't be contained within the other. 
     if ((x1 == x2) && (y1 == y2)) return boost::none;
     Point2D result0;
     Point2D result1;
     if (x1 == x2) {
          auto Y = (r1*r1-r2*r2-y1*y1+y2*y2)/(-2*(y1-y2));
          auto a = 1;
          auto b = -2*x1;
          auto c = x1*x1 +Y*Y - 2*Y*y2 + y2*y2 - r2*r2;
          auto delta = b*b -4*a*c;
          assert(!(delta < 0));
          auto sol1_x = (-b - std::sqrt(delta))/(2*a) ;
          auto sol2_x = (-b + std::sqrt(delta))/(2*a) ;
          result0.x = sol1_x;
          result0.y = Y;
          result1.x = sol2_x;
          result1.y = Y;
     } else if (y1 == y2) {
          auto X = (r1*r1-r2*r2-x1*x1+x2*x2)/(2*(x2-x1));
          auto a = 1;
          auto b = -2*y1;
          auto c = X*X - 2*X*x2 + x2*x2 + y1*y1 - r2*r2;
          auto delta = b*b -4*a*c;
          assert(!(delta < 0));
          auto sol1_y = (-b -std::sqrt(delta))/(2*a) ;
          auto sol2_y = (-b +std::sqrt(delta))/(2*a) ;
          result0.x = X;
          result0.y = sol1_y;
          result1.x = X;
          result1.y = sol2_y;
     } else {
          auto X = (r1*r1-r2*r2-x1*x1+x2*x2-y1*y1+y2*y2)/(2*(x2-x1));
          auto Y = (y2-y1)/(x2-x1);
          auto a = 1 + Y*Y;
          auto b = 2*Y*x2 - 2*X*Y - 2*y2;
          auto c = X*X - 2*X*x2 + x2*x2 + y2*y2 - r2*r2;
          auto delta = b*b -4*a*c;
          assert(!(delta < 0));
          auto sol1_y = (-b -std::sqrt(delta))/(2*a) ;
          auto sol2_y = (-b +std::sqrt(delta))/(2*a) ;
          auto sol1_x = X - Y*sol1_y;
          auto sol2_x = X - Y*sol2_y;
          result0.x = sol1_x;
          result0.y = sol1_y;
          result1.x = sol2_x;
          result1.y = sol2_y;
     }
     return std::make_tuple(result0, result1);
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

