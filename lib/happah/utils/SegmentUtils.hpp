// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <glm/glm.hpp>

#include "happah/Happah.hpp"

//TODO: Segment::Utils as in Plane::Utils
class SegmentUtils {
public:
     static Point2D intersect(const Point2D& p1, const Point2D& p2, hpreal y, const hpreal& epsilon = happah::EPSILON) {
          hpreal temp = p1.y - p2.y;
          if (glm::abs(temp) < epsilon) return Point2D(p1.x, y);
          else {
               hpreal alpha = (y - p2.y) / temp;
               return alpha * p1 + (1.0f - alpha) * p2;
          }
     }

     static Point2D intersect(const std::tuple<Point2D&, Point2D&>& segment, hpreal y) { return intersect(std::get<0>(segment), std::get<1>(segment), y); }

     //TODO: tuple<const Point2D&, const Point2D&>
     static bool isLeftOf(const std::tuple<Point2D&, Point2D&>& segment, const Point2D& point, const hpreal& epsilon = happah::EPSILON) { return isLeftOf(std::get<0>(segment), std::get<1>(segment), point, epsilon); }

     static bool isLeftOf(const Point2D& p1, const Point2D& p2, const Point2D& point, const hpreal& epsilon = happah::EPSILON) {
          if (p1.y < p2.y) return doIsLeftOf(p1, p2, point, epsilon);
          else return doIsLeftOf(p2, p1, point, epsilon);
     }

private:
     static bool doIsLeftOf(const Point2D& p1, const Point2D& p2, const Point2D& point, const hpreal& epsilon = happah::EPSILON) {
          Vector2D v = point - p1;
          hpreal temp = glm::dot(p2 - p1, v);
          if(glm::abs(temp) < epsilon) return v.x < 0.0;
          return temp < 0.0;
     }

};

