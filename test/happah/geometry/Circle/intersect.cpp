// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometry/Circle.hpp>
#include <cmath> 
#include <glm/gtc/constants.hpp>

int main() {
     using namespace happah;

     auto epsilon = EPSILON;
   
     //circles have the same center
     auto circle00 = make_circle({ 3, 4 }, 3);
     auto circle01 = make_circle({ 3, 4 }, 5);
     assert(intersect(circle00, circle01) == boost::none);

     //circles are the same
     auto circle02 = make_circle({ 3, 4 }, 5);
     auto circle03 = make_circle({ 3, 4 }, 5);
     assert(intersect(circle02, circle03) == boost::none);

     //circles are disjoint
     auto circle04 = make_circle({ 2, 2 }, 2);
     auto circle05 = make_circle({ 7, -3 }, 4);
     assert(intersect(circle04, circle05) == boost::none);

     //circles are disjoint and centers have same x coordinate
     auto circle06 = make_circle({ 7, 4 }, 2.5);
     auto circle07 = make_circle({ 7, -3 }, 4);
     assert(intersect(circle06, circle07) == boost::none);

     //circles are disjoint and centers have same y coordinate
     auto circle08 = make_circle({ -1, -3 }, 3);
     auto circle09 = make_circle({ 7, -3 }, 4);
     assert(intersect(circle08, circle09) == boost::none);
 
     //second circle contains first circle
     auto circle10 = make_circle({ 6, -2 }, 1);
     auto circle11 = make_circle({ 7, -3 }, 4);
     assert(intersect(circle10, circle11) == boost::none);

     //circles intersect with same x coordinate
     auto circle12 = make_circle({ 2, 2 }, 2);
     auto circle13 = make_circle({ 2, 6 }, 3);
     auto intersections1 = *intersect(circle12, circle13);
     auto intersection1 = boost::get<std::tuple<Point2D, Point2D> >(intersections1);
     assert(glm::length(std::get<0>(intersection1) - Point2D(2 - std::sqrt(135) / 8, 27.0 / 8.0)) < epsilon);
     assert(glm::length(std::get<1>(intersection1) - Point2D(2 + std::sqrt(135) / 8, 27.0 / 8.0)) < epsilon);

     //circles intersect with same y coordinate
     auto circle14 = make_circle({ 3, 1 }, 3);
     auto circle15 = make_circle({ -1, 1 }, 2);
     auto intersections2 = *intersect(circle14, circle15);
     auto intersection2 = boost::get<std::tuple<Point2D, Point2D> >(intersections2);
     assert(glm::length(std::get<0>(intersection2) - Point2D(3.0 / 8.0, 1 - std::sqrt(135) / 8)) < epsilon);
     assert(glm::length(std::get<1>(intersection2) - Point2D(3.0 / 8.0, 1 + std::sqrt(135) / 8)) < epsilon);

     //circles intersect
     auto circle16 = make_circle({ 1, 1 }, 1);
     auto circle17 = make_circle({ 3, 0 }, 2);     
     auto intersections3 = *intersect(circle16, circle17);
     auto intersection3 = boost::get<std::tuple<Point2D, Point2D> >(intersections3);
     assert(glm::length(std::get<0>(intersection3) - Point2D(1.8, 16.0 / 10)) < epsilon);
     assert(glm::length(std::get<1>(intersection3) - Point2D(1.0, 0.0)) < epsilon);

     //circles are tangent with same x coordinate
     auto circle18 = make_circle({ 3, 2 }, 3);
     auto circle19 = make_circle({ 3, -2 }, 1);
     auto intersections4 = *intersect(circle18, circle19);
     auto intersection4 = boost::get<Point2D>(intersections4);
     assert(glm::length(intersection4 - Point2D(3.0, -1.0)) < epsilon);

     //circles are tangent with same y coordinate
     auto circle20 = make_circle({ 1, -2 }, 2);
     auto circle21 = make_circle({-4, -2 }, 3);
     auto intersections5 = *intersect(circle20, circle21);
     auto intersection5 = boost::get<Point2D>(intersections5);
     assert(glm::length(intersection5 - Point2D(-1.0, -2.0)) < epsilon);

     //circles are tangent
     auto circle22 = make_circle({ 1, 1 }, std::sqrt(8));
     auto circle23 = make_circle({ 4, 4 }, std::sqrt(2));
     auto intersections6 = *intersect(circle22, circle23);
     auto intersection6 = boost::get<Point2D>(intersections6);
     assert(glm::length(intersection6 - Point2D(3.0, 3.0)) < epsilon);

     return 0;
}

