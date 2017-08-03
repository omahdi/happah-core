// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometries/Circle.hpp>
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
     auto intersections = *intersect(circle12, circle13);
     auto intersection = boost::get<std::tuple<Point2D, Point2D> >(intersections);
     assert(glm::length(std::get<0>(intersection) - Point2D(2 - std::sqrt(135) / 8, 27.0 / 8.0)) < epsilon);
     assert(glm::length(std::get<1>(intersection) - Point2D(2 + std::sqrt(135) / 8, 27.0 / 8.0)) < epsilon);

     /*//test when the circles intersect with the same ordinate
     center1.x = 3;
     center1.y = 1;
     C1 = make_circle(center1, 3);
     center2.x = -1;
     center2.y = 1;
     C2 = make_circle(center2,2);
     intersections = intersect(C1,C2);
     assert(std::get<0>(*intersections).x == (3.0/8.0));
     assert(std::abs(std::get<0>(*intersections).y - (1 - (sqrt(135)/8))) < epsilon);
     assert(std::get<1>(*intersections).x == (3.0/8.0));
     assert(std::abs(std::get<1>(*intersections).y - (1.0 + (std::sqrt(135.0)/8.0))) < epsilon);

    //test when the circles intersect
     center1.x = 1;
     center1.y = 1;
     C1 = make_circle(center1, 1);
     center2.x = 3;
     center2.y = 0;
     C2 = make_circle(center2,2);
     intersections = intersect(C1,C2);
     assert(std::get<0>(*intersections).x == 1.0);
     assert(std::get<0>(*intersections).y == 0.0);
     assert(std::abs(std::get<1>(*intersections).x - 1.8) < epsilon);
     assert(std::abs(std::get<1>(*intersections).y - (16.0/10.0)) < epsilon);

     //test when the circles are tangent with the same abscisse
     center1.x = 3;
     center1.y = 2;
     C1 = make_circle(center1, 3);
     center2.x = 3;
     center2.y = -2;
     C2 = make_circle(center2, 1);
     intersections = intersect(C1,C2);
     assert(std::get<0>(*intersections).x == 3);
     assert(std::get<0>(*intersections).y == -1);
     assert(std::get<1>(*intersections).x == 3);
     assert(std::get<1>(*intersections).y == -1);

     //test when the circles are tangent with the same ordinate
     center1.x = 1;
     center1.y = -2;
     C1 = make_circle(center1, 2);
     center2.x = -4;
     center2.y = -2;
     C2 = make_circle(center2, 3);
     intersections = intersect(C1,C2);
     assert(std::get<0>(*intersections).x == -1);
     assert(std::get<0>(*intersections).y == -2);
     assert(std::get<1>(*intersections).x == -1);
     assert(std::get<1>(*intersections).y == -2);

     //test when the circles are tangent
     center1.x = 1;
     center1.y = 1;
     C1 = make_circle(center1, sqrt(8));
     center2.x = 4;
     center2.y = 4;
     C2 = make_circle(center2, sqrt(2));
     intersections = intersect(C1,C2);
     assert(std::get<0>(*intersections).x == 3);
     assert(std::get<0>(*intersections).y == 3);
     assert(std::get<1>(*intersections).x == 3);
     assert(std::get<1>(*intersections).y == 3);*/

     return 0;
}

