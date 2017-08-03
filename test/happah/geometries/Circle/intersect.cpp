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
     auto intersections01 = *intersect(circle12, circle13);
     auto intersection01 = boost::get<std::tuple<Point2D, Point2D> >(intersections01);
     assert(glm::length(std::get<0>(intersection01) - Point2D(2 - std::sqrt(135) / 8, 27.0 / 8.0)) < epsilon);
     assert(glm::length(std::get<1>(intersection01) - Point2D(2 + std::sqrt(135) / 8, 27.0 / 8.0)) < epsilon);

     //circles intersect with same y coordinate
     auto circle14 = make_circle({ 3, 1 }, 3);
     auto circle15 = make_circle({ -1, 1 }, 2);
     auto intersections02 = *intersect(circle14, circle15);
     auto intersection02 = boost::get<std::tuple<Point2D, Point2D> >(intersections02);
     assert(glm::length(std::get<0>(intersection02) - Point2D(3.0 / 8.0, 1 - sqrt(135) / 8)) < epsilon);
     assert(glm::length(std::get<1>(intersection02) - Point2D(3.0 / 8.0, 1 + sqrt(135) / 8)) < epsilon);

    //circles intersect
     auto circle16 = make_circle({ 1, 1 }, 1);
     auto circle17 = make_circle({ 3, 0 }, 2);     
     auto intersections03 = *intersect(circle16, circle17);
     auto intersection03 = boost::get<std::tuple<Point2D, Point2D> >(intersections03);
     assert(glm::length(std::get<1>(intersection03) - Point2D(1.0, 0.0)) < epsilon);
     assert(glm::length(std::get<0>(intersection03) - Point2D(1.8, 16.0/10)) < epsilon);

     //circles are tangent with same x coordinate
     auto circle18 = make_circle({ 3, 2 }, 3);
     auto circle19 = make_circle({ 3, -2 }, 1);
     auto intersections04 = *intersect(circle18, circle19);
     auto intersection04 = boost::get<Point2D>(intersections04);
     assert(glm::length(intersection04 - Point2D(3.0, -1.0)) < epsilon);

     //circles are tangent with same y coordinate
     auto circle20 = make_circle({ 1, -2 }, 2);
     auto circle21 = make_circle({-4, -2 }, 3);
     auto intersections05 = *intersect(circle20, circle21);
     auto intersection05 = boost::get<Point2D>(intersections05);
     assert(glm::length(intersection05 - Point2D(-1.0, -2.0)) < epsilon);

     //circles are tangent
     auto circle22 = make_circle({ 1, 1 }, sqrt(8));
     auto circle23 = make_circle({ 4, 4 }, sqrt(2));
     auto intersections06 = *intersect(circle22, circle23);
     auto intersection06 = boost::get<Point2D>(intersections06);
     assert(glm::length(intersection06 - Point2D(3.0, 3.0)) < epsilon);

     return 0;
}

