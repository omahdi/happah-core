// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometries/Circle.hpp>
#include <cmath> 
#include <glm/gtc/constants.hpp>

int main() {
     using namespace happah;
     auto epsilon = pow(10, -6);
   
     //test when the two circles have the same center
     Point2D center1;
     center1.x = 3;
     center1.y = 4;
     Point2D center2;
     center2.x = 3;
     center2.y = 4;
     Circle C1 = make_circle(center1, 3);
     Circle C2 = make_circle(center2, 5);
     assert(intersect(C1, C2) == boost::none);

     //test when the two circles are the same
     C1 = make_circle(center1, 5);
     assert(intersect(C1, C2) == boost::none);

     //test when the two circles don't intersect
     center1.x = 2;
     center1.y = 2;
     C1 = make_circle(center1, 2);
     center2.x = 7;
     center2.y = -3;
     C2 = make_circle(center2, 4);
     assert(intersect(C1, C2) == boost::none);

     //test when the two circles don't intersect with the same absciss
     center1.x = 7;
     center1.y = 4;
     C1 = make_circle(center1, 2.5);
     assert(intersect(C1, C2) == boost::none);

     //test when the two circles don't intersect with the same ordinate
     center1.x = -1;
     center1.y = -3;
     C1 = make_circle(center1, 3);
     assert(intersect(C1, C2) == boost::none);
 
     //test when the two circles are within the other
     center1.x = 6;
     center1.y = -2;
     C1 = make_circle(center1, 1);
     assert(intersect(C1, C2) == boost::none);

     //test when the circles intersect with the same absciss
     center1.x = 2;
     center1.y = 2;
     C1 = make_circle(center1, 2);
     center2.x = 2;
     center2.y = 6;
     C2 = make_circle(center2,3);
     auto intersections = intersect(C1,C2);
     assert(std::abs(std::get<0>(*intersections).x - (2 - (sqrt(135)/8))) < epsilon);
     assert(std::get<0>(*intersections).y == (27.0/8.0));
     assert(std::abs(std::get<1>(*intersections).x - (2 + (sqrt(135)/8))) < epsilon);
     assert(std::get<1>(*intersections).y == (27.0/8.0));

     //test when the circles intersect with the same ordinate
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
     assert(std::get<1>(*intersections).y == 3);

     return 0;
}

