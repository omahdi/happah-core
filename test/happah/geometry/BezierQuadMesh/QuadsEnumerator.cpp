// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometry/BezierQuadMesh.hpp>

int main() {
     using namespace happah;

     /*1x1*/
     auto solution1 = Indices({ 0,1,3,2 });
     /*2x2*/
     auto solution2 = Indices({ 0,1,4,3, 1,2,5,4, 3,4,7,6, 4,5,8,7 });
     /*3x3*/
     auto solution3 = Indices({ 0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10, 8,9,13,12, 9,10,14,13, 10,11,15,14 });
     /*4x4*/
     auto solution4 = Indices({ 0,1,6,5, 1,2,7,6, 2,3,8,7, 3,4,9,8, 5,6,11,10, 6,7,12,11, 7,8,13,12, 8,9,14,13, 10,11,16,15, 11,12,17,16, 12,13,18,17, 13,14,19,18, 15,16,21,20, 16,17,22,21, 17,18,23,22, 18,19,24,23 });
     /*3x2*/
     auto solution5 = Indices({ 0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10 });
     /*3x4*/
     auto solution6 = Indices({ 0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10, 8,9,13,12, 9,10,14,13, 10,11,15,14, 12,13,17,16, 13,14,18,17, 14,15,19,18 });
     /*4x1*/
     auto solution7 = Indices({ 0,1,6,5, 1,2,7,6, 2,3,8,7, 3,4,9,8 });
     
     auto answer1 = expand(make_quads_enumerator(1, 1));
     auto answer2 = expand(make_quads_enumerator(2, 2));
     auto answer3 = expand(make_quads_enumerator(3, 3));
     auto answer4 = expand(make_quads_enumerator(4, 4));
     auto answer5 = expand(make_quads_enumerator(3, 2));
     auto answer6 = expand(make_quads_enumerator(3, 4));
     auto answer7 = expand(make_quads_enumerator(4, 1));
     
     assert(answer1 == solution1);
     assert(answer2 == solution2);
     assert(answer3 == solution3);
     assert(answer4 == solution4);
     assert(answer5 == solution5);
     assert(answer6 == solution6);
     assert(answer7 == solution7);
 
     return 0;
}
