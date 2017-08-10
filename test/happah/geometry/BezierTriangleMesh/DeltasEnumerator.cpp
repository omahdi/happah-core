// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometry/SurfaceSplineBEZ.hpp>

int main() {
     using namespace happah;

     auto solution0 = Indices();
     auto solution1 = Indices({ 0,1,2 });
     auto solution2 = Indices({ 0,1,3,  1,2,4,  3,4,5 });
     auto solution3 = Indices({ 0,1,4,  1,2,5,  2,3,6,  4,5,7,  5,6,8,  7,8,9 });
     auto solution4 = Indices({ 0,1,5,  1,2,6,  2,3,7,  3,4,8,  5,6,9,  6,7,10,  7,8,11,  9,10,12,  10,11,13,  12,13,14 });
     auto solution5 = Indices({ 0,1,6,  1,2,7,  2,3,8,  3,4,9,  4,5,10,  6,7,11,  7,8,12,  8,9,13,  9,10,14,  11,12,15,  12,13,16,  13,14,17,  15,16,18,  16,17,19,  18,19,20 });
     auto answer0 = expand(make_deltas_enumerator(0));
     auto answer1 = expand(make_deltas_enumerator(1));
     auto answer2 = expand(make_deltas_enumerator(2));
     auto answer3 = expand(make_deltas_enumerator(3));
     auto answer4 = expand(make_deltas_enumerator(4));
     auto answer5 = expand(make_deltas_enumerator(5));
     
     assert(answer0 == solution0);
     assert(answer1 == solution1);
     assert(answer2 == solution2);
     assert(answer3 == solution3);
     assert(answer4 == solution4);
     assert(answer5 == solution5);
 
     return 0;
}

