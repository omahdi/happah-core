// Copyright 2017
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometries/SurfaceSplineBEZ.hpp>

int main() {
     using namespace happah;

     auto solution0 = Indices();
     auto solution1 = Indices();
     auto solution2 = Indices({ 4,3,1 });
     auto solution3 = Indices({ 5,4,1, 6,5,2, 8,7,5 });
     auto solution4 = Indices({ 6,5,1, 7,6,2, 8,7,3, 10,9,6, 11,10,7, 13,12,10 });
     auto solution5 = Indices({ 7,6,1, 8,7,2, 9,8,3, 10,9,4, 12,11,7, 13,12,8, 14,13,9, 16,15,12, 17,16,13, 19,18,16 });
     auto answer0 = expand(make_nablas_enumerator(0));
     auto answer1 = expand(make_nablas_enumerator(1));
     auto answer2 = expand(make_nablas_enumerator(2));
     auto answer3 = expand(make_nablas_enumerator(3));
     auto answer4 = expand(make_nablas_enumerator(4));
     auto answer5 = expand(make_nablas_enumerator(5));
     
     assert(answer0 == solution0);
     assert(answer1 == solution1);
     assert(answer2 == solution2);
     assert(answer3 == solution3);
     assert(answer4 == solution4);
     assert(answer5 == solution5);
 
     return 0;
}