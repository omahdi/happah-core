// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometry/BezierTriangleMesh.hpp>

int main() {
     using namespace happah;

     auto solution0 = Indices({ 0, 1, 2, 3, 4 });
     auto solution1 = Indices({ 5, 6, 7, 8 });
     auto solution2 = Indices({ 9, 10, 11 });
     auto solution3 = Indices({ 12, 13 });
     auto solution4 = Indices({ 14 });
     auto solution5 = Indices({ 4, 8, 11, 13, 14 });
     auto solution6 = Indices({ 14, 12, 9, 5, 0 });
     auto answer0 = make(make_row_enumerator<0>(4, 0));
     auto answer1 = make(make_row_enumerator<0>(4, 1));
     auto answer2 = make(make_row_enumerator<0>(4, 2));
     auto answer3 = make(make_row_enumerator<0>(4, 3));
     auto answer4 = make(make_row_enumerator<0>(4, 4));
     auto answer5 = make(make_row_enumerator<1>(4, 0));
     auto answer6 = make(make_row_enumerator<2>(4, 0));

     assert(answer0 == solution0);
     assert(answer1 == solution1);
     assert(answer2 == solution2);
     assert(answer3 == solution3);
     assert(answer4 == solution4);
     assert(answer5 == solution5);
     assert(answer6 == solution6);
 
     return 0;
}


