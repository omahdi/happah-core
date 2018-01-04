// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometry/BezierTriangleMesh.hpp>

int main() {
     using namespace happah;
     
     auto controlPoints1 = std::vector<Point3D>({
          Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(2, 0, 0), Point3D(3, 0, 0), Point3D(4, 0, 0),
          Point3D(0.5, 1, 0), Point3D(1.5, 1, 0), Point3D(2.5, 1, 0), Point3D(3.5, 1, 0),
          Point3D(1, 2, 0), Point3D(2, 2, 0), Point3D(3, 2, 0),
          Point3D(1.5, 3, 0), Point3D(2.5, 3, 0),
          Point3D(2, 4, 0),
          Point3D(3.5, -1, 0), Point3D(2.5, -1, 0), Point3D(1.5, -1, 0), Point3D(0.5, -1, 0),
          Point3D(3, -2, 0), Point3D(2, -2, 0), Point3D(1, -2, 0),
          Point3D(2.5, -3, 0), Point3D(1.5, -3, 0),
          Point3D(2, -4, 0)
     });
     auto indices1 = Tuples<hpindex>(15, {
          hpindex(0), hpindex(1), hpindex(2), hpindex(3), hpindex(4),
          hpindex(5), hpindex(6), hpindex(7), hpindex(8),
          hpindex(9), hpindex(10), hpindex(11),
          hpindex(12), hpindex(13),
          hpindex(14),
          hpindex(4), hpindex(3), hpindex(2), hpindex(1), hpindex(0),
          hpindex(15), hpindex(16), hpindex(17), hpindex(18),
          hpindex(19), hpindex(20), hpindex(21),
          hpindex(22), hpindex(23),
          hpindex(24)
     });
     auto mesh1 = BezierTriangleMesh<Space3D, 4>(controlPoints1, indices1);
     
     auto controlPoints2 = std::vector<Point3D>({
          Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(2, 0, 0), Point3D(3, 0, 0), Point3D(4, 0, 0),
          Point3D(0.5, 1, 0), Point3D(1.5, 1, 0), Point3D(2.5, 1, 0), Point3D(3.5, 1, 0),
          Point3D(1, 2, 0), Point3D(2, 2, 0), Point3D(3, 2, 0),
          Point3D(1.5, 3, 0), Point3D(2.5, 3, 0),
          Point3D(2, 4, 0),
          Point3D(3.5, 0, -1), Point3D(2.5, 0, -1), Point3D(1.5, 0, -1), Point3D(0.5, 0, -1),
          Point3D(3, 0, -2), Point3D(2, 0, -2), Point3D(1, 0, -2),
          Point3D(2.5, 0, -3), Point3D(1.5, 0, -3),
          Point3D(2, 0, -4)
     });
     auto indices2 = Tuples<hpindex>(15, {
          hpindex(0), hpindex(1), hpindex(2), hpindex(3), hpindex(4),
          hpindex(5), hpindex(6), hpindex(7), hpindex(8),
          hpindex(9), hpindex(10), hpindex(11),
          hpindex(12), hpindex(13),
          hpindex(14),
          hpindex(4), hpindex(3), hpindex(2), hpindex(1), hpindex(0),
          hpindex(15), hpindex(16), hpindex(17), hpindex(18),
          hpindex(19), hpindex(20), hpindex(21),
          hpindex(22), hpindex(23),
          hpindex(24)
     });
     auto mesh2 = BezierTriangleMesh<Space3D, 4>(controlPoints2, indices2);

     auto answer1 = is_g1(mesh1, make_neighbors(mesh1), 0, TRIT0);
     auto answer2 = is_g1(mesh2, make_neighbors(mesh2), 0, TRIT0);
     
     assert(answer1);
     assert(!answer2);
 
     return 0;
}
