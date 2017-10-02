// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometry/BezierTriangleMesh.hpp>

int main() {
     using namespace happah;
     
     std::vector<Point3D> controlPoints1 = {
          Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(2, 0, 0), Point3D(3, 0, 0), Point3D(4, 0, 0),
          Point3D(0.5, 1, 0), Point3D(1.5, 1, 0), Point3D(2.5, 1, 0), Point3D(3.5, 1, 0),
          Point3D(1, 2, 0), Point3D(2, 2, 0), Point3D(3, 2, 0),
          Point3D(1.5, 3, 0), Point3D(2.5, 3, 0),
          Point3D(2, 4, 0),
          Point3D(3.5, -1, 0), Point3D(2.5, -1, 0), Point3D(1.5, -1, 0), Point3D(0.5, -1, 0),
          Point3D(3, -2, 0), Point3D(2, -2, 0), Point3D(1, -2, 0),
          Point3D(2.5, -3, 0), Point3D(1.5, -3, 0),
          Point3D(2, -4, 0)
     };
     auto indices1 = Indices();
     indices1.assign({
          hpuint(0), hpuint(1), hpuint(2), hpuint(3), hpuint(4),
          hpuint(5), hpuint(6), hpuint(7), hpuint(8),
          hpuint(9), hpuint(10), hpuint(11),
          hpuint(12), hpuint(13),
          hpuint(14),
          hpuint(4), hpuint(3), hpuint(2), hpuint(1), hpuint(0),
          hpuint(15), hpuint(16), hpuint(17), hpuint(18),
          hpuint(19), hpuint(20), hpuint(21),
          hpuint(22), hpuint(23),
          hpuint(24)
     });
     auto surface1 = BezierTriangleMesh<Space3D, 4>(controlPoints1, indices1);
     
     std::vector<Point3D> controlPoints2 = {
          Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(2, 0, 0), Point3D(3, 0, 0), Point3D(4, 0, 0),
          Point3D(0.5, 1, 0), Point3D(1.5, 1, 0), Point3D(2.5, 1, 0), Point3D(3.5, 1, 0),
          Point3D(1, 2, 0), Point3D(2, 2, 0), Point3D(3, 2, 0),
          Point3D(1.5, 3, 0), Point3D(2.5, 3, 0),
          Point3D(2, 4, 0),
          Point3D(3.5, 0, -1), Point3D(2.5, 0, -1), Point3D(1.5, 0, -1), Point3D(0.5, 0, -1),
          Point3D(3, 0, -2), Point3D(2, 0, -2), Point3D(1, 0, -2),
          Point3D(2.5, 0, -3), Point3D(1.5, 0, -3),
          Point3D(2, 0, -4)
     };
     auto indices2 = Indices();
     indices2.assign({
          hpuint(0), hpuint(1), hpuint(2), hpuint(3), hpuint(4),
          hpuint(5), hpuint(6), hpuint(7), hpuint(8),
          hpuint(9), hpuint(10), hpuint(11),
          hpuint(12), hpuint(13),
          hpuint(14),
          hpuint(4), hpuint(3), hpuint(2), hpuint(1), hpuint(0),
          hpuint(15), hpuint(16), hpuint(17), hpuint(18),
          hpuint(19), hpuint(20), hpuint(21),
          hpuint(22), hpuint(23),
          hpuint(24)
     });
     auto surface2 = BezierTriangleMesh<Space3D, 4>(controlPoints2, indices2);

     auto answer1 = is_g1(surface1, make_neighbors(surface1), 0, 0);
     auto answer2 = is_g1(surface2, make_neighbors(surface2), 0, 0);
     
     assert(answer1);
     assert(!answer2);
 
     return 0;
}