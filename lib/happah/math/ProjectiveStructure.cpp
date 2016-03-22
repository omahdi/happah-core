// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/ProjectiveStructure.h"

namespace happah {

/*const ProjectiveStructure3D ProjectiveStructure3D::DOUBLE_TORUS_PROJECTIVE_STRUCTURE = ProjectiveStructure3D(//1
     std::vector<Transition>({
          Transition(1, 1, -1),
          Transition(2, 2, -1)
     }),
     {{
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1
     }},
     {{
          1, 8, 7,
          0, 2, 9,
          3, 12, 1,
          2, 4, 13,
          5, 16, 3,
          4, 6, 17,
          7, 20, 5,
          6, 0, 21,
          23, 0, 22,
          1, 10, 11,
          9, 11, 23,
          10, 12, 9,
          11, 2, 21,
          3, 14, 20,
          13, 15, 16,
          14, 16, 18,
          15, 4, 14,
          5, 18, 19,
          17, 19, 15,
          18, 20, 17,
          19, 6, 13,
          7, 22, 12,
          21, 23, 8,
          22, 8, 10
     }},
     {{
          9, 11,
          10, 23,
          12, 21,
          13, 20,
          14, 16,
          15, 18,
          17, 19,
          22, 8
     }}
);*/

/*const ProjectiveStructure3D ProjectiveStructure3D::DOUBLE_TORUS_PROJECTIVE_STRUCTURE = ProjectiveStructure3D(//2
     std::vector<Transition>({
          Transition(1, 1, -1),
          Transition(2, 2, -1)
     }),
     {{
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          1, 0, 0,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1,
          0, 0, 1
     }},
     {{
          1, 8, 7,
          0, 2, 9,
          3, 12, 1,
          2, 4, 13,
          5, 16, 3,
          4, 6, 17,
          7, 20, 5,
          6, 0, 21,
          23, 0, 17,
          1, 10, 16,
          9, 11, 23,
          10, 12, 14,
          11, 2, 21,
          3, 14, 20,
          13, 15, 11,
          14, 16, 18,
          15, 4, 9,
          5, 18, 8,
          17, 19, 15,
          18, 20, 22,
          19, 6, 13,
          7, 22, 12,
          21, 23, 19,
          22, 8, 10
     }},
     {{
          9, 16,
          10, 23,
          11, 14,
          12, 21,
          13, 20,
          15, 18,
          17, 8,
          22, 19
     }}
);*/

const ProjectiveStructure3D ProjectiveStructure3D::DOUBLE_TORUS_PROJECTIVE_STRUCTURE = ProjectiveStructure3D(//3
     std::vector<Transition>({
          Transition(1, 1, -1),
          Transition(2, 2, -1)
     }),
     {{
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0,
          0, 1, 0
     }},
     {{
          5, 6, 1,
          0, 17, 2,
          1, 10, 3,
          2, 15, 4,
          3, 8, 5,
          4, 13, 0,
          11, 0, 7,
          6, 23, 8,
          7, 4, 9,
          8, 21, 10,
          9, 2, 11,
          10, 19, 6,
          17, 18, 13,
          12, 5, 14,
          13, 22, 15,
          14, 3, 16,
          15, 20, 17,
          16, 1, 12,
          23, 12, 19,
          18, 11, 20,
          19, 16, 21,
          20, 9, 22,
          21, 14, 23,
          22, 7, 18
     }},
     {{
          2, 10,
          3, 15,
          4, 8,
          5, 13,
          7, 23,
          9, 21,
          14, 22,
          16, 20
     }}
);

}//namespace happah

