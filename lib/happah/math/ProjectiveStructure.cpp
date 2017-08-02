// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/ProjectiveStructure.hpp"

namespace happah {

bool validate(const ProjectiveStructure& structure, hpreal epsilon) {
     return false;
     /*auto is_one = [&](auto a) { return glm::abs(1.0 - a) < epsilon; };
     auto is_zero = [&](auto a) { return glm::abs(a) < epsilon; };
     return !find_fan(neighbors, [&](auto p, auto i, auto fan) -> auto {
          auto A1 = hpvec3(0.0, 1.0, 0.0);
          auto A2 = hpvec3(0.0, 0.0, 1.0);
          while(fan) {
               auto q = 0u, j = 0u;
               std::tie(q, j) = *fan;
               auto transition = std::begin(transitions) + (9 * q + 3 * j);
               auto temp = A2;
               if(j == 0) {
                    A2 = transition[0] * A1 + transition[2] * A2;
                    A2.x += transition[1];
               } else if(j == 1) {
                    A2 = transition[0] * A2 + transition[1] * A1;
                    A2.x += transition[2];
               } else {
                    A2 = transition[1] * A2 + transition[2] * A1;
                    A2.x += transition[0];
               }
               A1 = temp;
               ++fan;
          }
          return !(is_zero(A1[0]) && is_one(A1[1]) && is_zero(A1[2]) && is_zero(A2[0]) && is_zero(A2[1]) && is_one(A2[2]));
     });*/
}

}//namespace happah

