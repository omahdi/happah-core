// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/TriangleMeshUtils.h"

#include <unordered_map>
#include <utility>

#include "happah/utils/visitors.h"

namespace happah {

hpuint TriangleMeshUtils::findTriangle(const std::vector<hpuint>& indices, hpuint v0, hpuint v1, hpuint& v2) {
     hpuint triangle = 0;
     for(auto i = indices.begin(), end = indices.end(); i != end; ++i) {
          hpuint i0 = *i;
          hpuint i1 = *(++i);
          hpuint i2 = *(++i);
          if(i0 == v0) {
               if(i1 == v1) {
                    v2 = i2;
                    return triangle;
               } else if(i2 == v1) {
                    v2 = i1;
                    return triangle;
               }
          } else if(i0 == v1) {
               if(i1 == v0) {
                    v2 = i2;
                    return triangle;
               } else if(i2 == v0) {
                    v2 = i1;
                    return triangle;
               }
          } else if((i1 == v0 && i2 == v1) || (i1 == v1 && i2 == v0)) {
               v2 = i0;
               return triangle;
          }
          ++triangle;
     }
     return happah::UNULL;
}

std::tuple<hpuint, hpuint, hpuint> TriangleMeshUtils::getNeighbors(const std::vector<hpuint>& indices, hpuint triangle) {
     hpuint n1 = happah::UNULL;
     hpuint n2 = happah::UNULL;
     hpuint n3 = happah::UNULL;
     hpuint temp = 3 * triangle;
     std::vector<hpuint>::const_iterator v = indices.cbegin() + temp;
     hpuint v1 = *v;
     hpuint v2 = *(++v);
     hpuint v3 = *(++v);

     hpuint t = -1;
     for(std::vector<hpuint>::const_iterator i = indices.cbegin(); i != indices.cend(); ++i) {
          ++t;
          if(t == triangle) {
               ++(++i);
               continue;
          }
          hpuint w1 = *i;
          hpuint w2 = *(++i);
          hpuint w3 = *(++i);
          if((v1 == w1 && v2 == w2) || (v1 == w2 && v2 == w1) || (v1 == w2 && v2 == w3) || (v1 == w3 && v2 == w2) || (v1 == w1 && v2 == w3) || (v1 == w3 && v2 == w1)) n1 = t;
          if((v2 == w1 && v3 == w2) || (v2 == w2 && v3 == w1) || (v2 == w2 && v3 == w3) || (v2 == w3 && v3 == w2) || (v2 == w1 && v3 == w3) || (v2 == w3 && v3 == w1)) n2 = t;
          if((v3 == w1 && v1 == w2) || (v3 == w2 && v1 == w1) || (v3 == w2 && v1 == w3) || (v3 == w3 && v1 == w2) || (v3 == w1 && v1 == w3) || (v3 == w3 && v1 == w1)) n3 = t;

          if(n1 != happah::UNULL && n2 != happah::UNULL && n3 != happah::UNULL) break;
     }

     return std::make_tuple(n1, n2, n3);
}

}//namespace happah

