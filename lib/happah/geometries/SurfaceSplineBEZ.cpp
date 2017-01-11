// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/SurfaceSplineBEZ.h"

namespace happah {

hpuint make_boundary_offset(hpuint degree, hpuint i, hpuint k) {
     switch(i) {
     case 0u: return k + 1u;
     case 1u: return (degree << 1) + ((k * ((degree << 1) - k - 1u)) >> 1);
     case 2u: return make_patch_size(degree) - 3u - ((k * (5u + k)) >> 1);
     }
}

hpuint make_interior_offset(hpuint degree, hpuint i) {
     auto delta = degree - 2u;
     auto end = degree - 3u;
     auto j = degree + 2u;//absolute index
     while(i > end) {
          --delta;
          end += delta;
          j += delta + 3u;
     }
     return j + i - (end - delta + 1u);
}

}//namespace happah

