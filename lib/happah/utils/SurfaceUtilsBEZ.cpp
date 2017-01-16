// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/utils/SurfaceUtilsBEZ.h"

namespace happah {

std::vector<hpreal> make_de_casteljau_matrix(hpuint degree, hpuint nSamples) {
     std::vector<hpreal> matrix;
     sample(nSamples, [&](hpreal u, hpreal v, hpreal w) {
          hpuint coefficient = 1;
          hpuint k = 0;
          hpreal wk = 1.0;
          while(k < degree) {
               hpuint limit = degree - k;
               hpuint i = limit;
               hpreal ui = MathUtils::pow(u, i);
               matrix.push_back(coefficient * ui * wk);//j=0
               coefficient *= i;
               --i;
               ui = MathUtils::pow(u, i);//NOTE: Here we could also do 'ui /= u' but to avoid division by zero, we recalculate u^i.
               hpuint j = 1;
               hpreal vj = v;
               while(j < limit) {
                    matrix.push_back(coefficient * ui * vj * wk);
                    coefficient *= i;
                    --i;
                    ui = MathUtils::pow(u, i);
                    ++j;
                    coefficient /= j;
                    vj *= v;
               }
               matrix.push_back(coefficient * vj * wk);//i=0
               coefficient *= limit;
               wk *= w;
               ++k;
               coefficient /= k;
          }
          matrix.push_back(wk);//k=degree
     });
     return std::move(matrix);
}

}//namespace happah

