// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

class MathUtils {
public:
     template<bool t_check = true>
     static hpuint binom(hpuint n, hpuint i) {
          if(t_check && i > n) return 0;
          if(i == n) return 1;
          hpuint result = 1;
          for(hpuint d = 1; d <= i; ++d) {
               result *= n--;
               result /= d;
          }
          return result;
     }

     template<bool t_check = true>
     static hpuint munom(hpuint n, hpuint i, hpuint j) {
          if(t_check && (i+j) > n) return 0;
          hpuint result = 1;
          for(hpuint d = 1; d <= i; ++d) {
               result *= n--;
               result /= d;
          }
          for(hpuint d = 1; d <= j; ++d) {
               result *= n--;
               result /= d;
          }
          return result;
     }

     static hpreal pow(hpreal x, hpuint i) {
          if(i == 0) return 1;
          if(i == 1) return x;
          hpreal tmp = pow(x, (i >> 1));
          if ((i & 1) == 0) return tmp * tmp;
          else return x * tmp * tmp;
     }

};

}//namespace happah

