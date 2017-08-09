// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

//DECLARATIONS

//Evaluate the Bernstein polynomial B^n_i.
template<bool safe = true>
hpreal bernstein(hpuint degree, hpuint i, hpreal t);

template<bool safe = true>
hpuint binom(hpuint n, hpuint i);

template<bool safe = true>
hpuint munom(hpuint n, hpuint i, hpuint j);

inline hpreal pow(hpreal x, hpuint i);

//DEFINITIONS

template<bool safe>
hpreal bernstein(hpuint degree, hpuint i, hpreal t) { return binom<safe>(degree, i) * pow(t, i) * pow(hpreal(1.0) - t, degree - i); }

template<bool safe>
hpuint binom(hpuint n, hpuint i) {
     if(safe && i > n) return 0;
     if(i == n) return 1;
     auto result = hpuint(1);
     for(auto d = hpuint(1); d <= i; ++d) {
          result *= n--;
          result /= d;
     }
     return result;
}

template<bool safe = true>
hpuint munom(hpuint n, hpuint i, hpuint j) {
     if(safe && (i+j) > n) return 0;
     auto result = hpuint(1);
     for(auto d = hpuint(1); d <= i; ++d) {
          result *= n--;
          result /= d;
     }
     for(auto d = hpuint(1); d <= j; ++d) {
          result *= n--;
          result /= d;
     }
     return result;
}

inline hpreal pow(hpreal x, hpuint i) {
     if(i == 0) return 1;
     if(i == 1) return x;
     auto tmp = hpreal(pow(x, (i >> 1)));
     if((i & 1) == 0) return tmp * tmp;
     else return x * tmp * tmp;
}

}//namespace happah

