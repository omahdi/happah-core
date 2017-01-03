// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/SurfaceSplineBEZ.h"

namespace happah {

namespace curves {

//Returns the coefficients for computing the inner control points of the Bezier representation elevated by one degree.
std::vector<hpreal> make_elevation_coefficients(hpuint degree) {
     std::vector<hpreal> coefficients;
     coefficients.reserve(degree << 1);
     auto alpha = hpreal(1.0 / (degree + 1));
     for(auto i = 1u; i <= degree; ++i) {
          coefficients.push_back(alpha * i);
          coefficients.push_back(alpha * (degree + 1 - i));
     }
     return coefficients;
}

}

namespace surfaces {

//Returns the coefficients for computing the inner control points of the Bezier representation elevated by one degree.
std::vector<hpreal> make_elevation_coefficients(hpuint degree) {
     std::vector<hpreal> coefficients;

     return coefficients;
}

}//namespace surfaces

}//namespace happah

