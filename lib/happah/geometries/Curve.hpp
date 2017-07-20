// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

namespace curves {

template<class Point, class Iterator>
std::vector<Point> elevate(hpuint degree, const Point& first, Iterator middle, const Point& last) {
     std::vector<Point> middle1;
     middle1.reserve(degree);
     auto alpha = hpreal(1.0 / (degree + 1));
     middle1.push_back(alpha * (first + hpreal(degree) * (*middle)));
     for(auto i = 2u; i < degree; ++i, ++middle) middle1.push_back(alpha * (hpreal(i) * (*middle) + hpreal(degree + 1 - i) * (*(middle + 1))));
     middle1.push_back(alpha * (hpreal(degree) * (*middle) + last));
     assert(middle1.size() == degree);
     return middle1;
}

}//namespace curves

}//namespace happah

