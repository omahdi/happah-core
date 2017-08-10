// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

/**
 * Sample parameter triangle uniformly and pass u,v,w to visitor.
 * @param[nSamples] Number of samples on one edge of the parameter triangle.
 * 
 * Example:
 *   sample(3, [] (hpreal u, hpreal v, hpreal w) {
 *        std::cout << '(' << u << ',' << v << ',' << w << ")\n";
 *   });
 */
//TODO: what if nSamples = 1?
template<class Visitor>
void sample(hpuint nSamples, Visitor&& visit) {
     hpuint degree = nSamples - 1;
     hpuint nPoints = nSamples * (nSamples + 1) >> 1;
     hpreal delta = 1.0 / degree;
     hpreal u = 1.0, v = 0.0, w = 0.0;
     hpuint rowLength = nSamples;
     hpuint limit = rowLength;
     hpreal ou = u;
     hpuint i = 0;
     while(i < nPoints) {
          visit(u, v, w);
          ++i;
          if(i == limit) {
               --rowLength;
               limit = i + rowLength;
               ou -= delta;
               u = ou;
               v = 0.0;
               w += delta;
          } else {
               u -= delta;
               v += delta;
          }
     }
}

}//namespace happah

