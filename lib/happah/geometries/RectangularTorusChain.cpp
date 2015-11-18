// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/RectangularTorusChain.h"

namespace happah {

//NOTE: nHandles must be greater than one.
RectangularTorusChain::RectangularTorusChain(hpuint nHandles, hpreal side, hpreal spacing, hpreal hole)
     : m_hole(hole), m_nHandles(nHandles), m_side(side), m_spacing(spacing) {
     assert(nHandles > 1);//TODO: throw exception
     assert(0.0 <= hole && hole <= 1.0);//TODO: throw exception
}

}//namespace happah

