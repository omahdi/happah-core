// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>

#include "happah/geometries/RectangularTorusRing.h"

namespace happah {

//NOTE: nHandles must be greater than one.
RectangularTorusRing::RectangularTorusRing(hpuint nHandles, hpreal side, hpreal spacing, hpreal hole)
     : m_hole(hole), m_nFaceVertices((nHandles << 3) + 1), m_nHandles(nHandles), m_side(side), m_spacing(spacing) {
     assert(nHandles > 1);//TODO: throw exception
     assert(0.0 <= hole && hole <= 1.0);//TODO: throw exception

     hpreal temp1 = M_PI / (hpreal) nHandles;
     hpreal sideSquared = side * side;
     hpreal temp2 = spacing * side * std::cos(temp1);
     hpreal aSquared = spacing * spacing + sideSquared + temp2 + temp2;
     hpreal temp3 = aSquared / (1.0 - std::cos(temp1 + temp1));
     hpreal baseRadiusSquared = temp3 / 2.0;
     m_baseRadius = std::sqrt(baseRadiusSquared);
     m_theta = std::acos(1.0 - sideSquared / temp3);

     hpreal temp4 = sideSquared + baseRadiusSquared + 2.0 * side * m_baseRadius * std::cos(m_theta / 2.0);
     m_r0 = (temp4 + baseRadiusSquared - sideSquared) / temp3;
     m_r1 = std::sqrt(temp4 / baseRadiusSquared - m_r0 * m_r0);
}

}//namespace happah

