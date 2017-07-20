// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

struct ConstantWeigher {
     using Weight = hpuint;
     static const Weight MAX_WEIGHT = -1;

     Weight weigh(hpuint v0, hpuint v1) const;

     Weight weigh(hpuint v0, hpuint v1, hpuint edge) const;

};//ConstantWeigher

}//namespace happah

