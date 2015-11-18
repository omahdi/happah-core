// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/weighers/ConstantWeigher.h"

namespace happah {

auto ConstantWeigher::weigh(hpuint v0, hpuint v1) const -> Weight { return 1; }

auto ConstantWeigher::weigh(hpuint v0, hpuint v1, hpuint edge) const -> Weight { return weigh(v0, v1); }

}//namespace happah

