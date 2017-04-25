// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/formats/off.h"

namespace happah {

namespace off {

hpuint make_vertex_size(const Header& header) {
     auto n = header.dimension;
     if(header.normal) n <<= 1;
     if(header.texture) n += 2;
     if(header.color) n += 4;
     return n;
}

}//namespace off

}//namespace happah

