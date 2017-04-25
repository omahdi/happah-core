// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.h"

namespace happah {

namespace off {

struct Header {
     bool color;
     hpuint dimension;
     hpuint nEdges;
     hpuint nFaces;
     bool normal;
     hpuint nVertices;
     bool texture;

};//Header

using Vertices = std::vector<hpreal>;

struct Faces {
     Indices vertices;
     std::vector<hpreal> colors;
     Indices indices;

};

struct Content {
     Header header;
     Vertices vertices;
     Faces faces;

};

}//namespace off

}//namespace happah

