// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cmath>

#include "happah/Happah.hpp"
#include "happah/geometries/Vertex.hpp"

class GeometryUtils {
public:
     static std::vector<VertexA2O1N>* sampleSineWave(hpuint xEdgeLength, hpuint yEdgeLength, hpuint nx, hpuint ny, hpreal amplitude);

};

