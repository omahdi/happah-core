// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/utils/GeometryUtils.h"

std::vector<VertexA2O1N>* GeometryUtils::sampleSineWave(hpuint xEdgeLength, hpuint yEdgeLength, hpuint nx, hpuint ny, hpreal amplitude) {
     std::vector<VertexA2O1N>* points = new std::vector<VertexA2O1N>();
     points->reserve(nx * ny);
     hpreal xDelta = xEdgeLength / hpreal(nx - 1);
     hpreal yDelta = yEdgeLength / hpreal(ny - 1);
     hpreal x = -0.5 * xEdgeLength;
     hpreal y = -0.5 * yEdgeLength;
     for(hpuint i = 0; i < nx; ++i) {
          hpreal temp = y;
          for(hpuint j = 0; j < ny; ++j) {
               Point1D ordinate(amplitude * glm::sin(y));
               Vector3D normal(0.0, amplitude * glm::cos(y), 1.0);
               points->push_back(VertexA2O1N(Point2D(x, y), ordinate, normal));
               y += yDelta;
          }
          x += xDelta;
          y = temp;
     }
     return points;
}
