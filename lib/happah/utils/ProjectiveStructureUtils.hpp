// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/Happah.hpp"
#include "happah/math/Space.hpp"

namespace happah {

class ProjectiveStructureUtils {
public:
     static std::pair<std::vector<Point3D>, std::vector<hpuint> > toTransitions(const std::vector<Point3D>& points, const std::vector<hpuint>& pointIndices, const std::vector<hpuint>& neighbors, hpreal epsilon = EPSILON) {
          std::vector<Point3D> transitions;
          std::vector<hpuint> transitionIndices;

          auto insert = [&] (const Point3D& p0, const Point3D& p1, const Point3D& p2, hpuint i0, hpuint i1, hpuint i2, hpuint n) {
               if(n == UNULL) {
                    transitionIndices.push_back(UNULL);
                    return;
               }
               auto q = pointIndices.cbegin() + 3 * n;
               hpuint j0 = *q;
               hpuint i3;
               if(j0 == i0) i3 = *(q+1);
               else if(j0 == i1) i3 = *(q+2);
               else i3 = j0;
               hpmat3x3 m(p1, p0, p2);
               Point3D transition = glm::inverse(m) * points[i3];
               auto result = std::find_if(transitions.begin(), transitions.end(), [&] (const Point3D& other) -> bool { return (std::abs(transition.x - other.x) < happah::EPSILON) && (std::abs(transition.y - other.y) < happah::EPSILON) && (std::abs(transition.z - other.z) < happah::EPSILON); });
               if(result != transitions.end()) transitionIndices.push_back(std::distance(transitions.begin(), result));
               else {
                    transitions.push_back(transition);
                    transitionIndices.push_back(transitions.size() - 1);
               }
          };

          auto p = pointIndices.cbegin();
          for(auto n = neighbors.cbegin(), end = neighbors.cend(); n != end; ++n) {
               hpuint n0 = *n;
               hpuint n1 = *(++n);
               hpuint n2 = *(++n);
               hpuint i0 = *p;
               hpuint i1 = *(++p);
               hpuint i2 = *(++p);
               ++p;
               const Point3D& p0 = points[i0];
               const Point3D& p1 = points[i1];
               const Point3D& p2 = points[i2];
               insert(p0, p1, p2, i0, i1, i2, n0);
               insert(p1, p2, p0, i1, i2, i0, n1);
               insert(p2, p0, p1, i2, i0, i1, n2);
          }

          return {std::move(transitions), std::move(transitionIndices)};
     }

};//ProjectiveStructureUtils

}//namespace happah

