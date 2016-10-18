// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.h"
#include "happah/utils/SurfaceUtilsBEZ.h"

namespace happah {

class SurfaceSplineUtilsBEZ {
public:
     enum class Mode { CONTROL_POINTS, NEIGHBORS };

     template<hpuint t_degree>
     static std::vector<hpuint> buildTriangleMeshIndices(const std::vector<hpuint>& controlPointIndices) {
          std::vector<hpuint> indices;
          indices.reserve(3 * controlPointIndices.size() * SurfaceUtilsBEZ::get_number_of_control_polygon_triangles<t_degree>::value / SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value);

          switch(t_degree) {
          case 1:
               for(auto i = controlPointIndices.cbegin(), end = controlPointIndices.cend(); i != end; ++i) {
                    hpuint c0 = *i;
                    hpuint c1 = *(++i);
                    hpuint c2 = *(++i);

                    hpuint temp[3] = { c0, c1, c2 };
                    indices.insert(indices.end(), temp, temp+3);
               }
               break;
          case 2:
               for(auto i = controlPointIndices.cbegin(), end = controlPointIndices.cend(); i != end; ++i) {
                    hpuint c0 = *i;
                    hpuint c1 = *(++i);
                    hpuint c2 = *(++i);
                    hpuint c3 = *(++i);
                    hpuint c4 = *(++i);
                    hpuint c5 = *(++i);

                    hpuint temp[12] = {
                         c3, c0, c1,
                         c4, c1, c2,
                         c1, c4, c3,
                         c5, c3, c4
                    };
                    indices.insert(indices.end(), temp, temp+12);
               }
               break;
          case 3:
               for(auto i = controlPointIndices.cbegin(), end = controlPointIndices.cend(); i != end; ++i) {
                    hpuint c0 = *i;
                    hpuint c1 = *(++i);
                    hpuint c2 = *(++i);
                    hpuint c3 = *(++i);
                    hpuint c4 = *(++i);
                    hpuint c5 = *(++i);
                    hpuint c6 = *(++i);
                    hpuint c7 = *(++i);
                    hpuint c8 = *(++i);
                    hpuint c9 = *(++i);

                    hpuint temp[27] = {
                         c1, c4, c0,
                         c2, c5, c1,
                         c6, c5, c2,
                         c6, c2, c3,
                         c1, c5, c4,
                         c7, c4, c5,
                         c8, c5, c6,
                         c5, c8, c7,
                         c9, c7, c8
                    };
                    indices.insert(indices.end(), temp, temp+27);
               }
               break;
          case 4:
               for(auto i = controlPointIndices.cbegin(), end = controlPointIndices.cend(); i != end; ++i) {
                    hpuint c0 = *i;
                    hpuint c1 = *(++i);
                    hpuint c2 = *(++i);
                    hpuint c3 = *(++i);
                    hpuint c4 = *(++i);
                    hpuint c5 = *(++i);
                    hpuint c6 = *(++i);
                    hpuint c7 = *(++i);
                    hpuint c8 = *(++i);
                    hpuint c9 = *(++i);
                    hpuint c10 = *(++i);
                    hpuint c11 = *(++i);
                    hpuint c12 = *(++i);
                    hpuint c13 = *(++i);
                    hpuint c14 = *(++i);
               
                    hpuint temp[48] = {
                         c5, c0, c1,
                         c6, c1, c2,
                         c7, c2, c3,
                         c8, c3, c4,
                         c9, c5, c6,
                         c10, c6, c7,
                         c11, c7, c8,
                         c12, c9, c10,
                         c13, c10, c11,
                         c14, c12, c13,
                         c6, c1, c5,
                         c10, c6, c9,
                         c13, c10, c12,
                         c7, c2, c6,
                         c11, c7, c10,
                         c8, c3, c7
                    };
                    indices.insert(indices.end(), temp, temp+48);
               }
               break;
          case 5:
               for(auto i = controlPointIndices.cbegin(), end = controlPointIndices.cend(); i != end; ++i) {
                    hpuint c0 = *i;
                    hpuint c1 = *(++i);
                    hpuint c2 = *(++i);
                    hpuint c3 = *(++i);
                    hpuint c4 = *(++i);
                    hpuint c5 = *(++i);
                    hpuint c6 = *(++i);
                    hpuint c7 = *(++i);
                    hpuint c8 = *(++i);
                    hpuint c9 = *(++i);
                    hpuint c10 = *(++i);
                    hpuint c11 = *(++i);
                    hpuint c12 = *(++i);
                    hpuint c13 = *(++i);
                    hpuint c14 = *(++i);
                    hpuint c15 = *(++i);
                    hpuint c16 = *(++i);
                    hpuint c17 = *(++i);
                    hpuint c18 = *(++i);
                    hpuint c19 = *(++i);
                    hpuint c20 = *(++i);

                    hpuint temp[75] = {
                         c0, c1, c6,
                         c1, c2, c7,
                         c2, c3, c8,
                         c3, c4, c9,
                         c4, c5, c10,
                         c1, c7, c6,
                         c2, c8, c7,
                         c3, c9, c8,
                         c4, c10, c9,
                         c6, c7, c11,
                         c7, c8, c12,
                         c8, c9, c13,
                         c9, c10, c14,
                         c7, c12, c11,
                         c8, c13, c12,
                         c9, c14, c13,
                         c11, c12, c15,
                         c12, c13, c16,
                         c13, c14, c17,
                         c12, c16, c15,
                         c13, c17, c16,
                         c15, c16, c18,
                         c16, c17, c19,
                         c16, c19, c18,
                         c18, c19, c20,
                    };
                    indices.insert(indices.end(), temp, temp+75);
               } 
               break;
          default: std::cerr << "ERROR: Not implemented yet.\n";//TODO
          }

          return std::move(indices);
     }

};//SurfaceSplineUtilsBEZ

}//namespace happah

