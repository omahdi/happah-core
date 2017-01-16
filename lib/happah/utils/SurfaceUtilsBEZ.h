// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/Triangle.h"
#include "happah/math/MathUtils.h"
#include "happah/math/Space.h"

namespace happah {

constexpr hpuint make_control_polygon_size(hpuint degree);

/**
 * @param[in] nSamples Number of times an edge of the parameter triangle should be sampled.  The entire triangle is sampled uniformly such that this parameter is respected.
 * @return Matrix whose rows are the Bernstein polynomials evaluated at some point u.  The matrix is returned row-major.  To evaluate a B\'ezier polynomial at the sampled u values given a vector of control points, simply compute the product of the matrix with the vector of control points.
 */
std::vector<hpreal> make_evaluation_matrix(hpuint degree, hpuint nSamples);

constexpr hpuint make_offset(hpuint degree, hpuint i0, hpuint i1, hpuint i2);

constexpr hpuint make_patch_size(hpuint degree);

class SurfaceUtilsBEZ {
public:
     template<hpuint t_degree>
     static std::vector<hpuint> buildTriangleMeshIndices() {
          switch(t_degree) {
          case 1: return { 0, 1, 2 };
          case 2: return {
                    3, 0, 1,
                    4, 1, 2,
                    1, 4, 3,
                    5, 3, 4
               };
          case 3: return {
                    1, 4, 0,
                    2, 5, 1,
                    6, 5, 2,
                    6, 2, 3,
                    1, 5, 4,
                    7, 4, 5,
                    8, 5, 6,
                    5, 8, 7,
                    9, 7, 8
               };
          case 4: return {
                    5, 0, 1,
                    6, 1, 2,
                    7, 2, 3,
                    8, 3, 4,
                    9, 5, 6,
                    10, 6, 7,
                    11, 7, 8,
                    12, 9, 10,
                    13, 10, 11,
                    14, 12, 13
                    //6, 1, 5,
                    //10, 6, 9,
                    //13, 10, 12
               };
          default: {
                    std::vector<hpuint> indices;
                    indices.reserve(3 * make_control_polygon_size(t_degree));

                    auto index = indices.begin();
                    hpuint i = 0;
                    hpuint ai = i+t_degree;
                    while(i < t_degree) {
                         *index = i;//TODO: shouldn't this be push_back
                         ++index;
                         *index = ++i;
                         ++index;
                         *index = ++ai;
                         ++index;
                    }
                    ++i;
                    hpuint d = t_degree;
                    while(d > 1) {
                         hpuint bi = i-d;
                         --d;
                         hpuint end = i+d;
                         ai = end;
                         while(i < end) {
                              hpuint oi = i;
                              ++i;

                              *index = oi;
                              ++index;
                              *index = i;
                              ++index;
                              *index = ++ai;
                              ++index;

                              *index = oi;
                              ++index;
                              *index = i;
                              ++index;
                              *index = bi++;
                              ++index;
                         }
                         ++i;
                    }//TODO: fix that vertices are in counterclockwise order
                    return indices;
                    //TODO: arrange provoking vertices correctly
                    //NOTE: If t_degree = 4, there are 16 triangle and 15 vertices which means we cannot map each vertex to a single triangle, which is necessary for flat shading (one of the vertices needs to a provoking vertex).  If t_degree = 5, there are 25 triangles but only 21 vertices.  The general solution must duplicate some vertices, namely, t_degree*t_degree-(t_degree+1)*(t_degree+2)/2=t_degree*(t_degree-3)/2-1 of them.
               }
          }
     }

};//SurfaceUtilsBEZ

constexpr hpuint make_control_polygon_size(hpuint degree) { return degree * degree; }

constexpr hpuint make_offset(hpuint degree, hpuint i0, hpuint i1, hpuint i2) { return make_patch_size(degree) - make_patch_size(degree - i2) + i1; }

constexpr hpuint make_patch_size(hpuint degree) { return (degree + 1) * (degree + 2) >> 1; }

}//namespace happah

