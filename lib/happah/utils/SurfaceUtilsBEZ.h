// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <vector>

#include "happah/Happah.h"
#include "happah/math/MathUtils.h"
#include "happah/math/Space.h"

namespace happah {

constexpr hpuint make_offset(hpuint degree, hpuint i0, hpuint i1, hpuint i2);

constexpr hpuint make_patch_size(hpuint degree);

class SurfaceUtilsBEZ {
public:
     template<hpuint t_degree>
     struct get_number_of_control_polygon_triangles : public std::integral_constant<hpuint, (t_degree * t_degree)> {};

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
                    indices.reserve(3 * get_number_of_control_polygon_triangles<t_degree>::value);

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

     /**
      * @param[in] nSamples Number of times an edge of the parameter triangle should be sampled.  The entire triangle is sampled uniformly such that this parameter is respected.
      * @return Matrix whose rows are the Bernstein polynomials evaluated at some point u.  The matrix is returned row-major.  To evaluate a B\'ezier polynomial at the sampled u values given a vector of control points, simply compute the product of the matrix with the vector of control points.
      */
     template<hpuint t_degree>
     static std::vector<hpreal> getEvaluationMatrix(hpuint nSamples) {
          std::vector<hpreal> matrix;
          sample(nSamples, [&] (hpreal u, hpreal v, hpreal w) {
               hpuint coefficient = 1;
               hpuint k = 0;
               hpreal wk = 1.0;
               while(k < t_degree) {
                    hpuint limit = t_degree - k;
                    hpuint i = limit;
                    hpreal ui = MathUtils::pow(u, i);
                    matrix.push_back(coefficient * ui * wk);//j=0
                    coefficient *= i;
                    --i;
                    ui = MathUtils::pow(u, i);//NOTE: Here we could also do 'ui /= u' but to avoid division by zero, we recalculate u^i.
                    hpuint j = 1;
                    hpreal vj = v;
                    while(j < limit) {
                         matrix.push_back(coefficient * ui * vj * wk);
                         coefficient *= i;
                         --i;
                         ui = MathUtils::pow(u, i);
                         ++j;
                         coefficient /= j;
                         vj *= v;
                    }
                    matrix.push_back(coefficient * vj * wk);//i=0
                    coefficient *= limit;
                    wk *= w;
                    ++k;
                    coefficient /= k;
               }
               matrix.push_back(wk);//k=degree
          });
          return std::move(matrix);
     }

     static hpuint getNumberOfControlPolygonTriangles(hpuint degree) { return degree * degree; }

     /**
      * Sample parameter triangle uniformly and pass u,v,w to visitor.
      * @param[nSamples] Number of samples on one edge of the parameter triangle.
      * 
      * Example:
      *   SurfaceUtilsBEZ::sample(3, [] (hpreal u, hpreal v, hpreal w) {
      *        std::cout << '(' << u << ',' << v << ',' << w << ")\n";
      *   });
      */
     //TODO: what if nSamples = 1?
     //TODO: move to Triangle class
     template<class Visitor>
     static void sample(hpuint nSamples, Visitor&& visit) {
          hpuint degree = nSamples - 1;
          hpuint nPoints = make_patch_size(degree);
          hpreal delta = 1.0 / degree;
          hpreal u = 1.0, v = 0.0, w = 0.0;
          hpuint rowLength = nSamples;
          hpuint limit = rowLength;
          hpreal ou = u;
          hpuint i = 0;
          while(i < nPoints) {
               visit(u, v, w);
               ++i;
               if(i == limit) {
                    --rowLength;
                    limit = i + rowLength;
                    ou -= delta;
                    u = ou;
                    v = 0.0;
                    w += delta;
               } else {
                    u -= delta;
                    v += delta;
               }
          }
     }

};//SurfaceUtilsBEZ

constexpr hpuint make_offset(hpuint degree, hpuint i0, hpuint i1, hpuint i2) { return make_patch_size(degree) - make_patch_size(degree - i2) + i1; }

constexpr hpuint make_patch_size(hpuint degree) { return (degree + 1) * (degree + 2) >> 1; }

}//namespace happah

