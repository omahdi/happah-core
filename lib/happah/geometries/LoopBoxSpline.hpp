// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>

#include "happah/geometries/TriangleMesh.hpp"

namespace happah {

template<class Vertex>
class LoopBoxSpline {
public:
     LoopBoxSpline(std::array<Vertex, 12> controlPoints)
          : m_controlPoints(std::move(controlPoints)) {}

     auto& getControlPoints() const { return m_controlPoints; }

private:
     std::array<Vertex, 12> m_controlPoints;

};//LoopBoxSpline

//    10  11
//   7   8   9
// 3   4   5   6
//   0   1   2
template<class Vertex>
TriangleMesh<Vertex> make_control_polygon(const LoopBoxSpline<Vertex>& surface) {
     auto& controlPoints = surface.getControlPoints();
     return make_triangle_mesh({ std::begin(controlPoints), std::end(controlPoints) }, { 0, 1, 4, 1, 2, 5, 3, 0, 4, 4, 1, 5, 5, 2, 6, 3, 4, 7, 4, 5, 8, 5, 6, 9, 7, 4, 8, 8, 5, 9, 7, 8, 10, 8, 9, 11, 10, 8, 11 });
}

//NOTE: Merge the rings of three vertices of valence six into a regular patch.
template<class Iterator, class T = typename std::iterator_traits<Iterator>::value_type>
static std::array<T, 12> make_loop_box_spline_control_points(T v0, Iterator begin0, T v1, Iterator begin1, T v2, Iterator begin2) {
     auto end0 = begin0 + 6;
     auto end1 = begin1 + 6;
     auto end2 = begin2 + 6;
     auto i0 = begin0;
     while(*i0 != v2) ++i0;
     ++i0;
     if(i0 == end0) i0 = begin0;
     auto a7 = *i0; 
     ++i0;
     if(i0 == end0) i0 = begin0;
     auto a3 = *i0;
     ++i0;
     if(i0 == end0) i0 = begin0;
     auto a0 = *i0;
     ++i0;
     if(i0 == end0) i0 = begin0;
     auto a1 = *i0;
     auto i1 = begin1;
     while(*i1 != a1) ++i1;
     ++i1;
     if(i1 == end1) i1 = begin1;
     auto a2 = *i1;
     ++i1;
     if(i1 == end1) i1 = begin1;
     auto a6 = *i1;
     ++i1;
     if(i1 == end1) i1 = begin1;
     auto a9 = *i1;
     auto i2 = begin2;
     while(*i2 != a9) ++i2;
     ++i2;
     if(i2 == end2) i2 = begin2;
     auto a11 = *i2;
     ++i2;
     if(i2 == end2) i2 = begin2;
     auto a10 = *i2;
     return { a0, a1, a2, a3, v0, v1, a6, a7, v2, a9, a10, a11 };
}

}//namespace happah

