// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <functional>

#include "happah/Happah.h"
#include "happah/utils/SurfaceUtilsBEZ.h"

namespace happah {

template<class Space, hpuint t_degree>
class SurfaceHEZ : public Geometry3D<Space> {
     using Point = typename Space::POINT;
     using ControlPoints = std::vector<Point>;

public:
     template<class Iterator>
     SurfaceHEZ(Iterator begin)
          : m_controlPoints(begin, begin + SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value) {}

     template<hpuint i0, hpuint i1, hpuint i2>
     const Point& getControlPoint() const { return m_controlPoints[SurfaceUtilsBEZ::get_index<t_degree, i0, i1, i2>::value]; }

     const Point& getControlPoint(hpuint i0, hpuint i1, hpuint i2) const { return m_controlPoints[SurfaceUtilsBEZ::template getIndex<t_degree>(i0, i1, i2)]; }

     std::tuple<const Point&, const Point&, const Point&> getCorners() const { return std::tie(m_controlPoints[0], m_controlPoints[t_degree], m_controlPoints.back()); }

     /*Point getPoint(const Point3D& p) const {
          //TODO: can this be done more efficiently?

          //NOTE: Here we use Cramer's rule.
          hpmat3x3 matA(m_p0, m_p1, m_p2);

          //NOTE: This will be zero if either p0, p1, or p2 is at the origin or if any two or three of them are collinear.  In that case, we have a degenerate parameter tetrahedron.
          hpreal detA = glm::determinant(matA);

          hpmat3x3 matU(p, m_p1, m_p2);
          hpmat3x3 matV(m_p0, p, m_p2);
          hpmat3x3 matW(m_p0, m_p1, p);

          hpreal u = glm::determinant(matU) / detA;
          hpreal v = glm::determinant(matV) / detA;
          hpreal w = glm::determinant(matW) / detA;
          return getPoint(u, v, w);
     }*/

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     template<hpuint i0, hpuint i1, hpuint i2>
     void setControlPoint(const Point& controlPoint) { m_controlPoints[SurfaceUtilsBEZ::get_index<t_degree, i0, i1, i2>::value] = controlPoint; }

     void setControlPoint(hpuint i0, hpuint i1, hpuint i2, const Point& controlPoint) { m_controlPoints[SurfaceUtilsBEZ::template getIndex<t_degree>(i0, i1, i2)] = controlPoint; }

private:
     //NOTE: Order is bn00 bn-110 bn-220 ... bn-101 bn-211 ... bn-202 bn-212 ... b00n.
     ControlPoints m_controlPoints;

};
template<class Space>
using ConstantSurfaceHEZ = SurfaceHEZ<Space, 0>;
template<class Space>
using CubicSurfaceHEZ = SurfaceHEZ<Space, 3>;
template<class Space>
using LinearSurfaceHEZ = SurfaceHEZ<Space, 1>;
template<class Space>
using QuadraticSurfaceHEZ = SurfaceHEZ<Space, 2>;
template<class Space>
using QuarticSurfaceHEZ = SurfaceHEZ<Space, 4>;

template<class Space, hpuint degree>
auto de_casteljau(const SurfaceHEZ<Space, degree>& surface, hpreal u, hpreal v, hpreal w) { return de_casteljau<degree>(std::begin(surface.getControlPoints()), u, v, w); }

}//namespace happah

