// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <glm/gtx/norm.hpp>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/Plane.h"
#include "happah/geometries/Ray.h"
#include "happah/geometries/Surface.h"
#include "happah/geometries/SurfaceHEZ.h"
#include "happah/geometries/SurfaceSplineBEZ.h"

namespace happah {

template<class Space, hpuint t_degree>
class SurfaceSplineHEZ : public Surface<Space> {
     using Point = typename Space::POINT;

public:
     using ControlPoints = std::vector<Point>;
     using Indices = std::vector<hpuint>;
     using ParameterPoints = std::vector<Point3D>;

     //NOTE: Here we do not use TriangleMesh to represent the parameter space but instead use a direct implementation because a TriangleMesh is a model but a parameter space is not a model.
     SurfaceSplineHEZ(ControlPoints controlPoints, Indices controlPointIndices, ParameterPoints parameterPoints, Indices parameterPointIndices) 
          : m_controlPointIndices{std::move(controlPointIndices)}, m_controlPoints{std::move(controlPoints)}, m_parameterPointIndices{std::move(parameterPointIndices)}, m_parameterPoints{std::move(parameterPoints)} {}

     const Indices& getControlPointIndices() const { return m_controlPointIndices; }

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     hpuint getNumberOfPatches() const {
          static auto nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
          return m_controlPointIndices.size() / nControlPoints;
     }

     const Indices& getParameterPointIndices() const { return m_parameterPointIndices; }

     const ParameterPoints& getParameterPoints() const { return m_parameterPoints; }

     SurfaceHEZ<Space, t_degree> getPatch(hpuint p) const {
          static auto nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
          std::vector<Point> controlPoints;
          auto c = deindex(m_controlPoints, m_controlPointIndices);
          for(auto i = c.cbegin() + p * nControlPoints, end = i + nControlPoints; i != end; ++i) controlPoints.push_back(*i);//TODO: use vector input iterator constructor?
          auto r = deindex(m_parameterPoints, m_parameterPointIndices).cbegin() + 3 * p;
          auto p0 = *r;
          auto p1 = *(++r);
          auto p2 = *(++r);
          return {p0, p1, p2, controlPoints};
     }

     std::tuple<const ControlPoints&, const Indices&> getPatches() const { return std::tie(m_controlPoints, m_controlPointIndices); }

     //NOTE: Replace sth patch with a linear piece whose corners are p0, p1, and p2.  Be aware the shared control points along the edge are also moved.
     void linearize(hpuint s, const Point& p0, const Point& p1, const Point& p2) {
          auto c = m_controlPointIndices.begin() + s * SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
          SurfaceUtilsBEZ::sample(t_degree + 1, [&] (hpreal u, hpreal v, hpreal w) {
               m_controlPoints[*c] = u * p0 + v * p1 + w * p2;
               ++c;
          });
     }

     void reparametrize(const std::vector<hpreal>& factors, const Indices& indices) {
          //TODO: get rid of std::pair in project and just use std::tuple for everything
          //TODO: Similar control points can be shared between different edges and/or patches.  Therefore, we cannot reuse the original control point indices but generate new ones.  We can optimize slightly; any control points shared on an edge between two adjacent patches can be shared in the resulting reparametrization.  In fact, for usability in the graphical interface, these should be the same if they are the same originally.  Thus, we need a data structure that tells us which edges are shared and update the generated control point indices accordingly.
          //TODO: microoptimizations; skip computation on shared edges

          auto n = m_controlPointIndices.size();

          ControlPoints controlPoints;
          controlPoints.reserve(n);

          auto j = m_controlPointIndices.cbegin();
          visit_triangles(factors, indices, [&](hpreal f0, hpreal f1, hpreal f2) {
               auto rowLength = t_degree + 1u;
               while(rowLength > 0u) {
                    hpreal factor = std::pow(f0, rowLength - 1u) * std::pow(f2, t_degree - rowLength + 1u);
                    auto end = j + rowLength;
                    while(j != end) {
                         controlPoints.push_back(m_controlPoints[*j] * factor);
                         factor /= f0;
                         factor *= f1;
                         ++j;
                    }
                    --rowLength;
               }
          });
          m_controlPoints.swap(controlPoints);
          
          m_controlPointIndices.clear();
          m_controlPointIndices.reserve(n);
          for(auto i = 0lu; i < n; ++i) m_controlPointIndices.push_back(i);
     }

     void operator+=(const SurfaceSplineHEZ<Space, t_degree>& surface) {
          //TODO: surface must have the same parameter space
          //TODO: maybe compare indices using simd?
          if(m_controlPointIndices == surface.m_controlPointIndices) for(auto i = m_controlPoints.begin(), j = surface.m_controlPoints.cbegin(), end = m_controlPoints.end(); i != end; ++i, ++j) *i += *j;
          else throw std::runtime_error("+= with different control point indices is not implemented yet.");
     }

     void operator*=(hpreal a) {
          for(auto& controlPoint : m_controlPoints) controlPoint *= a;
     }

     friend SurfaceSplineBEZ<Space, t_degree> restrict(SurfaceSplineHEZ<Space, t_degree> surface, const std::vector<hpreal>& factors, const std::vector<hpuint>& indices) {
          surface.reparametrize(factors, indices);
          return { std::move(surface.m_controlPoints), std::move(surface.m_controlPointIndices) };
     }

private:
     Indices m_controlPointIndices;
     ControlPoints m_controlPoints;
     Indices m_parameterPointIndices;
     ParameterPoints m_parameterPoints;

};//SurfaceSplineHEZ
template<class Space>
using ConstantSurfaceSplineHEZ = SurfaceSplineHEZ<Space, 0>;
template<class Space>
using CubicSurfaceSplineHEZ = SurfaceSplineHEZ<Space, 3>;
template<class Space>
using LinearSurfaceSplineHEZ = SurfaceSplineHEZ<Space, 1>;
template<class Space>
using QuadraticSurfaceSplineHEZ = SurfaceSplineHEZ<Space, 2>;
template<class Space>
using QuarticSurfaceSplineHEZ = SurfaceSplineHEZ<Space, 4>;

template<hpuint degree>
SurfaceSplineHEZ<Space4D, degree> operator*(const SurfaceSplineHEZ<Space1D, degree>& surface, const Point4D& p) {
     std::vector<Point4D> controlPoints;
     controlPoints.reserve(surface.getControlPoints().size());
     for(auto& controlPoint : surface.getControlPoints()) controlPoints.push_back(controlPoint.x * p);
     return SurfaceSplineHEZ<Space4D, degree>(controlPoints, surface.getControlPointIndices(), surface.getParameterPoints(), surface.getParameterPointIndices());
}

}//namespace happah

