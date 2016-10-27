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
     using ControlPoints = std::vector<Point>;

public:
     SurfaceSplineHEZ(ControlPoints controlPoints, Indices indices) 
          : m_indices{std::move(indices)}, m_controlPoints{std::move(controlPoints)}  {}

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     hpuint getNumberOfPatches() const { return m_indices.size() / SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value; }

     SurfaceHEZ<Space, t_degree> getPatch(hpuint p) const {
          static auto nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
          return { deindex(m_controlPoints, m_indices).begin() + p * nControlPoints };
     }

     std::tuple<const ControlPoints&, const Indices&> getPatches() const { return std::tie(m_controlPoints, m_indices); }

     //NOTE: Replace sth patch with a linear piece whose corners are p0, p1, and p2.  Be aware the shared control points along the edge are also moved.
     void linearize(hpuint s, const Point& p0, const Point& p1, const Point& p2) {
          static auto nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
          auto c = m_indices.begin() + s * nControlPoints;
          SurfaceUtilsBEZ::sample(t_degree + 1, [&] (hpreal u, hpreal v, hpreal w) {
               m_controlPoints[*c] = u * p0 + v * p1 + w * p2;
               ++c;
          });
     }

     void reparametrize(const std::vector<hpreal>& factors, const Indices& indices) {
          //TODO: get rid of std::pair in project and just use std::tuple for everything
          //TODO: Similar control points can be shared between different edges and/or patches.  Therefore, we cannot reuse the original control point indices but generate new ones.  We can optimize slightly; any control points shared on an edge between two adjacent patches can be shared in the resulting reparametrization.  In fact, for usability in the graphical interface, these should be the same if they are the same originally.  Thus, we need a data structure that tells us which edges are shared and update the generated control point indices accordingly.
          //TODO: microoptimizations; skip computation on shared edges

          auto n = m_indices.size();

          ControlPoints controlPoints;
          controlPoints.reserve(n);

          auto j = m_indices.cbegin();
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
          
          m_indices.clear();
          m_indices.reserve(n);
          for(auto i = 0lu; i < n; ++i) m_indices.push_back(i);
     }

     void operator+=(const SurfaceSplineHEZ<Space, t_degree>& surface) {
          //TODO: surface must have the same parameter space
          //TODO: maybe compare indices using simd?
          if(m_indices == surface.m_indices) for(auto i = m_controlPoints.begin(), j = surface.m_controlPoints.cbegin(), end = m_controlPoints.end(); i != end; ++i, ++j) *i += *j;
          else throw std::runtime_error("+= with different control point indices is not implemented yet.");
     }

     void operator*=(hpreal a) { for(auto& controlPoint : m_controlPoints) controlPoint *= a; }

     friend SurfaceSplineBEZ<Space, t_degree> restrict(SurfaceSplineHEZ<Space, t_degree> surface, const std::vector<hpreal>& factors, const std::vector<hpuint>& indices) {
          surface.reparametrize(factors, indices);
          return { std::move(surface.m_controlPoints), std::move(surface.m_indices) };
     }

private:
     Indices m_indices;
     ControlPoints m_controlPoints;

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
     return { controlPoints, std::get<1>(surface.getPatches()) };
}

template<class Space, hpuint degree, class Visitor>
void visit_patches(const SurfaceSplineHEZ<Space, degree>& surface, Visitor&& visit) {//TODO: visit_patches on object with getPatches method only?
     auto patches = surface.getPatches();
     auto temp = deindex(std::get<0>(patches), std::get<1>(patches));
     visit_patches<degree>(temp.begin(), temp.end(), std::forward<Visitor>(visit));
}

}//namespace happah

