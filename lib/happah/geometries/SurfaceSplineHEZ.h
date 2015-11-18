// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

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

     const Indices& getParameterPointIndices() const { return m_parameterPointIndices; }

     const ParameterPoints& getParameterPoints() const { return m_parameterPoints; }

     SurfaceHEZ<Space, t_degree> getPatch(hpuint p) const {
          static auto nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
          std::vector<Point> controlPoints;
          auto c = deindex(m_controlPoints, m_controlPointIndices);
          for(auto i = c.cbegin() + p * nControlPoints, end = i + nControlPoints; i != end; ++i) controlPoints.push_back(*i);
          auto r = deindex(m_parameterPoints, m_parameterPointIndices).cbegin() + 3 * p;
          auto p0 = *r;
          auto p1 = *(++r);
          auto p2 = *(++r);
          return {p0, p1, p2, controlPoints};
     }

     SurfaceSplineBEZ<Space, t_degree> restrict(const Plane& plane, hpreal epsilon = EPSILON) const {
          std::vector<Point2D> parameterPoints;
          parameterPoints.reserve(m_parameterPoints.size());
          std::vector<hpreal> factors;
          factors.reserve(m_parameterPoints.size());

          Ray3D ray{{0.0,0.0,0.0}, {1.0,0.0,0.0}};
          for(const Point3D& parameterPoint : m_parameterPoints) {
               ray.setDirection(glm::normalize(parameterPoint));
               if(auto t = plane.intersect(ray, epsilon)) {
                    auto intersection = ray.getPoint(*t);
                    factors.push_back(std::sqrt(glm::length2(intersection) / glm::length2(parameterPoint)));
                    parameterPoints.push_back(plane.project(intersection));
               } else throw std::runtime_error("Plane does not properly intersect all tetrahedra.");
          }

          ControlPoints controlPoints;
          Indices controlPointIndices;
          std::tie(controlPoints, controlPointIndices) = reparametrize(factors);
          return SurfaceSplineBEZ<Space, t_degree>(std::move(controlPoints), std::move(controlPointIndices), std::move(parameterPoints), m_parameterPointIndices); 
     }

     template<class Visitor>
     void sample(hpuint nSamples, Visitor visit) const {
          auto matrix = SurfaceUtilsBEZ::getEvaluationMatrix<t_degree>(nSamples);
          //TODO; skip multiple computations on common edges and eliminate common points in array
          //TODO: t_degree == 1,2,4, general
          //TODO: optimize calculations below; (u,v,w) do not have to be recalculated each time
          switch(t_degree) {
          case 1:

               break;
          case 2:
 
               break;
          case 3: {
               auto c = m_controlPointIndices.cbegin();
               for(auto i = m_parameterPointIndices.cbegin(), end = m_parameterPointIndices.cend(); i != end; ++i) {
                    auto& p0 = m_parameterPoints[*i];
                    auto& p1 = m_parameterPoints[*(++i)];
                    auto& p2 = m_parameterPoints[*(++i)];

                    auto& c0 = m_controlPoints[*c];
                    auto& c1 = m_controlPoints[*(++c)];
                    auto& c2 = m_controlPoints[*(++c)];
                    auto& c3 = m_controlPoints[*(++c)];
                    auto& c4 = m_controlPoints[*(++c)];
                    auto& c5 = m_controlPoints[*(++c)];
                    auto& c6 = m_controlPoints[*(++c)];
                    auto& c7 = m_controlPoints[*(++c)];
                    auto& c8 = m_controlPoints[*(++c)];
                    auto& c9 = m_controlPoints[*(++c)];
                    ++c;

                    //if(glm::length2(c0) < EPSILON && glm::length2(c1) < EPSILON && glm::length2(c2) < EPSILON && glm::length2(c3) < EPSILON && glm::length2(c4) < EPSILON && glm::length2(c5) < EPSILON && glm::length2(c6) < EPSILON && glm::length2(c7) < EPSILON && glm::length2(c8) < EPSILON && glm::length2(c9) < EPSILON) std::cout << "ERROR: Patch is ZERO!!!\n";
                    
                    //if(glm::length2(c0) < EPSILON && glm::length2(c1) < EPSILON && glm::length2(c2) < EPSILON && glm::length2(c3) < EPSILON) std::cout << "ERROR: Edge 0 is ZERO!!!\n";
                    
                    //if(glm::length2(c3) < EPSILON && glm::length2(c6) < EPSILON && glm::length2(c8) < EPSILON && glm::length2(c9) < EPSILON) std::cout << "ERROR: Edge 1 is ZERO!!!\n";
                    
                    //if(glm::length2(c0) < EPSILON && glm::length2(c4) < EPSILON && glm::length2(c7) < EPSILON && glm::length2(c9) < EPSILON) std::cout << "ERROR: Edge 2 is ZERO!!!\n";
                    
                    auto m = matrix.cbegin();
                    SurfaceUtilsBEZ::sample(nSamples, [&] (hpreal u, hpreal v, hpreal w) {
                         auto sample = *m * c0;
                         sample += *(++m) * c1;
                         sample += *(++m) * c2;
                         sample += *(++m) * c3;
                         sample += *(++m) * c4;
                         sample += *(++m) * c5;
                         sample += *(++m) * c6;
                         sample += *(++m) * c7;
                         sample += *(++m) * c8;
                         sample += *(++m) * c9;
                         ++m;
                         visit(u * p0 + v * p1 + w * p2, sample);
                    });
               }
               break;
          }
          case 4: {
               auto c = m_controlPointIndices.cbegin();
               for(auto i = m_parameterPointIndices.cbegin(), end = m_parameterPointIndices.cend(); i != end; ++i) {
                    auto& p0 = m_parameterPoints[*i];
                    auto& p1 = m_parameterPoints[*(++i)];
                    auto& p2 = m_parameterPoints[*(++i)];

                    auto& c0 = m_controlPoints[*c];
                    auto& c1 = m_controlPoints[*(++c)];
                    auto& c2 = m_controlPoints[*(++c)];
                    auto& c3 = m_controlPoints[*(++c)];
                    auto& c4 = m_controlPoints[*(++c)];
                    auto& c5 = m_controlPoints[*(++c)];
                    auto& c6 = m_controlPoints[*(++c)];
                    auto& c7 = m_controlPoints[*(++c)];
                    auto& c8 = m_controlPoints[*(++c)];
                    auto& c9 = m_controlPoints[*(++c)];
                    auto& c10 = m_controlPoints[*(++c)];
                    auto& c11 = m_controlPoints[*(++c)];
                    auto& c12 = m_controlPoints[*(++c)];
                    auto& c13 = m_controlPoints[*(++c)];
                    auto& c14 = m_controlPoints[*(++c)];
                    ++c;

                    auto m = matrix.cbegin();
                    SurfaceUtilsBEZ::sample(nSamples, [&] (hpreal u, hpreal v, hpreal w) {
                         auto sample = *m * c0;
                         sample += *(++m) * c1;
                         sample += *(++m) * c2;
                         sample += *(++m) * c3;
                         sample += *(++m) * c4;
                         sample += *(++m) * c5;
                         sample += *(++m) * c6;
                         sample += *(++m) * c7;
                         sample += *(++m) * c8;
                         sample += *(++m) * c9;
                         sample += *(++m) * c10;
                         sample += *(++m) * c11;
                         sample += *(++m) * c12;
                         sample += *(++m) * c13;
                         sample += *(++m) * c14;
                         ++m;
                         visit(u * p0 + v * p1 + w * p2, sample);
                    });
               }
               break;
          }
          default:

               break;
          }
     }

     //NOTE: Replace sth piece with a linear piece whose corners are p0, p1, and p2.  Be aware the shared control points along the edge are also moved.
     void linearize(hpuint s, const Point& p0, const Point& p1, const Point& p2) {
          auto c = m_controlPointIndices.begin() + s * SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
          SurfaceUtilsBEZ::sample(t_degree + 1, [&] (hpreal u, hpreal v, hpreal w) {
               m_controlPoints[*c] = u * p0 + v * p1 + w * p2;
               ++c;
          });
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

private:
     Indices m_controlPointIndices;
     ControlPoints m_controlPoints;
     Indices m_parameterPointIndices;
     ParameterPoints m_parameterPoints;

     std::tuple<ControlPoints, Indices> reparametrize(const std::vector<hpreal>& factors) const {
          //TODO: get rid of std::pair in project and just use std::tuple for everything
          //TODO: Similar control points can be shared between different edges and/or patches.  Therefore, we cannot reuse the original control point indices but generate new ones.  We can optimize slightly; any control points shared on an edge between two adjacent patches can be shared in the resulting reparametrization.  In fact, for usability in the graphical interface, these should be the same if they are the same originally.  Thus, we need a data structure that tells us which edges are shared and update the generated control point indices accordingly.
          //TODO: microoptimizations; skip computation on shared edges

          ControlPoints controlPoints;
          controlPoints.reserve(m_controlPointIndices.size());

          auto j = m_controlPointIndices.cbegin();
          for(auto i = m_parameterPointIndices.cbegin(), end = m_parameterPointIndices.cend(); i != end; ++i) {
               auto f0 = factors[*i];
               auto f1 = factors[*(++i)];
               auto f2 = factors[*(++i)];

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
          }
          
          Indices controlPointIndices;
          controlPointIndices.reserve(m_controlPointIndices.size());
          for(auto i = 0lu, end = m_controlPointIndices.size(); i < end; ++i) controlPointIndices.push_back(i);

          return std::make_tuple(std::move(controlPoints), std::move(controlPointIndices));
     }

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

