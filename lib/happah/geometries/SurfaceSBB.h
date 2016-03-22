// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <functional>

#include "happah/Happah.h"
#include "happah/geometries/Manifold3HEZ.h"
#include "happah/geometries/Surface.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/utils/SurfaceUtilsBEZ.h"

template<class Space, hpuint t_degree>
class SurfaceSBB : public Surface<Space> {
     typedef typename Space::POINT Point;

public:
     typedef SurfaceUtilsBEZ::ControlPoints<Space> ControlPoints;

     static const hpuint NUMBER_OF_CONTROL_POINTS = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
     static const hpuint NUMBER_OF_CONTROL_POLYGON_TRIANGLES = SurfaceUtilsBEZ::get_number_of_control_polygon_triangles<t_degree>::value;

     SurfaceSBB(const Manifold3HEZ<Space, t_degree>& manifold)
          : m_manifold(manifold) {
          auto ps = m_manifold.getParameterTriangle();
          m_p0 = Sphere::Utils::getAbscissa(std::get<0>(ps));
          m_p1 = Sphere::Utils::getAbscissa(std::get<1>(ps));
          m_p2 = Sphere::Utils::getAbscissa(std::get<2>(ps));
     }
     //TODO: getparameterVertex in surfaces
     //NOTE: Here we use spherical coordinates.
     //NOTE: There exist two triangles on the sphere associated with the given vertices.
     //TODO: rethink last note and behavior at ends of parameter space
     SurfaceSBB(const Point2D& p0, const Point2D& p1, const Point2D& p2)
          : SurfaceSBB(p0, p1, p2, ControlPoints()) {}

     SurfaceSBB(const Point2D& p0, const Point2D& p1, const Point2D& p2, const ControlPoints& controlPoints)
          : m_manifold(Sphere::Utils::getPoint(p0), Sphere::Utils::getPoint(p1), Sphere::Utils::getPoint(p2), controlPoints), m_p0(p0), m_p1(p1), m_p2(p2) {}

     template<hpuint t_i0, hpuint t_i1, hpuint t_i2>
     const Point& getControlPoint() const { return m_manifold.getControlPoint<t_i0, t_i1, t_i2>(); }
     const Point& getControlPoint(hpuint i0, hpuint i1, hpuint i2) const { return m_manifold.getControlPoint(i0, i1, i2); }
     std::tuple<const Point2D&, const Point2D&, const Point2D&> getParameterTriangle() const { return std::make_tuple(std::cref(m_p0), std::cref(m_p1), std::cref(m_p2)); }
     Point getPoint(const Point2D& p) const { return getPoint(p.x, p.y); }
     Point getPoint(hpreal u, hpreal v) const { return m_manifold.getPoint(Sphere::Utils::getPoint(u, v)); }
     template<hpuint t_i0, hpuint t_i1, hpuint t_i2>
     void setControlPoint(const Point& controlPoint) { m_manifold.setControlPoint<t_i0, t_i1, t_i2>(controlPoint); }
     void setControlPoint(hpuint i0, hpuint i1, hpuint i2, const Point& controlPoint) { m_manifold.setControlPoint(i0, i1, i2, controlPoint); }
     void setControlPoints(const ControlPoints& controlPoints) { m_manifold.setControlPoints(controlPoints); }
     //NOTE: nSamples is the number of samples on one edge of the parameter triangle.
     template<class Vertex = VertexP<Space>, typename = typename std::enable_if<(t_degree > 0)>::type>
     TriangleMesh<Vertex>* toTriangleMesh(hpuint nSamples = 100) const {
          hpuint pseudodegree = nSamples - 1;
          auto pt = m_manifold.getParameterTriangle();
          auto factory = [&](hpreal u, hpreal v, hpreal w) {
               Point3D abscissa = u * std::get<0>(pt) + v * std::get<1>(pt) + w * std::get<2>(pt);
               return build_vertex<Vertex>::exec(Sphere::Utils::getAbscissa(abscissa), m_manifold.getPoint(abscissa));
          };
          std::vector<Vertex>* vertices = SurfaceUtilsBEZ::sample<Vertex>(pseudodegree, factory);
          std::vector<hpuint>* indices = new std::vector<hpuint>();
          indices->reserve(3 * SurfaceUtilsBEZ::getNumberOfControlPolygonTriangles(pseudodegree));

          //TODO: move to SurfaceUtilsBEZ
          //TODO: arrange vertices in counterclockwise order
          //TODO: arrange provoking vertices correctly, if not too difficult
          hpuint index = 0, nPointsInRow = nSamples;
          for(hpuint row = 0; row < (nSamples - 1); ++row) {
               for(hpuint i = 0; i < (2 * nPointsInRow - 3); ++i) {
                    if(i % 2 == 0) {
                         indices->push_back(index);
                         indices->push_back(++index);
                         indices->push_back(index + nPointsInRow - 1);
                    } else {
                         indices->push_back(index);
                         indices->push_back(index + nPointsInRow);
                         indices->push_back(index + nPointsInRow - 1);
                    }
               }
               --nPointsInRow;
               ++index;
          }

          return new TriangleMesh<Vertex>(vertices, indices);
     }

private:
     Manifold3HEZ<Space, t_degree> m_manifold;
     Point2D m_p0, m_p1, m_p2;

     template<class Vertex, typename = void>
     struct build_vertex;

     template<class Vertex>
     struct build_vertex<Vertex, typename enable_if_absolute_vertex<Vertex>::type> {
          static Vertex exec(const Point2D& abscissa, const Point1D& ordinate) { return VertexFactory<Vertex>::build(Sphere::Utils::getPoint(abscissa, ordinate)); }
     };

     template<class Vertex>
     struct build_vertex<Vertex, typename enable_if_relative_vertex<Vertex>::type> {
          static Vertex exec(const Point2D& abscissa, const Point1D& ordinate) { return VertexFactory<Vertex>::build(abscissa, ordinate); }
     };

};
template<class Space>
using ConstantSurfaceSBB = SurfaceSBB<Space, 0>;
template<class Space>
using CubicSurfaceSBB = SurfaceSBB<Space, 3>;
template<class Space>
using LinearSurfaceSBB = SurfaceSBB<Space, 1>;
template<class Space>
using QuadraticSurfaceSBB = SurfaceSBB<Space, 2>;
template<class Space>
using QuarticSurfaceSBB = SurfaceSBB<Space, 4>;

