// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.h"
#include "happah/geometries/Surface.h"
#include "happah/geometries/SurfaceSplineBEZ.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/utils/SurfaceSubdividerBEZ.h"
#include "happah/utils/SurfaceUtilsBEZ.h"

namespace happah {

template<class Space, hpuint t_degree>
class SurfaceBEZ : public Surface<Space> {
     using Point = typename Space::POINT;

public:
     using ControlPoints = typename SurfaceUtilsBEZ::ControlPoints<Space>;

     SurfaceBEZ(Point2D p0, Point2D p1, Point2D p2)
          : SurfaceBEZ(getDefaultControlPoints(), std::move(p0), std::move(p1), std::move(p2)) {}

     SurfaceBEZ(ControlPoints controlPoints, Point2D p0, Point2D p1, Point2D p2)
          : m_controlPoints(std::move(controlPoints)), m_p0(std::move(p0)), m_p1(std::move(p1)), m_p2(std::move(p2)), m_a(m_p1.y - m_p2.y), m_b(m_p0.x - m_p2.x), m_c(m_p2.x - m_p1.x), m_d(m_p2.y - m_p0.y) {
          hpreal e = m_a * m_b - m_c * m_d;
          m_a /= e;
          m_b /= e;
          m_c /= e;
          m_d /= e;
     }

     template<hpuint t_i0, hpuint t_i1, hpuint t_i2>
     const Point& getControlPoint() const { return m_controlPoints[SurfaceUtilsBEZ::get_index<t_degree, t_i0, t_i1, t_i2>::value]; }

     const Point& getControlPoint(hpuint i0, hpuint i1, hpuint i2) const { return m_controlPoints[SurfaceUtilsBEZ::template getIndex<t_degree>(i0, i1, i2)]; }

     template<class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, bool t_includeParameterPoints = is_relative_vertex<Vertex>::value>
     TriangleMesh<Vertex> getControlPolygon(VertexFactory&& factory = VertexFactory()) const { return toTriangleMesh<Vertex, VertexFactory, t_includeParameterPoints>(0, std::forward<VertexFactory>(factory)); }

     std::tuple<const Point2D&, const Point2D&, const Point2D&> getParameterPoints() const { return std::make_tuple(std::cref(m_p0), std::cref(m_p1), std::cref(m_p2)); }

     Point getPoint(const Point2D& p) const {
          hpreal tx = p.x - m_p2.x;
          hpreal ty = p.y - m_p2.y;
          hpreal u = m_a * tx + m_c * ty;
          hpreal v = m_d * tx + m_b * ty;
          return getPoint(u, v);
     }

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     Point getPoint(hpreal u, hpreal v) const { return SurfaceUtilsBEZ::template evaluate<Space, t_degree>(u, v, 1.0 - u - v, m_controlPoints); }

     template<hpuint t_i0, hpuint t_i1, hpuint t_i2>
     void setControlPoint(Point controlPoint) { m_controlPoints[SurfaceUtilsBEZ::get_index<t_degree, t_i0, t_i1, t_i2>::value] = std::move(controlPoint); }

     void setControlPoint(hpuint i0, hpuint i1, hpuint i2, Point controlPoint) { m_controlPoints[SurfaceUtilsBEZ::template getIndex<t_degree>(i0, i1, i2)] = std::move(controlPoint); }

     void setControlPoints(ControlPoints controlPoints) { m_controlPoints = std::move(controlPoints); }

     template<bool t_includeParameterPoints = false>
     SurfaceSplineBEZ<Space, t_degree> subdivide(hpuint nSubdivisions) const {
          using Subdivider = SurfaceSubdividerBEZ<Space, t_degree>;
          Subdivider subdivider(m_controlPoints.cbegin());
          auto subdivided = subdivider.subdivide(nSubdivisions);
          if(t_includeParameterPoints) {
               auto temp = Subdivider::getParameterPoints(m_p0, m_p1, m_p2, nSubdivisions);
               return SurfaceSplineBEZ<Space, t_degree>(subdivided.first, subdivided.second, temp.first, temp.second);
          } else return SurfaceSplineBEZ<Space, t_degree>(subdivided.first, subdivided.second);
     }

     template<class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, bool t_includeParameterPoints = is_relative_vertex<Vertex>::value, typename = typename std::enable_if<(t_degree > 0)>::type>
     TriangleMesh<Vertex> toTriangleMesh(hpuint nSubdivisions = 4, VertexFactory&& factory = VertexFactory()) const {
          if(nSubdivisions > 0) {
               auto subdivided = subdivide<t_includeParameterPoints>(nSubdivisions);
               return subdivided.getControlPolygon<Vertex, VertexFactory, t_includeParameterPoints>(std::forward<VertexFactory>(factory));
          }
          TriangleMeshBuilder<Vertex, t_includeParameterPoints> builder(*this);
          return builder.build(std::forward<VertexFactory>(factory));
     }

private:
     //NOTE: Order is bn00 bn-110 bn-220 ... bn-101 bn-211 ... bn-202 bn-212 ... b00n.
     ControlPoints m_controlPoints;
     Point2D m_p0, m_p1, m_p2;//parameter triangle
     hpreal m_a, m_b, m_c, m_d;//cached calculations

     static ControlPoints getDefaultControlPoints() {
          ControlPoints points;
          points.resize(SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value);
          return std::move(points);
     }

     template<class Vertex, bool t_includeParameterPoints>
     struct TriangleMeshBuilder;

     template<class Vertex>
     struct TriangleMeshBuilder<Vertex, false> {
          
          TriangleMeshBuilder(const SurfaceBEZ& surface) : m_surface(surface) {}

          template<class VertexFactory>
          TriangleMesh<Vertex> build(VertexFactory&& factory) const {
               static_assert(std::is_base_of<Vertex, decltype(factory(Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");

               std::vector<Vertex> vertices;
               vertices.reserve(m_surface.m_controlPoints.size());

               for(const Point& controlPoint : m_surface.m_controlPoints)
                    vertices.push_back(factory(controlPoint));

               std::vector<hpuint> indices = SurfaceUtilsBEZ::template buildTriangleMeshIndices<t_degree>();

               return TriangleMesh<Vertex>(std::move(vertices), std::move(indices));
          }

     private:
          const SurfaceBEZ& m_surface;

     };//TriangleMeshBuilder

     template<class Vertex>
     struct TriangleMeshBuilder<Vertex, true> {

          TriangleMeshBuilder(const SurfaceBEZ<Space, t_degree>& surface) : m_surface(surface) {}

          template<class VertexFactory>
          TriangleMesh<Vertex> build(VertexFactory&& factory) const {
               static_assert(std::is_base_of<Vertex, decltype(factory(Point2D(0.0), Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");
               
               std::vector<Vertex> vertices;

               switch(t_degree) {
               case 1:
                    vertices = {
                         build<1, 0, 0>(factory),
                         build<0, 1, 0>(factory),
                         build<0, 0, 1>(factory)
                    };
                    break;
               case 2:
                    vertices = {
                         build<2, 1, 0>(factory),
                         build<1, 1, 0>(factory),
                         build<0, 2, 0>(factory),
                         build<1, 0, 1>(factory),
                         build<0, 1, 1>(factory),
                         build<0, 0, 2>(factory)
                    };
                    break;
               case 3:
                    vertices = {
                         build<3, 0, 0>(factory),
                         build<2, 1, 0>(factory),
                         build<1, 2, 0>(factory),
                         build<0, 3, 0>(factory),
                         build<2, 0, 1>(factory),
                         build<1, 1, 1>(factory),
                         build<0, 2, 1>(factory),
                         build<1, 0, 2>(factory),
                         build<0, 1, 2>(factory),
                         build<0, 0, 3>(factory)
                    };
                    break;
               case 4:
                    vertices = {
                         build<4, 0, 0>(factory),
                         build<3, 1, 0>(factory),
                         build<2, 2, 0>(factory),
                         build<1, 3, 0>(factory),
                         build<0, 4, 0>(factory),
                         build<3, 0, 1>(factory),
                         build<2, 1, 1>(factory),
                         build<1, 2, 1>(factory),
                         build<0, 3, 1>(factory),
                         build<2, 0, 2>(factory),
                         build<1, 1, 2>(factory),
                         build<0, 2, 2>(factory),
                         build<1, 0, 3>(factory),
                         build<0, 1, 3>(factory),
                         build<0, 0, 4>(factory)
                    };
                    break;
               case 5:
                    vertices = {
                         build<5, 0, 0>(factory),
                         build<4, 1, 0>(factory),
                         build<3, 2, 0>(factory),
                         build<2, 3, 0>(factory),
                         build<1, 4, 0>(factory),
                         build<0, 5, 0>(factory),
                         build<4, 0, 1>(factory),
                         build<3, 1, 1>(factory),
                         build<2, 2, 1>(factory),
                         build<1, 3, 1>(factory),
                         build<0, 4, 1>(factory),
                         build<3, 0, 2>(factory),
                         build<2, 1, 2>(factory),
                         build<1, 2, 2>(factory),
                         build<0, 3, 2>(factory),
                         build<2, 0, 3>(factory),
                         build<1, 1, 3>(factory),
                         build<0, 2, 3>(factory),
                         build<1, 0, 4>(factory),
                         build<0, 1, 4>(factory),
                         build<0, 0, 5>(factory)
                    };
                    break;
               default: std::cerr << "ERROR: Not implemented yet.\n";//TODO
               }

               std::vector<hpuint> indices = SurfaceUtilsBEZ::template buildTriangleMeshIndices<t_degree>();

               return TriangleMesh<Vertex>({new std::vector<Vertex>(std::move(vertices))}, {new std::vector<hpuint>(std::move(indices))});
          }

     private:
          const SurfaceBEZ& m_surface;

          template<hpuint t_i0, hpuint t_i1, hpuint t_i2, class VertexFactory>
          Vertex build(VertexFactory& factory) const { return factory(((hpreal)t_i0 / (hpreal)t_degree) * m_surface.m_p0 + ((hpreal)t_i1 / (hpreal)t_degree) * m_surface.m_p1 + ((hpreal) t_i2 / (hpreal)t_degree) * m_surface.m_p2, m_surface.m_controlPoints[SurfaceUtilsBEZ::get_index<t_degree, t_i0, t_i1, t_i2>::value]); }

     };//TriangleMeshBuilder

};//SurfaceBEZ
template<class Space>
using ConstantSurfaceBEZ = SurfaceBEZ<Space, 0>;
template<class Space>
using CubicSurfaceBEZ = SurfaceBEZ<Space, 3>;
template<class Space>
using LinearSurfaceBEZ = SurfaceBEZ<Space, 1>;
template<class Space>
using QuadraticSurfaceBEZ = SurfaceBEZ<Space, 2>;
template<class Space>
using QuarticSurfaceBEZ = SurfaceBEZ<Space, 4>;

}//namespace happah

