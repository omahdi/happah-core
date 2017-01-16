// Copyright 2015 - 2016
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
     using ControlPoints = std::vector<Point>;

public:
     SurfaceBEZ(ControlPoints controlPoints)
          : m_controlPoints(std::move(controlPoints)) {}

     template<hpuint i0, hpuint i1, hpuint i2>
     const Point& getControlPoint() const { return m_controlPoints[make_offset(t_degree, i0, i1, i2)]; }

     const Point& getControlPoint(hpuint i0, hpuint i1, hpuint i2) const { return m_controlPoints[make_offset(t_degree, i0, i1, i2)]; }

     /* m_a(m_p1.y - m_p2.y), m_b(m_p0.x - m_p2.x), m_c(m_p2.x - m_p1.x), m_d(m_p2.y - m_p0.y)
     hpreal e = m_a * m_b - m_c * m_d;
     m_a /= e;
     m_b /= e;
     m_c /= e;
     m_d /= e;
     Point getPoint(const Point2D& p) const {
          hpreal tx = p.x - m_p2.x;
          hpreal ty = p.y - m_p2.y;
          hpreal u = m_a * tx + m_c * ty;
          hpreal v = m_d * tx + m_b * ty;
          return getPoint(u, v);
     }*/

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     hpuint getNumberOfControlPoints() const { return m_controlPoints.size(); }

     template<hpuint i0, hpuint i1, hpuint i2>
     void setControlPoint(Point controlPoint) { m_controlPoints[make_offset(t_degree, i0, i1, i2)] = std::move(controlPoint); }

     void setControlPoint(hpuint i0, hpuint i1, hpuint i2, Point controlPoint) { m_controlPoints[make_offset(t_degree, i0, i1, i2)] = std::move(controlPoint); }

private:
     //NOTE: Order is bn00 bn-110 bn-220 ... bn-101 bn-211 ... bn-202 bn-212 ... b00n.
     ControlPoints m_controlPoints;

     /*template<class Vertex>
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

     };//TriangleMeshBuilder*/

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

template<class Space, hpuint degree>
auto de_casteljau(const SurfaceBEZ<Space, degree>& surface, hpreal u, hpreal v) { return de_casteljau<degree>(std::begin(surface.getControlPoints()), u, v, 1.0 - u - v); }

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceBEZ<Space, degree>& surface, hpuint nSubdivisions) {
     if(nSubdivisions == 0) return { surface.getControlPoints() };
     auto subdivider = SurfaceSubdividerBEZ<Space, degree>(surface.getControlPoints().begin());
     auto subdivided = subdivider.subdivide(nSubdivisions);
     return { std::move(std::get<0>(subdivided)), std::move(std::get<1>(subdivided)) };
}

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory()) {
     using Point = typename Space::POINT;
     static_assert(std::is_base_of<Vertex, decltype(factory(Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");

     std::vector<Vertex> vertices;
     vertices.reserve(make_patch_size(degree));
     for(auto& controlPoint : surface.getControlPoints()) vertices.push_back(factory(controlPoint));

     auto make_indices = [&]() -> auto {
          switch(degree) {
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
                    Indices indices;
                    indices.reserve(3 * make_control_polygon_size(degree));

                    auto index = indices.begin();
                    hpuint i = 0;
                    hpuint ai = i+degree;
                    while(i < degree) {
                         *index = i;//TODO: shouldn't this be push_back
                         ++index;
                         *index = ++i;
                         ++index;
                         *index = ++ai;
                         ++index;
                    }
                    ++i;
                    hpuint d = degree;
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
     };

     auto indices = make_indices();

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceBEZ<Space, degree>& surface, hpuint nSubdivisions, VertexFactory&& factory = VertexFactory()) {
     if(nSubdivisions > 0) return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(subdivide(surface, nSubdivisions), std::forward<VertexFactory>(factory));
     else return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, std::forward<VertexFactory>(factory));
}

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
TriangleMesh<Vertex> make_control_polygon(const SurfaceBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory()) { return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, std::forward<VertexFactory>(factory)); }

}//namespace happah

