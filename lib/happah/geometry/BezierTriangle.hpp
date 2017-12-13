// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/geometry/BezierTriangleMesh.hpp"
#include "happah/geometry/TriangleMesh.hpp"
#include "happah/util/BezierTriangleSubdivider.hpp"

namespace happah {

//DECLARATIONS

template<class Space, hpuint t_degree>
class BezierTriangle;

template<class Space, hpuint degree>
auto de_casteljau(const BezierTriangle<Space, degree>& surface, hpreal u, hpreal v);

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
TriangleMesh<Vertex> make_control_polygon(const BezierTriangle<Space, degree>& surface, VertexFactory&& build = VertexFactory());

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const BezierTriangle<Space, degree>& surface, VertexFactory&& build = VertexFactory());

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const BezierTriangle<Space, degree>& surface, hpuint nSubdivisions, VertexFactory&& build = VertexFactory());

template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> subdivide(const BezierTriangle<Space, degree>& surface, hpuint nSubdivisions);

//DEFINITIONS

template<class Space, hpuint t_degree>
class BezierTriangle {
     using Point = typename Space::POINT;

public:
     BezierTriangle(std::vector<Point> controlPoints)
          : m_controlPoints(std::move(controlPoints)) {}

     template<hpindex i0, hpindex i1, hpindex i2>
     auto& getControlPoint() const { return m_controlPoints[make_control_point_index(t_degree, i0, i1, i2)]; }

     auto& getControlPoint(hpindex i0, hpindex i1, hpindex i2) const { return m_controlPoints[make_control_point_index(t_degree, i0, i1, i2)]; }

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

     auto& getControlPoints() const { return m_controlPoints; }

     hpuint getNumberOfControlPoints() const { return m_controlPoints.size(); }

     template<hpindex i0, hpindex i1, hpindex i2>
     void setControlPoint(Point point) { m_controlPoints[make_control_point_index(t_degree, i0, i1, i2)] = std::move(point); }

     void setControlPoint(hpindex i0, hpindex i1, hpindex i2, Point controlPoint) { m_controlPoints[make_control_point_index(t_degree, i0, i1, i2)] = std::move(point); }

private:
     //NOTE: Order is bn00 bn-110 bn-220 ... bn-101 bn-211 ... bn-202 bn-212 ... b00n.
     std::vector<Point> m_controlPoints;

};//BezierTriangle
template<class Space>
using ConstantBezierTriangle = BezierTriangle<Space, 0>;
template<class Space>
using CubicBezierTriangle = BezierTriangle<Space, 3>;
template<class Space>
using LinearBezierTriangle = BezierTriangle<Space, 1>;
template<class Space>
using QuadraticBezierTriangle = BezierTriangle<Space, 2>;
template<class Space>
using QuarticBezierTriangle = BezierTriangle<Space, 4>;

template<class Space, hpuint degree>
auto de_casteljau(const BezierTriangle<Space, degree>& surface, hpreal u, hpreal v) { return de_casteljau<degree>(std::begin(surface.getControlPoints()), u, v, 1.0 - u - v); }

template<class Space, hpuint degree, class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_control_polygon(const BezierTriangle<Space, degree>& surface, VertexFactory&& build) { return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, std::forward<VertexFactory>(build)); }

template<class Space, hpuint degree, class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const BezierTriangle<Space, degree>& surface, VertexFactory&& build) {
     auto vertices = std::vector<Vertex>();

     vertices.reserve(make_patch_size(degree));
     for(auto& point : surface.getControlPoints()) vertices.push_back(build(point));

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

     return make_triangle_mesh(std::move(vertices), make_indices());
}

template<class Space, hpuint degree, class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const BezierTriangle<Space, degree>& surface, hpuint nSubdivisions, VertexFactory&& build) {
     if(nSubdivisions > 0) return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(subdivide(surface, nSubdivisions), std::forward<VertexFactory>(build));
     else return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, std::forward<VertexFactory>(build));
}

template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> subdivide(const BezierTriangle<Space, degree>& surface, hpuint nSubdivisions) {
     if(nSubdivisions == 0) return { surface.getControlPoints() };
     auto subdivider = BezierTriangleSubdivider<Space, degree>(surface.getControlPoints().begin());
     auto subdivided = subdivider.subdivide(nSubdivisions);
     return { std::move(std::get<0>(subdivided)), std::move(std::get<1>(subdivided)) };
}

}//namespace happah

