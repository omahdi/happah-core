// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <type_traits>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/Surface.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/utils/DeindexedArray.h"
#include "happah/utils/SurfaceSubdividerBEZ.h"
#include "happah/utils/SurfaceSplineUtilsBEZ.h"
#include "happah/utils/SurfaceUtilsBEZ.h"
#include "happah/utils/VertexFactory.h"

namespace happah {

template<class Space, hpuint t_degree>
class SurfaceSplineBEZ : public Surface<Space> {
     using Point = typename Space::POINT;

public:
     using ControlPoints = std::vector<Point>;
     using Indices = std::vector<hpuint>;
     using Mode = SurfaceSplineUtilsBEZ::Mode;
     using ParameterPoints = std::vector<Point2D>;

     SurfaceSplineBEZ() {}

     //NOTE: There is no automatic way of figuring out the neighborhood of patches given only the control points.
     SurfaceSplineBEZ(ControlPoints controlPoints, Indices controlPointIndices)
          : m_controlPointIndices{std::move(controlPointIndices)}, m_controlPoints{std::move(controlPoints)} {}

     SurfaceSplineBEZ(ControlPoints controlPoints, Indices controlPointIndices, ParameterPoints parameterPoints, Indices parameterPointIndices)
          : m_controlPointIndices{std::move(controlPointIndices)}, m_controlPoints{std::move(controlPoints)}, m_parameterPointIndices{std::move(parameterPointIndices)}, m_parameterPoints{std::move(parameterPoints)} {}

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     const ParameterPoints& getParameterPoints() const { return m_parameterPoints; }

     std::tuple<const ParameterPoints&, const Indices&> getParameterSpace() const { return std::tie(m_parameterPoints, m_parameterPointIndices); }

     std::tuple<const ControlPoints&, const Indices&> getPatches() const { return std::tie(m_controlPoints, m_controlPointIndices); }

     hpuint getNumberOfPatches() const { return m_controlPointIndices.size() / SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value; }

private:
     Indices m_controlPointIndices;
     ControlPoints m_controlPoints;
     Indices m_parameterPointIndices;
     ParameterPoints m_parameterPoints;

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const SurfaceSplineBEZ<Space, t_degree>& surface) {
          stream << surface.m_controlPointIndices << '\n';
          stream << surface.m_controlPoints << '\n';
          stream << surface.m_parameterPointIndices << '\n';
          stream << surface.m_parameterPoints;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, SurfaceSplineBEZ<Space, t_degree>& surface) {
          stream >> surface.m_controlPointIndices;
          stream >> surface.m_controlPoints;
          stream >> surface.m_parameterPointIndices;
          stream >> surface.m_parameterPoints;
          return stream;
     }

};//SurfaceSplineBEZ
template<class Space>
using ConstantSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 0>;
template<class Space>
using CubicSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 3>;
template<class Space>
using LinearSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 1>;
template<class Space>
using SexticSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 6>;
template<class Space>
using QuadraticSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 2>;
template<class Space>
using QuarticSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 4>;
template<class Space>
using QuinticSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 5>;

template<class Space, hpuint t_degree, class Visitor>
void visit_patches(const SurfaceSplineBEZ<Space, t_degree>& surface, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
     auto patches = surface.getPatches();
     auto temp = deindex(std::get<0>(patches), std::get<1>(patches));
     for(auto i = temp.begin(), end = temp.end(); i != end; i += nControlPoints) visit(i);
}

template<class Space, hpuint t_degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, t_degree>& surface, Visitor&& visit) {}

template<class Space, hpuint degree, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSamples, Visitor&& visit) {
     //TODO; skip multiple computations on common edges and eliminate common points in array
     //TODO: t_degree == 1,2,4, general
     //TODO: optimize calculations below; (u,v,w) do not have to be recalculated each time
     auto matrix = SurfaceUtilsBEZ::getEvaluationMatrix<degree>(nSamples);
     auto patches = surface.getPatches();
     auto c = deindex(std::get<0>(patches), std::get<1>(patches)).begin();
     auto domain = surface.getParameterSpace();
     visit_triangles(std::get<0>(domain), std::get<1>(domain), [&](const Point2D& p0, const Point2D& p1, const Point2D& p2) {
          auto m = matrix.cbegin();
          switch(degree) {
          case 1:

               break;
          case 2:

               break;
          case 3: {
               auto& c0 = *c;
               auto& c1 = *(++c);
               auto& c2 = *(++c);
               auto& c3 = *(++c);
               auto& c4 = *(++c);
               auto& c5 = *(++c);
               auto& c6 = *(++c);
               auto& c7 = *(++c);
               auto& c8 = *(++c);
               auto& c9 = *(++c);

               auto m = matrix.cbegin();
               SurfaceUtilsBEZ::sample(nSamples, [&](hpreal u, hpreal v, hpreal w) {
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
               break;
          }
          case 4: {//TODO: implementation in hez is exactly the same; refactor to one location
               auto& c0 = *c;
               auto& c1 = *(++c);
               auto& c2 = *(++c);
               auto& c3 = *(++c);
               auto& c4 = *(++c);
               auto& c5 = *(++c);
               auto& c6 = *(++c);
               auto& c7 = *(++c);
               auto& c8 = *(++c);
               auto& c9 = *(++c);
               auto& c10 = *(++c);
               auto& c11 = *(++c);
               auto& c12 = *(++c);
               auto& c13 = *(++c);
               auto& c14 = *(++c);

               SurfaceUtilsBEZ::sample(nSamples, [&](hpreal u, hpreal v, hpreal w) {
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
               break;
          }
          default:

               break;
          }
          ++c;
     });
}

//**********************************************************************************************************************************
//TODO: cleanup below
//**********************************************************************************************************************************

template<class Space, hpuint t_degree>
class SurfaceSplineSubdividerBEZ {
     static_assert(t_degree > 1, "Surface subdivision only makes sense for degree greater than one because a surface of degree one is planar.");
     using Point = typename Space::POINT;

public:
     using ControlPoints = std::vector<Point>;
     using Indices = std::vector<hpuint>;

     template<class Point>
     static std::pair<std::vector<Point>, std::vector<hpuint> > getParameterPoints(const std::vector<Point>& parameterPoints, const std::vector<hpuint>& parameterPointIndices, hpuint nSubdivisions) {
          std::vector<Point> points;
          std::vector<hpuint> indices;

          //points.reserve();//TODO
          //indices.reserve();

          //TODO: we can significantly reduce the number of parameter points by taking common points into consideration
          for(auto i = parameterPointIndices.cbegin(), end = parameterPointIndices.cend(); i != end; ++i) {
               const Point& p0 = parameterPoints[*i];
               const Point& p1 = parameterPoints[*(++i)];
               const Point& p2 = parameterPoints[*(++i)];

               auto temp = SurfaceSubdividerBEZ<Space, t_degree>::getParameterPoints(p0, p1, p2, nSubdivisions);
               hpuint offset = points.size();
               for(auto& j : temp.second) j += offset;
               std::move(temp.first.begin(), temp.first.end(), std::back_inserter(points));
               std::move(temp.second.begin(), temp.second.end(), std::back_inserter(indices));
          }

          return std::make_pair(std::move(points), std::move(indices));
     }

     /*
      * NOTE: Each patch is a tuple of seven items:
      *   1. An iterator pointing to the control points of the patch.  The control points are ordered left to right starting at the bottom row and ending at the top row, which contains one point.
      *   2. Three unsigned integers specifying the three neighboring patches or UNULL if there is no neighboring patch.
      *   3. Three booleans specifying whether the edge with the corresponding neighbors is shared, that is, the control points are the same.  If that is the case, the algorithm should output a subdivided spline that shares the control points on the shared edges as well. (TODO: implement this; the entire spline has to be subdivided at once to avoid multiple computations on common edges; the data structure and algorithm should be similar to the one in the surfacesubdividerbez class but taking into account shared edges)
      */
     SurfaceSplineSubdividerBEZ(const SurfaceSplineBEZ<Space, t_degree>& surface) 
          : m_nPatches(surface.getNumberOfPatches()) {
          m_subdividers.reserve(m_nPatches);
          auto push_patch = [&](auto begin) { m_subdividers.emplace_back(begin); };
          visit_patches(surface, push_patch);
     }

     std::pair<ControlPoints, Indices> subdivide(hpuint nSubdivisions) {
          assert(nSubdivisions > 0);

          ControlPoints points;
          Indices indices;

          points.reserve(m_nControlPoints * m_nPatches);
          indices.reserve(3 * m_nTriangles * m_nPatches);

          for(auto& subdivider : m_subdividers) {
               auto subdivided = subdivider.subdivide(nSubdivisions);
               hpuint offset = points.size();
               for(auto& i : subdivided.second) i += offset;
               std::move(subdivided.first.begin(), subdivided.first.end(), std::back_inserter(points));
               std::move(subdivided.second.begin(), subdivided.second.end(), std::back_inserter(indices));
          }

          return { std::move(points), std::move(indices) };
     }

private:
     static constexpr hpuint m_nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
     hpuint m_nPatches;
     static constexpr hpuint m_nTriangles = SurfaceUtilsBEZ::get_number_of_control_polygon_triangles<t_degree>::value;
     std::vector<SurfaceSubdividerBEZ<Space, t_degree> > m_subdividers;

};//SurfaceSplineSubdividerBEZ
template<class Space>
using CubicSurfaceSplineSubdividerBEZ = SurfaceSplineSubdividerBEZ<Space, 3>;
template<class Space>
using QuadraticSurfaceSplineSubdividerBEZ = SurfaceSplineSubdividerBEZ<Space, 2>;
template<class Space>
using QuarticSurfaceSplineSubdividerBEZ = SurfaceSplineSubdividerBEZ<Space, 4>;

template<class Space, hpuint degree, bool t_includeParameterPoints = false>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions) {
     //assert(!t_includeParameterPoints || (t_includeParameterPoints && m_parameterPoints.size() > 0));
     using Subdivider = SurfaceSplineSubdividerBEZ<Space, degree>;
     Subdivider subdivider(surface);
     auto subdivided = subdivider.subdivide(nSubdivisions);
     //if(t_includeParameterPoints) {//TODO
     //     auto temp = Subdivider::getParameterPoints(m_parameterPoints, m_parameterPointIndices, nSubdivisions);
     //     return SurfaceSplineBEZ(subdivided.first, subdivided.second, temp.first, temp.second);
     //} else 
     return { std::move(subdivided.first), std::move(subdivided.second) };
}

template<class Vertex, bool t_includeParameterPoints>
struct TriangleMeshBuilder;

template<class Vertex>
struct TriangleMeshBuilder<Vertex, false> {
     
     template<class Space, hpuint degree, class VertexFactory>
     static TriangleMesh<Vertex> build(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory) {
          using Point = typename Space::POINT;
          static_assert(std::is_base_of<Vertex, decltype(factory(Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");

          std::vector<Vertex> vertices;
          vertices.reserve(surface.getControlPoints().size());
          for(auto& point : surface.getControlPoints()) vertices.push_back(factory(point));

          auto indices = SurfaceSplineUtilsBEZ::template buildTriangleMeshIndices<degree>(std::get<1>(surface.getPatches()));

          return make_triangle_mesh(std::move(vertices), std::move(indices));
     }

};//TriangleMeshBuilder

/*template<class Vertex>
struct TriangleMeshBuilder<Vertex, true> {

     TriangleMeshBuilder(const SurfaceSplineBEZ& surface) : m_surface(surface) {}

     template<class VertexFactory>
     TriangleMesh<Vertex> build(VertexFactory&& factory) const {
          static_assert(std::is_base_of<Vertex, decltype(factory(Point2D(0.0), Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");
          
          assert(m_surface.m_parameterPointIndices.size() > 0);
          assert(m_surface.m_parameterPoints.size() > 0);

          //TODO: Similar control points can be shared between different edges and/or patches.  Therefore, when including the parameter points, we cannot reuse the original control point indices but generate new ones.  Here, we do not optimize and assume each resulting vertex is different.  However, we can optimize slightly; any control points shared on an edge between two adjacent patches can be shared in the triangle mesh.  In fact, for usability in the graphical interface, these should be the same if they are the same originally.  Thus, we need a data structure that tells us which edges are shared and update the generated control point indices accordingly.

          std::vector<Vertex> vertices;
          vertices.reserve(m_surface.m_controlPointIndices.size());

          auto c = m_surface.m_controlPointIndices.cbegin();
          for(auto i = m_surface.m_parameterPointIndices.cbegin(), end = m_surface.m_parameterPointIndices.cend(); i != end; ++i) {
               const Point2D& p0 = m_surface.m_parameterPoints[*i];
               const Point2D& p1 = m_surface.m_parameterPoints[*(++i)];
               const Point2D& p2 = m_surface.m_parameterPoints[*(++i)];

               constexpr hpreal nth = 1.0 / (hpreal) t_degree;
               for(hpint i2 = 0; i2 <= t_degree; ++i2) {
                    hpint i1 = 0;
                    hpint i0 = t_degree - i2;
                    while(i0 >= 0) {
                         vertices.push_back(factory(((hpreal)i0 * p0 + (hpreal)i1 * p1 + (hpreal)i2 * p2) * nth, m_surface.m_controlPoints[*c]));
                         ++i1;
                         --i0;
                         ++c;
                    }
               }
          }

          std::vector<hpuint> controlPointIndices;
          controlPointIndices.reserve(m_surface.m_controlPointIndices.size());

          for(hpint i = 0, end = m_surface.m_controlPointIndices.size(); i < end; ++i) controlPointIndices.push_back(i);
          std::vector<hpuint> indices = SurfaceSplineUtilsBEZ::template buildTriangleMeshIndices<t_degree>(controlPointIndices);

          return make_triangle_mesh(std::move(vertices), std::move(indices));
     }

private:
     const SurfaceSplineBEZ& m_surface;

};//TriangleMeshBuilder*/

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, bool t_includeParameterPoints = is_relative_vertex<Vertex>::value, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions = 4, VertexFactory&& factory = VertexFactory()) {
     if(nSubdivisions > 0) {
          auto subdivided = subdivide(surface, nSubdivisions);
          return TriangleMeshBuilder<Vertex, t_includeParameterPoints>::build(subdivided, std::forward<VertexFactory>(factory));
     } else return TriangleMeshBuilder<Vertex, t_includeParameterPoints>::build(surface, std::forward<VertexFactory>(factory));
}

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, bool t_includeParameterPoints = is_relative_vertex<Vertex>::value>
TriangleMesh<Vertex> make_control_polygon(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory()) { return make_triangle_mesh<Space, degree, Vertex, VertexFactory, t_includeParameterPoints>(surface, 0, std::forward<VertexFactory>(factory)); }

}//namespace happah

