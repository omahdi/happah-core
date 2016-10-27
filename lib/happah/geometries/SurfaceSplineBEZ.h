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
     using ControlPoints = std::vector<Point>;

public:
     SurfaceSplineBEZ() {}

     //NOTE: There is no automatic way of figuring out the neighborhood of patches given only the control points.
     SurfaceSplineBEZ(ControlPoints controlPoints, Indices controlPointIndices)
          : m_controlPointIndices{std::move(controlPointIndices)}, m_controlPoints{std::move(controlPoints)} {}

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     hpuint getNumberOfPatches() const { return m_controlPointIndices.size() / SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value; }

     std::tuple<const ControlPoints&, const Indices&> getPatches() const { return std::tie(m_controlPoints, m_controlPointIndices); }

private:
     Indices m_controlPointIndices;
     ControlPoints m_controlPoints;

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const SurfaceSplineBEZ<Space, t_degree>& surface) {
          stream << surface.m_controlPointIndices << '\n';//TODO: not good pushing '\n'
          stream << surface.m_controlPoints;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, SurfaceSplineBEZ<Space, t_degree>& surface) {
          stream >> surface.m_controlPointIndices;
          stream >> surface.m_controlPoints;
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

template<hpuint degree, class Iterator, class Visitor>
void visit_patches(Iterator begin, hpuint nPatches, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     for(auto i = begin, end = begin + nPatches * nControlPoints; i != end; i += nControlPoints) visit(i);
}

template<hpuint degree, class Iterator, class Visitor>
void visit_patches(Iterator begin, Iterator end, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     for(auto i = begin; i != end; i += nControlPoints) visit(i);
}

template<class Space, hpuint degree, class Visitor>
void visit_patches(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto patches = surface.getPatches();
     auto temp = deindex(std::get<0>(patches), std::get<1>(patches));
     visit_patches<degree>(temp.begin(), temp.end(), std::forward<Visitor>(visit));
}

template<class Space, hpuint t_degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, t_degree>& surface, Visitor&& visit) {}

template<hpuint degree, class ControlPointsIterator, class DomainPointsIterator, class Visitor>
void sample(ControlPointsIterator controlPoints, DomainPointsIterator domainPoints, hpuint nPatches, hpuint nSamples, Visitor&& visit) {
     //TODO; skip multiple computations on common edges and eliminate common points in array
     //TODO: t_degree == 1,2,4, general
     //TODO: optimize calculations below; (u,v,w) do not have to be recalculated each time
     auto matrix = SurfaceUtilsBEZ::getEvaluationMatrix<degree>(nSamples);
     visit_patches<degree>(controlPoints, nPatches, [&](auto c) {
          auto& p0 = *domainPoints;
          auto& p1 = *(++domainPoints);
          auto& p2 = *(++domainPoints);
          ++domainPoints;

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
                    visit(mix(p0, u, p1, v, p2, w), sample);
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
                    visit(mix(p0, u, p1, v, p2, w), sample);
               });
               break;
          }
          default:

               break;
          }
     });
}

template<class Space, hpuint degree, class T, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, std::tuple<const std::vector<T>&, const Indices&> domain, hpuint nSamples, Visitor&& visit) {
     auto patches = surface.getPatches();
     auto controlPoints = deindex(std::get<0>(patches), std::get<1>(patches)).begin();
     auto domainPoints = deindex(std::get<0>(domain), std::get<1>(domain)).begin();
     sample<degree>(controlPoints, domainPoints, surface.getNumberOfPatches(), nSamples, std::forward<Visitor>(visit));
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

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions) {
     if(nSubdivisions == 0) return surface;
     auto subdivider = SurfaceSplineSubdividerBEZ<Space, degree>(surface);
     auto subdivided = subdivider.subdivide(nSubdivisions);
     return { std::move(subdivided.first), std::move(subdivided.second) };
}

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory()) {
     using Point = typename Space::POINT;
     static_assert(std::is_base_of<Vertex, decltype(factory(Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");

     std::vector<Vertex> vertices;
     vertices.reserve(surface.getControlPoints().size());
     for(auto& point : surface.getControlPoints()) vertices.push_back(factory(point));

     auto indices = SurfaceSplineUtilsBEZ::template buildTriangleMeshIndices<degree>(std::get<1>(surface.getPatches()));

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions, VertexFactory&& factory = VertexFactory()) {
     if(nSubdivisions > 0) return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(subdivide(surface, nSubdivisions), std::forward<VertexFactory>(factory));
     else return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, std::forward<VertexFactory>(factory));
}

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
TriangleMesh<Vertex> make_control_polygon(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory()) { return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, std::forward<VertexFactory>(factory)); }

}//namespace happah

