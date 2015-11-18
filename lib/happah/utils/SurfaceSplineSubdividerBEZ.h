// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <iterator>
#include <vector>

#include "happah/utils/SurfaceSubdividerBEZ.h"

namespace happah {

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
     template<class Iterator>
     SurfaceSplineSubdividerBEZ(Iterator patches, hpuint nPatches) 
          : m_nPatches(nPatches) {
          m_subdividers.reserve(m_nPatches);
          for(hpuint i = 0; i < m_nPatches; ++i) {
               m_subdividers.push_back(SurfaceSubdividerBEZ<Space, t_degree>(*patches));
               ++patches;
          }
     }

     std::pair<ControlPoints, Indices> subdivide(hpuint nSubdivisions) {
          assert(nSubdivisions > 0);

          ControlPoints points;
          Indices indices;

          points.reserve(m_nControlPoints * m_nPatches);
          indices.reserve(3 * m_nTriangles * m_nPatches);

          for(SurfaceSubdividerBEZ<Space, t_degree>& subdivider : m_subdividers) {
               auto subdivided = subdivider.subdivide(nSubdivisions);
               hpuint offset = points.size();
               for(auto& i : subdivided.second) i += offset;
               std::move(subdivided.first.begin(), subdivided.first.end(), std::back_inserter(points));
               std::move(subdivided.second.begin(), subdivided.second.end(), std::back_inserter(indices));
          }

          return std::make_pair(std::move(points), std::move(indices));
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

}//namespace happah

