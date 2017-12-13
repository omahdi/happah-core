// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.11 - Hedwig Amberg    - Introduce BezierQuadMesh.

#pragma once

#include "happah/geometry/QuadMesh.hpp"

namespace happah {

//DECLARATIONS

template<class Space, hpuint t_degree0, hpuint t_degree1>
class BezierQuadMesh;

template<class Space, hpuint degree0, hpuint degree1>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(hpuint n);

template<class Space, hpuint degree0, hpuint degree1, class Point = typename Space::POINT>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(std::vector<Point> controlPoints);

template<class Space, hpuint degree0, hpuint degree1, class Point = typename Space::POINT>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(std::vector<Point> controlPoints, Indices indices);

template<class Space, hpuint degree0, hpuint degree1, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
QuadMesh<Vertex> make_control_polygon(const BezierQuadMesh<Space, degree0, degree1>& surface, VertexFactory&& build = VertexFactory());

inline constexpr hpuint make_control_polygon_size(hpuint degree0, hpuint degree1) { return degree0 * degree1; }

inline constexpr hpuint make_patch_size(hpuint degree0, hpuint degree1) { return (degree0 + 1) * (degree1 + 1); }

template<class Space, hpuint degree0, hpuint degree1, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, 
typename = typename std::enable_if<( (degree0 > 0) && (degree1 > 0) )>::type>
QuadMesh<Vertex> make_quad_mesh(const BezierQuadMesh<Space, degree0, degree1>& surface, hpuint nSubdivisions, VertexFactory&& factory = VertexFactory());

template<class Space, hpuint degree0, hpuint degree1>
auto size(const BezierQuadMesh<Space, degree0, degree1>& surface);

template<hpuint degree0, hpuint degree1, class Iterator, class Visitor>
void visit_patches(Iterator patches, hpuint nPatches, Visitor&& visit);

//DEFINITIONS

template<class Space, hpuint t_degree0, hpuint t_degree1>
class BezierQuadMesh {
     using Point = typename Space::POINT;

public:
     BezierQuadMesh() {}

     BezierQuadMesh(std::vector<Point> controlPoints, Indices indices)
          : m_controlPoints{std::move(controlPoints)}, m_indices{std::move(indices)} {}

     auto& getControlPoints() const { return m_controlPoints; }

     std::tuple<const std::vector<Point>&, const Indices&> getPatches() const { return std::tie(m_controlPoints, m_indices); }

private:
     std::vector<Point> m_controlPoints;
     Indices m_indices;

};//BezierQuadMesh

template<class Space, hpuint degree0, hpuint degree1>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(hpuint n) {
     using Point = typename Space::POINT;

     return make_bezier_triangle_mesh({ 1, Point(0) }, { n * make_patch_size(degree0, degree1), 0 });
}

template<class Space, hpuint degree0, hpuint degree1, class Point>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(std::vector<Point> controlPoints) {
     auto indices = Indices(controlPoints.size());//TODO: check correct number of control points

     std::iota(std::begin(indices), std::end(indices), 0);

     return make_bezier_quad_mesh(std::move(controlPoints), std::move(indices));
}

template<class Space, hpuint degree0, hpuint degree1, class Point>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(std::vector<Point> controlPoints, Indices indices) { return { std::move(controlPoints), std::move(indices) }; }

template<class Space, hpuint degree0, hpuint degree1, class Vertex, class VertexFactory>
QuadMesh<Vertex> make_control_polygon(const BezierQuadMesh<Space, degree0, degree1>& surface, VertexFactory&& build) {
     //return make_quad_mesh<Space, degree0, degree1, Vertex, VertexFactory>(surface, 0, std::forward<VertexFactory>(build));
     using Point = typename Space::POINT;

     auto vertices = std::vector<Vertex>();
     vertices.reserve(surface.getControlPoints().size());
     for(auto& point : surface.getControlPoints()) { vertices.push_back(build(point)); }

     auto indices = Indices();
     indices.reserve(4 * make_control_polygon_size(degree0, degree1) * size(surface));
     auto inserter = make_back_inserter(indices);
     visit_patches<degree0, degree1>(std::begin(std::get<1>(surface.getPatches())), size(surface), [&](auto patch) {
          visit_quads(degree0, degree1, patch, inserter);
     });

     return make_quad_mesh(std::move(vertices), std::move(indices));
}

template<class Space, hpuint degree0, hpuint degree1, class Vertex, class VertexFactory, typename>
QuadMesh<Vertex> make_quad_mesh(const BezierQuadMesh<Space, degree0, degree1>& surface, hpuint nSubdivisions, VertexFactory&& factory) {
     
     //TODO
     
}

template<class Space, hpuint degree0, hpuint degree1>
auto size(const BezierQuadMesh<Space, degree0, degree1>& surface) { return std::get<1>(surface.getPatches()).size() / make_patch_size(degree0, degree1); }

template<hpuint degree0, hpuint degree1, class Iterator, class Visitor>
void visit_patches(Iterator patches, hpuint nPatches, Visitor&& visit) {
     static constexpr auto patchSize = make_patch_size(degree0, degree1);
     for(auto i = patches, end = patches + nPatches * patchSize; i != end; i += patchSize) visit(i);
}

template<class Iterator, class Visitor>
void visit_quads(hpuint degree0, hpuint degree1, Iterator patch, Visitor&& visit) {
     hpuint rowSize = degree0 + 1;
     hpuint colSize = degree1 + 1;
     
     auto N = [&](hpuint col, hpuint row) { return *(patch + row * rowSize + col); };
     
     for(hpuint r = 0; r < colSize-1; ++r) {
          for(hpuint c = 0; c < rowSize-1; ++c) {
               visit(N(c, r), N(c+1, r), N(c+1, r+1), N(c, r+1));
          }
     }
}

}//namespace happah
