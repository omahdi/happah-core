// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.11 - Hedwig Amberg    - Introduce BezierQuadMesh.

#pragma once

#include <numeric>

#include "happah/geometry/QuadMesh.hpp"
#include "happah/util/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

template<class Space, hpuint t_degree0, hpuint t_degree1>
class BezierQuadMesh;

namespace bqm {

class QuadsEnumerator;

}//namespace bqm

template<class Space, hpuint degree0, hpuint degree1>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(hpuint n);

template<class Space, hpuint degree0, hpuint degree1, class Point = typename Space::POINT>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(std::vector<Point> controlPoints);

template<class Space, hpuint degree0, hpuint degree1, class Point = typename Space::POINT>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(std::vector<Point> controlPoints, Tuples<hpindex> indices);

template<class Space, hpuint degree0, hpuint degree1, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
QuadMesh<Vertex> make_control_polygon(const BezierQuadMesh<Space, degree0, degree1>& mesh, VertexFactory&& build = VertexFactory());

inline constexpr hpuint make_control_polygon_size(hpuint degree0, hpuint degree1);

inline constexpr hpuint make_patch_size(hpuint degree0, hpuint degree1);

template<class Space, hpuint degree0, hpuint degree1, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, 
typename = typename std::enable_if<( (degree0 > 0) && (degree1 > 0) )>::type>
QuadMesh<Vertex> make_quad_mesh(const BezierQuadMesh<Space, degree0, degree1>& mesh, hpuint nSubdivisions, VertexFactory&& factory = VertexFactory());

inline bqm::QuadsEnumerator make_quads_enumerator(hpuint degree0, hpuint degree1);

template<class Space, hpuint degree0, hpuint degree1>
auto size(const BezierQuadMesh<Space, degree0, degree1>& mesh);

//DEFINITIONS

template<class Space, hpuint t_degree0, hpuint t_degree1>
class BezierQuadMesh {
     using Point = typename Space::POINT;

public:
     BezierQuadMesh() {}

     BezierQuadMesh(std::vector<Point> controlPoints, Tuples<hpindex> indices)
          : m_controlPoints{std::move(controlPoints)}, m_indices{std::move(indices)} {}

     auto& getControlPoints() const { return m_controlPoints; }

     auto& getIndices() const { return m_indices; }

     hpuint getNumberOfControlPoints() const { return m_controlPoints.size(); }//TODO; there may be fewer control points on the actual surface

     hpuint getNumberOfPatches() const { return m_indices.size() / make_patch_size(t_degree0, t_degree1); }

private:
     std::vector<Point> m_controlPoints;
     Tuples<hpindex> m_indices;

};//BezierQuadMesh

namespace bqm {

class QuadsEnumerator {
public:
     QuadsEnumerator(hpuint degree0, hpuint degree1)
     : m_delta(degree0 + 1), m_end(m_delta * degree1), m_i(1) {}

     explicit operator bool() const { return m_i < m_end; }

     auto operator*() const { return std::make_tuple(m_i - 1, m_i, m_i + m_delta, m_i + m_delta - 1); }

     auto& operator++() {
          ++m_i;
          if(m_i % m_delta == 0) ++m_i;
          return *this;
     }

private:
     hpuint m_delta;
     hpuint m_end;
     hpindex m_i;

};//QuadsEnumerator

}//namespace bqm

template<class Space, hpuint degree0, hpuint degree1>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(hpuint n) {
     using Point = typename Space::POINT;

     auto length = make_patch_size(degree0, degree1);

     return { { 1, Point(0) }, { length, n * length, 0 } };
}

template<class Space, hpuint degree0, hpuint degree1, class Point>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(std::vector<Point> controlPoints) {
     auto length = make_patch_size(degree0, degree1);
     auto indices = Tuples<hpindex>(length, controlPoints.size());//TODO: check correct number of control points

     std::iota(std::begin(indices), std::end(indices), 0);

     return { std::move(controlPoints), std::move(indices) };
}

template<class Space, hpuint degree0, hpuint degree1, class Point>
BezierQuadMesh<Space, degree0, degree1> make_bezier_quad_mesh(std::vector<Point> controlPoints, Indices indices) { return { std::move(controlPoints), std::move(indices) }; }

template<class Space, hpuint degree0, hpuint degree1, class Vertex, class VertexFactory>
QuadMesh<Vertex> make_control_polygon(const BezierQuadMesh<Space, degree0, degree1>& mesh, VertexFactory&& build) {
     auto vertices = std::vector<Vertex>();
     auto indices = Quadruples<hpindex>();

     vertices.reserve(mesh.getNumberOfControlPoints());
     indices.reserve((make_control_polygon_size(degree0, degree1) << 2) * size(mesh));
     for(auto& point : mesh.getControlPoints()) vertices.push_back(build(point));
     visit(mesh.getIndices(), [&](auto patch) {
          visit(make_quads_enumerator(degree0, degree1), [&](auto i0, auto i1, auto i2, auto i3) {
               indices.insert(std::end(indices), { patch[i0], patch[i1], patch[i2], patch[i3] });
          });
     });
     return make_quad_mesh(std::move(vertices), std::move(indices));
}

inline constexpr hpuint make_control_polygon_size(hpuint degree0, hpuint degree1) { return degree0 * degree1; }

inline constexpr hpuint make_patch_size(hpuint degree0, hpuint degree1) { return (degree0 + 1) * (degree1 + 1); }

template<class Space, hpuint degree0, hpuint degree1, class Vertex, class VertexFactory, typename>
QuadMesh<Vertex> make_quad_mesh(const BezierQuadMesh<Space, degree0, degree1>& mesh, hpuint nSubdivisions, VertexFactory&& factory) {
     
     //TODO
     
}

inline bqm::QuadsEnumerator make_quads_enumerator(hpuint degree0, hpuint degree1) { return { degree0, degree1 }; }

template<class Space, hpuint degree0, hpuint degree1>
auto size(const BezierQuadMesh<Space, degree0, degree1>& mesh) { return mesh.getNumberOfPatches(); }

}//namespace happah

