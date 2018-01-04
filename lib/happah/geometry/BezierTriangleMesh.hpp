// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.08 - Hedwig Amberg    - added NablasEnumerator and paint_edges function.

//DICTIONARY
//   border: control points of a given patch where at least one common index is zero
//   corner: one of the three control points of a given patch which have two zero indices
//   end: one of the two control points on a given border of a given patch which have two zero indices
//   boundary: border excluding the ends
//   interior: patch excluding all borders

#pragma once

#include <array>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/iterator_range.hpp>
#include <cmath>
#include <experimental/filesystem>
#include <glm/gtc/constants.hpp>
#include <limits>
#include <numeric>
#include <sstream>
#include <stack>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "happah/Eigen.hpp"
#include "happah/Happah.hpp"
#include "happah/format/hph.hpp"
#include "happah/geometry/TriangleGraph.hpp"
#include "happah/geometry/TriangleMesh.hpp"
#include "happah/math/ProjectiveStructure.hpp"
#include "happah/util/BezierTriangleSubdivider.hpp"
#include "happah/util/ProxyArray.hpp"
#include "happah/util/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

template<class Space, hpuint t_degree>
class BezierTriangleMesh;

namespace btm {

class DeltasEnumerator;

class EdgeDiamondsEnumerator;

class NablasEnumerator;

template<hpindex t_ring>
class RingWalker;

template<hpindex t_ring>
class RingEnumerator;

template<hpindex t_i>
class RowEnumerator;

template<hpindex t_ring>
class VertexDiamondsEnumerator;

}//namespace btm

template<hpuint degree, class Iterator>
auto de_casteljau(Iterator patch, hpreal u, hpreal v, hpreal w);

template<class Space, hpuint degree>
auto de_casteljau(const BezierTriangleMesh<Space, degree>& mesh, hpuint p, hpreal u, hpreal v);

template<class Space, hpuint degree>
auto deindex(const BezierTriangleMesh<Space, degree>& mesh);

template<class Space, hpuint degree>
BezierTriangleMesh<Space, (degree + 1)> elevate(const BezierTriangleMesh<Space, degree>& mesh);

template<class Space, hpuint degree>
BezierTriangleMesh<Space, (degree + 1)> elevate(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors);

template<class NewSpace, class OldSpace, hpuint degree, class Transformer>
BezierTriangleMesh<NewSpace, degree> embed(const BezierTriangleMesh<OldSpace, degree>& mesh, Transformer&& transform);

template<class Space, hpuint degree>
bool is_c0(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors, hpindex p, trit i);

template<class Space, hpuint degree>
bool is_g1(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors, hpindex p, trit i, hpuint nSamples = 5, hpreal epsilon = EPSILON);

template<hpuint degree>
void make_bernstein_polynomials(const std::string& directory);

template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(hpuint n);

template<class Space, hpuint degree, class Point = typename Space::POINT>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(std::vector<Point> controlPoints);

template<class Space, hpuint degree, class Point = typename Space::POINT>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(std::vector<Point> controlPoints, Tuples<hpindex> indices);

//Convert a string representation in HPH format of a Bezier triangle mesh.
template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(const std::string& mesh);

//Import a Bezier triangle mesh stored in the given file in HPH format.
template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(const std::experimental::filesystem::path& mesh);

//Convert a closed(!) triangle graph into a quartic Bezier triangle mesh.
template<class Vertex>
auto make_bezier_triangle_mesh(const TriangleGraph<Vertex>& graph);

template<hpuint degree, class Iterator>
auto make_boundary(Iterator patch, trit i);

constexpr hpindex make_control_point_index(hpuint degree, hpindex i0, hpindex i1, hpindex i2);

//Return the index of the jth point on the ith boundary.
inline hpindex make_control_point_index(hpuint degree, trit i, hpindex j);

//Return the index of the ith point in the interior.
inline hpindex make_control_point_index(hpuint degree, hpindex i);

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
TriangleMesh<Vertex> make_control_polygon(const BezierTriangleMesh<Space, degree>& mesh, VertexFactory&& factory = VertexFactory());

constexpr hpuint make_control_polygon_size(hpuint degree);

/**
 * @param[in] nSamples Number of times an edge of the parameter triangle should be sampled.  The entire triangle is sampled uniformly such that this parameter is respected.
 * @return Matrix whose rows are the Bernstein polynomials evaluated at some point u.  The matrix is returned row-major.  To evaluate a B\'ezier polynomial at the sampled u values given a vector of control points, simply compute the product of the matrix with the vector of control points.
 */
std::vector<hpreal> make_de_casteljau_matrix(hpuint degree, hpuint nSamples);

inline btm::DeltasEnumerator make_deltas_enumerator(hpuint degree);

inline btm::EdgeDiamondsEnumerator make_diamonds_enumerator(hpuint degree, trit i, trit j);

template<class Space, hpuint degree>
auto make_diamonds_enumerator(const BezierTriangleMesh<Space, degree>& mesh, hpindex p, trit i, hpindex q, trit j);

template<hpindex ring>
inline btm::VertexDiamondsEnumerator<ring> make_diamonds_enumerator(hpuint degree, const Triples<hpindex>& neighbors, hpindex p, trit i);

template<hpindex ring, class Space, hpuint degree>
auto make_diamonds_enumerator(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors, hpindex p, trit i);

inline btm::NablasEnumerator make_nablas_enumerator(hpuint degree);

template<class Space, hpuint degree>
Triples<hpindex> make_neighbors(const BezierTriangleMesh<Space, degree>& mesh);

constexpr hpuint make_patch_size(hpuint degree);

template<hpindex ring = 1>
btm::RingEnumerator<ring> make_ring_enumerator(hpuint degree, const Triples<hpindex>& neighbors, hpindex p, trit i);

template<hpindex ring, class Space, hpuint degree>
auto make_ring_enumerator(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors, hpindex p, trit i);

template<hpindex i>
btm::RowEnumerator<i> make_row_enumerator(hpuint degree, hpindex row);

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const BezierTriangleMesh<Space, degree>& mesh, hpuint nSubdivisions, VertexFactory&& factory = VertexFactory());

std::tuple<std::vector<hpcolor>, std::vector<hpcolor> > paint_boundary_edges(hpuint degree, std::vector<hpcolor> Vcolors, std::vector<hpcolor> Ecolors, const hpcolor& color);

std::vector<hpcolor> paint_boundary_triangles(hpuint degree, std::vector<hpcolor> colors, const hpcolor& color0, const hpcolor& color1);

template<class Space, hpuint degree, class Visitor>
void sample(const BezierTriangleMesh<Space, degree>& mesh, hpuint nSamples, Visitor&& visit);

template<class Space, hpuint degree>
auto size(const BezierTriangleMesh<Space, degree>& mesh);

//Return a G1 surface that interpolates the positions and the tangents planes at the corners of the patches in the given surface.
template<hpuint degree>
BezierTriangleMesh<Space4D, degree> smooth(BezierTriangleMesh<Space4D, degree> mesh, const ProjectiveStructure& structure, hpreal epsilon = EPSILON);

template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> subdivide(const BezierTriangleMesh<Space, degree>& mesh, hpuint nSubdivisions);

template<class Space, hpuint degree, class Visitor>
void visit(const BezierTriangleMesh<Space, degree>& mesh, Visitor&& visit);

template<class Visitor>
void visit_bernstein_indices(hpuint degree, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_boundary(Iterator patch, trit i, Visitor&& visit);

//Visit triangles in control polygon schematically pointing up.
//template<class Iterator, class Visitor>
//void visit_deltas(hpuint degree, Iterator patch, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_ends(Iterator patch, trit i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_interior(Iterator patch, Visitor&& visit);

//Visit triangles in control polygon schematically pointing down.  The points are given in counterclockwise order; the first point is the top right point.
//template<class Iterator, class Visitor>
//void visit_nablas(hpuint degree, Iterator patch, Visitor&& visit);

//Visit the subring stopping at patch p.
template<class Visitor>
void visit_subring(btm::RingEnumerator<1> e, hpindex p, Visitor&& visit);

//Choose constant weight.
template<hpuint degree>
BezierTriangleMesh<Space4D, degree> weigh(const BezierTriangleMesh<Space3D, degree>& mesh);

//Choose weights such that diamonds along edges satisfy transition constraints and C1 constraints as much as possible.
template<hpuint degree>
BezierTriangleMesh<Space4D, degree> weigh(BezierTriangleMesh<Space4D, degree> mesh, const ProjectiveStructure& structure);

//DEFINITIONS

template<class Space, hpuint t_degree>
class BezierTriangleMesh {
     using Point = typename Space::POINT;

public:
     BezierTriangleMesh() {}

     BezierTriangleMesh(std::vector<Point> controlPoints, Tuples<hpindex> indices)
          : m_controlPoints{std::move(controlPoints)}, m_indices{std::move(indices)} {}

     auto& getControlPoint(hpindex p, hpindex i) const { return m_controlPoints[m_indices(p, i)]; }

     auto& getControlPoint(hpindex p, hpindex i) { return m_controlPoints[m_indices(p, i)]; }

     //Get corner control point.
     auto& getControlPoint(hpindex p, trit i) const {
          static constexpr hpuint o[3] = { 0u, t_degree, make_patch_size(t_degree) - 1u };

          return m_controlPoints[m_indices(p, o[i])];
     }

     //Get corner control point.
     auto& getControlPoint(hpindex p, trit i) {
          static constexpr hpuint o[3] = { 0u, t_degree, make_patch_size(t_degree) - 1u };

          return m_controlPoints[m_indices(p, o[i])];
     }

     //Get boundary control point.
     auto& getControlPoint(hpindex p, trit i, hpindex j) const {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");

          return m_controlPoints[m_indices(p, make_control_point_index(t_degree, i, j))];
     }

     //Get boundary control point.
     auto& getControlPoint(hpindex p, trit i, hpindex j) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");

          return m_controlPoints[m_indices(p, make_control_point_index(t_degree, i, j))];
     }

     auto& getControlPoints() const { return m_controlPoints; }

     auto& getControlPoints() { return m_controlPoints; }

     auto& getIndices() const { return m_indices; }

     auto& getIndices() { return m_indices; }

     hpuint getNumberOfControlPoints() const { return m_controlPoints.size(); }//TODO; there may be fewer control points on the actual surface

     hpuint getNumberOfPatches() const { return m_indices.size() / make_patch_size(t_degree); }

     auto getPatch(hpindex p) const { return deindex(*this)(p); }

     auto getPatch(hpindex p) { return deindex(*this)(p); }

     auto getPatches() const { return deindex(*this); }

     auto getPatches() { return deindex(*this); }

     //Set the ith boundary of the pth patch.
     template<class Iterator>
     void setBoundary(hpindex p, trit i, Iterator begin) {//TODO: rename setControlPoints, trit
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");

          auto n = m_controlPoints.size();

          visit_boundary<t_degree>(m_indices(p), i, [&](auto& i) { i = n++; });
          m_controlPoints.insert(std::end(m_controlPoints), begin, begin + (t_degree - 1));
     }

     //Set the ith boundary of the pth patch to the jth boundary of the qth patch.
     void setBoundary(hpindex p, trit i, hpindex q, trit j) {//TODO: rename setControlPoints, trit
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");

          auto boundary = make_boundary<t_degree>(m_indices(q), j);
          auto n = std::end(boundary);

          visit_boundary<t_degree>(m_indices(p), i, [&](auto& i) { i = *(--n); });
     }

     void setControlPoint(hpindex p, hpindex i, Point point) {
          m_indices(p, i) = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     //Set the kth point on the ith boundary of the pth patch.
     void setControlPoint(hpindex p, trit i, hpindex k, Point point) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");

          auto o = make_control_point_index(t_degree, i, k);

          m_indices(p, o) = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     //Set the kth point on the ith boundary of the pth patch and the point opposite to it on the jth boundary of the qth patch.
     void setControlPoint(hpindex p, trit i, hpindex k, hpindex q, trit j, Point point) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");

          auto o0 = make_control_point_index(t_degree, i, k);
          auto o1 = make_control_point_index(t_degree, j, t_degree - 2 - k);

          //std::cout << "boundary: " << i << ' ' << o0 << ' ' << j << ' ' << o1 << '\n';

          m_indices(p, o0) = m_controlPoints.size();
          m_indices(q, o1) = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     //Set the kth point on the ith boundary of the pth patch to the opposite point on the jth boundary of the qth patch.
     void setControlPoint(hpindex p, trit i, hpindex k, hpindex q, trit j) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");

          auto o0 = make_control_point_index(t_degree, i, k);
          auto o1 = make_control_point_index(t_degree, j, t_degree - 2 - k);

          m_indices(p, o0) = m_indices(q, o1);
     }

     void setControlPoint(hpindex p, trit i, Point point) {
          static_assert(t_degree > 0, "There is no corner in a constant.");

          static constexpr hpindex o[3] = { 0u, t_degree, make_patch_size(t_degree) - 1u };

          m_indices(p, o[i]) = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     void setControlPoint(hpindex p, trit i, hpindex q, trit j) {
          static_assert(t_degree > 0, "There is no corner in a constant.");

          static constexpr hpindex o[3] = { 0u, t_degree, make_patch_size(t_degree) - 1u };

          m_indices(p, o[i]) = m_indices(q, o[j]);
     }

     template<class Iterator>
     void setInterior(hpindex p, Iterator begin) {
          static_assert(t_degree > 2, "There is no interior in a constant, linear, or quadratic.");

          auto n = m_controlPoints.size();

          visit_interior<t_degree>(m_indices(p), [&](auto& i) { i = n++; });
          m_controlPoints.insert(std::end(m_controlPoints), begin, begin + (make_patch_size(t_degree) - 3 * t_degree)); 
     }

private:
     std::vector<Point> m_controlPoints;
     Tuples<hpindex> m_indices;

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const BezierTriangleMesh<Space, t_degree>& mesh) {
          using happah::format::hph::operator<<;

          stream << mesh.m_controlPoints << '\n';
          stream << mesh.m_indices;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, BezierTriangleMesh<Space, t_degree>& mesh) {
          using happah::format::hph::operator>>;

          stream >> mesh.m_controlPoints;
          stream >> mesh.m_indices;
          return stream;
     }

};//BezierTriangleMesh
template<class Space>
using ConstantBezierTriangleMesh = BezierTriangleMesh<Space, 0>;
template<class Space>
using CubicBezierTriangleMesh = BezierTriangleMesh<Space, 3>;
template<class Space>
using LinearBezierTriangleMesh = BezierTriangleMesh<Space, 1>;
template<class Space>
using SexticBezierTriangleMesh = BezierTriangleMesh<Space, 6>;
template<class Space>
using QuadraticBezierTriangleMesh = BezierTriangleMesh<Space, 2>;
template<class Space>
using QuarticBezierTriangleMesh = BezierTriangleMesh<Space, 4>;
template<class Space>
using QuinticBezierTriangleMesh = BezierTriangleMesh<Space, 5>;

namespace btm {

class DeltasEnumerator {
public:
     DeltasEnumerator(hpuint degree)
          : m_bottom(0u), m_delta(degree), m_end(degree) { if(degree == hpuint(0)) throw std::runtime_error("There are no deltas in a constant Bezier triangle."); }
     
     auto operator*() const { return std::make_tuple(m_bottom, m_bottom + 1, m_bottom + m_delta + 1); }
          
     explicit operator bool() const { return m_delta > 0; }
     
     auto& operator++() {
          ++m_bottom;
          if(m_bottom == m_end) {
               --m_delta;
               ++m_bottom;
               m_end = m_bottom + m_delta;
          }
          return *this;
     }

private:
     hpindex m_bottom;
     hpindex m_delta;
     hpindex m_end;

};//DeltasEnumerator
 
/*
 *   k3
 * k0  k2
 *   k1
 * i is with respect to top patch; j is with respect to bottom patch.
 * k0 and k2 are control point indices with respect to the ith patch.
 */
class EdgeDiamondsEnumerator {
public:
     EdgeDiamondsEnumerator(hpuint degree, trit i, trit j)
          : m_i(i), m_j(j) {
          if(m_i == 0u) {
               m_end = degree;
               m_k0 = 0u;
               m_k2 = 1u;
               m_k3 = degree + 1u;
          } else if(m_i == 1u) {
               m_delta0 = degree;
               m_end = make_patch_size(degree) - 1u;
               m_k0 = degree;
               m_k2 = degree << 1;
               m_k3 = degree - 1u;
          } else {
               m_delta0 = 2u;
               m_end = 0u;
               m_k0 = make_patch_size(degree) - 1u;
               m_k2 = m_k0 - 2u;
               m_k3 = m_k0 - 1u;
          }
          if(m_j == 0u) m_k1 = degree << 1;
          else if(m_j == 1u) {
               m_delta1 = 1u;
               m_k1 = make_patch_size(degree) - 3u;
          } else {
               m_delta1 = degree + 2u;
               m_k1 = 1u;
          }
     }

     explicit operator bool() const { return m_k0 != m_end; }

     auto operator*() const { return std::make_tuple(m_k0, m_k1, m_k2, m_k3); }

     auto& operator++() {
          if(m_i == 0u) {
               ++m_k0;
               ++m_k2;
               ++m_k3;
          } else if(m_i == 1u) {
               m_k0 += m_delta0;
               m_k3 += m_delta0;
               --m_delta0;
               m_k2 += m_delta0;
          } else {
               m_k0 -= m_delta0;
               ++m_delta0; 
               m_k2 -= m_delta0;
               m_k3 -= m_delta0;
          }
          if(m_j == 0u) --m_k1;
          else if(m_j == 1u) m_k1 -= ++m_delta1;
          else m_k1 += --m_delta1;
          return *this;
     }

private:
     hpuint m_delta0;
     hpuint m_delta1;
     hpuint m_end;
     trit m_i;
     trit m_j;
     hpindex m_k0;
     hpindex m_k1;
     hpindex m_k2;
     hpindex m_k3;

};//EdgeDiamondsEnumerator

class NablasEnumerator {
public:
     NablasEnumerator(hpuint degree)
          : m_bottom(1u), m_delta(degree), m_end(degree) { if(degree < hpuint(2)) throw std::runtime_error("There are no nablas in constant or linear Bezier triangles."); }

     auto operator*() const { return std::make_tuple(m_bottom + m_delta + 1, m_bottom + m_delta, m_bottom); }

     explicit operator bool() const { return m_delta > 1; }
     
     auto& operator++() {
          ++m_bottom;
          if(m_bottom == m_end) {
               --m_delta;
               ++m_bottom;
               m_end = m_bottom + m_delta;
               ++m_bottom;
          }
          return *this;
     }

private:
     hpindex m_bottom;
     hpindex m_delta;
     hpindex m_end;

};//NablasEnumerator

template<>
class RingWalker<1> {
public:
     RingWalker(hpuint degree, trm::SpokesWalker i)
          : m_i(std::move(i)), m_o{ 1u, degree << 1, make_patch_size(degree) - 3u } {}

     //NOTE: Equality of degrees is not confirmed.
     auto operator==(const RingWalker& walker) const { return m_i == walker.m_i; }

     auto operator!=(const RingWalker& walker) const { return !(*this == walker); }

     auto operator*() const {
          auto p = std::get<0>(*m_i);
          auto i = std::get<1>(*m_i);

          return std::make_tuple(p, m_o[i]);
     }

     auto& operator++() {
          ++m_i;
          return *this;
     }

     auto& operator--() {
          --m_i;
          return *this;
     }

     auto& operator+=(hpuint n) {
          while(n--) ++(*this);
          return *this;
     }

     auto& operator-=(hpuint n) {
          while(n--) --(*this);
          return *this;
     }

     auto operator+(hpuint n) const {
          auto copy = *this;
          return copy += n;
     }

     auto operator-(hpuint n) const {
          auto copy = *this;
          return copy -= n;
     }

private:
     trm::SpokesWalker m_i;
     const hpindex m_o[3];

};//RingWalker

template<>
class RingWalker<2> {
public:
     RingWalker(hpuint degree, trm::SpokesWalker i)
          : m_i(std::move(i)), m_o0{ 2, 3 * degree - 1, make_patch_size(degree) - 6 }, m_o1{ degree + 2, (degree << 1) - 1, make_patch_size(degree) - 5 } { m_o = m_o0.data(); }

     RingWalker(const RingWalker& walker)
          : m_i(walker.m_i), m_o0(walker.m_o0), m_o1(walker.m_o1) { m_o = (walker.m_o == walker.m_o0.data()) ? m_o0.data() : m_o1.data(); }

     RingWalker(RingWalker&& walker)
          : m_i(walker.m_i), m_o0(walker.m_o0), m_o1(walker.m_o1) { m_o = (walker.m_o == walker.m_o0.data()) ? m_o0.data() : m_o1.data(); }

     //NOTE: Equality of degrees is not confirmed.
     auto operator==(const RingWalker& walker) const { return m_i == walker.m_i && ((m_o == m_o0.data() && walker.m_o == walker.m_o0.data()) || (m_o == m_o1.data() && walker.m_o == walker.m_o1.data())); }

     auto operator!=(const RingWalker& walker) const { return !(*this == walker); }

     auto operator*() const {
          auto p = std::get<0>(*m_i);
          auto i = std::get<1>(*m_i);

          return std::make_tuple(p, m_o[i]);
     }

     auto& operator++() {
          if(m_o == m_o0.data()) m_o = m_o1.data();
          else {
               m_o = m_o0.data();
               ++m_i;
          }
          return *this;
     }

     auto& operator--() {
          if(m_o == m_o1.data()) m_o = m_o0.data();
          else {
               m_o = m_o1.data();
               --m_i;
          }
          return *this;
     }

     auto& operator+=(hpuint n) {
          while(n--) ++(*this);
          return *this;
     }

     auto& operator-=(hpuint n) {
          while(n--) --(*this);
          return *this;
     }

     auto operator+(hpuint n) const {
          auto copy = *this;
          return copy += n;
     }

     auto operator-(hpuint n) const {
          auto copy = *this;
          return copy -= n;
     }

private:
     trm::SpokesWalker m_i;
     const hpindex* m_o;
     const std::array<hpindex, 3> m_o0;
     const std::array<hpindex, 3> m_o1;

};//RingWalker

template<hpindex t_ring>
class RingEnumerator {
public:
     RingEnumerator(RingWalker<t_ring> i)
          : m_begin(i), m_i(std::move(i)) {}

     explicit operator bool() const { return m_i != m_begin; }

     auto operator*() const { return *m_i; }

     auto& operator++() {
          ++m_i;
          return *this;
     }

     auto& operator+=(hpuint n) {
          m_i += n;
          return *this;
     }

     auto operator+(hpuint n) const {
          auto copy = *this;
          return copy += n;
     }

private:
     RingWalker<t_ring> m_begin;
     RingWalker<t_ring> m_i;

};//RingEnumerator

template<>
class RowEnumerator<0> {
public:
     RowEnumerator(hpuint degree, hpindex row)
          : m_i(row * (degree + 1) - ((row * (row - 1)) >> 1)), m_end(m_i + degree + 1 - row) {}

     explicit operator bool() const { return m_i != m_end; }

     auto operator*() const { return m_i; }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     hpindex m_i;
     hpindex m_end;

};//RowEnumerator

template<>
class RowEnumerator<1> {
public:
     RowEnumerator(hpuint degree, hpindex row)
          : m_i(degree), m_delta(degree + 1) {}

     explicit operator bool() const { return m_delta > 0; }

     auto operator*() const { return m_i; }

     auto& operator++() {
          --m_delta;
          m_i += m_delta;
          return *this;
     }

private:
     hpuint m_delta;
     hpindex m_i;

};//RowEnumerator

template<>
class RowEnumerator<2> {
public:
     RowEnumerator(hpuint degree, hpuint row);

private:

};//RowEnumerator

template<>
class VertexDiamondsEnumerator<1> {
public:
     VertexDiamondsEnumerator(hpuint degree, trm::SpokesWalker i)
          : m_begin(degree, i), m_center(std::make_tuple(std::get<0>(*i), make_array(0u, degree, make_patch_size(degree) - 1u)[std::get<1>(*i)])), m_i(degree, std::move(i)) {
          --m_begin;
          --m_i;
     }

     explicit operator bool() const { return m_i != m_begin; }

     auto operator*() const {
          auto p0 = hpindex(0), p1 = hpindex(0), p2 = hpindex(0), p3 = hpindex(0), i0 = hpindex(0), i1 = hpindex(0), i2 = hpindex(0), i3 = hpindex(0);

          std::tie(p0, i0) = m_center;
          std::tie(p1, i1) = *m_i;
          std::tie(p2, i2) = *(m_i + 1);
          std::tie(p3, i3) = *(m_i + 2);
          return std::make_tuple(p0, i0, p1, i1, p2, i2, p3, i3);
     }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     RingWalker<1> m_begin;
     std::tuple<hpindex, hpindex> m_center;
     RingWalker<1> m_i;

};//VertexDiamondsEnumerator

template<>
class VertexDiamondsEnumerator<2> {
public:
     VertexDiamondsEnumerator(hpuint degree, trm::SpokesWalker i)
          : m_e({ degree, i }), m_i(degree, std::move(i)) { --m_i; }

     explicit operator bool() const { return bool(m_e); }

     auto operator*() const {
          auto p = hpindex(0), q = hpindex(0), i0 = hpindex(0), i1 = hpindex(0), i2 = hpindex(0), i3 = hpindex(0);

          std::tie(p, i0) = *m_e;
          std::tie(q, i1) = *m_i;
          std::tie(p, i2) = *(m_i + 1);
          std::tie(p, i3) = *(m_i + 2);
          return std::make_tuple(p, i0, q, i1, p, i2, p, i3);
     }

     auto& operator++() {
          ++m_e;
          ++m_i;
          ++m_i;
          return *this;
     }

private:
     RingEnumerator<1> m_e;
     RingWalker<2> m_i;

};//VertexDiamondsEnumerator

}//namespace btm

template<hpuint degree, class Iterator>
auto de_casteljau(Iterator patch, hpreal u, hpreal v, hpreal w) {
     using T = typename std::remove_const<typename std::remove_reference<decltype(*patch)>::type>::type;

     if(degree == 0u) return patch[0];

     auto points = std::array<T, make_patch_size(degree - 1u)>();

     auto do_de_casteljau = [&](auto i, auto patch) {
          auto p = std::begin(points) - 1;
          visit_deltas(make_deltas_enumerator(i), [&](auto b0, auto b1, auto b2) { (++p)[0] = u * patch[b0] + v * patch[b1] + w * patch[b2]; });
     };

     do_de_casteljau(degree, patch);
     if(degree > 1u) for(auto i : boost::irange(1u, degree) | boost::adaptors::reversed) do_de_casteljau(i, std::begin(points));
     return points[0];
}

template<class Space, hpuint degree>
auto de_casteljau(const BezierTriangleMesh<Space, degree>& mesh, hpuint p, hpreal u, hpreal v) { return de_casteljau<degree>(deindex(mesh)(p), u, v, 1.0f - u - v); }

template<class Space, hpuint degree>
auto deindex(const BezierTriangleMesh<Space, degree>& mesh) { return deindex(mesh.getControlPoints(), mesh.getIndices()); }

template<class Space, hpuint degree>
BezierTriangleMesh<Space, (degree + 1)> elevate(const BezierTriangleMesh<Space, degree>& mesh) { return elevate(mesh, make_neighbors(mesh)); }

template<class Space, hpuint degree>
BezierTriangleMesh<Space, (degree + 1)> elevate(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors) {
     using Point = typename Space::POINT;

     auto n = size(mesh);
     auto mesh1 = make_bezier_triangle_mesh<Space, (degree + 1)>(n);
     auto p = hpindex(0);

     auto elevate_boundary = [&](auto p, auto i) {
          auto patch = mesh.getPatch(p);
          auto boundary = make_boundary<degree>(patch, i);
          visit_ends<degree>(patch, i, [&](auto& corner0, auto& corner1) {
               auto alpha = hpreal(1.0) / hpreal(degree + 1);
               auto points = std::vector<Point>();
               auto middle = std::begin(boundary);

               points.reserve(degree);

               points.push_back(alpha * (corner0 + hpreal(degree) * middle[0]));
               for(auto i = hpuint(2); i < degree; ++i, ++middle) points.push_back(alpha * (hpreal(i) * middle[0] + hpreal(degree + 1 - i) * middle[1]));
               points.push_back(alpha * (hpreal(degree) * middle[0] + corner1));
               assert(points.size() == degree);

               mesh1.setBoundary(p, i, std::begin(points));
          });
     };

     visit(make_vertices_enumerator(neighbors), [&](auto p, auto i) {
          mesh1.setControlPoint(p, i, mesh.getControlPoint(p, i));
          visit(make_spokes_enumerator(neighbors, p, i), [&](auto q, auto j) { if(q < n) mesh1.setControlPoint(q, j, p, i); });
     });

     visit_edges(neighbors, [&](auto p, auto i) {
          elevate_boundary(p, i);
          auto q = neighbors(p, i);
          if(q >= n) return;
          auto j = make_offset(neighbors, q, p);
          if(is_c0(mesh, neighbors, p, i)) mesh1.setBoundary(q, j, p, i);
          else elevate_boundary(q, j);
     });

     if(degree < 2) return mesh1;

     visit(mesh, [&](auto patch) {
          auto alpha = hpreal(1.0) / hpreal(degree + 1);
          auto t0 = degree;
          auto j0 = hpuint(1), j1 = hpuint(-1), j2 = hpuint(0);
          auto e = transform(make_nablas_enumerator(degree), [&](auto i0, auto i1, auto i2) {
               --j0;
               ++j1;
               if(j0 == hpuint(0)) {
                    j0 = --t0;
                    j2 = degree - t0;
                    j1 = degree - j0 - j2 + hpuint(1);
               }
               return alpha * (hpreal(j0) * patch[i0] + hpreal(j1) * patch[i1] + hpreal(j2) * patch[i2]);
          });
          auto interior = make(e);

          mesh1.setInterior(p, std::begin(interior));
          ++p;
     });

     return mesh1;
}

template<class NewSpace, class OldSpace, hpuint degree, class Transformer>
BezierTriangleMesh<NewSpace, degree> embed(const BezierTriangleMesh<OldSpace, degree>& mesh, Transformer&& transform) {
     auto& oldPoints = mesh.getControlPoints();
     auto& indices = mesh.getIndices();
     auto newPoints = std::vector<typename NewSpace::POINT>(oldPoints.size());

     std::transform(std::begin(oldPoints), std::end(oldPoints), std::begin(newPoints), transform);
     return { std::move(newPoints), indices };
}

template<class Space, hpuint degree>
bool is_c0(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors, hpindex p, trit i) {
     auto& indices = mesh.getIndices();
     auto q = neighbors(p, i);
     auto j = make_offset(neighbors, q, p);

     auto k0 = 0u, k1 = 0u, l0 = 0u, l1 = 0u;
     visit_ends<degree>(indices(p), i, [&](auto k, auto l) { k0 = k; l0 = l; });
     visit_ends<degree>(indices(q), j, [&](auto k, auto l) { k1 = k; l1 = l; });
     if(k0 != l1 || l0 != k1) return false;
     auto boundary0 = make_boundary<degree>(indices(p), i); 
     auto boundary1 = make_boundary<degree>(indices(q), j);
     std::reverse(std::begin(boundary1), std::end(boundary1));
     return boundary0 == boundary1;
}

template<class Space, hpuint degree>
bool is_g1(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors, hpindex p, trit i, hpuint nSamples, hpreal epsilon) {
     auto q = neighbors(p, i);
     auto j = make_offset(neighbors, q, p);
     auto t = hpreal(0);
     auto delta = hpreal(1) / hpreal(nSamples - 1);
     auto vectors = std::array<Vector3D, 3 * degree>();
     auto v0 = std::begin(vectors), v1 = v0 + degree, v2 = v1 + degree;
     
     visit(make_diamonds_enumerator(degree, i, j), [&](auto k0, auto k1, auto k2, auto k3) {
          auto size = make_patch_size(degree);
          auto& b0 = mesh.getControlPoint(p, k0);
          auto& b1 = mesh.getControlPoint(q, k1);
          auto& b2 = mesh.getControlPoint(p, k2);
          auto& b3 = mesh.getControlPoint(p, k3);
          
          *(v0++) = (b2 - b0);
          *(v1++) = (b3 - b0);
          *(v2++) = (b1 - b0);
     });
     
     while(nSamples--) {
          auto temp = vectors;
          
          auto de_casteljau = [&](auto first, auto last) {
               while(first != last){
                    for(hpuint i = first; i < last; i++){
                         temp[i] = (1 - t) * temp[i] + t * temp[i+1];
                    }
                    last--;
               }
          };
          
          de_casteljau(0, degree - 1); //t0
          de_casteljau(degree, 2 * degree - 1); //t1
          de_casteljau(2 * degree, 3 * degree - 1); //t2
          
          auto t0 = vectors[0], t1 = vectors[degree], t2 = vectors[2 * degree];
          
          if(std::abs(glm::dot(t0, glm::cross(t1, t2))) > epsilon) { return false; }
          t += delta;
     }
     
     return true;
}

template<hpuint degree>
void make_bernstein_polynomials(const std::string& directory) {
     auto points = std::vector<Point3D>();
     points.reserve(make_patch_size(degree));
     auto plane = LinearBezierTriangleMesh<Space2D>({ Point2D(0.0, 0.0), Point2D(1.0, 0.0), Point2D(0.0, 1.0) }, { 0, 1, 2 });
     sample(plane, degree + 1u, [&](auto sample) { points.emplace_back(sample.x, sample.y, 0.0); });
     visit_bernstein_indices(degree, [&](auto i, auto j, auto k) {
          auto temp = points;
          temp[make_control_point_index(degree, i, j, k)].z = 1.0;
          auto mesh = BezierTriangleMesh<Space3D, degree>(temp);
          std::ostringstream path;
          path << directory << "/b" << i << j << k << ".ss" << degree << ".bz.3.hph";
          write(mesh, path.str());
     });
}

template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(hpuint n) {
     using Point = typename Space::POINT;

     auto length = make_patch_size(degree);

     return { { 1, Point(0) }, { length, n * length, 0 } };
}

template<class Space, hpuint degree, class Point>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(std::vector<Point> controlPoints) {
     auto length = make_patch_size(degree);
     auto indices = Tuples<hpindex>(length, controlPoints.size());//TODO: check correct number of control points

     std::iota(std::begin(indices), std::end(indices), 0);

     return { std::move(controlPoints), std::move(indices) };
}

template<class Space, hpuint degree, class Point>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(std::vector<Point> controlPoints, Tuples<hpindex> indices) { return { std::move(controlPoints), std::move(indices) }; }

template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(const std::string& mesh) { return format::hph::read<BezierTriangleMesh<Space, degree> >(mesh); }

template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> make_bezier_triangle_mesh(const std::experimental::filesystem::path& mesh) { return format::hph::read<BezierTriangleMesh<Space, degree> >(mesh); }

template<class Vertex>
auto make_bezier_triangle_mesh(const TriangleGraph<Vertex>& graph) {
     using Space = typename Vertex::SPACE;
     using Vector = typename Space::VECTOR;

     auto mesh = make_bezier_triangle_mesh<Space, 4u>(size(graph));

     auto set_boundary_point = [&](auto t, auto i, auto k, auto&& point) {
          auto u = make_neighbor_index(graph, t, i);
          auto j = make_neighbor_offset(graph, u, t);

          mesh.setControlPoint(t, i, k, u, j, point);
     };

     visit_diamonds(graph, [&](auto e, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
          auto t = make_triangle_index(e);
          auto i = make_edge_offset(e);

          set_boundary_point(t, i, 1, (hpreal(1.0) / hpreal(6.0)) * (hpreal(2.0) * vertex0.position + vertex1.position + hpreal(2.0) * vertex2.position + vertex3.position));
     });

     for(auto v : boost::irange(0u, graph.getNumberOfVertices())) {
          auto ring = make(make_ring_enumerator(graph, v));
          auto begin = std::begin(ring);
          auto end = std::end(ring);
          auto t = make_triangle_index(graph.getOutgoing(v));
          auto i = make_edge_offset(graph.getOutgoing(v));
          auto& center = graph.getVertex(v);
          auto fan = Indices();
          visit(make_spokes_enumerator(graph.getEdges(), graph.getOutgoing(v)), [&](auto e) {
               auto u = make_triangle_index(e);
               auto j = make_edge_offset(e);
               fan.push_back(u);
               fan.push_back(j);
          });
          auto valence = std::distance(begin, end);
          auto n = (valence << 1) - 4;

          auto set_interior_point = [&](auto t, auto i, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
               static constexpr hpindex o[3] = { 6, 7, 10 };

               auto point = (hpreal(1.0) / hpreal(24.0)) * (hpreal(10.0) * center.position + vertex0.position + hpreal(6.0) * vertex1.position + hpreal(6.0) * vertex2.position + vertex3.position);

               mesh.setControlPoint(t, o[i], point);
          };

          auto corner = center.position;
          if(valence == 6) {
               corner *= hpreal(6.0);
               for(auto& vertex : boost::make_iterator_range(begin, end)) corner += vertex.position;
               corner /= hpreal(12.0);

               auto make_boundary_point = [&](auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3, auto& vertex4) -> auto { return (hpreal(1.0) / hpreal(24.0)) * (hpreal(12.0) * center.position + vertex0.position + hpreal(3.0) * vertex1.position + hpreal(4.0) * vertex2.position + hpreal(3.0) * vertex3.position + vertex4.position); };

               set_boundary_point(fan[0], trit(fan[1]), 0, make_boundary_point(begin[4], begin[5], begin[0], begin[1], begin[2]));
               set_boundary_point(fan[2], trit(fan[3]), 0, make_boundary_point(begin[5], begin[0], begin[1], begin[2], begin[3]));
               set_boundary_point(fan[4], trit(fan[5]), 0, make_boundary_point(begin[0], begin[1], begin[2], begin[3], begin[4]));
               set_boundary_point(fan[6], trit(fan[7]), 0, make_boundary_point(begin[1], begin[2], begin[3], begin[4], begin[5]));
               set_boundary_point(fan[8], trit(fan[9]), 0, make_boundary_point(begin[2], begin[3], begin[4], begin[5], begin[0]));
               set_boundary_point(fan[10], trit(fan[11]), 0, make_boundary_point(begin[3], begin[4], begin[5], begin[0], begin[1]));
          } else {
               auto omega = hpreal(3.0) / hpreal(8.0) + std::cos(glm::two_pi<hpreal>() / valence) / hpreal(4.0);
               omega = (hpreal(5.0) / hpreal(8.0)) - omega * omega;
               omega = hpreal(3.0) * valence / (hpreal(8.0) * omega);
               corner *= omega;
               for(auto& vertex : boost::make_iterator_range(begin, end)) corner += vertex.position;
               corner /= (valence + omega);

               auto f = std::begin(fan);
               auto delta = glm::two_pi<hpreal>() / valence;
               for(auto middle = begin; middle != end; ++middle, f += 2) {
                    auto theta = hpreal(0.0);
                    auto tangent = Vector(0.0);
                    auto update_tangent = [&](auto& vertex) {
                         tangent += (hpreal(2.0) / valence) * std::cos(theta) * vertex.position;
                         theta += delta;
                    };
                    std::for_each(middle, end, update_tangent);
                    std::for_each(begin, middle, update_tangent);//TODO: instead of recomputing the tagent, simply rotate the first one
                    tangent = glm::normalize(tangent);
                    auto vector = mesh.getControlPoint(t, i, 1) - corner;
                    auto r = std::fmin(glm::length2(vector) / std::abs(glm::dot(tangent, vector)), glm::length(vector)) / hpreal(2.0);
                    set_boundary_point(f[0], trit(f[1]), 0, corner + r * tangent);
               }
          }
          mesh.setControlPoint(t, i, corner);
          visit_pairs(fan, [&](auto u, auto j) { mesh.setControlPoint(u, trit(j), t, trit(i)); });

          set_interior_point(fan[0], fan[1], begin[valence - 1], begin[0], begin[1], begin[2]);
          if(valence > 3) {
               auto middle = begin;
               visit_pairs(std::begin(fan) + 2, valence - 3, 2, [&](auto u, auto j) {
                    set_interior_point(u, j, middle[0], middle[1], middle[2], middle[3]);
                    ++middle;
               });
          }
          set_interior_point(fan[n], fan[n + 1], begin[valence - 3], begin[valence - 2], begin[valence - 1], begin[0]);
          set_interior_point(fan[n + 2], fan[n + 3], begin[valence - 2], begin[valence - 1], begin[0], begin[1]);
     }

     return mesh;
}

template<hpuint degree, class Iterator>
auto make_boundary(Iterator patch, trit i) {
     using T = typename std::remove_const<typename std::remove_reference<decltype(*patch)>::type>::type;

     auto boundary = std::vector<T>();

     boundary.reserve(degree - 1);
     visit_boundary<degree>(patch, i, make_back_inserter(boundary));
     return boundary;
}

constexpr hpindex make_control_point_index(hpuint degree, hpindex i0, hpindex i1, hpindex i2) { return make_patch_size(degree) - make_patch_size(degree - i2) + i1; }

inline hpindex make_control_point_index(hpuint degree, trit i, hpindex k) {
     switch(i) {
     case 0u: return k + 1u;
     case 1u: return (degree << 1) + ((k * ((degree << 1) - k - 1u)) >> 1);
     case 2u: return make_patch_size(degree) - 3u - ((k * (5u + k)) >> 1);
     }
}
     
inline hpindex make_control_point_index(hpuint degree, hpindex i) {
     auto delta = degree - 2u;
     auto end = degree - 3u;
     auto j = degree + 2u;//index

     while(i > end) {
          --delta;
          end += delta;
          j += delta + 3u;
     }
     return j + i - (end - delta + 1u);
}

template<class Space, hpuint degree, class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_control_polygon(const BezierTriangleMesh<Space, degree>& mesh, VertexFactory&& factory) { return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(mesh, 0, std::forward<VertexFactory>(factory)); }

constexpr hpuint make_control_polygon_size(hpuint degree) { return degree * degree; }

inline btm::DeltasEnumerator make_deltas_enumerator(hpuint degree) { return { degree }; };

inline btm::EdgeDiamondsEnumerator make_diamonds_enumerator(hpuint degree, trit i, trit j) { return { degree, i, j }; };

template<class Space, hpuint degree>
auto make_diamonds_enumerator(const BezierTriangleMesh<Space, degree>& mesh, hpindex p, trit i, hpindex q, trit j) { return transform(make_diamonds_enumerator(degree, i, j), [&](auto i0, auto i1, auto i2, auto i3) {
     auto& point0 = mesh.getControlPoint(p, i0);
     auto& point1 = mesh.getControlPoint(q, i1);
     auto& point2 = mesh.getControlPoint(p, i2);
     auto& point3 = mesh.getControlPoint(p, i3);

     return std::tie(point0, point1, point2, point3);
}); }

template<hpindex ring>
inline btm::VertexDiamondsEnumerator<ring> make_diamonds_enumerator(hpuint degree, const Triples<hpindex>& neighbors, hpindex p, trit i) { return { degree, make_spokes_walker(neighbors, p, i) }; }

template<hpindex ring, class Space, hpuint degree>
auto make_diamonds_enumerator(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors, hpindex p, trit i) { return transform(make_diamonds_enumerator<ring>(degree, neighbors, p, i), [&](auto p0, auto i0, auto p1, auto i1, auto p2, auto i2, auto p3, auto i3) {
     auto& point0 = mesh.getControlPoint(p0, i0);
     auto& point1 = mesh.getControlPoint(p1, i1);
     auto& point2 = mesh.getControlPoint(p2, i2);
     auto& point3 = mesh.getControlPoint(p3, i3);

     return std::tie(point0, point1, point2, point3);
}); }

template<class Space, hpuint degree>
Triples<hpindex> make_neighbors(const BezierTriangleMesh<Space, degree>& mesh) {
     static constexpr hpuint patchSize = make_patch_size(degree);

     auto& indices = mesh.getIndices();
     auto corners = Triples<hpindex>();

     corners.reserve(3 * size(mesh));
     visit(indices, [&](auto i) { corners.insert(std::end(corners), { i[0], i[degree], i[patchSize - 1] }); });

     return make_neighbors(corners);
}

inline btm::NablasEnumerator make_nablas_enumerator(hpuint degree) { return { degree }; };

constexpr hpuint make_patch_size(hpuint degree) { return (degree + 1) * (degree + 2) >> 1; }

template<hpindex ring>
btm::RingEnumerator<ring> make_ring_enumerator(hpuint degree, const Triples<hpindex>& neighbors, hpindex p, trit i) { return { { degree, { neighbors, p, i } } }; }

template<hpindex ring, class Space, hpuint degree>
auto make_ring_enumerator(const BezierTriangleMesh<Space, degree>& mesh, const Triples<hpindex>& neighbors, hpindex p, trit i) { return transform(make_ring_enumerator<ring>(degree, neighbors, p, i), [&](auto p, auto i) { return mesh.getControlPoint(p, i); }); }

template<hpindex i>
btm::RowEnumerator<i> make_row_enumerator(hpuint degree, hpindex row) { return { degree, row }; }

template<class Space, hpuint degree, class Vertex, class VertexFactory, typename>
TriangleMesh<Vertex> make_triangle_mesh(const BezierTriangleMesh<Space, degree>& mesh, hpuint nSubdivisions, VertexFactory&& factory) {
     using Point = typename Space::POINT;
     static_assert(degree > 0u, "Constant spline surfaces cannot be converted into triangle meshes.");
     static_assert(std::is_base_of<Vertex, decltype(factory(Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");

     auto do_make_triangle_mesh = [&](const BezierTriangleMesh<Space, degree>& mesh) -> auto {
          auto vertices = std::vector<Vertex>();
          vertices.reserve(mesh.getControlPoints().size());
          for(auto& point : mesh.getControlPoints()) vertices.push_back(factory(point));

          auto indices = Triples<hpindex>();
          auto inserter = make_back_inserter(indices);

          indices.reserve(3 * make_control_polygon_size(degree) * size(mesh));
          visit(mesh.getIndices(), [&](auto patch) {
               auto lambda = [&](auto i0, auto i1, auto i2) { return std::tie(patch[i0], patch[i1], patch[i2]); };

               visit(transform(make_deltas_enumerator(degree), lambda), inserter);
               visit(transform(make_nablas_enumerator(degree), lambda), inserter);
          });

          return make_triangle_mesh(std::move(vertices), std::move(indices));
     };

     if(nSubdivisions > 0) return do_make_triangle_mesh(subdivide(mesh, nSubdivisions));
     else return do_make_triangle_mesh(mesh);
}

template<class Space, hpuint degree, class Visitor>
void sample(const BezierTriangleMesh<Space, degree>& mesh, hpuint nSamples, Visitor&& visit) {
     using Point = typename Space::POINT;

     //TODO; skip multiple computations on common edges and eliminate common points in array
     auto matrix = make_de_casteljau_matrix(degree, nSamples);

     happah::visit(mesh, [&](auto patch) { visit(std::inner_product(std::begin(matrix), std::end(matrix), patch, Point(0))); });
}

template<class Space, hpuint degree>
auto size(const BezierTriangleMesh<Space, degree>& mesh) { return mesh.getNumberOfPatches(); }

template<hpuint degree>
BezierTriangleMesh<Space4D, degree> smooth(BezierTriangleMesh<Space4D, degree> mesh, const ProjectiveStructure& structure, hpreal epsilon) {
     static_assert(degree > 4u, "The first two rings of control points surrounding the corners of the patches are assumed to be disjoint.");

     auto mesh1 = make_bezier_triangle_mesh<Space4D, degree>(size(mesh));
     auto& neighbors = structure.getNeighbors();
     auto& transitions = structure.getTransitions();

     auto make_coefficients_1 = [&](auto p, auto i, auto valence, auto& center) {
          auto coefficients = std::vector<double>(valence * 7, 0.0);
          auto a0 = std::begin(coefficients) - 1;
          auto a1 = a0 + valence;
          auto a2 = a1 + valence;
          auto b0 = a2 + valence;
          auto b1 = b0 + valence;
          auto b2 = b1 + valence;
          auto b3 = b2 + valence;
          auto e = make_ring_enumerator<1>(mesh, neighbors, p, i);
          auto f = transform(make_spokes_enumerator(neighbors, p, i), [&](auto p, auto i) { return std::begin(transitions) + 3 * (3 * p + i); });

          auto push_back = [&](auto t0, auto t1, auto t2) {
               auto point = *e - hpreal(t0) * center;

               (++a0)[0] = t0;
               (++a1)[0] = t1;
               (++a2)[0] = t2;
               (++b0)[0] = point.x;
               (++b1)[0] = point.y;
               (++b2)[0] = point.z;
               (++b3)[0] = point.w;
               ++e;
          };

          push_back(0.0, 1.0, 0.0);
          push_back(0.0, 0.0, 1.0);
          ++f;

          while(e) {
               auto l = *f;
               auto t0 = l[2] + l[0] * a0[0] + l[1] * (a0 - 1)[0];
               auto t1 = l[0] * a1[0] + l[1] * (a1 - 1)[0];
               auto t2 = l[0] * a2[0] + l[1] * (a2 - 1)[0];

               push_back(t0, t1, t2);
               ++f;
          }

          return coefficients;
     };

     auto set_ring_1 = [&](auto p, auto i, auto valence, auto& center, auto& coefficients) {
          auto a0 = std::begin(coefficients) + 1;
          auto a1 = a0 + valence;
          auto a2 = a1 + valence;
          auto A = Eigen::Map<Eigen::MatrixX2d>(coefficients.data() + valence, valence, 2);
          auto b0 = Eigen::Map<Eigen::VectorXd>(coefficients.data() + 3 * valence, valence);
          auto b1 = Eigen::Map<Eigen::VectorXd>(coefficients.data() + 4 * valence, valence);
          auto b2 = Eigen::Map<Eigen::VectorXd>(coefficients.data() + 5 * valence, valence);
          auto b3 = Eigen::Map<Eigen::VectorXd>(coefficients.data() + 6 * valence, valence);
          auto temp = A.colPivHouseholderQr();
          auto x0 = temp.solve(b0);
          auto x1 = temp.solve(b1);
          auto x2 = temp.solve(b2);
          auto x3 = temp.solve(b3);
          auto q1 = Point4D(x0[0], x1[0], x2[0], x3[0]);
          auto q2 = Point4D(x0[1], x1[1], x2[1], x3[1]);
          auto e = make_spokes_enumerator(neighbors, p, i);

          auto set_boundary_point = [&](auto point) {
               auto q = std::get<0>(*e);
               auto j = std::get<1>(*e);
               auto r = neighbors(q, j);
               auto k = make_offset(neighbors, r, q);

               mesh1.setControlPoint(q, j, 0, r, k, point);
               ++e;
          };

          assert(q1.w > epsilon);
          assert(q2.w > epsilon);
          set_boundary_point(q1);
          set_boundary_point(q2);

          while(e) set_boundary_point(hpreal((++a0)[0]) * center + hpreal((++a1)[0]) * q1 + hpreal((++a2)[0]) * q2);
     };

     auto make_coefficients_2 = [&](auto p, auto i, auto valence) {
          auto l0 = transitions[3 * (3 * p + i)];
          //auto nRows = valence << 1;
          auto nRows = (valence << 1) - 1;
          auto coefficients = std::vector<double>(nRows * (valence + 4), 0.0);
          auto b0 = std::begin(coefficients) - 1;
          auto b1 = b0 + nRows;
          auto b2 = b1 + nRows;
          auto b3 = b2 + nRows;
          auto A = b3 + (nRows + 1);
          auto r0 = hpreal(0), r1 = hpreal(0), r2 = hpreal(0), r3 = hpreal(0);
          auto e = make_diamonds_enumerator<2>(mesh, neighbors, p, i);
          auto f = transform(make_spokes_enumerator(neighbors, p, i), [&](auto p, auto i) { return std::begin(transitions) + 3 * (3 * p + i); });
          auto n = 0;

          auto push_back_0 = [&](auto& point) {
               (++b0)[0] = point.x;
               (++b1)[0] = point.y;
               (++b2)[0] = point.z;
               (++b3)[0] = point.w;
               A[n * nRows] = 1.0;
               ++A;
          };

          auto push_back_1 = [&](auto& point0, auto& point1, auto& point2, auto& point3) {
               auto l = *f;
               auto a = A + 1;

               push_back_0(point2);
               r0 = l[2] * point0.x + l[1] * r0;
               r1 = l[2] * point0.y + l[1] * r1;
               r2 = l[2] * point0.z + l[1] * r2;
               r3 = l[2] * point0.w + l[1] * r3;
               (++b0)[0] = point3.x - r0;
               (++b1)[0] = point3.y - r1;
               (++b2)[0] = point3.z - r2;
               (++b3)[0] = point3.w - r3;
               repeat(n, [&]() {
                    a[0] = l[1] * (a - 2)[0];
                    a += nRows;
               });
               a[0] = l[0];
               ++A;
          };

          /*auto push_back_2 = [&](auto point0, auto& point1, auto point2, auto point3) {
               auto l = *f;
               auto a = A;

               (++b0)[0] = point2.x + (l[2] * point0.x + l[1] * r0) / l[0];
               (++b1)[0] = point2.y + (l[2] * point0.y + l[1] * r1) / l[0];
               (++b2)[0] = point2.z + (l[2] * point0.z + l[1] * r2) / l[0];
               (++b3)[0] = point2.w + (l[2] * point0.w + l[1] * r3) / l[0];
               repeat(n, [&]() {
                    a[0] = (-l[1] / l[0]) * (a - 1)[0];
                    a += nRows;
               });
               *A += hpreal(1) / l[0];
          };

          auto push_back_3 = [&](auto point0, auto& point1, auto point2, auto point3) {
               auto l = *f;
               auto a = A;

               (++b0)[0] = l[2] * point0.x + l[1] * r0;
               (++b1)[0] = l[2] * point0.y + l[1] * r1;
               (++b2)[0] = l[2] * point0.z + l[1] * r2;
               (++b3)[0] = l[2] * point0.w + l[1] * r3;
               repeat(n, [&]() {
                    a[0] = -l[1] * (a - 1)[0];
                    a += nRows;
               });
               *A += hpreal(1);
          };*/

          push_back_0(std::get<3>(*e));
          ++e;
          ++f;
          ++n;

          while(e) {
               apply(push_back_1, *e);
               ++e;
               ++f;
               ++n;
          }
          //if(std::abs(l0) > epsilon) apply(push_back_2, *e);
          //else apply(push_back_3, *e);

          return coefficients;
     };

     auto set_ring_2 = [&](auto p, auto i, auto valence, auto& coefficients) {
          auto l0 = transitions[3 * (3 * p + i)];
          //auto nRows = valence << 1;
          auto nRows = (valence << 1) - 1;
          auto A = Eigen::Map<Eigen::MatrixXd>(coefficients.data() + (nRows << 2), nRows, valence);
          auto b0 = Eigen::Map<Eigen::VectorXd>(coefficients.data(), nRows);
          auto b1 = Eigen::Map<Eigen::VectorXd>(coefficients.data() + nRows, nRows);
          auto b2 = Eigen::Map<Eigen::VectorXd>(coefficients.data() + (nRows << 1), nRows);
          auto b3 = Eigen::Map<Eigen::VectorXd>(coefficients.data() + 3 * nRows, nRows);
          auto temp = A.colPivHouseholderQr();
          auto x0 = temp.solve(b0);
          auto x1 = temp.solve(b1);
          auto x2 = temp.solve(b2);
          auto x3 = temp.solve(b3);
          auto e1 = make_ring_enumerator<1>(mesh1, neighbors, p, i);
          auto e2 = make_ring_enumerator<2>(degree, neighbors, p, i);
          auto f = make_spokes_enumerator(neighbors, p, i);
          auto n = 0;
          auto point1 = Point4D(x0[0], x1[0], x2[0], x3[0]);

          auto set_boundary_point = [&](auto q, auto j, auto& point) {
               auto r = neighbors(q, j);
               auto k = make_offset(neighbors, r, q);

               mesh1.setControlPoint(q, j, 1, r, k, point);
          };

          auto set_interior_point = [&](auto& point) {
               auto q = 0u, j = 0u;

               std::tie(q, j) = *(++e2);
               mesh1.setControlPoint(q, j, point);
          };

          assert(x3[0] > epsilon);
          set_interior_point(point1);
          ++e1;
          ++e2;
          ++f;
          ++n;

          while(e1) {
               auto point0 = *e1;
               auto point2 = Point4D(x0[n], x1[n], x2[n], x3[n]);
               auto q = std::get<0>(*f);
               auto j = std::get<1>(*f);
               auto l = std::begin(transitions) + 3 * (3 * q + j);
               auto point3 = l[0] * point2 + l[1] * point1 + l[2] * point0;

               assert(x3[n] > epsilon);
               set_boundary_point(q, j, point2);
               set_interior_point(point3);
               point1 = point3;
               ++e1;
               ++e2;
               ++f;
               ++n;
          }

          if(std::abs(l0) > epsilon) {
               auto l = std::begin(transitions) + 3 * (3 * p + i);
               auto point3 = Point4D(x0[0], x1[0], x2[0], x3[0]);
               auto point0 = *e1;
               auto point2 = (point3 - l[1] * point1 - l[2] * point0) / l[0];

               set_boundary_point(p, i, point2);
          } else set_boundary_point(p, i, mesh.getControlPoint(p, i, 1));
     };

     //TODO: try to optimize both rings at once
     //TODO: use linear least squares with non-negative constraint from octave to ensure weights are positive

     visit(make_vertices_enumerator(neighbors), [&](auto p, auto i) {
          auto valence = size(make_spokes_enumerator(neighbors, p, i));
          auto& center = mesh.getControlPoint(p, i);
          auto coefficients1 = make_coefficients_1(p, i, valence, center);
          auto coefficients2 = make_coefficients_2(p, i, valence);

          mesh1.setControlPoint(p, i, center);
          visit(make_spokes_enumerator(neighbors, p, i), [&](auto q, auto j) { mesh1.setControlPoint(q, j, p, i); });
          set_ring_1(p, i, valence, center, coefficients1);
          set_ring_2(p, i, valence, coefficients2);
     });

     assert(degree == 5);//TODO: update edges and copy interior points for degrees > 5

     if(degree == 5) visit_edges(neighbors, [&](auto p, auto i) {
          static constexpr hpindex o0[3] = { 2, 14, 15 };
          static constexpr hpindex o1[3] = { 8, 13, 12 };
          static constexpr hpindex o2[3] = { 3, 17, 11 };

          auto q = neighbors(p, i);
          auto j = make_offset(neighbors, q, p);
          auto& point0 = mesh1.getControlPoint(p, o0[i]);
          auto& point1 = mesh.getControlPoint(q, o1[j]);
          auto& point2 = mesh1.getControlPoint(p, o2[i]);
          auto& point3 = mesh.getControlPoint(p, o1[i]);
          auto l = std::begin(transitions) + 3 * (3 * p + i);
          auto x1 = (point1 + l[1] * (point3 - l[0] * point2 - l[2] * point0)) / (hpreal(1) + l[1] * l[1]);
          auto x3 = l[0] * point2 + l[1] * x1 + l[2] * point0;

          mesh1.setControlPoint(p, o1[i], x3);
          mesh1.setControlPoint(q, o1[j], x1);
     });

     visit_edges(neighbors, [&](auto p, auto i) {
          auto q = neighbors(p, i);
          auto j = make_offset(neighbors, q, p);
          auto l = std::begin(transitions) + 3 * (3 * p + i);

          visit(make_diamonds_enumerator(mesh1, p, i, q, j), [&](auto& point0, auto& point1, auto& point2, auto& point3) { assert(glm::length(point3 - (l[0] * point2 + l[1] * point1 + l[2] * point0)) < epsilon); });
     });

     return mesh1;
}

template<class Space, hpuint degree>
BezierTriangleMesh<Space, degree> subdivide(const BezierTriangleMesh<Space, degree>& mesh, hpuint nSubdivisions) {
     using Point = typename Space::POINT;

     if(nSubdivisions == 0) return mesh;

     auto nPatches = size(mesh);

     std::vector<BezierTriangleSubdivider<Space, degree> > subdividers;
     subdividers.reserve(nPatches);
     visit(mesh, [&](auto patch) { subdividers.emplace_back(patch); });

     std::vector<Point> points;
     auto indices = Tuples<hpindex>(make_patch_size(degree));

     points.reserve(nPatches * make_patch_size(degree));
     indices.reserve(3 * nPatches * make_control_polygon_size(degree));

     for(auto& subdivider : subdividers) {
          auto subdivided = subdivider.subdivide(nSubdivisions);
          auto offset = points.size();
          for(auto& i : std::get<1>(subdivided)) i += offset;
          std::move(std::begin(std::get<0>(subdivided)), std::end(std::get<0>(subdivided)), std::back_inserter(points));
          std::move(std::begin(std::get<1>(subdivided)), std::end(std::get<1>(subdivided)), std::back_inserter(indices));
     }

     return { std::move(points), std::move(indices) };
}

template<class Space, hpuint degree, class Visitor>
void visit(const BezierTriangleMesh<Space, degree>& mesh, Visitor&& visit) { happah::visit(deindex(mesh), std::forward<Visitor>(visit)); }

template<class Visitor>
void visit_bernstein_indices(hpuint degree, Visitor&& visit) {
     auto i = degree;
     auto k = 0u;
     while(i > 0u) {
          for(auto j : boost::irange(0u, i + 1u)) visit(degree - j - k, j, k);
          --i;
          ++k;
     }
     visit(i, 0u, k);
}

template<hpuint degree, class Iterator, class Visitor>
void visit_boundary(Iterator patch, trit i, Visitor&& visit) {
     if(i == 0u) for(auto end = patch + (degree - 1u); patch != end; ) visit(*(++patch));
     else if(i == 1u) {
          auto delta = degree;
          for(patch += degree << 1; delta > 1u; patch += --delta) visit(*patch);
     } else {
          auto delta = 2u;
          for(patch += make_patch_size(degree) - 3u; delta <= degree; patch -= ++delta) visit(*patch);
     }
}

/*template<class Iterator, class Visitor>
void visit_deltas(hpuint degree, Iterator patch, Visitor&& visit) {
     auto bottom = patch;
     auto top = patch + (degree + 1u);
     auto delta = degree;
     while(delta > 0u) {
          for(auto end = bottom + delta; bottom != end; ++bottom, ++top) visit(bottom[0], bottom[1], top[0]);
          --delta;
          ++bottom;
     }
}*/

template<hpuint degree, class Iterator, class Visitor>
void visit_ends(Iterator patch, trit i, Visitor&& visit) {
     static constexpr hpuint o[3] = { 0u, degree, make_patch_size(degree) - 1u };

     if(i == 0u) visit(patch[o[0]], patch[o[1]]);
     else if(i == 1u) visit(patch[o[1]], patch[o[2]]);
     else visit(patch[o[2]], patch[o[0]]);
}

template<hpuint degree, class Iterator, class Visitor>
void visit_interior(Iterator patch, Visitor&& visit) {
     static_assert(degree > 2, "There is no interior in a constant, linear, or quadratic.");
     patch += degree + 2u;
     auto delta = degree - 2u;
     while(delta > 0u) {
          for(auto end = patch + delta; patch != end; ++patch) visit(*patch);
          patch += 2u;
          --delta;
     }
}

/*template<class Iterator, class Visitor>
void visit_nablas(hpuint degree, Iterator patch, Visitor&& visit) {
     auto bottom = patch + 1;
     auto top = patch + (degree + 1u);
     auto delta = degree; // + 1 so that visit is a noop if degree is zero
     while(delta > 1u) {
          for(auto end = bottom + (delta - 1u); bottom != end; ++bottom, ++top) visit(top[1], top[0], bottom[0]);
          --delta;
          bottom += 2;
          ++top;
     }
}*/

template<class Visitor>
void visit_subring(btm::RingEnumerator<1> e, hpindex p, Visitor&& visit) {
     while(std::get<0>(*e) != p) {
          apply(visit, *e);
          ++e;
     }
     apply(visit, *e);
}

template<hpuint degree>
BezierTriangleMesh<Space4D, degree> weigh(const BezierTriangleMesh<Space3D, degree>& mesh) { return embed<Space4D>(mesh, [](const Point3D& point) { return Point4D(point.x, point.y, point.z, 1.0); }); }

template<hpuint degree>
BezierTriangleMesh<Space4D, degree> weigh(BezierTriangleMesh<Space4D, degree> mesh, const ProjectiveStructure& structure) {
     using Vector = Eigen::Matrix<hpreal, Eigen::Dynamic, 1>;

     static_assert(degree > 4u, "The first two rings of control points surrounding the corners of the patches are assumed to be disjoint.");

     auto& neighbors = structure.getNeighbors();
     auto& transitions = structure.getTransitions();

     auto make_weights_1 = [&](auto p, auto i, auto valence, auto& center) {
          auto a = std::vector<Eigen::Triplet<hpreal> >();
          auto b = Vector(Vector::Zero(valence << 2));
          auto n = hpuint(0);
          auto e = make_diamonds_enumerator<1>(mesh, neighbors, p, i);
          auto f = transform(make_spokes_enumerator(neighbors, p, i), [&](auto p, auto i) { return std::begin(transitions) + 3 * (3 * p + i); });

          auto push_back = [&](auto& point1, auto n1, auto& point2, auto n2, auto& point3, auto n3) {
               auto l = *f;

               a.emplace_back((n << 2) + 0, n3, point3.x);
               a.emplace_back((n << 2) + 1, n3, point3.y);
               a.emplace_back((n << 2) + 2, n3, point3.z);
               a.emplace_back((n << 2) + 3, n3, point3.w);
               a.emplace_back((n << 2) + 0, n2, -l[0] * point2.x);
               a.emplace_back((n << 2) + 1, n2, -l[0] * point2.y);
               a.emplace_back((n << 2) + 2, n2, -l[0] * point2.z);
               a.emplace_back((n << 2) + 3, n2, -l[0] * point2.w);
               a.emplace_back((n << 2) + 0, n1, -l[1] * point1.x);
               a.emplace_back((n << 2) + 1, n1, -l[1] * point1.y);
               a.emplace_back((n << 2) + 2, n1, -l[1] * point1.z);
               a.emplace_back((n << 2) + 3, n1, -l[1] * point1.w);
               b[(n << 2) + 0] = l[2] * center.x;
               b[(n << 2) + 1] = l[2] * center.y;
               b[(n << 2) + 2] = l[2] * center.z;
               b[(n << 2) + 3] = l[2] * center.w;
               ++f;
               ++n;
          };

          auto push_back_0 = [&](auto& point0, auto& point1, auto& point2, auto& point3) { push_back(point1, valence - 1, point2, 0, point3, 1); };

          auto push_back_1 = [&](auto& point0, auto& point1, auto& point2, auto& point3) { push_back(point1, n - 1, point2, n, point3, n + 1); };

          auto push_back_2 = [&](auto& point0, auto& point1, auto& point2, auto& point3) { push_back(point1, valence - 2, point2, valence - 1, point3, 0); };

          apply(push_back_0, *e);
          repeat(valence - 2, [&]() { apply(push_back_1, *(++e)); });
          apply(push_back_2, *(++e));

          return lsq::solve(make_sparse_matrix(valence << 2, valence, a), b);
     };

     auto make_weights_2 = [&](auto p, auto i, auto valence) {
          auto a = std::vector<Eigen::Triplet<hpreal> >();
          auto b = Vector(Vector::Zero(6 * valence));
          auto n = hpuint(0);
          auto e = make_diamonds_enumerator<2>(mesh, neighbors, p, i);
          auto f = transform(make_spokes_enumerator(neighbors, p, i), [&](auto p, auto i) { return std::begin(transitions) + 3 * (3 * p + i); });

          auto push_back = [&](auto& point0, auto& point1, auto n1, auto& point2, auto n2, auto& point3, auto n3) {
               auto l = *f;

               a.emplace_back((6 * n) + 0, n3, point3.x);
               a.emplace_back((6 * n) + 1, n3, point3.y);
               a.emplace_back((6 * n) + 2, n3, point3.z);
               a.emplace_back((6 * n) + 3, n3, point3.w);
               a.emplace_back((6 * n) + 0, n2, -l[0] * point2.x);
               a.emplace_back((6 * n) + 1, n2, -l[0] * point2.y);
               a.emplace_back((6 * n) + 2, n2, -l[0] * point2.z);
               a.emplace_back((6 * n) + 3, n2, -l[0] * point2.w);
               a.emplace_back((6 * n) + 0, n1, -l[1] * point1.x);
               a.emplace_back((6 * n) + 1, n1, -l[1] * point1.y);
               a.emplace_back((6 * n) + 2, n1, -l[1] * point1.z);
               a.emplace_back((6 * n) + 3, n1, -l[1] * point1.w);
               a.emplace_back((6 * n) + 4, n2, 1);
               a.emplace_back((6 * n) + 5, n3, 1);
               b[(6 * n) + 0] = l[2] * point0.x;
               b[(6 * n) + 1] = l[2] * point0.y;
               b[(6 * n) + 2] = l[2] * point0.z;
               b[(6 * n) + 3] = l[2] * point0.w;
               b[(6 * n) + 4] = 1;
               b[(6 * n) + 5] = 1;
               ++f;
               ++n;
          };

          auto push_back_0 = [&](auto& point0, auto& point1, auto& point2, auto& point3) { push_back(point0, point1, (valence << 1) - 1, point2, 0, point3, 1); };

          auto push_back_1 = [&](auto& point0, auto& point1, auto& point2, auto& point3) { push_back(point0, point1, (n << 1) - 1, point2, n << 1, point3, (n << 1) + 1); };

          apply(push_back_0, *e);
          while(++e) apply(push_back_1, *e);

          return lsq::solve(make_sparse_matrix(6 * valence, valence << 1, a), b);
     };

     visit(make_vertices_enumerator(neighbors), [&](auto p, auto i) {
          auto valence = size(make_spokes_enumerator(neighbors, p, i));
          auto& center = mesh.getControlPoint(p, i);
          auto weights1 = make_weights_1(p, i, valence, center);
          auto weights2 = make_weights_2(p, i, valence);
          auto k1 = hpuint(-1);
          auto k2 = hpuint(-1);

          visit(make_ring_enumerator<1>(degree, neighbors, p, i), [&](auto p, auto i) { mesh.getControlPoint(p, i) *= weights1[++k1]; assert(weights1[k1] > hpreal(0)); });
          visit(make_ring_enumerator<2>(degree, neighbors, p, i), [&](auto p, auto i) { mesh.getControlPoint(p, i) *= weights2[++k2]; assert(weights2[k2] > hpreal(0)); });
     });

     assert(degree == 5);//TODO: update weights on inner edge diamonds for degree > 5

     if(degree == 5) visit_edges(neighbors, [&](auto p, auto i) {
          static constexpr hpindex o0[3] = { 2, 14, 15 };
          static constexpr hpindex o1[3] = { 8, 13, 12 };
          static constexpr hpindex o2[3] = { 3, 17, 11 };

          auto q = neighbors(p, i);
          auto j = make_offset(neighbors, q, p);
          auto& point0 = mesh.getControlPoint(p, o0[i]);
          auto& point1 = mesh.getControlPoint(q, o1[j]);
          auto& point2 = mesh.getControlPoint(p, o2[i]);
          auto& point3 = mesh.getControlPoint(p, o1[i]);
          auto l = std::begin(transitions) + 3 * (3 * p + i);
          auto m = std::begin(transitions) + 3 * (3 * q + j);
          //alternative weights 1
          auto g0 = l[1] * glm::dot(point1, point3) / glm::dot(point3, point3);
          auto g1 = glm::dot(l[0] * point2 + l[2] * point0, point3) / glm::dot(point3, point3);
          auto l0 = m[1] * glm::dot(point1, point3) / glm::dot(point1, point1);
          auto l1 = glm::dot(m[0] * point0 + m[2] * point2, point1) / glm::dot(point1, point1);
          auto ws = glm::inverse(hpmat2x2(-g0, 1, 1, -l0)) * Point2D(g1, l1);
          point1 *= ws.x;
          point3 *= ws.y;
          //alternative weights 2
          /*auto w3 = glm::dot((l[0] * point2 + l[1] * point1 + l[2] * point0), point3) / glm::length2(point3);
          auto w1 = glm::dot((m[0] * point0 + m[1] * point3 + m[2] * point2), point1) / glm::length2(point1);
          point1 *= w1;
          point3 *= w3;*/
     });

     return mesh;
}

//WORKSPACE

namespace mdz {

//Returns quadratics of the form $$ax_jx_k+bx_j'+c=0$$, where the quadratic coefficients are stored in the first vector, the linear coefficients in the second, and the constant coefficients in the third.  The first integer (the 'i') in each entry of the three vectors identifies the constraint.
std::tuple<std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const Triples<hpindex>& neighbors);

template<class Space, hpuint degree>
std::tuple<std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const BezierTriangleMesh<Space, degree>& mesh) {
     auto neighbors = make_neighbors(mesh);
     return make_constraints(neighbors);
}

//Returns an objective of the form |Ax - b| that is to be minimized.  There are 3 * 9 * nPatches variables.
template<hpuint degree>
std::tuple<std::vector<hpijr>, std::vector<hpir> > make_objective(const BezierTriangleMesh<Space3D, degree>& mesh) {
     using Point = Point3D;

     auto patches = deindex(mesh);
     auto neighbors = make_neighbors(mesh);
     auto nEdges = 3 * size(mesh) / 2;
     auto irs = std::vector<hpir>();
     auto ijrs = std::vector<hpijr>();
     auto row = -1;

     // indexing of rho:
     //   x0 x1 x2
     //   x3 x4 x5
     //   x6 x7 x8

     auto insert = [&](auto& source, auto& target, auto offset) {
          irs.emplace_back(++row, target.x);
          ijrs.emplace_back(row, offset + 0, source.x);
          ijrs.emplace_back(row, offset + 1, source.y);
          ijrs.emplace_back(row, offset + 2, source.z);
          irs.emplace_back(++row, target.y);
          ijrs.emplace_back(row, offset + 3, source.x);
          ijrs.emplace_back(row, offset + 4, source.y);
          ijrs.emplace_back(row, offset + 5, source.z);
          irs.emplace_back(++row, target.z);
          ijrs.emplace_back(row, offset + 6, source.x);
          ijrs.emplace_back(row, offset + 7, source.y);
          ijrs.emplace_back(row, offset + 8, source.z);
     };

     irs.reserve(2 * 4 * 3 * nEdges);
     ijrs.reserve(2 * 4 * 9 * nEdges);

     visit_edges(neighbors, [&](auto p, auto i) {
          static constexpr hpuint o[3] = { 1u, 2u, 0u };
          static const trit o1[3] = { TRIT1, TRIT0, TRIT1 };

          auto q = neighbors(p, i);
          auto j = make_offset(neighbors, q, p);
          auto op = 27 * p + 9 * i;
          auto oq = 27 * q + 9 * j;
          auto b = std::array<Point, 3>();
          auto c = std::array<Point, 3>();
          auto& b0 = mesh.getControlPoint(q, j);
          auto& b1 = b[1];
          auto& b2 = b[0];
          auto& b3 = b[2];
          auto& c0 = c[1];
          auto& c1 = mesh.getControlPoint(p, i);
          auto& c2 = c[2];
          auto& c3 = c[0];
          auto tb = b.data() - 1;
          auto tc = c.data() - 1;

          //TODO: simplify this method using ring enumerator and write directly into b0, b1, b2, b3
          visit_subring(make_ring_enumerator(degree, neighbors, p, o[i]), neighbors(q, o1[j]), [&](auto r, auto k) { *(++tb) = patches(r)[k]; });
          visit_subring(make_ring_enumerator(degree, neighbors, q, o[j]), neighbors(p, o1[i]), [&](auto r, auto k) { *(++tc) = patches(r)[k]; });

          insert(c0, b0, op);
          insert(c1, b1, op);
          insert(c2, b2, op);
          insert(c3, b3, op);
          insert(b0, c0, oq);
          insert(b1, c1, oq);
          insert(b2, c2, oq);
          insert(b3, c3, oq);
     });

     return std::make_tuple(std::move(ijrs), std::move(irs));
}

template<hpuint degree>
std::vector<hpreal> make_transitions(const BezierTriangleMesh<Space3D, degree>& mesh, const std::vector<hpreal>& solution) {
     static_assert(degree > 4u, "The first two rings of control points surrounding the corners of the patches are assumed to be disjoint.");

     auto transitions = std::vector<hpreal>(9 * size(mesh));
     auto patches = deindex(mesh);
     auto neighbors = make_neighbors(mesh);

     auto insert = [&](auto b0, auto b1, auto b2, auto b3, auto p, auto i) {
          auto transition = glm::inverse(hpmat3x3(b0, b1, b2)) * b3;
          auto t = std::begin(transitions) + (9 * p  + 3 * i);
          t[0] = transition.x;
          t[1] = transition.y;
          t[2] = transition.z;
     };

     auto make_point = [&](auto p, auto i) {
          auto& corner = mesh.getControlPoint(p, i);
          auto rho = solution.data() + (27 * p + 9 * i);

          return Point3D(
               corner.x * rho[0] + corner.y * rho[1] + corner.z * rho[2],
               corner.x * rho[3] + corner.y * rho[4] + corner.z * rho[5],
               corner.x * rho[6] + corner.y * rho[7] + corner.z * rho[8]
          );
     };

     visit_edges(neighbors, [&](auto p, auto i) {
          static constexpr hpuint o0[3] = { 1u, 2u, 0u };
          static constexpr hpuint o1[3] = { 2u, 0u, 1u };
          
          auto q = neighbors(p, i);
          auto j = make_offset(neighbors, q, p);
          auto r = neighbors(p, o0[i]);
          auto k = make_offset(neighbors, r, p);
          auto& b0 = mesh.getControlPoint(p, trit(o0[i]));
          auto b1 = make_point(p, i);
          auto b2 = make_point(r, k);
          auto b3 = make_point(q, o1[j]);

          insert(b0, b1, b2, b3, p, i);
          insert(b1, b0, b3, b2, q, j);
     });

     return transitions;
}

}//namespace mdz

namespace phm {

//Returns cubics of the form $$ax_jx_kx_l+bx_j'x_k'+cx_j''+d=0$$, where the cubic coefficients are stored in the first vector, the quadratic coefficients in the second, the linear coefficients in the third, and the constant coefficients in the fourth.  The first integer (the 'i') in each entry of the four vectors identifies the constraint.
std::tuple<std::vector<hpijklr>, std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const Triples<hpindex>& neighbors);

template<class Space, hpuint degree>
std::tuple<std::vector<hpijklr>, std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const BezierTriangleMesh<Space, degree>& mesh) {
     auto neighbors = make_neighbors(mesh);
     return make_constraints(neighbors);
}

//Returns an objective of the form |Ax - b| that is to be minimized.  There are 3 * 12 * nPatches variables.
template<hpuint degree>
std::tuple<std::vector<hpijr>, std::vector<hpir> > make_objective(const BezierTriangleMesh<Space3D, degree>& mesh) {
     using Point = Point3D;
     auto patches = deindex(mesh);
     auto neighbors = make_neighbors(mesh);
     auto nEdges = 3 * size(mesh) / 2;
     auto irs = std::vector<hpir>();
     auto ijrs = std::vector<hpijr>();
     auto row = -1;

     // indexing of rho:
     //   x0 x1 x2
     //   x3 x4 x5
     //   x6 x7 x8
     // indexing of lambda:
     //   1 0 x9
     //   0 0 x10
     //   0 1 x11

     auto do_column = [&](auto offset, auto x, auto y, auto z) {
          ijrs.emplace_back(++row, offset + 0, 1.0);
          irs.emplace_back(row, x);
          ijrs.emplace_back(++row, offset + 1, 1.0);
          irs.emplace_back(row, y);
          ijrs.emplace_back(++row, offset + 2, 1.0);
          irs.emplace_back(row, z);
     };

     auto do_row = [&](auto offset, auto x, auto y, auto z, auto a) {
          ijrs.emplace_back(++row, offset + 0, x);
          ijrs.emplace_back(row, offset + 1, y);
          ijrs.emplace_back(row, offset + 2, z);
          irs.emplace_back(row, a);
     };

     auto insert = [&](auto offset, auto b2, auto b3, auto c2, auto c3) {
          // |rho - id|
          ijrs.emplace_back(++row, offset + 0, 1.0);
          irs.emplace_back(row, 1.0);
          ijrs.emplace_back(++row, offset + 1, 1.0);
          ijrs.emplace_back(++row, offset + 2, 1.0);
          ijrs.emplace_back(++row, offset + 3, 1.0);
          ijrs.emplace_back(++row, offset + 4, 1.0);
          irs.emplace_back(row, 1.0);
          ijrs.emplace_back(++row, offset + 5, 1.0);
          ijrs.emplace_back(++row, offset + 6, 1.0);
          ijrs.emplace_back(++row, offset + 7, 1.0);
          ijrs.emplace_back(++row, offset + 8, 1.0);
          irs.emplace_back(row, 1.0);

          // |rho(e) - e|
          do_row(offset + 0, 1.0, 1.0, 1.0, 1.0);
          do_row(offset + 3, 1.0, 1.0, 1.0, 1.0);
          do_row(offset + 6, 1.0, 1.0, 1.0, 1.0);

          // |rho(c3) - b3|
          do_row(offset + 0, c3.x, c3.y, c3.z, b3.x);
          do_row(offset + 3, c3.x, c3.y, c3.z, b3.y);
          do_row(offset + 6, c3.x, c3.y, c3.z, b3.z);

          // |lambda(e2) - b2|
          do_column(offset + 9, b2.x, b2.y, b2.z);

          // |lambda(e2) - c2|
          do_column(offset + 9, c2.x, c2.y, c2.z);
     };

     auto A = [](auto& b0, auto& b1, auto& b2, auto& b3) -> auto { return glm::inverse(hpmat3x3(b0, b1, b2)) * b3; };

     irs.reserve(2 * 15 * nEdges);
     ijrs.reserve(2 * 33 * nEdges);

     visit_edges(neighbors, [&](auto p, auto i) {
          static constexpr hpuint o[3] = { 1u, 2u, 0u };
          static constexpr hpindex o1[3] = { 2u, 0u, 1u };

          auto q = neighbors(p, i);
          auto j = make_offset(neighbors, q, p);
          auto b = std::array<Point, 3>();
          auto c = std::array<Point, 3>();
          auto& b0 = mesh.getControlPoint(q, j);
          auto& b1 = b[1];
          auto& b2 = b[0];
          auto& b3 = b[2];
          auto& c0 = c[1];
          auto& c1 = mesh.getControlPoint(p, i);
          auto& c2 = c[2];
          auto& c3 = c[0];
          auto tb = b.data() - 1;
          auto tc = c.data() - 1;

          //TODO: simplify this method using ring enumerator and write directly into b0, b1, b2, b3
          visit_subring(make_ring_enumerator(degree, neighbors, p, o[i]), neighbors(q, o1[j]), [&](auto r, auto k) { *(++tb) = patches(r)[k]; });
          visit_subring(make_ring_enumerator(degree, neighbors, q, o[j]), neighbors(p, o1[i]), [&](auto r, auto k) { *(++tc) = patches(r)[k]; });

          auto Ab2 = A(b0, b3, b1, b2);
          auto Ab3 = A(b1, b2, b0, b3);
          auto Ac2 = A(c0, c3, c1, c2);
          auto Ac3 = A(c1, c2, c0, c3);

          insert(36 * p + 12 * i, Ab2, Ab3, Ac2, Ac3);
          insert(36 * q + 12 * j, Ac3, Ac2, Ab3, Ab2);
     });

     return std::make_tuple(std::move(ijrs), std::move(irs));
}

std::vector<hpreal> make_transitions(const std::vector<hpreal>& solution);

}//namespace phm

}//namespace happah

