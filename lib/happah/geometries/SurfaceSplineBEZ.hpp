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
#include "happah/geometries/Curve.hpp"
#include "happah/geometries/TriangleGraph.hpp"
#include "happah/geometries/TriangleMesh.hpp"
#include "happah/utils/DeindexedArray.hpp"
#include "happah/utils/SurfaceSubdividerBEZ.hpp"
#include "happah/utils/VertexFactory.hpp"
#include "happah/utils/visitors.hpp"

namespace happah {

//DECLARATIONS

template<class Space, hpuint t_degree>
class SurfaceSplineBEZ;

namespace ssb {

class DeltasEnumerator;

class DiamondsEnumerator;

class NablasEnumerator;

template<class Iterator>
class PatchesEnumerator;

template<hpindex t_ring>
class RingEnumerator;

}//namespace ssb

template<hpuint degree, class Iterator>
auto de_casteljau(Iterator patch, hpreal u, hpreal v, hpreal w);

template<class Space, hpuint degree>
auto de_casteljau(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpreal u, hpreal v);

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, (degree + 1)> elevate(const SurfaceSplineBEZ<Space, degree>& surface);

template<class NewSpace, class OldSpace, hpuint degree, class Transformer>
SurfaceSplineBEZ<NewSpace, degree> embed(const SurfaceSplineBEZ<OldSpace, degree>& surface, Transformer&& transform);

template<hpuint degree, class Iterator>
auto& get_boundary_point(Iterator patch, hpindex i, hpindex k);

template<hpuint degree, class Iterator>
auto& get_boundary_point(Iterator patches, hpindex p, hpindex i, hpindex k);

template<class Space, hpuint degree>
auto& get_boundary_point(const SurfaceSplineBEZ<Space, degree>& surface, hpindex p, hpindex i, hpindex k);

template<hpuint degree, class Iterator>
auto& get_corner(Iterator patch, hpuint i);

template<hpuint degree, class Iterator>
auto& get_corner(Iterator patches, hpuint p, hpuint i);

template<hpuint degree, class Iterator>
Iterator get_patch(Iterator patches, hpuint p);

template<class Space, hpuint degree>
auto get_patch(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p);

template<class Space, hpuint degree>
bool is_c0(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i);

template<hpuint degree>
void make_bernstein_polynomials(const std::string& directory);

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_boundary(Iterator patch, hpuint i);

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_boundary(Iterator patches, hpuint p, hpuint i);

//Return the absolute offset of the jth point on the ith boundary.
hpuint make_boundary_offset(hpuint degree, hpuint i, hpuint j);

constexpr hpuint make_control_point_offset(hpuint degree, hpuint i0, hpuint i1, hpuint i2);

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
TriangleMesh<Vertex> make_control_polygon(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory());

constexpr hpuint make_control_polygon_size(hpuint degree);

template<class Iterator>
auto make_corners_enumerator(hpuint degree, Iterator begin, Iterator end);

/**
 * @param[in] nSamples Number of times an edge of the parameter triangle should be sampled.  The entire triangle is sampled uniformly such that this parameter is respected.
 * @return Matrix whose rows are the Bernstein polynomials evaluated at some point u.  The matrix is returned row-major.  To evaluate a B\'ezier polynomial at the sampled u values given a vector of control points, simply compute the product of the matrix with the vector of control points.
 */
std::vector<hpreal> make_de_casteljau_matrix(hpuint degree, hpuint nSamples);

inline ssb::DeltasEnumerator make_deltas_enumerator(hpuint degree);

template<class Transformer>
EnumeratorTransformer<ssb::DeltasEnumerator, Transformer> make_deltas_enumerator(hpuint degree, Transformer&& transform);

inline ssb::DiamondsEnumerator make_diamonds_enumerator(hpuint degree, hpuint i, hpuint j);

template<class Transformer>
EnumeratorTransformer<ssb::DiamondsEnumerator, Transformer> make_diamonds_enumerator(hpuint degree, hpuint i, hpuint j, Transformer&& transform);

//Return the absolute offset of the ith point in the interior.
hpuint make_interior_offset(hpuint degree, hpuint i);

inline ssb::NablasEnumerator make_nablas_enumerator(hpuint degree);

template<class Transformer>
EnumeratorTransformer<ssb::NablasEnumerator, Transformer> make_nablas_enumerator(hpuint degree, Transformer&& transform);

template<class Space, hpuint degree>
Indices make_neighbors(const SurfaceSplineBEZ<Space, degree>& surface);

template<class Iterator>
ssb::PatchesEnumerator<Iterator> make_patches_enumerator(hpuint degree, Iterator begin, Iterator end);

template<class Iterator, class Transformer>
EnumeratorTransformer<ssb::PatchesEnumerator<Iterator>, Transformer> make_patches_enumerator(hpuint degree, Iterator begin, Iterator end, Transformer&& transform);

constexpr hpuint make_patch_size(hpuint degree);

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_ring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i);

template<hpindex ring = 1>
ssb::RingEnumerator<ring> make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i);

template<hpindex ring, class Transformer>
EnumeratorTransformer<ssb::RingEnumerator<ring>, Transformer> make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i, Transformer&& transform);

template<hpindex ring, class Space, hpuint degree>
auto make_ring_enumerator(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i);

//Convert a string representation in HPH format of a spline surface into a spline surface.
template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> make_spline_surface(const std::string& surface);

//Import a spline surface stored in the given file in HPH format.
template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> make_spline_surface(const std::experimental::filesystem::path& surface);
     
//Convert a closed(!) triangle graph into a quartic polynomial spline.
template<class Vertex>
auto make_spline_surface(const TriangleGraph<Vertex>& graph);

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions, VertexFactory&& factory = VertexFactory());

std::tuple<std::vector<hpcolor>, std::vector<hpcolor> > paint_boundary_edges(hpuint degree, std::vector<hpcolor> Vcolors, std::vector<hpcolor> Ecolors, const hpcolor& color);

std::vector<hpcolor> paint_boundary_triangles(hpuint degree, std::vector<hpcolor> colors, const hpcolor& color0, const hpcolor& color1);

template<hpuint degree, class Iterator, class Visitor>
void sample(Iterator patches, hpuint nPatches, hpuint nSamples, Visitor&& visit);

template<hpuint degree, class ControlPointsIterator, class DomainPointsIterator, class Visitor>
void sample(ControlPointsIterator controlPoints, DomainPointsIterator domainPoints, hpuint nPatches, hpuint nSamples, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSamples, Visitor&& visit);

template<class Space, hpuint degree, class T, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, std::tuple<const std::vector<T>&, const Indices&> domain, hpuint nSamples, Visitor&& visit);

template<class Space, hpuint degree>
auto size(const SurfaceSplineBEZ<Space, degree>& surface);

//Return a G1 surface that interpolates the positions and the tangents planes at the corners of the patches in the given surface.
template<hpuint degree>
SurfaceSplineBEZ<Space4D, degree> smooth(const SurfaceSplineBEZ<Space4D, degree>& surface, const std::vector<hpreal>& transitions, hpreal epsilon = EPSILON);

template<hpuint degree>
SurfaceSplineBEZ<Space4D, degree> smooth(const SurfaceSplineBEZ<Space3D, degree>& surface, const std::vector<hpreal>& transitions, hpreal epsilon = EPSILON);

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions);

template<class Visitor>
void visit_bernstein_indices(hpuint degree, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_boundary(Iterator patch, hpuint i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_boundary(Iterator patches, hpuint p, hpuint i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_corner(Iterator patch, hpuint i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_corner(Iterator patches, hpuint p, hpuint i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_corners(Iterator patch, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_corners(Iterator patches, hpuint p, Visitor&& visit);

//Visit triangles in control polygon schematically pointing up.
//template<class Iterator, class Visitor>
//void visit_deltas(hpuint degree, Iterator patch, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_edges(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_ends(Iterator patch, hpuint i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_ends(Iterator patches, hpuint p, hpuint i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_interior(Iterator patch, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_interior(Iterator patches, hpuint p, Visitor&& visit);

//Visit triangles in control polygon schematically pointing down.  The points are given in counterclockwise order; the first point is the top right point.
template<class Iterator, class Visitor>
void visit_nablas(hpuint degree, Iterator patch, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_patch(Iterator patches, hpuint p, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_patch(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_patches(Iterator patches, hpuint nPatches, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_patches(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit);

template<hpindex ring, class Visitor>
void visit_ring(ssb::RingEnumerator<ring> e, Visitor&& visit);

template<hpindex ring, class Transformer, class Visitor>
void visit_ring(EnumeratorTransformer<ssb::RingEnumerator<ring>, Transformer> e, Visitor&& visit);

//Visit the subring stopping at patch p.
template<class Visitor>
void visit_subring(ssb::RingEnumerator<1> e, hpindex p, Visitor&& visit);

//DEFINITIONS

template<class Space, hpuint t_degree>
class SurfaceSplineBEZ {
     using Point = typename Space::POINT;
     using ControlPoints = std::vector<Point>;

public:
     SurfaceSplineBEZ() {}

     SurfaceSplineBEZ(hpuint n)
          : m_controlPoints(1, Point(0)), m_indices(n * make_patch_size(t_degree), 0) {}

     SurfaceSplineBEZ(ControlPoints controlPoints)
          : m_controlPoints(std::move(controlPoints)), m_indices(m_controlPoints.size()) { std::iota(std::begin(m_indices), std::end(m_indices), 0); }

     SurfaceSplineBEZ(ControlPoints controlPoints, Indices indices)
          : m_controlPoints{std::move(controlPoints)}, m_indices{std::move(indices)} {}

     auto& getControlPoints() const { return m_controlPoints; }

     std::tuple<const ControlPoints&, const Indices&> getPatches() const { return std::tie(m_controlPoints, m_indices); }

     //Set the ith boundary of the pth patch.
     template<class Iterator>
     void setBoundary(hpindex p, hpindex i, Iterator begin) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          auto n = m_controlPoints.size();
          visit_boundary<t_degree>(std::begin(m_indices), p, i, [&](auto& i) { i = n++; });
          m_controlPoints.insert(std::end(m_controlPoints), begin, begin + (t_degree - 1));
     }

     //Set the ith boundary of the pth patch to the jth boundary of the qth patch.
     void setBoundary(hpindex p, hpindex i, hpindex q, hpindex j) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          auto boundary = make_boundary<t_degree>(std::begin(m_indices), q, j);
          auto n = std::end(boundary);
          visit_boundary<t_degree>(std::begin(m_indices), p, i, [&](auto& i) { i = *(--n); });
     }

     //Set the kth point on the ith boundary of the pth patch.
     void setBoundaryPoint(hpindex p, hpindex i, hpindex k, Point point) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          i = make_boundary_offset(t_degree, i, k);
          get_patch<t_degree>(std::begin(m_indices), p)[i] = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     //Set the kth point on the ith boundary of the pth patch and the point opposite to it on the jth boundary of the qth patch.
     void setBoundaryPoint(hpindex p, hpindex i, hpindex k, hpindex q, hpindex j, Point point) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          i = make_boundary_offset(t_degree, i, k);
          j = make_boundary_offset(t_degree, j, t_degree - 2 - k);
          get_patch<t_degree>(std::begin(m_indices), p)[i] = m_controlPoints.size();
          get_patch<t_degree>(std::begin(m_indices), q)[j] = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     //Set the kth point on the ith boundary of the pth patch to the opposite point on the jth boundary of the qth patch.
     void setBoundaryPoint(hpindex p, hpindex i, hpindex k, hpindex q, hpindex j) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          i = make_boundary_offset(t_degree, i, k);
          j = make_boundary_offset(t_degree, j, t_degree - 2 - k);
          get_patch<t_degree>(std::begin(m_indices), p)[i] = get_patch<t_degree>(std::begin(m_indices), q)[j];
     }

     void setControlPoint(hpindex p, hpindex i, Point point) {
          get_patch<t_degree>(std::begin(m_indices), p)[i] = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     void setCorner(hpindex p, hpindex i, Point point) {
          static_assert(t_degree > 0, "There is no corner in a constant.");
          get_corner<t_degree>(std::begin(m_indices), p, i) = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     void setCorner(hpindex p, hpindex i, hpindex q, hpindex j) {
          static_assert(t_degree > 0, "There is no corner in a constant.");
          get_corner<t_degree>(std::begin(m_indices), p, i) = get_corner<t_degree>(std::begin(m_indices), q, j);
     }

     template<class Iterator>
     void setInterior(hpindex p, Iterator begin) {
          static_assert(t_degree > 2, "There is no interior in a constant, linear, or quadratic.");
          auto n = m_controlPoints.size();
          visit_interior<t_degree>(std::begin(m_indices), p, [&](auto& i) { i = n++; });
          m_controlPoints.insert(std::end(m_controlPoints), begin, begin + (make_patch_size(t_degree) - 3 * t_degree)); 
     }

     void setInteriorPoint(hpindex p, hpindex i, Point point) {
          static_assert(t_degree > 2, "There is no interior in a constant, linear, or quadratic.");
          static constexpr auto patchSize = make_patch_size(t_degree);
          i = make_interior_offset(t_degree, i);
          get_patch<t_degree>(std::begin(m_indices), p)[i] = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

private:
     ControlPoints m_controlPoints;
     Indices m_indices;

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const SurfaceSplineBEZ<Space, t_degree>& surface) {
          using happah::format::hph::operator<<;

          stream << surface.m_controlPoints << '\n';
          stream << surface.m_indices;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, SurfaceSplineBEZ<Space, t_degree>& surface) {
          using happah::format::hph::operator>>;

          stream >> surface.m_controlPoints;
          stream >> surface.m_indices;
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

namespace ssb {

class DeltasEnumerator {
public:
     DeltasEnumerator(hpuint degree)
          : m_bottom(0u), m_delta(degree), m_end(degree) {}
     
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
     hpuint m_bottom;
     hpuint m_delta;
     hpuint m_end;

};//DeltasEnumerator
 
/*
 *   k3
 * k0  k2
 *   k1
 * i is with respect to top patch; j is with respect to bottom patch.
 * k0 and k2 are control point indices with respect to the ith patch.
 */
class DiamondsEnumerator {
public:
     DiamondsEnumerator(hpuint degree, hpuint i, hpuint j)
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
     hpuint m_i;
     hpuint m_j;
     hpuint m_k0;
     hpuint m_k1;
     hpuint m_k2;
     hpuint m_k3;

};//DiamondsEnumerator

class NablasEnumerator {
public:
     NablasEnumerator(hpuint degree)
          : m_bottom(1u), m_delta(degree), m_end(degree) {}

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
     hpuint m_bottom;
     hpuint m_delta;
     hpuint m_end;

};//NablasEnumerator

template<class Iterator>
class PatchesEnumerator {
public:
     PatchesEnumerator(hpuint degree, Iterator begin, Iterator end)
          : m_begin(begin), m_end(end), m_patchSize(make_patch_size(degree)) {}

     explicit operator bool() const { return m_begin != m_end; }

     auto operator*() const { return m_begin; }

     auto& operator++() {
          m_begin += m_patchSize;
          return *this;
     }

private:
     Iterator m_begin;
     Iterator m_end;
     hpuint m_patchSize;

};//PatchesEnumerator

template<>
class RingEnumerator<1> {
public:
     RingEnumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i)
          : m_e(make_spokes_enumerator(neighbors, p, i)), m_o{ 1u, degree << 1, make_patch_size(degree) - 3u } {}

     explicit operator bool() const { return bool(m_e); }

     auto operator*() const {
          auto p = 0u, i = 0u;
          std::tie(p, i) = *m_e;
          return std::make_tuple(p, m_o[i]);
     }

     auto& operator++() {
          ++m_e;
          return *this;
     }

private:
     trm::SpokesEnumerator m_e;
     hpindex m_o[3];

};//RingEnumerator

template<>
class RingEnumerator<2> {
public:
     RingEnumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i)
          : m_e(make_spokes_enumerator(neighbors, p, i)), m_o0{ 2, 3 * degree - 1, make_patch_size(degree) - 6 }, m_o1{ degree + 2, (degree << 1) - 1, make_patch_size(degree) - 5 } { m_o = m_o0; }

     explicit operator bool() const { return bool(m_e); }

     auto operator*() const {
          auto p = 0u, i = 0u;
          std::tie(p, i) = *m_e;
          return std::make_tuple(p, m_o[i]);
     }

     auto& operator++() {
          if(m_o == m_o0) m_o = m_o1;
          else {
               m_o = m_o0;
               ++m_e;
          }
          return *this;
     }

private:
     trm::SpokesEnumerator m_e;
     hpindex* m_o;
     hpindex m_o0[3];
     hpindex m_o1[3];

};//RingEnumerator

}//namespace ssb

template<hpuint degree, class Iterator>
auto de_casteljau(Iterator patch, hpreal u, hpreal v, hpreal w) {
     using T = typename std::iterator_traits<Iterator>::value_type;

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
auto de_casteljau(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpreal u, hpreal v) { return de_casteljau<degree>(get_patch(surface, p), u, v, 1.0f - u - v); }

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, (degree + 1)> elevate(const SurfaceSplineBEZ<Space, degree>& surface) {
     SurfaceSplineBEZ<Space, (degree + 1)> surface1(size(surface));
     auto& indices = std::get<1>(surface.getPatches());
     auto neighbors = make_neighbors(surface);
     auto patches = deindex(surface.getPatches());

     visit_vertices(neighbors, [&](auto p, auto i) {
          surface1.setCorner(p, i, get_corner<degree>(std::begin(patches), p, i));
          visit_spokes(make_spokes_enumerator(neighbors, p, i), [&](auto q, auto j) { surface1.setCorner(q, j, p, i); });
     });

     auto elevate_boundary = [&](auto p, auto i) {
          auto patch = get_patch<degree>(std::begin(patches), p);
          auto boundary = make_boundary<degree>(patch, i);
          visit_ends<degree>(patch, i, [&](auto& corner0, auto& corner1) { surface1.setBoundary(p, i, std::begin(happah::curves::elevate(degree, corner0, std::begin(boundary), corner1))); });
     };
     visit_edges(neighbors, [&](auto p, auto i) {
          elevate_boundary(p, i);
          auto q = neighbors[3 * p + i];
          if(q == UNULL) return;
          auto j = make_neighbor_offset(neighbors, q, p);
          if(is_c0(surface, neighbors, p, i)) surface1.setBoundary(q, j, p, i);
          else elevate_boundary(q, j);
     });

     if(degree < 2) return surface1;

     auto p = 0u;
     visit_patches<degree>(std::begin(patches), size(surface), [&](auto patch) {
          using Point = typename Space::POINT;
          std::vector<Point> interior;
          interior.reserve(make_patch_size(degree + 1u) - 3u * (degree + 1u));
          auto alpha = hpreal(1.0 / (degree + 1u));
          auto t0 = degree;
          auto j0 = 0u, j1 = 0u, j2 = 0u;
          visit_nablas(degree, patch, [&](auto& b0, auto& b1, auto& b2) {
               if(j0 == 0u) {
                    j0 = --t0;
                    j2 = degree - t0 + 1u;
                    j1 = degree - j0 - j2 + 1u;
               }
               interior.push_back(alpha * (hpreal(j0) * b0 + hpreal(j1) * b1 + hpreal(j2) * b2));
               --j0;
               ++j1;
          });
          surface1.setInterior(p, std::begin(interior));
          ++p;
     });
     return surface1;
}

template<class NewSpace, class OldSpace, hpuint degree, class Transformer>
SurfaceSplineBEZ<NewSpace, degree> embed(const SurfaceSplineBEZ<OldSpace, degree>& surface, Transformer&& transform) {
     auto& oldPoints = std::get<0>(surface.getPatches());
     auto& indices = std::get<1>(surface.getPatches());
     auto newPoints = std::vector<typename NewSpace::POINT>(oldPoints.size());
     std::transform(std::begin(oldPoints), std::end(oldPoints), std::begin(newPoints), transform);
     return { std::move(newPoints), indices };
}

template<hpuint degree, class Iterator>
auto& get_boundary_point(Iterator patch, hpindex i, hpindex k) {
     static_assert(degree > 1, "There is no boundary in a constant or linear.");
     return patch[make_boundary_offset(degree, i, k)];
}

template<hpuint degree, class Iterator>
auto& get_boundary_point(Iterator patches, hpindex p, hpindex i, hpindex k) { return get_boundary_point<degree>(get_patch<degree>(patches, p), i, k); }

template<class Space, hpuint degree>
auto& get_boundary_point(const SurfaceSplineBEZ<Space, degree>& surface, hpindex p, hpindex i, hpindex k) {
     auto patches = deindex(surface.getPatches());
     return get_boundary_point<degree>(std::begin(patches), p, i, k);
}

template<hpuint degree, class Iterator>
auto& get_corner(Iterator patch, hpuint i) {
     static constexpr hpuint o[3] = { 0u, degree, make_patch_size(degree) - 1u };
     return patch[o[i]];
}

template<hpuint degree, class Iterator>
auto& get_corner(Iterator patches, hpuint p, hpuint i) { return get_corner<degree>(get_patch<degree>(patches, p), i); }

template<hpuint degree, class Iterator>
Iterator get_patch(Iterator patches, hpuint p) {
     static constexpr auto patchSize = make_patch_size(degree);
     return patches + p * patchSize;
}

template<class Space, hpuint degree>
auto get_patch(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p) {
     auto patches = deindex(surface.getPatches());
     return get_patch<degree>(std::begin(patches), p);
}

template<class Space, hpuint degree>
bool is_c0(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i) {
     auto& indices = std::get<1>(surface.getPatches());
     auto q = neighbors[3 * p + i];
     auto j = make_neighbor_offset(neighbors, q, p);

     auto k0 = 0u, k1 = 0u, l0 = 0u, l1 = 0u;
     visit_ends<degree>(std::begin(indices), p, i, [&](auto k, auto l) { k0 = k; l0 = l; });
     visit_ends<degree>(std::begin(indices), q, j, [&](auto k, auto l) { k1 = k; l1 = l; });
     if(k0 != l1 || l0 != k1) return false;
     auto boundary0 = make_boundary<degree>(std::begin(indices), p, i); 
     auto boundary1 = make_boundary<degree>(std::begin(indices), q, j);
     std::reverse(std::begin(boundary1), std::end(boundary1));
     return boundary0 == boundary1;
}

template<hpuint degree>
void make_bernstein_polynomials(const std::string& directory) {
     auto points = std::vector<Point3D>();
     points.reserve(make_patch_size(degree));
     auto plane = LinearSurfaceSplineBEZ<Space2D>({ Point2D(0.0, 0.0), Point2D(1.0, 0.0), Point2D(0.0, 1.0) }, { 0, 1, 2 });
     sample(plane, degree + 1u, [&](auto sample) { points.emplace_back(sample.x, sample.y, 0.0); });
     visit_bernstein_indices(degree, [&](auto i, auto j, auto k) {
          auto temp = points;
          temp[make_control_point_offset(degree, i, j, k)].z = 1.0;
          auto surface = SurfaceSplineBEZ<Space3D, degree>(temp);
          std::ostringstream path;
          path << directory << "/b" << i << j << k << ".ss" << degree << ".bz.3.hph";
          write(surface, path.str());
     });
}

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_boundary(Iterator patch, hpuint i) {
     auto boundary = std::vector<typename std::iterator_traits<Iterator>::value_type>();
     boundary.reserve(degree - 1);
     visit_boundary<degree>(patch, i, make_back_inserter(boundary));
     return boundary;
}

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_boundary(Iterator patches, hpuint p, hpuint i) { return make_boundary<degree>(get_patch<degree>(patches, p), i); }

constexpr hpuint make_control_point_offset(hpuint degree, hpuint i0, hpuint i1, hpuint i2) { return make_patch_size(degree) - make_patch_size(degree - i2) + i1; }

template<class Space, hpuint degree, class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_control_polygon(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory) { return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, 0, std::forward<VertexFactory>(factory)); }

constexpr hpuint make_control_polygon_size(hpuint degree) { return degree * degree; }

template<class Iterator>
auto make_corners_enumerator(hpuint degree, Iterator begin, Iterator end) { return make_patches_enumerator(degree, begin, end, [&](auto patch) { return std::tie(patch[0], patch[degree], patch[make_patch_size(degree) - 1]); }); }

inline ssb::DeltasEnumerator make_deltas_enumerator(hpuint degree) { return { degree }; };

template<class Transformer>
EnumeratorTransformer<ssb::DeltasEnumerator, Transformer> make_deltas_enumerator(hpuint degree, Transformer&& transform) { return { make_deltas_enumerator(degree), std::forward<Transformer>(transform) }; }

inline ssb::DiamondsEnumerator make_diamonds_enumerator(hpuint degree, hpuint i, hpuint j) { return { degree, i, j }; };

template<class Transformer>
EnumeratorTransformer<ssb::DiamondsEnumerator, Transformer> make_diamonds_enumerator(hpuint degree, hpuint i, hpuint j, Transformer&& transform) { return { make_diamonds_enumerator(degree, i, j), std::forward<Transformer>(transform) }; }

template<class Space, hpuint degree>
Indices make_neighbors(const SurfaceSplineBEZ<Space, degree>& surface) {
     auto& indices = std::get<1>(surface.getPatches());
     auto corners = expand(make_corners_enumerator(degree, std::begin(indices), std::end(indices)));
     return make_neighbors(corners);
}

inline ssb::NablasEnumerator make_nablas_enumerator(hpuint degree) { return { degree }; };

template<class Transformer>
EnumeratorTransformer<ssb::NablasEnumerator, Transformer> make_nablas_enumerator(hpuint degree, Transformer&& transform) { return { make_nablas_enumerator(degree), std::forward<Transformer>(transform) }; }

template<class Iterator>
ssb::PatchesEnumerator<Iterator> make_patches_enumerator(hpuint degree, Iterator begin, Iterator end) { return { degree, begin, end }; }

template<class Iterator, class Transformer>
EnumeratorTransformer<ssb::PatchesEnumerator<Iterator>, Transformer> make_patches_enumerator(hpuint degree, Iterator begin, Iterator end, Transformer&& transform) { return { make_patches_enumerator(degree, begin, end), std::forward<Transformer>(transform) }; }

constexpr hpuint make_patch_size(hpuint degree) { return (degree + 1) * (degree + 2) >> 1; }

template<hpuint degree, class Iterator, class T = typename std::iterator_traits<Iterator>::value_type>
std::vector<T> make_ring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i) {
     auto ring = std::vector<T>();
     visit_ring<degree>(patches, neighbors, p, i, make_back_inserter(ring));
     return ring;
}

template<hpindex ring>
ssb::RingEnumerator<ring> make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i) { return { degree, neighbors, p, i }; }

template<hpindex ring, class Transformer>
EnumeratorTransformer<ssb::RingEnumerator<ring>, Transformer> make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i, Transformer&& transform) { return { make_ring_enumerator<ring>(degree, neighbors, p, i), std::forward<Transformer>(transform) }; }

template<hpindex ring, class Space, hpuint degree>
auto make_ring_enumerator(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpindex p, hpindex i) {
     auto patches = deindex(surface.getPatches());
     return make_ring_enumerator<ring>(degree, neighbors, p, i, [&](auto p, auto i) { return get_patch<degree>(std::begin(patches), p)[i]; });
}

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> make_spline_surface(const std::string& surface) { return format::hph::read<SurfaceSplineBEZ<Space, degree> >(surface); }

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> make_spline_surface(const std::experimental::filesystem::path& surface) { return format::hph::read<SurfaceSplineBEZ<Space, degree> >(surface); }

template<class Vertex>
auto make_spline_surface(const TriangleGraph<Vertex>& graph) {
     using Space = typename Vertex::SPACE;
     using Vector = typename Space::VECTOR;

     auto surface = QuarticSurfaceSplineBEZ<Space>(graph.getNumberOfTriangles());

     auto set_boundary_point = [&](auto t, auto i, auto k, auto&& point) {
          auto u = make_neighbor_index(graph, t, i);
          auto j = make_neighbor_offset(graph, u, t);
          surface.setBoundaryPoint(t, i, k, u, j, point);
     };

     visit_diamonds(graph, [&](auto e, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
          auto t = make_triangle_index(e);
          auto i = make_edge_offset(e);
          set_boundary_point(t, i, 1, (hpreal(1.0) / hpreal(6.0)) * (hpreal(2.0) * vertex0.position + vertex1.position + hpreal(2.0) * vertex2.position + vertex3.position));
     });

     for(auto v : boost::irange(0u, graph.getNumberOfVertices())) {
          auto ring = make_ring(make_ring_enumerator(graph, v));
          auto begin = std::begin(ring);
          auto end = std::end(ring);
          auto t = make_triangle_index(graph.getOutgoing(v));
          auto i = make_edge_offset(graph.getOutgoing(v));
          auto& center = graph.getVertex(v);
          auto fan = Indices();
          visit_spokes(make_spokes_enumerator(graph.getEdges(), graph.getOutgoing(v)), [&](auto e) {
               auto u = make_triangle_index(e);
               auto j = make_edge_offset(e);
               fan.push_back(u);
               fan.push_back(j);
          });
          auto valence = std::distance(begin, end);

          auto corner = center.position;
          if(valence == 6) {
               corner *= hpreal(6.0);
               for(auto& vertex : boost::make_iterator_range(begin, end)) corner += vertex.position;
               corner /= hpreal(12.0);

               auto make_boundary_point = [&](auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3, auto& vertex4) -> auto { return (hpreal(1.0) / hpreal(24.0)) * (hpreal(12.0) * center.position + vertex0.position + hpreal(3.0) * vertex1.position + hpreal(4.0) * vertex2.position + hpreal(3.0) * vertex3.position + vertex4.position); };

               set_boundary_point(fan[0], fan[1], 0, make_boundary_point(begin[4], begin[5], begin[0], begin[1], begin[2]));
               set_boundary_point(fan[2], fan[3], 0, make_boundary_point(begin[5], begin[0], begin[1], begin[2], begin[3]));
               set_boundary_point(fan[4], fan[5], 0, make_boundary_point(begin[0], begin[1], begin[2], begin[3], begin[4]));
               set_boundary_point(fan[6], fan[7], 0, make_boundary_point(begin[1], begin[2], begin[3], begin[4], begin[5]));
               set_boundary_point(fan[8], fan[9], 0, make_boundary_point(begin[2], begin[3], begin[4], begin[5], begin[0]));
               set_boundary_point(fan[10], fan[11], 0, make_boundary_point(begin[3], begin[4], begin[5], begin[0], begin[1]));
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
                    auto vector = get_boundary_point(surface, t, i, 1) - corner;
                    auto r = std::fmin(glm::length2(vector) / std::abs(glm::dot(tangent, vector)), glm::length(vector)) / hpreal(2.0);
                    set_boundary_point(f[0], f[1], 0, corner + r * tangent);
               }
          }
          surface.setCorner(t, i, corner);
          visit_pairs(fan, [&](auto u, auto j) { surface.setCorner(u, j, t, i); });

          auto set_interior_point = [&](auto t, auto i, auto& vertex0, auto& vertex1, auto& vertex2, auto& vertex3) {
               auto point = (hpreal(1.0) / hpreal(24.0)) * (hpreal(10.0) * center.position + vertex0.position + hpreal(6.0) * vertex1.position + hpreal(6.0) * vertex2.position + vertex3.position);
               surface.setInteriorPoint(t, i, point);
          };

          set_interior_point(fan[0], fan[1], begin[valence - 1], begin[0], begin[1], begin[2]);
          if(valence > 3) {
               auto middle = begin;
               visit_pairs(std::begin(fan) + 2, valence - 3, 2, [&](auto u, auto j) {
                    set_interior_point(u, j, middle[0], middle[1], middle[2], middle[3]);
                    ++middle;
               });
          }
          auto n = (valence << 1) - 4;
          set_interior_point(fan[n], fan[n + 1], begin[valence - 3], begin[valence - 2], begin[valence - 1], begin[0]);
          set_interior_point(fan[n + 2], fan[n + 3], begin[valence - 2], begin[valence - 1], begin[0], begin[1]);
     }

     return surface;
}

template<class Space, hpuint degree, class Vertex, class VertexFactory, typename>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions, VertexFactory&& factory) {
     using Point = typename Space::POINT;
     static_assert(degree > 0u, "Constant spline surfaces cannot be converted into triangle meshes.");
     static_assert(std::is_base_of<Vertex, decltype(factory(Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");

     auto do_make_triangle_mesh = [&](const SurfaceSplineBEZ<Space, degree>& surface) -> auto {
          auto vertices = std::vector<Vertex>();
          vertices.reserve(surface.getControlPoints().size());
          for(auto& point : surface.getControlPoints()) vertices.push_back(factory(point));

          auto indices = Indices();
          indices.reserve(3 * make_control_polygon_size(degree) * size(surface));
          auto inserter = make_back_inserter(indices);
          visit_patches<degree>(std::begin(std::get<1>(surface.getPatches())), size(surface), [&](auto patch) {
               visit(make_deltas_enumerator(degree, [&](auto i0, auto i1, auto i2) { return std::tie(patch[i0], patch[i1], patch[i2]); }), inserter);
               visit_nablas(degree, patch, inserter);
          });

          return make_triangle_mesh(std::move(vertices), std::move(indices));
     };

     if(nSubdivisions > 0) return do_make_triangle_mesh(subdivide(surface, nSubdivisions));
     else return do_make_triangle_mesh(surface);
}

//TODO: move non-member functions with iterators into subnamespace so as not to conflict with implementations for curves, for example
template<hpuint degree, class Iterator, class Visitor>
void sample(Iterator patches, hpuint nPatches, hpuint nSamples, Visitor&& visit) {
     //TODO; skip multiple computations on common edges and eliminate common points in array
     auto matrix = make_de_casteljau_matrix(degree, nSamples);
     visit_patches<degree>(patches, nPatches, [&](auto patch) {
          static constexpr auto patchSize = make_patch_size(degree);
          auto end = patch + patchSize;
          for(auto m = std::begin(matrix), mend = std::end(matrix); m != mend; ++m) {
               auto temp = patch;
               auto sample = *m * *temp;
               while(++temp != end) sample += *(++m) * *temp;
               visit(sample);
          }
     });
}

template<hpuint degree, class ControlPointsIterator, class DomainPointsIterator, class Visitor>
void sample(ControlPointsIterator controlPoints, DomainPointsIterator domainPoints, hpuint nPatches, hpuint nSamples, Visitor&& visit) {
     //TODO; skip multiple computations on common edges and eliminate common points in array
     auto matrixd = make_de_casteljau_matrix(degree, nSamples);
     auto matrix1 = make_de_casteljau_matrix(1, nSamples);
     visit_patches<degree>(controlPoints, nPatches, [&](auto patch) {
          static constexpr auto patchSize = make_patch_size(degree);
          auto end = patch + patchSize;
          auto& p0 = *domainPoints;
          auto& p1 = *(++domainPoints);
          auto& p2 = *(++domainPoints);
          ++domainPoints;

          auto d = std::begin(matrix1);
          for(auto m = std::begin(matrixd), mend = std::end(matrixd); m != mend; ++m) {
               auto u = *d;
               auto v = *(++d);
               auto w = *(++d);
               ++d;

               auto temp = patch;
               auto sample = *m * *temp;
               while(++temp != end) sample += *(++m) * *temp;
               visit(mix(p0, u, p1, v, p2, w), sample);
          }
     });
}

template<class Space, hpuint degree, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSamples, Visitor&& visit) {
     auto patches = deindex(surface.getPatches());
     sample<degree>(std::begin(patches), size(surface), nSamples, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class T, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, std::tuple<const std::vector<T>&, const Indices&> domain, hpuint nSamples, Visitor&& visit) {
     auto controlPoints = deindex(surface.getPatches());
     auto domainPoints = deindex(domain);
     sample<degree>(std::begin(controlPoints), std::begin(domainPoints), size(surface), nSamples, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree>
auto size(const SurfaceSplineBEZ<Space, degree>& surface) { return std::get<1>(surface.getPatches()).size() / make_patch_size(degree); }

template<hpuint degree>
SurfaceSplineBEZ<Space4D, degree> smooth(const SurfaceSplineBEZ<Space4D, degree>& surface, const std::vector<hpreal>& transitions, hpreal epsilon) {
     static_assert(degree > 4u, "The first two rings of control points surrounding the corners of the patches are assumed to be disjoint.");
     
     auto neighbors = make_neighbors(surface);
     auto patches = deindex(surface.getPatches());
     auto surface1 = SurfaceSplineBEZ<Space4D, degree>(size(surface));

     auto get_transition = [&](auto p, auto i) {
          auto q = make_neighbor_index(neighbors, p, i);
          auto j = make_neighbor_offset(neighbors, q, p);
          return get_triplet(transitions, 3 * q + j);
     };

     auto make_coefficients_1 = [&](auto p, auto i, auto valence, auto& center) {
          auto coefficients = std::vector<double>(valence * 7, 0.0);
          auto a0 = std::begin(coefficients) - 1;
          auto a1 = a0 + valence;
          auto a2 = a1 + valence;
          auto b0 = a2 + valence;
          auto b1 = b0 + valence;
          auto b2 = b1 + valence;
          auto b3 = b2 + valence;
          auto e = make_ring_enumerator<1>(degree, neighbors, p, i, [&](auto q, auto j) { return get_patch<degree>(std::begin(patches), q)[j]; });
          auto f = make_spokes_enumerator(neighbors, p, i);

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
               auto l0 = hpreal(0.0), l1 = hpreal(0.0), l2 = hpreal(0.0);
               auto q = 0u, j = 0u;
               std::tie(q, j) = *f;
               std::tie(l0, l1, l2) = get_transition(q, j);
               auto t0 = l0 + l1 * a0[0] + l2 * (a0 - 1)[0];
               auto t1 = l1 * a1[0] + l2 * (a1 - 1)[0];
               auto t2 = l1 * a2[0] + l2 * (a2 - 1)[0];
               push_back(t0, t1, t2);
               ++f;
          }

          return coefficients;
     };

     auto set_ring_1 = [&](auto p, auto i, auto valence, auto& center, auto& coefficients) {
          auto a0 = std::begin(coefficients) - 1;
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
               auto q = 0u, j = 0u;
               std::tie(q, j) = *e;
               auto r = make_neighbor_index(neighbors, q, j);
               auto k = make_neighbor_offset(neighbors, r, q);
               surface1.setBoundaryPoint(q, j, 0, r, k, point);
               ++e;
          };

          set_boundary_point(q1);
          set_boundary_point(q2);

          while(e) set_boundary_point(hpreal((++a0)[0]) * center + hpreal((++a1)[0]) * q1 + hpreal((++a2)[0]) * q2);
     };

     auto make_coefficients_2 = [&](auto p, auto i, auto valence) {
          auto nRows = valence << 1;
          auto coefficients = std::vector<double>(nRows * (valence + 4), 0.0);
          auto b0 = std::begin(coefficients) - 1;
          auto b1 = b0 + nRows;
          auto b2 = b1 + nRows;
          auto b3 = b2 + nRows;
          auto A = b3 + (nRows + 1);
          auto r0 = hpreal(0.0), r1 = hpreal(0.0), r2 = hpreal(0.0), r3 = hpreal(0.0);
          auto e1 = make_ring_enumerator<1>(degree, neighbors, p, i, [&](auto q, auto j) { return get_patch<degree>(std::begin(patches), q)[j]; });
          auto e2 = make_ring_enumerator<2>(degree, neighbors, p, i, [&](auto q, auto j) { return get_patch<degree>(std::begin(patches), q)[j]; });
          auto f = make_spokes_enumerator(neighbors, p, i);
          auto n = 0;

          auto push_back_0 = [&](auto point) {
               (++b0)[0] = point.x;
               (++b1)[0] = point.y;
               (++b2)[0] = point.z;
               (++b3)[0] = point.w;
               A[n * nRows] = 1.0;
               ++A;
          };

          auto push_back_1 = [&](auto q1, auto p3, auto l0, auto l1, auto l2) {
               r0 = l0 * q1.x + l2 * r0;
               r1 = l0 * q1.y + l2 * r1;
               r2 = l0 * q1.z + l2 * r2;
               r3 = l0 * q1.w + l2 * r3;
               (++b0)[0] = p3.x - r0;
               (++b1)[0] = p3.y - r1;
               (++b2)[0] = p3.z - r2;
               (++b3)[0] = p3.w - r3;
               auto temp = A;
               repeat(n, [&]() {
                    temp[0] = l2 * (temp - 2)[0];
                    temp += nRows;
               });
               temp[0] = l1;
               ++A;
          };

          //std::tie(p, i) = *f;//TODO: is this necessary?

          auto q3 = *e1;
          auto p6 = *e2;
          auto p1 = *(++e2);
          push_back_0(p1);
          ++e1;
          ++e2;
          ++f;
          ++n;

          while(e1) {
               auto q1 = *e1;
               auto p2 = *e2;
               auto p3 = *(++e2);
               auto l0 = hpreal(0.0), l1 = hpreal(0.0), l2 = hpreal(0.0);
               auto q = 0u, j = 0u;
               std::tie(q, j) = *f;
               std::tie(l0, l1, l2) = get_transition(q, j);
               push_back_0(p2);
               push_back_1(q1, p3, l0, l1, l2);
               ++e1;
               ++e2;
               ++f;
               ++n;
          }

          auto l0 = hpreal(0.0), l1 = hpreal(0.0), l2 = hpreal(0.0);
          std::tie(l0, l1, l2) = get_transition(p, i);
          if(std::abs(l1) < epsilon) {
               (++b0)[0] = l0 * q3.x + l2 * r0;
               (++b1)[0] = l0 * q3.y + l2 * r1;
               (++b2)[0] = l0 * q3.z + l2 * r2;
               (++b3)[0] = l0 * q3.w + l2 * r3;
               auto temp = A;
               repeat(valence, [&]() {
                    temp[0] = -l2 * (temp - 1)[0];
                    temp += nRows;
               });
               A[0] += 1.0;
          } else {
               (++b0)[0] = p6.x + (l0 * q3.x + l2 * r0) / l1;
               (++b1)[0] = p6.y + (l0 * q3.y + l2 * r1) / l1;
               (++b2)[0] = p6.z + (l0 * q3.z + l2 * r2) / l1;
               (++b3)[0] = p6.w + (l0 * q3.w + l2 * r3) / l1;
               auto temp = A;
               repeat(valence, [&]() {
                    temp[0] = (-l2 / l1) * (temp - 1)[0];
                    temp += nRows;
               });
               A[0] += 1.0 / l1;
          }

          return coefficients;
     };

     auto set_ring_2 = [&](auto p, auto i, auto valence, auto& coefficients) {
          auto nRows = valence << 1;
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
          auto e1 = make_ring_enumerator<1>(degree, neighbors, p, i, [&](auto q, auto j) { return get_patch<degree>(std::begin(patches), q)[j]; });
          auto e2 = make_ring_enumerator<2>(degree, neighbors, p, i);
          auto f = make_spokes_enumerator(neighbors, p, i);
          auto n = 0;

          auto set_boundary_point = [&](auto q, auto j, auto point) {
               auto r = make_neighbor_index(neighbors, q, j);
               auto k = make_neighbor_offset(neighbors, r, q);
               surface1.setBoundaryPoint(q, j, 1, r, k, point);
          };

          auto set_interior_point = [&](auto point) {
               static constexpr hpuint patchSize = make_patch_size(degree);

               auto q = 0u, j = 0u;
               std::tie(q, j) = *(++e2);
               surface1.setControlPoint(q, j, point);
          };

          //std::tie(p, i) = *f;//TODO: is this necessary?

          auto qb = *e1;
          auto p1 = Point4D(x0[0], x1[0], x2[0], x3[0]);
          set_interior_point(p1);
          ++e1;
          ++e2;
          ++f;
          ++n;

          while(e1) {
               auto q1 = *e1;
               auto l0 = hpreal(0.0), l1 = hpreal(0.0), l2 = hpreal(0.0);
               auto q = 0u, j = 0u;
               std::tie(q, j) = *f;
               std::tie(l0, l1, l2) = get_transition(q, j);
               auto p2 = Point4D(x0[n], x1[n], x2[n], x3[n]);
               auto p3 = l0 * q1 + l1 * p2 + l2 * p1;
               set_boundary_point(q, j, p2);
               set_interior_point(p3);
               p1 = p3;
               ++e1;
               ++e2;
               ++f;
               ++n;
          }

          auto l0 = hpreal(0.0), l1 = hpreal(0.0), l2 = hpreal(0.0);
          std::tie(l0, l1, l2) = get_transition(p, i);
          if(std::abs(l1) < epsilon) set_boundary_point(p, i, get_boundary_point<degree>(std::begin(patches), p, i, 1));
          else {
               auto p3 = Point4D(x0[0], x1[0], x2[0], x3[0]);
               auto p2 = (p3 - l0 * qb - l2 * p1) / l1;
               set_boundary_point(p, i, p2);
          }
     };

     visit_vertices(neighbors, [&](auto p, auto i) {
          auto valence = make_valence(make_spokes_enumerator(neighbors, p, i));
          auto& center = get_corner<degree>(std::begin(patches), p, i);
          surface1.setCorner(p, i, center);
          visit_spokes(make_spokes_enumerator(neighbors, p, i), [&](auto q, auto j) { surface1.setCorner(q, j, p, i); });
          auto coefficients1 = make_coefficients_1(p, i, valence, center);
          set_ring_1(p, i, valence, center, coefficients1);
          auto coefficients2 = make_coefficients_2(p, i, valence);
          set_ring_2(p, i, valence, coefficients2);
     });

     assert(degree == 5);//TODO: update edges and copy interior points for degrees > 5

     visit_edges(neighbors, [&](auto p, auto i) {
          static constexpr hpindex o[3] = { 8, 13, 12 };

          auto q = make_neighbor_index(neighbors, p, i);
          auto j = make_neighbor_offset(neighbors, q, p);
          auto patch0 = get_patch<degree>(std::begin(patches), p);
          auto patch1 = get_patch<degree>(std::begin(patches), q);
          auto e = make_diamonds_enumerator(degree, i, j, [&](auto k0, auto k1, auto k2, auto k3) { return std::tie(patch0[k0], patch1[k1], patch0[k2], patch0[k3]); });
          auto diamond = *(++(++e));
          auto& p0 = std::get<0>(diamond);
          auto& p1 = std::get<1>(diamond);
          auto& p2 = std::get<2>(diamond);
          auto& p3 = std::get<3>(diamond);
          auto l0 = hpreal(0.0), l1 = hpreal(0.0), l2 = hpreal(0.0);
          std::tie(l0, l1, l2) = get_triplet(transitions, 3 * q + j);
          auto x1 = (p1 + l2 * (p3 - l0 * p0 - l1 * p2)) / (hpreal(1.0) + l2 * l2);
          auto x3 = l0 * p0 + l1 * p2 + l2 * x1;
          surface1.setControlPoint(p, o[i], x3);
          surface1.setControlPoint(q, o[j], x1);
     });

     return surface1;
}

template<hpuint degree>
SurfaceSplineBEZ<Space4D, degree> smooth(const SurfaceSplineBEZ<Space3D, degree>& surface, const std::vector<hpreal>& transitions, hpreal epsilon) {
     auto temp = embed<Space4D>(surface, [](const Point3D& point) { return Point4D(point.x, point.y, point.z, 1.0); });
     return smooth(temp, transitions, epsilon);
}

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions) {
     using Point = typename Space::POINT;

     if(nSubdivisions == 0) return surface;

     auto nPatches = size(surface);

     std::vector<SurfaceSubdividerBEZ<Space, degree> > subdividers;
     subdividers.reserve(nPatches);
     visit_patches(surface, [&](auto patch) { subdividers.emplace_back(patch); });

     std::vector<Point> points;
     Indices indices;

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
void visit_boundary(Iterator patch, hpuint i, Visitor&& visit) {
     if(i == 0u) for(auto end = patch + (degree - 1u); patch != end; ) visit(*(++patch));
     else if(i == 1u) {
          auto delta = degree;
          for(patch += degree << 1; delta > 1u; patch += --delta) visit(*patch);
     } else {
          auto delta = 2u;
          for(patch += make_patch_size(degree) - 3u; delta <= degree; patch -= ++delta) visit(*patch);
     }
}

template<hpuint degree, class Iterator, class Visitor>
void visit_boundary(Iterator patches, hpuint p, hpuint i, Visitor&& visit) { visit_boundary<degree>(get_patch<degree>(patches, p), i, std::forward<Visitor>(visit)); }

template<hpuint degree, class Iterator, class Visitor>
void visit_corner(Iterator patch, hpuint i, Visitor&& visit) { visit(get_corner<degree>(patch, i)); }

template<hpuint degree, class Iterator, class Visitor>
void visit_corner(Iterator patches, hpuint p, hpuint i, Visitor&& visit) { visit(get_corner<degree>(patches, p, i)); }

template<hpuint degree, class Iterator, class Visitor>
void visit_corners(Iterator patch, Visitor&& visit) { visit(get_corner<degree>(patch, 0), get_corner<degree>(patch, 1), get_corner<degree>(patch, 2)); }

template<hpuint degree, class Iterator, class Visitor>
void visit_corners(Iterator patches, hpuint p, Visitor&& visit) { visit_corners<degree>(get_patch<degree>(patches, p), std::forward<Visitor>(visit)); }

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

template<class Space, hpuint degree, class Visitor>
void visit_edges(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_edges(neighbors, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_ends(Iterator patch, hpuint i, Visitor&& visit) {
     visit_corners<degree>(patch, [&](auto& c0, auto& c1, auto& c2) {
          if(i == 0u) visit(c0, c1);
          else if(i == 1u) visit(c1, c2);
          else visit(c2, c0);
     });
}

template<hpuint degree, class Iterator, class Visitor>
void visit_ends(Iterator patches, hpuint p, hpuint i, Visitor&& visit) { visit_ends<degree>(get_patch<degree>(patches, p), i, std::forward<Visitor>(visit)); }

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

template<hpuint degree, class Iterator, class Visitor>
void visit_interior(Iterator patches, hpuint p, Visitor&& visit) { visit_interior<degree>(get_patch<degree>(patches, p), std::forward<Visitor>(visit)); }

template<class Iterator, class Visitor>
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
}

template<hpuint degree, class Iterator, class Visitor>
void visit_patch(Iterator patches, hpuint p, Visitor&& visit) { visit(get_patch<degree>(patches, p)); }

template<class Space, hpuint degree, class Visitor>
void visit_patch(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, Visitor&& visit) { visit(get_patch<degree>(surface, p)); }

template<hpuint degree, class Iterator, class Visitor>
void visit_patches(Iterator patches, hpuint nPatches, Visitor&& visit) {
     static constexpr auto patchSize = make_patch_size(degree);
     for(auto i = patches, end = patches + nPatches * patchSize; i != end; i += patchSize) visit(i);
}

template<class Space, hpuint degree, class Visitor>
void visit_patches(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto patches = deindex(surface.getPatches());
     visit_patches<degree>(std::begin(patches), size(surface), std::forward<Visitor>(visit));
}

template<hpindex ring, class Visitor>
void visit_ring(ssb::RingEnumerator<ring> e, Visitor&& visit) { do apply(visit, *e); while(++e); }

template<hpindex ring, class Transformer, class Visitor>
void visit_ring(EnumeratorTransformer<ssb::RingEnumerator<ring>, Transformer> e, Visitor&& visit) { do apply(visit, *e); while(++e); }

template<class Visitor>
void visit_subring(ssb::RingEnumerator<1> e, hpindex p, Visitor&& visit) {
     while(std::get<0>(*e) != p) {
          apply(visit, *e);
          ++e;
     }
     apply(visit, *e);
}

//WORKSPACE

namespace mdz {

//Returns quadratics of the form $$ax_jx_k+bx_j'+c=0$$, where the quadratic coefficients are stored in the first vector, the linear coefficients in the second, and the constant coefficients in the third.  The first integer (the 'i') in each entry of the three vectors identifies the constraint.
std::tuple<std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const Indices& neighbors);

template<class Space, hpuint degree>
std::tuple<std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const SurfaceSplineBEZ<Space, degree>& surface) {
     auto neighbors = make_neighbors(surface);
     return make_constraints(neighbors);
}

//Returns an objective of the form |Ax - b| that is to be minimized.  There are 3 * 9 * nPatches variables.
template<hpuint degree>
std::tuple<std::vector<hpijr>, std::vector<hpir> > make_objective(const SurfaceSplineBEZ<Space3D, degree>& surface) {
     using Point = Point3D;
     auto patches = deindex(surface.getPatches());
     auto neighbors = make_neighbors(surface);
     auto nEdges = 3 * size(surface) / 2;
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
          static constexpr hpuint o1[3] = { 2u, 0u, 1u };

          auto q = make_neighbor_index(neighbors, p, i);
          auto j = make_neighbor_offset(neighbors, q, p);
          auto op = 27 * p + 9 * i;
          auto oq = 27 * q + 9 * j;
          auto b = std::array<Point, 3>();
          auto c = std::array<Point, 3>();
          auto& b0 = get_corner<degree>(std::begin(patches), q, j);
          auto& b1 = b[1];
          auto& b2 = b[0];
          auto& b3 = b[2];
          auto& c0 = c[1];
          auto& c1 = get_corner<degree>(std::begin(patches), p, i);
          auto& c2 = c[2];
          auto& c3 = c[0];
          auto tb = b.data() - 1;
          auto tc = c.data() - 1;

          //TODO: simplify this method using ring enumerator and write directly into b0, b1, b2, b3
          visit_subring(make_ring_enumerator(degree, neighbors, p, o[i]), make_neighbor_index(neighbors, q, o1[j]), [&](auto r, auto k) { *(++tb) = get_patch<degree>(std::begin(patches), r)[k]; });
          visit_subring(make_ring_enumerator(degree, neighbors, q, o[j]), make_neighbor_index(neighbors, p, o1[i]), [&](auto r, auto k) { *(++tc) = get_patch<degree>(std::begin(patches), r)[k]; });

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
std::vector<hpreal> make_transitions(const SurfaceSplineBEZ<Space3D, degree>& surface, const std::vector<hpreal>& solution) {
     static_assert(degree > 4u, "The first two rings of control points surrounding the corners of the patches are assumed to be disjoint.");

     auto transitions = std::vector<hpreal>(9 * size(surface));
     auto patches = deindex(surface.getPatches());
     auto neighbors = make_neighbors(surface);

     auto insert = [&](auto b0, auto b1, auto b2, auto b3, auto p, auto i) {
          auto transition = glm::inverse(hpmat3x3(b0, b1, b2)) * b3;
          auto t = std::begin(transitions) + (9 * p  + 3 * i);
          t[0] = transition.x;
          t[1] = transition.y;
          t[2] = transition.z;
     };

     auto make_point = [&](auto p, auto i) {
          auto& corner = get_corner<degree>(std::begin(patches), p, i);
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
          
          auto q = make_neighbor_index(neighbors, p, i);
          auto j = make_neighbor_offset(neighbors, q, p);
          auto r = make_neighbor_index(neighbors, p, o0[i]);
          auto k = make_neighbor_offset(neighbors, r, p);
          auto& b0 = get_corner<degree>(std::begin(patches), p, o0[i]);
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
std::tuple<std::vector<hpijklr>, std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const Indices& neighbors);

template<class Space, hpuint degree>
std::tuple<std::vector<hpijklr>, std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const SurfaceSplineBEZ<Space, degree>& surface) {
     auto neighbors = make_neighbors(surface);
     return make_constraints(neighbors);
}

//Returns an objective of the form |Ax - b| that is to be minimized.  There are 3 * 12 * nPatches variables.
template<hpuint degree>
std::tuple<std::vector<hpijr>, std::vector<hpir> > make_objective(const SurfaceSplineBEZ<Space3D, degree>& surface) {
     using Point = Point3D;
     auto patches = deindex(surface.getPatches());
     auto neighbors = make_neighbors(surface);
     auto nEdges = 3 * size(surface) / 2;
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

          auto q = make_neighbor_index(neighbors, p, i);
          auto j = make_neighbor_offset(neighbors, q, p);
          auto b = std::array<Point, 3>();
          auto c = std::array<Point, 3>();
          auto& b0 = get_corner<degree>(std::begin(patches), q, j);
          auto& b1 = b[1];
          auto& b2 = b[0];
          auto& b3 = b[2];
          auto& c0 = c[1];
          auto& c1 = get_corner<degree>(std::begin(patches), p, i);
          auto& c2 = c[2];
          auto& c3 = c[0];
          auto tb = b.data() - 1;
          auto tc = c.data() - 1;

          //TODO: simplify this method using ring enumerator and write directly into b0, b1, b2, b3
          visit_subring(make_ring_enumerator(degree, neighbors, p, o[i]), make_neighbor_index(neighbors, q, o1[j]), [&](auto r, auto k) { *(++tb) = get_patch<degree>(std::begin(patches), r)[k]; });
          visit_subring(make_ring_enumerator(degree, neighbors, q, o[j]), make_neighbor_index(neighbors, p, o1[i]), [&](auto r, auto k) { *(++tc) = get_patch<degree>(std::begin(patches), r)[k]; });

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

