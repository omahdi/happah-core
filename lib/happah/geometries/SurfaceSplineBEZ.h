// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

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
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "happah/Eigen.h"
#include "happah/Happah.h"
#include "happah/io/readers/ReaderHPH.h"
#include "happah/io/writers/WriterHPH.h"
#include "happah/geometries/Curve.h"
#include "happah/geometries/Surface.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/geometries/TriangleMeshUtils.h"
#include "happah/utils/DeindexedArray.h"
#include "happah/utils/SurfaceSubdividerBEZ.h"
#include "happah/utils/SurfaceUtilsBEZ.h"
#include "happah/utils/VertexFactory.h"

namespace happah {

//DECLARATIONS

template<class Space, hpuint t_degree>
class SurfaceSplineBEZ;

template<class Transformer>
class DiamondsEnumerator;

namespace ssb {

template<class Transformer>
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

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
TriangleMesh<Vertex> make_control_polygon(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory());

template<class Transformer>
DiamondsEnumerator<Transformer> make_diamonds_enumerator(hpuint degree, hpuint i, hpuint j, Transformer&& transform);

template<int dummy = 0>
auto make_diamonds_enumerator(hpuint degree, hpuint i, hpuint j);

//Return the absolute offset of the ith point in the interior.
hpuint make_interior_offset(hpuint degree, hpuint i);

template<class Space, hpuint degree>
std::vector<hpuint> make_neighbors(const SurfaceSplineBEZ<Space, degree>& surface);

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_ring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i);

template<class Transformer>
ssb::RingEnumerator<Transformer> make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i, Transformer&& transform);

template<int dummy = 0>
auto make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i);

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> make_spline_surface(const std::string& path);

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions, VertexFactory&& factory = VertexFactory());

template<hpuint degree, class Iterator, class Visitor>
void sample(Iterator patches, hpuint nPatches, hpuint nSamples, Visitor&& visit);

template<hpuint degree, class ControlPointsIterator, class DomainPointsIterator, class Visitor>
void sample(ControlPointsIterator controlPoints, DomainPointsIterator domainPoints, hpuint nPatches, hpuint nSamples, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSamples, Visitor&& visit);

template<class Space, hpuint degree, class T, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, std::tuple<const std::vector<T>&, const Indices&> domain, hpuint nSamples, Visitor&& visit);

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions);

bool validate_projective_structure(const Indices& neighbors, const std::vector<hpreal>& transitions, hpreal epsilon = EPSILON);

template<class Space, hpuint degree>
bool validate_projective_structure(const SurfaceSplineBEZ<Space, degree>& surface, const std::vector<hpreal>& transitions, hpreal epsilon = EPSILON);

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
template<class Iterator, class Visitor>
void visit_deltas(hpuint degree, Iterator patch, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_edges(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_ends(Iterator patch, hpuint i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_ends(Iterator patches, hpuint p, hpuint i, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_fans(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit);

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

template<hpuint degree, class Iterator, class Visitor>
void visit_ring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint i, Visitor&& visit);

template<hpuint degree, class Iterator, class Visitor>
void visit_rings(Iterator patches, hpuint nPatches, const Indices& neighbors, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit);

//Visit the subring starting at patch p rotating counterclockwise and stopping at patch q.
template<hpuint degree, class Iterator, class Visitor>
void visit_subring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i, hpuint q, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint q, Visitor&& visit);

template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint q, Visitor&& visit);

//DEFINITIONS

template<class Space, hpuint t_degree>
class SurfaceSplineBEZ : public Surface<Space> {
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

     auto& getBoundaryPoint(hpuint p, hpuint i, hpuint k) const {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          static constexpr auto patchSize = make_patch_size(t_degree);
          i = make_boundary_offset(t_degree, i, k);
          return m_controlPoints[m_indices[p * patchSize + i]];
     }

     auto& getControlPoints() const { return m_controlPoints; }

     auto getNumberOfPatches() const { return m_indices.size() / make_patch_size(t_degree); }

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
     friend auto& operator<<(Stream& stream, const SurfaceSplineBEZ<Space, t_degree>& surface) {
          stream << surface.m_controlPoints << '\n';
          stream << surface.m_indices;
          return stream;
     }

     template<class Stream>
     friend auto& operator>>(Stream& stream, SurfaceSplineBEZ<Space, t_degree>& surface) {
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

/*
 *   k3
 * k0  k2
 *   k1
 * i is with respect to top patch; j is with respect to bottom patch.
 */
template<class Transformer>
class DiamondsEnumerator {
public:
     DiamondsEnumerator(hpuint degree, hpuint i, hpuint j, Transformer transform)
          : m_i(i), m_j(j), m_transform(std::move(transform)) {
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

     explicit operator bool() { return m_k0 != m_end; }

     auto operator*() { return m_transform(m_k0, m_k1, m_k2, m_k3); }

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
     Transformer m_transform;

};

namespace ssb {

template<class Transformer>
class RingEnumerator {
public:
     RingEnumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i, Transformer transform)
          : m_degree(degree), m_e(neighbors, p, i), m_flag(true), m_o0{ 1u, degree << 1, make_patch_size(degree) - 3u }, m_o1{ degree + 1u, degree - 1u, make_patch_size(degree) - 2u }, m_transform(std::move(transform)) {}

     explicit operator bool() { return m_flag; }

     auto operator*() {
          if(m_e) {
               auto p = 0u, i = 0u;
               std::tie(p, i) = *m_e;
               return m_transform(p, m_o0[i]);
          } else return m_transform(m_p, m_o1[m_i]);
     }

     auto& operator++() {
          std::tie(m_p, m_i) = *m_e;
          m_flag = !!m_e;
          ++m_e;
          m_flag &= !!m_e || (m_flag && !m_e && std::get<0>(*m_e) == UNULL);
          return *this;
     }

private:
     hpuint m_degree;
     FanEnumerator<Format::SIMPLE> m_e;
     bool m_flag;
     hpuint m_i;
     hpindex m_o0[3];
     hpindex m_o1[3];
     hpuint m_p;
     Transformer m_transform;

};//RingEnumerator

}//namespace ssb

template<hpuint degree, class Iterator>
auto de_casteljau(Iterator patch, hpreal u, hpreal v, hpreal w) {
     using T = typename std::iterator_traits<Iterator>::value_type;

     if(degree == 0u) return patch[0];

     auto points = std::array<T, make_patch_size(degree - 1u)>();

     auto do_de_casteljau = [&](auto i, auto patch) {
          auto p = std::begin(points) - 1;
          visit_deltas(i, patch, [&](auto& b0, auto& b1, auto& b2) { (++p)[0] = u * b0 + v * b1 + w * b2; });
     };

     do_de_casteljau(degree, patch);
     if(degree > 1u) for(auto i : boost::irange(1u, degree) | boost::adaptors::reversed) do_de_casteljau(i, std::begin(points));
     return points[0];
}

template<class Space, hpuint degree>
auto de_casteljau(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpreal u, hpreal v) { return de_casteljau<degree>(get_patch(surface, p), u, v, 1.0f - u - v); }

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, (degree + 1)> elevate(const SurfaceSplineBEZ<Space, degree>& surface) {
     SurfaceSplineBEZ<Space, (degree + 1)> surface1(surface.getNumberOfPatches());
     auto& indices = std::get<1>(surface.getPatches());
     auto neighbors = make_neighbors(surface);
     auto patches = deindex(surface.getPatches());

     visit_vertices(neighbors, [&](auto p, auto i) {
          surface1.setCorner(p, i, get_corner<degree>(std::begin(patches), p, i));
          visit_fan(neighbors, p, i, [&](auto q, auto j) { surface1.setCorner(q, j, p, i); });
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
     visit_patches<degree>(std::begin(patches), surface.getNumberOfPatches(), [&](auto patch) {
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
          temp[make_offset(degree, i, j, k)].z = 1.0;
          auto surface = SurfaceSplineBEZ<Space3D, degree>(temp);
          std::ostringstream path;
          path << directory << "/b" << i << j << k << ".ss" << degree << ".bz.3.hph";
          WriterHPH::write(surface, path.str().c_str());
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

template<class Space, hpuint degree, class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_control_polygon(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory) { return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, 0, std::forward<VertexFactory>(factory)); }

template<class Transformer>
DiamondsEnumerator<Transformer> make_diamonds_enumerator(hpuint degree, hpuint i, hpuint j, Transformer&& transform) { return { degree, i, j, transform }; }

template<int dummy>
auto make_diamonds_enumerator(hpuint degree, hpuint i, hpuint j) { return make_diamonds_enumerator(degree, i, j, [](auto k0, auto k1, auto k2, auto k3) { return std::make_tuple(k0, k1, k2, k3); }); }

template<class Space, hpuint degree>
std::vector<hpuint> make_neighbors(const SurfaceSplineBEZ<Space, degree>& surface) {
     auto indices = Indices();
     indices.reserve(3 * surface.getNumberOfPatches());
     auto inserter = make_back_inserter(indices);
     visit_patches<degree>(std::begin(std::get<1>(surface.getPatches())), surface.getNumberOfPatches(), [&](auto patch) { visit_corners<degree>(patch, inserter); });
     return make_neighbors(indices);
}

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_ring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i) {
     auto ring = std::vector<typename std::iterator_traits<Iterator>::value_type>();
     visit_ring<degree>(patches, neighbors, p, i, make_back_inserter(ring));
     return ring;
}

template<class Transformer>
ssb::RingEnumerator<Transformer> make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i, Transformer&& transform) { return { degree, neighbors, p, i, transform }; }

template<int dummy>
auto make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i) { return make_ring_enumerator(degree, neighbors, p, i, [&](auto q, auto j) { return std::make_tuple(q, j); }); }

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> make_spline_surface(const std::string& path) { return ReaderHPH::read<SurfaceSplineBEZ<Space, degree> >(path); }

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
          indices.reserve(3 * make_control_polygon_size(degree) * surface.getNumberOfPatches());
          auto inserter = make_back_inserter(indices);
          visit_patches<degree>(std::begin(std::get<1>(surface.getPatches())), surface.getNumberOfPatches(), [&](auto patch) {
               visit_deltas(degree, patch, inserter);
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
     sample<degree>(std::begin(patches), surface.getNumberOfPatches(), nSamples, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class T, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, std::tuple<const std::vector<T>&, const Indices&> domain, hpuint nSamples, Visitor&& visit) {
     auto controlPoints = deindex(surface.getPatches());
     auto domainPoints = deindex(domain);
     sample<degree>(std::begin(controlPoints), std::begin(domainPoints), surface.getNumberOfPatches(), nSamples, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions) {
     using Point = typename Space::POINT;

     if(nSubdivisions == 0) return surface;

     auto nPatches = surface.getNumberOfPatches();

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

template<class Space, hpuint degree>
bool validate_projective_structure(const SurfaceSplineBEZ<Space, degree>& surface, const std::vector<hpreal>& transitions, hpreal epsilon) { return is_valid_projective_structure(make_neighbors(surface), transitions, epsilon); }

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

template<class Iterator, class Visitor>
void visit_deltas(hpuint degree, Iterator patch, Visitor&& visit) {
     auto bottom = patch;
     auto top = patch + (degree + 1u);
     auto delta = degree;
     while(delta > 0u) {
          for(auto end = bottom + delta; bottom != end; ++bottom, ++top) visit(bottom[0], bottom[1], top[0]);
          --delta;
          ++bottom;
     }
}

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

template<class Space, hpuint degree, class Visitor>
void visit_fans(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_fans(neighbors, std::forward<Visitor>(visit));
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
     visit_patches<degree>(std::begin(patches), surface.getNumberOfPatches(), std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_ring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit) { for(auto e = make_ring_enumerator(degree, neighbors, p, i, [&](auto q, auto j) { return patches[make_patch_size(degree) * q + j]; }); e; ++e) visit(*e); }

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit) {
     auto patches = deindex(surface.getPatches());
     visit_ring<degree>(std::begin(patches), neighbors, p, i, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint i, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_ring(surface, neighbors, p, i, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_rings(Iterator patches, hpuint nPatches, const Indices& neighbors, Visitor&& visit) { visit_vertices(neighbors, [&](auto p, auto i) { visit(p, i, make_ring_enumerator(degree, neighbors, p, i, [&](auto q, auto j) { return patches[make_patch_size(degree) * q + j]; })); }); }

template<class Space, hpuint degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, Visitor&& visit) {
     auto patches = deindex(surface.getPatches());
     visit_rings<degree>(std::begin(patches), surface.getNumberOfPatches(), neighbors, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_rings(surface, neighbors, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_subring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i, hpuint q, Visitor&& visit) {
     hpuint k;
     visit_subfan(neighbors, p, i, q, [&](hpuint r, hpuint j) {
          k = j;
          visit_patch<degree>(patches, r, [&](auto patch) {
               static constexpr auto patchSize = make_patch_size(degree);
               if(j == 0) visit(*(patch + 1));
               else if(j == 1) visit(*(patch + + (degree << 1)));
               else visit(*(patch + (patchSize - 3)));
          });
     });
     visit_patch<degree>(patches, q, [&](auto patch) {
          static constexpr auto patchSize = make_patch_size(degree);
          if(k == 0) visit(*(patch + (degree + 1)));
          else if(k == 1) visit(*(patch + (degree - 1)));
          else visit(*(patch + (patchSize - 2)));
     });
}

template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint q, Visitor&& visit) {
     auto& indices = std::get<1>(surface.getPatches());
     auto patches = deindex(surface.getPatches());
     hpuint i;
     visit_triplet(neighbors, p, [&](hpuint n0, hpuint n1, hpuint n2) { i = (q == n0) ? 1 : (q == n1) ? 2 : (q == n2) ? 0 : UNULL; });
     if(i == UNULL) visit_corners<degree>(std::begin(indices), p, [&](hpuint v0, hpuint v1, hpuint v2) {
          visit_corners<degree>(std::begin(indices), q, [&](hpuint w0, hpuint w1, hpuint w2) { i = (v0 == w0 || v0 == w1 || v0 == w2) ? 0 : (v1 == w0 || v1 == w1 || v1 == w2) ? 1 : 2; });
     });
     visit_subring<degree>(std::begin(patches), neighbors, p, i, q, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint q, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_subring(surface, neighbors, p, q, std::forward<Visitor>(visit));
}

//WORKSPACE

namespace tpfssb {

template<hpuint degree>
std::tuple<std::vector<hpijr>, std::vector<hpir> > make_objective(const SurfaceSplineBEZ<Space3D, degree>& surface) {
     auto patches = surface.getPatches();
     auto& indices = std::get<1>(patches);
     auto& points = std::get<0>(patches);
     auto neighbors = make_neighbors(surface);

     std::unordered_map<hpuint, Point3D> coordinates;

     visit_rings<degree>(std::begin(indices), surface.getNumberOfPatches(), neighbors, [&](auto p, auto i, auto begin, auto end) {
          hpuint v;
          visit_corner<degree>(std::begin(indices), p, i, [&](auto& w) { v = w; });

          auto& c0 = points[v];
          auto& c1 = points[begin[0]];
          auto& c2 = points[begin[1]];

          auto A = glm::inverse(hpmat3x3(c0, c1, c2));

          for(auto& q : boost::make_iterator_range(begin, end)) coordinates[q] = A * points[q];
     });

     std::vector<hpir> hirs;
     std::vector<hpijr> hijrs;

     auto row = 0u;
     visit_edges(neighbors, [&](hpuint p, hpuint i) {
          hpuint b[3];
          hpuint c[3];
          auto q = make_neighbor_index(neighbors, p, i);
          auto j = make_neighbor_offset(neighbors, q, p);

          auto tb = b + 3;
          auto tc = c - 1;
          visit_subring<degree>(std::begin(indices), neighbors, q, (j == 0) ? 1 : (j == 1) ? 2 : 0, p, [&](auto k) { *(--tb) = k; });
          visit_subring<degree>(std::begin(indices), neighbors, p, (i == 0) ? 1 : (i == 1) ? 2 : 0, q, [&](auto k) { *(++tc) = k; });

          auto e0 = Point3D(1, 0, 0);
          auto& b0 = coordinates[b[0]];
          auto& b1 = coordinates[b[1]];
          auto& b2 = coordinates[b[2]];
          auto& c0 = coordinates[c[0]];
          auto& c1 = coordinates[c[1]];
          auto& c2 = coordinates[c[2]];

          // indexing of matrix elements:
          //   x0 x1 x2
          //   x3 x4 x5
          //   x6 x7 x8

          auto insert = [&](auto& source, auto& target, hpuint offset) {
               hirs.emplace_back(row, target.x);
               hijrs.emplace_back(row, offset + 0, source.x);
               hijrs.emplace_back(row, offset + 1, source.y);
               hijrs.emplace_back(row, offset + 2, source.z);
               ++row;
               hirs.emplace_back(row, target.y);
               hijrs.emplace_back(row, offset + 3, source.x);
               hijrs.emplace_back(row, offset + 4, source.y);
               hijrs.emplace_back(row, offset + 5, source.z);
               ++row;
               hirs.emplace_back(row, target.z);
               hijrs.emplace_back(row, offset + 6, source.x);
               hijrs.emplace_back(row, offset + 7, source.y);
               hijrs.emplace_back(row, offset + 8, source.z);
               ++row;
          };

          auto sp = 27 * p + 9 * i;
          insert(b0, c0, sp);
          insert(b1, e0, sp);
          insert(b2, c2, sp);
          insert(e0, c1, sp);

          auto sq = 27 * q + 9 * j;
          insert(c0, b0, sq);
          insert(c1, e0, sq);
          insert(c2, b2, sq);
          insert(e0, b1, sq);
     });

     return std::make_tuple(std::move(hijrs), std::move(hirs));
}

/**
 *  The constraints are quadratics of the form $$ax_jx_k+bx_j=c$$, where the quadratic coefficients are stored in the first vector, the linear coefficients in the second, and the constant coefficients in the third.  The first integer (the 'i') in each entry of the three vectors identifies the constraint equality.
 */

template<class Space, hpuint degree>
std::tuple<std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const SurfaceSplineBEZ<Space, degree>& surface) {
     std::vector<hpijkr> ijkrs;
     std::vector<hpijr> ijrs;
     std::vector<hpir> irs;

     hpuint c = 0; //Counter for the constraint number.
     auto neighbors = make_neighbors(surface);

     //Identity constraint: `p . q = id`, where `p` and `q` are opposite projections on the same edge.
     visit_edges(neighbors, [&](hpuint p, hpuint i) {
               //`q` is CW!! in this case
               auto q = *(neighbors.begin() + 3*p + i);

               hpuint j;
               visit_triplet(neighbors, q, [&](hpuint n0, hpuint n1, hpuint n2) { j = (p == n0) ? 0 : (p == n1) ? 1 : 2; });

               /**
                * 3x3 matrix entries indices
                * Encoding: every patch identifies 3 3x3 matrices.
                * Use the position of neighbor index to identify each of the sets of 9 entries
                * These are the indices into the huge vector to optimise */
               std::vector<hpuint> p_v(9);
               std::iota(p_v.begin(), p_v.end(), 27*p + (9*i)); //+ [0..9]

               std::vector<hpuint> q_v(9);
               std::iota(q_v.begin(), q_v.end(), 27*q + (9*j));

               //Iterate over the matrix p . q - id
               for(hpuint u = 0; u < 3; u++) {
                    for(hpuint v = 0; v < 3; v++) {
                         //Rolled out (index) matrix multiplication
                         ijkrs.emplace_back(c, p_v[3*u + 0], q_v[v + 3*0], 1);
                         ijkrs.emplace_back(c, p_v[3*u + 1], q_v[v + 3*1], 1);
                         ijkrs.emplace_back(c, p_v[3*u + 2], q_v[v + 3*2], 1);

                         if(u == v)
                              irs.emplace_back(c, -1);
                         c++;
                    }
               }

          });//visit_edges, identity constraint


     //Compatibility constraint: `p_1 . p_2 . p_3 = id`, where `p_i` are the three projections inside a given patch. Note that this can be rewritten as `p_1 . p_2 = p_3^(-1) `where the inverse denotes the projection running opposite the edge (in the code this is `q`).
     for(auto p = 0lu; p < surface.getNumberOfPatches(); p++) {
          //`q` is CW!! in this case
          auto q = *(neighbors.begin() + 3*p + 2); //Choose an arbitrary neighbor.

          hpuint j;
          visit_triplet(neighbors, q, [&](hpuint n0, hpuint n1, hpuint n2) { j = (p == n0) ? 0 : (p == n1) ? 1 : 2; });

          //3x3 matrix entries indices
          //Encoding: every patch identifies 3 3x3 matrices.
          //Use the position of neighbor index to identify each of the sets of 9 entries
          std::vector<hpuint> p1_v(9);
          std::iota(p1_v.begin(), p1_v.end(), 27*p + (9*0)); //+ [0..9]

          std::vector<hpuint> p2_v(9);
          std::iota(p2_v.begin(), p2_v.end(), 27*p + (9*1));

          std::vector<hpuint> q_v(9);
          std::iota(q_v.begin(), q_v.end(), 27*q + (9*j));

          //p_1 . p_2 - q
          for(hpuint u = 0; u < 3; u++) {
               for(hpuint v = 0; v < 3; v++) {
                    //Rolled out (index) matrix multiplication
                    ijkrs.emplace_back(c, p1_v[3*u + 0], p2_v[v + 3*0], 1);
                    ijkrs.emplace_back(c, p1_v[3*u + 1], p2_v[v + 3*1], 1);
                    ijkrs.emplace_back(c, p1_v[3*u + 2], p2_v[v + 3*2], 1);
                    ijrs.emplace_back(c, q_v[3*u + v], -1);
                    c++;
               }
          }
     }//for(numOfPatches), compatibility constraint

     return std::make_tuple(std::move(ijkrs), std::move(ijrs), std::move(irs));
}

}//namespace tpfssb

namespace phm {

//Returns cubics of the form $$ax_jx_kx_l+bx_j'x_k'+cx_j''+d=0$$, where the cubic coefficients are stored in the first vector, the quadratic coefficients in the second, the linear coefficients in the third, and the constant coefficients in the fourth.  The first integer (the 'i') in each entry of the four vectors identifies the constraint.
template<class Space, hpuint degree>
std::tuple<std::vector<hpijklr>, std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const SurfaceSplineBEZ<Space, degree>& surface) {
     auto neighbors = make_neighbors(surface);
     auto nPatches = surface.getNumberOfPatches();
     auto nEdges = 3 * nPatches / 2;
     auto hpirs = std::vector<hpir>();
     auto hpijrs = std::vector<hpijr>();
     auto hpijkrs = std::vector<hpijkr>();
     auto hpijklrs = std::vector<hpijklr>();
     auto row = -1;

     //TODO: reserve hpi...s

     // indexing of rho: (see make_objective)
     // indexing of lambda: (see make_objective)

     auto insert = [&](auto& x, auto& y) {
          auto do_row = [&](auto x0, auto x1, auto x2, auto y0, auto y3, auto y6) {
               hpijkrs.emplace_back(++row, y0, x1, 1.0);
               hpijkrs.emplace_back(row, y3, x0, 1.0);
               hpijklrs.emplace_back(row, y6, y[11], x0, 1.0);
               hpijklrs.emplace_back(row, y6, y[9], x1, 1.0);
               hpijklrs.emplace_back(row, y6, y[10], x2, 1.0);
          };

          auto do_column = [&](auto y0, auto y3, auto y6) {
               do_row(x[0], x[1], x[2], y0, y3, y6);
               do_row(x[3], x[4], x[5], y0, y3, y6);
               do_row(x[6], x[7], x[8], y0, y3, y6);
          };

          hpijrs.emplace_back(row + 1, y[9], -1.0);
          hpijrs.emplace_back(row + 2, y[10], -1.0);
          hpijrs.emplace_back(row + 3, y[11], -1.0);
          hpirs.emplace_back(row + 4, -1.0);
          hpirs.emplace_back(row + 9, -1.0);

          do_column(y[0], y[3], y[6]);
          do_column(y[1], y[4], y[7]);
          do_column(y[2], y[5], y[8]);
     };

     auto make_array = [](hpuint op) -> auto { return std::array<hpuint, 12>{ op + 6u, op + 7u, op + 8u, op, op + 1u, op + 2u, op + 3u, op + 4u, op + 5u }; };

     visit_edges(neighbors, [&](auto p, auto i) {
          auto q = make_neighbor_index(neighbors, p, i);
          auto j = make_neighbor_offset(neighbors, q, p);
          auto op = 36 * p + 12 * i;
          auto oq = 36 * q + 12 * j;
          auto x = std::array<hpuint, 12>();
          auto y = std::array<hpuint, 12>();

          std::iota(std::begin(x), std::end(x), op);
          std::iota(std::begin(y), std::end(y), oq);

          // lambda * lambda' = id
          hpijkrs.emplace_back(++row, x[11], y[10], 1.0);
          hpijrs.emplace_back(row, y[9], 1.0);
          hpijkrs.emplace_back(++row, x[9], y[10], 1.0);
          hpijrs.emplace_back(row, y[11], 1.0);
          hpijkrs.emplace_back(++row, x[10], y[10], 1.0);
          hpirs.emplace_back(row, -1.0);

          // rho * lambda' * rho' = lambda'
          insert(x, y);

          // rho' * lambda * rho = lambda
          insert(y, x);
     });

     // rho0 * rho1 * rho2 = id
     for(auto p : boost::irange(0lu, nPatches)) {
          auto op = 36 * p;
          auto x = make_array(op);
          auto y = make_array(op + 12);
          auto z = make_array(op + 24);

          auto do_column_x = [&](auto z0, auto y0, auto x0, auto x3, auto x6) {
               hpijklrs.emplace_back(row + 1, z0, y0, x0, 1.0);
               hpijklrs.emplace_back(row + 2, z0, y0, x3, 1.0);
               hpijklrs.emplace_back(row + 3, z0, y0, x6, 1.0);
          };

          auto do_column_y = [&](auto z0, auto y0, auto y3, auto y6) {
               do_column_x(z0, y0, x[0], x[3], x[6]);
               do_column_x(z0, y3, x[1], x[4], x[7]);
               do_column_x(z0, y6, x[2], x[5], x[8]);
          };
          
          auto do_column_z = [&](auto z0, auto z3, auto z6) {
               do_column_y(z0, y[0], y[3], y[6]);
               do_column_y(z3, y[1], y[4], y[7]);
               do_column_y(z6, y[2], y[5], y[8]);
          };

          do_column_z(z[0], z[3], z[6]);
          hpirs.emplace_back(row, -1); 
          row += 3;
          do_column_z(z[1], z[4], z[7]);
          hpirs.emplace_back(row + 1, -1); 
          row += 3;
          do_column_z(z[2], z[5], z[8]);
          hpirs.emplace_back(row + 2, -1); 
          row += 3;
     }

     return std::make_tuple(std::move(hpijklrs), std::move(hpijkrs), std::move(hpijrs), std::move(hpirs));
}

//Returns an objective of the form |Ax - b| that is to be minimized.  There are 3 * 12 * nPatches variables.
template<hpuint degree>
std::tuple<std::vector<hpijr>, std::vector<hpir> > make_objective(const SurfaceSplineBEZ<Space3D, degree>& surface) {
     using Point = Point3D;
     auto patches = deindex(surface.getPatches());
     auto neighbors = make_neighbors(surface);
     auto nEdges = 3 * surface.getNumberOfPatches() / 2;
     auto hirs = std::vector<hpir>();
     auto hijrs = std::vector<hpijr>();
     auto row = -1;

     hirs.reserve(2 * 33 * nEdges);
     hijrs.reserve(2 * 21 * nEdges);

     // indexing of rho:
     //   x0 x1 x2
     //   x3 x4 x5
     //   x6 x7 x8
     // indexing of lambda:
     //   1 0 x9
     //   0 0 x10
     //   0 1 x11

     auto do_column = [&](auto offset, auto x, auto y, auto z) {
          hijrs.emplace_back(++row, offset + 0, 1.0);
          hirs.emplace_back(row, x);
          hijrs.emplace_back(++row, offset + 1, 1.0);
          hirs.emplace_back(row, y);
          hijrs.emplace_back(++row, offset + 2, 1.0);
          hirs.emplace_back(row, z);
     };

     auto do_row = [&](auto offset, auto x, auto y, auto z, auto a) {
          hijrs.emplace_back(++row, offset + 0, x);
          hijrs.emplace_back(row, offset + 1, y);
          hijrs.emplace_back(row, offset + 2, z);
          hirs.emplace_back(row, a);
     };

     auto insert = [&](auto offset, auto b2, auto b3, auto c2, auto c3) {
          // |rho - id|
          hijrs.emplace_back(++row, offset + 0, 1.0);
          hirs.emplace_back(row, 1.0);
          hijrs.emplace_back(++row, offset + 1, 1.0);
          hijrs.emplace_back(++row, offset + 2, 1.0);
          hijrs.emplace_back(++row, offset + 3, 1.0);
          hijrs.emplace_back(++row, offset + 4, 1.0);
          hirs.emplace_back(row, 1.0);
          hijrs.emplace_back(++row, offset + 5, 1.0);
          hijrs.emplace_back(++row, offset + 6, 1.0);
          hijrs.emplace_back(++row, offset + 7, 1.0);
          hijrs.emplace_back(++row, offset + 8, 1.0);
          hirs.emplace_back(row, 1.0);

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

     visit_edges(neighbors, [&](auto p, auto i) {
          static constexpr hpuint o[3] = { 1u, 2u, 0u };

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

          visit_subring<degree>(std::begin(patches), neighbors, p, o[i], q, [&](auto point) { *(++tb) = point; });
          visit_subring<degree>(std::begin(patches), neighbors, q, o[j], p, [&](auto point) { *(++tc) = point; });

          auto Ab2 = A(b0, b3, b1, b2);
          auto Ab3 = A(b0, b1, b2, b3);
          auto Ac2 = A(c0, c3, c1, c2);
          auto Ac3 = A(c0, c1, c2, c3);

          insert(36 * p + 12 * i, Ab2, Ab3, Ac2, Ac3);
          insert(36 * q + 12 * j, Ac3, Ac2, Ab3, Ab2);
     });

     return std::make_tuple(std::move(hijrs), std::move(hirs));
}

template<hpuint degree>
void smooth(SurfaceSplineBEZ<Space4D, degree>& surface, const std::vector<hpreal>& transitions) {
     static_assert(degree > 1u, "Not implemented for linear surfaces.");

     auto neighbors = make_neighbors(surface);
     auto patches = deindex(surface.getPatches());
     auto surface1 = SurfaceSplineBEZ<Space4D, degree>(surface.getNumberOfPatches());

     auto make_coefficients = [&](auto p, auto i, auto valence, auto& center) {
          auto coefficients = std::vector<double>(valence * 7, 0.0);
          auto a0 = std::begin(coefficients) - 1;
          auto a1 = a0 + valence;
          auto a2 = a1 + valence;
          auto b0 = a2 + valence;
          auto b1 = b0 + valence;
          auto b2 = b1 + valence;
          auto b3 = b2 + valence;
          auto e = make_ring_enumerator(degree, neighbors, p, i, [&](auto q, auto j) { return get_patch<degree>(std::begin(patches), q)[j]; });
          auto f = make_fan_enumerator(neighbors, p, i);

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
               auto q = 0u, j = 0u;
               std::tie(q, j) = *f;
               visit_triplet(transitions, 3 * q + j, [&](auto& l0, auto& l1, auto& l2) {
                    auto t0 = l0 + l1 * a0[0] + l2 * (a0 - 1)[0];
                    auto t1 = l1 * a1[0] + l2 * (a1 - 1)[0];
                    auto t2 = l1 * a2[0] + l2 * (a2 - 1)[0];
                    push_back(t0, t1, t2);
               });
               ++f;
          }

          return coefficients;
     };

     auto set_ring = [&](auto p, auto i, auto valence, auto& center, auto& coefficients) {
          auto a0 = std::begin(coefficients) - 1;
          auto a1 = a0 + valence;
          auto a2 = a1 + valence;
          auto A = Eigen::Map<Eigen::MatrixX2d>(coefficients.data() + valence, valence, 2);
          auto temp = A.colPivHouseholderQr();
          auto x0 = temp.solve(Eigen::Map<Eigen::VectorXd>(coefficients.data() + 3 * valence, valence));
          auto x1 = temp.solve(Eigen::Map<Eigen::VectorXd>(coefficients.data() + 4 * valence, valence));
          auto x2 = temp.solve(Eigen::Map<Eigen::VectorXd>(coefficients.data() + 5 * valence, valence));
          auto x3 = temp.solve(Eigen::Map<Eigen::VectorXd>(coefficients.data() + 6 * valence, valence));

          auto q1 = Point4D(x0[0], x1[0], x2[0], x3[0]);
          auto q2 = Point4D(x0[1], x1[1], x2[1], x3[1]);

          auto e = make_fan_enumerator(neighbors, p, i);

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

     visit_vertices(neighbors, [&](auto p, auto i) {
          auto valence = make_valence(neighbors, p, i);
          auto& center = get_corner<degree>(std::begin(patches), p, i);
          auto coefficients = make_coefficients(p, i, valence, center);
          set_ring(p, i, valence, center, coefficients);
          surface1.setCorner(p, i, center);
          visit_fan(neighbors, p, i, [&](auto q, auto j) { surface1.setCorner(q, j, p, i); });
     });
}

}//namespace phm

}//namespace happah

