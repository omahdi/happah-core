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
#include <boost/range/irange.hpp>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/Curve.h"
#include "happah/geometries/Surface.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/geometries/TriangleMeshUtils.h"
#include "happah/utils/DeindexedArray.h"
#include "happah/utils/SurfaceSubdividerBEZ.h"
#include "happah/utils/SurfaceSplineUtilsBEZ.h"
#include "happah/utils/SurfaceUtilsBEZ.h"
#include "happah/utils/VertexFactory.h"

namespace happah {

//DECLARATIONS

template<class Space, hpuint t_degree>
class SurfaceSplineBEZ;

//visitors

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
//template<hpuint degree, class Iterator, class Visitor>
//void visit_deltas(Iterator patch, Visitor&& visit);

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
template<hpuint degree, class Iterator, class Visitor>
void visit_nablas(Iterator patch, Visitor&& visit);

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

//algorithms

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, (degree + 1)> elevate(const SurfaceSplineBEZ<Space, degree>& surface);

template<class Space, hpuint degree>
bool is_c0(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i);

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_boundary(Iterator patch, hpuint i);

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_boundary(Iterator patches, hpuint p, hpuint i);

//Return the absolute offset of the jth point on the ith boundary.
hpuint make_boundary_offset(hpuint degree, hpuint i, hpuint j);

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex> >
TriangleMesh<Vertex> make_control_polygon(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory());

//Return the absolute offset of the ith point in the interior.
hpuint make_interior_offset(hpuint degree, hpuint i);

template<class Space, hpuint degree>
std::vector<hpuint> make_neighbors(const SurfaceSplineBEZ<Space, degree>& surface);

constexpr hpuint make_patch_size(hpuint degree);

template<class Space, hpuint degree, class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, typename = typename std::enable_if<(degree > 0)>::type>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory = VertexFactory());

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
     void setBoundary(hpuint p, hpuint i, Iterator begin) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          auto n = m_controlPoints.size();
          visit_boundary<t_degree>(m_indices.begin(), p, i, [&](auto& i) { i = n++; });
          m_controlPoints.insert(m_controlPoints.end(), begin, begin + (t_degree - 1));
     }

     //Set the ith boundary of the pth patch to the jth boundary of the qth patch.
     void setBoundary(hpuint p, hpuint i, hpuint q, hpuint j) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          auto boundary = make_boundary<t_degree>(m_indices.begin(), q, j);
          auto n = boundary.end();
          visit_boundary<t_degree>(m_indices.begin(), p, i, [&](auto& i) { i = *(--n); });
     }

     //Set the kth point on the ith boundary of the pth patch.
     void setBoundaryPoint(hpuint p, hpuint i, hpuint k, Point point) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          static constexpr auto patchSize = make_patch_size(t_degree);
          i = make_boundary_offset(t_degree, i, k);
          m_indices[p * patchSize + i] = m_controlPoints.size();
          m_controlPoints.push_back(point);
     }

     //Set the kth point on the ith boundary of the pth patch to the opposite point on the jth boundary of the qth patch.
     void setBoundaryPoint(hpuint p, hpuint i, hpuint k, hpuint q, hpuint j) {
          static_assert(t_degree > 1, "There is no boundary in a constant or linear.");
          static constexpr auto patchSize = make_patch_size(t_degree);
          i = make_boundary_offset(t_degree, i, k);
          j = make_boundary_offset(t_degree, j, t_degree - 2 - k);
          m_indices[p * patchSize + i] = m_indices[q * patchSize + j];
     }

     void setCorner(hpuint p, hpuint i, Point point) {
          static_assert(t_degree > 0, "There is no corner in a constant.");
          visit_corner<t_degree>(m_indices.begin(), p, i, [&](auto& i) { i = m_controlPoints.size(); });
          m_controlPoints.push_back(point);
     }

     void setCorner(hpuint p, hpuint i, hpuint q, hpuint j) {
          static_assert(t_degree > 0, "There is no corner in a constant.");
          visit_corner<t_degree>(m_indices.begin(), p, i, [&](auto& k) { visit_corner<t_degree>(m_indices.begin(), q, j, [&](auto& l) { k = l; }); });
     }

     template<class Iterator>
     void setInterior(hpuint p, Iterator begin) {
          static_assert(t_degree > 2, "There is no interior in a constant, linear, or quadratic.");
          auto n = m_controlPoints.size();
          visit_interior<t_degree>(m_indices.begin(), p, [&](auto& i) { i = n++; });
          m_controlPoints.insert(m_controlPoints.end(), begin, begin + (make_patch_size(t_degree) - 3 * t_degree)); 
     }

     void setInteriorPoint(hpuint p, hpuint i, Point point) {
          static_assert(t_degree > 2, "There is no interior in a constant, linear, or quadratic.");
          static constexpr auto patchSize = make_patch_size(t_degree);
          i = make_interior_offset(t_degree, i);
          m_indices[p * patchSize + i] = m_controlPoints.size();
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

//visitors

template<hpuint degree, class Iterator, class Visitor>
void visit_boundary(Iterator patch, hpuint i, Visitor&& visit) {
     switch(i) {
     case 0u: {
          for(auto end = patch + (degree - 1u); patch != end; ) visit(*(++patch));
          break;
     } case 1u: {
          auto delta = degree;
          for(patch += degree << 1; delta > 1u; patch += --delta) visit(*patch);
          break;
     } case 2u: {
          auto delta = 2u;
          for(patch += make_patch_size(degree) - 3u; delta <= degree; patch -= ++delta) visit(*patch);
          break;
     }
     }
}

template<hpuint degree, class Iterator, class Visitor>
void visit_boundary(Iterator patches, hpuint p, hpuint i, Visitor&& visit) { visit_patch<degree>(patches, p, [&](auto patch) { visit_boundary<degree>(patch, i, std::forward<Visitor>(visit)); }); }

template<hpuint degree, class Iterator, class Visitor>
void visit_corner(Iterator patch, hpuint i, Visitor&& visit) {
     visit_corners<degree>(patch, [&](auto& c0, auto& c1, auto& c2) {
          if(i == 0u) visit(c0);
          else if(i == 1u) visit(c1);
          else visit(c2);
     });
}

template<hpuint degree, class Iterator, class Visitor>
void visit_corner(Iterator patches, hpuint p, hpuint i, Visitor&& visit) { visit_patch<degree>(patches, p, [&](auto patch) { visit_corner<degree>(patch, i, std::forward<Visitor>(visit)); }); }

template<hpuint degree, class Iterator, class Visitor>
void visit_corners(Iterator patch, Visitor&& visit) { visit(*patch, *(patch + degree), *(patch + (make_patch_size(degree) - 1))); }

template<hpuint degree, class Iterator, class Visitor>
void visit_corners(Iterator patches, hpuint p, Visitor&& visit) { visit_patch<degree>(patches, p, [&](auto patch) { visit_corners<degree>(patch, std::forward<Visitor>(visit)); }); }

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
void visit_ends(Iterator patches, hpuint p, hpuint i, Visitor&& visit) { visit_patch<degree>(patches, p, [&](auto patch) { visit_ends<degree>(patch, i, std::forward<Visitor>(visit)); }); }

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
void visit_interior(Iterator patches, hpuint p, Visitor&& visit) { visit_patch<degree>(patches, p, [&](auto patch) { visit_interior<degree>(patch, std::forward<Visitor>(visit)); }); }

template<hpuint degree, class Iterator, class Visitor>
void visit_nablas(Iterator patch, Visitor&& visit) {
     auto bottom = patch + 1;
     auto top = patch + (degree + 1u);
     auto delta = degree - 1u;
     while(delta > 0u) {
          for(auto end = bottom + delta; bottom != end; ++bottom) {
               auto& b2 = *bottom;
               auto& b1 = *top;
               auto& b0 = *(++top);
               visit(b0, b1, b2);
          }
          --delta;
          bottom += 2;
          ++top;
     }
}

template<hpuint degree, class Iterator, class Visitor>
void visit_patch(Iterator patches, hpuint p, Visitor&& visit) {
     static constexpr auto patchSize = make_patch_size(degree);
     visit(patches + p * patchSize);
}

template<class Space, hpuint degree, class Visitor>
void visit_patch(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, Visitor&& visit) {
     auto patches = deindex(surface.getPatches());
     visit_patch<degree>(patches.begin(), p, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_patches(Iterator patches, hpuint nPatches, Visitor&& visit) {
     static constexpr auto patchSize = make_patch_size(degree);
     for(auto i = patches, end = patches + nPatches * patchSize; i != end; i += patchSize) visit(i);
}

template<class Space, hpuint degree, class Visitor>
void visit_patches(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto patches = deindex(surface.getPatches());
     visit_patches<degree>(patches.begin(), surface.getNumberOfPatches(), std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_ring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit) {
     static constexpr auto patchSize = make_patch_size(degree);
     hpuint first = UNULL, last, n = 0u, k;

     visit_fan(neighbors, p, i, [&](hpuint q, hpuint j) {
          if(first == UNULL) first = q;
          last = q;
          k = j;
          ++n;
          visit_patch<degree>(patches, q, [&](auto patch) {
               if(j == 0) visit(*(patch + 1));
               else if(j == 1) visit(*(patch + (degree << 1)));
               else visit(*(patch + (patchSize - 3)));
          });
     });
     if(n < 3 || !is_neighbor(neighbors, first, last)) {
          visit_patch<degree>(patches, last, [&](auto patch) {
               if(k == 0) visit(*(patch + (degree + 1)));
               else if(k == 1) visit(*(patch + (degree - 1)));
               else visit(*(patch + (patchSize - 2)));
          });
     }
}

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit) {
     auto patches = deindex(surface.getPatches());
     visit_ring<degree>(patches.begin(), neighbors, p, i, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint i, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_ring(surface, neighbors, p, i, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_rings(Iterator patches, hpuint nPatches, const Indices& neighbors, Visitor&& visit) {
     using T = typename std::iterator_traits<Iterator>::value_type;
     boost::dynamic_bitset<> visited(neighbors.size(), false);

     auto do_visit_rings = [&](hpuint p, hpuint i) {
          std::vector<T> ring;
          visit_ring<degree>(patches, neighbors, p, i, [&](auto t) { ring.push_back(t); });
          visit(p, i, ring);
          visit_fan(neighbors, p, i, [&](hpuint q, hpuint j) {
               if(j == 0) visited[3 * q] = true;
               else if(j == 1) visited[3 * q + 1] = true;
               else visited[3 * q + 2] = true;
          });
     };

     for(auto p : boost::irange(0u, nPatches)) {
          if(!visited[3 * p]) do_visit_rings(p, 0);
          if(!visited[3 * p + 1]) do_visit_rings(p, 1);
          if(!visited[3 * p + 2]) do_visit_rings(p, 2);
     }
}

template<class Space, hpuint degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, Visitor&& visit) {
     auto patches = surface.getPatches();
     auto temp = deindex(std::get<0>(patches), std::get<1>(patches));
     visit_rings<degree>(temp.begin(), temp.end(), neighbors, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_rings(surface, neighbors, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_subring(Iterator patches, const Indices& neighbors, hpuint p, hpuint i, hpuint q, Visitor&& visit) {
     static constexpr auto patchSize = make_patch_size(degree);
     hpuint k;
     visit_subfan(neighbors, p, i, q, [&](hpuint r, hpuint j) {
          k = j;
          visit_patch<degree>(patches, r, [&](auto patch) {
               if(j == 0) visit(*(patch + 1));
               else if(j == 1) visit(*(patch + + (degree << 1)));
               else visit(*(patch + (patchSize - 3)));
          });
     });
     visit_patch<degree>(patches, q, [&](auto patch) {
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
     if(i == UNULL) visit_corners<degree>(indices.begin(), p, [&](hpuint v0, hpuint v1, hpuint v2) {
          visit_corners<degree>(indices.begin(), q, [&](hpuint w0, hpuint w1, hpuint w2) { i = (v0 == w0 || v0 == w1 || v0 == w2) ? 0 : (v1 == w0 || v1 == w1 || v1 == w2) ? 1 : 2; });
     });
     visit_subring<degree>(patches.begin(), neighbors, p, i, q, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint q, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_subring(surface, neighbors, p, q, std::forward<Visitor>(visit));
}

//algorithms

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, (degree + 1)> elevate(const SurfaceSplineBEZ<Space, degree>& surface) {
     SurfaceSplineBEZ<Space, (degree + 1)> surface1(surface.getNumberOfPatches());
     auto& indices = std::get<1>(surface.getPatches());
     auto neighbors = make_neighbors(surface);
     auto patches = deindex(surface.getPatches());

     visit_fans(neighbors, [&](auto p, auto i, auto begin, auto end) {
          visit_corner<degree>(patches.begin(), p, i, [&](auto& corner) { surface1.setCorner(p, i, corner); });
          visit_pairs(begin, std::distance(begin, end) / 2, 2, [&](auto q, auto j) { surface1.setCorner(q, j, p, i); });
     });

     auto elevate_boundary = [&](auto p, auto i) {
          visit_patch<degree>(patches.begin(), p, [&](auto patch) {
               auto boundary = make_boundary<degree>(patch, i);
               visit_ends<degree>(patch, i, [&](auto& corner0, auto& corner1) { surface1.setBoundary(p, i, happah::curves::elevate(degree, corner0, boundary.begin(), corner1).begin()); });
          });
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
     visit_patches<degree>(patches.begin(), surface.getNumberOfPatches(), [&](auto patch) {
          using Point = typename Space::POINT;
          std::vector<Point> interior;
          interior.reserve(make_patch_size(degree + 1u) - 3u * (degree + 1u));
          auto alpha = hpreal(1.0 / (degree + 1u));
          auto t0 = degree;
          auto j0 = 0u, j1 = 0u, j2 = 0u;
          visit_nablas<degree>(patch, [&](auto& b0, auto& b1, auto& b2) {
               if(j0 == 0u) {
                    j0 = --t0;
                    j2 = degree - t0 + 1u;
                    j1 = degree - j0 - j2 + 1u;
               }
               interior.push_back(alpha * (hpreal(j0) * b0 + hpreal(j1) * b1 + hpreal(j2) * b2));
               --j0;
               ++j1;
          });
          surface1.setInterior(p, interior.begin());
          ++p;
     });
     return surface1;
}

template<class Space, hpuint degree>
bool is_c0(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i) {
     auto& indices = std::get<1>(surface.getPatches());
     auto q = neighbors[3 * p + i];
     auto j = make_neighbor_offset(neighbors, q, p);

     auto k0 = 0u, k1 = 0u, l0 = 0u, l1 = 0u;
     visit_ends<degree>(indices.begin(), p, i, [&](auto k, auto l) { k0 = k; l0 = l; });
     visit_ends<degree>(indices.begin(), q, j, [&](auto k, auto l) { k1 = k; l1 = l; });
     if(k0 != l1 || l0 != k1) return false;
     auto boundary0 = make_boundary<degree>(indices.begin(), p, i); 
     auto boundary1 = make_boundary<degree>(indices.begin(), q, j);
     std::reverse(boundary1.begin(), boundary1.end());
     return boundary0 == boundary1;
}

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_boundary(Iterator patch, hpuint i) {
     std::vector<typename std::iterator_traits<Iterator>::value_type> boundary;
     boundary.reserve(degree - 1);
     visit_boundary<degree>(patch, i, [&](auto& value) { boundary.push_back(value); });
     return boundary;
}

template<hpuint degree, class Iterator>
std::vector<typename std::iterator_traits<Iterator>::value_type> make_boundary(Iterator patches, hpuint p, hpuint i) {
     std::vector<typename std::iterator_traits<Iterator>::value_type> boundary;
     visit_patch<degree>(patches, p, [&](auto patch) { boundary = make_boundary<degree>(patch, i); });
     return boundary;
}

template<class Space, hpuint degree, class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_control_polygon(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory) { return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, std::forward<VertexFactory>(factory)); }

template<class Space, hpuint degree>
std::vector<hpuint> make_neighbors(const SurfaceSplineBEZ<Space, degree>& surface) {
     Indices indices;
     visit_patches<degree>(std::get<1>(surface.getPatches()).begin(), surface.getNumberOfPatches(), [&](auto patch) {
          visit_corners<degree>(patch, [&](hpuint i0, hpuint i1, hpuint i2) {
               indices.push_back(i0);
               indices.push_back(i1);
               indices.push_back(i2);
          });
     });
     return make_neighbors(indices);
}

constexpr hpuint make_patch_size(hpuint degree) { return (degree + 1) * (degree + 2) >> 1; }

template<class Space, hpuint degree, class Vertex, class VertexFactory, typename>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, VertexFactory&& factory) {
     using Point = typename Space::POINT;
     static_assert(std::is_base_of<Vertex, decltype(factory(Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");

     std::vector<Vertex> vertices;
     vertices.reserve(surface.getControlPoints().size());
     for(auto& point : surface.getControlPoints()) vertices.push_back(factory(point));

     auto indices = SurfaceSplineUtilsBEZ::template buildTriangleMeshIndices<degree>(std::get<1>(surface.getPatches()));

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

template<class Space, hpuint degree, class Vertex, class VertexFactory, typename>
TriangleMesh<Vertex> make_triangle_mesh(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions, VertexFactory&& factory) {
     if(nSubdivisions > 0) return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(subdivide(surface, nSubdivisions), std::forward<VertexFactory>(factory));
     else return make_triangle_mesh<Space, degree, Vertex, VertexFactory>(surface, std::forward<VertexFactory>(factory));
}

//TODO: move non-member functions with iterators into subnamespace so as not to conflict with implementations for curves, for example
template<hpuint degree, class Iterator, class Visitor>
void sample(Iterator patches, hpuint nPatches, hpuint nSamples, Visitor&& visit) {
     //TODO; skip multiple computations on common edges and eliminate common points in array
     auto matrix = SurfaceUtilsBEZ::getEvaluationMatrix<degree>(nSamples);
     visit_patches<degree>(patches, nPatches, [&](auto patch) {
          static constexpr auto patchSize = make_patch_size(degree);
          auto end = patch + patchSize;
          for(auto m = matrix.begin(), mend = matrix.end(); m != mend; ++m) {
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
     auto matrixd = SurfaceUtilsBEZ::getEvaluationMatrix<degree>(nSamples);
     auto matrix1 = SurfaceUtilsBEZ::getEvaluationMatrix<1>(nSamples);
     visit_patches<degree>(controlPoints, nPatches, [&](auto patch) {
          static constexpr auto patchSize = make_patch_size(degree);
          auto end = patch + patchSize;
          auto& p0 = *domainPoints;
          auto& p1 = *(++domainPoints);
          auto& p2 = *(++domainPoints);
          ++domainPoints;

          auto d = matrix1.begin();
          for(auto m = matrixd.begin(), mend = matrixd.end(); m != mend; ++m) {
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
     sample<degree>(patches.begin(), surface.getNumberOfPatches(), nSamples, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class T, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, std::tuple<const std::vector<T>&, const Indices&> domain, hpuint nSamples, Visitor&& visit) {
     auto controlPoints = deindex(surface.getPatches());
     auto domainPoints = deindex(domain);
     sample<degree>(controlPoints.begin(), domainPoints.begin(), surface.getNumberOfPatches(), nSamples, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions) {
     using Point = typename Space::POINT;
     static constexpr hpuint nTriangles = SurfaceUtilsBEZ::get_number_of_control_polygon_triangles<degree>::value;

     if(nSubdivisions == 0) return surface;

     auto nPatches = surface.getNumberOfPatches();

     std::vector<SurfaceSubdividerBEZ<Space, degree> > subdividers;
     subdividers.reserve(nPatches);
     auto push_patch = [&](auto patch) { subdividers.emplace_back(patch); };
     visit_patches(surface, push_patch);

     std::vector<Point> points;
     Indices indices;

     points.reserve(make_patch_size(degree) * nPatches);
     indices.reserve(3 * nTriangles * nPatches);

     for(auto& subdivider : subdividers) {
          auto subdivided = subdivider.subdivide(nSubdivisions);
          auto offset = points.size();
          for(auto& i : std::get<0>(subdivided)) i += offset;
          std::move(begin(std::get<0>(subdivided)), end(std::get<0>(subdivided)), std::back_inserter(points));
          std::move(begin(std::get<1>(subdivided)), end(std::get<1>(subdivided)), std::back_inserter(indices));
     }

     return { std::move(points), std::move(indices) };
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

     visit_rings<degree>(indices.begin(), surface.getNumberOfPatches(), neighbors, [&](hpuint p, hpuint i, auto ring) {
          hpuint v;
          visit_corner<degree>(indices.begin(), p, [&](auto& w) { v = w; });

          auto& c0 = points[v];
          auto& c1 = points[ring[0]];
          auto& c2 = points[ring[1]];

          auto A = glm::inverse(hpmat3x3(c0, c1, c2));

          for(auto& q : ring) coordinates[q] = A * points[q];
     });

     std::vector<hpir> hirs;
     std::vector<hpijr> hijrs;

     auto row = 0u;
     visit_edges(neighbors, [&](hpuint p, hpuint i) {
          hpuint b[3];
          hpuint c[3];
          auto q = neighbors[3 * p + i];
          auto j = 0u;

          auto tb = b + 3;
          auto tc = c - 1;
          visit_triplet(neighbors, q, [&](hpuint n0, hpuint n1, hpuint n2) {
               j = (p == n0) ? 0 : (p == n1) ? 1 : 2;
               visit_subring<degree>(indices.begin(), neighbors, q, (p == n0) ? 1 : (p == n1) ? 2 : 0, p, [&](auto k) { *(--tb) = k; });
          });
          visit_subring<degree>(indices.begin(), neighbors, p, (i == 0) ? 1 : (i == 1) ? 2 : 0, q, [&](auto k) { *(++tc) = k; });

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

}//namespace happah

