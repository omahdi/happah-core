// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <type_traits>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/Surface.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/geometries/TriangleMeshUtils.h"
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

     SurfaceSplineBEZ(ControlPoints controlPoints, Indices indices)
          : m_controlPoints{std::move(controlPoints)}, m_indices{std::move(indices)} {}

     SurfaceSplineBEZ(ControlPoints controlPoints)
          : m_controlPoints{std::move(controlPoints)}, m_indices{m_controlPoints.size()} { std::iota(std::begin(m_indices), std::end(m_indices), 0); }

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     hpuint getNumberOfPatches() const { return m_indices.size() / SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value; }

     std::tuple<const ControlPoints&, const Indices&> getPatches() const { return std::tie(m_controlPoints, m_indices); }

private:
     ControlPoints m_controlPoints;
     Indices m_indices;

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const SurfaceSplineBEZ<Space, t_degree>& surface) {
          stream << surface.m_controlPoints << '\n';
          stream << surface.m_indices;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, SurfaceSplineBEZ<Space, t_degree>& surface) {
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
void visit_corners(Iterator begin, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     visit(*begin, *(begin + degree), *(begin + (nControlPoints - 1)));
}

template<class Space, hpuint degree, class Visitor>
void visit_edges(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_edges(neighbors, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_fans(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_fans(neighbors, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_patch(Iterator begin, hpuint p, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     visit(begin + p * nControlPoints, begin + ((p + 1) * nControlPoints));
}

template<class Space, hpuint degree, class Visitor>
void visit_patch(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, Visitor&& visit) {
     auto patches = surface.getPatches();
     auto temp = deindex(std::get<0>(patches), std::get<1>(patches));
     visit_patch<degree>(temp.begin(), p, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_patches(Iterator begin, hpuint nPatches, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     for(auto i = begin, end = begin + nPatches * nControlPoints; i != end; i += nControlPoints) visit(i, i + nControlPoints);
}

template<hpuint degree, class Iterator, class Visitor>
void visit_patches(Iterator begin, Iterator end, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     for(auto i = begin; i != end; i += nControlPoints) visit(i, i + nControlPoints);
}

template<class Space, hpuint degree, class Visitor>
void visit_patches(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto patches = surface.getPatches();
     auto temp = deindex(std::get<0>(patches), std::get<1>(patches));
     visit_patches<degree>(temp.begin(), temp.end(), std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_ring(Iterator begin, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     hpuint first = UNULL, last, n = 0u, k;

     visit_fan(neighbors, p, i, [&](hpuint q, hpuint j) {
          if(first == UNULL) first = q;
          last = q;
          k = j;
          ++n;
          if(j == 0) visit(*(begin + (q * nControlPoints + 1)));
          else if (j == 1) visit(*(begin + (q * nControlPoints + (degree << 1))));
          else visit(*(begin + ((q + 1) * nControlPoints - 3)));
     });
     if(n < 3 || !is_neighbor(neighbors, first, last)) {
          if(k == 0) visit(*(begin + (last * nControlPoints + degree + 1)));
          else if (k == 1) visit(*(begin + (last * nControlPoints + degree - 1)));
          else visit(*(begin + ((last + 1) * nControlPoints - 2)));
     }
}

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit) {
     auto patches = surface.getPatches();
     auto temp = deindex(std::get<0>(patches), std::get<1>(patches));
     visit_ring<degree>(temp.begin(), neighbors, p, i, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint i, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_ring(surface, neighbors, p, i, std::forward<Visitor>(visit));
}

template<hpuint degree, class Iterator, class Visitor>
void visit_rings(Iterator begin, Iterator end, const Indices& neighbors, Visitor&& visit) {
     using T = typename Iterator::value_type;
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     boost::dynamic_bitset<> visited(neighbors.size(), false);

     auto do_visit_rings = [&](hpuint p, hpuint i) {
          std::vector<T> ring;
          visit_ring<degree>(begin, neighbors, p, i, [&](auto t) { ring.push_back(t); });
          visit(p, i, ring);
          visit_fan(neighbors, p, i, [&](hpuint q, hpuint j) {
               if(j == 0) visited[3 * q] = true;
               else if(j == 1) visited[3 * q + 1] = true;
               else visited[3 * q + 2] = true;
          });
     };

     for(auto p : boost::irange(0l, std::distance(begin, end) / nControlPoints)) {
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

/**
 * Visit the subring starting at patch p rotating counterclockwise and stopping at patch q.
 */
template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint q, Visitor&& visit) {
	static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     auto patches = surface.getPatches();
     auto& indices = std::get<1>(patches);
     auto begin = deindex(std::get<0>(patches), indices).begin();
     hpuint i;
     visit_triplet(neighbors, p, [&](hpuint n0, hpuint n1, hpuint n2) { i = (q == n0) ? 1 : (q == n1) ? 2 : (q == n2) ? 0 : UNULL; });
     if(i == UNULL) visit_corners<degree>(indices.begin() + p * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) {
          visit_corners<degree>(indices.begin() + q * nControlPoints, [&](hpuint w0, hpuint w1, hpuint w2) { i = (v0 == w0 || v0 == w1 || v0 == w2) ? 0 : (v1 == w0 || v1 == w1 || v1 == w2) ? 1 : 2; });
     });

     hpuint k;
     visit_subfan(neighbors, p, i, q, [&](hpuint r, hpuint j) {
          k = j;
          if(j == 0) visit(*(begin + (r * nControlPoints + 1)));
          else if(j == 1) visit(*(begin + (r * nControlPoints + (degree << 1))));
          else visit(*(begin + ((r + 1) * nControlPoints - 3)));
     });
     if(k == 0) visit(*(begin + (q * nControlPoints + degree + 1)));
     else if (k == 1) visit(*(begin + (q * nControlPoints + degree - 1)));
     else visit(*(begin + ((q + 1) * nControlPoints - 2)));
}

template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint q, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_subring(surface, neighbors, p, q, std::forward<Visitor>(visit));
}

//algorithms

namespace tpfssb {

template<class Space, hpuint degree>
std::vector<hpijr> make_objective(const SurfaceSplineBEZ<Space, degree>& surface) {
     //TODO SM
     using Point = typename Space::POINT;
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     auto patches = surface.getPatches();
     auto& indices = std::get<1>(patches);
     auto& points = std::get<0>(patches);
     auto neighbors = make_neighbors(surface);
     auto begin = deindex(std::get<0>(patches), indices).begin();

     std::unordered_map<hpuint, Point> coordCPs;

     //Calculate the transition functions within every ring starting from a projective frame
     visit_rings(indices.begin(), indices.end(), neighbors, [&](hpuint p, hpuint i, auto ring) {
               //Assuming ring points go CCW
               //Get centre
               hpuint v;
               visit_corners<degree>(indices.begin() + p * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) { v = (i == 0) ? v0 : (i == 1) ? v1 : v2; });

               auto homogenise = [&](const Point&& p) {
                    auto z = p.z;
                    return Point(p.x/z, p.y/z, 1.0f);
               };

               Point c0 = homogenise(*(points.begin + v)); //Is this the right level of indirection???
               assert(ring.size() >= 2);
               //Make frame and build the projection matrix (we are going to use the inverse)
               Point c1 = homogenise(points[ring[0]]);
               Point c2 = homogenise(points[ring[1]]);

               hpmat3x3 mat(c0, c1, c2);
               hpmat3x3 invMat = glm::inverse(mat);

               //Now calculate the coordinates of the control points on the ring
               for(auto& p : ring) {
                    Point coords = invMat * points[p];
                    coordCPs[p] = homogenise(coords);
               }
          }); //visit_rings

     //Visit all the edges (i.e. pair of patches) to gather the transition functions
     std::vector<hpijr> triplets;

     visit_edges(neighbors, [&](hpuint p, hpuint i) {
               //`q` is CW!! in this case
               auto q = *(neighbors.begin() + 3*p + i);
               visit_subring(surface, q, p, [&](auto i0, auto i1, auto i2){ //Note the swapping of `q` and `p`!!!!
                         //Encode this 3x3 matrix as triplets (row, col, val)
                         //NOTE: This takes up a 9x9 space in the sparse matrix!!!
                         for(auto i = 0; i < 3; i++) {
                              triplets.emplace_back(9*q + 0 + i, 27*p + (3*i) + 0, coordCPs[i0].x);
                              triplets.emplace_back(9*q + 0 + i, 27*p + (3*i) + 1, coordCPs[i0].y);
                              triplets.emplace_back(9*q + 0 + i, 27*p + (3*i) + 2, coordCPs[i0].z);

                              triplets.emplace_back(9*q + 3 + i, 27*p + (3*i) + 0, coordCPs[i1].x);
                              triplets.emplace_back(9*q + 3 + i, 27*p + (3*i) + 1, coordCPs[i1].y);
                              triplets.emplace_back(9*q + 3 + i, 27*p + (3*i) + 2, coordCPs[i1].z);

                              triplets.emplace_back(9*q + 6 + i, 27*p + (3*i) + 0, coordCPs[i2].x);
                              triplets.emplace_back(9*q + 6 + i, 27*p + (3*i) + 1, coordCPs[i2].y);
                              triplets.emplace_back(9*q + 6 + i, 27*p + (3*i) + 2, coordCPs[i2].z);
                         }
                    });

               //The other side of the egde
               visit_subring(surface, p, q, [&](auto i0, auto i1, auto i2){
                         //Encode this 3x3 matrix as triplets (row, col, val)
                         for(auto i = 0; i < 3; i++) {
                              triplets.emplace_back(9*p + 0 + i, 27*q + (3*i) + 0, coordCPs[i0].x);
                              triplets.emplace_back(9*p + 0 + i, 27*q + (3*i) + 1, coordCPs[i0].y);
                              triplets.emplace_back(9*p + 0 + i, 27*q + (3*i) + 2, coordCPs[i0].z);

                              triplets.emplace_back(9*p + 3 + i, 27*q + (3*i) + 0, coordCPs[i1].x);
                              triplets.emplace_back(9*p + 3 + i, 27*q + (3*i) + 1, coordCPs[i1].y);
                              triplets.emplace_back(9*p + 3 + i, 27*q + (3*i) + 2, coordCPs[i1].z);

                              triplets.emplace_back(9*p + 6 + i, 27*q + (3*i) + 0, coordCPs[i2].x);
                              triplets.emplace_back(9*p + 6 + i, 27*q + (3*i) + 1, coordCPs[i2].y);
                              triplets.emplace_back(9*p + 6 + i, 27*q + (3*i) + 2, coordCPs[i2].z);
                         }
                    });
          });

     return triplets;
}

/**
 *  The constraints are quadratics of the form $$ax_jx_k+bx_j=c$$, where the quadratic coefficients are stored in the first vector, the linear coefficients in the second, and the constant coefficients in the third.  The first integer (the 'i') in each entry of the three vectors identifies the constraint.
 */

template<class Space, hpuint degree>
std::tuple<std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const SurfaceSplineBEZ<Space, degree>& surface) {
     std::vector<hpijkr> ijkrs;
     std::vector<hpijr> ijrs;
     std::vector<hpir> irs;

     return std::make_tuple(std::move(ijkrs), std::move(ijrs), std::move(irs));
}

}//namespace tpfssb

template<hpuint n, class Space, hpuint degree>
SurfaceSplineBEZ<Space, (degree + n)> elevate(const SurfaceSplineBEZ<Space, degree>& surface) {
     //TODO: TR specialize for n = 1
     if(n == 1) {
     } else {

     }
     return {};
}

template<class Space, hpuint degree>
std::vector<hpuint> make_neighbors(const SurfaceSplineBEZ<Space, degree>& surface) {
     Indices indices;
     visit_patches<degree>(std::get<1>(surface.getPatches()).begin(), surface.getNumberOfPatches(), [&](auto begin, auto end) {
          visit_corners<degree>(begin, [&](hpuint i0, hpuint i1, hpuint i2) {
               indices.push_back(i0);
               indices.push_back(i1);
               indices.push_back(i2);
          });
     });
     return make_neighbors(indices);
}

template<class Space, hpuint degree>
SurfaceSplineBEZ<Space, degree> subdivide(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSubdivisions) {
     using Point = typename Space::POINT;
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     static constexpr hpuint nTriangles = SurfaceUtilsBEZ::get_number_of_control_polygon_triangles<degree>::value;

     if(nSubdivisions == 0) return surface;

     auto nPatches = surface.getNumberOfPatches();

     std::vector<SurfaceSubdividerBEZ<Space, degree> > subdividers;
     subdividers.reserve(nPatches);
     auto push_patch = [&](auto begin, auto end) { subdividers.emplace_back(begin); };
     visit_patches(surface, push_patch);

     std::vector<Point> points;
     Indices indices;

     points.reserve(nControlPoints * nPatches);
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

//TODO: move non-member functions with iterators into subnamespace so as not to conflict with implementations for curves, for example
template<hpuint degree, class Iterator, class Visitor>
void sample(Iterator begin, hpuint nPatches, hpuint nSamples, Visitor&& visit) {
     //TODO; skip multiple computations on common edges and eliminate common points in array
     auto matrix = SurfaceUtilsBEZ::getEvaluationMatrix<degree>(nSamples);
     visit_patches<degree>(begin, nPatches, [&](auto begin, auto end) {
          for(auto m = matrix.begin(), mend = matrix.end(); m != mend; ++m) {
               auto temp = begin;
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
     visit_patches<degree>(controlPoints, nPatches, [&](auto begin, auto end) {
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

               auto temp = begin;
               auto sample = *m * *temp;
               while(++temp != end) sample += *(++m) * *temp;
               visit(mix(p0, u, p1, v, p2, w), sample);
          }
     });
}

template<class Space, hpuint degree, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, hpuint nSamples, Visitor&& visit) {
     auto patches = surface.getPatches();
     sample<degree>(deindex(std::get<0>(patches), std::get<1>(patches)).begin(), surface.getNumberOfPatches(), nSamples, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class T, class Visitor>
void sample(const SurfaceSplineBEZ<Space, degree>& surface, std::tuple<const std::vector<T>&, const Indices&> domain, hpuint nSamples, Visitor&& visit) {
     auto patches = surface.getPatches();
     auto controlPoints = deindex(std::get<0>(patches), std::get<1>(patches)).begin();
     auto domainPoints = deindex(std::get<0>(domain), std::get<1>(domain)).begin();
     sample<degree>(controlPoints, domainPoints, surface.getNumberOfPatches(), nSamples, std::forward<Visitor>(visit));
}

}//namespace happah
