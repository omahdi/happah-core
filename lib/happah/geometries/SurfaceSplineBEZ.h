// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <type_traits>
#include <vector>

#include "happah/Eigen.h"
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

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     auto patches = surface.getPatches();
     auto& indices = std::get<1>(patches);
     auto begin = deindex(std::get<0>(patches), indices).begin();
     hpuint v;
     visit_corners<degree>(indices.begin() + p * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) { v = (i == 0) ? v0 : (i == 1) ? v1 : v2; });
     hpuint first = UNULL, last, n = 0u;

     visit_fan(neighbors, p, i, [&](hpuint f) {
          if(first == UNULL) first = f;
          last = f;
          ++n;
          visit_corners<degree>(indices.begin() + f * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) {
               if(v == v0) visit(*(begin + (f * nControlPoints + 1)));
               else if (v == v1) visit(*(begin + (f * nControlPoints + (degree << 1))));
               else visit(*(begin + ((f + 1) * nControlPoints - 3)));
          });
     });
     if(n < 3 || !is_neighbor(neighbors, first, last)) {
          visit_corners<degree>(indices.begin() + last * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) {
               if(v == v0) visit(*(begin + (last * nControlPoints + degree + 1)));
               else if (v == v1) visit(*(begin + (last * nControlPoints + degree - 1)));
               else visit(*(begin + ((last + 1) * nControlPoints - 2)));
          });
     }
}

template<class Space, hpuint degree, class Visitor>
void visit_ring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint i, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_ring(surface, neighbors, p, i, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, Visitor&& visit) {
     using Point = typename Space::POINT;
     boost::dynamic_bitset<> visited(neighbors.size(), false);
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     auto& indices = std::get<1>(surface.getPatches());

     auto do_visit_rings = [&](hpuint p, hpuint i) {
          std::vector<Point> ring;
          visit_ring(surface, neighbors, p, i, [&](auto point) { ring.push_back(point); });
          visit(p, i, ring);
          hpuint v;
          visit_corners<degree>(indices.begin() + p * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) { v = (i == 0) ? v0 : (i == 1) ? v1 : v2; });
          visit_fan(neighbors, p, i, [&](hpuint f) {
               visit_corners<degree>(indices.begin() + f * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) {
                    if(v == v0) visited[3 * f] = true;
                    else if(v == v1) visited[3 * f + 1] = true;
                    else visited[3 * f + 2] = true;
               });
          });
     };

     for(auto p = 0u, end = surface.getNumberOfPatches(); p != end; ++p) {
          if(!visited[3 * p]) do_visit_rings(p, 0);
          if(!visited[3 * p + 1]) do_visit_rings(p, 1);
          if(!visited[3 * p + 2]) do_visit_rings(p, 2);
     }
}

template<class Space, hpuint degree, class Visitor>
void visit_rings(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_rings(surface, neighbors, std::forward<Visitor>(visit));
}


//TODO This indexed functions are only to make the tracking of the proj. coordinates easier. There has to be a better way.
template<class Space, hpuint degree, class Visitor>
void visit_ring_indexed(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint i, Visitor&& visit) {
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     auto patches = surface.getPatches();
     auto& indices = std::get<1>(patches);
     auto begin = deindex(std::get<0>(patches), indices).begin();
     hpuint v;
     visit_corners<degree>(indices.begin() + p * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) { v = (i == 0) ? v0 : (i == 1) ? v1 : v2; });
     hpuint first = UNULL, last, n = 0u;

     visit_fan(neighbors, p, i, [&](hpuint f) {
          if(first == UNULL) first = f;
          last = f;
          ++n;
          visit_corners<degree>(indices.begin() + f * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) {
		          if(v == v0) {
			          auto idx = f * nControlPoints + 1;
			          visit(idx, *(begin + idx));
		          } else if (v == v1) {
			          auto idx = f * nControlPoints + (degree << 1);
			          visit(idx, *(begin + idx));
		          } else {
			          auto idx = (f + 1) * nControlPoints - 3;
			          visit(idx, *(begin + idx));
		          }
          });
     });
     hpuint idx;
     if(n < 3 || !is_neighbor(neighbors, first, last)) {
          visit_corners<degree>(indices.begin() + last * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) {
               if(v == v0) {
	               idx = last * nControlPoints + degree + 1;
	               visit(idx, *(begin + idx));
               } else if (v == v1) {
	               idx = last * nControlPoints + degree - 1;
	               visit(idx, *(begin + idx));
               } else {
	               idx = (last + 1) * nControlPoints - 2;
	               visit(idx, *(begin + idx));
               }
          });
     }
}

template<class Space, hpuint degree, class Visitor>
void visit_rings_indexed(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, Visitor&& visit) {
     using Point = typename Space::POINT;
     boost::dynamic_bitset<> visited(neighbors.size(), false);
     static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     auto& indices = std::get<1>(surface.getPatches());

     auto do_visit_rings = [&](hpuint p, hpuint i) {
          std::vector<Point> ring;
          visit_ring_indexed(surface, neighbors, p, i, [&](auto idx, auto point) { ring.push_back(std::make_pair(idx, point)); });
          visit(p, i, ring);
          hpuint v;
          visit_corners<degree>(indices.begin() + p * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) { v = (i == 0) ? v0 : (i == 1) ? v1 : v2; });
          visit_fan(neighbors, p, i, [&](hpuint f) {
               visit_corners<degree>(indices.begin() + f * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) {
                    if(v == v0) visited[3 * f] = true;
                    else if(v == v1) visited[3 * f + 1] = true;
                    else visited[3 * f + 2] = true;
               });
          });
     };

     for(auto p = 0u, end = surface.getNumberOfPatches(); p != end; ++p) {
          if(!visited[3 * p]) do_visit_rings(p, 0);
          if(!visited[3 * p + 1]) do_visit_rings(p, 1);
          if(!visited[3 * p + 2]) do_visit_rings(p, 2);
     }
}

template<class Space, hpuint degree, class Visitor>
void visit_rings_indexed(const SurfaceSplineBEZ<Space, degree>& surface, Visitor&& visit) {
     auto neighbors = make_neighbors(surface);
     visit_rings_indexed(surface, neighbors, std::forward<Visitor>(visit));
}

template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, const Indices& neighbors, hpuint p, hpuint q, Visitor&& visit) {
     //TODO: SM
	//Starting at `p` and going CCW to `q`. The required corner should be clear this way.
	//`visit` expects the three control points of the half-diamond.

	static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     auto patches = surface.getPatches();
     auto& indices = std::get<1>(patches);
     auto begin = deindex(std::get<0>(patches), indices).begin();

     //Get the "direction" `i` to `q` to calculate which corner of `p` (CCW is the one to the left of the incoming idx (probably?))
     hpuint i;
     visit_triplet(neighbors, p, [&](hpuint n0, hpuint n1, hpuint n2) { i = (q == n0) ? 1 : (q == n1) ? 2 : 0; });
     hpuint i0, i1, i2;
     if(i == 0) { //Left corner
	     i0 = p * nControlPoints + 0;
	     i1 = p * nControlPoints + 1;
	     i2 = p * nControlPoints + (degree + 1);
     } else if(i == 1) { //Right corner
	     i0 = p * nControlPoints + (degree - 1);
	     i1 = p * nControlPoints + (degree + 0);
	     i2 = p * nControlPoints + (degree << 1);
     } else { //Top corner
	     i0 = (p + 1) * nControlPoints - 1;
	     i1 = (p + 1) * nControlPoints - 3;
	     i2 = (p + 1) * nControlPoints - 2;
     }
     hpuint v0, v1, v2;
     v0 = *(begin + i0);
     v1 = *(begin + i1);
     v2 = *(begin + i2);

     visit(v0, v1, v2);
}

template<class Space, hpuint degree, class Visitor>
void visit_subring(const SurfaceSplineBEZ<Space, degree>& surface, hpuint p, hpuint q, Visitor&& visit) {
     //TODO: SM
	auto neighbors = make_neighbors(surface);
     visit_subring(surface, neighbors, p, q, std::forward<Visitor>(visit));
}

//algorithms

template<hpuint n, class Space, hpuint degree>
SurfaceSplineBEZ<Space, (degree + n)> elevate(const SurfaceSplineBEZ<Space, degree>& surface) {
     //TODO: TR specialize for n = 1
     if(n == 1) {
     } else {

     }
     return {};
}

namespace tpfssb {

template<class Space, hpuint degree>
std::vector<Eigen::Triplet<hpreal> > make_objective_function(const SurfaceSplineBEZ<Space, degree>& surface) {
	//TODO SM
	using Point = typename Space::POINT;
	static constexpr hpuint nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<degree>::value;
     auto patches = surface.getPatches();
     auto& indices = std::get<1>(patches);
	auto& points = std::get<0>(patches);
     auto begin = deindex(std::get<0>(patches), indices).begin();

	std::unordered_map<hpuint, Point> coordCPs;

	//Calculate the transition functions within every ring starting from a projective frame
	visit_rings_indexed(surface, [&](hpuint p, hpuint i, std::pair<hpuint, Point> ringIdx) {
			//Assuming ring points go CCW
			//Get centre
			hpuint v;
			visit_corners<degree>(indices.begin() + p * nControlPoints, [&](hpuint v0, hpuint v1, hpuint v2) { v = (i == 0) ? v0 : (i == 1) ? v1 : v2; });

			auto homogenise = [&](const Point&& p) {
				auto z = p.z;
				return Point(p.x/z, p.y/z, 1.0f);
			};

			Point c0 = homogenise(*(points.begin + v)); //Is this the right level of indirection???
			assert(ringIdx.size() >= 2);
			//Make frame and build the projection matrix (we are going to use the inverse)
			Point c1 = homogenise(ringIdx[0].second);
			Point c2 = homogenise(ringIdx[1].second);

			hpmat3x3 mat(c0, c1, c2);
			hpmat3x3 invMat = glm::inverse(mat);

			//Now calculate the coordinates of the control points on the ring
			for(auto& p : ringIdx) {
				Point coords = invMat * p.second;
				coordCPs[p.first] = homogenise(coords);
			}
		}); //visit_rings
}//make_objective_function
} //tfssp

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
