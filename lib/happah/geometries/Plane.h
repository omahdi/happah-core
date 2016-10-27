// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/optional.hpp>
#include <memory>
#include <list>
#include <random>
#include <type_traits>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/Geometry.h"
#include "happah/geometries/QuadMesh.h"
#include "happah/geometries/Ray.h"
#include "happah/geometries/SegmentMesh.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/geometries/VertexCloud.h"
#include "happah/math/RigidAffineTransformation.h"
#include "happah/utils/VertexFactory.h"

namespace happah {

//TODO: render planar structures as texture on top of plane
//TODO: move struct Utils into [class]Utils to be consistent
class Plane : public Geometry2D<Space3D> {
public:
     struct Utils {
          class PointFactory {
          public:
               typedef Point3D PRODUCT;

               PointFactory(const Plane& plane);

               Point3D build(const Point2D& abscissa, const Point1D& ordinate) const;

          private:
               const Plane& m_plane;
               Vector3D m_xAxis;
               Vector3D m_yAxis;

          };

          template<class Vertex>
          class BestFittingPlaneNormalSetter {
               static_assert(is_relative_vertex<Vertex>::value && std::is_same<typename Vertex::SPACEA, Space2D>::value && std::is_same<typename Vertex::SPACEO, Space1D>::value && contains_normal<Vertex>::value, "A neighborhood visitor was designed to be used with relative vertices that have 2D abscissae, 1D ordinates, and normals.");
               typedef Plane::Utils::PointFactory PointFactory;

          public:
               BestFittingPlaneNormalSetter(const Plane& plane) 
                    : m_pointFactory(plane) {}
               void visit(Vertex& v, const Vertex& n1) const { v.normal = Vector3D(0.0, 0.0, 1.0); }
               void visit(Vertex& v, const Vertex& n1, const Vertex& n2) const { v.normal = Vector3D(0.0, 0.0, 1.0); }
               void visit(Vertex& v, const Vertex& n1, const Vertex& n2, const Vertex& n3) const { v.normal = Plane::Utils::getBestFittingPlaneNormal({getPoint(v), getPoint(n1), getPoint(n2), getPoint(n3)}); }
               void visit(Vertex& v, const Vertex& n1, const Vertex& n2, const Vertex& n3, const Vertex& n4, const Vertex& n5) const { v.normal = Plane::Utils::getBestFittingPlaneNormal({getPoint(v), getPoint(n1), getPoint(n2), getPoint(n3), getPoint(n4), getPoint(n5)}); }
               void visit(Vertex& v, const Vertex& n1, const Vertex& n2, const Vertex& n3, const Vertex& n4, const Vertex& n5, const Vertex& n6, const Vertex& n7, const Vertex& n8) const { v.normal = Plane::Utils::getBestFittingPlaneNormal({getPoint(v), getPoint(n1), getPoint(n2), getPoint(n3), getPoint(n4), getPoint(n5), getPoint(n6), getPoint(n7), getPoint(n8)}); }

          private:
               PointFactory m_pointFactory;

               Point3D getPoint(const Vertex& v) const { return m_pointFactory.build(v.abscissa, v.ordinate); }

          };

          static Vector3D getBestFittingPlaneNormal(const std::vector<Point3D>& points, hpreal epsilon = happah::EPSILON);

     };

     static boost::optional<hpreal> intersect(const Point3D& origin, const Vector3D& normal, const Ray3D& ray, hpreal epsilon = happah::EPSILON);

     Plane(Point3D origin, const Vector3D& normal, hpreal epsilon = happah::EPSILON);

     std::pair<Vector3D, Vector3D> getAxes() const;

     const Vector3D& getNormal() const;

     const Point3D& getOrigin() const;

     RigidAffineTransformation3D getPlaneTransformation(hpreal epsilon = happah::EPSILON) const;

     boost::optional<hpreal> intersect(const Ray3D& ray, hpreal epsilon = happah::EPSILON) const;

     //NOTE: Point must be in local coordinate system of the plane.
     Point2D project(const Point3D& point) const;

     template<class Visitor>
     void sample(hpreal xEdgeLength, hpreal yEdgeLength, hpuint nx, hpuint ny, Visitor&& visit) const {
          Vector3D xAxis, yAxis;
          std::tie(xAxis, yAxis) = getAxes();
          Space2D::sample<Space3D>(m_origin, xAxis, yAxis, xEdgeLength, yEdgeLength, nx, ny, visit);
     }

     void setNormal(const Vector3D& normal, hpreal epsilon = happah::EPSILON);

     void setOrigin(Point3D origin);

     template<class Vertex = typename QuadMesh3D::VERTEX, class VertexFactory = happah::VertexFactory<Vertex> >
     QuadMesh<Vertex> toQuadMesh(hpreal xEdgeLength = 1.0, hpreal yEdgeLength = 1.0, hpuint nx = 2, hpuint ny = 2, VertexFactory&& factory = VertexFactory()) const {
          std::vector<Vertex> vertices = getVertices(xEdgeLength, yEdgeLength, nx, ny, factory);
          std::vector<hpuint> indices = Space2D::getQuadIndices(nx, ny);
          return QuadMesh<Vertex>(std::move(vertices), std::move(indices));
     }

     template<class Vertex = typename SegmentMesh3D::VERTEX, class VertexFactory = happah::VertexFactory<Vertex> >
     SegmentMesh<Vertex> toSegmentMesh(hpreal xEdgeLength = 1.0, hpreal yEdgeLength = 1.0, hpuint nx = 2, hpuint ny = 2, Diagonals diagonals = NONE, VertexFactory&& factory = VertexFactory()) const {
          std::vector<Vertex> vertices = getVertices(xEdgeLength, yEdgeLength, nx, ny, factory);
          std::vector<hpuint> indices = Space2D::getSegmentIndices(nx, ny, diagonals);
          return SegmentMesh<Vertex>(std::move(vertices), std::move(indices));
     }

     template<class Vertex = typename TriangleMesh3D::VERTEX, class VertexFactory = happah::VertexFactory<Vertex> >
     TriangleMesh<Vertex> toTriangleMesh(hpreal xEdgeLength = 1.0, hpreal yEdgeLength = 1.0, hpuint nx = 3, hpuint ny = 3, VertexFactory&& factory = VertexFactory()) const {
          std::vector<Vertex> vertices = getVertices(xEdgeLength, yEdgeLength, nx, ny, factory);
          std::vector<hpuint> indices = Space2D::getTriangleIndices(nx, ny);
          return TriangleMesh<Vertex>(std::move(vertices), std::move(indices));
     }

     template<class Vertex = VertexP3, class VertexFactory = happah::VertexFactory<Vertex> >
     VertexCloud<Vertex> toVertexCloud(hpreal xEdgeLength = 1.0, hpreal yEdgeLength = 1.0, hpuint nx = 3, hpuint ny = 3, VertexFactory&& factory = VertexFactory()) const {
          std::vector<Vertex> vertices = getVertices(xEdgeLength, yEdgeLength, nx, ny, factory);
          return VertexCloud<Vertex>(std::move(vertices));
     }

     Ray3D unproject(const Point2D& point) const;

private:
     Vector3D m_normal;
     Point3D m_origin;

     static Vector3D check(const Vector3D& normal, hpreal epsilon = happah::EPSILON);

     template<class VertexFactory>
     std::vector<typename VertexFactory::PRODUCT> getVertices(hpreal xEdgeLength, hpreal yEdgeLength, hpuint nx, hpuint ny, VertexFactory& factory) const {
          static_assert(is_vertex<typename VertexFactory::PRODUCT, Space3D>::value, "The conversion methods can only be parameterized by a vertex in 3D space.");
          std::vector<typename VertexFactory::PRODUCT> vertices;
          vertices.reserve(nx * ny);
          sample(xEdgeLength, yEdgeLength, nx, ny, [&](const Point3D& position) { vertices.push_back(factory(position)); });
          return std::move(vertices);
     }

};   
typedef std::shared_ptr<Plane> Plane_ptr;

template<class P, typename = void>
struct is_plane : public std::false_type {};

template<class P>
struct is_plane<P, typename std::enable_if<std::is_base_of<Plane, P>::value>::type> : public std::true_type {};

boost::optional<hpreal> intersect(const Plane& plane, const Ray3D& ray, hpreal epsilon = EPSILON);

Point2D project(const Plane& plane, const Point3D& point);

template<class T>
std::vector<hpreal> make_factors(const std::vector<T>& ts, const Plane& plane, hpreal epsilon = EPSILON) {
     std::vector<hpreal> factors;
     factors.reserve(ts.size());

     for(auto& t : ts) {
          auto ray = make_ray(t);
          ray.normalize();
          if(auto lambda = intersect(plane, ray, epsilon)) {
               auto intersection = ray.getPoint(*lambda);
               factors.push_back(std::sqrt(length2(intersection) / length2(t)));
          } else throw std::runtime_error("Plane does not properly intersect all tetrahedra.");
     }

     return factors;
}

template<class T>
std::vector<Point2D> restrict(const std::vector<T>& ts, const Plane& plane, const std::vector<hpreal>& factors) {
     //TODO: instead of returning Point2Ds return Projection<T> or Plane::Projection<T> or something like that
     std::vector<Point2D> points;
     points.reserve(ts.size());

     auto f = factors.begin();
     for(auto& t : ts) {
          points.push_back(project(plane, mix(t, *f)));
          ++f;
     }

     return points;
}

}//namespace happah

