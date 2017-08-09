// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/optional.hpp>
#include <boost/range/irange.hpp>
#include <boost/variant.hpp>
#include <glm/gtc/constants.hpp>

#include "happah/Happah.hpp"
#include "happah/geometries/Ray.hpp"
#include "happah/geometries/TriangleMesh.hpp"
#include "happah/utils/VertexFactory.hpp"

namespace happah {

//DECLARATIONS

class Sphere;

//NOTE: Ray's direction must be normalized.
boost::optional<boost::variant<Point3D, std::tuple<Point3D, Point3D> > > intersect(const Sphere& sphere, const Ray3D& ray, hpreal epsilon = EPSILON);

inline Point2D make_abscissa(const Sphere& sphere, const Point3D& point);

inline Point3D make_point(const Sphere& sphere, hpreal latitude, hpreal longitude);

inline Sphere make_sphere();

inline Sphere make_sphere(Point3D center, hpreal radius);

template<class Vertex = VertexP3, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const Sphere& sphere, hpuint nLatitudes = 50, hpuint nLongitudes = 50, VertexFactory&& build = VertexFactory());

std::vector<Point3D> sample(const Sphere& sphere, hpuint nLatitudes, hpuint nLongitudes);

//DEFINITIONS

class Sphere {
public:
     Sphere(Point3D center, hpreal radius)
          : m_center(std::move(center)), m_radius(radius) {}

     auto& getCenter() const { return m_center; }

     auto getRadius() const { return m_radius; }

private:
     Point3D m_center;
     hpreal m_radius;

};//Sphere

inline Point2D make_abscissa(const Sphere& sphere, const Point3D& point) {
     auto& center = sphere.getCenter();
     return Point2D(std::acos((point.z - center.z) / glm::length(point - center)), std::atan2(point.y - center.y, point.x - center.x));
}

inline Point3D make_point(const Sphere& sphere, hpreal latitude, hpreal longitude) {
     auto temp = glm::sin(latitude);
     return sphere.getCenter() + sphere.getRadius() * Point3D(temp * glm::cos(longitude), temp * glm::sin(longitude), glm::cos(latitude));
}

inline Sphere make_sphere() { return { Point3D(0.0), hpreal(1.0) }; }

inline Sphere make_sphere(Point3D center, hpreal radius) { return { std::move(center), radius }; }

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const Sphere& sphere, hpuint nLatitudes, hpuint nLongitudes, VertexFactory&& build) {
     auto indices = Indices();
     auto points = sample(sphere, nLatitudes, nLongitudes);
     auto vertices = std::vector<Vertex>();
     auto n = hpuint(points.size() - 1);
     
     vertices.reserve(points.size());
     indices.reserve(6 * nLongitudes * nLatitudes);

     for(auto& point : points) vertices.push_back(build(point));

     for(auto j : boost::irange(hpuint(1), nLongitudes)) {
          indices.push_back(j);
          indices.push_back(0);
          indices.push_back(j + 1);
     }
     indices.push_back(nLongitudes);
     indices.push_back(0);
     indices.push_back(1);

     for(auto j : boost::irange(n - nLongitudes, nLatitudes * nLongitudes)) {
          indices.push_back(n);
          indices.push_back(j);
          indices.push_back(j + 1);
     }
     indices.push_back(n);
     indices.push_back(nLatitudes * nLongitudes);
     indices.push_back(n - nLongitudes);

     for(auto i : boost::irange(hpindex(1), n - nLongitudes)) {
          indices.push_back(i);
          indices.push_back(i + nLongitudes);
          indices.push_back(i + nLongitudes - 1);
          indices.push_back(i);
          indices.push_back(i + 1);
          indices.push_back(i + nLongitudes);
     }

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

//TODO: take advantage of sphere symmetry (8 quadrants)
std::vector<Point3D> sample(const Sphere& sphere, hpuint nLatitudes, hpuint nLongitudes) {
     auto& center = sphere.getCenter();
     auto radius = sphere.getRadius();
     auto points = std::vector<Point3D>();
     auto phi = glm::pi<hpreal>() / (nLatitudes + 1);
     auto theta = 2.0 * glm::pi<hpreal>() / nLongitudes;

     points.reserve(nLatitudes * nLongitudes + 2);

     points.emplace_back(center.x, center.y, center.z + radius);
     for(auto i : boost::irange(hpuint(1), nLatitudes + hpuint(1))) {
          auto alpha = glm::half_pi<hpreal>() - i * phi;
          auto beta = std::cos(alpha);
          auto z = radius * std::sin(alpha);
          for(auto j : boost::irange(hpuint(0), nLongitudes)) {
               auto gamma = j * theta;
               auto x = radius * beta * std::cos(gamma);
               auto y = radius * beta * std::sin(gamma);
               points.emplace_back(x + center.x, y + center.y, z + center.z);
          }
     }
     points.emplace_back(center.x, center.y, center.z - radius);

     return points;
}

}//namespace happah

