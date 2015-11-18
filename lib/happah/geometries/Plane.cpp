// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#define GLM_FORCE_RADIANS

#include <glm/gtx/quaternion.hpp>
#include <iostream>
#include <math.h>

#include "happah/geometries/Plane.h"

//TODO: typedef glm stuff into hph namespace
//TODO: create hph namespace?
//TODO: use ux=u_0 format for plane? keep normals normalized for distance function?

Plane::Plane(Point3D origin, const Vector3D& normal, hpreal epsilon) 
     : m_normal(check(normal, epsilon)), m_origin(std::move(origin)) {}

auto Plane::check(const Vector3D& normal, hpreal epsilon) -> Vector3D {
     if(glm::abs(normal.x) < epsilon && glm::abs(normal.y) < epsilon && glm::abs(normal.z) < epsilon)
          std::cerr << "Plane with normal == 0 generated!" << std::endl;
     return glm::normalize(normal);
}

auto Plane::getAxes() const -> std::pair<Vector3D, Vector3D> {
     Vector3D normal = glm::normalize(m_normal);
     Vector3D xAxis = Vector3D(1.f, 0.f, 0.f);
     if(xAxis == normal) xAxis = Vector3D(0.f, 0.f, 1.f);//TODO: to epsilon
     xAxis -= glm::dot(xAxis, normal) * normal;
     Vector3D yAxis = glm::normalize(glm::cross(normal, xAxis));
     return std::make_pair(xAxis, yAxis);
}

auto Plane::getNormal() const -> const Vector3D& { return m_normal; }

auto Plane::getOrigin() const -> const Point3D& { return m_origin; }

//NOTE: Returns the transformation that transforms the x-y plane to this plane.
//TODO: check that this transformation transforms the x-axis, y-axis to the x-axis, y-axis of this plane
RigidAffineTransformation3D Plane::getPlaneTransformation(hpreal epsilon) const {//TODO: rename to something better
     if(glm::abs(m_normal.z-1.0) < epsilon) return RigidAffineTransformation3D(hpmat3x3(1.0), m_origin);
     glm::quat rotation(glm::degrees(std::acos(m_normal.z)),glm::normalize(Vector3D(m_normal.z-m_normal.y, m_normal.x-m_normal.z, 0.0)));
     return RigidAffineTransformation3D(glm::toMat3(rotation), m_origin);
}

boost::optional<hpreal> Plane::intersect(const Ray3D& ray, hpreal epsilon) const { return Plane::intersect(m_origin, m_normal, ray, epsilon); }

//NOTE: Ray's direction must be normalized.
boost::optional<hpreal> Plane::intersect(const Point3D& origin, const Vector3D& normal, const Ray3D& ray, hpreal epsilon) {
     const Vector3D& direction = ray.getDirection(); 
     hpreal a = glm::dot(direction, normal);
     if(glm::abs(a) < epsilon) return boost::none;
     const Point3D& temp = ray.getOrigin();
     return (hpreal) ((glm::dot(origin, normal) - glm::dot(temp, normal)) / a);
}

Point2D Plane::project(const Point3D& point) const {
     //TODO: can this be done more cheaply?
     /*typedef typename RigidAffineTransformation3D::MATRIX Matrix;
     const Matrix& matrix = getPlaneTransformation().getMatrix();
     Vector3D u(matrix[0][0], matrix[1][0], matrix[2][0]);
     Vector3D v(matrix[0][1], matrix[1][1], matrix[2][1]);
     Point3D temp = point - m_origin;
     return Point2D(glm::dot(temp, u), glm::dot(temp, v));*/
     Vector3D xAxis, yAxis;
     std::tie(xAxis, yAxis) = getAxes();
     Point3D temp = point - m_origin;
     return Point2D(glm::dot(temp, xAxis), glm::dot(temp, yAxis));
}

void Plane::setNormal(const Vector3D& normal, hpreal epsilon) { m_normal = check(normal, epsilon); }

void Plane::setOrigin(Vector3D origin) { m_origin = std::move(origin); }

Ray3D Plane::unproject(const Point2D& point) const {
     typedef typename RigidAffineTransformation3D::MATRIX Matrix;
     const Matrix& matrix = getPlaneTransformation().getMatrix();
     Vector3D u(matrix[0][0], matrix[1][0], matrix[2][0]);
     Vector3D v(matrix[0][1], matrix[1][1], matrix[2][1]);
     return Ray3D(point.x * u + point.y * v + m_origin, m_normal);
}

Plane::Utils::PointFactory::PointFactory(const Plane& plane)
     : m_plane(plane) { std::tie(m_xAxis, m_yAxis) = m_plane.getAxes(); }

Point3D Plane::Utils::PointFactory::build(const Point2D& abscissa, const Point1D& ordinate) const { return m_plane.getOrigin() + abscissa.x * m_xAxis + abscissa.y * m_yAxis + ordinate.x * m_plane.getNormal(); }

//NOTE: pawned from http://missingbytes.blogspot.de/2012/06/fitting-plane-to-point-cloud.html
Vector3D Plane::Utils::getBestFittingPlaneNormal(const std::vector<Point3D>& points, hpreal epsilon) {
     static const hpuint ITERATIONS = 100;

     //TODO: make sure algorithm is correct
     hpreal sumxx = 0, sumxy = 0, sumxz = 0, sumyy = 0, sumyz = 0, sumzz = 0;
     for(Point3D p : points) {
          sumxx += p.x * p.x;
          sumxy += p.x * p.y;
          sumxz += p.x * p.z;
          sumyy += p.y * p.y;
          sumyz += p.y * p.z;
          sumzz += p.z * p.z;
     }
     hpmat3x3 m(sumxx, sumxy, sumxz, sumxy, sumyy, sumyz, sumxz, sumyz, sumzz);
     //TODO: check for det = 0;
     hpmat3x3 mInverse = glm::inverse(m);
     Vector3D eigenvector(1.0,1.0,1.0);
     for(hpuint i = 0; i < ITERATIONS; ++i) eigenvector = glm::normalize(mInverse * eigenvector);
     //TODO: why check z?
     if(glm::abs(eigenvector.z) < epsilon) return Vector3D(0.0, 0.0, -1.0);
     return eigenvector;
}

