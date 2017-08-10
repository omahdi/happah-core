// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <glm/gtx/quaternion.hpp>

#include "happah/geometry/Plane.hpp"

namespace happah {

//TODO: typedef glm stuff into hph namespace
//TODO: use ux=u_0 format for plane? keep normals normalized for distance function?

//NOTE: Ray's direction must be normalized.
boost::optional<hpreal> intersect(const Plane& plane, const Ray3D& ray, hpreal epsilon) {
     auto& normal = plane.getNormal();
     auto& direction = ray.getDirection(); 
     auto a = glm::dot(direction, normal);
     if(glm::abs(a) < epsilon) return boost::none;
     return hpreal((glm::dot(plane.getOrigin(), normal) - glm::dot(ray.getOrigin(), normal)) / a);
}

std::tuple<Vector3D, Vector3D> make_axes(const Plane& plane, hpreal epsilon) {
     auto normal = glm::normalize(plane.getNormal());
     auto xAxis = Vector3D(1.0, 0.0, 0.0);
     if(glm::length2(xAxis - normal) < epsilon) xAxis = Vector3D(0.0, 0.0, 1.0);
     xAxis -= glm::dot(xAxis, normal) * normal;
     auto yAxis = glm::normalize(glm::cross(normal, xAxis));
     return std::make_tuple(xAxis, yAxis);
}

//TODO: check that this transformation transforms the x-axis, y-axis to the x-axis, y-axis of this plane
std::tuple<hpmat3x3, hpvec3> make_transformation(const Plane& plane, hpreal epsilon) {
     auto& normal = plane.getNormal();
     auto& origin = plane.getOrigin();
     if(glm::abs(normal.z - hpreal(1.0)) < epsilon) return std::make_tuple(hpmat3x3(1.0), origin);
     auto rotation = glm::quat(glm::degrees(std::acos(normal.z)), glm::normalize(Vector3D(normal.z - normal.y, normal.x - normal.z, hpreal(0.0))));
     return std::make_tuple(glm::toMat3(rotation), origin);
}

Point2D project(const Plane& plane, const Point3D& point) {
     //TODO: can this be done more cheaply?
     /*typedef typename RigidAffineTransformation3D::MATRIX Matrix;
     auto matrix = std::get<0>(getPlaneTransformation());
     Vector3D u(matrix[0][0], matrix[1][0], matrix[2][0]);
     Vector3D v(matrix[0][1], matrix[1][1], matrix[2][1]);
     Point3D temp = point - m_origin;
     return Point2D(glm::dot(temp, u), glm::dot(temp, v));*/
     auto xAxis = Vector3D(), yAxis = Vector3D();
     auto temp = point - plane.getOrigin();
     std::tie(xAxis, yAxis) = make_axes(plane);
     return { glm::dot(temp, xAxis), glm::dot(temp, yAxis) };
}

Ray3D unproject(const Plane& plane, const Point2D& point) {
     auto matrix = std::get<0>(make_transformation(plane));
     auto u = Vector3D(matrix[0][0], matrix[1][0], matrix[2][0]);
     auto v = Vector3D(matrix[0][1], matrix[1][1], matrix[2][1]);
     return { point.x * u + point.y * v + plane.getOrigin(), plane.getNormal() };
}

}//namespace happah

