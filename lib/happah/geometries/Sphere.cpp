// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/Sphere.h"

namespace happah {

Sphere::Sphere(const Point3D& center, hpreal radius) : m_center(center), m_radius(radius) {}

Sphere::~Sphere() {}

Point2D Sphere::getAbscissa(const Point3D& position) const { return Sphere::Utils::getAbscissa(m_center, position); }

auto Sphere::getCenter() const -> const Point3D& { return m_center; }

Point3D Sphere::getPoint(const Point2D& abscissa) const { return m_center + m_radius * Sphere::Utils::getNormal(abscissa); }

Point3D Sphere::getPoint(const Point2D& abscissa, const Point1D& ordinate) const { return m_center + (m_radius + ordinate.x) * Sphere::Utils::getNormal(abscissa); }

hpreal Sphere::getRadius() const { return m_radius; }

bool Sphere::intersect(const Ray3D& ray, Point3D& intersection, hpreal epsilon) const { return Sphere::Utils::intersect(m_center, m_radius, ray, intersection, epsilon); }

hpuint Sphere::intersect(const Ray3D& ray, hpreal& foot, hpreal& delta, hpreal epsilon) const { return Sphere::Utils::intersect(m_center, m_radius, ray, foot, delta, epsilon); }

hpuint Sphere::intersect(const Ray3D& ray, Point3D& intersection1, Point3D& intersection2, hpreal epsilon) const { return Sphere::Utils::intersect(m_center, m_radius, ray, intersection1, intersection2, epsilon); }

PointCloud3D* Sphere::toPointCloud(hpuint nLatitudes, hpuint nLongitudes) const {
     typedef typename PointCloud3D::VERTEX Vertex;
     std::vector<Vertex>* vertices = sample<Vertex>(nLatitudes, nLongitudes);//TODO: cleanup
     auto pointCloud = new PointCloud3D(*vertices);
     delete vertices;
     return pointCloud;
}

Point2D Sphere::Utils::getAbscissa(const Point3D& position) { return Point2D(std::acos(position.z  / glm::length(position)), std::atan2(position.y, position.x)); }

Point2D Sphere::Utils::getAbscissa(const Point3D& center, const Point3D& position) { return Point2D(std::acos((position.z - center.z) / glm::length(position - center)), std::atan2(position.y - center.y, position.x - center.x)); }

template<class T>
T Sphere::Utils::getCoordinates(hpreal latitude, hpreal longitude) {
     hpreal temp = glm::sin(latitude);
     return T(temp*glm::cos(longitude), temp*glm::sin(longitude), glm::cos(latitude));
}

//NOTE: Resulting vector is normalized.
Vector3D Sphere::Utils::getNormal(hpreal latitude, hpreal longitude) { return Sphere::Utils::getCoordinates<Vector3D>(latitude, longitude); }

//NOTE: Resulting vector is normalized.
Vector3D Sphere::Utils::getNormal(const Point2D& abscissa) { return Sphere::Utils::getNormal(abscissa.x, abscissa.y); }

//NOTE: Resulting position is normalized.
Point3D Sphere::Utils::getPoint(hpreal latitude, hpreal longitude) { return Sphere::Utils::getCoordinates<Point3D>(latitude, longitude); }

//NOTE: Resulting position is normalized.
Point3D Sphere::Utils::getPoint(const Point2D& abscissa) { return Sphere::Utils::getPoint(abscissa.x, abscissa.y); }

Point3D Sphere::Utils::getPoint(const Point2D& abscissa, const Point1D& ordinate) { return (1.0f + ordinate.x) * Sphere::Utils::getNormal(abscissa); }

//NOTE: Ray's direction must be normalized.
hpuint Sphere::Utils::intersect(const Point3D& center, hpreal radius, const Ray3D& ray, hpreal& foot, hpreal& delta, hpreal epsilon) {
     const Point3D& origin = ray.getOrigin();
     const Vector3D& direction = ray.getDirection(); 
     Vector3D v = center - origin;
     foot = glm::dot(v, direction);
     hpreal rsquared = radius * radius;
     hpreal vsquared = glm::dot(v, v);
     bool originOutside = vsquared > rsquared;
     if(foot < 0 && originOutside) return 0;
     hpreal dsquared = vsquared - foot * foot;
     if(dsquared > rsquared) return 0;
     delta = glm::sqrt(rsquared - dsquared);
     if(originOutside) return 2;
     return 1;
}

//NOTE: Ray's direction must be normalized.
bool Sphere::Utils::intersect(const Point3D& center, hpreal radius, const Ray3D& ray, Point3D& intersection, hpreal epsilon) {
     hpreal foot, delta;
     hpuint result = Sphere::Utils::intersect(center, radius, ray, foot, delta, epsilon);
     hpreal t0;
     if(result > 0) {
           if(result == 1) t0 = foot + delta; 
           else t0 = foot - delta;
           intersection = ray.getPoint(t0);
           return true;
     }
     return false;
}

//NOTE: Ray's direction must be normalized.
hpuint Sphere::Utils::intersect(const Point3D& center, hpreal radius, const Ray3D& ray, Point3D& intersection1, Point3D& intersection2, hpreal epsilon) {
     hpreal foot, delta;
     hpuint result = Sphere::Utils::intersect(center, radius, ray, foot, delta, epsilon);
     if(result == 0) return 0;
     hpreal t0;
     if(result == 1) {
          t0 = foot + delta;     
          intersection1 = ray.getPoint(t0);
     } else {
          t0 = foot - delta;
          hpreal t1 = foot + delta;
          intersection1 = ray.getPoint(t0);
          intersection2 = ray.getPoint(t1);;
     }
     return result;
}

bool Sphere::Utils::isSphere(const Point3D& p1, const Point3D& p2, const Point3D& p3, hpreal radius, Point3D& center1, Point3D& center2) {
     //TODO
     return false;
}

}//namespace happah

