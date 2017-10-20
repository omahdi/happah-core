// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>
#include <vector>

#include "happah/Happah.hpp"
#include "happah/math/Space.hpp"

namespace happah {

//DECLARATIONS

class Circle;

//NOTE: In the case of two intersections, intersect returns the left intersection point first and the right second from the perspective at center0 looking at center1.
boost::optional<boost::variant<Point2D, std::tuple<Point2D, Point2D> > > intersect(const Circle& circle0, const Circle& circle1, hpreal epsilon = EPSILON);

inline Circle make_circle(Point2D center, hpreal radius);

std::vector<Point2D> make_regular_polygon(hpreal radius, hpuint n);

Circle poincare_to_euclidean(const Circle& circle);

template<class Visitor>
void sample(const Circle& circle, hpuint n, Visitor&& visit);

//DEFINITIONS

class Circle {
public:
     Circle(Point2D center, hpreal radius)
          : m_center(std::move(center)), m_radius(radius) {}

     const Point2D& getCenter() const { return m_center; }

     hpreal getRadius() const { return m_radius; }

private:
     Point2D m_center;
     hpreal m_radius;

};//Circle

inline Circle make_circle(Point2D center, hpreal radius) { return { std::move(center), radius }; }

template<class Visitor>
void sample(const Circle& circle, hpuint n, Visitor&& visit) {
     auto rotation = glm::angleAxis(glm::two_pi<hpreal>() / n, Vector3D(0, 0, 1));
     auto point = Point2D(hpreal(0), circle.getRadius());

     visit(point);
     while(--n) {
          point = Point2D(glm::rotate(rotation, Point3D(point, 0)));
          visit(point);
     }
}

}//namespace happah

