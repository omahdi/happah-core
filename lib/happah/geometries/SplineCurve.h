// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <algorithm>
#include <vector>

namespace happah {

template<class Space, hpuint t_degree>
class SplineCurve : public Geometry1D<Space> {
     using Point = typename Space::POINT;

public:
     SplineCurve() {}

     SplineCurve(std::vector<Point> controlPoints, std::vector<hpreal> knots)
          : m_controlPoints{std::move(controlPoints)}, m_knots{std::move(knots)} {}

     auto& getControlPoints() const { return m_controlPoints; }

     auto& getKnots() const { return m_knots; }

private:
     std::vector<Point> m_controlPoints;
     std::vector<hpreal> m_knots;

     template<class Stream>
     friend auto& operator<<(Stream& stream, const SplineCurve<Space, t_degree>& curve) {
          stream << curve.m_controlPoints << '\n';
          stream << curve.m_knots;
          return stream;
     }

     template<class Stream>
     friend auto& operator>>(Stream& stream, SplineCurve<Space, t_degree>& curve) {
          stream >> curve.m_controlPoints;
          stream >> curve.m_knots;
          return stream;
     }

};//SplineCurve
template<class Space>
using ConstantSplineCurve = SplineCurve<Space, 0>;
template<class Space>
using CubicSplineCurve = SplineCurve<Space, 3>;
template<class Space>
using LinearSplineCurve = SplineCurve<Space, 1>;
template<class Space>
using SexticSplineCurve = SplineCurve<Space, 6>;
template<class Space>
using QuadraticSplineCurve = SplineCurve<Space, 2>;
template<class Space>
using QuarticSplineCurve = SplineCurve<Space, 4>;
template<class Space>
using QuinticSplineCurve = SplineCurve<Space, 5>;

template<class Space, hpuint degree>
SplineCurve<Space, degree> insert_knot(const SplineCurve<Space, degree>& curve, hpreal a) {
     using Point = typename Space::POINT;

     auto& controlPoints0 = curve.getControlPoints();
     auto& knots0 = curve.getKnots();
     auto temp = std::upper_bound(knots0.begin(), knots0.end(), a);
     auto n = std::distance(knots0.begin(), temp);
     auto c = controlPoints0.begin() + (n - degree);

     auto controlPoints1 = std::vector<Point>();
     controlPoints1.reserve(controlPoints0.size() + 1);
     controlPoints1.insert(controlPoints1.end(), controlPoints0.begin(), c + 1);
     auto l = temp - degree;
     auto r = temp;
     for(auto end = c + degree; c != end; ++c, ++l, ++r) {
          auto alpha = (a - *l) / (*r - *l);
          controlPoints1.push_back((1 - alpha) * (*c) + alpha * (*(c + 1)));
     }
     controlPoints1.insert(controlPoints1.end(), c, controlPoints0.end());

     auto knots1 = std::vector<hpreal>();
     knots1.reserve(knots0.size() + 1);
     knots1.insert(knots1.end(), knots0.begin(), temp);
     knots1.push_back(a);
     knots1.insert(knots1.end(), temp, knots0.end());

     return { controlPoints1, knots1 };
}

}//namespace happah

