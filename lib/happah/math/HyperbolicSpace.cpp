// Copyright 2017
//   Obada Mahdi    -                                       - omahdi@gmail.com
//   Pawel Herman   - Karlsruhe Institute of Technology     - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/HyperbolicSpace.hpp"

namespace happah {

std::vector<Point2D> make_convex_polygon(const std::vector<hpreal>& angles, hpreal epsilon) {
     if(std::accumulate(std::begin(angles), std::end(angles), hpreal(0)) > (angles.size() - 2) * glm::pi<hpreal>())
          throw std::runtime_error("Condition for existence of polygon with given angles is not satisfied.");

     auto vertices = std::vector<Point2D>();
     auto cache = std::vector<hpreal>();
     
     //NOTE: We reparameterize g by cosh(x) and g' by cosh(x)^2.
     auto g = [&](auto x) { return std::accumulate(std::begin(cache), std::end(cache), -glm::pi<hpreal>(), [x](auto s, auto a) { return s + std::asin(a / x); }); };
     auto h = [&](auto x) { return std::accumulate(std::begin(cache), std::end(cache), hpreal(0), [x](auto s, auto a) { return s + (1.0 / std::sqrt(x / (a * a) - 1.0)); }); };

     vertices.reserve(angles.size());
     cache.reserve(angles.size());
     for(auto& angle : angles) cache.push_back(std::cos(angle / hpreal(2)));

     auto t = hpreal(1);
     auto u = std::cosh(t);
     auto y = g(u);
     while(std::abs(y) > epsilon) {
          t -= y / (-std::tanh(t) * h(u * u));
          u = std::cosh(t);
          y = g(u);
     }

     auto ct = std::cosh(t);
     auto st = std::sinh(t);
     auto sum = hpreal(0);
     for(auto& angle : angles) {
          auto alpha = std::asin(std::cos(angle / hpreal(2)) / ct);
          auto r = st / std::sin(angle);
          r /= (std::sqrt(r * r + 1.0) + 1.0);
          sum += alpha;
          vertices.emplace_back(r * std::cos(sum), r * std::sin(sum));
          sum += alpha;
     }

     assert(std::abs(sum - glm::two_pi<hpreal>()) < epsilon);

     return vertices;
}

}//namespace happah

