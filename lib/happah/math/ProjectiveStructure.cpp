// Copyright 2017
//   Obada Mahdi    -                                       - omahdi@gmail.com
//   Pawel Herman   - Karlsruhe Institute of Technology     - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/ProjectiveStructure.hpp"

namespace happah {

std::vector<Point2D> make_convex_polygon(const std::vector<hpreal>& angles, hpreal epsilon) {
     if(std::accumulate(std::begin(angles), std::end(angles), hpreal(0)) > (angles.size() - 2) * glm::pi<hpreal>())
          throw std::runtime_error("Condition for existence of polygon with given angles is not satisfied.");

     auto vertices = std::vector<Point2D>();
     auto cache = std::vector<hpreal>();
     
     //NOTE: We reparameterize g by cosh(x) and g' by cosh(x)^2.
     auto g = [&](auto x) { return std::accumulate(std::begin(cache), std::end(cache), -glm::pi<hpreal>(), [x](auto s, auto a) { return s + std::asin(a / x); }); };
     auto h = [&](auto x) { return std::accumulate(std::begin(cache), std::end(cache), hpreal(0), [x](auto s, auto a) { return s + (a / std::sqrt(x - a * a)); }); };

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
          auto r = st / std::sin(angle / hpreal(2));
          r /= (std::sqrt(r * r + 1.0) + 1.0);
          sum += alpha;
          vertices.emplace_back(r * std::cos(sum), r * std::sin(sum));
          sum += alpha;
     }

     assert(std::abs(sum - glm::two_pi<hpreal>()) < epsilon);

     return vertices;
}

hpreal validate(const ProjectiveStructure& structure) {
     return EPSILON;
     /*auto is_one = [&](auto a) { return glm::abs(1.0 - a) < epsilon; };
     auto is_zero = [&](auto a) { return glm::abs(a) < epsilon; };
     return !find_fan(neighbors, [&](auto p, auto i, auto fan) -> auto {
          auto A1 = hpvec3(0.0, 1.0, 0.0);
          auto A2 = hpvec3(0.0, 0.0, 1.0);
          while(fan) {
               auto q = 0u, j = 0u;
               std::tie(q, j) = *fan;
               auto transition = std::begin(transitions) + (9 * q + 3 * j);
               auto temp = A2;
               if(j == 0) {
                    A2 = transition[0] * A1 + transition[2] * A2;
                    A2.x += transition[1];
               } else if(j == 1) {
                    A2 = transition[0] * A2 + transition[1] * A1;
                    A2.x += transition[2];
               } else {
                    A2 = transition[1] * A2 + transition[2] * A1;
                    A2.x += transition[0];
               }
               A1 = temp;
               ++fan;
          }
          return !(is_zero(A1[0]) && is_one(A1[1]) && is_zero(A1[2]) && is_zero(A2[0]) && is_zero(A2[1]) && is_one(A2[2]));
     });*/
}

}//namespace happah

