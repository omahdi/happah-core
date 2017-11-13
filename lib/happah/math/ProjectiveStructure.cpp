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
     auto h = [&](auto x) { return std::accumulate(std::begin(cache), std::end(cache), hpreal(0), [x](auto s, auto a) { return s + (a / std::sqrt(x + a * a)); }); };

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
          auto r = std::tanh(std::asinh(st / std::sin(angle / hpreal(2))) / hpreal(2));

          sum += alpha;
          vertices.emplace_back(r * std::cos(sum), r * std::sin(sum));
          sum += alpha;
     }

     assert(std::abs(sum - glm::two_pi<hpreal>()) < epsilon);

     return vertices;
}

std::vector<Point2D> make_convex_polygon(const Indices& valences, hpreal epsilon) {
     auto angles = std::vector<hpreal>();
     
     angles.reserve(valences.size());
     for(auto& valence : valences) angles.push_back(glm::two_pi<hpreal>() / hpreal(valence));

     return make_convex_polygon(angles);//TODO: avoid intermediate angles array by working on valences directly
}

std::vector<Point3D> make_sun(const Indices& valences, const Indices& pairings) {
     auto n = hpindex(valences.size());
     auto sun = std::vector<Point3D>();
     auto polygon = make_convex_polygon(valences);
     auto d = std::asinh(std::sinh(hpreal(2) * std::atanh(glm::length(polygon[0]))) * std::sin(glm::pi<hpreal>() / hpreal(valences[0])));
     auto r = std::tanh(d);
     auto c = std::cosh(d);
     auto factor = (hpreal(2) * r) / (hpreal(1) + r * r);
     auto sum = hpreal(0);
     auto cache = std::vector<hpreal>();

     auto lambda = [&](auto& x, auto& y, auto& z) {
          auto matrix = glm::inverse(hpmat3x3(x, y, z));
          cache.push_back(matrix[2][0]);
          cache.push_back(matrix[2][1]);
          z *= -matrix[2][2];
          x.z = std::numeric_limits<hpreal>::max();
     };

     sun.reserve(n << 1);
     cache.reserve(n << 1);
     for(auto& point : polygon) sun.emplace_back((hpreal(2) / (hpreal(1) + glm::length2(point))) * point, hpreal(1));
     for(auto& valence : valences) {
          sum += hpreal(2) * std::asin(std::cos(glm::pi<hpreal>() / valence) / c);
          sun.emplace_back(factor * std::cos(sum), factor * std::sin(sum), hpreal(1));
     }
     for(auto i = std::begin(sun), j = i + n, end = j - 1; i != end; ++i, ++j) lambda(*i, *(i + 1), *j);
     sun.front().z = hpreal(1);
     lambda(sun[n - 1], sun.front(), sun.back());
     sun.front().z = std::numeric_limits<hpreal>::max();

     for(auto i = hpindex(0); i != n; ++i) {
          if(sun[i].z != std::numeric_limits<hpreal>::max()) continue;
          auto begin = i;
          auto weight = hpreal(1);
          do {
               weight /= cache[i << 1];
               i = pairings[i];
               weight *= cache[(i << 1) + 1];
               assert(weight > 0);
               if(++i == n) i = hpindex(0);
               auto& point = sun[i];
               assert(point.z == std::numeric_limits<hpreal>::max());
               point.x *= weight;
               point.y *= weight;
               point.z = weight;
          } while(i != begin);
     }

     return sun;
}

hpreal validate(const ProjectiveStructure& structure) {
     auto epsilon = hpreal(0);
     auto& neighbors = structure.getNeighbors();
     auto& transitions = structure.getTransitions();

     visit(make_vertices_enumerator(neighbors), [&](auto t, auto i) {
          auto p1 = Point3D(1, 0, 1);
          auto p2 = Point3D(0, 1, 1);
          visit_spokes(make_spokes_enumerator(neighbors, t, i), [&](auto t, auto i) {
               auto transition = std::begin(transitions) + (9 * t + 3 * i);
               auto p3 = transition[0] * p2 + transition[1] * p1;
               p3.z += transition[2];
               p1 = p2;
               p2 = p3;
          });
          epsilon = std::max(epsilon, std::abs(p2.x));
          epsilon = std::max(epsilon, std::abs(hpreal(1) - p2.y));
          epsilon = std::max(epsilon, std::abs(hpreal(1) - p2.z));
     });

     return epsilon;
}

}//namespace happah

