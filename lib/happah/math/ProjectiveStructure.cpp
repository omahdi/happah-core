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

ProjectiveStructure make_projective_structure(const Indices& valences, const Indices& pairings) {
     auto w = hpreal(0);
     auto sun = std::vector<Point2D>();
     std::tie(sun, w) = make_sun(valences);
     auto n = hpindex(valences.size());
     auto neighbors = Indices();
     auto transitions = std::vector<hpreal>();
     auto center = hpvec3(0, 0, 1);
     auto i = hpindex(-1);
     auto x = std::begin(sun);
     auto y = x + n;
     auto end = x + (n - 2);

     auto lambda0 = [&](auto& p0, auto& p1, auto& p2, auto& q0) {
          auto transition0 = glm::inverse(hpmat3x3(Point3D(p1, 1), Point3D(w * q0, w), Point3D(p0, 1)))[2];
          auto transition1 = glm::inverse(hpmat3x3(center, Point3D(p2, 1), Point3D(p1, 1))) * Point3D(p0, 1);

          transitions.insert(std::end(transitions), {
               transition0.x, transition0.y, transition0.z,
               transition1.x, transition1.y, transition1.z,
               -transition1.z / transition1.y, hpreal(1) / transition1.y, -transition1.x / transition1.y
          });
     };
     auto lambda1 = [&](auto& p0, auto& p1, auto& p2, auto& q0) {
          auto transition0 = glm::inverse(hpmat3x3(Point3D(p1, 1), Point3D(w * q0, w), Point3D(p0, 1)))[2];
          auto transition1 = glm::inverse(hpmat3x3(center, Point3D(p2, 1), Point3D(p1, 1))) * Point3D(p0, 1);

          transitions.insert(std::end(transitions), {
               transition0.x, transition0.y, transition0.z,
               transition1.x, transition1.y, transition1.z
          });
          transitions[0] = -transition1.z / transition1.y;
          transitions[1] = hpreal(1) / transition1.y;
          transitions[2] = -transition1.x / transition1.y;
     };

     neighbors.reserve(3 * n);
     transitions.reserve(9 * n);
     for(auto j : pairings) {
          neighbors.push_back(i);
          neighbors.push_back(j);
          neighbors.push_back(i + 2);
          ++i;
     }
     neighbors.front() = hpindex(n - 1);
     neighbors.back() = hpindex(0);

     transitions.push_back(std::numeric_limits<hpreal>::max());
     transitions.push_back(std::numeric_limits<hpreal>::max());
     transitions.push_back(std::numeric_limits<hpreal>::max());
     while(x != end) {
          lambda0(x[0], x[1], x[2], y[0]);
          ++x;
          ++y;
     }
     lambda0(x[0], x[1], sun.front(), y[0]);
     lambda1(x[1], sun.front(), sun[1], y[1]);

     /*auto test = [&](auto transition0, auto transition1) {
          std::cout << transition1[1] * transition0[0] + transition1[2] << '\n';
          std::cout << transition1[1] * transition0[1] << '\n';
          std::cout << transition1[1] * transition0[2] + transition1[0] << '\n' << '\n';
     };

     boost::dynamic_bitset<> done(valences.size(), false);
     for(auto j : pairings) {
          if(done[j]) continue;
          auto begin = j;
          do {
               auto i = pairings[j];
               auto transition0 = std::begin(transitions) + (9 * i + 3);
               auto transition1 = std::begin(transitions) + (9 * j + 3);
               test(transition0, transition1);
               done[j] = true;
               i = j + 1;
               if(i == valences.size()) i = hpindex(0);
               j = pairings[i];
          } while(j != begin);
     }*/

     return make_projective_structure(std::move(neighbors), std::move(transitions));
}

std::tuple<std::vector<Point2D>, hpreal> make_sun(const Indices& valences) {
     auto n = hpindex(valences.size());
     auto sun = make_convex_polygon(valences);
     auto d = std::asinh(std::sinh(hpreal(2) * std::atanh(glm::length(sun[0]))) * std::sin(glm::pi<hpreal>() / hpreal(valences[0])));
     auto r0 = std::tanh(d / hpreal(2));
     auto r1 = std::tanh(d);
     auto c = std::cosh(d);
     auto l0 = (hpreal(2) * r0) / (hpreal(1) + r0 * r0);
     auto l1 = (hpreal(2) * r1) / (hpreal(1) + r1 * r1);
     auto sum = hpreal(0);

     sun.reserve(n << 1);
     for(auto& point : sun) point *= hpreal(2) / (hpreal(1) + glm::length2(point));
     for(auto& valence : valences) {
          sum += hpreal(2) * std::asin(std::cos(glm::pi<hpreal>() / valence) / c);
          sun.emplace_back(l1 * std::cos(sum), l1 * std::sin(sum));
     }

     return std::make_tuple(std::move(sun), l0 / (l1 - l0));
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

