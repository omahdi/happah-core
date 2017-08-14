// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/CirclePackingMetric.hpp"

namespace happah {

CirclePackingMetric make_circle_packing_metric(std::vector<hpreal> weights, Indices indices, const Indices& neighbors, const Indices& border, hpreal epsilon) {
     auto radii = std::vector<hpreal>(indices.size(), hpreal(1.0));

     auto is_border = [&](auto t, auto i) -> bool {
          auto result = false;
          visit(make_spokes_enumerator(neighbors, t, i), [&](auto t, auto i) { result |= std::binary_search(std::begin(border), std::end(border), 3 * t + i); });
          return result;
     };

     auto max = hpreal(0);
     auto counter = hpuint(0);
     do {
          max = std::numeric_limits<hpreal>::min();
          visit_vertices(neighbors, [&](auto t, auto i) {
               static constexpr hpuint o0[3] = { 1, 2, 0 };
               static constexpr hpuint o1[3] = { 2, 0, 1 };

               if(is_border(t, i)) return;
               auto sum = hpreal(0);
               auto& r0 = radii[indices[3 * t + i]];
               visit(make_spokes_enumerator(neighbors, t, i), [&](auto t, auto i) {
                    auto r1 = radii[indices[3 * t + o1[i]]];
                    auto r2 = radii[indices[3 * t + o0[i]]];
                    auto l0 = std::acosh(std::cosh(r0) * std::cosh(r2));
                    auto l1 = std::acosh(std::cosh(r0) * std::cosh(r1));
                    auto l2 = std::acosh(std::cosh(r1) * std::cosh(r2));
                    sum += std::acos((std::cosh(l0) * std::cosh(l1) - std::cosh(l2)) / (std::sinh(l0) * std::sinh(l1)));
               });
               auto temp = sum - glm::two_pi<hpreal>();
               if(std::abs(temp) > epsilon) {
                    if(std::abs(temp) > max) max = std::abs(temp);
                    if(temp > 0) r0 += epsilon * r0;
                    else r0 -= epsilon * r0;
               }
               assert(r0 > hpreal(0));
          });
          if(counter == 100) {
               std::cout << max << std::endl;
               counter = 0;
          } else ++counter;
     } while(max > epsilon);

     return make_circle_packing_metric(std::move(radii), std::move(weights), std::move(indices));
}

}//namespace happah

