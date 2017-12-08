// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/CirclePacking.hpp"

namespace happah {

hpreal angle_sum(const CirclePacking& packing, const Indices& neighbors, hpindex t, hptrit i) {
     auto sum = hpreal(0);
     auto r0 = packing.getRadius(t, i);

     visit(make_spokes_enumerator(neighbors, t, i), [&](auto t, auto i) {
          static constexpr hpuint o0[3] = { 1, 2, 0 };
          static constexpr hpuint o1[3] = { 2, 0, 1 };

          auto r1 = packing.getRadius(t, hptrit(o1[i]));
          auto r2 = packing.getRadius(t, hptrit(o0[i]));
          auto l0 = std::acosh(std::cosh(r0) * std::cosh(r2));
          auto l1 = std::acosh(std::cosh(r0) * std::cosh(r1));
          auto l2 = std::acosh(std::cosh(r1) * std::cosh(r2));
          sum += std::acos((std::cosh(l0) * std::cosh(l1) - std::cosh(l2)) / (std::sinh(l0) * std::sinh(l1)));
     });

     return sum;
}

CirclePacking make_circle_packing(std::vector<hpreal> weights, Indices indices, const Indices& neighbors, hpreal epsilon) {
     auto radii = std::vector<hpreal>(indices.size(), hpreal(1.0));
     auto packing = make_circle_packing(std::move(radii), std::move(weights), std::move(indices));

     /*auto is_border = [&](auto t, auto i) -> bool {
          auto result = false;
          visit(make_spokes_enumerator(neighbors, t, i), [&](auto t, auto i) { result |= std::binary_search(std::begin(border), std::end(border), 3 * t + i); });
          return result;
     };*/

     auto max = hpreal(0);
     auto counter = hpuint(0);
     do {
          max = std::numeric_limits<hpreal>::min();
          visit(make_vertices_enumerator(neighbors), [&](auto t, auto i) {
               //if(is_border(t, i)) return;
               auto sum = angle_sum(packing, neighbors, t, i);
               auto temp = std::abs(sum - glm::two_pi<hpreal>());
               if(temp > epsilon) {
                    if(temp > max) max = temp;
                    auto r = packing.getRadius(t, i);
                    packing.setRadius(t, i, (sum > glm::two_pi<hpreal>()) ? (hpreal(1.0) + epsilon) * r : (hpreal(1.0) - epsilon) * r);
               }
          });
          if(counter == 100) {
               std::cout << max << std::endl;
               counter = 0;
          } else ++counter;
     } while(max > epsilon);

     return packing;
}

hpreal validate(const CirclePacking& packing, const Indices& neighbors) {
     auto max  = std::numeric_limits<hpreal>::min();

     visit(make_vertices_enumerator(neighbors), [&](auto t, auto i) {
          auto sum = angle_sum(packing, neighbors, t, i);
          auto temp = std::abs(sum - glm::two_pi<hpreal>());

          if(temp > max) max = temp;
     });

     return max;
}

}//namespace happah

