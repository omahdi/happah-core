// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometry/QuadMesh.hpp"

namespace happah {

Quartets<hpindex> make_neighbors(const Quartets<hpindex>& indices) {
     auto map = make_map(0);
     auto neighbors = Quartets<hpindex>();
     auto q = hpindex(0);
     
     auto cache = [&](auto va, auto vb) {
          auto key = std::make_pair(va, vb);
          auto i = map.find(key);

          if(i == map.end()) map[key] = std::make_pair(q, std::numeric_limits<hpindex>::max());
          else i->second.second = q;
     };
     
     auto move = [&](auto va, auto vb) {
          auto value = map[{ va, vb }];

          if(value.first == q) neighbors.push_back(value.second);
          else neighbors.push_back(value.first);
     };

     neighbors.reserve(indices.size());

     q = hpindex(0);
     visit(indices, [&](auto v0, auto v1, auto v2, auto v3) {
          cache(v0, v1);
          cache(v1, v2);
          cache(v2, v3);
          cache(v3, v0);
          ++q;
     });

     q = hpindex(0);
     visit(indices, [&](auto v0, auto v1, auto v2, auto v3) {
          move(v0, v1);
          move(v1, v2);
          move(v2, v3);
          move(v3, v0);
          ++q;
     });

     return neighbors;
}

}//namespace happah

