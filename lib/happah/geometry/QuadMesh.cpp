// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometry/QuadMesh.hpp"

namespace happah {

Quadruples<quax> make_neighbors(const Quadruples<hpindex>& indices) {
     auto map = make_map<std::pair<quax, quax> >(0);
     auto neighbors = Quadruples<quax>();
     auto q = hpindex(0);
     
     auto cache = [&](auto va, auto vb, auto i) {
          auto key = std::make_pair(va, vb);
          auto temp = map.find(key);

          if(temp == map.end()) map[key] = std::make_pair(quax(q, i), quax());
          else temp->second.second = quax(q, i);
     };
     
     auto move = [&](auto va, auto vb) {
          auto value = map[{ va, vb }];

          if(value.first.getQuadruple() == q) neighbors.push_back(value.second);
          else neighbors.push_back(value.first);
     };

     neighbors.reserve(indices.size());

     q = hpindex(0);
     visit(indices, [&](auto v0, auto v1, auto v2, auto v3) {
          cache(v0, v1, QUAT0);
          cache(v1, v2, QUAT1);
          cache(v2, v3, QUAT2);
          cache(v3, v0, QUAT3);
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

