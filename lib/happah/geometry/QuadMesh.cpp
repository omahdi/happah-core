// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometry/QuadMesh.hpp"
#include "happah/util/visitors.hpp"

namespace happah {

Indices make_neighbors(quads, const Indices& indices) {
     auto map = make_map(0);
     auto neighbors = Indices();
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
     visit_quartets(indices, [&](auto v0, auto v1, auto v2, auto v3) {
          cache(v0, v1);
          cache(v1, v2);
          cache(v2, v3);
          cache(v3, v0);
          ++q;
     });

     q = hpindex(0);
     visit_quartets(indices, [&](auto v0, auto v1, auto v2, auto v3) {
          move(v0, v1);
          move(v1, v2);
          move(v2, v3);
          move(v3, v0);
          ++q;
     });

     return neighbors;
}

quat make_quad_neighbor_offset(const Indices& neighbors, hpindex q, hpindex r) {
     auto n = std::begin(neighbors) + (q << 2);
     
     if(r == n[0]) return QUAT0;
     if(r == n[1]) return QUAT1;
     if(r == n[2]) return QUAT2;
     assert(r == n[3]);
     return QUAT3;
}

quat make_quad_vertex_offset(const Indices& indices, hpindex q, hpindex v) {
     auto i = std::begin(indices) + (q << 2);

     if(v == i[0]) return QUAT0;
     if(v == i[1]) return QUAT1;
     if(v == i[2]) return QUAT2;
     assert(v == i[3]);
     return QUAT3;
}

}//namespace happah
