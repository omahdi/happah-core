// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometry/TriangleMesh.hpp"

namespace happah {

Triplets<hpindex> make_neighbors(const Triplets<hpindex>& indices) {
     auto map = make_map<std::pair<hpindex, hpindex> >(0);
     auto neighbors = Triplets<hpindex>();
     auto t = hpindex(0);

     auto cache = [&](auto va, auto vb) {
          auto key = std::make_pair(va, vb);
          auto i = map.find(key);

          if(i == map.end()) map[key] = std::make_pair(t, std::numeric_limits<hpindex>::max());
          else i->second.second = t;
     };

     auto move = [&](auto va, auto vb) {
          auto value = map[{ va, vb }];

          if(value.first == t) neighbors.push_back(value.second);
          else neighbors.push_back(value.first);
     };

     neighbors.reserve(indices.size());

     t = hpindex(0);
     visit(indices, [&](auto v0, auto v1, auto v2) {
          cache(v0, v1);
          cache(v1, v2);
          cache(v2, v0);
          ++t;
     });

     t = hpindex(0);
     visit(indices, [&](auto v0, auto v1, auto v2) {
          move(v0, v1);
          move(v1, v2);
          move(v2, v0);
          ++t;
     });

     return neighbors;
}

Triplets<hpindex> seal(Triplets<hpindex> neighbors) {
     auto nTriangles = neighbors.size() / 3;
     auto m = nTriangles;
     auto t = hpindex(0);

     auto do_seal = [&](auto t, auto& i) {
          static const trit o[3] = { TRIT1, TRIT2, TRIT0 };

          neighbors.push_back(t);
          auto walker0 = make_spokes_walker(neighbors, t, i);
          while(std::get<0>(*walker0) < nTriangles) ++walker0;
          neighbors.push_back(std::get<0>(*walker0));
          auto walker1 = make_spokes_walker(neighbors, t, o[i]);
          while(std::get<0>(*walker1) < nTriangles) --walker1;
          neighbors.push_back(std::get<0>(*walker1));
     };

     for(auto& neighbor : neighbors) if(neighbor == std::numeric_limits<hpindex>::max()) neighbor = m++;
     neighbors.reserve(neighbors.size() + 3 * (m - nTriangles));

     visit(neighbors, [&](auto n0, auto n1, auto n2) {
          if(n0 >= nTriangles) do_seal(t, TRIT0);
          if(n1 >= nTriangles) do_seal(t, TRIT1);
          if(n2 >= nTriangles) do_seal(t, TRIT2);
          ++t;
     });

     return neighbors;
}

hpuint size(trm::SpokesEnumerator e) {
     auto valence = hpuint(0);

     do ++valence; while(++e);
     return valence;
}

}//namespace happah

