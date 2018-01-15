// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometry/TriangleMesh.hpp"

namespace happah {

Triples<trix> make_neighbors(const Triples<hpindex>& indices) {
     auto map = make_map<std::pair<trix, trix> >(0);
     auto neighbors = Triples<trix>();
     auto t = hpindex(0);

     auto cache = [&](auto va, auto vb, auto i) {
          auto key = std::make_pair(va, vb);
          auto temp = map.find(key);

          if(temp == map.end()) map[key] = std::make_pair(trix(t, i), trix());
          else temp->second.second = trix(t, i);
     };

     auto move = [&](auto va, auto vb) {
          auto value = map[{ va, vb }];

          if(value.first.getTriple() == t) neighbors.push_back(value.second);
          else neighbors.push_back(value.first);
     };

     neighbors.reserve(indices.size());

     t = hpindex(0);
     visit(indices, [&](auto v0, auto v1, auto v2) {
          cache(v0, v1, TRIT0);
          cache(v1, v2, TRIT1);
          cache(v2, v0, TRIT2);
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

Triples<trix> seal(Triples<trix> neighbors) {
     auto nTriangles = neighbors.size() / 3;
     auto m = nTriangles;
     auto t = hpindex(0);

     auto do_seal = [&](auto x) {
          neighbors.push_back(x);
          auto walker0 = make_spokes_walker(neighbors, x);
          while(std::get<0>(*walker0) < nTriangles) ++walker0;
          neighbors.emplace_back(std::get<0>(*walker0), std::get<1>(*walker0));
          auto walker1 = make_spokes_walker(neighbors, x.getNext());
          while(std::get<0>(*walker1) < nTriangles) --walker1;
          neighbors.emplace_back(std::get<0>(*walker1), std::get<1>(*walker1));
     };

     for(auto& neighbor : neighbors) if(neighbor == std::numeric_limits<hpindex>::max()) neighbor = trix(m++, TRIT0);
     neighbors.reserve(neighbors.size() + 3 * (m - nTriangles));

     visit(neighbors, [&](auto n0, auto n1, auto n2) {
          if(n0.getTriple() >= nTriangles) do_seal(trix(t, TRIT0));
          if(n1.getTriple() >= nTriangles) do_seal(trix(t, TRIT1));
          if(n2.getTriple() >= nTriangles) do_seal(trix(t, TRIT2));
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

