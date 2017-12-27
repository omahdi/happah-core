// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <unordered_map>

#include "happah/geometry/QuadMesh.hpp"
#include "happah/util/visitors.hpp"

namespace happah {

Indices make_neighbors(quads, const Indices& indices) {
     
     using Key = std::pair<hpuint, hpuint>;
     using Value = std::pair<hpuint, hpuint>;

     auto getHash = [](const Key& k) -> uint64_t {
          int32_t d = k.first - k.second;
          int32_t min = k.second + (d & d >> 31);
          int32_t max = k.first - (d & d >> 31);
          return ((uint64_t)max << 32 | min);
     };
     
     auto isKeysEqual = [](const Key& k1, const Key& k2) {
          return (k1.first == k2.first && k1.second == k2.second) || (k1.first == k2.second && k1.second == k2.first);
     };

     using Map = std::unordered_map<Key, Value, decltype(getHash), decltype(isKeysEqual)>;

     auto map = Map(0, getHash, isKeysEqual);
     auto neighbors = Indices();
     neighbors.reserve(indices.size());
     auto q = hpindex(0);
     
     auto cache = [&](hpuint va, hpuint vb) {
          auto key = Key(va, vb);
          auto i = map.find(key);
          if(i == map.end()) map[key] = Value(q, std::numeric_limits<hpindex>::max());
          else i->second.second = q;
     };
     
     auto move = [&](hpuint va, hpuint vb) {
          auto value = map[{ va, vb }];
          if(value.first == q) neighbors.push_back(value.second);
          else neighbors.push_back(value.first);
     };

     visit_quartets(indices, [&](hpuint v0, hpuint v1, hpuint v2, hpuint v3) {
          cache(v0, v1);
          cache(v1, v2);
          cache(v2, v3);
          cache(v3, v0);
          ++q;
     });

     q = hpindex(0);

     visit_quartets(indices, [&](hpuint v0, hpuint v1, hpuint v2, hpuint v3) {
          move(v0, v1);
          move(v1, v2);
          move(v2, v3);
          move(v3, v0);
          ++q;
     });

     return neighbors;
}

}//namespace happah
