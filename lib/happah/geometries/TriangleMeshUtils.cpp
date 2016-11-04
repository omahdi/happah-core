// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/TriangleMeshUtils.h"

#include <unordered_map>
#include <utility>

#include "happah/utils/visitors.h"

namespace happah {

hpuint TriangleMeshUtils::findTriangle(const std::vector<hpuint>& indices, hpuint v0, hpuint v1, hpuint& v2) {
     hpuint triangle = 0;
     for(auto i = indices.begin(), end = indices.end(); i != end; ++i) {
          hpuint i0 = *i;
          hpuint i1 = *(++i);
          hpuint i2 = *(++i);
          if(i0 == v0) {
               if(i1 == v1) {
                    v2 = i2;
                    return triangle;
               } else if(i2 == v1) {
                    v2 = i1;
                    return triangle;
               }
          } else if(i0 == v1) {
               if(i1 == v0) {
                    v2 = i2;
                    return triangle;
               } else if(i2 == v0) {
                    v2 = i1;
                    return triangle;
               }
          } else if((i1 == v0 && i2 == v1) || (i1 == v1 && i2 == v0)) {
               v2 = i0;
               return triangle;
          }
          ++triangle;
     }
     return happah::UNULL;
}

std::vector<hpuint> TriangleMeshUtils::getNeighbors(const std::vector<hpuint>& indices) {
     using Key = std::pair<hpuint, hpuint>;
     using Value = std::pair<hpuint, hpuint>;

     auto getHash = [](const Key& k) -> uint64_t {
          int32_t d = k.first - k.second;
          int32_t min = k.second + (d & d >> 31);
          int32_t max = k.first - (d & d >> 31);
          return ((uint64_t)max << 32 | min);
     };

     auto isKeysEqual = [](const Key& k1, const Key& k2) { return (k1.first == k2.first && k1.second == k2.second) || (k1.first == k2.second && k1.second == k2.first); };

     using Map = std::unordered_map<Key, Value, decltype(getHash), decltype(isKeysEqual)>;

     std::vector<hpuint> neighbors;
     neighbors.reserve(indices.size());

     Map map(0, getHash, isKeysEqual);
     hpuint triangle;

     auto cache = [&](hpuint va, hpuint vb) {
          Key key(va, vb);
          auto i = map.find(key);
          if(i == map.end()) map[key] = Value(triangle, UNULL);
          else i->second.second = triangle;
     };

     triangle = 0;
     visit_triplets(indices, [&](hpuint v0, hpuint v1, hpuint v2) {
          cache(v0, v1);
          cache(v0, v2);
          cache(v1, v2);
          ++triangle;
     });

     auto move = [&](hpuint va, hpuint vb) {
          auto value = map[{ va, vb }];
          if(value.first == triangle) neighbors.push_back(value.second);
          else neighbors.push_back(value.first);
     };

     triangle = 0;
     visit_triplets(indices, [&](hpuint v0, hpuint v1, hpuint v2) {
          move(v0, v1);
          move(v1, v2);
          move(v2, v0);
          ++triangle;
     });

     return neighbors;
}

std::tuple<hpuint, hpuint, hpuint> TriangleMeshUtils::getNeighbors(const std::vector<hpuint>& indices, hpuint triangle) {
     hpuint n1 = happah::UNULL;
     hpuint n2 = happah::UNULL;
     hpuint n3 = happah::UNULL;
     hpuint temp = 3 * triangle;
     std::vector<hpuint>::const_iterator v = indices.cbegin() + temp;
     hpuint v1 = *v;
     hpuint v2 = *(++v);
     hpuint v3 = *(++v);

     hpuint t = -1;
     for(std::vector<hpuint>::const_iterator i = indices.cbegin(); i != indices.cend(); ++i) {
          ++t;
          if(t == triangle) {
               ++(++i);
               continue;
          }
          hpuint w1 = *i;
          hpuint w2 = *(++i);
          hpuint w3 = *(++i);
          if((v1 == w1 && v2 == w2) || (v1 == w2 && v2 == w1) || (v1 == w2 && v2 == w3) || (v1 == w3 && v2 == w2) || (v1 == w1 && v2 == w3) || (v1 == w3 && v2 == w1)) n1 = t;
          if((v2 == w1 && v3 == w2) || (v2 == w2 && v3 == w1) || (v2 == w2 && v3 == w3) || (v2 == w3 && v3 == w2) || (v2 == w1 && v3 == w3) || (v2 == w3 && v3 == w1)) n2 = t;
          if((v3 == w1 && v1 == w2) || (v3 == w2 && v1 == w1) || (v3 == w2 && v1 == w3) || (v3 == w3 && v1 == w2) || (v3 == w1 && v1 == w3) || (v3 == w3 && v1 == w1)) n3 = t;

          if(n1 != happah::UNULL && n2 != happah::UNULL && n3 != happah::UNULL) break;
     }

     return std::make_tuple(n1, n2, n3);
}

}//namespace happah

