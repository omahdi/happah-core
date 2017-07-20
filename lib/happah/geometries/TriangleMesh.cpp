// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <unordered_map>

#include "happah/geometries/TriangleMesh.hpp"
#include "happah/utils/visitors.hpp"

namespace happah {

bool is_neighbor(const Indices& neighbors, hpuint t, hpuint u) {
     bool result;
     visit_triplet(neighbors, t, [&](hpuint n0, hpuint n1, hpuint n2) { result = (u == n0) || (u == n1) || (u == n2); });
     return result;
}

hpindex make_edge_offset(hpindex e) { return e - 3 * make_triangle_index(e); }

Indices make_fan(trm::FanEnumerator e) {
     auto fan = Indices();
     do fan.push_back(*e); while(++e);
     return fan;
}

trm::FanEnumerator make_fan_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

hpindex make_neighbor_index(const Indices& neighbors, hpuint t, hpuint i) { return neighbors[3 * t + i]; }

hpuint make_neighbor_offset(const Indices& neighbors, hpuint t, hpuint u) {
     auto n = std::begin(neighbors) + 3 * t;
     return (u == n[0]) ? 0 : (u == n[1]) ? 1 : 2;
}

Indices make_neighbors(const Indices& indices) {
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

     auto map = Map(0, getHash, isKeysEqual);
     auto neighbors = Indices();
     auto triangle = hpindex(0);

     auto cache = [&](hpuint va, hpuint vb) {
          auto key = Key(va, vb);
          auto i = map.find(key);
          if(i == map.end()) map[key] = Value(triangle, std::numeric_limits<hpindex>::max());
          else i->second.second = triangle;
     };

     auto move = [&](hpuint va, hpuint vb) {
          auto value = map[{ va, vb }];
          if(value.first == triangle) neighbors.push_back(value.second);
          else neighbors.push_back(value.first);
     };

     neighbors.reserve(indices.size());

     triangle = 0;
     visit_triplets(indices, [&](hpuint v0, hpuint v1, hpuint v2) {
          cache(v0, v1);
          cache(v1, v2);
          cache(v2, v0);
          ++triangle;
     });

     triangle = 0;
     visit_triplets(indices, [&](hpuint v0, hpuint v1, hpuint v2) {
          move(v0, v1);
          move(v1, v2);
          move(v2, v0);
          ++triangle;
     });

     return neighbors;
}

trm::RingEnumerator make_ring_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

trm::SpokesEnumerator make_spokes_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

trm::SpokesWalker make_spokes_walker(const Indices& neighbors, hpindex t, hpindex i) { return { neighbors, t, i }; }

hpindex make_triangle_index(hpindex e) { return e / 3; }

hpindex make_triangle_index(const Indices& indices, hpindex v) { return std::distance(std::begin(indices), std::find(std::begin(indices), std::end(indices), v)) / 3; }

hpuint make_valence(trm::FanEnumerator e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

hpuint make_valence(trm::RingEnumerator e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

hpuint make_valence(trm::SpokesEnumerator e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

hpindex make_vertex_offset(const Indices& indices, hpindex t, hpindex v) {
     auto i = std::begin(indices) + 3 * t;
     return (v == i[0]) ? 0 : (v == i[1]) ? 1 : 2;
}

trm::VerticesEnumerator make_vertices_enumerator(const Indices& neighbors) { return { neighbors }; }

Indices seal(Indices neighbors) {
     auto nTriangles = neighbors.size() / 3;
     auto m = nTriangles;
     auto t = hpindex(0);

     auto do_seal = [&](auto t, auto i) {
          static constexpr hpindex o[3] = { 1, 2, 0 };

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

     auto i = std::begin(neighbors) - 1;
     auto end = std::end(neighbors) - 1;

     while(i != end) {
          auto n0 = *(++i);
          auto n1 = *(++i);
          auto n2 = *(++i);

          if(n0 >= nTriangles) do_seal(t, 0);
          if(n1 >= nTriangles) do_seal(t, 1);
          if(n2 >= nTriangles) do_seal(t, 2);

          ++t;
     }

     return neighbors;
}

}//namespace happah

