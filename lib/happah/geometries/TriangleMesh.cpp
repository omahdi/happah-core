// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/TriangleMesh.h"

namespace happah {

boost::optional<hpuint> find_in_ring(const std::vector<Edge>& edges, hpuint begin, hpuint v) { return find_if_in_spokes(edges, begin, [&](const Edge& edge) { return edge.vertex == v; }); }

bool is_neighbor(const Indices& neighbors, hpuint t, hpuint n) {
     bool result;
     visit_triplet(neighbors, t, [&](hpuint n0, hpuint n1, hpuint n2) { result = (n == n0) || (n == n1) || (n == n2); });
     return result;
}

std::vector<Edge> make_edges(const std::vector<hpuint>& indices) {
     std::vector<Edge> edges;

     auto nEdges = indices.size();//NOTE: The number of edges is >= to 3x the number of triangles; the number is greater if the mesh is not closed, that is, it has a border.
     edges.reserve(nEdges);

     using Key = std::pair<hpuint, hpuint>;
     using Value = hpuint;

     auto getHash = [](const Key& k) -> uint64_t {
          int32_t d = k.first - k.second;
          int32_t min = k.second + (d & d >> 31);
          int32_t max = k.first - (d & d >> 31);
          return ((uint64_t)max << 32 | min);
     };

     auto isKeysEqual = [](const Key& k1, const Key& k2) { return (k1.first == k2.first && k1.second == k2.second) || (k1.first == k2.second && k1.second == k2.first); };

     using Map = std::unordered_map<Key, Value, decltype(getHash), decltype(isKeysEqual)>;
     Map map(nEdges, getHash, isKeysEqual);

     auto e = 0u;

     auto push_edge = [&](hpuint va, hpuint vb, hpuint next, hpuint previous) {
          auto key = Key(va, vb);
          auto i = map.find(key);
          if(i == map.end()) {
               map[key] = e;
               edges.emplace_back(vb, next, UNULL, previous);
          } else {
               auto opposite = (*i).second;
               edges[opposite].opposite = e;
               edges.emplace_back(vb, next, opposite, previous);
          }
     };

     for(auto i = std::begin(indices), end = std::end(indices); i != end; ++i) {
          auto v0 = *i;
          auto v1 = *(++i);
          auto v2 = *(++i);

          auto e0 = e;
          auto e1 = e + 1;
          auto e2 = e + 2;

          push_edge(v0, v1, e1, e2);
          ++e;
          push_edge(v1, v2, e2, e0);
          ++e;
          push_edge(v2, v0, e0, e1);
          ++e;
     }

     assert(edges.size() == indices.size());

     auto i = std::find_if(std::begin(edges), std::end(edges), [](const Edge& edge) { return edge.opposite == UNULL; });
     if(i != edges.end()) {
          e = std::distance(edges.begin(), i);
          auto begin = e;
          auto next = indices.size();
          auto previous = next - 2;
          do {
               (*i).opposite = next;
               i = edges.begin() + (*i).previous;
               edges.emplace_back((*i).vertex, ++next, e, ++previous);
               e = edges[(*i).opposite].previous;
               i = edges.begin() + e;
          } while(e != begin);
          edges[indices.size()].previous = edges.size() - 1;
          edges.back().next = indices.size();
     }

     return edges;
}//make_edges


hpuint make_neighbor_offset(const Indices& neighbors, hpuint p, hpuint q) {
     auto n = neighbors.begin() + 3 * p;
     return (q == *n) ? 0 : (q == *(n + 1)) ? 1 : 2;
}

std::vector<hpuint> make_neighbors(const Indices& indices) {
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
          cache(v1, v2);
          cache(v2, v0);
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

}//namespace happah

