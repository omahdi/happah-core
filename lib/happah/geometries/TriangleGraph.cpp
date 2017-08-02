// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <unordered_map>

#include "happah/geometries/TriangleGraph.hpp"
#include "happah/utils/visitors.hpp"

namespace happah {

Indices cut(const std::vector<Edge>& edges) {
     auto cache = boost::dynamic_bitset<>(edges.size(), false);
     auto range = std::mt19937();

     cache[0] = true;
     cache[1] = true;
     cache[2] = true;
     range.seed(std::random_device()());

     return cut(edges, 0, [&](auto& neighbors) {
          //for(auto e : boost::irange(0u, hpindex(mesh.getEdges().size())))
          //     if(neighbors[e << 1] != std::numeric_limits<hpindex>::max() && neighbors[mesh.getEdge(e).opposite << 1] == std::numeric_limits<hpindex>::max()) return e;
          if(cache.count() == 0u) return std::numeric_limits<hpindex>::max();
          auto distribution = std::uniform_int_distribution<std::mt19937::result_type>(1, cache.count());
          auto n = distribution(range);
          auto e = cache.find_first();
          while(--n) e = cache.find_next(e);
          auto& edge = edges[edges[e].opposite];
          auto e0 = edges[edge.previous].opposite;
          auto e1 = edges[edge.next].opposite;
          auto b0 = neighbors[e0 << 1] == std::numeric_limits<hpindex>::max();
          auto b1 = neighbors[e1 << 1] == std::numeric_limits<hpindex>::max();
          if(b0) cache[edge.previous] = true;
          else cache[e0] = false;
          if(b1) cache[edge.next] = true;
          else cache[e1] = false;
          cache[e] = false;
          return hpindex(e);
     });
}

std::vector<Edge> make_edges(const Indices& indices) {
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

     auto push_edge = [&](hpuint va, hpuint vb, hpuint next, hpuint previous) {
          auto key = Key(va, vb);
          auto i = map.find(key);
          if(i == map.end()) {
               map[key] = edges.size();
               edges.emplace_back(vb, next, std::numeric_limits<hpuint>::max(), previous);
          } else {
               auto opposite = (*i).second;
               edges[opposite].opposite = edges.size();
               edges.emplace_back(vb, next, opposite, previous);
          }
     };

     visit_triplets(indices, [&](auto v0, auto v1, auto v2) {
          auto e0 = edges.size();
          auto e1 = e0 + 1;
          auto e2 = e0 + 2;

          push_edge(v0, v1, e1, e2);
          push_edge(v1, v2, e2, e0);
          push_edge(v2, v0, e0, e1);
     });

     assert(edges.size() == indices.size());

     auto i = std::find_if(std::begin(edges), std::end(edges), [](auto& edge) { return edge.opposite == std::numeric_limits<hpuint>::max(); });
     while(i != (std::begin(edges) + indices.size())) {
          auto e = std::distance(std::begin(edges), i);
          auto f = edges.size();
          auto begin = e;
          auto next = f;
          auto previous = next - 2;
          while(true) {
               auto opposite = e;
               edges[e].opposite = next;
               e = edges[e].previous;
               edges.emplace_back(edges[e].vertex, ++next, opposite, ++previous);
               if(e == begin) break;
               while(edges[e].opposite != std::numeric_limits<hpuint>::max()) {
                    e = edges[edges[e].opposite].previous;
                    if(e == begin) goto exit;
               }
          }
          exit:
          edges[f].previous = edges.size() - 1;
          edges.back().next = f;
          i = std::find_if(std::begin(edges), std::begin(edges) + indices.size(), [](auto& edge) { return edge.opposite == std::numeric_limits<hpuint>::max(); });
     }

     return edges;
}//make_edges

Indices make_fan(trg::FanEnumerator e) {
     auto fan = Indices();
     do fan.push_back(*e); while(++e);
     return fan;
}

hpuint make_neighbor_offset(const std::vector<Edge>& edges, hpuint t, hpuint u) {
     auto& edge = edges[3 * t];
     if(make_triangle_index(edge.opposite) == u) return 0;
     else if(make_triangle_index(edges[edge.next].opposite) == u) return 1;
     else return 2;
}

Indices make_neighbors(const std::vector<Edge>& edges, hpuint nTriangles) {
     auto neighbors = Indices();
     neighbors.reserve(3 * nTriangles);

     for(auto e = std::begin(edges), end = e + 3 * nTriangles; e != end; ++e) {
          auto n = (*e).opposite / 3;
          if(n >= nTriangles) neighbors.push_back(std::numeric_limits<hpuint>::max());
          else neighbors.push_back(n);
     }

     return neighbors;
}

hpuint make_valence(trg::FanEnumerator e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

hpuint make_valence(trg::RingEnumerator e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

hpuint make_valence(trg::SpokesEnumerator e) {
     auto valence = 0u;
     do ++valence; while(++e);
     return valence;
}

}//namespace happah

