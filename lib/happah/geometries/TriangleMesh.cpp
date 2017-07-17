// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <unordered_map>

#include "happah/geometries/TriangleMesh.h"

namespace happah {

namespace trg {

FanEnumerator make_fan_enumerator(const std::vector<Edge>& edges, hpuint e) { return { { edges, e } }; }

RingEnumerator make_ring_enumerator(const std::vector<Edge>& edges, hpuint e) { return { { edges, e } }; }

SpokesEnumerator make_spokes_enumerator(const std::vector<Edge>& edges, hpuint e) { return { { edges, e } }; }

SpokesWalker make_spokes_walker(const std::vector<Edge>& edges, hpindex e) { return { edges, e }; }

}//namespace trg

namespace trm {

FanEnumerator make_fan_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

RingEnumerator make_ring_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

SpokesEnumerator make_spokes_enumerator(const Indices& neighbors, hpuint t, hpuint i) { return { { neighbors, t, i } }; }

SpokesWalker make_spokes_walker(const Indices& neighbors, hpindex t, hpindex i) { return { neighbors, t, i }; }

VerticesEnumerator make_vertices_enumerator(const Indices& neighbors) { return { neighbors }; }

}//namespace trm

bool is_neighbor(const Indices& neighbors, hpuint t, hpuint u) {
     bool result;
     visit_triplet(neighbors, t, [&](hpuint n0, hpuint n1, hpuint n2) { result = (u == n0) || (u == n1) || (u == n2); });
     return result;
}

hpindex make_edge_index(const Edge& edge) { return 3 * make_triangle_index(edge) + make_edge_offset(edge); }

hpindex make_edge_offset(hpindex e) { return e - 3 * make_triangle_index(e); }

hpindex make_edge_offset(const Edge& edge) { return 3 - edge.next - edge.previous + 6 * make_triangle_index(edge); }

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

Indices make_fan(trm::FanEnumerator e) {
     auto fan = Indices();
     do fan.push_back(*e); while(++e);
     return fan;
}

Indices make_fan(trg::FanEnumerator e) {
     auto fan = Indices();
     do fan.push_back(*e); while(++e);
     return fan;
}

hpindex make_neighbor_index(const Indices& neighbors, hpuint t, hpuint i) { return neighbors[3 * t + i]; }

hpindex make_neighbor_index(const std::vector<Edge>& edges, hpuint t, hpuint i) { return make_triangle_index(edges[3 * t + i].opposite); }

hpuint make_neighbor_offset(const Indices& neighbors, hpuint t, hpuint u) {
     auto n = std::begin(neighbors) + 3 * t;
     return (u == n[0]) ? 0 : (u == n[1]) ? 1 : 2;
}

hpuint make_neighbor_offset(const std::vector<Edge>& edges, hpuint t, hpuint u) {
     auto& edge = edges[3 * t];
     if(make_triangle_index(edge.opposite) == u) return 0;
     else if(make_triangle_index(edges[edge.next].opposite) == u) return 1;
     else return 2;
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

     std::vector<hpuint> neighbors;
     neighbors.reserve(indices.size());

     Map map(0, getHash, isKeysEqual);
     hpuint triangle;

     auto cache = [&](hpuint va, hpuint vb) {
          Key key(va, vb);
          auto i = map.find(key);
          if(i == map.end()) map[key] = Value(triangle, std::numeric_limits<hpuint>::max());
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

hpindex make_triangle_index(hpindex e) { return e / 3; }

hpindex make_triangle_index(const Indices& indices, hpindex v) { return std::distance(std::begin(indices), std::find(std::begin(indices), std::end(indices), v)) / 3; }

hpindex make_triangle_index(const Edge& edge) { return make_triangle_index(edge.next); }

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

hpindex make_vertex_offset(const Indices& indices, hpindex t, hpindex v) {
     auto i = std::begin(indices) + 3 * t;
     return (v == i[0]) ? 0 : (v == i[1]) ? 1 : 2;
}

}//namespace happah

