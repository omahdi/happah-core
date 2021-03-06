// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <unordered_map>

#include "happah/geometry/TriangleGraph.hpp"

namespace happah {

Indices cut(const Triples<Edge>& edges) {
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
          auto& edge = edges[edges[e].getOpposite()];
          auto e0 = edges[edge.getPrevious()].getOpposite();
          auto e1 = edges[edge.getNext()].getOpposite();
          auto b0 = neighbors[e0 << 1] == std::numeric_limits<hpindex>::max();
          auto b1 = neighbors[e1 << 1] == std::numeric_limits<hpindex>::max();
          if(b0) cache[edge.getPrevious()] = true;
          else cache[e0] = false;
          if(b1) cache[edge.getNext()] = true;
          else cache[e1] = false;
          cache[e] = false;
          return hpindex(e);
     });
}

Triples<Edge> make_edges(const Triples<hpindex>& indices) {
     auto edges = Triples<Edge>();
     auto nEdges = indices.size();
     auto map = make_map<hpindex>(nEdges);

     auto push_edge = [&](auto va, auto vb, auto id) {
          auto key = std::make_pair(va, vb);
          auto i = map.find(key);

          if(i == std::end(map)) {
               map[key] = edges.size();
               edges.emplace_back(id, vb);
          } else {
               auto& edge = edges[(*i).second];

               edge.setOpposite(id);
               edges.emplace_back(id, edge.getId(), vb);
          }
     };

     edges.reserve(nEdges);

     visit(indices, [&](auto v0, auto v1, auto v2) {
          auto t = hpindex(edges.size() / 3);

          push_edge(v0, v1, trix(t, TRIT0));
          push_edge(v1, v2, trix(t, TRIT1));
          push_edge(v2, v0, trix(t, TRIT2));
     });

     assert(edges.size() == indices.size());
     //assert(std::find_if(std::begin(edges), std::end(edges), [](auto& edge) { return edge.opposite == std::numeric_limits<hpuint>::max(); }) == std::end(edges));
     //TODO
     /*auto i = std::find_if(std::begin(edges), std::end(edges), [](auto& edge) { return edge.opposite == std::numeric_limits<hpuint>::max(); });
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
               edges.emplace_back(edges[e].vertex, ++next, opposite, ++previous, 0);
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
     }*/

     return edges;
}

Triples<trix> make_neighbors(const Triples<Edge>& edges, hpuint nTriangles) {
     auto neighbors = Triples<trix>();

     neighbors.reserve(3 * nTriangles);
     for(auto e = std::begin(edges), end = e + 3 * nTriangles; e != end; ++e) neighbors.emplace_back((*e).getOpposite());

     return neighbors;
}

std::vector<Point2D> parametrize(const Indices& lengths, const std::vector<Point3D>& polyline) {
     auto points = std::vector<Point2D>();
     auto j = std::begin(polyline);
     auto* point0 = &polyline[0];

     auto do_parametrize = [&](auto& point0, auto& point1, auto n) {
          auto delta = (hpreal(1.0) / hpreal(n + 1)) * (point1 - point0);
          auto point = point0;

          points.emplace_back(point.x / point.z, point.y / point.z);
          while(n--) {
               point += delta;
               points.emplace_back(point.x / point.z, point.y / point.z);
          }
     };

     for(auto i : boost::make_iterator_range(std::begin(lengths), std::end(lengths) - 1)) {
          auto& point1 = *(++j);

          do_parametrize(*point0, point1, i);
          point0 = &point1;
     }
     do_parametrize(polyline.back(), polyline[0], lengths.back());

     return points;
}

hpuint size(trg::SpokesEnumerator e) {
     auto valence = 0u;

     do ++valence; while(++e);
     return valence;
}

}//namespace happah

