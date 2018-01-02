// Copyright 2017
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

//11.2017 - Hedwig Amberg     - Initial commit.
//11.2017 - Pawel Herman      - Fix bug that skipped triangles with no edges on cut but with all three vertices on cut.

#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/NutChain.hpp>

int main() {
     using namespace happah;
     
     auto chain = NutChain(2, 1, 2, 1, 0.5);
     auto mesh = make_triangle_mesh<VertexP3>(chain);
     auto graph = make_triangle_graph(mesh);
     auto& edges = graph.getEdges();
     auto n = 10;
     
     while(n--) {
          auto path = undegenerate(graph, trim(graph, cut(graph)));
          auto indices = std::get<1>(analyze(graph, path));
          auto cache0 = Triplets<hpindex>(3 * size(graph), std::numeric_limits<hpindex>::max());
          auto cache1 = Triplets<hpindex>(3 * size(graph), std::numeric_limits<hpindex>::max());
          auto b = hpindex(0);
          auto i = hpindex(-1);

          auto lambda = [&](auto i0, auto i1, auto i2) {
               if(i0 == std::numeric_limits<hpindex>::max()) return;
               if(i1 == std::numeric_limits<hpindex>::max()) return;
               if(i2 == std::numeric_limits<hpindex>::max()) return;
               assert(i0 != i1 || i1 != i2);
          };

          for(auto e : path) {
               if(++i == indices[b]) {
                    auto walker0 = make_spokes_walker(edges, edges[e].opposite);
                    do cache0[*(--walker0)] = b; while(std::find(std::begin(path), std::end(path), *walker0) == std::end(path));
                    if(++b == indices.size()) b = hpindex(0);
                    auto walker1 = make_spokes_walker(edges, edges[e].opposite);
                    do cache1[*(--walker1)] = b; while(std::find(std::begin(path), std::end(path), *walker1) == std::end(path));
               } else {
                    auto walker = make_spokes_walker(edges, edges[e].opposite);
                    do {
                         --walker;
                         cache0[*walker] = b;
                         cache1[*walker] = b;
                    } while(std::find(std::begin(path), std::end(path), *walker) == std::end(path));
               }
          }

          visit(cache0, lambda);
          visit(cache1, lambda);
     }
}
