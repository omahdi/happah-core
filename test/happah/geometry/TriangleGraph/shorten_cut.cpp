// Copyright 2017
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

//11.2017 - Hedwig Amberg     - Initial commit.
//11.2017 - Pawel Herman      - Fix bug that skipped triangles with no edges on cut but with all three vertices on cut.

#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/NutChain.hpp>
#include <happah/util/visitors.hpp>

int main() {
     using namespace happah;
     
     auto chain = NutChain(2, 1, 2, 1, 0.5);
     auto mesh = make_triangle_mesh<VertexP3>(chain);
     auto graph = make_triangle_graph(mesh);
     auto& edges = graph.getEdges();
     auto n = 10;
     
     while(n--) {
          auto path = shorten_cut(trim(graph, cut(graph)), graph);
          auto indices = std::get<1>(analyze(graph, path));
          auto cache = std::vector<Indices>(graph.getNumberOfVertices());
          auto b = hpindex(0);
          auto i = hpindex(-1);

          for(auto e : path) {
               auto& temp = cache[edges[e].vertex];

               if(++i == indices[b]) {
                    temp.push_back(b);
                    if(++b == indices.size()) b = hpindex(0);
               }
               temp.push_back(b);
          }

          visit_triplets(edges, [&](auto& edge0, auto& edge1, auto& edge2) {
               auto& t0 = cache[edge0.vertex];
               auto& t1 = cache[edge1.vertex];
               auto& t2 = cache[edge2.vertex];

               for(auto x : t0) for(auto y : t1) for(auto z : t2) assert(x != y || y != z);
          });
     }
}
