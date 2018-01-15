// Copyright 2017
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

//11.2017 - Hedwig Amberg     - Initial commit.

#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/NutChain.hpp>

int main() {
     using namespace happah;
     
     auto chain = NutChain(2, 1, 2, 1, 0.5);
     auto mesh = make_triangle_mesh<VertexP3>(chain);
     auto graph = make_triangle_graph(mesh);
     auto& edges = graph.getEdges();
     auto n = 100;
     
     while(n--) {
          auto path = trim(graph, cut(graph));
          auto analysis = analyze(graph, path);
          auto& valences = std::get<0>(analysis);
          auto& indices = std::get<1>(analysis);
          auto& pairings = std::get<2>(analysis);
          auto cache = Indices(graph.getNumberOfEdges(), std::numeric_limits<hpindex>::max());
          auto b = size(indices) - 1;
          auto f = path[indices[0]];
          auto v = std::begin(valences);
          auto i = std::begin(indices);
          auto j = hpindex(-1);

          assert(size(valences) == size(indices));
          for(auto e : path) {
               auto n = hpuint(0);

               visit(make_spokes_enumerator(edges, edges[e].getOpposite()), [&](auto e) { if(std::find(std::begin(path), std::end(path), e) != std::end(path)) ++n; });
               if(i != std::end(indices) && ++j == *i) {
                    assert(n == *v);
                    ++v;
                    ++i;
               } else assert(n == 2);
               cache[e] = b;
               if(e == f) {
                    if(++b == size(indices)) b = hpindex(0);
                    f = (b + 1 == size(indices)) ? std::numeric_limits<hpindex>::max() : path[indices[b + 1]];
               }
          }
          
          for(auto e : path) assert(cache[edges[e].getOpposite()] == pairings[cache[e]]);
     }
}
