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
          
          auto analyzed = analyze(graph, path);
          auto valences = std::get<0>(analyzed);
          auto indices = std::get<1>(analyzed);
          auto pairings = std::get<2>(analyzed);
          
          auto v = Indices();
          auto i = Indices();
          for(hpuint j = 0; j < size(path); ++j) {
               auto e = path[j];
               auto n = hpuint(0);
               visit(make_spokes_enumerator(edges, edges[e].opposite), [&](auto e) { if(std::find(std::begin(path), std::end(path), e) != std::end(path)) ++n; });
               if(n > hpuint(2)){
                    i.push_back(j);
                    v.push_back(n);
               }
          }
          
          assert(v == valences);
          assert(i == indices);
          
          auto branch = Indices(graph.getNumberOfEdges(), std::numeric_limits<hpindex>::max());
          auto b = size(indices) - 1;
          auto next = 0;
          for(auto e : path){
               branch[e] = b;
               if(e == path[indices[next]]){
                    b = next;
                    next = (next+1) % size(indices);
               }
          }
          for(auto e : path){
               assert(branch[edges[e].opposite] == pairings[branch[e]]);
          }
          
     }
}