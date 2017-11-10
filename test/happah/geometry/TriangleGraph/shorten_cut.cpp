// Copyright 2017
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/NutChain.hpp>

int main() {
     using namespace happah;
     
     auto problematic = [&](Indices& cut, const TriangleGraph<VertexP3>& graph){
          auto& edges = graph.getEdges();
          auto analyzed = analyze(graph, cut);
          Indices& indices = std::get<1>(analyzed);

          /*collect all branches a vertex is contained in*/
          auto branches = std::map<hpuint, std::vector<hpuint> >();
          hpuint branch = 0;
          for(hpuint i = 0; i < size(cut); ++i){
               hpuint v = edges[cut[i]].vertex;
               if(i == indices[branch]){
                    auto b = branches[v];
                    b.push_back(branch);
                    branches[v] = b;
                    branch = (branch + 1) % size(indices);
               }
               auto b = branches[v];
               b.push_back(branch);
               branches[v] = b;
          }

          /*check if a full triangle is contained in one branch*/
          for(auto e : cut){
               hpuint v0 = edges[e].vertex;
               hpuint v1 = edges[edges[e].next].vertex;
               hpuint v2 = edges[edges[edges[e].next].next].vertex;
               for(auto b0 : branches[v0]){
                    for(auto b1 : branches[v1]){
                         for(auto b2 : branches[v2]){
                              if( (b0 == b1) && (b1 == b2) ){
                                   return true;
                              }
                         }
                    }
               }
          }
          return false;
     };
     
     auto chain = NutChain(2, 1, 2, 1, 0.5);
     auto mesh = make_triangle_mesh<VertexP3>(chain);
     const auto graph = make_triangle_graph(mesh);
     
     hpuint tests = 10;
     
     while(tests > 0){
          
          auto cut_path = trim(graph, cut(graph));
          auto length = size(cut_path);

          if(problematic(cut_path, graph)){
               
               --tests;
               auto new_path = shorten_cut(cut_path, graph);

               assert(size(new_path) < length);
               assert(!problematic(new_path, graph));
          }
     }
}
