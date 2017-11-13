// Copyright 2017
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

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

          std::cout << "building cache" << std::endl;

          auto b = hpindex(0);
          for(auto i = std::begin(indices) + 1, end = std::end(indices); i != end; ++i) {
               for(auto e : boost::make_iterator_range(std::begin(path) + *(i - 1), std::begin(path) + *i)) cache[edges[e].vertex].push_back(b);
               cache[edges[path[*i]].vertex].push_back(b);
               ++b;
          }
          for(auto e : boost::make_iterator_range(std::begin(path) + indices.back(), std::end(path))) cache[edges[e].vertex].push_back(b);
          for(auto e : boost::make_iterator_range(std::begin(path), std::begin(path) + indices[0])) cache[edges[e].vertex].push_back(b);
          cache[edges[path[indices[0]]].vertex].push_back(b);

          std::cout << "testing triangles" << std::endl;

          visit_triplets(edges, [&](auto& edge0, auto& edge1, auto& edge2) {
               auto& t0 = cache[edge0.vertex];
               auto& t1 = cache[edge1.vertex];
               auto& t2 = cache[edge2.vertex];

               for(auto x : t0) for(auto y : t1) for(auto z : t2) assert(x != y || y != z);
          });
     }
}
