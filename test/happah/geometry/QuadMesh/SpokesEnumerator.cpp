// Copyright 2017
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

//12.2017 - Hedwig Amberg     - Initial commit.

#include <happah/geometry/QuadMesh.hpp>
#include <happah/geometry/RectangularCuboid.hpp>
#include <happah/util/visitors.hpp>

int main() {
     using namespace happah;
     
     auto rec = RectangularCuboid(1, 2, 3);
     auto mesh = make_quad_mesh<VertexP3>(rec);
     auto neighbors = make_neighbors(mesh);
     
     std::vector<std::vector<std::pair<hpuint, hpuint> >> solutions = {
          { {0,3}, {2,0}, {3,1} },
          { {0,2}, {1,0}, {2,1} },
          { {0,1}, {4,0}, {1,1} },
          { {0,0}, {3,0}, {4,1} },
          { {2,3}, {5,0}, {3,2} },
          { {1,3}, {5,1}, {2,2} },
          { {1,2}, {4,3}, {5,2} },
          { {3,3}, {5,3}, {4,2} } };

     hpuint max = 8;
     for(hpuint v = 0; v < max; ++v) {
          auto spokes_enum = make_spokes_enumerator(mesh, neighbors, hpindex(v));

          hpuint j = 0;
          do {
               assert( solutions[v][j].first == std::get<0>(*spokes_enum) && solutions[v][j].second == std::get<1>(*spokes_enum) );
               ++j;
          } while(++spokes_enum);
     }
}
