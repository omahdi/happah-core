// Copyright 2017
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

//12.2017 - Hedwig Amberg     - Initial commit.

#include <boost/range/irange.hpp>

#include <happah/geometry/QuadMesh.hpp>
#include <happah/geometry/RectangularCuboid.hpp>

int main() {
     using namespace happah;
     
     auto cuboid = RectangularCuboid(1, 2, 3);
     auto mesh = make_quad_mesh<VertexP3>(cuboid);
     auto neighbors = make_neighbors(mesh);
     auto solutions = std::vector<quax>({
          quax(0, QUAT3), quax(2, QUAT0), quax(3, QUAT1),
          quax(0, QUAT2), quax(1, QUAT0), quax(2, QUAT1),
          quax(0, QUAT1), quax(4, QUAT0), quax(1, QUAT1),
          quax(0, QUAT0), quax(3, QUAT0), quax(4, QUAT1),
          quax(2, QUAT3), quax(5, QUAT0), quax(3, QUAT2),
          quax(1, QUAT3), quax(5, QUAT1), quax(2, QUAT2),
          quax(1, QUAT2), quax(4, QUAT3), quax(5, QUAT2),
          quax(3, QUAT3), quax(5, QUAT3), quax(4, QUAT2)
     });
     auto s = std::begin(solutions) - 1;

     for(auto v : boost::irange(hpindex(0), mesh.getNumberOfVertices())) visit(make_spokes_enumerator(mesh, neighbors, v), [&](auto x) { assert(*(++s) == x); });
}

