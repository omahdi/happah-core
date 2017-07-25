// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <memory>
#include <vector>

#include "happah/geometries/Model.hpp"

namespace happah {

template<class Vertex>
class VertexCloud : public Model<Vertex> {
public:
     VertexCloud(std::vector<Vertex> vertices) 
          : Model<Vertex>(std::move(vertices)) {}

     virtual ~VertexCloud() {}

};

typedef VertexCloud<VertexP2> PointCloud2D;
typedef VertexCloud<VertexP3> PointCloud3D;

}//namespace happah

