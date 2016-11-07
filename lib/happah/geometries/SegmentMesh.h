// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <memory>

#include "happah/geometries/Geometry.h"
#include "happah/geometries/Mesh.h"
#include "happah/math/Space.h"

namespace happah {

//NOTE: The technically most accurate name of this class would be LineSegmentMesh but this is shortened to SegmentMesh.
template<class Vertex>
class SegmentMesh : public Geometry1D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;

public:
     using Vertices = typename Mesh<Vertex>::Vertices;

     SegmentMesh(Vertices vertices, Indices indices)
          : Geometry1D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)) {}

};
typedef SegmentMesh<VertexP<Space2D> > SegmentMesh2D;
typedef std::shared_ptr<SegmentMesh2D> SegmentMesh2D_ptr;
typedef SegmentMesh<VertexPN<Space3D> > SegmentMesh3D;
typedef std::shared_ptr<SegmentMesh3D> SegmentMesh3D_ptr;

template<class M, class Space = typename M::SPACE, class Vertex = typename M::VERTEX>
struct is_segment_mesh {
     static const bool value = std::is_base_of<SegmentMesh<Vertex>, M>::value && std::is_base_of<typename M::SPACE, Space>::value;
};

}//namespace happah

