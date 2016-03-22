// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/geometries/Geometry.h"
#include "happah/geometries/Loops.h"

template<class Vertex>
class SegmentLoops : public Geometry1D<typename Vertex::SPACE>, public Loops<Vertex> {
     using Space = typename Vertex::SPACE;

public:
     using IndicesArrays = typename Loops<Vertex>::IndicesArrays; 
     using Vertices = typename Loops<Vertex>::Vertices;

     SegmentLoops(Vertices vertices, IndicesArrays loops)
          : Geometry1D<Space>(), Loops<Vertex>(std::move(vertices), std::move(loops)) {}

     virtual ~SegmentLoops() {}

};

