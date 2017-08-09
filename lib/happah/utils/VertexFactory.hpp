// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/geometries/Vertex.hpp"

namespace happah {

template<class Vertex>
class VertexFactory;

template<class Space>
class VertexFactory<VertexP<Space> > {
     using Point = typename Space::POINT;
     using Vector = typename Space::VECTOR;
     using Vertex = VertexP<Space>;

public:
     Vertex operator()(Point position) const { return {std::move(position)}; }

     Vertex operator()(Point position, Vector normal) const { return {std::move(position)}; }

};

template<class Space>
class VertexFactory<VertexPC<Space> > {
     using Point = typename Space::POINT;
     using Vector = typename Space::VECTOR;
     using Vertex = VertexPC<Space>;

public:
     Vertex operator()(Point position) const { return {std::move(position)}; }

     Vertex operator()(Point position, hpcolor color) const { return {std::move(position), std::move(color)}; }

};

template<class Space>
class VertexFactory<VertexPN<Space> > {
     using Point = typename Space::POINT;
     using Vector = typename Space::VECTOR;
     using Vertex = VertexPN<Space>;

public:
     Vertex operator()(Point position) const { return {std::move(position), Vector(0.0)}; }

     Vertex operator()(Point position, Vector normal) const { return {std::move(position), std::move(normal)}; }

};

}//namespace happah

