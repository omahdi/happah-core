// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/math/Space.hpp"

/**
 * @section DESCRIPTION
 *
 * A vertex is a zero-dimensional geometry with attributes including position, 
 * normal, color, or texture coordinates.  A vertex always has a position.
 */

namespace happah {

template<class Space>
class VertexP {
     using Point = typename Space::POINT;

public:
     using SPACE = Space;

     Point position{ 0.0 };

     VertexP() = default;
     VertexP(Point position) : position(std::move(position)) {}

};//VertexP

template<class Space>
class VertexPC {
     using Point = typename Space::POINT;

public:
     using SPACE = Space;

     Point position{ 0.0 };
     hpcolor color{ 0.0, 0.0, 0.0, 1.0 };

     VertexPC() = default;
     VertexPC(Point position) : position(std::move(position)) {}
     VertexPC(Point position, hpcolor color) : position(std::move(position)), color(std::move(color)) {}

};//VertexPC

template<class Space>
class VertexPN {
     using Point = typename Space::POINT;
     using Vector = typename Space::VECTOR;

public:
     using SPACE = Space;

     Point position{ 0.0 };
     Vector normal{ 0.0 };

     VertexPN() = default;
     VertexPN(Point position) : position(std::move(position)) {}
     VertexPN(Point position, Vector normal) : position(std::move(position)), normal(std::move(normal)) {}

};//VertexPN

template<class Space>
class VertexPNC {
     using Point = typename Space::POINT;
     using Vector = typename Space::VECTOR;

public:
     using SPACE = Space;

     Point position{ 0.0 };
     Vector normal{ 0.0 };
     hpcolor color{ 0.0, 0.0, 0.0, 1.0 };

     VertexPNC() = default;
     VertexPNC(Point position) : position(std::move(position)) {}
     VertexPNC(Point position, Vector normal) : position(std::move(position)), normal(std::move(normal)) {}
     VertexPNC(Point position, hpcolor color) : position(std::move(position)), color(std::move(color)) {}
     VertexPNC(Point position, Vector normal, hpcolor color) : position(std::move(position)), normal(std::move(normal)), color(std::move(color)) {}

};//VertexPNC

using VertexP2 = VertexP<Space2D>;
using VertexP3 = VertexP<Space3D>;
using VertexP4 = VertexP<Space4D>;
using VertexP3C = VertexPC<Space3D>;
using VertexP3N = VertexPN<Space3D>;

template<class Vertex, typename = void>
struct contains_position : public std::false_type {};

template<class Vertex>
struct contains_position<Vertex, typename std::enable_if<std::is_base_of<VertexP<typename Vertex::SPACE>, Vertex>::value || std::is_base_of<VertexPN<typename Vertex::SPACE>, Vertex>::value>::type> : public std::true_type {};

template<class Vertex>
struct contains_color : public std::integral_constant<bool, std::is_base_of<VertexPC<typename Vertex::SPACE>, Vertex>::value || std::is_base_of<VertexPNC<typename Vertex::SPACE>, Vertex>::value> {};

template<class Vertex>
struct contains_normal : public std::integral_constant<bool, std::is_base_of<VertexPN<typename Vertex::SPACE>, Vertex>::value> {};

namespace detail {

template<class Vertex>
struct make_vertex;

template<>
struct make_vertex<VertexP2> {
     template<class Iterator>
     static VertexP2 call(Iterator begin) { return VertexP2(Point2D(begin[0], begin[1])); }
};

template<>
struct make_vertex<VertexP3> {
     template<class Iterator>
     static VertexP3 call(Iterator begin) { return VertexP3(Point3D(begin[0], begin[1], begin[2])); }

};

template<>
struct make_vertex<VertexP3N> {
     template<class Iterator>
     static VertexP3N call(Iterator begin) { return VertexP3N(Point3D(begin[0], begin[1], begin[2]), Vector3D(begin[3], begin[4], begin[5])); }

};

template<>
struct make_vertex<VertexP3C> {
     template<class Iterator>
     static VertexP3C call(Iterator begin) { return VertexP3C(Point3D(begin[0], begin[1], begin[2]), hpcolor(begin[3], begin[4], begin[5], begin[6])); }

};

}//namespace detail

template<class Vertex, class Iterator>
Vertex make_vertex(Iterator begin) { return detail::make_vertex<Vertex>::call(begin); }

}//namespace happah

