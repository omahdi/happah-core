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
 * normal, color, or texture coordinates.  A vertex always has a position or
 * an abscissa and an ordinate.
 */

namespace happah {

template<class Space>
class Vertex {
     static_assert(is_space<Space>::value, "A vertex can only be parameterized by a space.");

public:
     using SPACE = Space;

protected:
     Vertex() {}

};

template<class Space>
class VertexP : public Vertex<Space> {
     typedef typename Space::POINT Point;

public:
     Point position;     

     VertexP() : position(0.0) {}
     VertexP(Point position) : position(std::move(position)) {}
     
};
typedef VertexP<Space2D> VertexP2;
typedef VertexP<Space3D> VertexP3;
typedef VertexP<Space4D> VertexP4;

template<class Space>
class VertexPC : public VertexP<Space> {
     typedef typename Space::POINT Point;

public:
     hpcolor color {0.0,0.0,0.0,1.0};
     
     VertexPC() : VertexP<Space>() {}
     VertexPC(Point position) : VertexP<Space>(std::move(position)) {}
     VertexPC(Point position, hpcolor color) : VertexP<Space>(std::move(position)), color(std::move(color)) {}

};
typedef VertexPC<Space3D> VertexP3C;

template<class Space>
class VertexPN : public Vertex<Space> {
     typedef typename Space::POINT Point;
     typedef typename Space::VECTOR Vector;

public:   
     Point position;     
     Vector normal;

     VertexPN() : position(0.0), normal(0.0) {}
     VertexPN(Point position) : position(std::move(position)), normal(0.0) {}
     VertexPN(Point position, Vector normal) : position(std::move(position)), normal(std::move(normal)) {}

};
typedef VertexPN<Space3D> VertexP3N;

template<class Space>
class VertexPNC : public VertexPN<Space> {
     typedef typename Space::POINT Point;
     typedef typename Space::VECTOR Vector;

public:
     hpcolor color;

     VertexPNC() : color(0.0,0.0,0.0,1.0) {}
     VertexPNC(Point position, Vector normal, hpcolor color) : VertexPN<Space>(std::move(position), std::move(normal)), color(std::move(color)) {}

};

template<class Vertex, typename = void>
struct contains_position : public std::false_type {};

template<class Vertex>
struct contains_position<Vertex, typename std::enable_if<std::is_base_of<VertexP<typename Vertex::SPACE>, Vertex>::value || std::is_base_of<VertexPN<typename Vertex::SPACE>, Vertex>::value>::type> : public std::true_type {};

template<class V>
struct contains_color : public std::integral_constant<bool, std::is_base_of<VertexPC<typename V::SPACE>, V>::value || std::is_base_of<VertexPNC<typename V::SPACE>, V>::value> {};

template<class Vertex>
struct contains_normal : public std::integral_constant<bool, std::is_base_of<VertexPN<typename Vertex::SPACE>, Vertex>::value> {};

template<class Vertex>
struct do_make_vertex;

template<>
struct do_make_vertex<VertexP3> {
     template<class Iterator>
     static VertexP3 call(Iterator begin) { return VertexP3(Point3D(begin[0], begin[1], begin[2])); }

};

template<>
struct do_make_vertex<VertexP3N> {
     template<class Iterator>
     static VertexP3N call(Iterator begin) { return VertexP3N(Point3D(begin[0], begin[1], begin[2]), Vector3D(begin[3], begin[4], begin[5])); }

};

template<>
struct do_make_vertex<VertexP3C> {
     template<class Iterator>
     static VertexP3C call(Iterator begin) { return VertexP3C(Point3D(begin[0], begin[1], begin[2]), hpcolor(begin[3], begin[4], begin[5], begin[6])); }

};
template<class Vertex, class Iterator>
Vertex make_vertex(Iterator begin) { return do_make_vertex<Vertex>::call(begin); }

}//namespace happah

