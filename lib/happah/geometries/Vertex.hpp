// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/geometries/Geometry.hpp"
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
class Vertex : public Geometry0D<Space> {
     static_assert(is_space<Space>::value, "A vertex can only be parameterized by a space.");

protected:
     Vertex() {}

};

template<class SpaceA, class SpaceO>
class VertexAO : public Vertex<typename get_combined_space<SpaceA, SpaceO>::SPACE> {
     static_assert(is_space<SpaceA>::value && is_space<SpaceO>::value, "An abscissa and ordinate vertex can only be parameterized by two spaces.");
     typedef typename SpaceA::POINT Abscissa;
     typedef typename SpaceO::POINT Ordinate;

public:
     typedef SpaceA SPACEA;
     typedef SpaceO SPACEO;

     Abscissa abscissa;
     Ordinate ordinate;

     VertexAO() : abscissa{0.0}, ordinate{0.0} {}
     VertexAO(Abscissa abscissa, Ordinate ordinate) : abscissa(std::move(abscissa)), ordinate(std::move(ordinate)) {}

};
typedef VertexAO<Space2D, Space1D> VertexA2O1;

template<class SpaceA, class SpaceO>
class VertexAON : public VertexAO<SpaceA, SpaceO> {
     typedef typename SpaceA::POINT Abscissa;
     typedef typename SpaceO::POINT Ordinate;
     typedef typename VertexAO<SpaceA, SpaceO>::SPACE::VECTOR Vector;

public:   
     Vector normal;

     VertexAON() : normal{0.0} {}
     VertexAON(Abscissa abscissa, Ordinate ordinate) : VertexAO<SpaceA, SpaceO>(std::move(abscissa), std::move(ordinate)), normal(0.0) {}
     VertexAON(Abscissa abscissa, Ordinate ordinate, Vector normal) : VertexAO<SpaceA, SpaceO>(std::move(abscissa), std::move(ordinate)), normal(std::move(normal)) {}

};
typedef VertexAON<Space2D, Space1D> VertexA2O1N;

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
struct contains_abscissa_color : public std::false_type {};

template<class V, typename = void>
struct contains_abscissa_ordinate : public std::false_type {};

template<class V>
struct contains_abscissa_ordinate<V, typename std::enable_if<std::is_base_of<VertexAO<typename V::SPACEA, typename V::SPACEO>, V>::value>::type> : public std::true_type {};

template<class Vertex, typename = void>
struct contains_position : public std::false_type {};

template<class Vertex>
struct contains_position<Vertex, typename std::enable_if<std::is_base_of<VertexP<typename Vertex::SPACE>, Vertex>::value || std::is_base_of<VertexPN<typename Vertex::SPACE>, Vertex>::value>::type> : public std::true_type {};

template<class V>
struct contains_color : public std::integral_constant<bool, std::is_base_of<VertexPC<typename V::SPACE>, V>::value || std::is_base_of<VertexPNC<typename V::SPACE>, V>::value> {};

template<class Vertex, typename = void>
struct contains_ordinate_color : public std::false_type {};

template<class Vertex, typename = void>
struct get_abscissa_color_offset : public std::integral_constant<size_t, 0> {};

//TODO: check if we can use std::integral_constant in other places
template<class Vertex, typename = void>
struct get_abscissa_dimension : public std::integral_constant<hpuint, 0> {};

template<class Vertex>
struct get_abscissa_dimension<Vertex, typename std::enable_if<std::is_base_of<VertexAO<typename Vertex::SPACEA, typename Vertex::SPACEO>, Vertex>::value>::type> : public std::integral_constant<hpuint, Vertex::SPACEA::DIMENSION> {};

template<class Vertex, class Space = typename Vertex::SPACE>
struct get_color_offset : public std::integral_constant<size_t, 0> {};

template<class Space>
struct get_color_offset<VertexPC<Space> > : public std::integral_constant<size_t, sizeof(typename Space::POINT)> {};

template<class Space>
struct get_color_offset<VertexPNC<Space> > : public std::integral_constant<size_t, sizeof(typename Space::POINT) + sizeof(typename Space::VECTOR)> {};

template<class Vertex, class Space = typename Vertex::SPACE>
struct get_normal_offset : public std::integral_constant<size_t, 0> {};

template<class Space>
struct get_normal_offset<VertexPN<Space> > : public std::integral_constant<size_t, sizeof(typename Space::POINT)> {};

template<class SpaceA, class SpaceO>
struct get_normal_offset<VertexAON<SpaceA, SpaceO> > : public std::integral_constant<size_t, sizeof(typename SpaceA::POINT) + sizeof(typename SpaceO::POINT)> {};

template<class Space>
struct get_normal_offset<VertexPNC<Space> > : public std::integral_constant<size_t, sizeof(typename Space::POINT)> {};

template<class Vertex, typename = void>
struct get_ordinate_color_offset : public std::integral_constant<size_t, 0> {};

template<class Vertex, typename = void>
struct get_ordinate_offset : public std::integral_constant<size_t, 0> {};

//TODO: is_same taking array of template template arguments and checking if it is one of them
template<class Vertex>
struct get_ordinate_offset<Vertex, typename std::enable_if<std::is_same<VertexAO<typename Vertex::SPACEA, typename Vertex::SPACEO>, Vertex>::value || std::is_same<VertexAON<typename Vertex::SPACEA, typename Vertex::SPACEO>, Vertex>::value>::type> : public std::integral_constant<size_t, sizeof(typename Vertex::SPACEA::POINT)> {};

template<class V, class Space = typename V::SPACE>
struct is_vertex : public std::integral_constant<bool, std::is_base_of<Vertex<Space>, V>::value> {};

template<class V, class Space = typename V::SPACE>
struct enable_if_vertex : public std::enable_if<is_vertex<V, Space>::value, V> {};

//NOTE: The location of the vertex is given by a position (i.e. coordinates) with respect to some basis of the corresponding space.
template<class Vertex, class Space = typename Vertex::SPACE>
struct is_absolute_vertex : public std::integral_constant<bool, is_vertex<Vertex, Space>::value && contains_position<Vertex>::value> {};

//NOTE: The location of the vertex is given by an abscissa and ordinate with respect to some n-dimensional structure that serves as a base.  This structure can be a manifold or a vector space.  If the base is a vector space, then the interpretation of the abscissa and ordinate of the vertex coincides with the interpretation of a vertex positioned absolutely.  The ordinate is plotted on a normal of the underlying base.
template<class Vertex>
struct is_relative_vertex : public std::integral_constant<bool, is_vertex<Vertex>::value && contains_abscissa_ordinate<Vertex>::value> {};

template<class Geometry1, class Geometry2, typename = void>
struct is_relativizable : public std::false_type {};

template<class Vertex, class Geometry>
struct is_relativizable<Vertex, Geometry, typename std::enable_if<is_relative_vertex<Vertex>::value && is_geometry<Geometry>::value && Vertex::SPACEA::DIMENSION == Geometry::DIMENSION>::type> : public std::true_type {};

template<class Vertex, class Space = typename Vertex::SPACE>
struct enable_if_absolute_vertex : public std::enable_if<is_absolute_vertex<Vertex, Space>::value> {};

template<class Vertex>
struct enable_if_relative_vertex : public std::enable_if<is_relative_vertex<Vertex>::value> {};

template<class Vertex, typename = void>
struct contains_normal;

template<class Vertex>
struct contains_normal<Vertex, typename enable_if_absolute_vertex<Vertex>::type> : public std::integral_constant<bool, std::is_base_of<VertexPN<typename Vertex::SPACE>, Vertex>::value> {};

template<class Vertex>
struct contains_normal<Vertex, typename enable_if_relative_vertex<Vertex>::type> : public std::integral_constant<bool, std::is_base_of<VertexAON<typename Vertex::SPACEA, typename Vertex::SPACEO>, Vertex>::value> {};

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

