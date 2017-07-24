// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/math/Space.hpp"

namespace happah {

template<class Space>
struct Geometry {
     static_assert(is_space<Space>::value, "A geometry can only be parameterized by a space.");
     using SPACE = Space;

protected:
     Geometry() {}
     
};

template<class Space>
struct Geometry0D : public Geometry<Space> {
     static constexpr hpuint DIMENSION = 0;

protected:
     Geometry0D() {}

};//points

template<class Space>
struct Geometry1D : public Geometry<Space> {
     static constexpr hpuint DIMENSION = 1;

protected:
     Geometry1D() {}

};//curves

template<class Space>
struct Geometry2D : public Geometry<Space> {
     static constexpr hpuint DIMENSION = 2;

protected:
     Geometry2D() {}

};//surfaces

template<class Space>
struct Geometry3D : public Geometry<Space> {
     static constexpr hpuint DIMENSION = 3;

protected:
     Geometry3D() {}

};//three manifolds

template<class G, class Space = typename G::SPACE>
struct is_geometry : std::integral_constant<bool, std::is_base_of<Geometry<Space>, G>::value> {};

template<class Geometry, class Space = typename Geometry::SPACE> 
struct enable_if_geometry : std::enable_if<is_geometry<Geometry, Space>::value> {};

}//namespace happah

