// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

class Space1D {
public:
     using POINT = Point1D;
     using VECTOR = Vector1D;
     static constexpr hpuint DIMENSION = 1;

};//Space1D

class Space2D {
public:
     using POINT = Point2D;
     using VECTOR = Vector2D;
     static constexpr hpuint DIMENSION = 2;

};//Space2D

class Space3D {
public:
     using POINT = Point3D;
     using VECTOR = Vector3D;
     static constexpr hpuint DIMENSION = 3;

};//Space3D

class Space4D {
public:
     using POINT = Point4D;
     using VECTOR = Vector4D;
     static constexpr hpuint DIMENSION = 4;

};//Space4D

}//namespace happah

