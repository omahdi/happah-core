// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <string>

#include "happah/Happah.h"

namespace happah {

template<hpindex t_format>
class Writer;

template<class Stream>
Stream& operator<<(Stream& stream, const hpvec2& v) { 
     stream << v.x << ' ' << v.y;
     return stream; 
}

template<class Stream>
Stream& operator<<(Stream& stream, const hpvec3& v) { 
     stream << v.x << ' ' << v.y << ' ' << v.z;
     return stream; 
}

template<class Stream>
Stream& operator<<(Stream& stream, const hpvec4& v) { 
     stream << v.x << ' ' << v.y << ' ' << v.z << ' ' << v.w;
     return stream; 
}

template<class Stream>
Stream& operator<<(Stream& stream, const hpmat3x3& m) { 
     stream << m[0] << ' ' << m[1] << ' ' << m[2];
     return stream; 
}

template<class Stream>
Stream& operator<<(Stream& stream, const VertexP3& vertex) {
     stream << vertex.position;
     return stream;
}

template<class Stream>
Stream& operator<<(Stream& stream, const VertexP3N& vertex) {
     stream << vertex.position << ' ' << vertex.normal;
     return stream;
}

}//namespace happah
