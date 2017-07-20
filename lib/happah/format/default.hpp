// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/geometries/Vertex.hpp"

namespace happah {

namespace format {

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
Stream& operator<<(Stream& stream, const hpucolor& v) {
     stream << v.r << ' ' << v.g << ' ' << v.b << ' ' << v.a;
     return stream;
}

template<class Stream, class Space>
Stream& operator<<(Stream& stream, const VertexP<Space>& vertex) {
     stream << vertex.position;
     return stream;
}

template<class Stream, class Space>
Stream& operator<<(Stream& stream, const VertexPN<Space>& vertex) {
     stream << vertex.position << ' ' << vertex.normal;
     return stream;
}

template<class Stream, class Space>
Stream& operator<<(Stream& stream, const VertexPC<Space>& vertex) {
     stream << vertex.position << ' ' << vertex.color;
     return stream;
}

template<class Stream, class Space>
Stream& operator<<(Stream& stream, const VertexPNC<Space>& vertex) {
     stream << vertex.position << ' ' << vertex.normal << ' ' << vertex.color;
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, hpvec2& point) {
     stream >> point.x;
     stream >> point.y;
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, hpvec3& point) {
     stream >> point.x;
     stream >> point.y;
     stream >> point.z;
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, hpvec4& point) {
     stream >> point.x;
     stream >> point.y;
     stream >> point.z;
     stream >> point.w;
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, VertexP3& vertex) {
     stream >> vertex.position;
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, VertexP3N& vertex) {
     stream >> vertex.position;
     stream >> vertex.normal;
     return stream;
}

}//namespace format

}//namespace happah

