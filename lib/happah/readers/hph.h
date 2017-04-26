// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <string>
#include <vector>

#include "happah/math/Space.h"
#include "happah/geometries/Vertex.h"

namespace happah {

template<class Stream>
Stream& operator>>(Stream& stream, Point2D& point) {
     stream >> point.x;
     stream >> point.y;
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, Point3D& point) {
     stream >> point.x;
     stream >> point.y;
     stream >> point.z;
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, Point4D& point) {
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

namespace hph {

using happah::operator>>;

template<class T>
static T read(const std::string& path) {
     auto stream = std::ifstream();
     T t;

     stream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
     stream.open(path);
     stream >> t;
     return t;
}

template<class Stream>
Stream& operator>>(Stream& stream, boost::dynamic_bitset<>& bits) {
     size_t nBits;
     stream >> nBits;
     bits.clear();
     bits.resize(nBits, false);
     for(auto i = 0lu; i < nBits; ++i) {
          char bit;
          stream >> bit;
          if(bit == '1') bits[i] = true;
     }
     return stream;
}

template<class Stream, class T>
Stream& operator>>(Stream& stream, std::vector<T>& ts) {
     size_t nTs;
     stream >> nTs;
     ts.clear();
     ts.reserve(nTs);
     for(auto i = 0lu; i < nTs; ++i) {
          T t;
          (stream) >> t;
          ts.push_back(t);
     }
     return stream;
}

}//namespace hph

}//namespace happah

