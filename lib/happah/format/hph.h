// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <string>
#include <vector>

#include "happah/format/default.h"
#include "happah/math/Space.h"
#include "happah/geometries/Vertex.h"
#include "happah/writers/Writer.h"

namespace happah {

namespace format {

namespace hph {

using happah::format::operator<<;
using happah::format::operator>>;

template<class T>
static T read(const std::string& path);

template<class T>
static void write(const T& t, const std::string& path);

template<class Stream>
Stream& operator<<(Stream& stream, const boost::dynamic_bitset<>& bits) {
     stream << bits.size() << ' ';
     for(auto i = 0lu, end = bits.size(); i < end; ++i) stream << ((bits[i]) ? '1' : '0');
     return stream;
}

template<class Stream, class T>
Stream& operator<<(Stream& stream, const std::vector<T>& ts) {
     stream << ts.size();
     for(auto& t : ts) {
          stream << ' ';
          stream << t;
     }
     return stream;
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

template<class T>
static T read(const std::string& path) {
     auto stream = std::ifstream();
     T t;

     stream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
     stream.open(path);
     stream >> t;
     return t;
}

template<class T>
static void write(const T& t, const std::string& path) {
     auto stream = std::ofstream();

     stream.exceptions(std::ofstream::failbit | std::ofstream::badbit);
     stream.open(path);
     stream << t;
}

}//namespace hph

}//namespace format

template<class T>
void write(const T& t, const std::string& path) { format::hph::write(t, path); }

}//namespace happah

