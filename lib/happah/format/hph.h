// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "happah/format/default.h"
#include "happah/math/Space.h"
#include "happah/geometries/Vertex.h"

namespace happah {

namespace format {

namespace hph {

using happah::format::operator<<;
using happah::format::operator>>;

template<class T>
T read(const std::string& source);

template<class T>
T read(const std::experimental::filesystem::path& source);

template<class T>
void write(const T& t, const std::experimental::filesystem::path& target);

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
T read(const std::string& source) {
     auto stream = std::stringstream(source);
     auto t = T();

     stream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
     stream >> t;
     return t;
}

template<class T>
T read(const std::experimental::filesystem::path& source) {
     auto stream = std::ifstream();
     auto t = T();

     stream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
     stream.open(source);
     stream >> t;
     return t;
}

template<class T>
void write(const T& t, const std::experimental::filesystem::path& target) {
     auto stream = std::ofstream();

     stream.exceptions(std::ofstream::failbit | std::ofstream::badbit);
     stream.open(target);
     stream << t;
}

}//namespace hph

}//namespace format

template<class T>
void write(const T& t, const std::experimental::filesystem::path& target) { format::hph::write(t, target); }

}//namespace happah

