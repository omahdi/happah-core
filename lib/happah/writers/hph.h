// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <string>
#include <vector>

#include "happah/writers/Writer.h"

namespace happah {

namespace hph {

using happah::operator<<;

template<class T>
static void write(const T& t, const std::string& path) {
     auto stream = std::ofstream();

     stream.exceptions(std::ofstream::failbit | std::ofstream::badbit);
     stream.open(path);
     stream << t;
}

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

}//namespace hph

template<class T>
void write(const T& t, const std::string& path) { hph::write(t, path); }

}//namespace happah

