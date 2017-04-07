// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <vector>

#include "happah/math/Space.h"
#include "happah/geometries/Vertex.h"
#include "happah/writers/Writer.h"

namespace happah {

constexpr hpindex HPH = 1000;

template<>
class Writer<HPH> : public std::ofstream {
public:
     Writer(const std::string& path);

     ~Writer();

     template<class T>
     static void write(const T& t, const std::string& path) {
          Writer<HPH> writer(path);
          writer << t;
     }

};//Writer<HPH>

Writer<HPH>& operator<<(Writer<HPH>& writer, const boost::dynamic_bitset<>& bits);

template<class T>
Writer<HPH>& operator<<(Writer<HPH>& writer, const std::vector<T>& ts) {
     writer << ts.size();
     for(auto& t : ts) {
          writer << ' ';
          writer << t;
     }
     return writer;
}

template<class T>
void write(const T& t, const std::string& path) { Writer<HPH>::write(t, path); }

}//namespace happah

