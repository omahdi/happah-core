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

namespace happah {

class ReaderHPH : public std::ifstream {
public:
     ReaderHPH(const char* path);

     ~ReaderHPH();

     template<class T>
     static T read(const std::string& path) {
          T t;
          ReaderHPH reader(path.c_str());
          reader >> t;
          return t;
     }

};//ReaderHPH

ReaderHPH& operator>>(ReaderHPH& reader, boost::dynamic_bitset<>& bits);

template<class T>
ReaderHPH& operator>>(ReaderHPH& reader, std::vector<T>& ts) {
     size_t nTs;
     reader >> nTs;
     ts.clear();
     ts.reserve(nTs);
     for(auto i = 0lu; i < nTs; ++i) {
          T t;
          (reader) >> t;
          ts.push_back(t);
     }
     return reader;
}

ReaderHPH& operator>>(ReaderHPH& reader, Point2D& point);

ReaderHPH& operator>>(ReaderHPH& reader, Point3D& point);

ReaderHPH& operator>>(ReaderHPH& reader, Point4D& point);

}//namespace happah

