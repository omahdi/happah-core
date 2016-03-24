/**********************************************************************************
 * Copyright (c) 2012-2015  See the COPYRIGHT-HOLDERS file in the top-level
 * directory of this distribution and at http://github.com/happah-graphics/happah.
 *
 * This file is part of Happah. It is subject to the license terms in the LICENSE
 * file found in the top-level directory of this distribution and at
 * http://github.com/happah-graphics/happah. No part of Happah, including this
 * file, may be copied, modified, propagated, or distributed except according to
 * the terms contained in the LICENSE file.
 **********************************************************************************/

#pragma once

#include <fstream>

#include "happah/math/HexagonDecomposition.h"

namespace happah {

class ReaderHPH : public std::ifstream {
public:
     ReaderHPH(const char* path)
          : std::ifstream(path) {
          if(this->fail()) {
               this->close();
               throw std::runtime_error("Failed to open file.");
          }
     }

     ~ReaderHPH() { this->close(); }

     template<class T>
     static void read(T& t, const char* path) {
          ReaderHPH reader(path);
          reader >> t;
     }

};//ReaderHPH

//TODO: to cpp
ReaderHPH& operator>>(ReaderHPH& reader, boost::dynamic_bitset<>& bits) {
     size_t nBits;
     reader >> nBits;
     bits.clear();
     bits.resize(nBits);
     for(auto i = 0lu; i < nBits; ++i) {
          char bit;
          reader >> bit;
          if(bit == '0') bits[i] = false;
          else bits[i] = true;
     }
     return reader;
}

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

ReaderHPH& operator>>(ReaderHPH& reader, Point3D& point) {
     reader >> point.x;
     reader >> point.y;
     reader >> point.z;
     return reader;
}

ReaderHPH& operator>>(ReaderHPH& reader, Point4D& point) {
     reader >> point.x;
     reader >> point.y;
     reader >> point.z;
     reader >> point.w;
     return reader;
}

}//namespace happah

