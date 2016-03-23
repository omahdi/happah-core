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

class WriterHPH : public std::ofstream {
public:
     WriterHPH(const char* path)
          : std::ofstream(path) {
          if(this->fail()) {
               this->close();
               throw std::runtime_error("Failed to open file.");
          }
     }

     ~WriterHPH() { this->close(); }

     template<class T>
     static void write(const T& t, const std::string& path) { write(t, path.c_str()); }

     template<class T>
     static void write(const T& t, const char* path) {
          WriterHPH writer(path);
          writer << t;
     }

};//WriterHPH

WriterHPH& operator<<(WriterHPH& writer, const boost::dynamic_bitset<>& bits) {
     writer << bits.size() << ' ';
     for(auto i = 0lu, end = bits.size(); i < end; ++i) writer << ((bits[i]) ? '1' : '0');
     writer << '\n' << '\n';
     return writer;
}

WriterHPH& operator<<(WriterHPH& writer, const Point3D& point) {
     writer << point.x << ' ' << point.y << ' ' << point.z << '\n';
     return writer;
}

WriterHPH& operator<<(WriterHPH& writer, const Point4D& point) {
     writer << point.x << ' ' << point.y << ' ' << point.z << ' ' << point.w << '\n';
     return writer;
}

//TODO: to cpp
template<class T>
WriterHPH& operator<<(WriterHPH& writer, const std::vector<T>& ts) {
     writer << ts.size() << '\n';
     for(auto& t : ts) writer << t << ' ';
     writer << '\n' << '\n';
     return writer;
}

}//namespace happah

