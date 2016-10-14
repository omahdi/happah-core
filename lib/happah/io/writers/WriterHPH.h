// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <vector>

#include "happah/math/Space.h"

namespace happah {

class WriterHPH : public std::ofstream {
public:
     WriterHPH(const char* path);

     ~WriterHPH();

     template<class T>
     static void write(const T& t, const std::string& path) { write(t, path.c_str()); }

     template<class T>
     static void write(const T& t, const char* path) {
          WriterHPH writer(path);
          writer << t;
     }

};//WriterHPH

WriterHPH& operator<<(WriterHPH& writer, const boost::dynamic_bitset<>& bits);

WriterHPH& operator<<(WriterHPH& writer, const Point3D& point);

WriterHPH& operator<<(WriterHPH& writer, const Point4D& point);

template<class T>
WriterHPH& operator<<(WriterHPH& writer, const std::vector<T>& ts) {
     writer << ts.size();
     for(auto& t : ts) writer << ' ' << t;
     return writer;
}

}//namespace happah

