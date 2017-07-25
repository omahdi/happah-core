// Copyright 2017
//   Obada Mahdi - Karlsruhe Institute of Technology - omahdi@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Based on hph.h

#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <type_traits>

#include "happah/format/default.hpp"

namespace happah {

namespace format {

namespace xyz {

template<class T>
static T read(const std::string& path);

template<class T>
static void write(const T& t, const std::string& path);

template<class Stream, class Vertex>
Stream& operator<<(Stream& stream, const Vertex& v) {
     using happah::format::operator<<;
     return stream << v.position;
}

template<class Stream, class T>
Stream& operator<<(Stream& stream, const std::vector<T>& ts) {
     using happah::format::operator<<;
     for(auto& t : ts)
          stream << t << '\n';
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, hpvec2& v) {
     using happah::format::operator>>;
     std::string coords;
     if (!!std::getline(stream, coords)) {
          std::istringstream s(coords);
          s >> v.x >> v.y;
     }
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, hpvec3& v) {
     using happah::format::operator>>;
     std::string coords;
     if (!!std::getline(stream, coords)) {
          std::istringstream s(coords);
          s >> v.x >> v.y >> v.z;
     }
     return stream;
}

template<class Stream>
Stream& operator>>(Stream& stream, hpvec4& v) {
     using happah::format::operator>>;
     std::string coords;
     if (!!std::getline(stream, coords)) {
          std::istringstream s(coords);
          s >> v.x >> v.y >> v.z >> v.w;
     }
     return stream;
}

template<class Stream, class Vertex>
Stream& operator>>(Stream& stream, Vertex& v) {
     decltype(v.position) t;
     stream >> t;
     if (!stream.eof() && !stream.fail())
          v.position = t;
     return stream;
}

template<class Stream, class T>
Stream& operator>>(Stream& stream, std::vector<T>& ts) {
     ts.clear();
     while (!stream.eof()) {
          T t;
          stream >> t;
          if (stream.eof() || stream.fail())
               break;
          ts.push_back(t);
     }
     return stream;
}

template<class T>
static T read(const std::string& path) {
     auto stream = std::ifstream();
     T t;

     stream.exceptions(std::ifstream::badbit);
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

}//namespace xyz

}//namespace format

}//namespace happah

