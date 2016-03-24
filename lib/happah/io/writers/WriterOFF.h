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
#include <string>

#include "happah/geometries/TriangleMesh.h"

namespace happah {

//TODO: Triangle indices as std::vector<TriangleIndices>* or std::vector<hpuint3>*
class WriterOFF {
public:
     template<class Vertex>
     static bool write(const TriangleMesh<Vertex>& triangleMesh, const std::string& path) { return write(triangleMesh, path.c_str()); }

     template<class Vertex>
     static bool write(const TriangleMesh<Vertex>& triangleMesh, const char* path) {
          std::ofstream out(path);
          if(out.fail()) {
               out.close();
               return false;
          }
          write(triangleMesh, out);
          out.close();
          return true;
     }

     template<class Stream, class Vertex>
     static void write(const TriangleMesh<Vertex>& triangleMesh, Stream& out) {//TODO: include format?
          static_assert(is_absolute_vertex<Vertex>::value, "This write method can only be parameterized by absolute vertices.");

          write<Stream, Vertex>(out);
          auto& indices = triangleMesh.getIndices();
          auto& vertices = triangleMesh.getVertices();
          out << vertices.size() << ' ' << indices.size() / 3 << " 0\n\n";
          write<Stream, Vertex>(vertices, out);
          out << '\n';
          for(auto i = indices.cbegin(), end = indices.cend(); i != end; ++i) {
               out << "3 " << *i;
               ++i;
               out << ' ' << *i;
               ++i;
               out << ' ' << *i << '\n';
          }
          out << '\n';
     }

private:
     template<class Vertex, typename = void>
     struct do_contains_color;

     template<class Vertex>
     struct do_contains_color<Vertex, typename enable_if_absolute_vertex<Vertex>::type> : public contains_color<Vertex> {};

     template<class Vertex>
     struct do_contains_color<Vertex, typename enable_if_relative_vertex<Vertex>::type> : public contains_ordinate_color<Vertex> {};

     template<class Vertex, typename = void>
     struct write_normal {
          template<class Stream>
          static void exec(const Vertex& v, Stream& out) {}

     };

     template<class Vertex>
     struct write_normal<Vertex, typename std::enable_if<contains_normal<Vertex>::value>::type> {
          template<class Stream>
          static void exec(const Vertex& v, Stream& out) { out << ' ' << v.normal; }

     };

     template<class Stream, class Vertex>
     static void write(Stream& out) {
          if(do_contains_color<Vertex>::value) out << 'C';
          if(contains_normal<Vertex>::value) out << 'N';
          if(Vertex::SPACE::DIMENSION == 4) out << "4OFF\n";
          else if(Vertex::SPACE::DIMENSION == 3) out << "OFF\n";
          else {
               out << "nOFF\n";
               out << Vertex::SPACE::DIMENSION << '\n';
          }
     }

     template<class Stream, class Vertex>
     static void write(const typename Model<Vertex>::Vertices& vertices, Stream& out) {
          for(auto i = vertices.cbegin(), end = vertices.cend(); i != end; ++i) {
               auto& v = *i;
               out << v.position;
               write_normal<Vertex>::exec(v, out);
               //TODO: color
               out << '\n';
          }
     }

};

}

