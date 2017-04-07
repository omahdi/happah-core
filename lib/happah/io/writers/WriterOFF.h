// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <fstream>
#include <string>

#include "happah/geometries/TriangleMesh.h"
#include "happah/io/writers/Writer.h"

namespace happah {

constexpr hpindex OFF = 1001;

template<>
class Writer<OFF> : public std::ofstream {
public:
     Writer(const std::string& path);

     ~Writer();

     template<class Vertex, Format format>
     static void write(const TriangleMesh<Vertex, format>& mesh, const std::string& path) {
          static_assert(is_absolute_vertex<Vertex>::value, "This write method can only be parameterized by absolute vertices.");

          Writer<OFF> writer(path);
          write<Vertex>(writer);
          auto& vertices = mesh.getVertices();
          writer << size(vertices) << ' ' << size(mesh) << " 0\n\n";
          write(writer, vertices);
          writer << '\n';
          visit_triplets(mesh.getIndices(), [&](auto i0, auto i1, auto i2) { writer << "3 " << i0 << ' ' << i1 << ' ' << i2 << '\n'; });
          writer << '\n';
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

     template<class Vertex, class Stream>
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

     template<class Vertex, class Stream>
     static void write(Stream& out, const std::vector<Vertex>& vertices) {
          for(auto& vertex : vertices) {
               out << vertex.position;
               write_normal<Vertex>::exec(vertex, out);
               //TODO: color
               out << '\n';
          }
     }

};//Writer<OFF>

}//namespace happah

