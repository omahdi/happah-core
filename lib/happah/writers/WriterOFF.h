// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <fstream>
#include <string>

#include "happah/geometries/TriangleMesh.h"
#include "happah/writers/Writer.h"

namespace happah {

constexpr hpindex OFF = 1001;

struct Header {
     bool color;
     hpuint dimension;
     hpuint nFaces;
     bool normal;
     hpuint nVertices;
};

template<class Vertex>
Header make_header(hpuint nFaces, hpuint nVertices) {
     auto header = Header();
     header.color = contains_color<Vertex>::value;
     header.dimension = Vertex::SPACE::DIMENSION;
     header.nFaces = nFaces;
     header.normal = contains_normal<Vertex>::value;
     header.nVertices = nVertices;
     return header;
}

template<>
class Writer<OFF> : public std::ofstream {
public:
     Writer(const std::string& path);

     ~Writer();

     template<class Vertex, Format format>
     static void write(const TriangleMesh<Vertex, format>& mesh, const std::string& path) {
          static_assert(is_absolute_vertex<Vertex>::value, "This write method can only be parameterized by absolute vertices.");

          Writer<OFF> writer(path);
          auto& indices = mesh.getIndices();
          auto& vertices = mesh.getVertices();

          writer << make_header<Vertex>(size(vertices), size(mesh)) << "\n\n";
          writer << vertices << "\n\n";
          visit_triplets(indices, [&](auto i0, auto i1, auto i2) { writer << "3 " << i0 << ' ' << i1 << ' ' << i2 << '\n'; });
     }

};//Writer<OFF>

template<class Stream>
Stream& operator<<(Stream& stream, const Header& header) {
     if(header.color) stream << 'C';
     if(header.normal) stream << 'N';
     if(header.dimension == 4) stream << "4OFF\n";
     else if(header.dimension == 3) stream << "OFF\n";
     else {
          stream << "nOFF\n";
          stream << header.dimension << '\n';
     }
     stream << header.nVertices << ' ' << header.nFaces << " 0";
}

template<class Vertex>
Writer<OFF>& operator<<(Writer<OFF>& writer, const std::vector<Vertex>& vertices) {
     for(auto& vertex : vertices) writer << vertex << '\n';
     return writer;
}

}//namespace happah

