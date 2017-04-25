// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <fstream>

#include "happah/formats/off.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/writers/Writer.h"

namespace happah {

namespace off {

using happah::operator<<;

template<class Vertex, Format format>
void write(const TriangleMesh<Vertex, format>& mesh, const std::string& path) {
     auto& vertices = mesh.getVertices();
     auto& indices = mesh.getIndices();
     auto stream = std::ofstream();

     stream.exceptions(std::ofstream::failbit);
     stream.open(path);
     stream << make_header<Vertex>(size(mesh), size(vertices)) << "\n\n";
     stream << std::fixed;
     stream << vertices << "\n\n";
     visit_triplets(indices, [&](auto i0, auto i1, auto i2) { stream << "3 " << i0 << ' ' << i1 << ' ' << i2 << '\n'; });
}

template<class Stream>
Stream& operator<<(Stream& stream, const Header& header) {
     if(header.color) stream << 'C';
     if(header.normal) stream << 'N';
     if(header.dimension == 3) stream << "OFF\n";
     else if(header.dimension == 4) stream << "4OFF\n";
     else {
          stream << "nOFF\n";
          stream << header.dimension << '\n';
     }
     stream << header.nVertices << ' ' << header.nFaces << " 0";
     return stream;
}

template<class Stream, class T>
Stream& operator<<(Stream& stream, const std::vector<T>& ts) {
     for(auto t = std::begin(ts), end = std::end(ts) - 1; t != end; ++t) {
          stream << *t;
          stream << '\n';
     }
     stream << ts.back();
     return stream;
}

}//namespace off

}//namespace happah

