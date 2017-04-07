// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/io/writers/WriterHPH.h"

namespace happah {

WriterHPH::WriterHPH(const char* path)
     : std::ofstream(path) {
     if(this->fail()) {
          this->close();
          throw std::runtime_error("Failed to open file.");
     }
}

WriterHPH::~WriterHPH() { this->close(); }

WriterHPH& operator<<(WriterHPH& writer, const boost::dynamic_bitset<>& bits) {
     writer << bits.size() << ' ';
     for(auto i = 0lu, end = bits.size(); i < end; ++i) writer << ((bits[i]) ? '1' : '0');
     return writer;
}

WriterHPH& operator<<(WriterHPH& writer, const Point3D& point) {
     writer << point.x << ' ' << point.y << ' ' << point.z;
     return writer;
}

WriterHPH& operator<<(WriterHPH& writer, const Point4D& point) {
     writer << point.x << ' ' << point.y << ' ' << point.z << ' ' << point.w;
     return writer;
}

WriterHPH& operator<<(WriterHPH& writer, const VertexP3& vertex) {
     writer << vertex.position;
     return writer;
}

WriterHPH& operator<<(WriterHPH& writer, const VertexP3N& vertex) {
     writer << vertex.position << ' ' << vertex.normal;
     return writer;
}

}//namespace happah

