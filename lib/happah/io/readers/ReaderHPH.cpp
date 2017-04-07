// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/io/readers/ReaderHPH.h"

namespace happah {

ReaderHPH::ReaderHPH(const char* path)
     : std::ifstream(path) {
     if(this->fail()) {
          this->close();
          throw std::runtime_error("Failed to open file.");
     }
}

ReaderHPH::~ReaderHPH() { this->close(); }

ReaderHPH& operator>>(ReaderHPH& reader, boost::dynamic_bitset<>& bits) {
     size_t nBits;
     reader >> nBits;
     bits.clear();
     bits.resize(nBits, false);
     for(auto i = 0lu; i < nBits; ++i) {
          char bit;
          reader >> bit;
          if(bit == '1') bits[i] = true;
     }
     return reader;
}

ReaderHPH& operator>>(ReaderHPH& reader, Point2D& point) {
     reader >> point.x;
     reader >> point.y;
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

ReaderHPH& operator>>(ReaderHPH& reader, VertexP3& vertex) {
     reader >> vertex.position;
     return reader;
}

ReaderHPH& operator>>(ReaderHPH& reader, VertexP3N& vertex) {
     reader >> vertex.position;
     reader >> vertex.normal;
     return reader;
}

}//namespace happah


