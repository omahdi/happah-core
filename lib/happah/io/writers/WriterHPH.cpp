// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/io/writers/WriterHPH.h"

namespace happah {

Writer<HPH>::Writer(const std::string& path)
     : std::ofstream(path.c_str()) {
     if(this->fail()) {
          this->close();
          throw std::runtime_error("Failed to open file.");
     }
}

Writer<HPH>::~Writer() { this->close(); }

Writer<HPH>& operator<<(Writer<HPH>& writer, const boost::dynamic_bitset<>& bits) {
     writer << bits.size() << ' ';
     for(auto i = 0lu, end = bits.size(); i < end; ++i) writer << ((bits[i]) ? '1' : '0');
     return writer;
}

Writer<HPH>& operator<<(Writer<HPH>& writer, const VertexP3& vertex) {
     writer << vertex.position;
     return writer;
}

Writer<HPH>& operator<<(Writer<HPH>& writer, const VertexP3N& vertex) {
     writer << vertex.position << ' ' << vertex.normal;
     return writer;
}

}//namespace happah

