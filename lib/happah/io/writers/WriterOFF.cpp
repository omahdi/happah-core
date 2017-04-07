// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/io/writers/WriterOFF.h"

namespace happah {

Writer<OFF>::Writer(const std::string& path)
     : std::ofstream(path.c_str()) {
     if(this->fail()) {
          this->close();
          throw std::runtime_error("Failed to open file.");
     }
}

Writer<OFF>::~Writer() { this->close(); }

}//namespace happah

