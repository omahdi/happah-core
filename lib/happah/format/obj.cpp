// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/iostreams/device/mapped_file.hpp>

#include "happah/format/obj.hpp"

namespace happah {

namespace format {

namespace obj {

Content read(const std::string& path) {
     using boost::iostreams::mapped_file;
     mapped_file file(path, mapped_file::readonly);
     auto begin = file.const_data();
     auto end = begin + file.size();
     return read(begin, end);
}

}//namespace obj

}//namespace format

}//namespace happah


