// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>

#include "happah/Happah.h"
#include "happah/format/hph.h"

namespace happah {

Indices::iterator defrag(Indices& indices) {
     auto j = std::begin(indices) - 1;

     for(auto i = std::begin(indices), end = std::end(indices); i != end; ++i) {
          if(*i == std::numeric_limits<hpindex>::max()) continue;
          *(++j) = *i;
     }

     return j + 1;
}

Indices make_indices(const std::string& path) { return format::hph::read<Indices>(path); }

std::vector<hpreal> make_reals(const std::string& path) { return format::hph::read<std::vector<hpreal> >(path); }

std::string slurp(const std::string& path) {
     auto file = std::ifstream(path);
     return { (std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()) };
}

}//namespace happah

