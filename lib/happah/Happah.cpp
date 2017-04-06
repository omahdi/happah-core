// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/Happah.h"
#include "happah/io/readers/ReaderHPH.h"

namespace happah {

Indices make_indices(const std::string& path) { return ReaderHPH::read<Indices>(path); }

std::vector<hpreal> make_reals(const std::string& path) { return ReaderHPH::read<std::vector<hpreal> >(path); }

}//namespace happah

