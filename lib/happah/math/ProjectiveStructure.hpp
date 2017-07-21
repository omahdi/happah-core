// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <vector>

#include "happah/Happah.hpp"

namespace happah {

class ProjectiveStructure;

ProjectiveStructure make_projective_structure(Indices neighbors, std::vector<hpreal> transitions);

class ProjectiveStructure {
public:
     ProjectiveStructure(Indices neighbors, std::vector<hpreal> transitions);

     const Indices& getNeighbors() const;

     const std::vector<hpreal>& getTransitions() const;

private:
     Indices m_neighbors;
     std::vector<hpreal> m_transitions;

};//ProjectiveStructure

}//namespace happah

