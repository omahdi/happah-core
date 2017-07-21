// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/ProjectiveStructure.hpp"

namespace happah {

ProjectiveStructure::ProjectiveStructure(Indices neighbors, std::vector<hpreal> transitions)
     : m_neighbors(std::move(neighbors)), m_transitions(std::move(transitions)) {}

const Indices& ProjectiveStructure::getNeighbors() const { return m_neighbors; }

const std::vector<hpreal>& ProjectiveStructure::getTransitions() const { return m_transitions; }

ProjectiveStructure make_projective_structure(Indices neighbors, std::vector<hpreal> transitions) { return { std::move(neighbors), std::move(transitions) }; }

}//namespace happah

