// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/CirclePackingMetric.hpp"

namespace happah {

CirclePackingMetric::CirclePackingMetric(std::vector<hpreal> radii, std::vector<hpreal> weights)
     : m_radii(std::move(radii)), m_weights(std::move(weights)) {}

const std::vector<hpreal>& CirclePackingMetric::getRadii() const { return m_radii; }

const std::vector<hpreal>& CirclePackingMetric::getWeights() const { return m_weights; }

CirclePackingMetric make_circle_packing_metric(std::vector<hpreal>& radii, std::vector<hpreal>& weights) { return { radii, weights }; }

}//namespace happah

