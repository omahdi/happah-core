// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/math/MathUtils.hpp"

class CurveUtilsBEZ {
     /**
      * Evaluate the Bernstein polynomial B^n_i.
      */
     template<hpuint t_degree, bool t_check = true>
     static hpreal evaluate(hpreal t, hpuint i) { return evaluate<t_degree, t_check>(t, 1.0 - t, i, t_degree - i); }

     template<hpuint t_degree, bool t_check = true>
     static hpreal evaluate(hpreal t, hpreal oneMinusT, hpuint i, hpuint degreeMinusI) { return MathUtils::binom<t_check>(t_degree, i) * MathUtils::pow(t, i) * MathUtils::pow(oneMinusT, degreeMinusI); }

};

