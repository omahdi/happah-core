// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/Triangle.h"
#include "happah/math/MathUtils.h"
#include "happah/math/Space.h"

namespace happah {

constexpr hpuint make_control_polygon_size(hpuint degree);

/**
 * @param[in] nSamples Number of times an edge of the parameter triangle should be sampled.  The entire triangle is sampled uniformly such that this parameter is respected.
 * @return Matrix whose rows are the Bernstein polynomials evaluated at some point u.  The matrix is returned row-major.  To evaluate a B\'ezier polynomial at the sampled u values given a vector of control points, simply compute the product of the matrix with the vector of control points.
 */
std::vector<hpreal> make_evaluation_matrix(hpuint degree, hpuint nSamples);

constexpr hpuint make_offset(hpuint degree, hpuint i0, hpuint i1, hpuint i2);

constexpr hpuint make_patch_size(hpuint degree);

constexpr hpuint make_control_polygon_size(hpuint degree) { return degree * degree; }

constexpr hpuint make_offset(hpuint degree, hpuint i0, hpuint i1, hpuint i2) { return make_patch_size(degree) - make_patch_size(degree - i2) + i1; }

constexpr hpuint make_patch_size(hpuint degree) { return (degree + 1) * (degree + 2) >> 1; }

}//namespace happah

