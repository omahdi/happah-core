#pragma once

#include <array>

#include "happah/math/BezierPatch.h"

namespace happah {

namespace loop {

template <typename Vertex>
using LoopPatch = std::array<Vertex, 12>;

namespace detail {
    extern const LoopPatch<QuarticBezierPatch<int>> bases;
}

/**
 * Computes the Bézier coefficients of a quartic box spline patch.
 * Order of the input vertices & basis functions
 *     10--11
 *    / \ / \
 *   7---8---9
 *  / \ / \ / \
 * 3---4---5---6
 *  \ / \ / \ /
 *   0---1---2
 * Order of the patch triangle vertices
 *       1
 *      / \
 *     0---2
 */
template <typename Vertex>
QuarticBezierPatch<Vertex> loop_to_bezier(const LoopPatch<Vertex>& vertices) {
    QuarticBezierPatch<Vertex> result{};
    const int nCoefficients = detail::bases[0].size();
    const int nBases = vertices.size();
    for (int j = 0; j < nCoefficients; ++j) {
        for (int i = 0; i < nBases; ++i) {
            result[j] += vertices[i] * detail::bases[i][j] / 24.0;
        }
    }
    return result;
}

} // namespace loop

} // namespace happah
