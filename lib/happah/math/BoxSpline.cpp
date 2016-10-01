#include "BoxSpline.h"

#include <algorithm>
#include <tuple>

namespace happah {

namespace loop {

namespace {
using Basis = QuarticBezierPatch<int>;

// HACK to allow constexpr access to std::array
template<typename ControlPoints = decltype(bezier_patch_ctrl_points<4>()),
         typename ControlPointIndices = std::array<std::array<int, 5>, 5>>
/* Compute the order of coefficients we have to use.
 * Thus, the implementation is order-independent. */
constexpr QuarticBezierPatch<int> coefficient_permutation() {
    ControlPoints ctrl_pts = bezier_patch_ctrl_points<4>();
    ControlPointIndices c{};
    for (int i = 0; i < 15; ++i) {
        c[ctrl_pts[i].r][ctrl_pts[i].s] = i;
    }
    return {c[4][0],
            c[3][1], c[3][0],
            c[2][2], c[2][1], c[2][0],
            c[1][3], c[1][2], c[1][1], c[1][0],
            c[0][4], c[0][3], c[0][2], c[0][1], c[0][0],
    };
}

// HACK to allow constexpr access to std::array
template<typename Array = Basis>
/*
 * This method is used for reordering the Bézier coefficients.
 *
 * Order of the parameters
 *             s
 *           10
 *         6 11
 *       3 7 12
 *     1 4 8 13
 *   0 2 5 9 14
 *  r          t
 */
constexpr Basis basis(int b400, int b310, int b301, int b220, int b211, int b202, int b130, int b121,
                      int b112, int b103, int b040, int b031, int b022, int b013, int b004) {

    Array result {{b400, b310, b301, b220, b211, b202, b130, b121, b112, b103, b040, b031, b022, b013, b004}};
    Array result_permuted{};
    auto perm = coefficient_permutation();
    for (int i = 0; i < 15; ++i) {
        result_permuted[perm[i]] = result[i];
    }
    const auto& r = result_permuted;
    return Basis{r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7],
                 r[8], r[9],r[10],r[11],r[12],r[13],r[14]};
}
/* Box spline
 *                 0---0---0---0---0---0---0---0---0
 *                / \             / \             / \
 *               0   0   0   0   0   0   0   0   0   0
 *              /     \         /     \         /     \
 *             0   0   0   0   0   0   0   0   0   0   0
 *            /         \     /         \     /         \
 *           0   0   0   0   1   1   1   1   0   0   0   0
 *          /             \ /             \ /             \
 *         0---0---0---1---2---3---4---3---2---1---0---0---0
 *        / \             / \             / \             / \
 *       0   0   0   1   3   4   6   6   4   3   1   0   0   0
 *      /     \         /     \         /     \         /     \
 *     0   0   0   1   4   6   8  10   8   6   4   1   0   0   0
 *    /         \     /         \     /         \     /         \
 *   0   0   0   1   3   6  10  12  12  10   6   3   1   0   0   0
 *  /             \ /             \ /             \ /             \
 * 0---0---0---0---2---4---8--12--12--12---8---4---2---0---0---0---0
 *  \             / \             / \             / \             /
 *   0   0   0   1   3   6  10  12  12  10   6   3   1   0   0   0
 *    \         /     \         /     \         /     \         /
 *     0   0   0   1   4   6   8  10   8   6   4   1   0   0   0
 *      \     /         \     /         \     /         \     /
 *       0   0   0   1   3   4   6   6   4   3   1   0   0   0
 *        \ /             \ /             \ /             \ /
 *         0---0---0---1---2---3---4---3---2---1---0---0---0
 *          \             / \             / \             /
 *           0   0   0   0   1   1   1   1   0   0   0   0
 *            \         /     \         /     \         /
 *             0   0   0   0   0   0   0   0   0   0   0
 *              \     /         \     /         \     /
 *               0   0   0   0   0   0   0   0   0   0
 *                \ /             \ /             \ /
 *                 0---0---0---0---0---0---0---0---0
 *
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
} // namespace

namespace detail {

extern const LoopPatch<Basis> bases = {
    //     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
    basis( 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), // 0
    basis( 2, 1, 3, 0, 1, 4, 0, 0, 1, 3, 0, 0, 0, 1, 2), // 1
    basis( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2), // 2
    basis( 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), // 3
    basis(12,12,12, 8,10, 8, 4, 6, 6, 4, 2, 3, 4, 3, 2), // 4
    basis( 2, 3, 4, 4, 6, 8, 3, 6,10,12, 2, 4, 8,12,12), // 5
    basis( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2), // 6
    basis( 2, 3, 1, 4, 1, 0, 3, 1, 0, 0, 2, 1, 0, 0, 0), // 7
    basis( 2, 4, 3, 8, 6, 4,12,10, 6, 3,12,12, 8, 4, 2), // 8
    basis( 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 3, 4, 3, 2), // 9
    basis( 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0), // 10
    basis( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0), // 11
};

} // namespace detail

} // namespace loop

} // namespace happah
