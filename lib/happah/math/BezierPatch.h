#pragma once

#include <array>

#include "happah/math/MathUtils.h"

namespace happah {

constexpr int bezier_patch_num_ctrl_points(int degree) {
    return (degree + 1) * (degree + 2) / 2;
}

template <typename Vertex, int Degree>
using BezierPatch = std::array<Vertex, bezier_patch_num_ctrl_points(Degree)>;

template <typename Vertex>
using ConstantBezierPatch = BezierPatch<Vertex, 0>;

template <typename Vertex>
using LinearBezierPatch = BezierPatch<Vertex, 1>;

template <typename Vertex>
using QuadraticBezierPatch = BezierPatch<Vertex, 2>;

template <typename Vertex>
using CubicBezierPatch = BezierPatch<Vertex, 3>;

template <typename Vertex>
using QuarticBezierPatch = BezierPatch<Vertex, 4>;

struct BezierAbscissa {
    int r, s, t;
};

template <int Degree>
constexpr BezierPatch<BezierAbscissa, Degree> bezier_patch_ctrl_points() {
    BezierPatch<BezierAbscissa, Degree> ctrl_points{};
    int i = 0;
    for (int t = 0; t <= Degree; ++t) {
        for (int s = 0; s <= Degree - t; ++s) {
            int r = Degree - s - t;
            ctrl_points[i] = {r, s, t};
            ++i;
        }
    }
    return ctrl_points;
}

template <typename Vertex, int Degree, typename Scalar = float>
Vertex bezier_patch_evaluate(const BezierPatch<Vertex, Degree>& ctrl_pts, Scalar u, Scalar v, Scalar w) {
    Vertex result{};
    auto ctrl_pts_loc = bezier_patch_ctrl_points<Degree>();
    for (int i = 0; i < ctrl_pts.size(); ++i) {
        auto coord = ctrl_pts_loc[i];
        // ctrl point * (n over i,j,k) * u^i * v^j * w^t
        result += ctrl_pts[i]
                * MathUtils::munom(Degree, coord.r, coord.s)
                * MathUtils::pow(u, coord.r)
                * MathUtils::pow(v, coord.s)
                * MathUtils::pow(w, coord.t);
    }
    return result;
}
} // namespace happah
