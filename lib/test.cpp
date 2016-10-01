#include <iostream>
#include <iomanip>
#include <algorithm>

#include "happah/math/BezierPatch.h"
#include "happah/math/BoxSpline.h"
#include "happah/math/MathUtils.h"

void assert_equal(double expected, double real) {
    if (std::abs(expected - real) > happah::EPSILON) {
        std::cout << "Expected " << expected << ", but got " << real << std::endl;
    }
}

void test_linear_precision() {
    auto coords = happah::bezier_patch_ctrl_points<4>();
    happah::QuarticBezierPatch<double> patch{};
    // simple check: linear precision
    for (int i = 0; i < 15; ++i) {
        patch[i] = (coords[i].r + coords[i].s) / 4.0;
    }
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            double u = i / 10.0;
            double v = j / 10.0;
            double w = 1 - u - v;
            double expected = (u + v);
            double real = happah::bezier_patch_evaluate<double, 4, double>(patch, u, v, w);
            assert_equal(expected, real);
        }
    }
}

void test_linear_precision_loop() {
    /* input values:
     *     5---7
     *    / \ / \
     *   3---5---7
     *  / \ / \ / \
     * 1---3---5---7
     *  \ / \ / \ /
     *   1---3---5
     * -> f(x,y) = x + y;
     * Triangle: (2, 1) (4, 1) (3, 2)
     */
    happah::loop::LoopPatch<double> patch = {1, 3, 5,
                                       1, 3, 5, 7,
                                          3, 5, 7,
                                             5, 7};
    auto bezier = happah::loop::loop_to_bezier(patch);
    for (int i = 0; i <= 10; ++i) {
        for (int j = 0; j <= 10; ++j) {
            double u = i / 10.0;
            double v = j / 10.0;
            double w = 1 - u - v;
            double x = 2 * u + 3 * v + 4 * w;
            double y = 1 * u + 2 * v + 1 * w;
            double expected = (x + y);
            double real = happah::bezier_patch_evaluate<double, 4, double>(bezier, u, v, w);
            assert_equal(expected, real);
        }
    }
}

int main() {
    std::cout << "Bézier linear precision" << std::endl;
    test_linear_precision();
    std::cout << "Box spline linear precision" << std::endl;
    test_linear_precision_loop();

}
