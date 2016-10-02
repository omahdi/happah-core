#include <iostream>
#include <iomanip>
#include <algorithm>

#include "happah/math/BezierPatch.h"
#include "happah/math/BoxSpline.h"
#include "happah/math/MathUtils.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/io/writers/WriterOFF.h"

using namespace happah;

/*TriangleMesh<VertexP3N, Format::DIRECTED_EDGE>*/ auto build_simple_torus() {
    std::vector<VertexP3N> vertices{};
    std::vector<hpuint> indices{};
    auto insert_verts = [&](int x, int y) {
        vertices.emplace_back(glm::vec3(x, y, 1));
        vertices.emplace_back(glm::vec3(x, y, 0));
        vertices.emplace_back(2.f * glm::vec3(x, y, 0));
        vertices.emplace_back(2.f * glm::vec3(x, y, 0.5f));
    };
    insert_verts(1, 0);
    insert_verts(0, 1);
    insert_verts(-1, 0);
    insert_verts(0, -1);

    auto insert_faces = [&](hpuint a, hpuint b, hpuint c, hpuint d) {
        indices.insert(indices.end(), {a, b, c, a, c, d});
    };
    auto index = [](hpuint x, hpuint y) {
        return x + 4 * y;
    };
    for (hpuint i = 0; i < 4; ++i) {
        for (hpuint j = 0; j < 4; ++j) {
            int x = (i + 1) % 4;
            int y = (j + 1) % 4;
            insert_faces(index(i, j), index(i, y), index(x, y), index(x, j));
        }
    }
    return TriangleMesh<VertexP3N, Format::DIRECTED_EDGE>(vertices, indices);
}

void assert_equal(double expected, double real) {
    if (std::abs(expected - real) > happah::EPSILON) {
        std::cout << "Expected " << expected << ", but got " << real << std::endl;
    }
}

void test_linear_precision() {
    auto coords = happah::bezier_patch_ctrl_points<4>();
    QuarticBezierPatch<double> patch{};
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
            double real = bezier_patch_evaluate<double, 4, double>(patch, u, v, w);
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
    loop::LoopPatch<double> patch = {1, 3, 5,
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
            double real = bezier_patch_evaluate<double, 4, double>(bezier, u, v, w);
            assert_equal(expected, real);
        }
    }
}

VertexP3N& operator+=(VertexP3N& fst, VertexP3N snd) {
    fst.position += snd.position;
    return fst;
}

VertexP3N operator*(VertexP3N fst, double snd) {
    return VertexP3N(fst.position * float(snd));
}

VertexP3N operator/(VertexP3N fst, double snd) {
    return VertexP3N(fst.position / float(snd));
}

int main() {
    std::cout << "Bézier linear precision" << std::endl;
    test_linear_precision();
    std::cout << "Box spline linear precision" << std::endl;
    test_linear_precision_loop();
    auto mesh = build_simple_torus();
    // TODO so far, WriterOFF doesn't seem to work with Format::...
    // WriterOFF::write(mesh, "out.off");
    std::vector<VertexP3N> vertices{};
    std::vector<hpuint> out_indices{};
    int i = 0;
    const auto& indices = mesh.getIndices();
    for (int k = 0; k < indices.size(); k += 3) {
        auto patch = loop::extract_regular_patch(mesh, indices[k], indices[k + 1], indices[k + 2]);
        auto bezier = loop::loop_to_bezier(patch);
        const double steps = 10;
        // TODO extract method for later use
        for (int t = 0; t <= steps; ++t) {
            for (int s = 0; s <= steps - t; ++s) {
                int r = steps - s - t;
                double u = r / steps;
                double v = s / steps;
                double w = t / steps;

                vertices.emplace_back(bezier_patch_evaluate<VertexP3N, 4, double>(bezier, u, v, w));
                if (s > 0) {
                    out_indices.insert(out_indices.end(), {i - 1, i, i + r + s});
                    if (t > 0) {
                        out_indices.insert(out_indices.end(), {i, i - 1, i - (r + s + 2)});
                    }
                }
                ++i;
            }
        }
    }
    WriterOFF::write(TriangleMesh3D(vertices, out_indices), "out2.off");
}
