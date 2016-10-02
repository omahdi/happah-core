#pragma once

#include <array>

#include "happah/math/BezierPatch.h"
#include "happah/geometries/TriangleMesh.h"

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

template<class Vertex>
static LoopPatch<Vertex> extract_regular_patch(const TriangleMesh<Vertex, Format::DIRECTED_EDGE>& mesh, hpuint u, hpuint v, hpuint w) {
    // locate edge uv
    auto it = mesh.template cbegin<View::VERTEX, Mode::EDGES>(u);
    assert(mesh.getDegree(u) == 6);
    assert(mesh.getDegree(v) == 6);
    assert(mesh.getDegree(w) == 6);
    auto begin = it;
    do {
        if ((*it).first.vertex == v) break;
        ++it;
    } while(it != begin);
    assert(it != begin || (*begin).first.vertex == v); // if we are regular, we should never wrap around

    LoopPatch<Vertex> result{};


    {
        /* Index order (internal/output):
         *     7---6        10--11
         *    / \ / \       / \ / \
         *   1---4---5     7---8---9
         *  / \ / \ / \   / \ / \ / \
         * 2---0---8--11 3---4---5---6
         *  \ / \ / \ /   \ / \ / \ /
         *   3---9--10     0---1---2
         * for
         *       v
         *      / \
         *     u---w
         */
        //                        0  1  2  3  4  5  6  7  8  9 10 11
        const LoopPatch<int> perm{4, 7, 3, 0, 8, 9,11,10, 5, 1, 2, 6};
        using edge_t = decltype((*it).first);
        // navigation helpers
        // next edge along the face
        const auto f_next = [&](const edge_t& edge) { return mesh.getEdge(edge.next); };
        // opposite edge
        const auto opp    = [&](const edge_t& edge) { return mesh.getEdge(edge.opposite); };
        // previous edge around the vertex (clockwise)
        const auto v_prev = [&](const edge_t& edge) { return f_next(opp(edge)); };
        // vertex from edge
        const auto vertex = [&](const edge_t& edge) { return mesh.getVertex(edge.vertex); };
        int i = 0;
        auto edge = (*it).first;
        assert(edge.vertex == v && opp(edge).vertex == u);
        for (int x = 0; x < 3; ++x) {
            // first the origin vertex
            result[perm[i++]] = vertex(opp(edge));
            auto temp = edge;
            // then the first 3 vertices outside the triangle in counter-clockwise order
            for(int y = 0; y < 3; ++y) {
                temp = v_prev(temp);
                result[perm[i++]] = vertex(temp);
            }
            edge = f_next(edge);
        }
    }
    return result;
}

} // namespace loop

} // namespace happah
