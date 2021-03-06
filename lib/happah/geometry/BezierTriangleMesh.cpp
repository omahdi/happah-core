// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.08 - Hedwig Amberg    - Added NablasEnumerator and paint_edges function.

#include "happah/geometry/BezierTriangleMesh.hpp"
#include "happah/geometry/Triangle.hpp"
#include "happah/math/functions.hpp"

namespace happah {

std::vector<hpreal> make_de_casteljau_matrix(hpuint degree, hpuint nSamples) {
     std::vector<hpreal> matrix;
     sample(nSamples, [&](hpreal u, hpreal v, hpreal w) {
          hpuint coefficient = 1;
          hpuint k = 0;
          hpreal wk = 1.0;
          while(k < degree) {
               hpuint limit = degree - k;
               hpuint i = limit;
               hpreal ui = pow(u, i);
               matrix.push_back(coefficient * ui * wk);//j=0
               coefficient *= i;
               --i;
               ui = pow(u, i);//NOTE: Here we could also do 'ui /= u' but to avoid division by zero, we recalculate u^i.
               hpuint j = 1;
               hpreal vj = v;
               while(j < limit) {
                    matrix.push_back(coefficient * ui * vj * wk);
                    coefficient *= i;
                    --i;
                    ui = pow(u, i);
                    ++j;
                    coefficient /= j;
                    vj *= v;
               }
               matrix.push_back(coefficient * vj * wk);//i=0
               coefficient *= limit;
               wk *= w;
               ++k;
               coefficient /= k;
          }
          matrix.push_back(wk);//k=degree
     });
     return std::move(matrix);
}

namespace mdz {

std::tuple<std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const Triples<trix>& neighbors) {
     auto nPatches = hpuint(neighbors.size() / 3);
     auto irs = std::vector<hpir>();
     auto ijrs = std::vector<hpijr>();
     auto ijkrs = std::vector<hpijkr>();
     auto row = -1;

     // indexing of rho: (see make_objective)

     auto insert = [&](auto x0, auto x3, auto x6, auto offset) {
          ijkrs.emplace_back(++row, offset + 0, x0, 1.0);
          ijkrs.emplace_back(  row, offset + 1, x3, 1.0);
          ijkrs.emplace_back(  row, offset + 2, x6, 1.0);
          ijkrs.emplace_back(++row, offset + 3, x0, 1.0);
          ijkrs.emplace_back(  row, offset + 4, x3, 1.0);
          ijkrs.emplace_back(  row, offset + 5, x6, 1.0);
          ijkrs.emplace_back(++row, offset + 6, x0, 1.0);
          ijkrs.emplace_back(  row, offset + 7, x3, 1.0);
          ijkrs.emplace_back(  row, offset + 8, x6, 1.0);
     };

     // pi * qj = id
     visit(make_edges_enumerator(neighbors), [&](auto x) {
          auto op = 9 * x;
          auto oq = 9 * neighbors[x];

          irs.emplace_back(row + 1, -1.0);
          insert(oq + 0, oq + 3, oq + 6, op);
          irs.emplace_back(row + 2, -1.0);
          insert(oq + 1, oq + 4, oq + 7, op);
          irs.emplace_back(row + 3, -1.0);
          insert(oq + 2, oq + 5, oq + 8, op);
     });

     // p0 * p1 = qj
     for(auto p : boost::irange(hpuint(0), nPatches)) {
          auto x = neighbors(p, TRIT2);
          auto q = x.getTriple();
          auto j = x.getOffset();
          auto op0 = 27 * p;
          auto op1 = 27 * p + 9;
          auto oq = 27 * q + 9 * j;

          auto do_column = [&](auto o0, auto o3, auto o6) {
               ijrs.emplace_back(row + 1, oq + o0, -1.0);
               ijrs.emplace_back(row + 2, oq + o3, -1.0);
               ijrs.emplace_back(row + 3, oq + o6, -1.0);
               insert(op1 + o0, op1 + o3, op1 + o6, op0);
          };

          do_column(0, 3, 6);
          do_column(1, 4, 7);
          do_column(2, 5, 8);
     }

     return std::make_tuple(std::move(ijkrs), std::move(ijrs), std::move(irs));
}

}//namespace mdz

namespace phm {

std::tuple<std::vector<hpijklr>, std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const Triples<trix>& neighbors) {
     //NOTE: There are 36p variables and 27p constraints, where p is the number of patches.
     auto nPatches = hpuint(neighbors.size() / 3);
     auto irs = std::vector<hpir>();
     auto ijrs = std::vector<hpijr>();
     auto ijkrs = std::vector<hpijkr>();
     auto ijklrs = std::vector<hpijklr>();
     auto row = -1;

     //TODO: reserve hpi...s

     // indexing of rho: (see make_objective)
     // indexing of lambda: (see make_objective)

     visit(make_edges_enumerator(neighbors), [&](auto e) {
          auto op = 12 * e;
          auto oq = 12 * neighbors[e];
          auto x = std::array<hpuint, 12>();
          auto y = std::array<hpuint, 12>();

          auto do_row = [&](auto x0, auto x1, auto x2, auto y0, auto y3, auto y6) {
               ijkrs.emplace_back(++row, y0, x1, 1.0);
               ijkrs.emplace_back(row, y3, x0, 1.0);
               ijklrs.emplace_back(row, y6, y[11], x0, 1.0);
               ijklrs.emplace_back(row, y6, y[9], x1, 1.0);
               ijklrs.emplace_back(row, y6, y[10], x2, 1.0);
          };

          auto do_column = [&](auto y0, auto y3, auto y6) {
               do_row(x[0], x[1], x[2], y0, y3, y6);
               do_row(x[3], x[4], x[5], y0, y3, y6);
               do_row(x[6], x[7], x[8], y0, y3, y6);
          };

          std::iota(std::begin(x), std::end(x), op);
          std::iota(std::begin(y), std::end(y), oq);

          // lambda * lambda' = id
          ijkrs.emplace_back(++row, x[11], y[10], 1.0);
          ijrs.emplace_back(row, y[9], 1.0);
          ijkrs.emplace_back(++row, x[9], y[10], 1.0);
          ijrs.emplace_back(row, y[11], 1.0);
          ijkrs.emplace_back(++row, x[10], y[10], 1.0);
          irs.emplace_back(row, -1.0);

          // rho * lambda' * rho' = lambda'
          irs.emplace_back(row + 2, -1.0);
          irs.emplace_back(row + 4, -1.0);
          ijrs.emplace_back(row + 7, y[11], -1.0);
          ijrs.emplace_back(row + 8, y[9], -1.0);
          ijrs.emplace_back(row + 9, y[10], -1.0);
          do_column(y[0], y[3], y[6]);
          do_column(y[1], y[4], y[7]);
          do_column(y[2], y[5], y[8]);
     });

     auto make_array = [](hpuint op) -> auto { return std::array<hpuint, 9>{ op + 6u, op + 7u, op + 8u, op, op + 1u, op + 2u, op + 3u, op + 4u, op + 5u }; };

     // rho0 * rho1 * rho2 = id
     for(auto p : boost::irange(hpuint(0), nPatches)) {
          auto op = 36 * p;
          auto x = make_array(op);
          auto y = make_array(op + 12);
          auto z = make_array(op + 24);

          auto do_column_x = [&](auto z0, auto y0, auto x0, auto x3, auto x6) {
               ijklrs.emplace_back(row + 1, z0, y0, x0, 1.0);
               ijklrs.emplace_back(row + 2, z0, y0, x3, 1.0);
               ijklrs.emplace_back(row + 3, z0, y0, x6, 1.0);
          };

          auto do_column_y = [&](auto z0, auto y0, auto y3, auto y6) {
               do_column_x(z0, y0, x[0], x[3], x[6]);
               do_column_x(z0, y3, x[1], x[4], x[7]);
               do_column_x(z0, y6, x[2], x[5], x[8]);
          };
          
          auto do_column_z = [&](auto z0, auto z3, auto z6) {
               do_column_y(z0, y[0], y[3], y[6]);
               do_column_y(z3, y[1], y[4], y[7]);
               do_column_y(z6, y[2], y[5], y[8]);
          };

          do_column_z(z[0], z[3], z[6]);
          irs.emplace_back(row + 1, -1);
          row += 3;
          do_column_z(z[1], z[4], z[7]);
          irs.emplace_back(row + 2, -1);
          row += 3;
          do_column_z(z[2], z[5], z[8]);
          irs.emplace_back(row + 3, -1);
          row += 3;
     }

     return std::make_tuple(std::move(ijklrs), std::move(ijkrs), std::move(ijrs), std::move(irs));
}

std::vector<hpreal> make_transitions(const std::vector<hpreal>& solution) {
     auto transitions = std::vector<hpreal>();

     transitions.reserve(solution.size() << 2);
     for(auto i = std::begin(solution) + 9, end = std::end(solution) + 9; i != end; i += 12) transitions.insert(transitions.end(), i, i + 3);

     return transitions;
}

}//namespace phm

}//namespace happah

