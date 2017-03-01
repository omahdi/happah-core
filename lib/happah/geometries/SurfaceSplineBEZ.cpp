// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/SurfaceSplineBEZ.h"

namespace happah {

hpuint make_boundary_offset(hpuint degree, hpuint i, hpuint k) {
     switch(i) {
     case 0u: return k + 1u;
     case 1u: return (degree << 1) + ((k * ((degree << 1) - k - 1u)) >> 1);
     case 2u: return make_patch_size(degree) - 3u - ((k * (5u + k)) >> 1);
     }
}

hpuint make_interior_offset(hpuint degree, hpuint i) {
     auto delta = degree - 2u;
     auto end = degree - 3u;
     auto j = degree + 2u;//absolute index
     while(i > end) {
          --delta;
          end += delta;
          j += delta + 3u;
     }
     return j + i - (end - delta + 1u);
}

bool validate_projective_structure(const Indices& neighbors, const std::vector<hpreal>& transitions, hpreal epsilon) {
     auto is_one = [&](auto a) { return glm::abs(1.0 - a) < epsilon; };
     auto is_zero = [&](auto a) { return glm::abs(a) < epsilon; };
     return !find_fan(neighbors, [&](auto p, auto i, auto fan) -> auto {
          auto A1 = hpvec3(0.0, 1.0, 0.0);
          auto A2 = hpvec3(0.0, 0.0, 1.0);
          while(fan) {
               auto q = 0u, j = 0u;
               std::tie(q, j) = *fan;
               auto transition = std::begin(transitions) + (9 * q + 3 * j);
               auto temp = A2;
               if(j == 0) {
                    A2 = transition[0] * A1 + transition[2] * A2;
                    A2.x += transition[1];
               } else if(j == 1) {
                    A2 = transition[0] * A2 + transition[1] * A1;
                    A2.x += transition[2];
               } else {
                    A2 = transition[1] * A2 + transition[2] * A1;
                    A2.x += transition[0];
               }
               A1 = temp;
               ++fan;
          }
          return !(is_zero(A1[0]) && is_one(A1[1]) && is_zero(A1[2]) && is_zero(A2[0]) && is_zero(A2[1]) && is_one(A2[2]));
     });
}

namespace phm {

//Returns cubics of the form $$ax_jx_kx_l+bx_j'x_k'+cx_j''+d=0$$, where the cubic coefficients are stored in the first vector, the quadratic coefficients in the second, the linear coefficients in the third, and the constant coefficients in the fourth.  The first integer (the 'i') in each entry of the four vectors identifies the constraint.
std::tuple<std::vector<hpijklr>, std::vector<hpijkr>, std::vector<hpijr>, std::vector<hpir> > make_constraints(const Indices& neighbors) {
     //NOTE: There are 36p variables and 27p constraints, where p is the number of patches.
     auto nPatches = neighbors.size() / 3;
     auto hpirs = std::vector<hpir>();
     auto hpijrs = std::vector<hpijr>();
     auto hpijkrs = std::vector<hpijkr>();
     auto hpijklrs = std::vector<hpijklr>();
     auto row = -1;

     //TODO: reserve hpi...s

     // indexing of rho: (see make_objective)
     // indexing of lambda: (see make_objective)

     visit_edges(neighbors, [&](auto p, auto i) {
          auto q = make_neighbor_index(neighbors, p, i);
          auto j = make_neighbor_offset(neighbors, q, p);
          auto op = 36 * p + 12 * i;
          auto oq = 36 * q + 12 * j;
          auto x = std::array<hpuint, 12>();
          auto y = std::array<hpuint, 12>();

          auto do_row = [&](auto x0, auto x1, auto x2, auto y0, auto y3, auto y6) {
               hpijkrs.emplace_back(++row, y0, x1, 1.0);
               hpijkrs.emplace_back(row, y3, x0, 1.0);
               hpijklrs.emplace_back(row, y6, y[11], x0, 1.0);
               hpijklrs.emplace_back(row, y6, y[9], x1, 1.0);
               hpijklrs.emplace_back(row, y6, y[10], x2, 1.0);
          };

          auto do_column = [&](auto y0, auto y3, auto y6) {
               do_row(x[0], x[1], x[2], y0, y3, y6);
               do_row(x[3], x[4], x[5], y0, y3, y6);
               do_row(x[6], x[7], x[8], y0, y3, y6);
          };

          std::iota(std::begin(x), std::end(x), op);
          std::iota(std::begin(y), std::end(y), oq);

          // lambda * lambda' = id
          hpijkrs.emplace_back(++row, x[11], y[10], 1.0);
          hpijrs.emplace_back(row, y[9], 1.0);
          hpijkrs.emplace_back(++row, x[9], y[10], 1.0);
          hpijrs.emplace_back(row, y[11], 1.0);
          hpijkrs.emplace_back(++row, x[10], y[10], 1.0);
          hpirs.emplace_back(row, -1.0);

          // rho * lambda' * rho' = lambda'
          hpijrs.emplace_back(row + 1, y[9], -1.0);
          hpijrs.emplace_back(row + 2, y[10], -1.0);
          hpijrs.emplace_back(row + 3, y[11], -1.0);
          hpirs.emplace_back(row + 4, -1.0);
          hpirs.emplace_back(row + 9, -1.0);
          do_column(y[0], y[3], y[6]);
          do_column(y[1], y[4], y[7]);
          do_column(y[2], y[5], y[8]);
     });

     auto make_array = [](hpuint op) -> auto { return std::array<hpuint, 9>{ op + 6u, op + 7u, op + 8u, op, op + 1u, op + 2u, op + 3u, op + 4u, op + 5u }; };

     // rho0 * rho1 * rho2 = id
     for(auto p : boost::irange(0lu, nPatches)) {
          auto op = 36 * p;
          auto x = make_array(op);
          auto y = make_array(op + 12);
          auto z = make_array(op + 24);

          auto do_column_x = [&](auto z0, auto y0, auto x0, auto x3, auto x6) {
               hpijklrs.emplace_back(row + 1, z0, y0, x0, 1.0);
               hpijklrs.emplace_back(row + 2, z0, y0, x3, 1.0);
               hpijklrs.emplace_back(row + 3, z0, y0, x6, 1.0);
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
          hpirs.emplace_back(row + 1, -1);
          row += 3;
          do_column_z(z[1], z[4], z[7]);
          hpirs.emplace_back(row + 2, -1);
          row += 3;
          do_column_z(z[2], z[5], z[8]);
          hpirs.emplace_back(row + 3, -1);
          row += 3;
     }

     return std::make_tuple(std::move(hpijklrs), std::move(hpijkrs), std::move(hpijrs), std::move(hpirs));
}

}//namespace phm

}//namespace happah

