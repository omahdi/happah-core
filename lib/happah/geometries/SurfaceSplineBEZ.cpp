// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/SurfaceSplineBEZ.h"

namespace happah {

namespace ssb {

RingEnumerator::RingEnumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i)
     : m_degree(degree), m_e(neighbors, p, i), m_flag(true) {}

RingEnumerator::operator bool() { return m_flag; }

std::tuple<hpuint, hpuint> RingEnumerator::operator*() {
     auto p = 0u, i = 0u;
     std::tie(p, i) = *m_e;
     if(m_e) return std::make_tuple(p, (i == 0) ? 1 : (i == 1) ? (m_degree << 1) : (make_patch_size(m_degree) - 3));
     else return std::make_tuple(m_p, (m_i == 0) ? (m_degree + 1) : (m_i == 1) ? (m_degree - 1) : (make_patch_size(m_degree) - 2));
}

RingEnumerator& RingEnumerator::operator++() {
     std::tie(m_p, m_i) = *m_e;
     m_flag = !!m_e;
     ++m_e;
     m_flag &= !!m_e || (m_flag && !m_e && std::get<0>(*m_e) == UNULL);
     return *this;
}

}//namespace ssb

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

ssb::RingEnumerator make_ring_enumerator(hpuint degree, const Indices& neighbors, hpuint p, hpuint i) { return { degree, neighbors, p, i }; }

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

}//namespace happah

