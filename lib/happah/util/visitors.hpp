// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

template<class Iterator>
auto get_triplet(Iterator begin, hpuint i, hpuint stride) {
     begin += stride * i;
     return std::tie(begin[0], begin[1], begin[2]);
}

template<class T>
auto get_triplet(const std::vector<T>& ts, hpuint i) { return get_triplet(std::begin(ts), i, 3); }

template<class Iterator, class Visitor>
void visit_pairs(Iterator begin, hpuint n, hpuint stride, Visitor&& visit) {
     repeat(n, [&]() {
          visit(begin[0], begin[1]);
          begin += stride;
     });
}

template<class T, class Visitor>
void visit_pairs(const std::vector<T>& ts, Visitor&& visit) { visit_pairs(std::begin(ts), ts.size() / 2, 2, std::forward<Visitor>(visit)); }

template<class Iterator, class Visitor>
void visit_quartets(Iterator begin, hpuint n, hpuint stride, Visitor&& visit) {
     repeat(n, [&]() {
          visit(begin[0], begin[1], begin[2], begin[3]);
          begin += stride;
     });
}

template<class T, class Visitor>
void visit_quartets(const std::vector<T>& ts, Visitor&& visit) { visit_quartets(std::begin(ts), ts.size() / 4, 4, std::forward<Visitor>(visit)); }


template<class Visitor>
void visit_rings(const Indices& offsets, const Indices& indices, Visitor&& visit) {
     auto i = std::begin(indices);
     for(auto o = std::begin(offsets), end = std::end(offsets) - 1; o != end; ++o) visit(i + o[0], i + o[1]);
}

template<class Visitor>
void visit_rings(const Indices& centers, const Indices& offsets, const Indices& indices, Visitor&& visit) {
     auto o = std::begin(offsets);
     auto i = std::begin(indices);
     for(auto center : centers) {
          visit(center, i + o[0], i + o[1]);
          ++o;
     }
}

template<class Iterator, class Visitor>
void visit_triplets(Iterator begin, hpuint n, hpuint stride, Visitor&& visit) {
     repeat(n, [&]() {
          visit(begin[0], begin[1], begin[2]);
          begin += stride;
     });
}

template<class T, class Visitor>
void visit_triplets(const std::vector<T>& ts, Visitor&& visit) { visit_triplets(std::begin(ts), ts.size() / 3, 3, std::forward<Visitor>(visit)); }

template<class Iterator, class Visitor>
void visit_triplet(Iterator begin, hpuint i, hpuint stride, Visitor&& visit) {
     auto temp = get_triplet(begin, i, stride);
     visit(std::get<0>(temp), std::get<1>(temp), std::get<2>(temp));
}

template<class T, class Visitor>
void visit_triplet(const std::vector<T>& ts, hpuint i, Visitor&& visit) { visit_triplet(std::begin(ts), i, 3, std::forward<Visitor>(visit)); }

}//namespace happah

