// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

template<class Iterator, class Visitor>
void visit_pairs(Iterator begin, hpuint n, hpuint stride, Visitor&& visit) {
     repeat(n, [&]() {
          visit(begin[0], begin[1]);
          begin += stride;
     });
}

template<class T, class Visitor>
void visit_pairs(const std::vector<T>& ts, Visitor&& visit) { visit_pairs(std::begin(ts), ts.size() >> 1, 2, std::forward<Visitor>(visit)); }

template<class Iterator, class Visitor>
void visit_triples(Iterator begin, hpuint n, hpuint stride, Visitor&& visit) {
     repeat(n, [&]() {
          visit(begin[0], begin[1], begin[2]);
          begin += stride;
     });
}

template<class T, class Visitor>
void visit_triples(const std::vector<T>& ts, Visitor&& visit) { visit_triples(std::begin(ts), ts.size() / 3, 3, std::forward<Visitor>(visit)); }

}//namespace happah

