// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.h"

namespace happah {

template<class Iterator, class Visitor>
void visit_pairs(Iterator begin, hpuint n, hpuint stride, Visitor&& visit) {
     while(n > 0) {
          visit(*begin, *(begin + 1));
          begin += stride;
          --n;
     }
}

template<class T, class Visitor>
void visit_pairs(const std::vector<T>& ts, Visitor&& visit) { visit_pairs(ts.begin(), ts.size() / 2, 2, std::forward<Visitor>(visit)); }

template<class Visitor>
void visit_rings(const Indices& offsets, const Indices& indices, Visitor&& visit) {
     auto o = offsets.begin();
     auto i = indices.begin();
     for(auto o = offsets.begin(), end = offsets.end(); o + 1 != end; ++o) visit(i + *o, i + *(o + 1));
}

template<class Visitor>
void visit_rings(const Indices& centers, const Indices& offsets, const Indices& indices, Visitor&& visit) {
     auto o = offsets.begin();
     auto i = indices.begin();
     for(auto center : centers) {
          auto ob = o;
          ++o;
          visit(center, i + *ob, i + *o);
     }
}

template<class Iterator, class Visitor>
void visit_triplets(Iterator begin, hpuint n, hpuint stride, Visitor&& visit) {
     while(n > 0) {
          visit(*begin, *(begin + 1), *(begin + 2));
          begin += stride;
          --n;
     }
}

template<class T, class Visitor>
void visit_triplets(const std::vector<T>& ts, Visitor&& visit) { visit_triplets(ts.begin(), ts.size() / 3, 3, std::forward<Visitor>(visit)); }

template<class Iterator, class Visitor>
void visit_triplet(Iterator begin, hpuint i, hpuint stride, Visitor&& visit) {
     begin += stride * i;
     visit(*begin, *(begin + 1), *(begin + 2));
}

template<class T, class Visitor>
void visit_triplet(const std::vector<T>& ts, hpuint i, Visitor&& visit) { visit_triplet(ts.begin(), i, 3, std::forward<Visitor>(visit)); }

}//namespace happah

