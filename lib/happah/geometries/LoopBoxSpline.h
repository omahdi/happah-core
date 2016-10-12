// Copyright 2015-2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <glm/glm.hpp>
#include <vector>

#include "TriangleMesh.h"

class LoopBoxSpline {
public:
     //NOTE: Merge the rings of three vertices of valence six into a regular patch.
     template<class Iterator0, class Iterator1, class Iterator2>
     static std::array<unsigned int, 12> merge(unsigned int v0, Iterator0 begin0, Iterator0 end0, unsigned int v1, Iterator1 begin1, Iterator1 end1, unsigned int v2, Iterator2 begin2, Iterator2 end2) {
          auto i0 = begin0;
          while(*i0 != v2) ++i0;
          ++i0;
          if(i0 == end0) i0 = begin0;
          auto a7 = *i0; 
          ++i0;
          if(i0 == end0) i0 = begin0;
          auto a3 = *i0;
          ++i0;
          if(i0 == end0) i0 = begin0;
          auto a0 = *i0;
          ++i0;
          if(i0 == end0) i0 = begin0;
          auto a1 = *i0;
          auto i1 = begin1;
          while(*i1 != a1) ++i1;
          ++i1;
          if(i1 == end1) i1 = begin1;
          auto a2 = *i1;
          ++i1;
          if(i1 == end1) i1 = begin1;
          auto a6 = *i1;
          ++i1;
          if(i1 == end1) i1 = begin1;
          auto a9 = *i1;
          auto i2 = begin2;
          while(*i2 != a9) ++i2;
          ++i2;
          if(i2 == end2) i2 = begin2;
          auto a11 = *i2;
          ++i2;
          if(i2 == end2) i2 = begin2;
          auto a10 = *i2;
          return {a0, a1, a2, a3, v0, v1, a6, a7, v2, a9, a10, a11};
     }

     LoopBoxSpline(std::array<glm::vec3, 12> controlPoints);

     const std::array<glm::vec3, 12>& getControlPoints() const;

     //    10  11
     //   7   8   9
     // 3   4   5   6
     //   0   1   2
     TriangleMesh getControlPolygon() const;

private:
     std::array<glm::vec3, 12> m_controlPoints;

};

template<class Iterator, class Visitor>
void visit_pairs(Iterator begin, Iterator end, Visitor visit) {
     auto i = begin;
     auto i0 = i;
     ++i;
     do {
          visit(*i0, *i);
          i0 = i;
          ++i;
     } while(i != end);
     visit(*i0, *begin);
}

template<class Visitor>
void visit_rings(const Indices& offsets, const Indices& indices, Visitor visit) {
     auto o = offsets.begin();
     auto i = indices.begin();
     for(auto o = offsets.begin(), end = offsets.end(); o + 1 != end; ++o) visit(i + *o, i + *(o + 1));
}

template<class Visitor>
void visit_rings(const Indices& centers, const Indices& offsets, const Indices& indices, Visitor visit) {
     auto o = offsets.begin();
     auto i = indices.begin();
     for(auto center : centers) {
          auto ob = o;
          ++o;
          visit(center, i + *ob, i + *o);
     }
}

template<class Iterator, class Visitor>
void visit_triplets(Iterator begin, Iterator end, Visitor visit) {
     auto i = begin;
     auto i0 = i;
     auto i1 = ++i;
     ++i;
     do {
          visit(*i0, *i1, *i);
          i0 = i1;
          i1 = i;
          ++i;
     } while(i != end);
     //if(std::distance(begin, end) == 3) return; TODO
     visit(*i0, *i1, *begin);
     i0 = i1;
     i1 = begin;
     visit(*i0, *i1, *(++begin));
}

