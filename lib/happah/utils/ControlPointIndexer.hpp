// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"

namespace happah {

namespace prv {

template<hpuint t_degree>
class ControlPointIndexerBaseBEZ {
public:
     static constexpr hpuint CLOCKWISE = 1;
     static constexpr hpuint COUNTER_CLOCKWISE = 2;

protected:
     template<hpuint t_direction, hpuint t_edge, hpuint t_row>
     static void test() {
          static_assert(t_direction == CLOCKWISE || t_direction == COUNTER_CLOCKWISE, "The direction has to either be clockwise or counter-clockwise.");
          static_assert(t_edge < 3, "A parameter triangle has at most three edges.");
          static_assert(t_row < (t_degree + 1), "There are only (degree+1) rows in a Bezier configuration.");
     }

     template<hpuint t_direction, hpuint t_edge>
     static void test(hpuint row) {
          static_assert(t_direction == CLOCKWISE || t_direction == COUNTER_CLOCKWISE, "The direction has to either be clockwise or counter-clockwise.");
          static_assert(t_edge < 3, "A parameter triangle has at most three edges.");
          assert(row < (t_degree + 1));
     }

};//ControlPointIndexerBaseBEZ

}//namespace prv

template<hpuint t_degree>
class ControlPointIndexerBEZ;//TODO: general case
using CubicControlPointIndexerBEZ = ControlPointIndexerBEZ<3>;
using LinearControlPointIndexerBEZ = ControlPointIndexerBEZ<1>;
using QuadraticControlPointIndexerBEZ = ControlPointIndexerBEZ<2>;
using QuarticControlPointIndexerBEZ = ControlPointIndexerBEZ<4>;

template<>
class ControlPointIndexerBEZ<3> : public prv::ControlPointIndexerBaseBEZ<3> {
public:
     typedef const hpuint* Iterator;

     template<hpuint t_direction, hpuint t_edge, hpuint t_row>
     static Iterator getIndices() {
          test<t_direction, t_edge, t_row>();
          if(t_direction == CLOCKWISE) {
               switch(t_row) {
               case 0: return ROW0 + (t_edge << 2);
               case 1: return ROW1 + (3 * t_edge);
               case 2: return ROW2 + (t_edge << 1);
               case 3: return ROW3 + t_edge;
               };
          } else {
               switch(t_row) {
               case 0: return ROW0 + 12 + (t_edge << 2);
               case 1: return ROW1 + 9 + (3 * t_edge);
               case 2: return ROW2 + 6 + (t_edge << 1);
               case 3: return ROW3 + t_edge;
               };
          }
     }

     template<hpuint t_direction, hpuint t_edge>
     static Iterator getIndices(hpuint row) {
          test<t_direction, t_edge>(row);
          if(t_direction == CLOCKWISE) {
               switch(row) {
               case 0: return ROW0 + (t_edge << 2);
               case 1: return ROW1 + (3 * t_edge);
               case 2: return ROW2 + (t_edge << 1);
               case 3: return ROW3 + t_edge;
               };
          } else {
               switch(row) {
               case 0: return ROW0 + 12 + (t_edge << 2);
               case 1: return ROW1 + 9 + (3 * t_edge);
               case 2: return ROW2 + 6 + (t_edge << 1);
               case 3: return ROW3 + t_edge;
               };
          }
     }

private:
     static const hpuint ROW0[24];
     static const hpuint ROW1[18];
     static const hpuint ROW2[12];
     static const hpuint ROW3[3];

};//CubicControlPointIndexerBEZ

template<>
class ControlPointIndexerBEZ<1> : public prv::ControlPointIndexerBaseBEZ<1> {
public:
     typedef const hpuint* Iterator;

     template<hpuint t_direction, hpuint t_edge, hpuint t_row>
     static Iterator getIndices() {
          test<t_direction, t_edge, t_row>();
          if(t_direction == CLOCKWISE) {
               switch(t_row) {
               case 0: return ROW0 + (t_edge << 1);
               case 1: return ROW1 + t_edge;
               };
          } else {
               switch(t_row) {
               case 0: return ROW0 + 6 + (t_edge << 1);
               case 1: return ROW1 + t_edge;
               };
          }
     }

     template<hpuint t_direction, hpuint t_edge>
     static Iterator getIndices(hpuint row) {
          test<t_direction, t_edge>(row);
          if(t_direction == CLOCKWISE) {
               switch(row) {
               case 0: return ROW0 + (t_edge << 1);
               case 1: return ROW1 + t_edge;
               };
          } else {
               switch(row) {
               case 0: return ROW0 + 6 + (t_edge << 1);
               case 1: return ROW1 + t_edge;
               };
          }
     }

private:
     static const hpuint ROW0[12];
     static const hpuint ROW1[3];

};//LinearControlPointIndexerBEZ

template<>
class ControlPointIndexerBEZ<2> : public prv::ControlPointIndexerBaseBEZ<2> {
public:
     typedef const hpuint* Iterator;

     template<hpuint t_direction, hpuint t_edge, hpuint t_row>
     static Iterator getIndices() {
          test<t_direction, t_edge, t_row>();
          if(t_direction == CLOCKWISE) {
               switch(t_row) {
               case 0: return ROW0 + (3 * t_edge);
               case 1: return ROW1 + (t_edge << 1);
               case 2: return ROW2 + t_edge;
               };
          } else {
               switch(t_row) {
               case 0: return ROW0 + 9 + (3 * t_edge);
               case 1: return ROW1 + 6 + (t_edge << 1);
               case 2: return ROW2 + t_edge;
               };
          }
     }

     template<hpuint t_direction, hpuint t_edge>
     static Iterator getIndices(hpuint row) {
          test<t_direction, t_edge>(row);
          if(t_direction == CLOCKWISE) {
               switch(row) {
               case 0: return ROW0 + (3 * t_edge);
               case 1: return ROW1 + (t_edge << 1);
               case 2: return ROW2 + t_edge;
               };
          } else {
               switch(row) {
               case 0: return ROW0 + 9 + (3 * t_edge);
               case 1: return ROW1 + 6 + (t_edge << 1);
               case 2: return ROW2 + t_edge;
               };
          }
     }

private:
     static const hpuint ROW0[18];
     static const hpuint ROW1[12];
     static const hpuint ROW2[3];

};//QuadraticControlPointIndexerBEZ

template<>
class ControlPointIndexerBEZ<4> : public prv::ControlPointIndexerBaseBEZ<4> {
public:
     using Iterator = const hpuint*;

     template<hpuint t_direction, hpuint t_edge, hpuint t_row>
     static Iterator getIndices() {
          test<t_direction, t_edge, t_row>();
          if(t_direction == CLOCKWISE) {
               switch(t_row) {
               case 0: return ROW0 + (5 * t_edge);
               case 1: return ROW1 + (t_edge << 2);
               case 2: return ROW2 + (3 * t_edge);
               case 3: return ROW3 + (t_edge << 1);
               case 4: return ROW4 + t_edge;
               };
          } else {
               switch(t_row) {
               case 0: return ROW0 + 15 + (5 * t_edge);
               case 1: return ROW1 + 12 + (t_edge << 2);
               case 2: return ROW2 + 9 + (3 * t_edge);
               case 3: return ROW3 + 6 + (t_edge << 1);
               case 4: return ROW4 + t_edge;
               };
          }
     }

     template<hpuint t_direction, hpuint t_edge>
     static Iterator getIndices(hpuint row) {
          test<t_direction, t_edge>(row);
          if(t_direction == CLOCKWISE) {
               switch(row) {
               case 0: return ROW0 + (5 * t_edge);
               case 1: return ROW1 + (t_edge << 2);
               case 2: return ROW2 + (3 * t_edge);
               case 3: return ROW3 + (t_edge << 1);
               case 4: return ROW4 + t_edge;
               };
          } else {
               switch(row) {
               case 0: return ROW0 + 15 + (5 * t_edge);
               case 1: return ROW1 + 12 + (t_edge << 2);
               case 2: return ROW2 + 9 + (3 * t_edge);
               case 3: return ROW3 + 6 + (t_edge << 1);
               case 4: return ROW4 + t_edge;
               };
          }
     }

private:
     static const hpuint ROW0[30];
     static const hpuint ROW1[24];
     static const hpuint ROW2[18];
     static const hpuint ROW3[12];
     static const hpuint ROW4[3];

};//QuarticControlPointIndexerBEZ

}//namespace happah

