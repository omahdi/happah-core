// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/utils/ControlPointIndexer.hpp"

namespace happah {

const hpuint CubicControlPointIndexerBEZ::ROW0[] = {3, 2, 1, 0, 9, 8, 6, 3, 0, 4, 7, 9, 0, 1, 2, 3, 3, 6, 8, 9, 9, 7, 4, 0};
const hpuint CubicControlPointIndexerBEZ::ROW1[] = {6, 5, 4, 7, 5, 2, 1, 5, 8, 4, 5, 6, 2, 5, 7, 8, 5, 1};
const hpuint CubicControlPointIndexerBEZ::ROW2[] = {8, 7, 4, 1, 2, 6, 7, 8, 1, 4, 6, 2};
const hpuint CubicControlPointIndexerBEZ::ROW3[] = {9, 0, 3};

const hpuint LinearControlPointIndexerBEZ::ROW0[] = {1, 0, 2, 1, 0, 2, 0, 1, 1, 2, 2, 0};
const hpuint LinearControlPointIndexerBEZ::ROW1[] = {2, 0, 1};

const hpuint QuadraticControlPointIndexerBEZ::ROW0[] = {2, 1, 0, 5, 4, 2, 0, 3, 5, 0, 1, 2, 2, 4, 5, 5, 3, 0};
const hpuint QuadraticControlPointIndexerBEZ::ROW1[] = {4, 3, 3, 1, 1, 4, 3, 4, 1, 3, 4, 1};
const hpuint QuadraticControlPointIndexerBEZ::ROW2[] = {5, 0, 2};

const hpuint QuarticControlPointIndexerBEZ::ROW0[] = {4, 3, 2, 1, 0, 14, 13, 11, 8, 4, 0, 5, 9, 12, 14, 0, 1, 2, 3, 4, 4, 8, 11, 13, 14, 14, 12, 9, 5, 0};
const hpuint QuarticControlPointIndexerBEZ::ROW1[] = {8, 7, 6, 5, 12, 10, 7, 3, 1, 6, 10, 13, 5, 6, 7, 8, 3, 7, 10, 12, 13, 10, 6, 1};
const hpuint QuarticControlPointIndexerBEZ::ROW2[] = {11, 10, 9, 9, 6, 2, 2, 7, 11, 9, 10, 11, 2, 6, 9, 11, 7, 2};
const hpuint QuarticControlPointIndexerBEZ::ROW3[] = {13, 12, 5, 1, 3, 8, 12, 13, 1, 5, 8, 3};
const hpuint QuarticControlPointIndexerBEZ::ROW4[] = {14, 0, 4};

}//namespace happah

