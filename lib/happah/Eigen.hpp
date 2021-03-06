// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#if !defined(__clang__)
#include <quadmath.h>
#endif
#include <vector>

namespace std {

__float128 sqrt(__float128 f);
__float128 abs(__float128 f);

}//namespace std

#include <Eigen/Dense>
#include <Eigen/Sparse>

template<class T>
Eigen::SparseMatrix<T> make_sparse_matrix(int m, int n, const std::vector<Eigen::Triplet<T> >& triplets) {
     auto matrix = Eigen::SparseMatrix<T>(m, n);

     matrix.setFromTriplets(std::begin(triplets), std::end(triplets));
     matrix.makeCompressed();

     return matrix;
}

namespace lsq {

// Minimize x - a with constraints Bx = b.
// The solution is a - B^t(BB^t)^{-1}(Ba - b).
template<class T>
auto solve(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a, const std::vector<Eigen::Triplet<T> >& B, const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) {
     auto solver = Eigen::SparseLU<Eigen::SparseMatrix<T> >();
     auto temp = make_sparse_matrix(B);
     auto S = temp * temp.transpose();
     auto eye = Eigen::SparseMatrix<T>(S.size());

     solver.compute(S);
     assert(solver.info() == Eigen::Success);
     eye.setIdentity();

     return a - temp.transpose() * solver.solve(eye) * (temp * a - b);
}

// Minimize Ax - b.
template<class T>
auto solve(const Eigen::SparseMatrix<T>& A, const Eigen::Matrix<T, Eigen::Dynamic, 1>& b) {
     Eigen::SparseQR<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int> > solver;
     solver.analyzePattern(A);
     solver.factorize(A);
     return Eigen::Matrix<T, Eigen::Dynamic, 1>(solver.solve(b));
}

}//namespace lsq

