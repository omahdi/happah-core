// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <lpsolve/lp_lib.h>

#include "happah/Eigen.h"
#include "happah/Happah.h"
#include "happah/utils/ControlPointIndexer.h"
#include "happah/utils/SurfaceUtilsBEZ.h"

namespace happah {

//NOTE: C0 continuity is guaranteed.
template<class Iterator, hpuint t_degree>
class SurfaceSplineConstrainerBEZ {
     static_assert(t_degree > 0, "A constant spline is not implemented.");
     using Indexer = ControlPointIndexerBEZ<t_degree>;
     using real = double;//__float128;//double//TODO: make class template parameter
     using Transition = hpvec3;
     using Triplet = Eigen::Triplet<real>;

public:
     using Matrix = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>;
     using PermutationMatrix = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>;
     using SparseMatrix = Eigen::SparseMatrix<real>;
     using Vector = Eigen::Matrix<real, Eigen::Dynamic, 1>;

     SurfaceSplineConstrainerBEZ(Iterator begin, Iterator end, hpreal epsilon = EPSILON)
          : m_begin(begin), m_end(end), m_nTriangles(std::distance(begin, end)) {
          //TODO: factor out control point indexing into special class
          m_indices.resize(m_nTriangles * NUMBER_OF_CONTROL_POINTS, -1);
          m_nControlPoints = 0;
          boost::dynamic_bitset<> done(m_nTriangles);
          hpuint index = 0;
          auto c = m_indices.begin();
          for(Iterator i = m_begin; i != m_end; ++i) {
               auto v = *i;
               auto n = std::get<0>(v);
               hpuint n0 = std::get<0>(n);
               hpuint n1 = std::get<1>(n);
               hpuint n2 = std::get<2>(n);
               if(n0 != UNULL && done[n0]) {
                    if(n1 != UNULL && done[n1]) {
                         if(n2 != UNULL && done[n2]) doIndexing<false, false, false>(done, c, index, n0, n1, n2, epsilon);
                         else doIndexing<false, false, true>(done, c, index, n0, n1, n2, epsilon);
                    } else {
                         if(n2 != UNULL && done[n2]) doIndexing<false, true, false>(done, c, index, n0, n1, n2, epsilon);
                         else doIndexing<false, true, true>(done, c, index, n0, n1, n2, epsilon);
                    }
               } else {
                    if(n1 != UNULL && done[n1]) {
                         if(n2 != UNULL && done[n2]) doIndexing<true, false, false>(done, c, index, n0, n1, n2, epsilon);
                         else doIndexing<true, false, true>(done, c, index, n0, n1, n2, epsilon);
                    } else {
                         if(n2 != UNULL && done[n2]) doIndexing<true, true, false>(done, c, index, n0, n1, n2, epsilon);
                         else doIndexing<true, true, true>(done, c, index, n0, n1, n2, epsilon);
                    }
               }
               done[index] = true;
               ++index;
          }

          m_zeroed.resize(m_nControlPoints, false);
     }

     void constrain(hpuint continuity) {
          assert(continuity > 0);

          std::vector<hpuint> indices = getLinearSystemIndices();
          hpuint nRows;
          std::vector<Triplet> triplets;
          if(continuity == 1) {
               if(m_zeroed.count() > 0) {
                    TripletBuilderC1 builder(*this, indices);
                    auto result = builder.build();
                    triplets = std::move(result.first);
                    nRows = result.second;
               } else {
                    TripletBuilderC1 builder(*this);
                    auto result = builder.build();
                    triplets = std::move(result.first);
                    nRows = result.second;
               }
          } else if(continuity == 2) {
               std::cerr << "TODO\n";//TODO
          } else {
               std::cerr << "TODO\n";//TODO
          }

          m_matrix = SparseMatrix(nRows, m_nControlPoints - m_zeroed.count());
          m_matrix.setFromTriplets(triplets.begin(), triplets.end());

          std::cout << "INFO: There are " << m_nControlPoints << " control points.\n";
          std::cout << "INFO: There are " << m_zeroed.count() << " control points that are set to zero.\n";
          std::cout << "INFO: There are " << nRows << " continuity conditions.\n";
          std::cout << "INFO: Constraint matrix is " << m_matrix.rows() << 'x' << m_matrix.cols() << ".\n";

          m_matrix = m_matrix.transpose();
          m_matrix.makeCompressed();
     }

     /**
      * Samples the surface and checks if there exist solutions such that all the samples are positive.
      *
      * @param[in] nSamples number of samples on one edge of parameter triangle, the triangle is then sampled uniformly respecting this value
      */
     //NOTE: The epsilon needs to be relatively large because the solver returns the trivial solution when epsilon is small.  This is probably due to numerical errors.
     bool existsPositiveSolution(hpuint nSamples = 5, double epsilon = 0.001) {
          auto nVariables = m_dimension + SurfaceUtilsBEZ::getNumberOfControlPoints(nSamples - 1) * m_nTriangles;
          auto lp = make_lp(0, nVariables);
          if(lp == nullptr) throw std::runtime_error("Failed to initialize the linear program.");
          double factors[nVariables];
          double objective[nVariables];
          memset(objective, 0, sizeof(double) * nVariables);
          std::fill(objective + m_dimension, objective + nVariables, 1.0);
          if(!set_obj_fn(lp, objective - 1)) {
               delete_lp(lp);
               throw std::runtime_error("Failed to set the objective function.");
          }
          set_maxim(lp);//maximize objective
          for(auto i = 1u; i <= m_dimension; ++i) if(!set_unbounded(lp, i)) {
               delete_lp(lp);
               throw std::runtime_error("Failed to set " + std::to_string(i) + "th variable as unbounded.");
          }
          for(auto i = m_dimension + 1; i <= nVariables; ++i) if(!set_binary(lp, i, TRUE)) {
               delete_lp(lp);
               throw std::runtime_error("Failed to set " + std::to_string(i) + "th variable as binary.");
          }
          set_add_rowmode(lp, TRUE);
          if(!add_constraint(lp, objective - 1, GE, 1)) {
               delete_lp(lp);
               throw std::runtime_error("Failed to set 'at least one' constraint.");
          }
          try {
               if(m_zeroed.count() > 0) addPositiveSamplesConstraints<true>(lp, nSamples, epsilon);
               else addPositiveSamplesConstraints<false>(lp, nSamples, epsilon);
          } catch(...) {
               delete_lp(lp);
               throw;
          }
          set_add_rowmode(lp, FALSE);
          set_verbose(lp, CRITICAL);
          //write_LP(lp, stdout);
          int code;
          if((code = solve(lp)) == OPTIMAL) {
               get_variables(lp, factors);
               bool zero = true;
               for(auto f = factors, end = f + m_dimension; f != end; ++f) {
                    std::cout << *f << '\n';
                    zero &= *f < epsilon;
               }
               if(zero) std::cout << "WARN: The optimal solution is zero, which is probably due to numerical errors.  T)ry a larger epsilon.\n";
          }
          delete_lp(lp);
          return code == OPTIMAL;
     }

     const Matrix& getBasis() const { return m_q2; }

     boost::optional<std::vector<real> > getBoundedSolution(double minimum, double maximum, double epsilon = EPSILON) const {
          if(maximum < epsilon || minimum > -epsilon) return doGetBoundedSolutionSameSign(minimum, maximum, epsilon);
          else return doGetBoundedSolutionDifferentSigns(minimum, maximum, epsilon);
     }

     /**
      * &param[values] set of (control point index, value)
      */
     boost::optional<std::pair<std::vector<Point1D>, std::vector<hpuint> > > getControlPoints(const std::vector<std::pair<hpuint, double> >& values, double epsilon = EPSILON) const {
          if(auto s = isSolution(values, epsilon)) return getControlPoints(*s);
          else return boost::none;
     }

     //NOTE: Coefficients should be const but Eigen::Map does not accept pointer to const.
     std::pair<std::vector<Point1D>, std::vector<hpuint> > getControlPoints(std::vector<real>& coefficients, double epsilon = EPSILON) const {
          assert(coefficients.size() == m_dimension);//TODO: exception
          Eigen::Map<Vector> temp(coefficients.data(), m_dimension);
          return getControlPoints(m_q2 * temp);
     }

     std::pair<std::vector<Point1D>, std::vector<hpuint> > getControlPoints(const Vector& solution) const {
          assert(solution.rows() == (m_nControlPoints - m_zeroed.count()));//TODO: exception
          std::vector<Point1D> points;
          points.reserve(m_nControlPoints);
          hpuint i = 0;
          hpuint j = 0;//index nonzero control point
          while(i < m_nControlPoints) {
               if(m_zeroed[i]) points.emplace_back(0.0);
               else {
                    points.emplace_back(solution(j));
                    ++j;
               }
               ++i;
          }
          return std::make_pair(std::move(points), m_indices);
     }

     hpuint getDimension() const { return m_dimension; }

     const std::vector<hpuint>& getIndices() const { return m_indices; }

     Vector getMinimalCurvatureSolution(const Vector3D& v, hpuint i, double epsilon = EPSILON) const {
          std::vector<Triplet> triplets;
          hpuint nRows = 0;
          auto a0 = v.x * v.x;
          auto a1 = 2.0 * v.x * v.y;
          auto a2 = v.y * v.y;
          auto a3 = 2.0 * v.x * v.z;
          auto a4 = 2.0 * v.y * v.z;
          auto a5 = v.z * v.z;
          std::vector<hpuint> indices = getLinearSystemIndices();
          auto c = m_indices.begin(), end = m_indices.end();
          while(c != end) {
               auto r0 = c;
               auto r1 = r0 + t_degree + 1;
               auto r2 = r1 + t_degree;
               auto stride = t_degree;
               c += NUMBER_OF_CONTROL_POINTS;
               while(r2 != c) {
                    auto re = r2 + (--stride);
                    while(r2 != re) {
                         bool atLeastOne = false;
                         if(!m_zeroed[*r0]) {
                              triplets.emplace_back(nRows, indices[*r0], a0);
                              atLeastOne = true;
                         }
                         ++r0;
                         if(!m_zeroed[*r0]) {
                              triplets.emplace_back(nRows, indices[*r0], a1);
                              atLeastOne = true;
                         }
                         if(!m_zeroed[*(r0+1)]) {
                              triplets.emplace_back(nRows, indices[*(r0+1)], a2);
                              atLeastOne = true;
                         }
                         if(!m_zeroed[*r1]) {
                              triplets.emplace_back(nRows, indices[*r1], a3);
                              atLeastOne = true;
                         }
                         ++r1;
                         if(!m_zeroed[*r1]) {
                              triplets.emplace_back(nRows, indices[*r1], a4);
                              atLeastOne = true;
                         }
                         if(!m_zeroed[*r2]) {
                              triplets.emplace_back(nRows, indices[*r2], a5);
                              atLeastOne = true;
                         }
                         if(atLeastOne) ++nRows;
                         ++r2;
                    }
                    ++r1;
                    r0 += 2;
               }
          }
          auto C = SparseMatrix(nRows, m_nControlPoints - m_zeroed.count());
          C.setFromTriplets(triplets.begin(), triplets.end());
          C.makeCompressed();
          auto q2 = m_q2;
          q2.row(indices[i]) = Vector::Zero(m_q2.cols());
          Matrix A = C * q2;
          Vector b = C.col(indices[i]) * -1.0;
          Vector x = A.colPivHouseholderQr().solve(b);
          Vector solution = m_q2 * x;
          solution[indices[i]] = 1.0;
          return solution;
     }

     hpuint getNumberOfTriangles() const { return m_nTriangles; }

     std::vector<hpuint> getSupport() const {
          std::vector<hpuint> support;
          support.reserve(m_indices.size() - m_zeroed.count());
          for(hpuint i = 0, end = m_indices.size(); i < end; ++i) if(!m_zeroed[m_indices[i]]) support.push_back(i);
          return support;
     }

     bool isSolution(const Vector& solution, double epsilon = EPSILON) const {
          Vector temp = m_matrix.transpose() * solution;
          return temp.isZero(epsilon);
     }

     /**
      * @param[values] set of (control point index, value)
      */
     boost::optional<Vector> isSolution(const std::vector<std::pair<hpuint, real> >& values, double epsilon = EPSILON) const {
          assert(values.size() == m_dimension);//TODO: expand to allow any number of values and handle all cases
          std::vector<hpuint> indices = getLinearSystemIndices();
          auto v = values.cbegin();
          auto end = values.cend();
          hpuint i = 0;
          Matrix A(m_dimension, m_dimension);
          Vector b(m_dimension);
          if(m_zeroed.count() == 0) std::cerr << "ERROR: The case when there are no zeroed control points is not implemented.\n";//TODO
          while(i < m_dimension && v != end) {//TODO: if no zeroed control points, do not use indices
               auto value = *v;
               if(m_zeroed[value.first]) return boost::none;//TODO: throw error
               A.row(i) = m_q2.row(indices[value.first]);
               b(i) = value.second;
               ++i;
               ++v;
          }
          auto lu = A.fullPivLu();
          if(lu.rank() == m_dimension) {
               auto x = lu.solve(b);
               return (Vector)(m_q2 * x);
          } else std::cout << "ERROR: Control points are not independent.\n";//TODO: check if values in image of matrix
          return boost::none;
     }

     bool solve(double epsilon = EPSILON) {
          using Ordering = Eigen::COLAMDOrdering<int>;

          Eigen::SparseQR<SparseMatrix, Ordering> solver;
          solver.setPivotThreshold(epsilon);
          solver.compute(m_matrix);

          if(solver.info() != Eigen::Success) {
               std::cerr << "ERROR: Failed to compute QR decomposition.\n";//TODO: throw error
               return false;
          }

          auto rank = solver.rank();

          m_dimension = m_nControlPoints - m_zeroed.count() - rank;
          std::cout << "INFO: Solution space has dimension " << m_dimension << ".\n";

          auto r = solver.matrixR();
          m_r1 = r.topLeftCorner(rank, rank);
          m_r2 = r.topRightCorner(rank, m_dimension);

          /*std::cout << "INFO: R is " << r.rows() << 'x' << r.cols() << ".\n";
          std::cout << "INFO: R1 is " << m_r1.rows() << 'x' << m_r1.cols() << ".\n";
          std::cout << "INFO: R2 is " << m_r2.rows() << 'x' << m_r2.cols() << ".\n";*/
          
          SparseMatrix q;
          q = solver.matrixQ();
          m_q2 = q.rightCols(m_dimension);

          //NOTE: Sanity check.  Make sure A*Q2 = 0.
          Matrix temp = m_matrix.transpose() * m_q2;
          assert(temp.isZero(epsilon));

          /*std::cout << "INFO: Q is " << q.rows() << 'x' << q.cols() << ".\n";
          std::cout << "INFO: Q2 is " << m_q2.rows() << 'x' << m_q2.cols() << ".\n";*/

          m_p = solver.colsPermutation();

          return true;
     }

     void unzero() { m_zeroed.reset(); }

     void zero(hpuint triangle) { for(auto c = m_indices.cbegin() + triangle * NUMBER_OF_CONTROL_POINTS, end = c + NUMBER_OF_CONTROL_POINTS; c != end; ++c) m_zeroed[*c] = true; }//zero control points of triangle

     void zero(const std::vector<hpuint>& triangles) { for(hpuint triangle : triangles) zero(triangle); }

     void zero(hpuint triangle, hpuint neighbor, hpuint rows) {//zero rows rows of control points in triangle on edge shared by triangle and neighbor
          assert(rows < (t_degree + 1));
          auto v = *(m_begin + triangle);
          auto c = m_indices.cbegin() + triangle * NUMBER_OF_CONTROL_POINTS;
          if(v.n0 == neighbor) zero<0>(c, rows);
          else if(v.n1 == neighbor) zero<1>(c, rows);
          else {
               assert(v.n2 == neighbor);
               zero<2>(c, rows);
          }
     }

     void zero(const std::vector<std::pair<hpuint, hpuint> >& pairs, hpuint rows) {
          for(auto& pair : pairs) zero(pair.first, pair.second, rows);
     }

     static constexpr hpuint NUMBER_OF_CONTROL_POINTS = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;//TODO: private
private:

     const Iterator m_begin;
     hpuint m_dimension;
     const Iterator m_end;
     std::vector<hpuint> m_indices;//indices of the control points in the projective structure mesh
     hpuint m_nControlPoints;
     const hpuint m_nTriangles;
     SparseMatrix m_matrix;
     PermutationMatrix m_p;
     Matrix m_q2;
     Matrix m_r1;
     Matrix m_r2;
     boost::dynamic_bitset<> m_zeroed;

     static int solve(lprec* lp) {
          int code;
          switch(code = ::solve(lp)) {
          case OPTIMAL:
               std::cout << "INFO: Found an optimal solution.\n";
               break;
          case INFEASIBLE:
               std::cout << "INFO: The model is infeasible.\n";
               break;
          case UNBOUNDED:
               std::cout << "INFO: The model is unbounded.\n";
               break;
          case DEGENERATE:
               std::cout << "INFO: The model is degenerate.\n";
               break;
          default:
               std::cerr << "ERROR: Solver returned " << code << ".\n";
               break;
          }
          return code;
     }

     template<class TripletBuilderImpl>
     class TripletBuilder {
     public:
          TripletBuilder(const SurfaceSplineConstrainerBEZ& constrainer)
               : m_constrainer(constrainer), m_indices(nullptr) {}

          TripletBuilder(const SurfaceSplineConstrainerBEZ& constrainer, const std::vector<hpuint>& indices)
               : m_constrainer(constrainer), m_indices(&indices) {}

          std::pair<std::vector<Triplet>, hpuint> build() {
               m_row = 0;
               m_triplets.clear();
               boost::dynamic_bitset<> done(m_constrainer.m_nTriangles);
               auto index = 0u;
               if(m_indices == nullptr) {
                    for(auto i = m_constrainer.m_begin; i != m_constrainer.m_end; ++i) {
                         auto v = *i;
                         auto n = std::get<0>(v);
                         auto n0 = std::get<0>(n);
                         auto n1 = std::get<1>(n);
                         auto n2 = std::get<2>(n);
                         auto t = std::get<1>(v);
                         auto& t0 = std::get<0>(t);
                         auto& t1 = std::get<1>(t);
                         auto& t2 = std::get<2>(t);
                         if(n0 != UNULL && !done[n0]) static_cast<TripletBuilderImpl*>(this)->template build<0, false>(t0, index, n0);
                         if(n1 != UNULL && !done[n1]) static_cast<TripletBuilderImpl*>(this)->template build<1, false>(t1, index, n1);
                         if(n2 != UNULL && !done[n2]) static_cast<TripletBuilderImpl*>(this)->template build<2, false>(t2, index, n2);
                         done[index] = true;
                         ++index;
                    }
               } else {
                    for(auto i = m_constrainer.m_begin; i != m_constrainer.m_end; ++i) {
                         auto v = *i;
                         auto n = std::get<0>(v);
                         auto n0 = std::get<0>(n);
                         auto n1 = std::get<1>(n);
                         auto n2 = std::get<2>(n);
                         auto t = std::get<1>(v);
                         auto& t0 = std::get<0>(t);
                         auto& t1 = std::get<1>(t);
                         auto& t2 = std::get<2>(t);
                         if(n0 != UNULL && !done[n0]) static_cast<TripletBuilderImpl*>(this)->template build<0, true>(t0, index, n0);
                         if(n1 != UNULL && !done[n1]) static_cast<TripletBuilderImpl*>(this)->template build<1, true>(t1, index, n1);
                         if(n2 != UNULL && !done[n2]) static_cast<TripletBuilderImpl*>(this)->template build<2, true>(t2, index, n2);
                         done[index] = true;
                         ++index;
                    }
               }
               return std::make_pair(std::move(m_triplets), m_row);
          }

     protected:
          const SurfaceSplineConstrainerBEZ& m_constrainer;
          const std::vector<hpuint>* m_indices;
          hpuint m_row;
          std::vector<Triplet> m_triplets;

     };//TripletBuilder

     class TripletBuilderC1 : public TripletBuilder<TripletBuilderC1> {
          friend class TripletBuilder<TripletBuilderC1>;

     public:
          TripletBuilderC1(const SurfaceSplineConstrainerBEZ& constrainer)
               : TripletBuilder<TripletBuilderC1>(constrainer) {}

          TripletBuilderC1(const SurfaceSplineConstrainerBEZ& constrainer, const std::vector<hpuint>& indices)
               : TripletBuilder<TripletBuilderC1>(constrainer, indices) {}

          using TripletBuilder<TripletBuilderC1>::build;

     protected:
          template<hpuint t_e0, bool t_zeroed>
          void build(const Transition& transition, hpuint i0, hpuint i1) {
               auto c0 = this->m_constrainer.m_indices.cbegin() + NUMBER_OF_CONTROL_POINTS * i0;
               auto c1 = this->m_constrainer.m_indices.cbegin() + NUMBER_OF_CONTROL_POINTS * i1;
               auto r00 = Indexer::template getIndices<Indexer::CLOCKWISE, t_e0, 0>();
               auto r01 = Indexer::template getIndices<Indexer::CLOCKWISE, t_e0, 1>();
               auto r11 = this->m_constrainer.template getIndices<Indexer::COUNTER_CLOCKWISE, 1>(i1, i0);

               auto end = r00 + (t_degree + 1);
               auto cx = *r00;
               ++r00;
               do {
                    auto cy = *r00;
                    auto cz = *r01;
                    auto cr = *r11;
                    if(t_zeroed) {
                         auto dx = *(c0 + cx);
                         auto dy = *(c0 + cy);
                         auto dz = *(c0 + cz);
                         auto dr = *(c1 + cr);
                         bool atLeastOne = false;
                         if(!this->m_constrainer.m_zeroed[dx]) {
                              this->m_triplets.emplace_back(this->m_row, (*this->m_indices)[dx], transition.x);
                              atLeastOne = true;
                         }
                         if(!this->m_constrainer.m_zeroed[dy]) {
                              this->m_triplets.emplace_back(this->m_row, (*this->m_indices)[dy], transition.y);
                              atLeastOne = true;
                         }
                         if(!this->m_constrainer.m_zeroed[dz]) {
                              this->m_triplets.emplace_back(this->m_row, (*this->m_indices)[dz], transition.z);
                              atLeastOne = true;
                         }
                         if(!this->m_constrainer.m_zeroed[dr]) {
                              this->m_triplets.emplace_back(this->m_row, (*this->m_indices)[dr], -1);
                              atLeastOne = true;
                         }
                         if(atLeastOne) ++this->m_row;
                    } else {
                         this->m_triplets.emplace_back(this->m_row, *(c0 + cx), transition.x);
                         this->m_triplets.emplace_back(this->m_row, *(c0 + cy), transition.y);
                         this->m_triplets.emplace_back(this->m_row, *(c0 + cz), transition.z);
                         this->m_triplets.emplace_back(this->m_row, *(c1 + cr), -1);
                         ++this->m_row;
                    }
                    cx = cy;
                    ++r00;
                    ++r01;
                    ++r11;
               } while(r00 != end);
          }
          
     };//TripletBuilderC1

     /*class TripletBuilderC2 : public TripletBuilder<TripletBuilderC2> {
          friend class TripletBuilder<TripletBuilderC2>;

     public:
          TripletBuilderC2(const SurfaceSplineConstrainerBEZ& constrainer)
               : TripletBuilder<TripletBuilderC2>(constrainer) {}

          TripletBuilderC2(const SurfaceSplineConstrainerBEZ& constrainer, const std::vector<hpuint>& indices)
               : TripletBuilder<TripletBuilderC2>(constrainer, indices) {}

          using TripletBuilder<TripletBuilderC2>::build;

     protected:
          template<hpuint t_e0, bool t_zeroed>
          void build(const Transition& transition, hpuint i0, hpuint i1) {
               auto c0 = this->m_constrainer.m_indices.cbegin() + NUMBER_OF_CONTROL_POINTS * i0;
               auto c1 = this->m_constrainer.m_indices.cbegin() + NUMBER_OF_CONTROL_POINTS * i1;

          }

     };//TripletBuilderC2*/

     //class TripletBuilderCk

     template<bool zeroed>
     void addPositiveSamplesConstraints(lprec* lp, hpuint nSamples, double epsilon) {
          const double nolispe = -1.0 / epsilon;
          auto nVariables = m_dimension + SurfaceUtilsBEZ::getNumberOfControlPoints(nSamples - 1) * m_nTriangles;
          auto matrix = SurfaceUtilsBEZ::getEvaluationMatrix<t_degree>(nSamples);
          std::vector<hpuint> indices;
          if(zeroed) indices = getLinearSystemIndices();
          //TODO: eliminate common points on edges
          auto sample = m_dimension;
          double factors[nVariables];
          for(auto tb = m_indices.cbegin(), te = tb + NUMBER_OF_CONTROL_POINTS, ts = m_indices.cend(); tb != ts; tb = te, te += NUMBER_OF_CONTROL_POINTS) {
               memset(factors, 0, sizeof(double) * nVariables);
               auto t = tb;
               bool nonzero = false;
               for(auto m = matrix.cbegin(), end = matrix.cend(); m != end; ++m) {
                    if(!zeroed || !m_zeroed[*t]) {
                         Vector r = m_q2.row((zeroed) ? indices[*t] : *t);
                         if(!r.isZero(epsilon) && std::abs(*m) > epsilon) {
                              auto temp = (*m) * r;
                              auto f = factors;
                              for(auto i = 0u; i < m_dimension; ++i, ++f) *f += temp[i];
                              nonzero = true;
                         }
                    }
                    ++t;
                    if(t == te) {
                         //NOTE: Need to shift back by one because the solver is one-based.
                         if(nonzero) {
                              if(!add_constraint(lp, factors - 1, GE, 0.0)) throw std::runtime_error("Failed to add positive sample constraint.");
                              factors[sample] = nolispe;
                              if(!add_constraint(lp, factors - 1, GE, nolispe + epsilon)) throw std::runtime_error("Failed to add strictly positive sample constraint.");
                         } else {
                              factors[sample] = 1.0;
                              if(!add_constraint(lp, factors - 1, EQ, 0.0)) throw std::runtime_error("Failed to add zero sample constraint.");
                         }
                         ++sample;
                         memset(factors, 0, sizeof(double) * nVariables);
                         t = tb;
                         nonzero = false;
                    }
               }
          }
     }

     boost::optional<std::vector<real> > doGetBoundedSolutionDifferentSigns(double minimum, double maximum, double epsilon = EPSILON) const {
          const hpuint nControlPoints = m_nControlPoints - m_zeroed.count();
          auto nVariables = m_dimension + 3 * nControlPoints;
          auto lp = make_lp(0, nVariables);
          if(lp == NULL) return boost::none;
          std::vector<real> values;
          std::vector<double> factors(nVariables, 0.0);
          std::vector<std::vector<double> > cache;
          double objective[nVariables];
          memset(objective, 0, sizeof(double) * nVariables);
          std::fill(objective + m_dimension, objective + nVariables - nControlPoints, 1.0);
          if(!set_obj_fn(lp, objective - 1)) goto cleanup;
          set_maxim(lp);//maximize objective
          for(hpuint i = 1; i <= m_dimension; ++i) {
               if(!set_unbounded(lp, i)) {
                    std::cerr << "ERROR: Failed to set " << i << "th variable as unbounded.\n";
                    goto cleanup;
               }
          }
          for(hpuint i = (nControlPoints << 1) + 1; i <= nVariables; ++i) {
               if(!set_binary(lp, i, TRUE)) {
                    std::cerr << "ERROR: Failed to set " << i << "th variable as binary.\n";
                    goto cleanup;
               }
          }
          set_add_rowmode(lp, TRUE);
          for(hpuint i = 0, end = nControlPoints; i < end; ++i) {
               bool zero = true;
               for(hpuint j = 0; j < m_dimension; ++j) {
                    factors[j] = m_q2(i, j);
                    zero &= factors[j] < epsilon;
               }
               hpuint i0 = m_dimension + (i << 1);
               hpuint i1 = i0 + 1;
               hpuint i2 = m_dimension + (nControlPoints << 1) + i;
               if(zero || std::find_if(cache.cbegin(), cache.cend(), [&](const std::vector<double>& other) -> bool {
                    for(hpuint i = 0; i < m_dimension; ++i) if(std::abs(other[i] - factors[i]) > epsilon) return false;
                    return true;
               }) != cache.cend()) {
                    std::fill(factors.begin(), factors.begin() + m_dimension, 0.0);//set to zero just to be sure
                    factors[i0] = 1.0;
                    if(!add_constraint(lp, factors.data() - 1, EQ, 0.0)) {
                         std::cerr << "ERROR: Failed to set +" << i << "th control point to zero.\n";
                         goto cleanup;
                    }
                    factors[i0] = 0.0;
                    factors[i1] = 1.0;
                    if(!add_constraint(lp, factors.data() - 1, EQ, 0.0)) {
                         std::cerr << "ERROR: Failed to set -" << i << "th control point to zero.\n";
                         goto cleanup;
                    }
                    factors[i1] = 0.0;
               } else {
                    cache.push_back(factors); 
                    factors[i0] = -1.0;
                    factors[i1] = 1.0;
                    if(!add_constraint(lp, factors.data() - 1, EQ, 0.0)) {
                         std::cerr << "ERROR: Failed to set " << i << "th control point constraint.\n";
                         goto cleanup;
                    }
                    std::fill(factors.begin(), factors.begin() + m_dimension, 0.0);
                    factors[i0] = 1.0;
                    factors[i1] = 0.0;
                    factors[i2] = -maximum;
                    if(!add_constraint(lp, factors.data() - 1, LE, 0.0)) {
                         std::cerr << "ERROR: Failed to set " << i << "th maximum constraint.\n";
                         goto cleanup;
                    }
                    factors[i0] = 0.0;
                    factors[i1] = 1.0;
                    factors[i2] = std::abs(minimum);
                    if(!add_constraint(lp, factors.data() - 1, LE, std::abs(minimum))) {
                         std::cerr << "ERROR: Failed to set " << i << "th minimum constraint.\n";
                         goto cleanup;
                    }
                    factors[i1] = 0.0;
                    factors[i2] = 0.0;
               }
          }
          set_add_rowmode(lp, FALSE);
          set_verbose(lp, CRITICAL);
          //write_LP(lp, stdout);
          if(solve(lp) == OPTIMAL) {
               double temp[nVariables];
               values.reserve(m_dimension);
               get_variables(lp, temp);
               for(auto f = temp, end = f + m_dimension; f != end; ++f) {
                    std::cout << *f << '\n';
                    values.push_back(*f);
               }
          }
          cleanup: delete_lp(lp);
          if(values.size() > 0) return values;
          else return boost::none;
     }

     boost::optional<std::vector<real> > doGetBoundedSolutionSameSign(double minimum, double maximum, double epsilon = EPSILON) const {
          const hpuint nControlPoints = m_nControlPoints - m_zeroed.count();
          hpuint nVariables = m_dimension + nControlPoints;
          lprec* lp = make_lp(0, nVariables);
          if(lp == NULL) return boost::none;
          std::vector<real> values;
          double factors[nVariables];
          memset(factors, 0, sizeof(double) * nVariables);
          double objective[nVariables];
          memset(objective, 0, sizeof(double) * nVariables);
          std::fill(objective + m_dimension, objective + nVariables, 1.0);
          if(!set_obj_fn(lp, objective - 1)) goto cleanup;
          if(minimum > -epsilon) set_maxim(lp);//maximize objective if minimum >= 0
          for(hpuint i = 1; i <= nVariables; ++i) {
               if(!set_unbounded(lp, i)) {
                    std::cerr << "ERROR: Failed to set " << i << "th variable as unbounded.\n";
                    goto cleanup;
               }
          }
          set_add_rowmode(lp, TRUE);
          for(hpuint i = 0, end = nControlPoints; i < end; ++i) {
               bool zero = true;
               for(hpuint j = 0; j < m_dimension; ++j) {
                    factors[j] = m_q2(i, j);
                    zero &= factors[j] < epsilon;
               }
               if(zero) {
                    memset(factors, 0, sizeof(double) * m_dimension);//set to zero just to be sure
                    factors[m_dimension + i] = 1.0;
                    if(!add_constraint(lp, factors - 1, EQ, 0.0)) {
                         std::cerr << "ERROR: Failed to set " << i << "th control point to zero.\n";
                         goto cleanup;
                    }
                    factors[m_dimension + i] = 0.0;
               } else {
                    factors[m_dimension + i] = -1.0;
                    if(!add_constraint(lp, factors - 1, EQ, 0.0)) {
                         std::cerr << "ERROR: Failed to set " << i << "th control point constraint.\n";
                         goto cleanup;
                    }
                    memset(factors, 0, sizeof(double) * m_dimension);
                    factors[m_dimension + i] = 1.0;
                    if(!add_constraint(lp, factors - 1, LE, maximum)) {
                         std::cerr << "ERROR: Failed to set " << i << "th maximum constraint.\n";
                         goto cleanup;
                    }
                    if(!add_constraint(lp, factors - 1, GE, minimum)) {
                         std::cerr << "ERROR: Failed to set " << i << "th minimum constraint.\n";
                         goto cleanup;
                    }
                    factors[m_dimension + i] = 0.0;
               }
          }
          set_add_rowmode(lp, FALSE);
          set_verbose(lp, CRITICAL);
          //write_LP(lp, stdout);
          if(solve(lp) == OPTIMAL) {
               values.reserve(m_dimension);
               get_variables(lp, factors);
               for(auto f = factors, end = f + m_dimension; f != end; ++f) {
                    std::cout << *f << '\n';
                    values.push_back(*f);
               }
          }
          cleanup: delete_lp(lp);
          if(values.size() > 0) return values;
          else return boost::none;
     }

     template<bool t_e0, bool t_e1, bool t_e2>
     void doIndexing(const boost::dynamic_bitset<>& done, std::vector<hpuint>::iterator& c, hpuint index, hpuint n0, hpuint n1, hpuint n2, hpreal epsilon) {
          typename Indexer::Iterator f0;
          typename Indexer::Iterator f1;
          typename Indexer::Iterator f2;
          std::vector<hpuint>::const_iterator i0;
          std::vector<hpuint>::const_iterator i1;
          std::vector<hpuint>::const_iterator i2;

          if(!t_e0) {
               f0 = getIndices<Indexer::CLOCKWISE, 0>(n0, index);
               i0 = m_indices.cbegin() + NUMBER_OF_CONTROL_POINTS * n0;
          }
          if(!t_e1) {
               f1 = getIndices<Indexer::CLOCKWISE, 0>(n1, index);
               i1 = m_indices.cbegin() + NUMBER_OF_CONTROL_POINTS * n1;
          }
          if(!t_e2) {
               f2 = getIndices<Indexer::COUNTER_CLOCKWISE, 0>(n2, index);
               i2 = m_indices.cbegin() + NUMBER_OF_CONTROL_POINTS * n2;
          }

          //first row
          if(t_e0 && t_e2) {
               hpuint d = findIndex(done, index, n0, n2, epsilon);
               if(d == -1u) *c = m_nControlPoints++;
               else *c = d;
          } else if(t_e2) *c = *(i0 + *f0);
          else *c = *(i2 + *f2);
          if(t_e0) {
               hpuint limit = m_nControlPoints + (t_degree - 1);
               while(m_nControlPoints < limit) *(++c) = m_nControlPoints++;
          } else {
               auto limit = f0 + (t_degree - 1);
               while(f0 != limit) *(++c) = *(i0 + *(++f0));
          }
          if(t_e0 && t_e1) {
               hpuint d = findIndex(done, index, n1, n0, epsilon);
               if(d == -1u) *(++c) = m_nControlPoints++;
               else *(++c) = d;
          } else if(t_e1) *(++c) = *(i0 + *(++f0));
          else *(++c) = *(i1 + *f1);
          //middle rows
          hpuint rowLength = t_degree - 1;
          while(rowLength > 0) {
               if(t_e2) *(++c) = m_nControlPoints++;
               else *(++c) = *(i2 + *(++f2));
               --rowLength;
               hpuint limit = m_nControlPoints + rowLength;
               while(m_nControlPoints < limit) *(++c) = m_nControlPoints++;
               if(t_e1) *(++c) = m_nControlPoints++;
               else *(++c) = *(i1 + *(++f1));
          }
          //last row
          if(t_e2 && t_e1) {
               hpuint d = findIndex(done, index, n2, n1, epsilon);
               if(d == -1u) *(++c) = m_nControlPoints++;
               else *(++c) = d;
          } else if(t_e1) *(++c) = *(i2 + *(++f2));
          else *(++c) = *(i1 + *(++f1));
          ++c;
     }

     hpuint findIndex(const boost::dynamic_bitset<>& done, hpuint index, hpuint neighborCW, hpuint neighborCCW, hpreal epsilon) {
          hpuint p = index;
          hpuint n = neighborCW;
          if(neighborCW != neighborCCW) {
               while(n != UNULL && (n != index || (n == index && p != neighborCCW)) && !done[n]) {
                    auto nv = *(m_begin + n);
                    auto nn = std::get<0>(nv);
                    hpuint n0 = std::get<0>(nn);
                    hpuint n1 = std::get<1>(nn);
                    hpuint n2 = std::get<2>(nn);
                    if(n0 == p) {
                         p = n;
                         n = n1;
                    } else if(n1 == p) {
                         p = n;
                         n = n2;
                    } else {
                         assert(n2 == p);
                         p = n;
                         n = n0;
                    }
               }
          }
          if(n == UNULL) {//go counterclockwise if previous loop reached the border
               p = index;
               n = neighborCCW;
               while(n != UNULL && (n != index || (n == index && p != neighborCW)) && !done[n]) {
                    auto nv = *(m_begin + n);
                    auto nn = std::get<0>(nv);
                    hpuint n0 = std::get<0>(nn);
                    hpuint n1 = std::get<1>(nn);
                    hpuint n2 = std::get<2>(nn);
                    if(n0 == p) {
                         p = n;
                         n = n2;
                    } else if(n1 == p) {
                         p = n;
                         n = n0;
                    } else {
                         assert(n2 == p);
                         p = n;
                         n = n1;
                    }
               }
          }
          if(n != UNULL && done[n]) {
               auto i = m_indices.cbegin() + NUMBER_OF_CONTROL_POINTS * n;
               auto v = *(m_begin + n);
               auto nn = std::get<0>(v);
               if(std::get<0>(nn) == p) return *(i + t_degree);
               else if(std::get<1>(nn) == p) return *(i + (NUMBER_OF_CONTROL_POINTS - 1));
               else {
                    assert(std::get<2>(nn) == p);
                    return *i;
               }
          }
          return -1;
     }

     template<hpuint t_direction, hpuint t_row>
     typename Indexer::Iterator getIndices(hpuint index, hpuint neighbor) const {
          auto v = *(m_begin + index);
          auto n = std::get<0>(v);
          if(std::get<0>(n) == neighbor) return Indexer::template getIndices<t_direction, 0, t_row>();
          else if(std::get<1>(n) == neighbor) return Indexer::template getIndices<t_direction, 1, t_row>();
          else {
               assert(std::get<2>(n) == neighbor);
               return Indexer::template getIndices<t_direction, 2, t_row>();
          }
     }

     std::vector<hpuint> getLinearSystemIndices() const {
          std::vector<hpuint> indices;//real indices with zeroed control points removed
          indices.reserve(m_nControlPoints);
          if(m_zeroed.count() > 0) {
               hpuint index = -1;
               for(hpuint i = 0; i < m_nControlPoints; ++i) {
                    if(m_zeroed[i]) indices.push_back(-1);
                    else indices.push_back(++index);
               }
          } else for(hpuint i = 0; i < m_nControlPoints; ++i) indices.push_back(i);
          return indices;
     }

     template<hpuint t_edge>
     void zero(std::vector<hpuint>::const_iterator c, hpuint rows) {
          for(hpuint i = 0; i < rows; ++i) {
               auto r = Indexer::template getIndices<Indexer::CLOCKWISE, t_edge>(i);
               auto end = r + (t_degree + 1 - i);
               while(r != end) {
                    m_zeroed[*(c + *r)] = true;
                    ++r;
               }
          }
     }

};//SurfaceSplineConstrainerBEZ
template<class Iterator>
using CubicSurfaceSplineConstrainerBEZ = SurfaceSplineConstrainerBEZ<Iterator, 3>;
template<class Iterator>
using LinearSurfaceSplineConstrainerBEZ = SurfaceSplineConstrainerBEZ<Iterator, 1>;
template<class Iterator>
using QuadraticSurfaceSplineConstrainerBEZ = SurfaceSplineConstrainerBEZ<Iterator, 2>;
template<class Iterator>
using QuarticSurfaceSplineConstrainerBEZ = SurfaceSplineConstrainerBEZ<Iterator, 4>;

template<hpuint degree, class T>
static auto make_constrainer(const T& t) -> SurfaceSplineConstrainerBEZ<decltype(t.cbegin()), degree> { return SurfaceSplineConstrainerBEZ<decltype(t.cbegin()), degree>(t.cbegin(), t.cend()); }

}//namespace happah

