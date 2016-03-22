// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

#include "happah/geometries/SurfaceSplineHEZ.h"
#include "happah/math/ProjectiveStructure.h"
#include "happah/math/TriangleRefinementScheme.h"
#include "happah/utils/Arrays.h"
#include "happah/utils/SurfaceSplineConstrainerBEZ.h"

namespace happah {

class BasisBuilder {
     using IndicesArrays = Arrays<hpuint>;

public:
     template<hpuint degree>
     using Constrainer = SurfaceSplineConstrainerBEZ<ProjectiveStructure3D::default_const_iterator, degree>;
     enum class Type { VERSION1, VERSION2, VERSION3 };

private:
     using Mode = mode::ProjectiveStructure3D;
     using View = view::ProjectiveStructure3D;

     template<int t_dummy, Type t_type>
     class Iterator;

     template<int t_dummy>
     class Iterator<t_dummy, Type::VERSION1> {
          using Fans = ProjectiveStructure3D::const_iterator<View::VERTICES, Mode::TETRAHEDRA>;

     public:
          using value_type = std::vector<hpuint>;

          static Iterator cbegin(const BasisBuilder& basisBuilder) { return Iterator(basisBuilder.getProjectiveStructure().cbegin<View::VERTICES, Mode::TETRAHEDRA>(), basisBuilder.getTriangleRefinementScheme()); }

          static Iterator cend(const BasisBuilder& basisBuilder) { return Iterator(basisBuilder.getProjectiveStructure().cend<View::VERTICES, Mode::TETRAHEDRA>(), basisBuilder.getTriangleRefinementScheme()); }

          Iterator(Fans i, const TriangleRefinementScheme& scheme)
               : m_i(i), m_nTriangles(scheme.getNumberOfTriangles()) {}

          Iterator& operator++() {
               ++m_i;
               return *this;
          }

          Iterator operator++(int) {
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          bool operator==(const Iterator& iterator) const { return iterator.m_i == m_i; }

          bool operator!=(const Iterator& iterator) const { return !(*this == iterator); }

          value_type operator*() const {
               std::vector<hpuint> support;
               for(auto i : *m_i) for(hpuint f = m_nTriangles * i, end = f + m_nTriangles; f < end; ++f) support.push_back(f);
               std::sort(support.begin(), support.end());
               return support;
          }

     private:
          Fans m_i;
          hpuint m_nTriangles;

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, Type::VERSION2> {
          using Fans = std::vector<std::vector<hpuint> >;
          using Rings = ProjectiveStructure3D::const_iterator<View::VERTICES, Mode::VERTICES>;

     public:
          using value_type = std::vector<hpuint>;

          static Iterator cbegin(const BasisBuilder& basisBuilder) {
               Fans fans;
               auto& refinedProjectiveStructure = basisBuilder.getRefinedProjectiveStructure();
               for(auto i = refinedProjectiveStructure.cbegin<View::VERTICES, Mode::TETRAHEDRA>(), end = refinedProjectiveStructure.cend<View::VERTICES, Mode::TETRAHEDRA>(); i != end; ++i) { fans.push_back(*i); }
               return Iterator(refinedProjectiveStructure.cbegin<View::VERTICES, Mode::VERTICES>(), std::move(fans)); 
          }

          static Iterator cend(const BasisBuilder& basisBuilder) { return Iterator(basisBuilder.getRefinedProjectiveStructure().cend<View::VERTICES, Mode::VERTICES>()); }

          Iterator(Rings ring, Fans fans)
               : m_fans(std::move(fans)), m_fan(m_fans.cbegin()), m_ring(ring) {}

          Iterator(Rings ring)
               : m_ring(ring) {}

          Iterator& operator++() {
               ++m_fan;
               ++m_ring;
               return *this;
          }

          Iterator operator++(int) {
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          bool operator==(const Iterator& iterator) const { return iterator.m_ring == m_ring; }

          bool operator!=(const Iterator& iterator) const { return !(*this == iterator); }

          value_type operator*() const {
               std::vector<hpuint> support;
               for(auto& f : *m_fan) support.push_back(f);
               for(auto r : *m_ring) {
                    auto& fan = m_fans[r];
                    for(auto f : fan) support.push_back(f);
               }
               std::sort(support.begin(), support.end());
               return support;
          }

     private:
          Fans m_fans;
          Fans::const_iterator m_fan;
          Rings m_ring;

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, Type::VERSION3> {
          using Fans = std::vector<std::vector<hpuint> >;
          using Vertices = ProjectiveStructure3D::const_iterator<View::TETRAHEDRA, Mode::VERTICES>;

     public:
          using value_type = std::vector<hpuint>;

          static Iterator cbegin(const BasisBuilder& basisBuilder) {
               Fans fans;
               auto& refinedProjectiveStructure = basisBuilder.getRefinedProjectiveStructure();
               for(auto i = refinedProjectiveStructure.cbegin<View::VERTICES, Mode::TETRAHEDRA>(), end = refinedProjectiveStructure.cend<View::VERTICES, Mode::TETRAHEDRA>(); i != end; ++i) { fans.push_back(*i); }
               return Iterator(refinedProjectiveStructure.cbegin<View::TETRAHEDRA, Mode::VERTICES>(), std::move(fans)); 
          }

          static Iterator cend(const BasisBuilder& basisBuilder) { return Iterator(basisBuilder.getRefinedProjectiveStructure().cend<View::TETRAHEDRA, Mode::VERTICES>()); }

          Iterator(Vertices vertices, Fans fans)
               : m_fans(std::move(fans)), m_fan(m_fans.cbegin()), m_vertices(vertices) {}

          Iterator(Vertices vertices)
               : m_vertices(vertices) {}

          Iterator& operator++() {
               ++m_fan;
               ++m_vertices;
               return *this;
          }

          Iterator operator++(int) {
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          bool operator==(const Iterator& iterator) const { return iterator.m_vertices == m_vertices; }

          bool operator!=(const Iterator& iterator) const { return !(*this == iterator); }

          value_type operator*() const {
               std::vector<hpuint> support;
               hpuint v0, v1, v2;
               std::tie(v0, v1, v2) = *m_vertices;
               auto& fan0 = m_fans[v0];
               for(auto f : fan0) support.push_back(f);
               auto& fan1 = m_fans[v1];
               for(auto f : fan1) support.push_back(f);
               auto& fan2 = m_fans[v2];
               for(auto f : fan2) support.push_back(f);
               std::sort(support.begin(), support.end());
               return support;
          }

     private:
          Fans m_fans;
          Fans::const_iterator m_fan;
          Vertices m_vertices;

     };//Iterator

public:
     BasisBuilder(const ProjectiveStructure3D& projectiveStructure, hpuint tetrahedron, Point3D p0, Point3D p1, Point3D p2, const TriangleRefinementScheme& scheme, hpreal epsilon = EPSILON);

     template<hpuint degree>
     std::vector<SurfaceSplineHEZ<Space1D, degree> > build(hpuint continuity, const IndicesArrays& supports, hpreal epsilon = EPSILON) {
          auto factory = [&](Constrainer<degree>& constrainer) -> std::pair<std::vector<Point1D>, std::vector<hpuint> > {
               hpreal minimum = -1;
               //if(constrainer.existsPositiveSolution()) minimum = 0;
               auto values = constrainer.getBoundedSolution(minimum, 1);
               if(!values) throw std::runtime_error("Failed to find bounded solution.");
               return constrainer.getControlPoints(*values);
          };
          return build<degree>(continuity, supports, factory, epsilon);
     }

     template<hpuint degree>
     std::vector<SurfaceSplineHEZ<Space1D, degree> > build(hpuint continuity, hpreal epsilon = EPSILON) { 
          auto supports = this->template getSupports<Type::VERSION1>();
          return build<degree>(continuity, supports, epsilon);
     }

     //NOTE: Solutions should be const but Eigen refuses to take const data.  See the constrainer's implementation.
     template<hpuint degree>
     std::vector<SurfaceSplineHEZ<Space1D, degree> > build(hpuint continuity, std::vector<std::vector<double> >& solutions, const IndicesArrays& supports, hpreal epsilon = EPSILON) {
          auto s = solutions.begin();
          auto factory = [&](Constrainer<degree>& constrainer) -> std::pair<std::vector<Point1D>, std::vector<hpuint> > { return constrainer.getControlPoints(*(s++)); };
          return build<degree>(continuity, supports, factory, epsilon);
     }

     //NOTE: Solutions should be const but Eigen refuses to take const data.  See the constrainer's implementation.
     template<hpuint degree>//TODO: std::vector<std::vector<double> >& to RealsArrays = Arrays<double>
     std::vector<SurfaceSplineHEZ<Space1D, degree> > build(hpuint continuity, std::vector<std::vector<double> >& solutions, hpreal epsilon = EPSILON) {
          auto supports = this->template getSupports<Type::VERSION1>();
          return build<degree>(continuity, solutions, supports, epsilon);
     }

     template<hpuint degree, class ControlPointFactory>
     std::vector<SurfaceSplineHEZ<Space1D, degree> > build(hpuint continuity, const IndicesArrays& supports, ControlPointFactory factory, hpreal epsilon = EPSILON) {
          std::vector<SurfaceSplineHEZ<Space1D, degree> > surfaces;
          auto constrainer = make_constrainer<degree>(m_refinedProjectiveStructure);
          for(auto begin = supports.begin(), end = supports.end(); begin != end; ++begin) {
               auto support = *begin;
               auto s = support.first;
               auto nTriangles = constrainer.getNumberOfTriangles();
               hpuint i;
               std::cout << "INFO: Support is";
               for(i = 0; i < nTriangles; ++i) {
                    if(i == *s) {
                         std::cout << ' ' << *s;
                         while(i == *s && s != support.second) ++s;//skip duplicates
                         if(s == support.second) {
                              ++i;
                              break;
                         }
                    } else constrainer.zero(i);
               }
               while(i < nTriangles) {
                    constrainer.zero(i);
                    ++i;
               }
               std::cout << ".\n";
               constrainer.constrain(continuity);
               constrainer.solve(epsilon);

               //if(constrainer.existsPositiveSolution()) std::cout << "INFO: There may exist a positive solution.\n";
               //else std::cout << "INFO: No positive solution exists.\n";

               auto& basis = constrainer.getBasis();
               for(auto c = 0l, cols = basis.cols(); c < cols; ++c) {
                    bool positive = true;
                    for(auto r = 0l, rows = basis.rows(); r < rows; ++r) if(!(positive = (basis(r,c) > -epsilon))) break;
                    if(positive) {
                         std::cout << "INFO: " << c << "th basis vector is positive.\n";
                         //for(auto r = 0l, rows = basis.rows(); r < rows; ++r) std::cout << basis(r, c) << '\n';
                         //std::cout << '\n';
                         //auto controlPoints = constrainer.getControlPoints(basis.col(c));
                         //surfaces.push_back(SurfaceSplineHEZ<Space1D, degree>(std::move(controlPoints.first), std::move(controlPoints.second), parameterPoints, mesh.getIndices()));
                    }
               }

               /*Plane plane{{1.0,0.0,0.0}, {1.0,0.0,0.0}};
               for(auto i = 0l, end = basis.cols(); i < end; ++i) {
                    auto controlPoints = constrainer.getControlPoints(basis.col(i));
                    hpreal max = 0.0;
                    for(auto c : controlPoints.first) if(std::abs(c.x) > max) max = std::abs(c.x);
                    for(auto j = controlPoints.first.begin(), end = controlPoints.first.end(); j != end; ++j) *j /= max;
                    SurfaceSplineHEZ<Space1D, degree> surface(std::move(controlPoints.first), std::move(controlPoints.second), parameterPoints, mesh.getIndices());
                    if(auto s = surface.restrict(plane)) {
                         TriangleMesh<VertexP3N> mesh((*s).template toTriangleMesh<VertexP3N, VertexFactory<VertexP3N>, true>());
                         WriterOFF::write(mesh, "/home/pawel/Workspace/col" + std::to_string(i) + ".off");
                    } else std::cerr << "ERROR: Reparametrization failed.\n";
               }*/

               try {
                    auto controlPoints = factory(constrainer);
                    surfaces.push_back(SurfaceSplineHEZ<Space1D, degree>(std::move(controlPoints.first), std::move(controlPoints.second), m_parameterPoints, m_parameterPointIndices));
               } catch(const std::runtime_error& e) { std::cerr << "ERROR: " << e.what() << '\n'; }

               constrainer.unzero();
          }

          return surfaces;
     }

     const ProjectiveStructure3D& getProjectiveStructure() const;

     const ProjectiveStructure3D& getRefinedProjectiveStructure() const;

     const TriangleRefinementScheme& getTriangleRefinementScheme() const;

     template<Type type>
     IndicesArrays getSupports() const {
          IndicesArrays supports;
          for(auto begin = Iterator<0, type>::cbegin(*this), end = Iterator<0, type>::cend(*this); begin != end; ++begin) {
               auto support = *begin;
               supports.push_back(support.begin(), support.end());
          }
          return supports;
          //return std::make_tuple(Iterator<0, type>::cbegin(*this), Iterator<0, type>::cend(*this));//TODO: more efficient to return iterators
     }

     void setTriangleRefinementScheme(const TriangleRefinementScheme& scheme, hpreal epsilon = EPSILON);

private:
     Point3D m_p0, m_p1, m_p2;
     std::vector<hpuint> m_parameterPointIndices;
     std::vector<Point3D> m_parameterPoints;
     const ProjectiveStructure3D& m_projectiveStructure;
     std::reference_wrapper<const TriangleRefinementScheme> m_scheme;
     hpuint m_tetrahedron;
     ProjectiveStructure3D m_refinedProjectiveStructure;

     ProjectiveStructure3D refinedProjectiveStructure(hpreal epsilon = EPSILON);

};//BasisBuilder

}//namespace happah

