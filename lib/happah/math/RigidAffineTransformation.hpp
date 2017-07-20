// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <memory>
#include <type_traits>

#include "happah/Happah.hpp"
#include "happah/math/Space.hpp"

namespace RigidAffineTransformationTraits {
     template<class Space>
     struct conversion;
}

template<class Space>
class RigidAffineTransformation {
     static_assert(is_space<Space>::value, "A rigid affine transformation can only be parameterized by a space.");
     typedef typename std::enable_if<Space::DIMENSION == 2 || Space::DIMENSION == 3, typename std::conditional<Space::DIMENSION == 2, hpmat2x2, hpmat3x3>::type>::type Matrix;
     typedef typename Space::VECTOR Vector;

public:
     typedef Matrix MATRIX;
     typedef Space SPACE;

     static const RigidAffineTransformation<Space> IDENTITY;

     static RigidAffineTransformation<Space3D> to3D(const RigidAffineTransformation<Space2D>& transformation) { 
          typedef typename RigidAffineTransformation<Space3D>::MATRIX Matrix;
          typedef typename RigidAffineTransformation<Space3D>::SPACE::VECTOR Vector;
          return RigidAffineTransformation<Space3D>(Matrix(transformation.m_matrix[0][0],transformation.m_matrix[0][1],0.0,transformation.m_matrix[1][0],transformation.m_matrix[1][1],0.0,0.0,0.0,0.0),Vector(transformation.m_translation,0.0)); }

     RigidAffineTransformation()
          : m_matrix(1.0), m_translation(0.0) {}
     //TODO: change matrix to rotation quaternion; otherwise cannot ensure matrix is rigid
     RigidAffineTransformation(const Matrix& matrix, const Vector& translation)
          : m_matrix(matrix), m_translation(translation) {}
     RigidAffineTransformation(const RigidAffineTransformation& transformation)
          : m_matrix(transformation.m_matrix), m_translation(transformation.m_translation) {}
     ~RigidAffineTransformation() {}

     RigidAffineTransformation<Space> compose(const RigidAffineTransformation<Space>& transformation) const { return RigidAffineTransformation<Space>(m_matrix * transformation.m_matrix, m_matrix * transformation.m_translation + m_translation); }
     const Matrix& getMatrix() const { return m_matrix; }
     const Vector& getTranslation() const { return m_translation; }
     RigidAffineTransformation<Space> inverse() const {
          Matrix matrix = glm::inverse(m_matrix);
          return RigidAffineTransformation<Space>(matrix, -(matrix * m_translation));
     }
     void setMatrix(const Matrix& matrix) { m_matrix = matrix; }
     void setTranslation(const Vector& translation) { m_translation = translation; }
     hpmat4x4 toMatrix4x4() const { return RigidAffineTransformationTraits::conversion<Space>::toMatrix4x4(*this); }

     RigidAffineTransformation<Space>& operator=(const RigidAffineTransformation<Space>& other) {
          if(this != &other) {
               m_matrix = other.m_matrix;
               m_translation = other.m_translation;
          }
          return *this;
     }

     friend class RigidAffineTransformationTraits::conversion<Space>;

protected:
     Matrix m_matrix;
     Vector m_translation;

};
typedef RigidAffineTransformation<Space2D> RigidAffineTransformation2D;
typedef std::shared_ptr<RigidAffineTransformation2D> RigidAffineTransformation2D_ptr;
typedef RigidAffineTransformation<Space3D> RigidAffineTransformation3D;
typedef std::shared_ptr<RigidAffineTransformation3D> RigidAffineTransformation3D_ptr;

template<class Space>
const RigidAffineTransformation<Space> RigidAffineTransformation<Space>::IDENTITY(Matrix(1.0), Vector(0.0));

namespace RigidAffineTransformationTraits {
     template<>
     struct conversion<Space3D> {
          static hpmat4x4 toMatrix4x4(const RigidAffineTransformation<Space3D>& transformation) {
               return hpmat4x4(transformation.m_matrix[0][0],transformation.m_matrix[0][1],transformation.m_matrix[0][2],0.0f,transformation.m_matrix[1][0],transformation.m_matrix[1][1],transformation.m_matrix[1][2],0.0f,transformation.m_matrix[2][0],transformation.m_matrix[2][1],transformation.m_matrix[2][2],0.0f,transformation.m_translation.x,transformation.m_translation.y,transformation.m_translation.z,1.0f);
          }
     };

     //NOTE: There is intentionally no definition of conversion in 2D space because the conversion to a 4x4 matrix is not well-defined in 2D space.
}

