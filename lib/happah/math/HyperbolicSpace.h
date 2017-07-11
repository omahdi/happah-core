/// \file HyperbolicSpace.h
/// \brief Methods for computations and coordinate transformations in the
/// various models of hyperbolic space
///
/// \author Obada Mahdi <omahdi@gmail.com>
/// \copyright Copyright 2017 Obada Mahdi <omahdi@gmail.com>
/// \copyright Distributed under the Boost Software License, Version 1.0.
/// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#ifndef HYPERBOLIC_SPACE_H
#define HYPERBOLIC_SPACE_H

#include <cmath>
#include <tuple>

namespace happah {

namespace detail {
inline constexpr double hyperbolic_EPS() noexcept { return 1e-6; }
}

/// Convert hyperboloid coordinates to conformal disk coordinates.
inline hpvec2 hyp_HtoC(hpvec3 _co) { return hpvec2(_co.x/(1+_co.z), _co.y/(1+_co.z)); }
/// Convert conformal disk coordinates to hyperboloid coordinates.
inline hpvec3 hyp_CtoH(hpvec2 _co) { return hpvec3(2*_co.x, 2*_co.y, 1 + glm::length2(_co)) / hpreal(1 - glm::length2(_co)); }
/// Convert hyperboloid coordinates to projective disk coordinates.
inline hpvec2 hyp_HtoP(hpvec3 _co) { return hpvec2(_co.x/_co.z, _co.y/_co.z); }
/// Convert projective disk coordinates to hyperboloid coordinates.
inline hpvec3 hyp_PtoH(hpvec2 _co) { return hpvec3(_co.x, _co.y, 1) / hpreal(glm::sqrt(1 - glm::length2(_co))); }
/// Convert conformal disk coordinates to projective disk coordinates.
inline hpvec2 hyp_CtoP(hpvec2 _co) { return hpvec2::value_type(2)*_co / (1 + glm::length2(_co)); }
/// Convert projective disk coordinates to conformal disk coordinates.
inline hpvec2 hyp_PtoC(hpvec2 _co) { return _co / hpvec2::value_type(1 + glm::sqrt(1 - glm::length2(_co))); }
//inline hpvec2 hyp_PtoC(hpvec2 _co) { const auto len2 = glm::length2(_co); return ((1.0 - glm::sqrt(1.0 - len2)) / len2)*_co; }

/**
 * \brief Compute parameters of circle perpendicular to the unit disk and running
 * through points \p p and \p q, which contains the geodesic through both
 * points in the conformal disk model.
 *
 * \param[in] p,q hpvec2
 *
 * \return A \c std::tuple giving center and radius of the circle containing
 * the geodesic arc connecting \p p and \p q.
 *
 * \note Method from http://cs.unm.edu/~joel/NonEuclid/model.html#AppendixA
 * - Let \f$(c_x,c_y)\f$ be the center of the circle with radius \f$r\f$ that
 *   contains the geodesic arc in the conformal disk model connecting
 *   \f$p\f$ and \f$q\f$
 * - Pythagoras:
 *   \f{eqnarray}
        (p_x-c_x)^2 + (p_y-c_y)^2 &=& r^2 \label{eq:1} \\
        (q_x-c_x)^2 + (q_y-c_y)^2 &=& r^2 \label{eq:2}
      \f}
 * - The unit circle and circle \f$C\f$ around \f$(c_x,c_y)\f$ with radius
 *   \f$r\f$ are necessarily orthogonal (hyperbolic geodesics). Thus the
 *   radii from their respective centers to any intersection point form a
 *   right angle, and the line segment connecting their centers the
 *   hypotenuse of a right triangle, again with Pythagoras' theorem giving
 *   \f{equation}
       c_x^2 + c_y^2 = r^2 + 1. \label{eq:3}
     \f}
 * - Combining the formulas:
 *   \f{eqnarray}
       2 p_x c_x + 2 p_y c_y &=& p_x^2 + p_y^2 + 1 \label{eq:4}\\
       2 q_x c_x + 2 q_y c_y &=& q_x^2 + q_y^2 + 1 \label{eq:5}
     \f}
 * - Using Cramer's rule, solve for \f$c_x\f$ and \f$c_y\f$ satisfying both
 *   (4) and (5). This yields
 *   \f{eqnarray}
       c_x &=& (2 (-p_y(q_x^2 + q_y^2 + 1) + q_y(p_x^2 + p_y^2 + 1))) / (4 (p_x q_y - q_x p_y)) \nonumber\\
           &=& (   -p_y(q_x^2 + q_y^2 + 1) + q_y(p_x^2 + p_y^2 + 1))  / (2 (p_x q_y - q_x p_y)) \\
       c_y &=& (2 (-q_x(p_x^2 + p_y^2 + 1) + p_x(q_x^2 + q_y^2 + 1))) / (4 (p_x q_y - q_x p_y)) \nonumber\\
           &=& (   -q_x(p_x^2 + p_y^2 + 1) + p_x(q_x^2 + q_y^2 + 1))  / (2 (p_x q_y - q_x p_y))
     \f}
 *   and the radius can be obtained by solving equation (3) for \f$r\f$.
 */
inline std::tuple<hpvec2, double> hyp_Cgeo(hpvec2 p, hpvec2 q) {
     constexpr auto EPS = detail::hyperbolic_EPS();
     const auto det = p.x*q.y - q.x*p.y;
// special case: p or q are on a diameter, i.e. vectors p, q collinear
     if (glm::abs(det) < EPS)
          return {hpvec2(0.0, 0.0), 0.0};
     const auto lp1 = glm::length2(p) + 1.0, lq1 = glm::length2(q) + 1.0;
     const hpvec2 c {(q.y*lp1 - p.y*lq1) / (2.0*det), (p.x*lq1 - q.x*lp1) / (2.0*det)};
     const hpreal r = glm::sqrt(glm::length2(c) - 1);
     return std::make_tuple(c, r);
}

/// Computes the hyperbolic length of the radius of a circle circumscribing a
/// the regular p-gon of a (p,q)-tesselation (q regular p-gons meeting at each
/// vertex) in the conformal disk model.
///
/// This function uses a formula presented in [1].
///
/// \note [1] Coxeter, "The Trigonometry of Hyperbolic Tesselations", 1997.
/// https://cms.math.ca/cmb/v40/p158,
/// http://dx.doi.org/10.4153/CMB-1997-019-0
inline double hyp_CregularTesselationRadius(int p, int q) {
     assert((p-2)*(q-2) > 4);      // there is no regular tesselation otherwise
     const double s = std::sin(M_PI / p), c = std::cos(M_PI / q);
     return std::acosh(c/s);
}

}    // namespace happah

#endif // #ifdef HYPERBOLIC_SPACE_H
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
