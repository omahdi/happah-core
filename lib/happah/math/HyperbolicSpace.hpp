/// \file HyperbolicSpace.hpp
/// \brief Methods for computations and coordinate transformations in the
/// various models of hyperbolic space
///
/// \author Obada Mahdi <omahdi@gmail.com>
/// \copyright Copyright 2017 Obada Mahdi <omahdi@gmail.com>
/// \copyright Distributed under the Boost Software License, Version 1.0.
/// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <cmath>
#include <numeric>
#include <tuple>
#include <glm/gtx/norm.hpp>

#include <happah/Happah.hpp>

namespace happah {

namespace detail {
inline constexpr double hyperbolic_EPS() noexcept { return 1e-9; }
}

/// Convert hyperboloid coordinates to conformal disk coordinates.
inline hpvec2 hyp_HtoC(hpvec3 _co) { return hpvec2(_co.x/(1+_co.z), _co.y/(1+_co.z)); }
/// Convert conformal disk coordinates to hyperboloid coordinates.
inline hpvec3 hyp_CtoH(hpvec2 _co) { return hpvec3(2*_co.x, 2*_co.y, 1 + glm::length2(_co)) / hpvec3::value_type(1 - glm::length2(_co)); }
/// Convert hyperboloid coordinates to projective disk coordinates.
inline hpvec2 hyp_HtoP(hpvec3 _co) { return hpvec2(_co.x/_co.z, _co.y/_co.z); }
/// Convert projective disk coordinates to hyperboloid coordinates.
inline hpvec3 hyp_PtoH(hpvec2 _co) { return hpvec3(_co.x, _co.y, 1) / hpvec3::value_type(glm::sqrt(1 - glm::length2(_co))); }
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
inline std::tuple<hpvec2, double> hyp_geo_C(hpvec2 p, hpvec2 q) {
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

/// Computes the radius of a circle circumscribing a regular p-gon, centered
/// at the origin, of a (p,q)-tesselation (q regular p-gons meeting at each
/// vertex) in the conformal disk model.
///
/// This function was originally based on ideas presented in [1], which uses a
/// construction that builds a tesselation based on right-angled triangles
/// with acute angles \f$\pi/p\f$ and \f$\pi/q\f$. Since two angles suffice in
/// order to determine the length of the hypothenuse in hyperbolic geometry,
/// we can readily obtain the desired circumradius as
/// \f$\cosh R_{hyp} = \cot(\pi/p)\cot(\pi/q)\f$.
///
/// Note that this is the hyperbolic length, not the Euclidean length in the
/// Poincaré model; the desired length then is \f$tanh(R/2)\f$. By taking
/// advantage of the identity
/// \f$\tanh(\Arcosh(x)/2) = \sqrt{x-1}/\sqrt{x+1}\f$ and further
/// simplification, we arrive at the formula
/// \f$R_c = \sqrt{\cos(\pi/p + \pi/q) / \cos(\pi/p - \pi/q)}\f$.
///
/// \note [1] Coxeter, "The Trigonometry of Hyperbolic Tesselations", 1997.
/// https://cms.math.ca/cmb/v40/p158,
/// http://dx.doi.org/10.4153/CMB-1997-019-0
inline double hyp_tesselation_circumradius_C(int p, int q) {
     assert((p-2)*(q-2) > 4);      // there is no regular tesselation otherwise
     const double A = M_PI/p, B = M_PI/q;
     return std::sqrt(std::cos(A+B)/std::cos(A-B));
}

/// Constructs a convex hyperbolic polygon with given interior angles
/// \p thetas and returns the coordinates of its corner vertices in the
/// conformal disk model.
///
/// The polygon is constructed from quadrilateral pieces that have one corner
/// at the origin, the opposite one matching a polygon corner with a
/// prescribed angle, and two right-angled corners on polygon sides. This
/// method is based on the proof of Theorem 7.16.2 in [1, p. 155f].
///
/// \note [1]: Beardon, "The Geometry of Discrete Groups"
inline std::vector<hpvec2> hyp_polygon_from_angles_C(const std::vector<hpreal>& thetas) {
     using std::begin;
     using std::end;
// be more tolerant for the sum of angles than for intermediate computations
     constexpr double eps = 1e-10;
     constexpr double sum_eps = 1e-6;
     const auto sum_thetas = std::accumulate(begin(thetas), end(thetas), 0.0, [](auto s, auto th) { return s+th; });
     if (sum_thetas - (thetas.size()-2)*M_PI > 0)
          throw std::runtime_error("hyp_polygon_from_angles_C(): There is no hyperbolic polygon with given angles.");

     double t = 1.0;
     const auto n = thetas.size();
// precompute sine and cosine of (theta_i/2) for i=1,...,n
// used as lookup tables in the evaluation of the function g
     struct sincos {
          sincos() = default;
          sincos(double s, double c) : sin(s), cos(c) { }
          double sin {0.0}, cos {1.0};
     };
     std::vector<sincos> halfthetas;
     halfthetas.reserve(n);
     for (auto th : thetas)
          halfthetas.emplace_back(std::sin(double(th)/2), std::cos(double(th)/2));

// "goal function" g(x): sum of alpha_i for i=1,...,n minus pi, given interior
// angles (thetas) and hyperbolic length x of the two sides of the
// quadrilateral patches incident with the origin.
// Note: In Beardon's proof, g is defined without subtracting pi; we define it
// in a way such that g(t) = 0 for the desired solution t.
     auto eval_g = [&] (auto x) {
          const auto ch = std::cosh(x);
          return std::accumulate(begin(halfthetas), end(halfthetas), -M_PI,
               [ch] (auto s, auto ht) { return s + std::asin(ht.cos / ch); });
     };
// first derivative of g(x)
// Recall:
//   d/dx arcsin(x) = 1/sqrt(1-x^2)
// Substitute x = x(t) =  cos(theta_k/2) / cosh(t) (chain rule):
//   d/dx arcsin(x(t)) = (1/sqrt(1-x^2)) * x'(t)
// with x'(t) = d/dx x(t) = -cos(theta_k/2)*sinh(t)/cosh^2(t)
//            = -tanh(t)*(cos(theta_k)/cosh(t))
     auto eval_gprime = [&] (auto x) {
          const auto ch = std::cosh(x);
          const auto ch2 = ch*ch;
          return -std::tanh(x) * std::accumulate(begin(halfthetas), end(halfthetas), 0.0,
               [ch, ch2] (auto s, auto ht) { return s + (1.0 / std::sqrt(ch2/(ht.cos*ht.cos) - 1.0)); });
     };

// Newton's method to find root of g, which is a continuous and
// decreasing function for t > 0.
     auto g = eval_g(t);
     unsigned loopmax = 100;
     while (std::abs(g) > eps && loopmax > 0) {
          const auto gp = eval_gprime(t);
          t -= g / gp;
          g = eval_g(t);
          loopmax--;
     }
     if (loopmax == 0)
          throw std::runtime_error("hyp_polygon_from_angles_C(): maximum number of iterations reached.");
     const double ch = std::cosh(t), sh = std::sinh(t);
     double last_alpha = 0.0, alpha = 0.0, sum_alpha = 0.0;
     std::vector<hpvec2> vertices;
     vertices.reserve(n);
     for (unsigned k = 0; k < n; k++) {
          const auto& ht_k = halfthetas[k];
          last_alpha = alpha;
          alpha = std::asin(ht_k.cos / ch);
          const auto cosech_r = ht_k.sin / sh;
          const auto rk = 1.0 / (std::sqrt(1 + cosech_r*cosech_r) + cosech_r);
          sum_alpha += last_alpha+alpha;
          vertices.emplace_back(rk*std::cos(sum_alpha), rk*std::sin(sum_alpha));
     }
     sum_alpha += alpha;
// Note: this assertion may fail, even though the computations are
// theoretically correct, due to rounding errors. This problem is especially
// prominent when using single-precision computations.
     //std::cout << "sum_alpha - 2*PI = " << (sum_alpha - 2*M_PI) << "\n";
     assert(std::abs(sum_alpha - 2.0*M_PI) < sum_eps);
     return vertices;
}
}    // namespace happah
