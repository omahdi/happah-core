#include "happah/math/nlopt.h"

using namespace happah;

TNLPImpl::TNLPImpl(vector<hpijklr> cc, vector<hpijkr> cq, vector<hpijr> cl, vector<hpir> ck, vector<hpijr> in_A, const vector<hpir>& in_b, Index n, Index m, Index dim_b)
    : cCub{std::move(cc)},
      cQuad{std::move(cq)},
      cLin{std::move(cl)},
      cConst{std::move(ck)},
      A{dim_b, n}, b{dim_b},
      AtA{n, n},
      residuum{dim_b},
      grad_f{n}, g{m},
      jac_g{m, n}, hess{n, n},
      n{n}, m{m}
{
    for (auto& e : in_A) {
        A.value(e.i, e.j) += e.r;
    }
    for (auto& e : in_b) {
        b[e.i] += e.r;
    }

    // compute A^t A
    auto index = A.index();
    for (Index i = 0; i < A.num_rows; ++i) {
        const auto& cols = index[i];
        for (auto& j : cols) {
            for (auto& k : cols) {
                AtA.value(j, k) += A.value(i, j) * A.value(i, k);
            }
        }
    }

    // initialize non-zeros
    vector<Number> zero(max(n, m), 0);
    recompute_temp(zero.data());
    recompute_hess(zero.data(), zero.data(), 0);
    nz_jac = jac_g.rows.size();
    nz_hes = hess.rows.size();
}

void TNLPImpl::recompute_temp(const Number* x) {
    residuum = A.apply(x);
    residuum -= b;
    grad_f = 2 * A.apply_transposed(residuum);

    // set to zero
    g *= 0;
    jac_g.zero();

    // constant part
    for (auto& d : cConst) {
        g[d.i] += d.r;
    }

    // linear part
    for (auto& c : cLin) {
        g[c.i] += c.r * x[c.j];
        jac_g.value(c.i, c.j) += c.r;
    }

    // quadratic part
    for (auto& b : cQuad) {
        g[b.i] += b.r * x[b.j] * x[b.k];
        jac_g.value(b.i, b.j) += b.r * x[b.k];
        jac_g.value(b.i, b.k) += b.r * x[b.j];
    }

    // cubic part
    for (auto& a : cCub) {
        g[a.i] += a.r * x[a.j] * x[a.k] * x[a.l];
        jac_g.value(a.i, a.j) += a.r * x[a.k] * x[a.l];
        jac_g.value(a.i, a.k) += a.r * x[a.i] * x[a.l];
        jac_g.value(a.i, a.l) += a.r * x[a.i] * x[a.k];
    }
}

void TNLPImpl::recompute_hess(const Number *x, const Number *lambda, Number obj_factor) {
    int nz = AtA.rows.size();
    for (int i = 0; i < nz; ++i) {
        hess.value(AtA.rows[i], AtA.cols[i]) = obj_factor * 2 * AtA.values[i];
    }

    for (auto& b : cQuad) {
        hess.value(b.j, b.k) += lambda[b.i] * b.r;
    }

    for (auto& a : cCub) {
        hess.value(a.j, a.k) += lambda[a.i] * a.r * x[a.l];
        hess.value(a.k, a.j) += lambda[a.i] * a.r * x[a.l];
        hess.value(a.j, a.l) += lambda[a.i] * a.r * x[a.k];
        hess.value(a.l, a.j) += lambda[a.i] * a.r * x[a.k];
        hess.value(a.k, a.l) += lambda[a.i] * a.r * x[a.j];
        hess.value(a.l, a.k) += lambda[a.i] * a.r * x[a.j];
    }
}

bool TNLPImpl::get_nlp_info(Index &out_n, Index &out_m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
    out_n = n;
    out_m = m;
    nnz_jac_g = nz_jac;
    nnz_h_lag = nz_hes;
    index_style = C_STYLE;

    return true;
}

bool TNLPImpl::get_bounds_info(Index in_n, Number *x_l, Number *x_u, Index in_m, Number *g_l, Number *g_u) {
    assert(in_n == n);
    assert(in_m == m);
    // unbounded solution
    for (Index i = 0; i < n; ++i) {
        x_l[i] = -1e90;
        x_u[i] =  1e90;
    }
    // equality constraints
    for (Index j = 0; j < m; ++j) {
        g_l[j] = 0;
        g_u[j] = 0;
    }

    return true;
}

bool TNLPImpl::get_starting_point(Index in_n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index in_m, bool init_lambda, Number *lambda) {
    assert(in_n == n);
    assert(in_m == m);
    assert(!init_z);
    assert(!init_lambda);

    if (init_x) {
        // starting point 0
        for (Index i = 0; i < n; ++i) {
            x[i] = 0;
        }
    }

    return true;
}

bool TNLPImpl::eval_f(Index in_n, const Number *x, bool new_x, Number &obj_value) {
    assert(in_n == n);

    if (new_x) {
        recompute_temp(x);
    }

    obj_value = residuum.squaredNorm();

    return true;
}

bool TNLPImpl::eval_grad_f(Index in_n, const Number *x, bool new_x, Number *out_grad_f) {
    assert(in_n == n);

    if (new_x) {
        recompute_temp(x);
    }

    std::copy_n(grad_f.data(), n, out_grad_f);

    return true;
}

bool TNLPImpl::eval_g(Index in_n, const Number *x, bool new_x, Index in_m, Number *out_g) {
    assert(in_n == n);
    assert(in_m == m);

    if (new_x) {
        recompute_temp(x);
    }

    std::copy_n(g.data(), m, out_g);

    return true;
}

bool TNLPImpl::eval_jac_g(Index in_n, const Number *x, bool new_x, Index in_m, Index nele_jac, Index *iRow, Index *jCol, Number *values) {
    assert(in_n == n);
    assert(in_m == m);
    assert(nele_jac == nz_jac);

    if (new_x) {
        recompute_temp(x);
    }

    if (iRow) std::copy(jac_g.rows.begin(), jac_g.rows.end(), iRow);
    if (jCol) std::copy(jac_g.cols.begin(), jac_g.cols.end(), jCol);
    if (values) std::copy(jac_g.values.begin(), jac_g.values.end(), values);

    return true;
}

bool TNLPImpl::eval_h(Index in_n, const Number *x, bool new_x, Number obj_factor, Index in_m, const Number *lambda, bool new_lambda, Index nele_hess, Index *iRow, Index *jCol, Number *values) {
    assert(in_n == n);
    assert(in_m == m);
    assert(nele_hess == nz_hes);

    if (new_x) {
        recompute_temp(x);
    }

    if (new_x || new_lambda) {
        recompute_hess(x, lambda, obj_factor);
    }

    if (iRow) std::copy(hess.rows.begin(), hess.rows.end(), iRow);
    if (jCol) std::copy(hess.cols.begin(), hess.cols.end(), jCol);
    if (values) std::copy(hess.values.begin(), hess.values.end(), values);

    return true;
}

void TNLPImpl::finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L, const Number *z_U, Index m, const Number *g, const Number *lambda, Number obj_value, const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq) {

}
