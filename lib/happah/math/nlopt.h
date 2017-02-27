#ifndef NLOPT_H
#define NLOPT_H

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <coin/IpTNLP.hpp>
#include <vector>
#include <unordered_map>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <boost/functional/hash.hpp>

#include "happah/Happah.h"

namespace happah {

using namespace Ipopt;

using std::vector;
using std::max;
using SparseMatrix = Eigen::SparseMatrix<Number, Eigen::ColMajor, Index>;
using Vector = Eigen::Matrix<Number, 1, Eigen::Dynamic>;
using Entry = std::pair<Index, Index>;
using EntryMap = std::unordered_map<Entry, Index, boost::hash<Entry>>;

struct MySparseMatrix {
    vector<Number> values;
    vector<Index> rows;
    vector<Index> cols;
    EntryMap map;
    Index num_rows;
    Index num_cols;

    Number& value(Index row, Index col) {
        assert(row >= 0 && row < num_rows);
        assert(col >= 0 && col < num_cols);
        Entry e{row, col};
        auto it = map.find(e);
        if (it == map.end()) {
            it = map.insert(make_pair(e, values.size())).first;
            rows.push_back(row);
            cols.push_back(col);
            values.push_back(0);
        }
        return values[(*it).second];
    }

    void zero() {
        for (auto& v : values) v = 0;
    }

    vector<vector<Index>> index() const {
        vector<vector<Index>> index{num_rows, vector<Index>{}};
        for (size_t i = 0; i < rows.size(); ++i) {
            index[rows[i]].push_back(cols[i]);
        }
        return index;
    }

    Vector apply(const Number* x) {
        Vector result{num_rows};
        for (Index i = 0; i < rows.size(); ++i) {
            result[rows[i]] += values[i] * x[cols[i]];
        }
        return result;
    }

    Vector apply_transposed(const Vector& b) {
        assert(b.size() == num_rows);
        Vector result{num_cols};
        for (Index i = 0; i < rows.size(); ++i) {
            result[cols[i]] += values[i] * b[rows[i]];
        }
        return result;
    }

    MySparseMatrix(Index num_rows, Index num_cols)
        : num_rows{num_rows}, num_cols{num_cols} {}
};

struct MySymmetricSparseMatrix : MySparseMatrix {
private:
    Number temp;

public:
    MySymmetricSparseMatrix(Index rows, Index cols) : MySparseMatrix{rows, cols} {
        assert(rows == cols);
    }

    Number& value(Index row, Index col) {
        if (row > col) return temp;
        return MySparseMatrix::value(row, col);
    }
};

class TNLPImpl : public TNLP  {
private:
    // Meta information for the NLP:
    Index n; // Number of variables
    Index m; // Number of constraints
    Index nz_jac; // Number of non-zero entries of the Jacobian
    Index nz_hes; // Number of non-zero entries of the Hessian of the Lagrangian

    // Cubic constraints
    vector<hpijklr> cCub;
    vector<hpijkr>  cQuad;
    vector<hpijr>   cLin;
    vector<hpir>    cConst;

    // Objective function || Ax-b ||
    MySparseMatrix A;
    Vector b;

    // Auxillary values
    MySymmetricSparseMatrix AtA; // A^tA for computation of the Hessian

    // Temporary values
    Vector residuum; // Ax - b
    Vector grad_f;
    Vector g;
    MySparseMatrix jac_g;
    MySymmetricSparseMatrix hess;

    void recompute_temp(const Number* x);
    void recompute_hess(const Number* x, const Number* lambda, Number obj_factor);

public:
    TNLPImpl(vector<hpijklr> cc, vector<hpijkr> cq, vector<hpijr> cl, vector<hpir> ck, vector<hpijr> in_A, const vector<hpir>& in_b, Index n, Index m, Index dim_b);
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style);
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u);
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda);
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values);
    virtual bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values);
    virtual void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq);
};

}
#endif // NLOPT_H
