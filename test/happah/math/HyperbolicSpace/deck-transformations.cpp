#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <happah/Happah.hpp>
#include <happah/math/HyperbolicSpace.hpp>

// {{{ ---- Logging
namespace utils {
void _log_debug(const std::string& msg) {
     std::cerr << msg << "\n";
}
void _log_error(const std::string& msg) {
// Always flush stderr (via endl)
     std::cerr << msg << std::endl;
}
void _log_output(const std::string& msg) {
     std::cout << msg << "\n";
}
void _log_flush() {
     std::cout << std::flush;
}
}
// }}} ---- Logging
// {{{ ---- Namespace imports
using namespace utils;
using namespace happah;
// }}} ---- Namespace imports
// {{{ ---- Global variables
constexpr double EPS = 1e-9;
// }}} ---- Global variables
// {{{ ---- Test assertions
unsigned g_testcount = 0;
unsigned g_testfail = 0;
#define ASSERT(e)          test_assert(__FILE__, __LINE__, (e), #e);
#define ASSERT_MSG(e, m)   test_assert(__FILE__, __LINE__, (e), (m));
#define ASSERT_EQ(a, b)    test_assert(__FILE__, __LINE__, (a == b), #a " != " #b);
#define SASSERT(e)        test_assert<true>(__FILE__, __LINE__, (e), #e);
#define SASSERT_MSG(e, m) test_assert<true>(__FILE__, __LINE__, (e), (m));
#define SASSERT_EQ(a, b)  test_assert<true>(__FILE__, __LINE__, (a == b), #a " != " #b);
template<bool _strict = false>
void test_assert(std::string fname, unsigned long lineno, bool expr, std::string msg = "") {
     using std::to_string;
     g_testcount++;
     if (!expr) {
          g_testfail++;
          std::string m = fname+':'+to_string(lineno)+
               std::string(": ")+(msg.size() ? msg : std::string("Failed test assertion"));
          utils::_log_error(m);
          if (_strict)
               throw std::runtime_error(m);
     }
}
// }}} ---- Test assertions

template<class Matrix>
bool is_identity(Matrix&& mat) {
// Note: glm-functionality for obtaining the size of matrices/vectors at
// compile-time is not reliable at the moment:
//
// A version prior to version 0.9.9.x introduced "metaprogramming helpers,",
// which are ``constexpr'' values providing the dimensions of vectors and
// matrices. This functionality was available by defining
// "GLM_META_PROG_HELPERS" prior to including glm-related headers. A later
// change, however, removed these constants again, in favor of "constexpr"
// instance functions "length()" of matrix and vector types. We cannot,
// however, make use of these ``constexpr'' functions in glm version 0.9.7,
// for instance, because there is no suitable constexpr constructor and they
// cannot be instantiated at runtime.
//
     //constexpr auto num_cols = std::declval<typename Matrix::row_type>().length();
     //constexpr auto num_rows = std::declval<typename Matrix::col_type>().length();
     const auto num_cols = mat.length();
     const auto num_rows = mat[0].length();
     for (unsigned i = 0; i < num_rows; i++) {
          for (unsigned j = 0; j < num_cols; j++) {
               if ( !(std::abs(mat[j][i] - ((i == j) ? 1.0 : 0.0)) < EPS) )
                    return false;
          }
     }
     return true;
}

template<class Matrix>
void print_matrix(Matrix&& mat) {
     const auto num_cols = mat.length();
     const auto num_rows = mat[0].length();
     std::cout << std::fixed << std::setprecision(6);
     for (unsigned i = 0; i < num_rows; i++) {
          std::cout << "[";
          for (unsigned j = 0; j < num_cols; j++)
               std::cout << ' ' << std::setw(12) << mat[j][i];
          std::cout << " ]\n";
     }
}

void test_compute_simple() {
     const hpvec2 p1(0.2, 0.4), q1(-0.3, 0.5), p2(-0.7, -0.4), q2(0.4, -0.2);
     const auto tr21 {hyp_decktrans_P(p2, q2, p1, q1)};
     const auto tr12 {hyp_decktrans_P(p1, q1, p2, q2)};
     const auto test_id = tr12*tr21;
     print_matrix(test_id);
     ASSERT_MSG(is_identity(test_id), "testing whether reciprocal deck transformations combine to the indentity");
}

int main() {
     try {
          test_compute_simple();
// TODO: more elaborate tests
     } catch(const std::exception& err) {
          utils::_log_error(std::string("Caught exception: ")+std::string(err.what()));
     }
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
