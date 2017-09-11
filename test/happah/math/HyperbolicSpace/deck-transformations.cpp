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

void test_compute_simple() {
     hpvec2 p1(0.2, 0.4), q1(-0.3, 0.5), p2(-0.7, -0.4), q2(0.4, -0.2);

     const auto tr21 {hyp_decktrans_P(p2, q2, p1, q1)};
     const auto tr12 {hyp_decktrans_P(p1, q1, p2, q2)};
     //std::cout << "p1=(" << p1.x << "," << p1.y << ")\n";
     //std::cout << "q1=(" << q1.x << "," << q1.y << ")\n";
     //std::cout << "p2=(" << p2.x << "," << p2.y << ")\n";
     //std::cout << "q2=(" << q2.x << "," << q2.y << ")\n";
     //std::cout << "deck tr21ansformation:\n";
     //std::cout << "[" << tr21[0][0] << " " << tr21[1][0] << " " << tr21[2][0] << "]\n";
     //std::cout << "[" << tr21[0][1] << " " << tr21[1][1] << " " << tr21[2][1] << "]\n";
     //std::cout << "[" << tr21[0][2] << " " << tr21[1][2] << " " << tr21[2][2] << "]\n";
     const auto test_id = tr12*tr21;
     for (unsigned i = 0; i < 3; i++) {
          for (unsigned j = 0; j < 3; j++) {
               ASSERT_MSG(std::abs(test_id[j][i] - ((i == j) ? 1.0 : 0.0)) < EPS,
                    "testing whether reciprocal deck transformations combine to the indentity");
          }
     }
}

int main() {
     try {
          test_compute_simple();
     } catch(const std::exception& err) {
          utils::_log_error(std::string("Caught exception: ")+std::string(err.what()));
     }
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
