#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <happah/Happah.hpp>
#include <happah/geometry/Vertex.hpp>
#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/NutChain.hpp>
#include <happah/geometry/DiskEmbedding.hpp>
// only for debugging:
#include <happah/format.hpp>

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
constexpr double EPS = 1e-6;
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

template<bool _conformal>
TriangleMesh<VertexP2> regular_polygon(unsigned p, unsigned q) {
     using std::begin;
     using std::end;
     const double phi = (2*M_PI) / p;        // angle at center vertex
     const double theta = (2*M_PI) / q;      // interior angle
     const double cr = _conformal ? hyp_tesselation_circumradius_C(p, q) : hyp_tesselation_circumradius_P(p, q);
     std::vector<hpvec2> vpos;
     for (unsigned k = 0; k < p; k++)
          vpos.emplace_back(hpvec2(cr*std::cos((phi/2) + k*phi), cr*std::sin((phi/2) + k*phi)));
     std::vector<VertexP2> vertices;
     vertices.reserve(p+1);
     std::vector<hpindex> indices;
     indices.reserve(p*3);
     vertices.emplace_back(hpvec2(0, 0));
     for (auto& vp : vpos)
          vertices.emplace_back(std::move(vp));
     for (unsigned i = 1; i <= p; i++) {
          indices.emplace_back(0);
          indices.emplace_back(i);
          indices.emplace_back(1+(i % p));
     }
     return make_triangle_mesh(vertices, indices);
}

TriangleMesh<VertexP2> regular_polygon_C(unsigned p, unsigned q) {
     return regular_polygon<true>(p, q);
}

TriangleMesh<VertexP2> regular_polygon_P(unsigned p, unsigned q) {
     return regular_polygon<false>(p, q);
}

void test_regular_8gon() {
     const unsigned p = 8, q = 8;
     const auto fp_mesh {regular_polygon_P(p, q)};
     //format::off::write(fp_mesh, "fp-8_8.off");
     NutChain doubletorus {2, 2.0, 1.0, 1.0, 0.5};
     const auto dt_mesh {make_triangle_mesh<VertexP3>(doubletorus)};
     format::off::write(dt_mesh, "dt-2nut.off");
}

int main() {
     try {
          test_regular_8gon();

     } catch(const std::exception& err) {
          utils::_log_error(std::string("Caught exception: ")+std::string(err.what()));
     }
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
