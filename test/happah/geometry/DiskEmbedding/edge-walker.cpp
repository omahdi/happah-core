#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <boost/format.hpp>
#include <happah/Happah.hpp>
#include <happah/geometry/Vertex.hpp>
#include <happah/geometry/TriangleMesh.hpp>
#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/DiskEmbedding.hpp>

#define TEST_DEBUG_OUTPUT

// {{{ ---- Logging
namespace utils {
void _log_debug(const std::string& msg) {
#ifdef TEST_DEBUG_OUTPUT
     std::cerr << msg << "\n";
#endif
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
#define ASSERT(e)         test_assert(__FILE__, __LINE__, (e), #e);
#define ASSERT_MSG(e, m)  test_assert(__FILE__, __LINE__, (e), (m));
#define ASSERT_EQ(a, b)   test_assert(__FILE__, __LINE__, (a == b), #a " != " #b);
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

template<class _Vertex>
auto
build_tetrahedron(double _side = 1.0) {
     using namespace happah;
     if (_side <= 0)
          throw std::runtime_error("build_tetrahedron(): side length must be positive");
     return happah::make_triangle_mesh<_Vertex>(
          std::vector<_Vertex>{{{0.0, 0.0, 0.0}}, {{_side, 0.0, 0.0}}, {{0.0, 0.0, _side}}, {{0.0, _side, 0.0}}},
          std::vector<hpindex>{0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2}
     );
}

void test_edge_walker_basic() {
     using namespace happah;
     using namespace utils;
     auto mesh {make_triangle_graph(build_tetrahedron<VertexP<Space3D>>())};
     auto walker = [&mesh]() { return make_edge_walker(mesh); };

     for (unsigned i = 0; i < 4; i++) {
          auto w {walker().uv(0, 1) + i};
          _log_debug(str(boost::format("edge w [#%d]: %d e, %d opp, %d u, %d v, %d next_e, %d prev_e")
               % i % w.e() % w.opp() % w.u() % w.v() % (w+1).e() % (w-1).e()));
     }
     // Sanity check: we assume that the first face is (0, 1, 2) and edges 0, 1
     // and 2 correspond to pairs (0, 1), (1, 2), and (2, 0).
     {
          const auto& ind {make_indices(mesh)};
          ASSERT_EQ(ind[0], 0u);
          ASSERT_EQ(ind[1], 1u);
          ASSERT_EQ(ind[2], 2u);
          ASSERT_EQ(0u, *make_edge_index(mesh, 0, 1));
          ASSERT_EQ(1u, *make_edge_index(mesh, 1, 2));
          ASSERT_EQ(2u, *make_edge_index(mesh, 2, 0));
     }

     const auto w {walker().uv(0, 1)};
     {
          auto walker_copy_construct {w};
          ASSERT_EQ(w, walker_copy_construct);
          auto walker_copy_assign = w;
          ASSERT_EQ(w, walker_copy_assign);
     }
     {
          ASSERT_EQ(w.opp(), (-w).e());
          ASSERT_EQ(w, -(-w));
          ASSERT_EQ(w, w+3);
          ASSERT_EQ(walker().e(1), walker().uv(1, 2));
     }
     {
          auto x {walker().e(0)};
          ASSERT_EQ(x.e(), 0u); ASSERT_EQ(x.u(), 0u); ASSERT_EQ(x.v(), 1u); ASSERT_EQ(x.opp(), 8u);
          ++x;
          ASSERT_EQ(x.e(), 1u); ASSERT_EQ(x.u(), 1u); ASSERT_EQ(x.v(), 2u); ASSERT_EQ(x.opp(), 11u);
          ++x;
          ASSERT_EQ(x.e(), 2u); ASSERT_EQ(x.u(), 2u); ASSERT_EQ(x.v(), 0u); ASSERT_EQ(x.opp(), 3u);
          ++x; ASSERT_EQ(x.e(), 0u);
          --x; ASSERT_EQ(x.e(), 2u);
          --x; ASSERT_EQ(x.e(), 1u);
          --x; ASSERT_EQ(x.e(), 0u);
     }
     {
          auto x {walker().uv(1, 2)};
          ASSERT_EQ(x.e(), 1u); ASSERT_EQ(x.u(), 1u); ASSERT_EQ(x.v(), 2u);
          auto y {x};
          y.flip();
          ASSERT_EQ(-x, y);
          y.flip();
          ASSERT_EQ(x, y);
     }
}

int main() {
     try {
          test_edge_walker_basic();
     } catch(const std::exception& err) {
          utils::_log_error(std::string("Caught exception: ")+std::string(err.what()));
     }
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
