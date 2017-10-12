#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <happah/Happah.hpp>
#include <happah/format/off.hpp>

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

template<class _Vertex>
auto
build_minitorus(unsigned _segments, double _height = 1.0, double _r_outer = -1, double _r_inner = -1) {
     using namespace happah;
     using std::sin;
     using std::cos;
     using Point = decltype(std::declval<_Vertex>().position);
     if (_height <= 0)
          throw std::runtime_error("build_minitorus(): height must be positive");
     const double height  = _height;
     // use some sensible defaults
     const double r_outer = _r_outer > 0 ? _r_outer : 2.0*_height;
     const double r_inner = _r_inner > 0 ? _r_inner : _r_outer / std::sqrt(2.0);
     if (r_inner >= r_outer)
          throw std::runtime_error("build_minitorus(): inner radius must be smaller than outer radius");
     const double phi = 2*M_PI / _segments;
     std::vector<_Vertex> vertices;
     std::vector<hpindex> indices;
     if (_segments < 3)
          throw std::runtime_error("build_minitorus(): need at least three segments");
     vertices.reserve(3*_segments);
     indices.reserve(6*_segments);
     for (hpindex s = 0; s < _segments; s++) {
          vertices.emplace_back(Point{r_outer*cos(s*phi), r_outer*sin(s*phi), -height/2.0});
          vertices.emplace_back(Point{r_outer*cos(s*phi), r_outer*sin(s*phi),  height/2.0});
          vertices.emplace_back(Point{r_inner*cos(s*phi), r_inner*sin(s*phi),  0.0});
     }
     const auto num_verts = vertices.size();
     for (hpindex s = 0; s < _segments; s++) {
          const auto offset = 3*s;
          auto push_triplet = [&indices, num_verts, offset](auto u, auto v, auto w) {
               indices.emplace_back((offset + u) % num_verts);
               indices.emplace_back((offset + v) % num_verts);
               indices.emplace_back((offset + w) % num_verts);
          };
          push_triplet(0, 4, 3);
          push_triplet(0, 1, 4);
          push_triplet(1, 5, 4);
          push_triplet(1, 2, 5);
          push_triplet(2, 3, 5);
          push_triplet(2, 0, 3);
     }
     return happah::make_triangle_mesh<_Vertex>(std::move(vertices), std::move(indices));
}

int main() {
     try {
          //auto raw_mesh {format::off::read("minitorus.off")};
          //auto mesh {make_triangle_mesh<VertexP3>(raw_mesh)};
          auto mesh {build_minitorus<VertexP3>(4, 0.5, 1.0, 0.66)};
          format::off::write(mesh, "mt-new-2.off");
          const auto& verts {mesh.getVertices()};
          bool table_format = false;
          auto print_coords = [&table_format, &verts] (auto v0, auto v1, auto v2) {
               auto print_coords = [&table_format] (const auto& v) {
                    if (!table_format)
                         std::cout << " (" << v.position.x << "," << v.position.y << "," << v.position.z << ")";
                    else
                         std::cout << v.position.x << " " << v.position.y << " " << v.position.z << "\n";
               };
               print_coords(verts[v0]);
               print_coords(verts[v1]);
               print_coords(verts[v2]);
               if (!table_format)
                    std::cout << "\n";
          };
          std::cout << "---- TikZ coordinate format ----\n";
          visit_triplets(mesh.getIndices(), print_coords);
          std::cout << "---- table format ----\n";
          table_format = true;
          visit_triplets(mesh.getIndices(), print_coords);
     } catch(const std::exception& err) {
          utils::_log_error(std::string("Caught exception: ")+std::string(err.what()));
     }
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
