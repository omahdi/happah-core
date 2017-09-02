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
// TODO: More elaborate logging
     std::cerr << msg << "\n";
}
void _log_error(const std::string& msg) {
// TODO: More elaborate logging
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
#define ASSERT(e)        test_assert(__FILE__, __LINE__, (e), #e);
#define ASSERT_MSG(e, m) test_assert(__FILE__, __LINE__, (e), (m));
#define ASSERT_EQ(a, b)  test_assert(__FILE__, __LINE__, (a == b), #a " != " #b);
void test_assert(std::string fname, unsigned long lineno, bool expr, std::string msg = "") {
     using std::to_string;
     g_testcount++;
     if (!expr) {
          g_testfail++;
          utils::_log_error(fname+':'+to_string(lineno)+std::string(": ")+(msg.size() ? msg : std::string("Failed test assertion")));
     }
}
// }}} ---- Test assertions

/// Produces JSON output of polygon as triangle fan, implicitly adding the
/// center vertex at the origin.
void write_json(const std::vector<hpvec2>& vertices) {
     std::cout << "{\n  \"vertices\": [\n";
     std::cout << "    [" << 0.0 << "," << 0.0 << "],\n";
     for (unsigned i = 0, n = vertices.size(); i < n; i++)
          std::cout << "    [" << vertices[i].x << "," << vertices[i].y << (i+1 == n ? "]\n" : "],\n");
     std::cout << "  ],\n";
     std::cout << "  \"faces\": [\n";
     for (unsigned i = 1, n = vertices.size(); i <= n; i++)
          std::cout << "    [0," << i << "," << (1+(i%n)) << (i == n ? "]\n" : "],\n");
     std::cout << "  ]\n}\n";
}

inline double r2d(auto r) { return r*(180.0/M_PI); }
inline double d2r(auto d) { return d*(M_PI/180.0); }

void verify_interior_angles(const std::vector<hpvec2>& vertices, const std::vector<hpreal> thetas) {
// make_tangent: build Euclidean tangent vector at p to circular arc
// representing the geodesic from p to q (not normalized)
     auto make_tangent = [] (const auto& p, const auto& q) {
          using std::get;
          auto geo {hyp_geo_C(p, q)};
          auto dx = q.x - p.x, dy = q.y - p.y;
          if (get<1>(geo) == 0.0)
               return hpvec2(dx, dy);
// compute tangent vector at p
          hpvec2 t {-(p.y - get<0>(geo).y), (p.x - get<0>(geo).x)};
// flip if tangent points into wrong direction
          if (t.x*dx + t.y*dy < 0.0)
               return -t;
          else
               return t;
     };
     for (unsigned k = 0, n = vertices.size(); k < n; k++) {
          const auto& p = vertices[(k+n-1) % n];
          const auto& v = vertices[k];
          const auto& q = vertices[(k+1) % n];
          const auto tp {make_tangent(v, p)}, tq {make_tangent(v, q)};
          const auto dot = glm::dot(tp, tq);
          const double theta {std::acos(dot / std::sqrt(glm::length2(tp)*glm::length2(tq)))};
          //std::cout << "thetas[" << k << "] = " << r2d(thetas[k]) << "; computed theta = " << r2d(theta) << "\n";
          const double err = theta - thetas[k];
          ASSERT(std::abs(err) < EPS);
     }
}

void test_schema_regular(unsigned p, unsigned q) {
// No regular hyperbolic polygon exists unless (p-2)*(q-2) > 4.
     if ((p-2)*(q-2) <= 4)
          throw std::runtime_error("test_schema_regular(): invalid tesselation parameters");
     const double cr = hyp_tesselation_circumradius_C(p, q);
     std::vector<hpvec2> vertices;
     vertices.reserve(p);
     const double phi = (2*M_PI) / p;        // angle at center vertex
     const double theta = (2*M_PI) / q;      // interior angle
     for (unsigned k = 0; k < p; k++)
          vertices.emplace_back(hpvec2(cr*std::cos((phi/2) + k*phi), cr*std::sin((phi/2) + k*phi)));
     std::vector<hpreal> thetas(p, theta);
     verify_interior_angles(vertices, thetas);
// Generate mesh with prescribed angles, which should result in the same set
// of vertices (within a reasonable error range).
     auto gen_vertices {hyp_polygon_from_angles_C(thetas)};
     for (unsigned k = 0; k < p; k++) {
          //std::cout << "vertex #" << k << ": vertices (" << vertices[k].x << "," << vertices[k].y << ") vs gen_vertices (" << gen_vertices[k].x << "," << gen_vertices[k].y << ")\n";
          //std::cout << "vertex #" << k << " error: " << glm::length(vertices[k]-gen_vertices[k]) << "\n";
          ASSERT(glm::length(vertices[k]-gen_vertices[k]) < EPS);
     }
}

void test_schema_4x3_1x4() {
     std::vector<hpreal> thetas = {
          2*M_PI/3, 2*M_PI/3, 2*M_PI/4, 2*M_PI/4, 2*M_PI/3, 2*M_PI/3, 2*M_PI/3, 2*M_PI/3,
          2*M_PI/3, 2*M_PI/3, 2*M_PI/4, 2*M_PI/4, 2*M_PI/3, 2*M_PI/3, 2*M_PI/3, 2*M_PI/3
     };
     const auto vertices {hyp_polygon_from_angles_C(thetas)};
     ASSERT_EQ(vertices.size(), thetas.size());
// dump json mesh
     //std::cout << std::fixed << std::setprecision(5);
     //write_json(vertices);
// Compute interior angles and compare against thetas
     verify_interior_angles(vertices, thetas);
}

void test_schema_random() {
     // TODO: generate random thetas
}

int main() {
     test_schema_regular(8, 8);
     test_schema_regular(5, 4);
     test_schema_regular(18, 3);
     test_schema_4x3_1x4();
     test_schema_random();
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
