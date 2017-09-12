#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <experimental/filesystem>
#include <happah/Happah.hpp>
#include <happah/geometry/Vertex.hpp>
#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/NutChain.hpp>
#include <happah/geometry/DiskEmbedding.hpp>
// only for debugging:
#include <happah/format.hpp>

#define PRODUCE_TEST_OUTPUT

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
namespace fs = std::experimental::filesystem;
// }}} ---- Namespace imports
// {{{ ---- Global variables
constexpr double EPS = 1e-6;
// }}} ---- Global variables
// {{{ -- Measuring time
using Duration = std::chrono::microseconds;

template<class TimePoint>
inline auto roundtime(TimePoint& _t) {
     const auto start = _t;
     _t = std::chrono::high_resolution_clock::now();
     return std::chrono::duration_cast<Duration>(_t - start);
}

template<class TimePoint, class Func>
inline auto timecmd(TimePoint& _t, Func&& _cmd) {
     const auto start = _t;
     _cmd();
     _t = std::chrono::high_resolution_clock::now();
     return std::chrono::duration_cast<Duration>(_t - start);
}

template<class Func>
inline auto timecmd(Func&& _cmd) {
     auto start = std::chrono::high_resolution_clock::now();
     return timecmd(start, std::forward<Func>(_cmd));
}

void show_time(std::string _name, Duration _elapsed) {
     constexpr double scale = 1000.0*Duration::period::num/Duration::period::den;
     std::stringstream s;
     s << std::fixed << std::setprecision(3);
     s << "..."  << _elapsed.count()*scale << "ms [" << _name << "]";
     utils::_log_output(s.str());
}

template<class TimePoint>
auto reset_timer(TimePoint& _t) {
     return _t = std::chrono::high_resolution_clock::now();
}
// }}} -- Measuring time
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
     using std::to_string;

     const unsigned p = 8, q = 8;
     const auto fp_mesh {regular_polygon_P(p, q)};
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(fp_mesh, "fp-8_8.off");
#endif
     // TODO: actual tests, verify against manually computed data for a small
     // mesh, e.g. verify transition maps generated for a regular polygon.
}

void test_double_nutchain() {
     using std::to_string;

     auto laptime = std::chrono::high_resolution_clock::now();
     auto t_rst = [&laptime] () { return reset_timer(laptime); };
     auto t_log = [&laptime] (auto name) { show_time(name, roundtime(laptime)); };

     NutChain doubletorus {2, 2.0, 1.0, 1.0, 0.5};
     const auto dt_graph {make_triangle_graph(make_triangle_mesh<VertexP3>(doubletorus))};
     t_log("generating nut chain");

#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(dt_graph, "dt-2nut.off");
#endif
// Manually constructing a cut, see comment in happah/geometry/NutChain.hpp
// regarding the order of generated vertices.
//
//   L       TOP       R     BOTTOM     
// .___._____________.___._____________.
// |   |             |   |             |
// |   |    .___.    |   |    .___.    |
// |   |    |   |    |   |    |   |    |
// |   |    |___|    |   |    |___|    |
// |   |             |   |             |
// *+++*+++++++++++++*+++*+++++++++++++*
// |   |             |   |             |
// |   |      x      |   |      x      |
// |   |             |   |             |
// *+++*+++++++++++++A+++D+++++++++++++*
// |   |      E      |   |             |
// |   |    F___B    |   |    C___.    |
// |   | G  |   |  K |   |    |   |    |
// |   |    H___J    |   |    |___|    |
// |   |      I      |   |             |
// |___|_____________|___|_____________|
//
// (TODO)
//     const auto cut_edges = trim(dt_graph, cut(dt_graph));
//#ifdef PRODUCE_TEST_OUTPUT
//     format::hph::write(cut_edges, fs::path("dt-cut.hph"));
//#endif
     // random cut produced with the above procedure: add here as static data
     const std::string cut_data_hph = "56 131 129 50 48 30 350 366 242 323 321 308 306 305 231 250 275 273 257 339 259 296 290 239 215 300 301 315 316 319 320 262 263 370 89 77 75 147 154 124 127 136 137 46 82 94 386 387 394 54 55 67 164 158 134 132 141";
     const auto cut {format::hph::read<std::vector<hpindex>>(cut_data_hph)};
     std::vector<VertexP3C> cverts;
     cverts.reserve(dt_graph.getNumberOfVertices());
     for (const auto& v : dt_graph.getVertices())
          cverts.emplace_back(v.position, hpcolor(0.4, 0.2, 0.2, 0.5));
     for (auto ei : cut)
          cverts[dt_graph.getEdge(ei).vertex].color = hpcolor(1.0, 0.0, 0.4, 1.0);
     t_rst();
     auto cut_graph {cut_graph_from_edges(dt_graph, cut)};
     t_log("cut_graph_from_edges()");
     utils::_log_output("  Cut graph with "+ to_string(segment_count(cut_graph)) + " and " + to_string(branch_node_count(cut_graph)) + " branch nodes");
     for (unsigned k = 0, n = segment_count(cut_graph); k < n; k++) {
          utils::_log_output("  node #"+to_string(k)+" of degree "+to_string(branch_node_degree(cut_graph, k)));
          const auto segment = cut_segment(cut_graph, k);
          cverts[dt_graph.getEdge(*(segment.end()-1)).vertex].color = hpcolor(0.0, 1.0, 0.4, 1.0);
     }
//
// Remove ``chords'' to prevent degenerate triangles in the boundary mapping.
     t_rst();
     auto has_chords = false; //remove_chords(cut_graph, dt_graph);
     utils::_log_error("Warning: not calling broken remove_chords()");
     t_log("remove_chords()");
     if (has_chords)
          utils::_log_output("Detected and removed chords in cut segments.");
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(make_triangle_mesh(cverts, make_indices(dt_graph)), "dtnut-cut.off");
#endif
//
// Build mesh of fundamental region with a single center vertex.
     t_rst();
     auto fp_mesh {make_fundamental_domain(cut_graph)};
     t_log("make_fundamental_domain()");
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(fp_mesh, "dtnut-fp.off");
#endif
//
// Remap disk boundary to boundary our fundamental region.
// - build disk mesh with duplicated vertices on the boundary
     t_rst();
     auto disk_result {cut_to_disk(cut_graph, dt_graph)};
     t_log("cut_to_disk()");
     auto& disk_mesh = std::get<0>(disk_result);
// - obtain cached info on boundary edges (used by make_projective_structure)
     const auto& boundary_edge_info = std::get<1>(disk_result);
// - equidistant distribution along boundary of fundamental region
//   FIXME: paired sides are treated separately, which works for now because
//   of the equidistant distribution!
     t_rst();
     set_boundary_constraints(disk_mesh, boundary_edge_info,
          [&fp_mesh, &cut_graph](auto ei, auto side, auto offset) {
               const auto t = double(offset) / segment_length(cut_graph, side);
               const auto& v0_ref = fp_mesh.getVertex(side, 1);
               const auto& v1_ref = fp_mesh.getVertex(side, 2);
               // project onto hyperboloid before interpolating?
               //const auto v0 {hyp_PtoH(hpvec2(v0_ref.position.x, v0_ref.position.y))};
               //const auto v1 {hyp_PtoH(hpvec2(v1_ref.position.x, v1_ref.position.y))};
               const auto v0 {hpvec2(v0_ref.position.x, v0_ref.position.y)};
               const auto v1 {hpvec2(v1_ref.position.x, v1_ref.position.y)};
               const auto x = (1.0-t)*v0.x + t*v1.x, y = (1.0-t)*v0.y + t*v1.y, z = (1.0-t)*v0.z + t*v1.z;
               return Point2D{x/z, y/z};
          });
     t_log("set_boundary_constraints()");
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(disk_mesh, "dtnut-disk-raw.off");
#endif
// - compute Tutte embedding based on mean-value coordinates
     t_rst();
     compute_disk_embedding(disk_mesh, make_mv_coord_generator(dt_graph));
     t_log("compute_disk_embedding()");
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(disk_mesh, "dtnut-disk.off");
#endif
// Compute projective transition maps
     t_rst();
     auto disk_ps = make_projective_structure(cut_graph, boundary_edge_info, disk_mesh, fp_mesh);
     t_log("make_projective_structure()");
#ifdef PRODUCE_TEST_OUTPUT
     // TODO: define operator<<() for ProjectiveStructure
     //format::hph::write(disk_ps, fs::path("dtnut-ps.hph"));
#endif
     auto sorted_cut {cut_edges(cut_graph)};
     std::sort(std::begin(sorted_cut), std::end(sorted_cut));
     const hpindex seed_triangle = 0;
     const auto seed_t {disk_mesh.getTriangle(seed_triangle)};
     const auto& seed_v0_pos = std::get<0>(seed_t).position;
     const auto& seed_v1_pos = std::get<1>(seed_t).position;
     const auto& seed_v2_pos = std::get<2>(seed_t).position;
     auto recons_mesh {make_triangle_mesh(disk_ps, sorted_cut, seed_triangle,
          Point3D(seed_v0_pos.x, seed_v0_pos.y, 1.0),
          Point3D(seed_v1_pos.x, seed_v1_pos.y, 1.0),
          Point3D(seed_v2_pos.x, seed_v2_pos.y, 1.0),
          VertexFactory<VertexP3>())};
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(recons_mesh, "dtnut-recons.off");
#endif
     for (auto& v : recons_mesh.getVertices()) {
          v.position.x /= v.position.z;
          v.position.y /= v.position.z;
          v.position.z = 1.0;
     }
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(recons_mesh, "dtnut-recons-proj.off");
#endif
}

void test_minitorus() {
     using std::to_string;

     auto laptime = std::chrono::high_resolution_clock::now();
     auto t_rst = [&laptime] () { return reset_timer(laptime); };
     auto t_log = [&laptime] (auto name) { show_time(name, roundtime(laptime)); };

     auto raw_mesh {format::off::read("minitorus.off")};
     const auto dt_graph {make_triangle_graph<VertexP3>(make_triangle_mesh<VertexP3>(raw_mesh))};
     t_log("reading minitorus");

     // random cut produced with the above procedure: add here as static data
     std::vector<hpindex> cut_vertices {0,1,2,0,2,5,8,11,2};
     std::vector<hpindex> cut_offsets {0, 4, 9};
     std::vector<hpindex> cut;
     //TODO cut graph from paths
     for (unsigned k = 0, n = cut_offsets.size()-1; k < n; k++) {
          const hpindex first = cut_offsets[k], last = cut_offsets[k+1];
          const auto len = last-first;
          auto last_v = cut_vertices[first];
          for (unsigned j = first+1; j < last; j++) {
               const auto maybe_e = make_edge_index(dt_graph, last_v, cut_vertices[j == last-1 ? first : j]);
               ASSERT_MSG(!!maybe_e, "consecutive vertices in cut must be connected by an edge");
               last_v = cut_vertices[j];
               cut.emplace_back(*maybe_e);
          }
     }
     std::vector<VertexP3C> cverts;
     cverts.reserve(dt_graph.getNumberOfVertices());
     for (const auto& v : dt_graph.getVertices())
          cverts.emplace_back(v.position, hpcolor(0.4, 0.2, 0.2, 0.5));
     for (auto ei : cut)
          cverts[dt_graph.getEdge(ei).vertex].color = hpcolor(1.0, 0.0, 0.4, 1.0);
     t_rst();
     auto cut_graph {cut_graph_from_edges(dt_graph, cut)};
     t_log("cut_graph_from_edges()");
     utils::_log_output("  Cut graph with "+ to_string(segment_count(cut_graph)) + " and " + to_string(branch_node_count(cut_graph)) + " branch nodes");
     for (unsigned k = 0, n = segment_count(cut_graph); k < n; k++) {
          utils::_log_output("  node #"+to_string(k)+" of degree "+to_string(branch_node_degree(cut_graph, k)));
          const auto segment = cut_segment(cut_graph, k);
          cverts[dt_graph.getEdge(*(segment.end()-1)).vertex].color = hpcolor(0.0, 1.0, 0.4, 1.0);
     }
//
// Remove ``chords'' to prevent degenerate triangles in the boundary mapping.
     t_rst();
     auto has_chords = false; //remove_chords(cut_graph, dt_graph);
     utils::_log_error("Warning: not calling broken remove_chords()");
     t_log("remove_chords()");
     if (has_chords)
          utils::_log_output("Detected and removed chords in cut segments.");
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(make_triangle_mesh(cverts, make_indices(dt_graph)), "minit-cut.off");
#endif
//
// Build mesh of fundamental region with a single center vertex.
     t_rst();
     auto fp_mesh {make_fundamental_domain(cut_graph)};
     t_log("make_fundamental_domain()");
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(fp_mesh, "minit-fp.off");
#endif
//
// Remap disk boundary to boundary our fundamental region.
// - build disk mesh with duplicated vertices on the boundary
     t_rst();
     auto disk_result {cut_to_disk(cut_graph, dt_graph)};
     t_log("cut_to_disk()");
     auto& disk_mesh = std::get<0>(disk_result);
// - obtain cached info on boundary edges (used by make_projective_structure)
     const auto& boundary_edge_info = std::get<1>(disk_result);
// - equidistant distribution along boundary of fundamental region
//   FIXME: paired sides are treated separately, which works for now because
//   of the equidistant distribution!
     t_rst();
     set_boundary_constraints(disk_mesh, boundary_edge_info,
          [&fp_mesh, &cut_graph](auto ei, auto side, auto offset) {
               const auto t = double(offset) / segment_length(cut_graph, side);
               const auto& v0_ref = fp_mesh.getVertex(side, 1);
               const auto& v1_ref = fp_mesh.getVertex(side, 2);
               // project onto hyperboloid before interpolating?
               //const auto v0 {hyp_PtoH(hpvec2(v0_ref.position.x, v0_ref.position.y))};
               //const auto v1 {hyp_PtoH(hpvec2(v1_ref.position.x, v1_ref.position.y))};
               const auto v0 {hpvec2(v0_ref.position.x, v0_ref.position.y)};
               const auto v1 {hpvec2(v1_ref.position.x, v1_ref.position.y)};
               const auto x = (1.0-t)*v0.x + t*v1.x, y = (1.0-t)*v0.y + t*v1.y, z = (1.0-t)*v0.z + t*v1.z;
               return Point2D{x/z, y/z};
          });
     t_log("set_boundary_constraints()");
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(disk_mesh, "minit-disk-raw.off");
#endif
// - compute Tutte embedding based on mean-value coordinates
     t_rst();
     compute_disk_embedding(disk_mesh, make_mv_coord_generator(dt_graph));
     t_log("compute_disk_embedding()");
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(disk_mesh, "minit-disk.off");
#endif
// Compute projective transition maps
     t_rst();
     auto disk_ps = make_projective_structure(cut_graph, boundary_edge_info, disk_mesh, fp_mesh);
     t_log("make_projective_structure()");
#ifdef PRODUCE_TEST_OUTPUT
     // TODO: define operator<<() for ProjectiveStructure
     //format::hph::write(disk_ps, fs::path("minit-ps.hph"));
#endif
     auto sorted_cut {cut_edges(cut_graph)};
     std::sort(std::begin(sorted_cut), std::end(sorted_cut));
     const hpindex seed_triangle = 0;
     const auto seed_t {disk_mesh.getTriangle(seed_triangle)};
     const auto& seed_v0_pos = std::get<0>(seed_t).position;
     const auto& seed_v1_pos = std::get<1>(seed_t).position;
     const auto& seed_v2_pos = std::get<2>(seed_t).position;
     auto recons_mesh {make_triangle_mesh(disk_ps, sorted_cut, seed_triangle,
          Point3D(seed_v0_pos.x, seed_v0_pos.y, 1.0),
          Point3D(seed_v1_pos.x, seed_v1_pos.y, 1.0),
          Point3D(seed_v2_pos.x, seed_v2_pos.y, 1.0),
          VertexFactory<VertexP3>())};
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(recons_mesh, "minit-recons.off");
#endif
     for (auto& v : recons_mesh.getVertices()) {
          v.position.x /= v.position.z;
          v.position.y /= v.position.z;
          v.position.z = 1.0;
     }
#ifdef PRODUCE_TEST_OUTPUT
     format::off::write(recons_mesh, "minit-recons-proj.off");
#endif
}

int main() {
     try {
          test_regular_8gon();
          test_double_nutchain();
          test_minitorus();
     } catch(const std::exception& err) {
          utils::_log_error(std::string("Caught exception: ")+std::string(err.what()));
     }
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
