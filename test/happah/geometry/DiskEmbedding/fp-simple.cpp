#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <experimental/filesystem>
#include <boost/dynamic_bitset.hpp>
#include <happah/Happah.hpp>
#include <happah/geometry/Vertex.hpp>
#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/NutChain.hpp>
#include <happah/geometry/DiskEmbedding.hpp>
// only for debugging:
#include <happah/format.hpp>

#define PRODUCE_TEST_OUTPUT
#define SHOW_BRANCH_NODES

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
using namespace std::string_literals;
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

// write_off(): Wrappers for format::off::write(); purpose is to simplify
// switching to 3D coordinate output for easy viewing in MeshLab.
template<class Vertex, class = std::enable_if_t<std::is_base_of<VertexP<Space2D>, Vertex>::value>>
std::vector<VertexP3> vert2to3(const std::vector<Vertex>& vertices, char = 0) {
     std::vector<VertexP3> result;
     result.reserve(vertices.size());
     for (const auto& v : vertices)
          result.emplace_back(Space3D::POINT(v.position.x, v.position.y, 1.0));
     return result;
}
template<class Vertex, class = std::enable_if_t<std::is_base_of<VertexPC<Space2D>, Vertex>::value>>
std::vector<VertexPC<Space3D>> vert2to3(const std::vector<Vertex>& vertices, int = 0) {
     std::vector<VertexPC<Space3D>> result;
     result.reserve(vertices.size());
     for (const auto& v : vertices)
          result.emplace_back(Space3D::POINT(v.position.x, v.position.y, 1.0), v.color);
     return result;
}
void write_off(const TriangleMesh<VertexP3>& mesh, const std::string& filename) {
     format::off::write(mesh, filename);
}
void write_off(const TriangleMesh<VertexP3C>& mesh, const std::string& filename) {
     format::off::write(mesh, filename);
}
void write_off(const TriangleGraph<VertexP3>& mesh, const std::string& filename) {
     format::off::write(mesh, filename);
}
void write_off(const TriangleGraph<VertexP3C>& mesh, const std::string& filename) {
     format::off::write(mesh, filename);
}
template<class Vertex, class = std::enable_if_t<std::is_base_of<VertexP<Space2D>, Vertex>::value>>
void write_off(const TriangleMesh<Vertex>& mesh, const std::string& filename, char = 0) {
     auto vertices {vert2to3(mesh.getVertices())};
     format::off::write(make_triangle_mesh(vertices, mesh.getIndices()), filename);
}
template<class Vertex, class = std::enable_if_t<std::is_base_of<VertexPC<Space2D>, Vertex>::value>>
void write_off(const TriangleMesh<Vertex>& mesh, const std::string& filename, int = 0) {
     auto vertices {vert2to3(mesh.getVertices())};
     format::off::write(make_triangle_mesh(vertices, mesh.getIndices()), filename);
}
template<class Vertex, class = std::enable_if_t<std::is_base_of<VertexP<Space2D>, Vertex>::value>>
void write_off(const TriangleGraph<Vertex>& mesh, const std::string& filename, char = 0) {
     auto vertices {vert2to3(mesh.getVertices())};
     format::off::write(make_triangle_mesh(vertices, make_indices(mesh)), filename);
}
template<class Vertex, class = std::enable_if_t<std::is_base_of<VertexPC<Space2D>, Vertex>::value>>
void write_off(const TriangleGraph<Vertex>& mesh, const std::string& filename, int = 0) {
     auto vertices {vert2to3(mesh.getVertices())};
     format::off::write(make_triangle_mesh(vertices, make_indices(mesh)), filename);
}

void test_regular_8gon() {
     using std::to_string;

     const unsigned p = 8, q = 8;
     const auto fp_mesh {regular_polygon_P(p, q)};
#ifdef PRODUCE_TEST_OUTPUT
     write_off(fp_mesh, "fp-8_8.off");
#endif
     // TODO: actual tests, verify against manually computed data for a small
     // mesh, e.g. verify transition maps generated for a regular polygon.
}

template<class Vertex>
void test_embedding(const TriangleGraph<Vertex>& src_mesh, const CutGraph& cut_graph, const TriangleMesh<VertexP2>& fp_mesh, const std::string& file_prefix) {
     using std::to_string;

     auto laptime = std::chrono::high_resolution_clock::now();
     auto t_rst = [&laptime] () { return reset_timer(laptime); };
     auto t_log = [&laptime] (auto name) { show_time(name, roundtime(laptime)); };

     utils::_log_output("test_embedding(file_prefix="+file_prefix+")");
// Remap disk boundary to boundary our fundamental region.
// - build disk mesh with duplicated vertices on the boundary
     t_rst();
     auto disk_result {cut_to_disk(cut_graph, src_mesh)};
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
#ifdef TRANSITION_MAPS_USE_HYPERBOLOID
               const auto v0 {hyp_PtoH(hpvec2(v0_ref.position.x, v0_ref.position.y))};
               const auto v1 {hyp_PtoH(hpvec2(v1_ref.position.x, v1_ref.position.y))};
#else
               const auto v0 {hpvec3(v0_ref.position.x, v0_ref.position.y, 1.0)};
               const auto v1 {hpvec3(v1_ref.position.x, v1_ref.position.y, 1.0)};
#endif
               const auto x = (1.0-t)*v0.x + t*v1.x, y = (1.0-t)*v0.y + t*v1.y, z = (1.0-t)*v0.z + t*v1.z;
               return Point2D{x/z, y/z};
          });
     t_log("set_boundary_constraints()");
#ifdef PRODUCE_TEST_OUTPUT
     write_off(disk_mesh, file_prefix+"-disk-raw.off");
#endif
// - compute Tutte embedding based on mean-value coordinates
     t_rst();
     compute_disk_embedding(disk_mesh, make_mv_coord_generator(src_mesh));
     t_log("compute_disk_embedding()");
#ifdef PRODUCE_TEST_OUTPUT
     write_off(disk_mesh, file_prefix+"-disk.off");
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
#ifdef TRANSITION_MAPS_USE_HYPERBOLOID
          hyp_PtoH(seed_v0_pos), hyp_PtoH(seed_v1_pos), hyp_PtoH(seed_v2_pos),
#else
          Point3D(seed_v0_pos.x, seed_v0_pos.y, 1.0),
          Point3D(seed_v1_pos.x, seed_v1_pos.y, 1.0),
          Point3D(seed_v2_pos.x, seed_v2_pos.y, 1.0),
#endif
          VertexFactory<VertexP3>())};
#ifdef PRODUCE_TEST_OUTPUT
     write_off(recons_mesh, file_prefix+"-recons.off");
#endif
#if defined(PRODUCE_TEST_OUTPUT) && defined(TRANSITION_MAPS_USE_HYPERBOLOID)
     for (auto& v : recons_mesh.getVertices()) {
          v.position.x /= v.position.z;
          v.position.y /= v.position.z;
          v.position.z = 1.0;
     }
     write_off(recons_mesh, file_prefix+"-recons-proj.off");
#endif
#ifdef PRODUCE_TEST_OUTPUT
     {
          std::vector<VertexP3> vertices;
          vertices.reserve(disk_mesh.getNumberOfVertices());
          for (const auto &v : disk_mesh.getVertices())
               vertices.emplace_back(hyp_PtoH(v.position));
          write_off(make_triangle_mesh(vertices, make_indices(disk_mesh)), file_prefix+"-disk-hyp.off");
     }
#endif
}

void test_double_nutchain() {
     using std::to_string;

     auto laptime = std::chrono::high_resolution_clock::now();
     auto t_rst = [&laptime] () { return reset_timer(laptime); };
     auto t_log = [&laptime] (auto name) { show_time(name, roundtime(laptime)); };
     const std::string file_prefix = "dtnut";

     NutChain doubletorus {2, 2.0, 1.0, 1.0, 0.5};
     const auto src_mesh {make_triangle_graph(make_triangle_mesh<VertexP3>(doubletorus))};
     t_log("generating nut chain");

#ifdef PRODUCE_TEST_OUTPUT
     write_off(src_mesh, "dt-2nut.off");
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
     const auto the_cut = trim(src_mesh, cut(src_mesh));
#ifdef PRODUCE_TEST_OUTPUT
     format::hph::write(the_cut, fs::path("dt-cut.hph"));
#endif
     // random cut produced with the above procedure: add here as static data
     //const std::string cut_data_hph = "56 131 129 50 48 30 350 366 242 323 321 308 306 305 231 250 275 273 257 339 259 296 290 239 215 300 301 315 316 319 320 262 263 370 89 77 75 147 154 124 127 136 137 46 82 94 386 387 394 54 55 67 164 158 134 132 141";
     //const auto the_cut {format::hph::read<std::vector<hpindex>>(cut_data_hph)};
     std::vector<VertexP3C> cverts;
     cverts.reserve(src_mesh.getNumberOfVertices());
     for (const auto& v : src_mesh.getVertices())
          cverts.emplace_back(v.position, hpcolor(0.4, 0.2, 0.2, 0.5));
     for (auto ei : the_cut)
          cverts[src_mesh.getEdge(ei).vertex].color = hpcolor(1.0, 0.0, 0.4, 1.0);
     t_rst();
     auto cut_graph {cut_graph_from_edges(src_mesh, the_cut)};
     t_log("cut_graph_from_edges()");
// Remove ``chords'' to prevent degenerate triangles in the boundary mapping.
     t_rst();
     const auto has_chords = remove_chords(cut_graph, src_mesh);
     t_log("remove_chords() [has_chords="s + to_string(has_chords) + "]"s);
     if (has_chords)
          utils::_log_output("Detected and removed chords in cut segments.");
     utils::_log_output("  Cut graph with "+ to_string(segment_count(cut_graph)) + " and " + to_string(branch_node_count(cut_graph)) + " branch nodes");
#ifdef SHOW_BRANCH_NODES
     for (unsigned k = 0, n = segment_count(cut_graph); k < n; k++) {
          utils::_log_output("  node #"+to_string(k)+" of degree "+to_string(branch_node_degree(cut_graph, k)));
          const auto segment = cut_segment(cut_graph, k);
          cverts[src_mesh.getEdge(*(segment.end()-1)).vertex].color = hpcolor(0.0, 1.0, 0.4, 1.0);
     }
#endif
#ifdef PRODUCE_TEST_OUTPUT
     write_off(make_triangle_mesh(cverts, make_indices(src_mesh)), file_prefix+"-cut.off");
#endif
//
// Build mesh of fundamental region with a single center vertex.
     t_rst();
     auto fp_mesh {make_fundamental_domain(cut_graph)};
     t_log("make_fundamental_domain()");
#ifdef PRODUCE_TEST_OUTPUT
     write_off(fp_mesh, file_prefix+"-fp.off");
#endif
     test_embedding(src_mesh, cut_graph, fp_mesh, file_prefix);
}

void test_minitorus() {
     using std::to_string;

     auto laptime = std::chrono::high_resolution_clock::now();
     auto t_rst = [&laptime] () { return reset_timer(laptime); };
     auto t_log = [&laptime] (auto name) { show_time(name, roundtime(laptime)); };
     const std::string file_prefix = "minit"s;

     auto raw_mesh {format::off::read("minitorus.off")};
     const auto src_mesh {make_triangle_graph<VertexP3>(make_triangle_mesh<VertexP3>(raw_mesh))};
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
               const auto maybe_e = make_edge_index(src_mesh, last_v, cut_vertices[j == last-1 ? first : j]);
               ASSERT_MSG(!!maybe_e, "consecutive vertices in cut must be connected by an edge");
               last_v = cut_vertices[j];
               cut.emplace_back(*maybe_e);
          }
     }
     std::vector<VertexP3C> cverts;
     cverts.reserve(src_mesh.getNumberOfVertices());
     for (const auto& v : src_mesh.getVertices())
          cverts.emplace_back(v.position, hpcolor(0.4, 0.2, 0.2, 0.5));
     for (auto ei : cut)
          cverts[src_mesh.getEdge(ei).vertex].color = hpcolor(1.0, 0.0, 0.4, 1.0);
     t_rst();
     auto cut_graph {cut_graph_from_edges(src_mesh, cut)};
     t_log("cut_graph_from_edges()");
     utils::_log_output("  Cut graph with "+ to_string(segment_count(cut_graph)) + " and " + to_string(branch_node_count(cut_graph)) + " branch nodes");
#ifdef SHOW_BRANCH_NODES
     for (unsigned k = 0, n = segment_count(cut_graph); k < n; k++) {
          utils::_log_output("  node #"+to_string(k)+" of degree "+to_string(branch_node_degree(cut_graph, k)));
          const auto segment = cut_segment(cut_graph, k);
          cverts[src_mesh.getEdge(*(segment.end()-1)).vertex].color = hpcolor(0.0, 1.0, 0.4, 1.0);
     }
#endif
//
// Remove ``chords'' to prevent degenerate triangles in the boundary mapping.
     t_rst();
     auto has_chords = remove_chords(cut_graph, src_mesh);
     t_log("remove_chords()");
     if (has_chords)
          utils::_log_output("Detected and removed chords in cut segments.");
#ifdef PRODUCE_TEST_OUTPUT
     write_off(make_triangle_mesh(cverts, make_indices(src_mesh)), file_prefix+"-cut.off");
#endif
//
// Build mesh of fundamental region with a single center vertex.
// Note: A torus does not admit a hyperbolic structure. There is no hyperbolic
// rectangle with four right angles, as its area would vanish. To test at
// least the embedding procedure, we create a rectangle with interior angles
// less than pi/2.
     t_rst();
     //const auto fp_corners {hyp_polygon_from_angles_P(std::vector<double>(4, (85.0/180.0)*M_PI))};
     const auto fp_corners {make_convex_polygon(std::vector<hpreal>(4, (85.0/180.0)*M_PI))};
     const std::vector<VertexP2> fp_vertices {{Point2D(0.0, 0.0),
          hyp_CtoP(fp_corners[0]), hyp_CtoP(fp_corners[1]),
          hyp_CtoP(fp_corners[2]), hyp_CtoP(fp_corners[3])}};
     //fp_vertices.reserve(4);
     //for (const auto& p : fp_corners)
     //     fp_vertices.emplace_back(p);
     auto fp_mesh {make_triangle_mesh(fp_vertices, std::vector<hpindex>{0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1})};
     t_log("make_fundamental_domain()");
#ifdef PRODUCE_TEST_OUTPUT
     write_off(fp_mesh, file_prefix+"-fp.off");
#endif
     test_embedding(src_mesh, cut_graph, fp_mesh, file_prefix);
}

void test_double_torus() {
     using std::to_string;

     auto laptime = std::chrono::high_resolution_clock::now();
     auto t_rst = [&laptime] () { return reset_timer(laptime); };
     auto t_log = [&laptime] (auto name) { show_time(name, roundtime(laptime)); };
     const std::string file_prefix = "dtorus"s;

     auto raw_mesh {format::off::read("double-torus.off")};
     const auto src_mesh {make_triangle_graph<VertexP3>(make_triangle_mesh<VertexP3>(raw_mesh))};
     t_log("reading double torus");

#if 0
     std::vector<hpindex> cut;
     cut.reserve(cut_edge_flags.size());
     for (hpindex ei = 0, n = src_mesh.getNumberOfEdges(); ei < n; ei++)
          if (cut_edge_flags[ei])
               cut.emplace_back(ei);
     std::vector<VertexP3C> cverts;
     cverts.reserve(src_mesh.getNumberOfVertices());
     for (const auto& v : src_mesh.getVertices())
          cverts.emplace_back(v.position, hpcolor(0.4, 0.2, 0.2, 0.5));
     for (auto ei : cut)
          cverts[src_mesh.getEdge(ei).vertex].color = hpcolor(1.0, 0.0, 0.4, 1.0);
     t_rst();
     auto cut_graph {cut_graph_from_edges(src_mesh, cut)};
#else
     auto cut_graph {format::hph::read<CutGraph>(fs::path("dt-cut-graph.hph"))};
     const auto cut {cut_edges(cut_graph)};
     std::vector<VertexP3C> cverts;
     cverts.reserve(src_mesh.getNumberOfVertices());
     for (const auto& v : src_mesh.getVertices())
          cverts.emplace_back(v.position, hpcolor(0.4, 0.2, 0.2, 0.5));
     for (auto ei : cut)
          cverts[src_mesh.getEdge(ei).vertex].color = hpcolor(1.0, 0.0, 0.4, 1.0);
#endif
     t_log("cut_graph_from_edges()");
     utils::_log_output("  Cut graph with "+ to_string(segment_count(cut_graph)) + " and " + to_string(branch_node_count(cut_graph)) + " branch nodes");
#ifdef SHOW_BRANCH_NODES
     for (unsigned k = 0, n = segment_count(cut_graph); k < n; k++) {
          utils::_log_output("  node #"+to_string(k)+" of degree "+to_string(branch_node_degree(cut_graph, k)));
          const auto segment = cut_segment(cut_graph, k);
          cverts[src_mesh.getEdge(*(segment.end()-1)).vertex].color = hpcolor(0.0, 1.0, 0.4, 1.0);
     }
#endif
//
// Remove ``chords'' to prevent degenerate triangles in the boundary mapping.
     t_rst();
     const auto has_chords = remove_chords(cut_graph, src_mesh);
     t_log("remove_chords()");
     if (has_chords)
          utils::_log_output("Detected and removed chords in cut segments.");
#ifdef PRODUCE_TEST_OUTPUT
     write_off(make_triangle_mesh(cverts, make_indices(src_mesh)), file_prefix+"-cut.off");
#endif
//
// Build mesh of fundamental region with a single center vertex.
     t_rst();
     auto fp_mesh {make_fundamental_domain(cut_graph)};
     t_log("make_fundamental_domain()");
#ifdef PRODUCE_TEST_OUTPUT
     write_off(fp_mesh, file_prefix+"-fp.off");
#endif
     test_embedding(src_mesh, cut_graph, fp_mesh, file_prefix);
}

int main() {
     try {
          std::cout << "---- test_regular_8gon ----\n";
          test_regular_8gon();
     } catch(const std::exception& err) {
          utils::_log_error("Caught exception: "s + err.what());
          g_testfail++;
     }
     try {
          std::cout << "---- test_minitorus ----\n";
          test_minitorus();
     } catch(const std::exception& err) {
          utils::_log_error("Caught exception: "s + err.what());
          g_testfail++;
     }
     try {
          std::cout << "---- test_double_nutchain ----\n";
          test_double_nutchain();
     } catch(const std::exception& err) {
          utils::_log_error("Caught exception: "s + err.what());
          g_testfail++;
     }
     try {
          std::cout << "---- test_double_torus ----\n";
          test_double_torus();
     } catch(const std::exception& err) {
          utils::_log_error("Caught exception: "s + err.what());
          g_testfail++;
     }
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
