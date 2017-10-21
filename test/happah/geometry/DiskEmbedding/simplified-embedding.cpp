#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <queue>

#include <glm/glm.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <happah/Happah.hpp>
#include <happah/util/VertexFactory.hpp>
#include <happah/geometry/NutChain.hpp>
#include <happah/geometry/TriangleGraph.hpp>
#include <happah/geometry/DiskEmbedding.hpp>
#include <happah/format.hpp>

#ifndef LOG_DEBUG
#define LOG_DEBUG { (void)0; }
#endif

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
// {{{ ---- Measuring time
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
// }}} ---- Measuring time

namespace test_embedding {
/// Builds coefficient matrices for computing Tutte embeddings.
///
/// The \p disk_mesh argument provides adjacency information, and lookup
/// tables \p v_inner and \p v_boundary are used to map vertex indices of the
/// disk mesh to contiguous integer intervals \f$[0..n)\f$ and \f$[0..m)\f$,
/// where \f$n\f$ denotes the number of interior vertices and \f$m\f$ the
/// number of boundary vertices, respectively; they are used as corresponding
/// row and column indices in the computed matrices. Barycentric coordinates
/// like mean-value coordinates are computed by calling \p coord_builder,
/// provided by the caller.
///
/// The function handle \p coord_builder is called with signature
/// \code{.cpp}
///   double coord_builder(hpindex ei, hpindex v, hpindex w);
/// \endcode
/// where \p v and \p w are vertex indices of \p disk_mesh corresponding to
/// edge \p ei.  The return value is expected to be the computed
/// barycentric/convex-combination "coordinate" of \p v with respect to \p w
/// (but need not be normalized; normalization is carried out automatically).
///
/// \returns A tuple of two matrices: the first matrix is quadratic and
/// represents coefficients for all inner vertices; the second one is
/// rectangular with the same number of rows as the first (number of inner
/// vertices) and one column per boundary vertex. Coefficients correspond to
/// the \e normalized coordinates; this matrix can then be used to compute the
/// right-hand vector when solving for a particular coordinate component.
template<class Vertex, class Coords>
auto
compute_embedding_coeff( // {{{1
     const TriangleGraph<Vertex>& disk_mesh,
     std::unordered_map<hpindex, hpindex> v_inner,
     std::unordered_map<hpindex, hpindex> v_boundary,
     Coords&& coord_builder
) {
     using InnerMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
     using BoundaryMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
     const auto num_inner = v_inner.size(), num_boundary = v_boundary.size();
// Count inner and boundary neighbors for each inner vertex
     using Eigen::ArrayXi;
     ArrayXi deg_inner {ArrayXi::Constant(num_inner, 1)};   // one on the diagonal
     ArrayXi deg_boundary {ArrayXi::Zero(num_inner)};
     for (const auto& iv : v_inner) {
          visit_spokes(make_spokes_enumerator(disk_mesh.getEdges(), disk_mesh.getOutgoing(iv.first)),
               [&] (auto ei) {
                    const auto v = disk_mesh.getEdge(ei).vertex;
                    if (v_boundary.count(v) == 0)
                         deg_inner(iv.second)++;
                    else
                         deg_boundary(iv.second)++;
               });
     }
// Use degrees to accurately reserve memory for sparse matrices.
     InnerMatrix coeff_inner(num_inner, num_inner);
     coeff_inner.reserve(deg_inner);
     BoundaryMatrix coeff_boundary(num_inner, num_boundary);
     coeff_boundary.reserve(deg_boundary);
     auto col_less = [](auto a, auto b) { return a.col() < b.col(); };

// Create reverse lookup table mapping row indices back to inner vertex IDs.
     std::vector<hpindex> v_inner_rev(num_inner, 0);
     for (const auto& iv : v_inner) {
          assert(iv.second < num_inner);
          v_inner_rev[iv.second] = iv.first;
     }
// Loop over inner vertices again and populate sparse matrix coefficients
// - inner neighbors contribute coefficients for coeff_inner, which is the
//   matrix describing the linear system to be solved for the unknown coords
// - boundary neighbors contribute coefficients for coeff_boundary, which is
//   used to compute the right-hand side for each coordinate component
     using Triplet = Eigen::Triplet<double>;
     using std::begin;
     using std::end;
     const auto& disk_vertices = disk_mesh.getVertices();
     for (int i = 0; i < int(num_inner); i++) {
          auto v = v_inner_rev[i];
          std::vector<Triplet> inner_row, boundary_row;
          inner_row.reserve(deg_inner(i));
          boundary_row.reserve(deg_boundary(i));
          double row_sum = 0.0;
          visit_spokes(make_spokes_enumerator(disk_mesh.getEdges(), disk_mesh.getOutgoing(v)),
               [&] (auto ei) {
                    const auto w = disk_mesh.getEdge(ei).vertex;
// vw_coeff are the mean-value coordinates computed for the vertex pair (v,w)
                    double vw_coeff = coord_builder(ei, v, w);
                    row_sum += vw_coeff;
                    if (v_boundary.count(w) == 0)
                         inner_row.emplace_back(i, v_inner[w], -vw_coeff);
                    else
                         boundary_row.emplace_back(i, v_boundary[w], vw_coeff);
               });
          inner_row.emplace_back(i, i, 1.0);
          std::sort(begin(inner_row), end(inner_row), col_less);
          std::sort(begin(boundary_row), end(boundary_row), col_less);
          for (const auto t : inner_row) {
               if (t.col() != i)
                    coeff_inner.insert(t.row(), t.col()) = t.value() / row_sum;
               else
                    coeff_inner.insert(t.row(), t.col()) = t.value();      // insert 1.0
          }
          for (const auto t : boundary_row)
               coeff_boundary.insert(t.row(), t.col())  = t.value() / row_sum;
     }
     coeff_inner.makeCompressed();
     coeff_boundary.makeCompressed();
     return std::make_tuple(std::move(coeff_inner), std::move(coeff_boundary));
}
// }}}1 compute_embedding_coeff()

/// Compute a generalized Tutte embedding for \p disk_mesh, with the boundary
/// identified by edges with indices greater than
/// <tt>3*disk_mesh.getNumberOfTriangles()</tt>.
///
/// \param[in] disk_mesh disk topology mesh with fixed boundary
/// \param[in] coord_builder callable object for computing coordinates
///
/// The \c coord_builder() is called with an edge ID followed by the
/// corresponding ordered pair <tt>(v,w)</tt> of vertex indices and is
/// expected to return the coordinate of \c v with respect to its neighbor
/// \c w. Note that the latter vertex indices are with respect to
/// \p disk_mesh. We let the caller deal with the details of figuring out how
/// to compute coordinates and how to relate edge IDs to vertices in the
/// original (closed) mesh.
template<class Vertex, class Coords>
void
compute_disk_embedding(TriangleGraph<Vertex>& disk_mesh, Coords&& coord_builder) {   // {{{1
     using namespace std::string_literals;

     auto laptime = std::chrono::high_resolution_clock::now();
     auto t_rst = [&laptime] () { return laptime = std::chrono::high_resolution_clock::now(); };
     auto t_log = [&laptime] (auto name) { show_time(name, roundtime(laptime)); };

     auto& disk_vertices {disk_mesh.getVertices()};
     const auto num_verts = disk_vertices.size();
     const auto num_faces = disk_mesh.getNumberOfTriangles();
     const auto num_edges = disk_mesh.getEdges().size();
     std::unordered_map<hpindex, hpindex> v_inner, v_boundary;
     for (hpindex ei = 3*num_faces; ei < num_edges; ei++) {
          const auto& e = disk_mesh.getEdge(ei);
          v_boundary.emplace(e.vertex, v_boundary.size());
     }
     for (hpindex vi = 0; vi < num_verts; vi++)
          if (v_boundary.count(vi) == 0)
               v_inner.emplace(vi, v_inner.size());
     const auto num_boundary = v_boundary.size();

     t_log("compute_disk_embedding:prepare");
     t_rst();
     using std::get;
     auto coeff_mat {compute_embedding_coeff(disk_mesh, v_inner, v_boundary, std::forward<Coords>(coord_builder))};
     //Eigen::MatrixX2d boundary_coords {Eigen::MatrixX2d::Zero(num_boundary, 2)};
     Eigen::MatrixX2d boundary_coords {num_boundary, 2};
     for (const auto& bv : v_boundary) {
          const auto& vpos {disk_vertices[bv.first].position};
          boundary_coords(bv.second, 0) = vpos.x;
          boundary_coords(bv.second, 1) = vpos.y;
     }
     Eigen::MatrixX2d rhs {get<1>(coeff_mat)*boundary_coords};
     t_log("compute_disk_embedding:coefficients");
     t_rst();
// Create solver with column-major storage (required!)
     Eigen::SparseLU<Eigen::SparseMatrix<double>> tutte_solver;
     tutte_solver.compute(get<0>(coeff_mat));
     t_log("...compute_disk_embedding:compute_step");
     if (tutte_solver.info() != Eigen::Success)
          throw std::runtime_error("tutte_solver reported error computing LU decomposition: "s +
               tutte_solver.lastErrorMessage());
     Eigen::MatrixX2d coords {tutte_solver.solve(rhs)};
     t_log("...compute_disk_embedding:solve_step");
     if (tutte_solver.info() != Eigen::Success)
          throw std::runtime_error("tutte_solver reported error trying to solve for coordinates"s +
               tutte_solver.lastErrorMessage());
// Transfer computed coordinates of inner vertices
     for (const auto& iv : v_inner) {
          auto& v {disk_vertices[iv.first]};
          v.position.x = coords(iv.second, 0);
          v.position.y = coords(iv.second, 1);
     }
     t_log("...compute_disk_embedding:setup_output_coords");
}
// }}}1 compute_disk_embedding()
}    // namespace test_embedding

/// Variant of cut() that picks next edge by hop distance.
/// Note: not fit for meshes with border!
Indices hopdist_cut(const std::vector<Edge>& edges, hpindex t0 = 0) { // {{{1
     struct edge_info {
          hpindex ei, dist;
          // Note: invert meaning to turn max-heap into min-heap
          bool operator<(const edge_info& other) const { return dist > other.dist; }
     };
     auto cache = boost::dynamic_bitset<>(edges.size(), false);
     std::priority_queue<edge_info> pq;
     assert(t0 < edges.size() / 3);
     for (hpindex ei = 3*t0, last_ei = 3*t0+2; ei <= last_ei; ei++) {
          pq.push(edge_info{ei, 0});
          cache[ei] = true;
     }
     return basic_cut(edges, t0, [&](auto& neighbors) {
          //for(auto e : boost::irange(0u, hpindex(mesh.getEdges().size())))
          //     if(neighbors[e << 1] != std::numeric_limits<hpindex>::max() && neighbors[mesh.getEdge(e).opposite << 1] == std::numeric_limits<hpindex>::max()) return e;
          edge_info info;
          hpindex e;
          while (!pq.empty()) {
               info = pq.top();
               e = info.ei;
               if (cache[e])
                    break;
               pq.pop();
          }
          if (pq.empty())
               return std::numeric_limits<hpindex>::max();
          pq.pop();
          auto& edge = edges[edges[e].opposite];
          auto e0 = edges[edge.previous].opposite;
          auto e1 = edges[edge.next].opposite;
          auto b0 = neighbors[e0 << 1] == std::numeric_limits<hpindex>::max();
          auto b1 = neighbors[e1 << 1] == std::numeric_limits<hpindex>::max();
          if (b0) {
               pq.push(edge_info{edge.previous, info.dist+1});
               cache[edge.previous] = true;
          } else
               cache[e0] = false;
          if (b1) {
               pq.push(edge_info{edge.next, info.dist+1});
               cache[edge.next] = true;
          } else
               cache[e1] = false;
          cache[e] = false;
          return hpindex(e);
     });
}    // }}}1

/// Variant of cut() that picks next edge based on smallest curvature.
/// Note: not fit for meshes with border!
template<class Vertex>
Indices curvdist_cut(const TriangleGraph<Vertex>& graph, hpindex t0 = 0) { // {{{1
     const auto& edges = graph.getEdges();
     struct edge_info {
          hpindex ei;
          hpreal dist;
          unsigned hops;
          // Note: invert meaning to turn max-heap into min-heap
          bool operator<(const edge_info& other) const {
               //return (dist - other.dist) < 0.1 || hops > other.hops;
               //return dist - other.dist > 0;
               return dist - other.dist < 0;
          }
     };
     auto cache = boost::dynamic_bitset<>(edges.size(), false);
     std::priority_queue<edge_info> pq;
     assert(t0 < edges.size() / 3);
     for (hpindex ei = 3*t0, last_ei = 3*t0+2; ei <= last_ei; ei++) {
          pq.push(edge_info{ei, hpreal{0}, 0});
          cache[ei] = true;
     }
     // computes angle in radians at p0 between (p1-p0) and (p2-p0)
     const auto edge_angle = [](hpvec3 p0, hpvec3 p1, hpvec3 p2) {
          hpvec3 u {p1-p0}, v {p2-p0};
          return std::acos(glm::dot(u, v) / (glm::length(u)*glm::length(v)));
     };
     const auto triarea = [](hpvec3 p0, hpvec3 p1, hpvec3 p2) {
          return 0.5*glm::length(glm::cross(p1-p0, p2-p0));
     };
     const auto heron = [](hpvec3 p0, hpvec3 p1, hpvec3 p2) {
          hpvec3 u {p1-p0}, v {p2-p1}, w {p0-p2};
          hpreal lu {glm::length(u)}, lv {glm::length(v)}, lw {glm::length(w)};
          hpreal s = 0.5*(lu+lv+lw);
          return std::sqrt(s*(s-lu)*(s-lv)*(s-lw));
     };
     const auto gaussc = [&](hpindex v0) {
          const auto pos_v0 = graph.getVertex(v0).position;
          hpreal sum_area = 0.0, sum_angles = 0.0;
          visit_spokes(make_spokes_enumerator(edges, graph.getOutgoing(v0)),
               [&] (auto ei) {
                    const auto& e {graph.getEdge(ei)};
                    const auto& pos_v1 = graph.getVertex(e.vertex).position;
                    const auto& pos_v2 = graph.getVertex(graph.getEdge(e.next).vertex).position;
                    sum_area += triarea(pos_v0, pos_v1, pos_v2);
                    sum_angles += edge_angle(pos_v0, pos_v1, pos_v2);
               });
          return 3*(2.0*M_PI - sum_angles) / sum_area;
          //return 2.0*M_PI - sum_angles;      // only deficit
     };
     const hpindex num_vertices = graph.getNumberOfVertices();
     const hpindex num_triangles = graph.getNumberOfTriangles();
     std::vector<hpreal> angular_deficit;
     angular_deficit.reserve(num_vertices);
     for (hpindex vi = 0; vi < num_vertices; vi++)
          angular_deficit.push_back(gaussc(vi));
     std::vector<hpreal> avg_deficit;
     avg_deficit.reserve(num_triangles);
     // Simple heuristic for judging the "flatness" of a triangle in terms of
     // the angle deficit of its defining vertices.
     for (hpindex ti = 0; ti < num_triangles; ti++) {
          const auto& e0 = graph.getEdge(3*ti + 0);
          const auto& e1 = graph.getEdge(3*ti + 1);
          const auto& e2 = graph.getEdge(3*ti + 2);
          avg_deficit.push_back((angular_deficit[e2.vertex] + angular_deficit[e0.vertex] + angular_deficit[e1.vertex])/3);
          //avg_deficit.push_back(std::max(std::max(angular_deficit[e2.vertex], angular_deficit[e0.vertex]), angular_deficit[e1.vertex]));
          //avg_deficit.push_back(std::min(std::min(angular_deficit[e2.vertex], angular_deficit[e0.vertex]), angular_deficit[e1.vertex]));
     }

     return basic_cut(edges, t0, [&](auto& neighbors) {
          //for(auto e : boost::irange(0u, hpindex(mesh.getEdges().size())))
          //     if(neighbors[e << 1] != std::numeric_limits<hpindex>::max() && neighbors[mesh.getEdge(e).opposite << 1] == std::numeric_limits<hpindex>::max()) return e;
          edge_info info;
          hpindex e;
          while (!pq.empty()) {
               info = pq.top();
               e = info.ei;
               if (cache[e])
                    break;
               pq.pop();
          }
          if (pq.empty())
               return std::numeric_limits<hpindex>::max();
          pq.pop();
          auto& edge = edges[edges[e].opposite];
          auto e0 = edges[edge.previous].opposite;
          auto e1 = edges[edge.next].opposite;
          auto b0 = neighbors[e0 << 1] == std::numeric_limits<hpindex>::max();
          auto b1 = neighbors[e1 << 1] == std::numeric_limits<hpindex>::max();
          if (b0) {
               auto next_tri = make_triangle_index(graph.getEdge(edge.previous).opposite);
               pq.push(edge_info{edge.previous, avg_deficit[next_tri], info.hops+1});
               cache[edge.previous] = true;
          } else
               cache[e0] = false;
          if (b1) {
               auto next_tri = make_triangle_index(graph.getEdge(edge.next).opposite);
               pq.push(edge_info{edge.next, avg_deficit[next_tri], info.hops+1});
               cache[edge.next] = true;
          } else
               cache[e1] = false;
          cache[e] = false;
          return hpindex(e);
     });
}    // }}}1

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

void test_minitorus_embedding() { // {{{1
     const auto raw_mesh {format::off::read("mt-new-2.off")};
     const auto mesh {make_triangle_graph(make_triangle_mesh<VertexP3>(raw_mesh))};
     //const auto the_cut {trim(nut_mesh, hopdist_cut(nut_mesh.getEdges()))};
     auto the_cut {format::hph::read<std::vector<hpindex>>(p("minit-cut-edges.hph"))};
     auto cut_graph {cut_graph_from_edges(mesh, the_cut)};
     remove_chords(cut_graph, mesh);
     {
          auto xedges = cut_edges(cut_graph);
          format::hph::write(xedges, p("mt-cut-edges.hph"));
     }
     auto disk_result {cut_to_disk(cut_graph, mesh, CutGraph::VertexMapping::CONTIGUOUS_BOUNDARY)};
     //auto disk_result {cut_to_disk(cut_graph, mesh, CutGraph::VertexMapping::KEEP_INNER)};
     auto& disk = std::get<0>(disk_result);
     //auto& boundary_info = std::get<1>(disk_result);
     std::cout << "---- minitorus disk graph ----\n";
     auto walker = make_edge_walker(disk);
     //visit_triplets(begin(disk.getEdges()), disk.getNumberOfTriangles(), 3,
     //  [&walker] (const auto& e0, const auto& e1, const auto& e2) {
     //  });
     for (unsigned i = 3*disk.getNumberOfTriangles(), num_edges = disk.getEdges().size(); i < num_edges; i++) {
          walker.e(i);
          std::cout << "edge[" << i << "]: " << walker.u() << " -> " << walker.v() << "\n";
     }

     using namespace std::string_literals;
     using std::to_string;
     const auto num_faces = disk.getNumberOfTriangles();
     const auto num_edges = disk.getEdges().size();
     // Map boundary vertices onto unit circle
     const auto bi_start = 3*num_faces;
     const auto b_length = num_edges - bi_start;
     for (hpindex ei = bi_start; ei < num_edges; ei++) {
          //hpindex bc_last;       // last edge index of current boundary component
          //for (bc_last = ei; bc_last < num_edges && bc_last+1 == disk.getEdge(bc_last).next; bc_last++);
          const auto& e = disk.getEdge(ei);
          auto& v = disk.getVertex(e.vertex);
          const double t = double(ei-bi_start) / b_length;
          v.position.x = std::cos(2.0*M_PI * t);
          v.position.y = std::sin(2.0*M_PI * t);
     }
     // We only use the egde ID, which is not affected by the cutting
     // procedure, only vertex IDs are.
     auto edge_walker {make_edge_walker(mesh)};
     auto mv_coord_builder = [num_faces, &edge_walker, &mesh] (auto ei, auto vi, auto wi) {
          if (ei >= 3*num_faces)
               throw std::runtime_error("coord_builder: unexpected boundary edge #"s + to_string(ei));
          //std::cout << "coord_builder(" << ei << ", " << vi << ", " << wi << ")\n";
          const auto edge {edge_walker().e(ei)};
          const auto v {mesh.getVertex(edge.u())}, w {mesh.getVertex(edge.v())},
               vl {mesh.getVertex(edge().next().v())}, vr {mesh.getVertex(edge().flip().next().v())};
          const auto vec_vw = w.position-v.position, vec_vl = vl.position-v.position, vec_vr = vr.position-v.position;
          const double len_vw = glm::length(vec_vw);
          const double alpha_vw = std::acos(glm::dot(vec_vw, vec_vr) / (len_vw * glm::length(vec_vr))),
               beta_wv  = std::acos(glm::dot(vec_vw, vec_vl) / (len_vw * glm::length(vec_vl)));
          return (std::tan(alpha_vw/2) + std::tan(beta_wv/2))/ glm::length(vec_vw);
     };
     test_embedding::compute_disk_embedding(disk, mv_coord_builder);
     //auto mv_builder = make_mv_coord_generator(mesh);
     //compute_disk_embedding(disk, mv_builder);
     format::off::write(disk, "mt-new-2-disk.off");
     write_off(disk, "mt-new-2-disk.3.off");

     auto uniform_coord_builder = [] (auto ei, auto vi, auto wi) { return 1.0; };
     test_embedding::compute_disk_embedding(disk, uniform_coord_builder);
     format::off::write(disk, "mt-new-2-unidisk.off");
     write_off(disk, "mt-new-2-unidisk.3.off");

     for (hpindex ei = bi_start, count=0; ei < num_edges; ei++, count++) {
          const auto& e = disk.getEdge(ei);
          auto& v = disk.getVertex(e.vertex);
          const double t = double(ei-bi_start) / b_length;
          v.position.x = std::cos(2.0*M_PI * t);
          v.position.y = std::sin(2.0*M_PI * t);
     }

     test_embedding::compute_disk_embedding(disk, uniform_coord_builder);
     format::off::write(disk, "mt-new-2-unitri.off");
     write_off(disk, "mt-new-2-unitri.3.off");
}    // }}}1

void test_mesh_embedding(const std::string& filename, const std::string& output_dir = ".", const std::string& cut_method = "curv") { // {{{1
     using namespace std::string_literals;
     auto laptime = std::chrono::high_resolution_clock::now();
     auto t_rst = [&laptime] () { return laptime = std::chrono::high_resolution_clock::now(); };
     auto t_log = [&laptime] (auto name) { show_time(name, roundtime(laptime)); };
     auto output_path {p(output_dir)};

     std::cout << "---- test_mesh_embedding(filename=\"" + filename + "\") ----\n";
// flip switch to test with "double-torus.off"
#if 0
     const auto double_nut {NutChain(2, 1.5, 1.0, 1.0, 0.5)};
     const auto nut_mesh {make_triangle_graph(make_triangle_mesh<VertexP3>(double_nut, VertexFactory<VertexP3>()))};
#else
     const auto raw_mesh {format::off::read(filename)};
     const auto nut_mesh_simple {make_triangle_mesh<VertexP3>(raw_mesh)};
     const auto nut_mesh {make_triangle_graph(nut_mesh_simple)};
#endif
     hpuint nv = nut_mesh.getNumberOfVertices();
     hpuint ne = nut_mesh.getNumberOfEdges();
     hpuint nf = nut_mesh.getNumberOfTriangles();
     std::cout << nv << " vertices\n";
     std::cout << ne << " edges\n";
     std::cout << nf << " faces\n";
     std::cout << "genus " << ((2.0 - nv + ne/2.0 - nf)/2.0) << "\n";
     t_rst();
     Indices gen_cut_edges;
     if (cut_method == "curv"s)
          gen_cut_edges = curvdist_cut(nut_mesh);
     else
          gen_cut_edges = hopdist_cut(nut_mesh.getEdges());
     const auto the_cut {trim(nut_mesh, gen_cut_edges)};
     t_log(cut_method+"dist_cut"s);
     format::hph::write(the_cut, output_path / p("the-cut-raw.hph"));
     //const auto the_cut {trim(nut_mesh, cut(nut_mesh))};
     auto cut_graph {cut_graph_from_edges(nut_mesh, the_cut)};
     remove_chords(cut_graph, nut_mesh);
     format::hph::write(cut_edges(cut_graph), output_path / p("the-cut-nochords.hph"));
     //auto nut_disk_result {cut_to_disk(cut_graph, nut_mesh, CutGraph::VertexMapping::CONTIGUOUS_BOUNDARY)};
//
// Caution: KEEP_INNER is required in order to get a meaningful representation
// of the original mesh with disk topology (i.e. restoring the original mesh
// in 3D, but with boundary along the cut). This representatino is used for
// checkerboard rendering, along with the "uv coordinates" written as .xyz
     auto nut_disk_result {cut_to_disk(cut_graph, nut_mesh, CutGraph::VertexMapping::KEEP_INNER)};
     auto& nut_disk = std::get<0>(nut_disk_result);
     //auto& boundary_info = std::get<1>(nut_disk_result);
     std::cout << size(cut_edges(cut_graph)) << " cut edges\n";

     using namespace std::string_literals;
     using std::to_string;
     const auto num_faces = nut_disk.getNumberOfTriangles();
     const auto num_edges = nut_disk.getEdges().size();
     // Map boundary vertices onto unit circle
     const auto bi_start = 3*num_faces;
     const auto b_length = num_edges - bi_start;
     std::cout << b_length << " boundary length (should be equal to number of cut edges)\n";
     for (hpindex ei = bi_start; ei < num_edges; ei++) {
          //hpindex bc_last;       // last edge index of current boundary component
          //for (bc_last = ei; bc_last < num_edges && bc_last+1 == nut_disk.getEdge(bc_last).next; bc_last++);
          const auto& e = nut_disk.getEdge(ei);
          auto& v = nut_disk.getVertex(e.vertex);
          const double t = double(ei-bi_start) / b_length;
          v.position.x = std::cos(2.0*M_PI * t);
          v.position.y = std::sin(2.0*M_PI * t);
     }
     format::off::write(nut_disk, output_path / p("the-disk-raw.off"));
     // We only use the egde ID, which is not affected by the cutting
     // procedure, only vertex IDs are.
     auto edge_walker {make_edge_walker(nut_mesh)};
     auto mv_coord_builder = [num_faces, &edge_walker, &nut_mesh] (auto ei, auto vi, auto wi) {
          if (ei >= 3*num_faces)
               throw std::runtime_error("coord_builder: unexpected boundary edge #"s + to_string(ei));
          //std::cout << "coord_builder(" << ei << ", " << vi << ", " << wi << ")\n";
          const auto edge {edge_walker().e(ei)};
          const auto v {nut_mesh.getVertex(edge.u())}, w {nut_mesh.getVertex(edge.v())},
               vl {nut_mesh.getVertex(edge().next().v())}, vr {nut_mesh.getVertex(edge().flip().next().v())};
          const auto vec_vw = w.position-v.position, vec_vl = vl.position-v.position, vec_vr = vr.position-v.position;
          const double len_vw = glm::length(vec_vw);
          const double alpha_vw = std::acos(glm::dot(vec_vw, vec_vr) / (len_vw * glm::length(vec_vr))),
               beta_wv  = std::acos(glm::dot(vec_vw, vec_vl) / (len_vw * glm::length(vec_vl)));
          return (std::tan(alpha_vw/2) + std::tan(beta_wv/2))/ glm::length(vec_vw);
     };
     t_rst();
     test_embedding::compute_disk_embedding(nut_disk, mv_coord_builder);
     t_log("compute_disk_embedding");
     //auto mv_builder = make_mv_coord_generator(nut_mesh);
     //compute_disk_embedding(nut_disk, mv_builder);
     format::off::write(nut_disk, output_path / p("the-disk.off"));
     write_off(nut_disk, output_path / p("the-disk.3.off"));
     {
          std::vector<VertexP3> disktopo_verts;
          disktopo_verts.reserve(nut_disk.getVertices().size());
          std::vector<hpvec2> uv_coords;
          uv_coords.reserve(nut_disk.getVertices().size());
          for (const auto& v : nut_disk.getVertices()) {
               disktopo_verts.emplace_back(nut_mesh.getVertex(v.org_id).position);
               uv_coords.emplace_back(v.position.x, v.position.y);
          }
          auto cut_mesh {make_triangle_mesh<VertexP3>(disktopo_verts, make_indices(nut_disk))};
          format::off::write(cut_mesh, output_path / p("the-cutmesh.off"));
          format::xyz::write(uv_coords, output_path / p("the-uv-coords.xyz"));
     }

#if 0
     std::vector<Point2D> boundary_verts;
     boundary_verts.reserve(the_cut.size());
     boundary_verts.emplace(1.0, 0.0);
     std::unordered_set<hpindex> cut_verts(the_cut.size());
     cut_verts.emplace(edge_walker.e(the_cut[0]).u());
     for (unsigned ei = 0, cut_length = the_cut.size(); ei < cut_length-1; ei++) {
          double t = double(ei+1) / double(cut_length);
          boundary_verts.emplace(std::cos(2.0*M_PI * t), std::sin(2.0*M_PI * t));
          cut_verts.emplace(edge_walker.e(the_cut[ei]).v());
     }
     t_rst();
     auto embed_verts {embed(nut_mesh, the_cut, boundary_verts)};
     t_log("embed");
     auto& nut_disk_verts = nut_disk.getVertices();
     for (hpindex vi = 0, num_verts = nut_mesh.getNumberOfVertices(); vi < num_verts; vi++) {
          if (cut_verts.count(vi) == 0)
               nut_disk_verts[vi].position = embed_verts[vi];
     }
     format::off::write(nut_disk, output_path / p("the-disk-embed.off"));
     write_off(nut_disk, output_path / p("the-disk-embed.3.off"));
#endif
}    // }}}1

int main(int argc, const char** argv) {
     using namespace std::string_literals;

     try {
          test_minitorus_embedding();
          std::string mesh_filename = "double-torus.off"s;
          if (argc >= 2)
               mesh_filename = std::string(argv[1]);
          std::string output_dir {"."};
          if (argc >= 3)
               output_dir = std::string(argv[2]);
          std::string cut_method = "hop"s;
          if (argc >= 4)
               cut_method = std::string(argv[3]);
          test_mesh_embedding(mesh_filename, output_dir, cut_method);
     } catch(const std::exception& err) {
          utils::_log_error(std::string("Caught exception: ")+std::string(err.what()));
          g_testfail++;
     }
     return (g_testfail == 0) ? 0 : 1;
}
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
