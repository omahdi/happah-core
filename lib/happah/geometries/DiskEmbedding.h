/// \file DiskEmbedding.h
/// \brief Methods for cutting surface meshes into disks, producing Tutte
/// embeddings and computing projective structures from them.
///
/// \author Obada Mahdi <omahdi@gmail.com>
/// \copyright Copyright 2017 Obada Mahdi <omahdi@gmail.com>
/// Distributed under the Boost Software License, Version 1.0.
/// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#ifndef DISKEMBEDDING_H
#define DISKEMBEDDING_H

#include <algorithm>
#include <memory>
#include <stack>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <happah/geometries/TriangleMesh.h>

/// Default \c LOG_DEBUG macro, disables debug messages
#ifndef LOG_DEBUG
#define DISKEMBEDDING_H__LOG_DEBUG
#define LOG_DEBUG(...)   { void(0); }
#endif

/// namespace happah
namespace happah {
/// namespace obi
namespace obi {
// {{{ -- namespace detail
namespace detail {
void vertex_transfer_color(...) { }
template<class TV, class SV>
auto
vertex_transfer_color(TV& _t, const SV& _s) -> std::enable_if_t<sizeof(_t.color) && sizeof(_s.color)> {
     _t.color = _s.color;
}
}    // namespace detail
// }}} -- namespace detail

/// An iterator-like class for navigating along edges of a (triangle) mesh.
// {{{ -- EdgeWalker: navigating along edges
template<class _Mesh>
class EdgeWalker;

/// Specialization of EdgeWalker for directed-edge triangle meshes
template<class _Vertex>
class EdgeWalker<TriangleMesh<_Vertex, Format::DIRECTED_EDGE>> {
     static constexpr hpindex NIL_INDEX {std::numeric_limits<hpindex>::max()};
public:
     using Vertex = _Vertex;
     static constexpr Format MeshFormat = Format::DIRECTED_EDGE;
     using Mesh = TriangleMesh<Vertex, MeshFormat>;

     using value_type = Edge;
     using reference_type = const Edge&;

private:
     const Mesh* m_mesh {nullptr};
     hpindex m_ei {NIL_INDEX};
     const Edge* m_e {nullptr};

// Internal helper: check if initialized with mesh
     inline void check_mesh() const {
          if (m_mesh == nullptr)
               throw std::runtime_error("uninitialized EdgeWalker (not associated with a mesh)");
     }
     inline void check_edge() const {
          if (m_ei == NIL_INDEX)
               throw std::runtime_error("uninitialized EdgeWalker (not pointing to any edge)");
     }
     inline void check_init() const {
          check_mesh();
          check_edge();
     }
     inline void _set_ei(hpindex _ei) noexcept {
          m_e = &m_mesh->getEdge(_ei);
          m_ei = _ei;
     }
     inline EdgeWalker& _next() { _set_ei(m_e->next); return *this; }
     inline EdgeWalker& _prev() { _set_ei(m_e->previous); return *this; }
     inline EdgeWalker& _flip() { _set_ei(m_e->opposite); return *this; }
     hpindex _u() const   { return m_mesh->getEdge(m_e->opposite).vertex; }
     hpindex _v() const   { return m_e->vertex; }
     hpindex _opp() const { return m_e->opposite; }

public:
     EdgeWalker() { }
     EdgeWalker(const EdgeWalker& _x) : m_mesh{_x.m_mesh}, m_ei{_x.m_ei}, m_e{_x.m_e} { }
     EdgeWalker(const Mesh& _mesh) : m_mesh{&_mesh} { }
     EdgeWalker(const Mesh& _mesh, hpindex _ei) : m_mesh{&_mesh} { e(_ei); }
     ~EdgeWalker() { }

/// Binary op: equal
     bool operator==(const EdgeWalker& _x) const {
          return (this->m_ei == _x.m_ei) && (this->m_mesh == _x.m_mesh);
     }
/// Binary op: not equal
     bool operator!=(const EdgeWalker& _x) const {
          return !(*this == _x);
     }

/// Return source vertex
     hpindex u() const { check_init(); return _u(); }
/// Return target vertex
     hpindex v() const { check_init(); return _v(); }
/// Return pair (u,v) of edge source and target vertices
     std::pair<hpindex, hpindex> uv() const { check_init(); return std::make_pair(_u(), _v()); }
/// Return edge index
     hpindex e() const { check_init(); return m_ei; }
/// Return opposite edge index
     hpindex opp() const { check_init(); return _opp(); }
/// Dereference to Edge
     reference_type operator*() const { check_init(); return *m_e; }
/// Return copy
     EdgeWalker operator()() const { return {*this}; }

/// Assign edge index \p _ei
     EdgeWalker& e(hpindex _ei) {
          check_mesh();
          assert(_ei < m_mesh->getEdges().size());
          //   throw std::invalid_argument("EdgeWalker: edge index out of bounds");
          _set_ei(_ei);
          return *this;
     }
/// Assign edge corresponding to vertex pair \p (u,v)
     EdgeWalker& uv(hpindex _u, hpindex _v) {
          check_mesh();
          assert(_u < m_mesh->getVertices().size());
          assert(_v < m_mesh->getVertices().size());
          auto ei = make_edge_index(*m_mesh, _u, _v);
          if (!ei)
               throw std::invalid_argument("EdgeWalker: no directed edge corresponds to given vertex pair");
          _set_ei(*ei);
          return *this;
     }
/// Move forward
     EdgeWalker& next() { check_init(); return _next(); }
/// Move backward
     EdgeWalker& prev() { check_init(); return _prev(); }
/// Move to opposite
     EdgeWalker& flip() { check_init(); return _flip(); }
/// Move forward; same as next()
     EdgeWalker& operator++() { return next(); }
/// Move backward; same as prev()
     EdgeWalker& operator--() { return prev(); }
/// Move forward (post-increment)
     EdgeWalker operator++(int) { auto result {*this}; next(); return result; }
/// Move backward (post-decrement)
     EdgeWalker operator--(int) { auto result {*this}; prev(); return result; }
/// Move forward by offset
     EdgeWalker& operator+=(int _offset) {
          if (_offset < 0)
               _offset = 3 - (-_offset % 3);
          switch (_offset % 3) {
               case 1: next(); break;
               case 2: prev(); break;
          }
          return *this;
     }
/// Move backward by offset
     EdgeWalker& operator-=(int _offset) { return (*this) += -_offset; }
/// Move backward by offset (returns copy)
     EdgeWalker operator+(int _offset) const { return EdgeWalker{*this} += _offset; }
/// Move backward by offset (returns copy)
     EdgeWalker operator-(int _offset) const { return EdgeWalker{*this} -= _offset; }
/// Move to opposite (returns copy)
     EdgeWalker operator-() const { return EdgeWalker{*this}.flip(); }
};

/// Convenience function for constructing an \c EdgeWalker
template<class _V, Format _fmt>
EdgeWalker<TriangleMesh<_V, _fmt>>
make_edge_walker(const TriangleMesh<_V, _fmt>& _mesh) { return {_mesh}; }

template<class _V, Format _fmt>
EdgeWalker<TriangleMesh<_V, _fmt>>
make_edge_walker(const TriangleMesh<_V, _fmt>& _mesh, hpindex _ei) { return {_mesh, _ei}; }
// }}} -- EdgeWalker: navigating along edges

/// 2D vertex for embeddings with disk topology.
///
/// In addition to a regular 2D vertex with absolute position, this type
/// stores additional information about side pairings and the relation to the
/// original (closed-topology) mesh.
/// FIXME: temporarily changing to 3D vertex; make_projective_structure
/// currently depends on it
/// FIXME: color vs. no color? Base vertex should probably be configurable
class DiskVertex : public VertexPC<Space3D> {    // {{{
public:
     using Point = Point2D;   ///< type for storing the vertex position

private:
     using BaseVertex = VertexPC<Space3D>;
     using BasePoint = typename std::decay_t<decltype(std::declval<BaseVertex>().position)>;
     static constexpr hpindex IS_INNER_VALUE {std::numeric_limits<hpindex>::max()};

     static inline BasePoint _build_point(Point2D _p) { return {_p.x, _p.y, 1.0}; }
     static inline BasePoint _build_point(Point3D _p) { return _p; }

public:
/// Vertex index in the original mesh.
     hpindex org_id {0};
/// Vertex index of the following boundary vertex. Only meaningful if
/// is_inner() is false, undefined otherwise (an implementation-defined magic
/// value is used as indicator).
     hpindex next_v {IS_INNER_VALUE};
/// Default-constructs a vertex at the origin, marked as inner vertex.
     DiskVertex() : BaseVertex(_build_point(Point2D{0.0, 0.0})) {}
/// Constructs a vertex at \p _position, marked as inner vertex.
     DiskVertex(Point _position) : BaseVertex(_build_point(_position)) {}
/// Constructs a vertex at \p _position, setting its original index to \p
/// _org_id and optionally marking it as vertex on the boundary, with vertex
/// \p _next_v following on the boundary.
     DiskVertex(Point _position, hpindex _org_id, hpindex _next_v = IS_INNER_VALUE)
          : BaseVertex(_build_point(_position)), org_id {_org_id}, next_v {_next_v} { }
/// Returns whether this vertex lies on the border.
     bool is_inner() const noexcept { return next_v == IS_INNER_VALUE; }
/// Mark vertex as inner (non-boundary) vertex.
     void set_inner() noexcept { next_v = IS_INNER_VALUE; }

     template<class _Mesh>
     friend bool is_inner_vertex(const _Mesh&, const DiskVertex& _v) {
          return _v.is_inner();
     }
};   // }}} class DiskVertex

using DiskMesh = TriangleMesh<DiskVertex, Format::DIRECTED_EDGE>;

/// Additional information for boundary edges required by
/// make_projective_structure.
struct BoundaryEdgeInfo {     // {{{
/// Indicates which side of the fundamental polygon an edge lies on.
     hpindex side {0};
/// Offset relative to the start of segment \p side in cut circuit
     hpindex offset {0};
/// Index of the paired edge (opposite half-edge).
     hpindex opposite {0};
     BoundaryEdgeInfo() { }
     BoundaryEdgeInfo(hpindex _s, hpindex _o, hpindex _opp) : side {_s}, offset {_o}, opposite {_opp} { }
     BoundaryEdgeInfo(const BoundaryEdgeInfo&) = default;
     BoundaryEdgeInfo(BoundaryEdgeInfo&&) = default;
     ~BoundaryEdgeInfo() { }
};   // }}} struct BoundaryEdgeInfo

/// Representation of a cut graph as a topologically sorted list of half-edges
///
/// Edges are sorted according to the order of a walk along the Eulerian
/// circuit that will become the boundary of the resulting disk after cutting.
/// Edges are stored as a compact "array of arrays", similar to
/// \a happah::Arrays<hpindex> (but storing start indices instead of lengths,
/// and adding an extra entry at the end to allow easy computation of lengths
/// without distinction for the last item).
class CutGraph {    // {{{
public:
     using BoundaryInfo = std::unordered_map<hpindex, BoundaryEdgeInfo>;

/// Constants for controlling how vertex indices are generated or remapped:
/// - \p CONTIGUOUS_BOUNDARY reorders indices such that all inner vertices are
///   assigned indices 0 to n-1, followed by boundary vertices numbered from
///   n to n+b-1, where n+b is the total number of vertices and b the number
///   of vertices on the disk boundary, respectively
/// - \p KEEP_INNER leaves indices of all inner vertices unchanged and
///   assigns new indices for duplicated vertices on the boundary as needed
     enum class VertexMapping {
          CONTIGUOUS_BOUNDARY,
          KEEP_INNER
     };

private:
// FIXME make circuit_type and segments_type public?
     using circuit_type = std::vector<hpindex>;
     using segments_type = std::vector<hpindex>;
/// List of (half-)edge indices, topologically sorted into an Eulerian circuit
/// of the cut graph. When built using cut_graph_from_edges() or
/// cut_graph_from_paths(), the first edge is guaranteed to originate at a
/// branch node.
     circuit_type m_circuit;
/// Offsets into m_circuit of first edge of each segment; an additional entry
/// at index segment_count() is included to simplify computations of segment
/// lenghts or index ranges without distinction of cases for the last segment.
     segments_type m_segments;

// Check if cut graph is usable (non-empty)
     void check() const {
          if (m_segments.size() == 0 || m_circuit.size() == 0)
               throw std::runtime_error("CutGraph: not holding with any cut segments, possibly uninitialized");
     }

/// Simple static range begin()/end() support representing a range of indices
/// corresponding to a cut segment ("cut path").
///
/// This uses the \c const_iterator of the underlying data type
/// std::vector<hpindex>.
     class segment_range_const {
     public:
          using iterator = circuit_type::const_iterator;
          using value_type = circuit_type::value_type;

          segment_range_const(iterator _begin,  hpindex _si, hpindex _ei) : m_begin{_begin+_si}, m_end{_begin+_ei} { }
          segment_range_const(const segment_range_const&) = default;
          segment_range_const(segment_range_const&&) = default;
          ~segment_range_const() = default;

          iterator begin() const { return m_begin; }
          iterator end() const { return m_end; }
          iterator cbegin() const { return begin(); }
          iterator cend() const { return end(); }

     private:
          const iterator m_begin, m_end;
     };

/// Simple static range begin()/end() support representing a range of indices
/// corresponding to a cut segment ("cut path").
///
/// This uses the \c iterator of the underlying data type
/// std::vector<hpindex>.
     class segment_range {
     public:
          using iterator = circuit_type::iterator;
          using value_type = circuit_type::value_type;

          segment_range(iterator _begin, hpindex _si, hpindex _ei) : m_begin{_begin+_si}, m_end{_begin+_ei} { }
          segment_range(const segment_range&) = default;
          segment_range(segment_range&&) = default;
          ~segment_range() = default;

          iterator begin() const { return m_begin; }
          iterator end() const { return m_end; }

     private:
          const iterator m_begin, m_end;
     };

public:
     CutGraph() { }
     CutGraph(const CutGraph& _g) : m_circuit{_g.m_circuit}, m_segments{_g.m_segments} { }
     CutGraph(CutGraph&& _g) : m_circuit{std::move(_g.m_circuit)}, m_segments{std::move(_g.m_segments)} { }
     ~CutGraph() { }

/// Returns a const reference to the underlying vector of edge indices.
///
/// \sa segment_index(), segment_length()
     friend const auto& cut_edges(const CutGraph& cut_graph) { return cut_graph.m_circuit; }
/// Returns a reference to the underlying vector of edge indices.
///
/// \sa segment_index(), segment_length()
     friend auto& cut_edges(CutGraph& cut_graph) { return cut_graph.m_circuit; }
/// Returns the number of cut paths (connected chains of degree 2).
     friend hpindex segment_count(const CutGraph& cut_graph) noexcept {
          cut_graph.check();
          return cut_graph.m_segments.size()-1;
     }
/// Returns an index into the vector returned by cut_edges() of the start of
/// the <tt>_i</tt>-th cut segment in \p cut_graph.
///
/// \sa segment_length(), cut_edges()
     friend auto segment_index(const CutGraph& cut_graph, hpindex _i) {
          cut_graph.check();
          assert(_i < segment_count(cut_graph));
          return cut_graph.m_segments[_i];
     }
/// Returns the number of edges of the <tt>_i</tt>-th cut segment in
/// \p cut_graph.
///
/// \sa segment_length(), cut_edges()
     friend std::size_t segment_length(const CutGraph& cut_graph, hpindex _i) noexcept {
          cut_graph.check();
          assert(_i < segment_count(cut_graph));
          return cut_graph.m_segments[_i+1] - cut_graph.m_segments[_i];
     }
/// Returns a range object suitable for iterating over the const-qualified
/// edge indices of the <tt>_i</tt>-th cut segment of \p cut_graph
/// (range-for-loop compatible).
     friend segment_range_const cut_segment(const CutGraph& cut_graph, hpindex _i) {
          using std::cbegin;
          using std::cend;
          const hpindex start = segment_index(cut_graph, _i),  length = segment_length(cut_graph, _i);
          return {std::cbegin(cut_edges(cut_graph)), start, start+length};
     }
/// Returns a range object suitable for iterating over the edge indices of the
/// <tt>_i</tt>-th cut segment of \p cut_graph (range-for-loop compatible).
     friend segment_range cut_segment(CutGraph& cut_graph, hpindex _i) {
          using std::begin;
          using std::end;
          const hpindex start = segment_index(cut_graph, _i), length = segment_length(cut_graph, _i);
          return {std::begin(cut_edges(cut_graph)), start, start+length};
     }

// see definition for documentation
     template<class Mesh>
     friend CutGraph cut_graph_from_edges(const Mesh& source_mesh, const std::vector<hpindex>& cut);
// see definition for documentation
     template<class Mesh>
     friend CutGraph cut_graph_from_paths(const Mesh& source_mesh, const Arrays<hpindex>& paths);
// see definition for documentation
     template<class Mesh>
     friend bool remove_chords(CutGraph& cut_graph, const Mesh& source_mesh);

     template<class S>
     friend S& operator<<(S& _s, const CutGraph& _v) {
          return _s << _v.m_circuit << _v.m_segments;
     }
     template<class S>
     friend S& operator>>(S& _s, CutGraph& _v) {
          return _s >> _v.m_circuit >> _v.m_segments;
     }

/// Convert cut graph into Arrays<hpindex>
     operator Arrays<hpindex>() const {
          check();
          using std::cbegin;
          using std::cend;
          Arrays<hpindex> result;
          result.reserve(segment_count(*this), m_circuit.size());
          for (unsigned i = 0; i < segment_count(*this); i++) {
               auto range = cut_segment(*this, i);
               result.push_back(cbegin(range), cend(range));
          }
          return result;
     }
};   // }}} class CutGraph

/// Build a CutGraph from an (unordered) list of edge indices.
template<class SourceMesh>
CutGraph
cut_graph_from_edges(const SourceMesh& source_mesh, const std::vector<hpindex>& cut) {    // {{{
     using std::begin;
     using std::end;
     if (cut.size() == 0)
          throw std::invalid_argument("An empty graph does not have a circuit.");
     CutGraph cut_graph;
     auto& circuit {cut_graph.m_circuit};
     auto& segments {cut_graph.m_segments};

// Boolean flags for quick O(1) testing if edge is in cut. Linear-time setup,
// linear space (1 bit per edge).
//
// Note: The input specifies only *one* directed edge; the cut has to consider
// both directions.
// TODO: properly handle boundary edges, which will have only *one* direction
// in the cut!
     auto edge_walker = make_edge_walker(source_mesh);
     std::vector<bool> is_cut_edge(3*source_mesh.getNumberOfTriangles());
     for (auto ei : cut) {
          is_cut_edge[ei] = true;
          is_cut_edge[edge_walker.e(ei).flip().e()] = true;
     }
// Walk edges in a left-first fashion, effectively rearranging `cut` in the
// order we visit them.
     std::unordered_map<hpindex, unsigned int> vertex_count(cut.size()/2);
     //std::unordered_set<hpindex> branch_nodes;
     std::vector<hpindex> branch_nodes;
     auto start_edge = cut[0];
     auto edge = edge_walker().e(start_edge);
     auto current_v = edge.u();         // track source vertex of current edge
     circuit.clear();
     circuit.reserve(cut.size()*2);
     do {
// Note: default-constructed vertex_count value is initialized to 0:
          auto branch_node_it {std::lower_bound(begin(branch_nodes), end(branch_nodes), current_v)};
          LOG_DEBUG(6, "  branch_node_it at position %d/%d", std::distance(begin(branch_nodes), branch_node_it), branch_nodes.size());
          if (++vertex_count[current_v] > 2 && (branch_node_it == end(branch_nodes) || *branch_node_it != current_v))
               branch_nodes.insert(branch_node_it, current_v);
          circuit.emplace_back(edge.e());
// Visit spokes of next vertex, starting with edge next to \p edge, which is
// the clockwise first outgoing edge of current_v following \p edge.
          auto last_edge = edge;
          LOG_DEBUG(6, "  edge %d = (%d, %d), opposite = %d", edge.e(), edge.u(), edge.v(), edge.opp());
          do {
               edge.next();
               if (is_cut_edge[edge.e()]) {
                    current_v = last_edge.v();
                    break;
               }
          } while (edge.flip() != last_edge);
          if (edge == last_edge)
               throw std::invalid_argument("unable to find closed circuit");
               //utils::log_error("Error: Unable to find successor for edge %d at vertex #%d", ei, current_v);
     } while (edge.e() != start_edge);

     const auto num_segments = std::accumulate(branch_nodes.cbegin(), branch_nodes.cend(), 0u,
          [&vertex_count](auto s, auto v) { return s+vertex_count[v]; });
     segments.clear();
     segments.reserve(num_segments+1);
     LOG_DEBUG(3, "- found %d cut nodes: %s", branch_nodes.size(), utils::str(branch_nodes));
     LOG_DEBUG(3, "- computing positions of cut nodes in %u-gon", num_segments);
     unsigned int first_branch_offset = 0;
     for (hpindex i = 0; i < circuit.size(); i++) {
          auto v = edge_walker.e(circuit[i]).u();
          if (std::count(begin(branch_nodes), end(branch_nodes), v) > 0) {
               if (segments.size() == 0)
                    first_branch_offset = i;
               if (!(segments.size() < num_segments))
                    throw std::out_of_range("assertion segments.size() < num_segments failed");
               segments.emplace_back(i-first_branch_offset);
          }
     }
     segments.emplace_back(circuit.size());
// Rotate m_circuit such that first element is an edge emanating from the
// first branch node seen.
     std::rotate(circuit.begin(), circuit.begin()+first_branch_offset, circuit.end());
     return cut_graph;
}    // }}} cut_graph_from_edges()

/// Turns a list of path segments into a list of edges.
template<class SourceMesh>
auto
edges_from_paths(const SourceMesh& mesh, const Arrays<hpindex>& paths) {   // {{{
     std::vector<hpindex> cut;

     for (auto path : paths) {
          if (std::distance(path.first, path.second) < 2)
               throw std::invalid_argument("invalid cut path: expected sequence of at least two vertices");
          const auto last_node = path.second-1;
          for (auto it = path.first; it != last_node; it++) {
               auto e = make_edge_index(mesh, it[0], it[1]);
               if (!e)
                    throw std::invalid_argument("invalid cut path: no edge connecting consecutive pair of vertices");
               cut.emplace_back(*e);
          }
     }
     return cut;
}    // }}}

/// Build a CutGraph from a list of paths, i.e. sequences of vertex indices.
template<class SourceMesh>
CutGraph
cut_graph_from_paths(const SourceMesh& mesh, const Arrays<hpindex>& paths) {    // {{{
     return cut_graph_from_edges(mesh, edges_from_paths(mesh, paths));
}    // }}} cut_graph_from_paths

/// Builds an array of indices mapping the index of a cut path segment to the
/// index of its paired counterpart.
///
/// \note This is an example where it would be much easier if we stored vertex
/// indices instead of edge indices, because we would not have to deal with a
/// mesh and a (possibly expensive, in the Format::SIMPLE case)
/// make_edge_walker() and detect opposing edges by looking at the 1-ring.
/// \note Complexity:
/// - space: linear in the number of segments
/// - time: quadratic time in the number of segments (FIXME: more careful
/// analysis?
template<class Mesh>
auto segment_pairings(const CutGraph& cut_graph, const Mesh& mesh) {  // {{{
     const auto num_segments = segment_count(cut_graph);
     const auto& circuit {cut_edges(cut_graph)};
     const auto circuit_length {circuit.size()};
     const hpindex NIL_VALUE = std::numeric_limits<hpindex>::max();
     std::vector<hpindex> pairing(num_segments, NIL_VALUE);
     auto edge_walker = make_edge_walker(mesh);
     hpindex num_pairs = 0;
     for (hpindex side = 0; side < num_segments; side++) {
          if (pairing[side] != NIL_VALUE)
               continue;
// w: position in m_circuit of first outgoing edge of segment following _side
          auto w = segment_index(cut_graph, (side+1) % num_segments);
// w_prev: take one step back from w and find opposite edge
          auto w_prev = w == 0 ? circuit_length-1 : w-1;
          auto out_ei = edge_walker.e(circuit[w_prev]).opp();
          for (hpindex i = 0; i < num_segments; i++) {
               if (circuit[segment_index(cut_graph, i)] != out_ei)
                    continue;
               pairing[side] = i;
               pairing[i] = side;
               num_pairs++;
               break;
          }
     }
     assert(2*num_pairs == num_segments);
     return pairing;
};   // }}} segment_pairings()

/// Build array of neighbors of triangle fan representing our fundamental
/// region circumscribed by \p cut_graph.
///
/// \param[in] cut_graph
/// \param[in] source_mesh
///
/// \returns An array of length 3*cut_graph.segment_count() holding the indices of
/// faces adjacent to face j at positions 3*j..3*j+2. We follow the convention
/// that adjacent triangles are given in counter-clockwise order. Assuming
/// that each face is defined as the vertex triple (v0, v1, v2) and v0 is the
/// common vertex of the macro fan, neighbors of triangle i in a fan of size n
/// are given in the order ((i-1) mod n, opposite(i), (i+1) mod n)).
template<class SourceMesh>
std::vector<hpindex>
make_neighbors(const CutGraph& cut_graph, const SourceMesh& mesh) { // {{{
     const auto num_segments = segment_count(cut_graph);
     std::vector<hpindex> neighbors(3*num_segments);
     auto pairings {segment_pairings(cut_graph, mesh)};

     for (auto i = 0u; i < num_segments; i++) {
          hpindex paired_pos = pairings[i];
          neighbors[3*i + 0] = (i+num_segments-1) % num_segments;
          neighbors[3*i + 1] = (paired_pos+num_segments) % num_segments;
          neighbors[3*i + 2] = (i+1) % num_segments;
     }
     return neighbors;
}    // }}} make_neighbors()

/// Replaces each cut segment of a given cut graph by with a shortest path
/// along the its linear graph, augmented by additional edges from the
/// underlying mesh, keeping its homotopy class intact.
///
/// The goal of this procedure is to ensure that all maximal linear subgraphs
/// of \p cut_graph (<em>cut segments</em>) are free of \e chords, by which we
/// mean edges from the original mesh that connect two non-consecutive nodes
/// of such subgraphs.
///
/// \note TODO: Complexity?
template<class Mesh>
bool
remove_chords(CutGraph& cut_graph, const Mesh& mesh) { // {{{
     using std::begin;
     using std::end;
     const auto num_segments = segment_count(cut_graph);
     auto& circuit {cut_edges(cut_graph)};
     std::vector<bool> processed(num_segments);
     LOG_DEBUG(5, "remove_chords(): processing %d sides", num_segments);
     auto edge_walker {make_edge_walker(mesh)};
     const auto pairings {segment_pairings(cut_graph, mesh)};

     struct path_node {
          std::size_t distance {0};
          hpindex visited : 1, prev : (8*sizeof(hpindex)-1);
          hpindex rev_edge {0};
          path_node() : visited {0}, prev {0} { }
     };
// FIXME: no sanity check
     const auto sorted_insert = [](auto& _set, auto v) {
          using std::begin;
          using std::end;
          auto it {std::lower_bound(begin(_set), end(_set), v)};
          if (it == end(_set) || *it != v)
               _set.emplace(it, v);
     };
// get_range_containing(): Let {(v, w(i)): i=1,..,|N(v)|} be the
// (counter-clockwise) set of ordered outgoing edges of node v, where N(v)
// denotes the 1-ring of v and the order is determined by the oriented
// triangulation. Let b(1), ..., b(r) be indices in 1,...,|N(v)| such that
// (v, w(b(k))) is a cut edge, for k in the range 1,...,r.
//
// Given a branch edge ei, let k be the index with ei = (v, w(b(k)). Since
// branch nodes have minimum outgoing degree 3, we can assume w.l.o.g. that
// k=2 by shifting the the indices appropriately.
// The following method then returns the set
//   {(v, w(j)): b(1) < j < b(3)}.
// This set is non-empty because b(1) < b(2) < b(3) (as stated above, there
// are at least three outgoing edges that lie on the cut).
     const auto get_range_containing = [&](auto v, auto ei) {
// 1. Collect both incoming and outgoing edges of branch nodes
          std::vector<hpindex> branch_edges;
          branch_edges.reserve(num_segments);     // simple heuristic
          auto walker {edge_walker()};
          for (hpindex i = 0; i < num_segments; i++) {
               const auto range {cut_segment(cut_graph, i)};
               walker.e(*begin(cut_segment(cut_graph, i)));
               if (walker.u() != v)
                    continue;
               sorted_insert(branch_edges, walker.e());
               sorted_insert(branch_edges, walker.opp());
               const auto prev_range {cut_segment(cut_graph, i == 0 ? num_segments-1 : i-1)};
               sorted_insert(branch_edges, walker.e(*(end(prev_range)-1)).e());
               sorted_insert(branch_edges, walker.opp());
          }
          branch_edges.shrink_to_fit();
          const auto edges_begin = begin(branch_edges), edges_end = end(branch_edges);
// 2. Visit spokes in both counter-clockwise and clockwise order, starting with
//    ei, until hitting another edge in branch_edge
          std::vector<hpindex> edge_range;
          edge_range.reserve(8);        // FIXME reasonable heuristic?
          edge_range.emplace_back(ei);
          walker.e(ei);
          while (std::find(edges_begin, edges_end, walker.prev().flip().e()) == edges_end)
               edge_range.emplace_back(walker.e());
          assert(walker.e() != ei);     // branch nodes have out-degree >= 3 in cut graph
          walker.e(ei);
          while (std::find(edges_begin, edges_end, walker.flip().next().e()) == edges_end)
               edge_range.emplace_back(walker.e());
          return edge_range;
     };
     bool did_remove = false;
     for (std::size_t s_index = 0; s_index < num_segments; s_index++) {
          if (processed[s_index]) {
               LOG_DEBUG(6, "- segment %d already processed, skipping", s_index);
               continue;
          }
// Run Dijkstra's algorithm restricted on a path segment.  Care has to be
// taken that we do not accidentally reverse the orientation of a loop: We
// compute "allowed edges" for leaving and entering the first and the last
// vertex of the current path segment, which are edges in the same range
// between edges used by other cuts.
          const auto other_index = pairings[s_index];
          auto segment {cut_segment(cut_graph, s_index)};
          const auto s_begin = segment.begin(), s_end = segment.end();
          const auto s_length = std::distance(segment.begin(), segment.end());
          LOG_DEBUG(5, "- removing chords in segments %d and %d of length %d", s_index, other_index, s_length);
          std::vector<path_node> path_nodes(s_length+1);
          std::unordered_map<hpindex, hpindex> nodeset(s_length+1);
// Note: If the path segment forms a loop, the following loop leaves nodeset
// with a record of the 2nd position at the end.  Changing this behavior will
// break the visit_spokes() handler.
// Note: "emplace()" does not assign if key already exists; solved in C++17
// with "insert_or_assign()".
          const auto start_vertex = edge_walker.e(*s_begin).u();
          nodeset.emplace(start_vertex, 0);
          for (auto j = 0u; j < s_length; j++)
               nodeset[edge_walker.e(s_begin[j]).v()] = j+1;
          const auto end_vertex = edge_walker.e(*(s_end-1)).v();
          LOG_DEBUG(5, "  source, target: #%4d -> #%4d", start_vertex, end_vertex);
// Partition edges at start and end node; edge_walker still on last edge
          const auto
               allowed_out = get_range_containing(start_vertex, *s_begin),
               allowed_in = get_range_containing(end_vertex, edge_walker.opp());
          LOG_DEBUG(5, "  allowed_out = %s", allowed_out);
          LOG_DEBUG(5, "  allowed_in  = %s", allowed_in);
// Set up start node
          hpindex next_node = 0;
          path_nodes[next_node].prev = next_node;      // redundant, but repeat for clarity
          while (next_node != s_length) {
               const auto current_node = next_node;
// relax edges that have their target vertex on the current path segment and
// point *forward*, i.e. reach a higher index.
               const auto relax_begin = s_begin[current_node];
               edge_walker.e(relax_begin);
               do {
                    if (nodeset.count(edge_walker.v()) == 0)
                         continue;
                    LOG_DEBUG(9, "  relaxing edge %d (to node #%d)", edge_walker.e(), edge_walker.v());
// refuse to relax edge from start in reverse direction
                    if (current_node == 0 && std::count(allowed_out.cbegin(), allowed_out.cend(), edge_walker.e()) == 0)
                         continue;
                    const auto pos = nodeset[edge_walker.v()];
                    LOG_DEBUG(9, "    testing possible shortcut to node %d (-> pos %d in segment; edge %5d, opposite %5d)", edge_walker.v(), pos, edge_walker.e(), edge_walker.opp());
                    if (edge_walker.v() == end_vertex && std::count(allowed_in.cbegin(), allowed_in.cend(), edge_walker.opp()) == 0)
                         continue;
                    const auto d = path_nodes[next_node].distance;
                    if (pos <= next_node || (path_nodes[pos].visited && (d+1 >= path_nodes[pos].distance)))
                         continue;
                    path_nodes[pos].visited = 1;
                    path_nodes[pos].distance = d+1;
                    path_nodes[pos].prev = next_node;
                    path_nodes[pos].rev_edge = edge_walker.opp();
                    next_node = pos;
               } while (edge_walker.flip().next().e() != relax_begin);
               if (next_node == current_node) {
                    throw std::runtime_error("remove_chords(): could not find next edge in cut segment");
                    //throw std::runtime_error(utils::logfmt("could not find
                    //edge following vertex #%d (segment index %d; original
                    //edge is #%d)", _get_source_at(s_begin+next_node), next_node, get_edge_id_at(segment_start+next_node)));
               }
          }
          processed[other_index] = true;
          const int delta = s_length - path_nodes[next_node].distance;
// Sanity check (shouldn't happen)
          //THROW_IF(0, std::runtime_error, delta < 0, "Internal error: path cannot grow longer");
          assert(delta >= 0 && "Internal error: path cannot grow longer");
          if (delta == 0)
               continue;      // no improvement
          LOG_DEBUG(3, "  segment %d shortened by %d edges", s_index, delta);
          //throw std::runtime_error("Untested case where polygon sides have chords");
          did_remove = true;
// FIXME untested!
// Trace back over shortest path, updating circuit.
          auto fp = s_begin + path_nodes[s_length].distance;
          auto rp = begin(circuit) + segment_index(cut_graph, other_index);
          for (; next_node != 0; next_node = path_nodes[next_node].prev) {
               *(--fp) = edge_walker.e(path_nodes[next_node].rev_edge).opp();
               *(rp++) = path_nodes[next_node].rev_edge;
          }
// m_circuit layout: [...i..j..i+1..//..o..k..o+1...]
//   [i,i+1): i-th segment before; [i..j): after taking shortcuts
//   [o,o+1): o-th segment before; [o..k): after taking shortcuts
// Ranges [j, i+1) and [k, o+1) are to be erased; erase range starting at k
// first, keeping indices for first intact.
          circuit.erase(rp, begin(circuit)+segment_index(cut_graph, other_index < num_segments-1 ? other_index+1 : circuit.size()));
          circuit.erase(s_end-delta, s_end);
// Adjust segment indices: every range up to and including "other_index" moved
// by "delta"; everything following "other_index" by twice that.
          for (unsigned j = s_index+1; j <= other_index; j++)
               cut_graph.m_segments[j] -= delta;
          for (unsigned j = other_index+1; j < cut_graph.m_segments.size(); j++)
               cut_graph.m_segments[j] -= 2*delta;
     }
     return did_remove;
}    // }}} remove_chords()

/// Compute transition matrices for macro tetrahedra corresponding to a
/// tesselation {3; p,6,6}, where p is the valence of the center vertex.
std::vector<hpreal>
compute_transitions_p_3(std::vector<hpindex> neighbors) {   // {{{
     const auto num_segments = neighbors.size() / 3;
     std::vector<hpreal> transitions(9*num_segments);
     const hpreal coeff_c = std::cos(2.0*M_PI / num_segments);     // c = \cos(2\pi/n)
     const hpvec3 mat_B     {2.0-2.0*coeff_c, -1.0, 2.0*coeff_c};
     const hpvec3 mat_B_inv {mat_B[2], mat_B[1], mat_B[0]};
     const hpvec3 mat_U     {1.0 / (2.0 - 2.0*coeff_c), -1.0, 1.0 / (2.0 - 2.0*coeff_c)};

     const auto copyvec = [](auto real_it, const hpvec3& vec) {
          *(real_it++) = vec.x; *(real_it++) = vec.y; *(real_it++) = vec.z;
     };
     for (unsigned int i = 0; i < num_segments; i++) {
          // did_rev: clockwise; did_fwd: counter-clockwise
          bool did_rev {false}, did_fwd {false}, did_flip {false};
          for (unsigned int k = 0; k <= 2; k++) {
               if (i == ((neighbors[3*i+k]+1) % num_segments) && ((i+1) % num_segments) == neighbors[3*i+((k+2) % 3)]) {
                    if (did_rev)
                         throw std::invalid_argument("clockwise rotation already set earlier");
                    copyvec(transitions.begin() + 9*i + 3*k, mat_B_inv);
                    did_rev = true;
               } else if (((i+1) % num_segments) == neighbors[3*i+k] && i == ((neighbors[3*i+((k+1) % 3)]+1) % num_segments)) {
                    if (did_fwd)
                         throw std::invalid_argument("counter-clockwise rotation already set earlier");
                    copyvec(transitions.begin() + 9*i + 3*k, mat_B);
                    did_fwd = true;
               } else {
                    if (did_flip)
                         throw std::invalid_argument("reflection already set earlier");
                    copyvec(transitions.begin() + 9*i + 3*k, mat_U);
                    did_flip = true;
               }
          }
     }
     return transitions;
}    // }}}

/// Build mesh with disk topology by cutting \p source_mesh along
/// \p cut_graph.
///
/// \param[in] cut_graph
/// \param[in] source_mesh
///
/// This involves remapping vertex indices depending on the requested
/// \p mapping_mode:
/// - \p KEEP_INNER will only duplicate and assign new indices to vertices on
///   the boundary as needed, leaving existing indices intact.
/// - \p CONTIGUOUS_BOUNDARY will reorder vertices such that all inner
///   vertices come first, followed by those on the boundary.
///
/// \note \b TODO
/// - add mesh_format argument to allow selection of TriangleMesh format?
/// \param[in] mesh_format instance of an integral constant, evaluated at
///   compile time to choose between Format::DIRECTED_EDGE or Format::SIMPLE
///     FormatValue&& mesh_format = {}
///     class FormatValue = std::integral_constant<Format, Format::DIRECTED_EDGE>>
/// - complexity analysis; this is the most expensive operation
template<class Mesh>
std::tuple<DiskMesh, typename CutGraph::BoundaryInfo>
cut_to_disk(   // {{{
     const CutGraph& cut_graph,
     const Mesh& source_mesh,
     const CutGraph::VertexMapping mapping_mode = CutGraph::VertexMapping::CONTIGUOUS_BOUNDARY
) {
     const bool contiguous_mode = mapping_mode == CutGraph::VertexMapping::CONTIGUOUS_BOUNDARY;
     const auto& circuit = cut_edges(cut_graph);
     const auto circuit_length = circuit.size();
     CutGraph::BoundaryInfo boundary_info(circuit_length);
     const auto& mesh_indices {source_mesh.getIndices()};
     const auto& mesh_vertices {source_mesh.getVertices()};
// Start with copy of source_mesh face indices; its size does not change, even
// though technically we will have one additional face (the generally
// non-triangular "outer" face), which is not represented explicitly.
     std::vector<hpindex> disk_indices(std::begin(mesh_indices), std::end(mesh_indices));
// We can update vertex indices in a single sweep over the cut edges and one
// additional sweep over vertices to set up vertex data. If mode is set to
// KEEP_INNER, we can directly update indices. For CONTIGUOUS_BOUNDARY, we
// need to know the number of inner vertices so we can set up the first new_id
// accordingly. Alternatively, we compute the number of extra_vertices,
// allowing us to pre-allocate disk_vertices to the right size, and start
// allocation new vertex indices starting from the end.
// The following formula for extra_vertices follows from combining Euler's
// formula for the source_mesh and that for the expected disk, taking into
// account that all faces are triangulated. Note that the term
//   (mesh_indices.size()/6) - mesh_vertices.size()
// is equal to 2g-2, where g is the genus of source_mesh. This also follows
// from Euler's formula, taking into account that the number of faces of the
// triangulated closed source_mesh is mesh_indices.size()/3 and the number of
// (undirected) edges is mesh_indices.size()/2.
     assert((circuit_length % 2) == 0);
     assert((mesh_indices.size() % 6) == 0);
     const int extra_vertices = 1 + (circuit_length/2) + (mesh_indices.size()/6) - mesh_vertices.size();
     assert(extra_vertices > 0);
     LOG_DEBUG(3, "  mesh stats: %d faces, %d vertices ==> genus %d",
          mesh_indices.size()/3, mesh_vertices.size(), (2 + mesh_indices.size()/6 - mesh_vertices.size())/2);
     LOG_DEBUG(3, "  %d source vertices, %d extra vertices ==> %d disk vertices total",
          mesh_vertices.size(), extra_vertices, mesh_vertices.size()+extra_vertices);
// Initialize each vertex with a magic value for org_id in order to be able to
// detect uninitialized inner vertices later (Note: marked is_inner by default)
     const auto num_disk_vertices = mesh_vertices.size() + extra_vertices;
     std::vector<DiskVertex> disk_vertices(num_disk_vertices, {{}, std::numeric_limits<hpindex>::max()});
// v_mapped keeps track of whether a particular source vertex index has
// already been seen while walking along the cut.
     std::unordered_set<hpindex> v_mapped(circuit_length/2);
     const hpindex first_v = mesh_vertices.size() + extra_vertices - 1;
     hpindex new_id = first_v;
     LOG_DEBUG(5, "  first new ID: %d", first_v);
// prev_v holds the previous allocated vertex ID for updating its next_v
// member; it is not used in the first iteration, because we do not know the
// new ID of the last vertex along the cut yet.
     hpindex prev_v = 0;
     unsigned current_side = 0;
     auto prev = make_edge_walker(source_mesh);
     for (hpindex bi = 0, current_segment_index = segment_index(cut_graph, current_side); bi < circuit_length; bi++) {
          if (bi == current_segment_index+segment_length(cut_graph, current_side)) {
               current_segment_index = bi;
               ++current_side;
          }
          auto current_ei = circuit[bi];
// Initialize edge walker "prev" with predecessor of current_ei in circuit
// and store the vertex ID it is pointing at in current_v.
          auto current_v  = prev.e(circuit[bi > 0 ? bi-1 : circuit_length-1]).v();
          hpindex new_v = 0;
          bool mapped = v_mapped.count(current_v) > 0;
          if (!mapped)
               v_mapped.emplace(current_v);
          LOG_DEBUG(7, "  current_ei <- %d <- circuit[%d], side=%d, current_v=%d (mapped=%s)", current_ei, bi, current_side, current_v, mapped);
// Assign new ID if current_v has been seen before, or if a contiguous mapping
// of boundary vertices is requested.
          if (mapped || contiguous_mode)
               new_v = new_id--;
          else {
               new_v = current_v;
          }
          LOG_DEBUG(7, "    new_v <- %d", new_v);
          disk_vertices[new_v].org_id = current_v;
          if (bi != 0)        // last vertex is updated after the loop
               disk_vertices[prev_v].next_v = new_v;
          prev_v = new_v;
// Visit incoming edges of current_v in clockwise order between prev
// (inclusive) and the opposite of current_ei (exclusive) and replace
// current_v with its new index in each corresponding face. This is only
// required if new_v actually differs from the original current_v. However, in
// COTIGUOUS_BOUNDARY mode, we need to distinguish between old and new
// indices, and we mark new ones by adding disk_vertices.size() to prevent
// ambiguities. Otherwise, there may be some new boundary vertex ID matching
// an old inner vertex ID, making it impossible to decide further down which
// inner indices need to be remapped.
          if (new_v != current_v || contiguous_mode) {
               auto updated_index = contiguous_mode ? new_v + num_disk_vertices : new_v;
               LOG_DEBUG(7, "  - updating index new_v=%d -> updated_index=%d", new_v, updated_index);
               do {
                    ++prev;
                    const auto t = make_triangle_index(prev.e()), i = make_edge_offset(prev.e());
                    LOG_DEBUG(9, "    edge %d -> triangle %d, offset %d (recorded index: %d)", prev.e(), t, i, disk_indices[3*t+i]);
                    disk_indices[3*t + i] = updated_index;
               } while (prev.flip().opp() != current_ei);
          } else
               prev.e(current_ei).flip();    // match while condition above
// Note: prev is now expected to point to the opposite of current_ei.
          LOG_DEBUG(7, "  - boundary_info[%d] <- (side=%d, offset=%d, opposite=%d)", current_ei, current_side, bi-current_segment_index, prev.e());
          boundary_info.emplace(std::piecewise_construct,
               std::forward_as_tuple(current_ei),
               std::forward_as_tuple(current_side, bi - current_segment_index, prev.e()));
     }
     LOG_DEBUG(5, "  fixing disk_vertices[%d].next_v = %d", prev_v, first_v);
     disk_vertices[prev_v].next_v = first_v;

// All boundary vertices have been set up at this point. Now iterate over all
// faces again and set up any remaining inner vertices. For inner vertices have
// kept their IDs (KEEP_INNER mode), we only need to sweep once over the
// actual vertices and adjust their org_id properties. In CONTIGUOUS_BOUNDARY
// mode, we have to visit all faces and adjust IDs there, too.
     if (!contiguous_mode) {
          LOG_DEBUG(4, "  fixing org_id properties of remaining inner vertices");
          for (hpindex vi = 0; vi < num_disk_vertices; vi++) {
               auto& ref_v = disk_vertices[vi];
               if (ref_v.org_id == std::numeric_limits<hpindex>::max())
                    ref_v.org_id = vi;
          }
     } else {
// Reserve memory for new_id+1 buckets, the number of inner vertices (started
// with first_v equal to the total number of vertices - 1, and decremented for
// each new vertex on the boundary).
          LOG_DEBUG(4, "  remapping inner vertices");
          std::unordered_map<hpindex, hpindex> inner_map(new_id+1);
          for (auto& ti : disk_indices) {
               if (ti >= num_disk_vertices) {
                    LOG_DEBUG(7, "    fixing index %d -> %d", ti, ti - num_disk_vertices);
                    ti -= num_disk_vertices;
               } else {
                    auto result = inner_map.emplace(ti, new_id);
                    LOG_DEBUG(7, "    org_id <- %d, ti <- new_id=%d", ti, result.first->second);
                    if (result.second) {
                         LOG_DEBUG(7, "    - assigning new ID %d (was: %d)", new_id, ti);
                         disk_vertices[new_id].org_id = ti;
                         --new_id;
                    }
                    ti = result.first->second;
               }
          }
     }
// Color transfer
     for (auto& v : disk_vertices)
          detail::vertex_transfer_color(v, source_mesh.getVertex(v.org_id));
     return std::make_tuple(
          make_triangle_mesh<DiskVertex, Format::DIRECTED_EDGE>(disk_vertices, disk_indices),
          std::move(boundary_info)
     );
}    // }}} cut_to_disk()

/// Fix positions of boundary if \p disk_mesh according to \p coord_builder.
///
/// Additional information about boundary edges is provided in
/// \p boundary_info, which provides cached information about where on which
/// side of the fundamental region a particular edge is located.
template<class DiskMesh, class Coords>
void
set_boundary_constraints(     // {{{
     DiskMesh& disk_mesh,
     const CutGraph::BoundaryInfo& boundary_info,
     Coords&& coord_builder
) {
     LOG_DEBUG(3, "-- set_boundary_constraints(): boundary of length %d", boundary_info.size());
     auto edge = make_edge_walker(disk_mesh);
     for (const auto& info : boundary_info) {
          edge.e(info.first);
          auto& ref_v = disk_mesh.getVertex(edge.u());
          LOG_DEBUG(4, "  edge %d [=(%d,%d)]: side=%d, offset=%d", edge.e(), edge.u(), edge.v(), info.second.side, info.second.offset);
          const auto coord = coord_builder(info.first, info.second.side, info.second.offset);
          ref_v.position.x = coord.x;
          ref_v.position.y = coord.y;
     }
}    // }}}

/// Compute mean-value coordinate vertex \p vi with respect to neighbor \p wi.
///
/// This computation is the straight-forward application of the formula found
/// on page 10 in in the SIGGRAPH Asia 2008 Course Notes [1] or in the survey
/// on surface parameterization by Floater and Hormann [2].
///
/// \note [1] SIGGRAPH Asia 2008 Course Notes, Chapter 2.3,
/// \note [2] Michael S. Floater, Kai Hormann. Surface Parameterization: a
///   Tutorial and Survey. Advances in Multiresolution for Geometric Modeling
///   (2005), pp. 157--186.
///
/// \b TODO
/// - numeric stability of this method of computation? Are area-based
///   computations any better than those relying on evaluating tan()?
/// - are there more efficient ways of batch-computing mean-value coordinates
///   instead of computing them one-by-one?
template<class _Mesh>
class MeanValueCoordinateGenerator {    // {{{
private:
     const _Mesh& m_mesh;
     const EdgeWalker<_Mesh> m_edge_walker;
public:
     using Mesh = _Mesh;
     MeanValueCoordinateGenerator(const Mesh& _mesh) : m_mesh {_mesh}, m_edge_walker {_mesh} { }
     ~MeanValueCoordinateGenerator() { }

     double operator()(hpindex vi, hpindex wi) const noexcept {
          const auto edge {m_edge_walker().uv(vi, wi)};
          const auto v {m_mesh.getVertex(vi)}, w {m_mesh.getVertex(wi)},
               vl {m_mesh.getVertex(edge().next().v())}, vr {m_mesh.getVertex(edge().flip().next().v())};
          const auto vec_vw = w.position-v.position, vec_vl = vl.position-v.position, vec_vr = vr.position-v.position;
          const double len_vw = glm::length(vec_vw);
          const double alpha_vw = std::acos(glm::dot(vec_vw, vec_vr) / (len_vw * glm::length(vec_vr))),
               beta_wv  = std::acos(glm::dot(vec_vw, vec_vl) / (len_vw * glm::length(vec_vl)));
          return (std::tan(alpha_vw/2) + std::tan(beta_wv/2))/ glm::length(vec_vw);
     }
};   // }}} class MeanValueCoordinateGenerator

/// Convenience method for constructing a MeanValueCoordinateGenerator
template<class _Mesh>
MeanValueCoordinateGenerator<_Mesh>
make_mv_coord_generator(const _Mesh& _mesh) { return {_mesh}; }

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
///   double coord_builder(hpindex v, hpindex w);
/// \endcode
/// where \p v and \p w are vertex indices in the (uncut) source mesh,
/// obtained from extra vertex data stored with vertices in \p disk_mesh. The
/// return value is expected to be the computed barycentric/convex-combination
/// "coordinate" of \p v with respect to \p w (but need not be normalized;
/// normalization is carried out automatically).
///
/// Vertices are expected to carry the following extra member data:
/// - org_id [hpindex]: the vertex index in the original mesh
/// - free function bool \c is_inner_vertex(), which takes a mesh reference
///   and vertex data reference and returns true if it refers to an inner
///   vertex, and false for vertices on the boundary
///
/// \returns A tuple of two matrices: the first matrix is quadratic and
/// represents coefficients for all inner vertices; the second one is
/// rectangular with the same number of rows as the first (number of inner
/// vertices) and one column per boundary vertex. Coefficients correspond to
/// the \e normalized coordinates; this matrix can then be used to compute the
/// right-hand vector when solving for a particular coordinate component.
///
/// \note Note: \c is_inner_vertex() is not used in the current
/// implementation, because the same information can be derived from testing
/// for membership in \c v_boundary.
template<class DiskMesh, class Coords>
auto
compute_embedding_coeff( // {{{
     const DiskMesh& disk_mesh,
     std::unordered_map<hpindex, hpindex> v_inner,
     std::unordered_map<hpindex, hpindex> v_boundary,
     Coords&& coord_builder
) {
     using InnerMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
     using BoundaryMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
     const auto num_inner = v_inner.size(), num_boundary = v_boundary.size();
     LOG_DEBUG(4, "compute_embedding_coeff(): %d inner and %d boundary nodes", num_inner, num_boundary);
     LOG_DEBUG(5, "compute_embedding_coeff(): %d vertices, %d indices", disk_mesh.getVertices().size(), disk_mesh.getIndices().size());
// Count inner and boundary neighbors for each inner vertex
     using Eigen::ArrayXi;
     ArrayXi deg_inner {ArrayXi::Constant(num_inner, 1)};   // one on the diagonal
     ArrayXi deg_boundary {ArrayXi::Zero(num_inner)};
     for (const auto& iv : v_inner) {
          visit_spokes(disk_mesh, iv.first,
               [i=iv.second, &v_boundary, &deg_inner, &deg_boundary](const auto& e) {
                    if (v_boundary.count(e.vertex) == 0)
                         deg_inner(i)++;
                    else
                         deg_boundary(i)++;
               });
     }
// Use degrees to accurately reserve memory for sparse matrices.
     LOG_DEBUG(5, "- reserving memory for coeff_inner(%d, %d)", num_inner, num_inner);
     InnerMatrix coeff_inner(num_inner, num_inner);
     coeff_inner.reserve(deg_inner);
     LOG_DEBUG(5, "- reserving memory for coeff_boundary(%d, %d)", num_inner, num_boundary);
     BoundaryMatrix coeff_boundary(num_inner, num_boundary);
     coeff_boundary.reserve(deg_boundary);
     auto col_less = [](auto a, auto b) { return a.col() < b.col(); };

// Create reverse lookup table mapping row indices to inner vertex IDs.
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
          LOG_DEBUG(6, "  - inner node [%d]: -> disk vertex #%d", i, v);
          std::vector<Triplet> inner_row, boundary_row;
          inner_row.reserve(deg_inner(i));
          boundary_row.reserve(deg_boundary(i));
          double row_sum = 0.0;
          const auto& ref_v = disk_vertices[v];
          visit_spokes(disk_mesh, v,
               [&](const auto& e) {
                    auto w = e.vertex;
                    const auto& ref_w = disk_vertices[w];
// vw_coeff are the mean-value coordinates computed for the vertex pair (v,w)
// in the original mesh
                    LOG_DEBUG(6, "  coord_builder(%d, %d)", ref_v.org_id, ref_w.org_id);
                    double vw_coeff = coord_builder(ref_v.org_id, ref_w.org_id);
                    row_sum += vw_coeff;
                    if (v_boundary.count(w) == 0) {
                         LOG_DEBUG(7, "    to inner node #%d (vertex #%d)", v_inner[w], w);
                         inner_row.emplace_back(i, v_inner[w], -vw_coeff);
                    } else {
                         LOG_DEBUG(7, "    to boundary node #%d (vertex #%d)", v_boundary[w], w);
                         boundary_row.emplace_back(i, v_boundary[w], vw_coeff);
                    }
               });
          inner_row.emplace_back(i, i, 1.0);
          std::sort(begin(inner_row), end(inner_row), col_less);
          std::sort(begin(boundary_row), end(boundary_row), col_less);
          for (const auto t : inner_row) {
               if (t.col() != i) {
                    LOG_DEBUG(9, "    - coeff_inner insert:    (%d, %d) <- %.4f / %.4f", t.row(), t.col(), t.value(), row_sum);
                    coeff_inner.insert(t.row(), t.col()) = t.value() / row_sum;
               } else {
                    LOG_DEBUG(9, "    - coeff_inner insert:    (%d, %d) <- %.4f", t.row(), t.col(), t.value());
                    coeff_inner.insert(t.row(), t.col()) = t.value();      // insert 1.0
               }
          }
          for (const auto t : boundary_row) {
               LOG_DEBUG(9, "    - coeff_boundary insert: (%d, %d) <- %.4f / %.4f", t.row(), t.col(), t.value(), row_sum);
               coeff_boundary.insert(t.row(), t.col())  = t.value() / row_sum;
          }
     }
     coeff_inner.makeCompressed();
     coeff_boundary.makeCompressed();
     return std::make_tuple(std::move(coeff_inner), std::move(coeff_boundary));
}
// }}} compute_embedding_coeff()

/// Compute Tutte embedding for \p disk_mesh, with the boundary identified by
/// additional information provided in the vertex data.
///
/// \param[in] disk_mesh disk topology mesh with fixed boundary
/// \param[in] coord_builder callable object for computing coordinates
///
/// Extra information expected to be stored with each vertex:
/// - org_id [hpindex]: the vertex index in the original mesh
/// - free function bool is_inner_vertex(), which takes a mesh reference and
///   vertex data reference and returns true if it refers to an inner vertex,
///   and false for vertices on the boundary
///
/// This interface is provided by DiskVertex.
///
/// TODO:
/// - circumvent v_inner/v_boundary maps if we are guaranteed a disk mesh with
///   a mapping according to CutGraph::VertexMapping::CONTIGUOUS_BOUNDARY
/// - re-think about how we tell compute_disk_embedding which (barycentric)
///   coordinates to use, and where we fix the boundary
/// - what about methods with "free-floating" boundaries, how/when do they
///   work? There are some pointers in the related work section of e.g. Choi
///   et al. [1], who reference a method developed by Hormann et al. [2],
///   which according to table 1, p. 5 in [1] is bijective with free boundary.
///
/// \note [1] Gary Pui-Tung Choi, Lok Ming Lui. A linear formulation for disk
///   conformal parameterization of simply-connected open surfaces, 2017.
/// \note [2] Hormann, K., Greiner, G.. MIPS: An efficient global
///   parametrization method. Curve Surf. Des. 153– 162 (2000)
template<class DiskMesh, class Coords>
void
compute_disk_embedding(DiskMesh& disk_mesh, Coords&& coord_builder) {   // {{{
// Build maps of inner and boundary vertices.
     auto& disk_vertices {disk_mesh.getVertices()};
     std::unordered_map<hpindex, hpindex> v_inner, v_boundary;
     for (hpindex vi = 0, num_verts = disk_vertices.size(); vi < num_verts; vi++) {
          const auto& v {disk_vertices[vi]};
          if (is_inner_vertex(disk_mesh, v))
               v_inner.emplace(vi, v_inner.size());
          else
               v_boundary.emplace(vi, v_boundary.size());
     }
     const auto num_boundary = v_boundary.size();

// Build matrix coefficients
// - compute matrix coefficients based on disk_mesh, which represents the
//   disk topology after cutting open the source mesh along a cut graph
// - the second matrix returned is used to compute right-hand sides for the
//   linear system of equations to be solved: every row corresponds to an inner
//   vertex, and each column to a vertex on the boundary; coefficients at (i,j)
//   are the normalized mean-value coordinates of inner vertex i with regard to
//   vertices of its 1-neighborhood, for boundary vertex j.
     using std::get;
     auto coeff_mat {compute_embedding_coeff(disk_mesh, v_inner, v_boundary, std::forward<Coords>(coord_builder))};
     LOG_DEBUG(3, "  building coordinate columns for boundary vertices");
     Eigen::MatrixX2d boundary_coords(num_boundary, 2);
     for (const auto& bv : v_boundary) {
          const auto& vpos {disk_vertices[bv.first].position};
          boundary_coords(bv.second, 0) = vpos.x;
          boundary_coords(bv.second, 1) = vpos.y;
     }
     LOG_DEBUG(3, "  computing right-hand sides for x and y coordinates");
     Eigen::MatrixX2d rhs {get<1>(coeff_mat)*boundary_coords};
// Create solver with column-major storage (required!)
     Eigen::SparseLU<Eigen::SparseMatrix<double>> tutte_solver;
     tutte_solver.compute(get<0>(coeff_mat));
     if (tutte_solver.info() != Eigen::Success)
          throw std::runtime_error("tutte_solver reported error computing LU decomposition");
     Eigen::MatrixX2d coords {tutte_solver.solve(rhs)};
     if (tutte_solver.info() != Eigen::Success)
          throw std::runtime_error("tutte_solver reported error trying to solve for coordinates");
// Transfer computed coordinates of inner vertices
     for (const auto& iv : v_inner) {
          auto& v {disk_vertices[iv.first]};
          v.position.x = coords(iv.second, 0);
          v.position.y = coords(iv.second, 1);
     }
}
// }}} compute_disk_embedding()

/// Computes the projective structure corresponding to an embedding of a mesh
/// into the Klein disk model.
///
/// \param[in] disk_mesh 2D mesh with disk topology, representing the original
///   mesh after cutting it open. This mesh is expected to have exactly one
///   boundary, the cut circuit, and extra vertex data as in DiskVertex
/// \param[in] boundary_info pairing information for boundary edges, generated by
///   cut_to_disk()
/// \param[in] fp_mesh triangle mesh of fundamental region
/// \param[in] fp_neighbors
/// \param[in] fp_transisions transition maps corresponding to fp_neighbors
///
/// We require the disk_mesh to be an embedding without any degenerate
/// triangles to ensure that transition maps are well-defined.
///
/// TODO
/// - save some computations and visit every (undirected) edge only once
template<class DiskMesh, class FundamentalMesh>
std::vector<hpreal>
make_projective_structure(    // {{{
     const DiskMesh& disk_mesh,
     const CutGraph::BoundaryInfo& boundary_info,
     const FundamentalMesh& fp_mesh,
     const std::vector<hpuint> fp_neighbors,
     const std::vector<hpreal>& fp_transitions
) {
     using FrameMatrix = glm::highp_dmat3;
// Iterate over all faces of disk_mesh and compute transition maps for each
// (half-)edge.  This is easy for non-boundary edges, where we transition
// between faces in the interior of the fundamental region. Whenever we step
// over its boundary, we need to find the proper preimage in the copy of the
// fundamental region in the universal cover.
//
// We use the "macro" mesh fp_mesh and maps in fp_transitions (elements of the
// Fuchsian group G) in order to take a point of the fundamental region to
// another in its G-orbit and compute the corresponding coordinate transform.
     const auto& indices {disk_mesh.getIndices()};
     const auto& macro_vertices {fp_mesh.getVertices()};
     const auto& macro_indices {fp_mesh.getIndices()};
     auto macro_frame = [&macro_vertices, &macro_indices](auto _side) {
          return FrameMatrix(
               macro_vertices[macro_indices[3*_side+0]].position,
               macro_vertices[macro_indices[3*_side+1]].position,
               macro_vertices[macro_indices[3*_side+2]].position);
     };
//
// Each triplet in `indices` corresponds to a single face.  We store a compact
// representation of transition maps (one column of the corresponding matrix)
// per ordered pair of adjacent triangles, which uniquely identifies one
// half-edge, so we store 3 scalars/half-edge.
     std::vector<hpreal> transitions(3*indices.size());
// Pre-compute frames of each triangle flipped over the boundary.
     const auto num_segments = fp_mesh.getNumberOfTriangles();
     std::vector<FrameMatrix> flipped_frame;
     flipped_frame.reserve(num_segments);
// Transition matrix `tr` below corresponds to R*U*S, where R and S are
// permutations of rows and columns, respectively; R swaps the first and third
// rows of U, while S swaps the first and second column. This arrangement is
// compatible with frames given in the order (v0,v1,v2) for the known frame,
// and (w,v2,v1) for the other: We use the appropriate transition map to
// express w in terms of the frame (v1,v0,v2).
     for (unsigned side = 0; side < num_segments; side++) {
          const auto neigh = fp_neighbors[3*side+1];
          FrameMatrix tr(
               {fp_transitions[9*neigh + (3*1 + 2)], 0.0, 0.0},
               {fp_transitions[9*neigh + (3*1 + 1)], 0.0, 1.0},
               {fp_transitions[9*neigh + (3*1 + 0)], 1.0, 0.0});
          flipped_frame.push_back(macro_frame(side) * tr);
     }
// compute_map(): (v0,v1) is the common edge; computes column of transition
// map that sends v2 onto vertex w of the adjacent triangle (v1,v0,w).
// Note: vertices v0, v1, v2 are given in counter-clockwise order.
// TODO: reverse disk projection and find intersections with base points of
// our macro tetrahedra.
     auto edge_walker {make_edge_walker(disk_mesh)};
     auto compute_map = [&](auto v0, auto v1, auto v2, auto tr) {
          FrameMatrix frame;
          //LOG_DEBUG(3, "  compute_map: frame (%d, %d, %d) -> (%d, %d, %d)", v1, v0, v2, v1, w, v0);
          const auto& ref_v0 = disk_mesh.getVertex(v0);
          const auto& ref_v1 = disk_mesh.getVertex(v1);
          auto maybe_e = make_edge_index(disk_mesh, v0, v1);
          assert(!!maybe_e);
          auto common_ei = *maybe_e;
// Simple case: Directly use coordinates if we stay within the same "copy" of
// our fundamental region, which is the case if the edge (v0,v1) is not on the
// boundary (i.e. not a cut edge). Note that it is not enough to just test if
// the common edge has both endpoints on the boundary: The edge might connect
// vertices on two different sides of the fundamental polygon.
          if (!boundary_info.count(common_ei)) {
               auto w = edge_walker().uv(v1, v0).next().v();
               frame = FrameMatrix(ref_v1.position, disk_mesh.getVertex(w).position, ref_v0.position);
          } else {
// Boundary case: We have to switch to a transformation of our fundamental
// region to find the proper preimage coordinates of w.
// Finding vertex w: The index of the "opposite" e' of edge (v0,v1) (on the
// paired side) is conveniently stored in the vertex data as rev_v0.opposite.
// Then w is obtained as the target of the edge following e'.
               const auto& edge_info = boundary_info.at(common_ei);
               assert(boundary_info.count(edge_info.opposite) > 0);
               auto w = edge_walker().e(edge_info.opposite).next().v();
               const auto s_here = edge_info.side;
               const auto s_neigh = boundary_info.at(edge_info.opposite).side;
               // express w in frame coords of s_neigh, flip over boundary, evaluate
               auto w_coords {glm::inverse(macro_frame(s_neigh)) * disk_mesh.getVertex(w).position};
               frame = FrameMatrix(ref_v1.position, flipped_frame[s_here] * w_coords, ref_v0.position);
               //LOG_DEBUG(3, "  - mapping over boundary edge: s_here=%d, s_neigh=%d", s_here, s_neigh);
          }
          auto tr_v2 {glm::inverse(frame) * disk_mesh.getVertex(v2).position};
          for (unsigned i = 0; i <= 2; i++)
               tr[i] = tr_v2[i];
          //LOG_DEBUG(1, "  - %s [sum=%.4f] coords of %s in frame %s",
          //     strmatrix(tr_v2), (tr_v2(0)+tr_v2(1)+tr_v2(2)),
          //     strmatrix(vector_from_vertex(disk_vertices, v2)), strmatrix(frame));
     };
// Iterate over triplets of indices, each representing a single face.
     unsigned triangle = 0;
     visit_triplets(indices, [&](auto v0, auto v1, auto v2) {
          auto tr_it = std::begin(transitions) + 9*triangle;
          //LOG_DEBUG(1, "- [#%4d] triplet (%d, %d, %d)", triangle, v0, v1, v2);
          compute_map(v0, v1, v2, tr_it + 0);
          compute_map(v1, v2, v0, tr_it + 3);
          compute_map(v2, v0, v1, tr_it + 6);
          triangle++;
     });
     return transitions;
}
// }}} make_projective_structure()
}    // namespace obi
}    // namespace happah
#ifdef DISKEMBEDDING_H__LOG_DEBUG
#undef DISKEMBEDDING_H__LOG_DEBUG
#undef LOG_DEBUG
#endif
#endif // #ifdef DISKEMBEDDING_H
// vim:ai:bs=2:fo=croq:expandtab:ts=5:sw=5:sbr=+++\ :lbr:bri:wrap
