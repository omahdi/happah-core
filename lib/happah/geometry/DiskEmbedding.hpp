/// \file DiskEmbedding.hpp
/// \brief Methods for cutting surface meshes into disks, producing Tutte
/// embeddings and computing projective structures from them.
///
/// \author Obada Mahdi <omahdi@gmail.com>
/// \copyright Copyright 2017 Obada Mahdi <omahdi@gmail.com>
/// \copyright Distributed under the Boost Software License, Version 1.0.
/// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <happah/format/hph.hpp>
#include <happah/math/HyperbolicSpace.hpp>
#include <happah/math/ProjectiveStructure.hpp>
#include <happah/geometry/Vertex.hpp>
#include <happah/geometry/TriangleMesh.hpp>
#include <happah/geometry/TriangleGraph.hpp>
#include <happah/util/visitors.hpp>

#include <boost/serialization/nvp.hpp>

//#define TRANSITION_MAPS_USE_TETRAHEDRAL_PROJECTION
//#define GREEDY_REMOVE_CHORDS
#ifndef MPS_DUMP_TRANSITIONS
#define _MPS_DUMP_TRANSITIONS false
#else
#define _MPS_DUMP_TRANSITIONS true
#endif

/// Default \c LOG_DEBUG macro, disables debug messages
#ifndef LOG_DEBUG
#define DISKEMBEDDING_HPP__LOG_DEBUG
#define LOG_DEBUG(...)   { void(0); }
#endif

/// namespace happah
namespace happah {
// {{{ -- namespace detail
namespace detail {
using happah::make_indices;
template<class Vertex>
Indices make_indices(const TriangleMesh<Vertex>& graph) { return graph.getIndices(); }

void vertex_transfer_color(...) { }
template<class TV, class SV>
auto
vertex_transfer_color(TV& _t, const SV& _s) -> std::enable_if_t<sizeof(_t.color) && sizeof(_s.color)> {
     _t.color = _s.color;
}
}    // namespace detail
// }}}1 -- namespace detail

/// An iterator-like class for navigating along edges of a (triangle) mesh.
// {{{ -- EdgeWalker: navigating along edges
template<class _Mesh>
class EdgeWalker;

/// Specialization of EdgeWalker for directed-edge triangle meshes
template<class _Vertex>
class EdgeWalker<TriangleGraph<_Vertex>> {
     static constexpr hpindex NIL_INDEX {std::numeric_limits<hpindex>::max()};
public:
     using Vertex = _Vertex;
     using Mesh = TriangleGraph<Vertex>;

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
     void _set_ei(hpindex _ei) noexcept {
          m_e = &m_mesh->getEdge(_ei);
          m_ei = _ei;
     }
     hpindex _u() const   { return m_mesh->getEdge(m_e->opposite).vertex; }
     hpindex _v() const   { return m_e->vertex; }
     hpindex _opp() const { return m_e->opposite; }
     hpindex _next() const   { return m_e->next; }
     hpindex _prev() const   { return m_e->previous; }
     EdgeWalker& _move_opp() { _set_ei(_opp()); return *this; }
     EdgeWalker& _move_next() { _set_ei(_next()); return *this; }
     EdgeWalker& _move_prev() { _set_ei(_prev()); return *this; }

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
     EdgeWalker& next() { check_init(); return _move_next(); }
/// Move backward
     EdgeWalker& prev() { check_init(); return _move_prev(); }
/// Move to opposite
     EdgeWalker& flip() { check_init(); return _move_opp(); }
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
template<class _V>
EdgeWalker<TriangleGraph<_V>>
make_edge_walker(const TriangleGraph<_V>& _mesh) { return {_mesh}; }

template<class _V>
EdgeWalker<TriangleGraph<_V>>
make_edge_walker(const TriangleGraph<_V>& _mesh, hpindex _ei) { return {_mesh, _ei}; }
// }}}1 -- EdgeWalker: navigating along edges

/// 2D vertex for embeddings with disk topology.
///
/// In addition to a regular 2D vertex with absolute position, this type
/// stores additional information about side pairings and the relation to the
/// original (closed-topology) mesh.
/// FIXME: temporarily changing to 3D vertex; make_projective_structure
/// currently depends on it
/// FIXME: color vs. no color? Base vertex should probably be configurable
class DiskVertex : public VertexPC<Space2D> {    // {{{1
public:
     using Point = Point2D;   ///< type for storing the vertex position

private:
     using BaseVertex = VertexPC<Space2D>;
     using BasePoint = typename std::decay_t<decltype(std::declval<BaseVertex>().position)>;
     static constexpr hpindex IS_INNER_VALUE {std::numeric_limits<hpindex>::max()};

     //static inline BasePoint _build_point(Point2D _p) { return {_p.x, _p.y, 1.0}; }
     //static inline BasePoint _build_point(Point3D _p) { return _p; }
     static inline BasePoint _build_point(Point2D _p) { return _p; }
     static inline BasePoint _build_point(Point3D _p) { return {_p.x, _p.y}; }

public:
/// Vertex index in the original mesh.
     hpindex org_id {0};
/// Vertex index of the following boundary vertex. Only meaningful if
/// is_inner() is false, undefined otherwise (an implementation-defined magic
/// value is used as indicator).
     hpindex next_v {IS_INNER_VALUE};
/// Default-constructs a vertex at the origin, marked as inner vertex.
     DiskVertex() : BaseVertex(_build_point(Point2D{0.0, 0.0})) {}
     ~DiskVertex() = default;
     DiskVertex(const DiskVertex&) = default;
     DiskVertex(DiskVertex&&) = default;
     DiskVertex& operator=(const DiskVertex&) = default;
     DiskVertex& operator=(DiskVertex&&) = default;
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
};   // }}}1 class DiskVertex

using DiskMesh = TriangleGraph<DiskVertex>;

/// Additional information for boundary edges required by
/// make_projective_structure.
struct BoundaryEdgeInfo {     // {{{1
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
};   // }}}1 struct BoundaryEdgeInfo

/// Representation of a cut graph as a topologically sorted list of half-edges
///
/// Edges are sorted according to the order of a walk along the Eulerian
/// circuit that will become the boundary of the resulting disk after cutting.
/// Edges are stored as a compact "array of arrays," implemented as a flat
/// array of data with an additional array storing index pointers, including
/// an extra entry at the end to allow easy computation of lengths without
/// distinction for the last item).
class CutGraph {    // {{{1
     static const unsigned long _classversion = 0;
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
     using pairings_type = std::vector<hpindex>;

     struct BranchNodeInfo {
          static const unsigned long _classversion = 0;
          hpindex next, prev;
          unsigned degree {0};
          hpindex u, vo, vi;

          template<class S>
          friend S& operator<<(S& _s, const BranchNodeInfo& _v) {
               using ::happah::format::hph::operator<<;
               _s << _v.next << ' ' << _v.prev << ' ' << _v.degree << ' ' << _v.u << ' ' << _v.vo << ' ' << _v.vi;
               return _s;
          }
          template<class S>
          friend S& operator>>(S& _s, BranchNodeInfo& _v) {
               using ::happah::format::hph::operator>>;
               _s >> _v.next >> _v.prev >> _v.degree >> _v.u >> _v.vo >> _v.vi;
               return _s;
          }
          template<class S>
          void save(S& _s, unsigned long version) const {
               using boost::serialization::make_nvp;
               _s << make_nvp("next", next);
               _s << make_nvp("prev", prev);
               _s << make_nvp("degree", degree);
               _s << make_nvp("u", u);
               _s << make_nvp("vi", vi);
               _s << make_nvp("vo", vo);
          }
          template<class S>
          void load(S& _s, unsigned long version) {
               using boost::serialization::make_nvp;
               _s >> make_nvp("next", next);
               _s >> make_nvp("prev", prev);
               _s >> make_nvp("degree", degree);
               _s >> make_nvp("u", u);
               _s >> make_nvp("vi", vi);
               _s >> make_nvp("vo", vo);
          }
     };
/// List of (half-)edge indices, topologically sorted into an Eulerian circuit
/// of the cut graph. When built using cut_graph_from_edges() or
/// cut_graph_from_paths(), the first edge is guaranteed to originate at a
/// branch node.
     circuit_type m_circuit;
/// Offsets into m_circuit of first edge of each segment; an additional entry
/// at index segment_count() is included to simplify computations of segment
/// lenghts or index ranges without distinction of cases for the last segment.
     segments_type m_segments;
/// Array of side pairings: \c m_pairings[i] is the index of the side
/// associated with side \c i. It is a permutation of the range
/// <tt>[0..segment_count(*this)-1]</tt>.
     pairings_type m_pairings;
/// Maps a segment index to an instance of \c BranchNodeInfo, which records
/// the (undirected) \p degree of a branch node in the cut graph, the index
/// into \c m_node_info of the \p next one corresponding to the same
/// vertex with index \p u, and the neighbor \p v at the outgoing edge.
     std::vector<BranchNodeInfo> m_node_info;
/// The topological genus of the source mesh.
     unsigned m_genus {0};
/// The number of boundary components of the source mesh.
/// FIXME: not implemented
     unsigned m_boundary_components {0};

// Check if cut graph is usable (non-empty)
     void check() const {
#if !defined(NDEBUG)
          if (m_segments.size() == 0 || m_circuit.size() == 0)
               throw std::runtime_error("CutGraph: not holding with any cut segments, possibly uninitialized");
#endif
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
     CutGraph(const CutGraph& _g) :
          m_circuit{_g.m_circuit},
          m_segments{_g.m_segments},
          m_pairings{_g.m_pairings},
          m_node_info{_g.m_node_info},
          m_genus{_g.m_genus},
          m_boundary_components{_g.m_boundary_components}
     { }
     CutGraph(CutGraph&& _g) :
          m_circuit{std::move(_g.m_circuit)},
          m_segments{std::move(_g.m_segments)},
          m_pairings{std::move(_g.m_pairings)},
          m_node_info{std::move(_g.m_node_info)},
          m_genus{_g.m_genus},
          m_boundary_components{_g.m_boundary_components}
     { }
     ~CutGraph() { }
     CutGraph& operator=(const CutGraph& _g) {
          m_circuit = _g.m_circuit;
          m_segments = _g.m_segments;
          m_pairings = _g.m_pairings;
          m_node_info = _g.m_node_info;
          m_genus = _g.m_genus;
          m_boundary_components = _g.m_boundary_components;
     }
     CutGraph& operator=(CutGraph&& _g) {
          m_circuit = std::move(_g.m_circuit);
          m_segments = std::move(_g.m_segments);
          m_pairings = std::move(_g.m_pairings);
          m_node_info = std::move(_g.m_node_info);
          m_genus = std::move(_g.m_genus);
          m_boundary_components = std::move(_g.m_boundary_components);
     }

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
/// Returns the number of branch nodes in the cut graph.
     friend hpindex branch_node_count(const CutGraph& cut_graph) noexcept {
          return 1 + segment_count(cut_graph)/2 - 2*cut_graph.m_genus;
     }
/// Given a polygonal fundamental domain \f$D\f$ corresponding to a cut, the
/// branch node degree corresponds to the number of translates of \f$D\f$
/// meeting at this vertex.
     friend hpindex branch_node_degree(const CutGraph& cut_graph, hpindex _i) noexcept {
          cut_graph.check();
          assert(_i < segment_count(cut_graph));
          return cut_graph.m_node_info[_i].degree;
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
/// Returns the index of the side paired with side <tt>_i</tt> in
/// \p cut_graph.
     friend auto paired_side(const CutGraph& cut_graph, hpindex _i) noexcept {
          cut_graph.check();
          assert(_i < segment_count(cut_graph));
          return cut_graph.m_pairings[_i];
     }

// see definition for documentation
     template<class Mesh>
     friend CutGraph cut_graph_from_edges(const Mesh& source_mesh, const std::vector<hpindex>& cut);
// see definition for documentation
     template<class Mesh>
     friend bool remove_chords(CutGraph& cut_graph, const Mesh& source_mesh);

     template<class S>
     friend S& operator<<(S& _s, const CutGraph& _v) {
          using ::happah::format::hph::operator<<;
          _s << _v.m_circuit << ' ';
          _s << _v.m_segments << ' ';
          _s << _v.m_pairings << ' ';
          _s << _v.m_node_info << ' ';
          _s << _v.m_genus;
          return _s;
     }
     template<class S>
     friend S& operator>>(S& _s, CutGraph& _v) {
          using ::happah::format::hph::operator>>;
          _s >> _v.m_circuit;
          _s >> _v.m_segments;
          _s >> _v.m_pairings;
          _s >> _v.m_node_info;
          _s >> _v.m_genus;
          return _s;
     }

     template<class S>
     void save(S& _s, unsigned long version) const {
          using boost::serialization::make_nvp;
          _s << make_nvp("circuit", m_circuit);
          _s << make_nvp("segments", m_segments);
          _s << make_nvp("pairings", m_pairings);
          _s << make_nvp("node_info", m_node_info);
          _s << make_nvp("genus", m_genus);
     }
     template<class S>
     void load(S& _s, unsigned long version) {
          using boost::serialization::make_nvp;
          _s >> make_nvp("circuit", m_circuit);
          _s >> make_nvp("segments", m_segments);
          _s >> make_nvp("pairings", m_pairings);
          _s >> make_nvp("node_info", m_node_info);
          _s >> make_nvp("genus", m_genus);
     }
};   // }}}1 class CutGraph

/// Build a CutGraph from an (unordered) list of edge indices.
template<class SourceMesh>
CutGraph
cut_graph_from_edges(const SourceMesh& source_mesh, const std::vector<hpindex>& cut) {    // {{{1
     using std::begin;
     using std::end;
     if (cut.size() == 0)
          throw std::invalid_argument("An empty graph does not have a circuit.");
     CutGraph cut_graph;
     decltype(cut_graph.m_circuit) circuit;
     decltype(cut_graph.m_segments) segments;
     constexpr auto NIL_VALUE = std::numeric_limits<std::decay_t<decltype(cut_graph.m_pairings[0])>>::max();

// Note: for closed meshes only!
//
// Determine topological genus based on Euler's formula for closed graphs.
// TODO: handle meshes with boundary
// FIXME: num_vertices needs to be the count of vertices actually referenced
// by the num_triangles triangles. See note in Triangle{Graph,Mesh}.hpp
     const auto num_vertices = source_mesh.getNumberOfVertices();
     const auto num_triangles = source_mesh.getNumberOfTriangles();
     if ((num_triangles % 2) != 0)
          throw std::runtime_error("Number of triangles is expected to be even for closed meshes without boundary.");
     if (((num_triangles/2 - num_vertices) % 2) != 0)
          throw std::runtime_error("``num_triangles/2 - num_vertices'' must be even.");
     cut_graph.m_genus = 1 + (num_triangles/2 - num_vertices) / 2;

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
     } while (edge.e() != start_edge);
     circuit.shrink_to_fit();

     const auto num_segments = std::accumulate(begin(branch_nodes), end(branch_nodes), 0u,
          [&vertex_count](auto s, auto v) { return s+vertex_count[v]; });
     segments.reserve(num_segments+1);

// Consistency checks based on Euler's formula
     if ((num_segments % 2) != 0)
          throw std::runtime_error("Invalid cut graph: expected an even number of cut segments");
     if (num_segments < 4*cut_graph.m_genus)
          throw std::runtime_error("Invalid cut graph: too few cuts for computed genus of source mesh");
     if (branch_nodes.size() != 1 + num_segments/2 - 2*cut_graph.m_genus)
          throw std::runtime_error("Invalid cut graph: number of branch nodes inconsistent with Euler's formula");

     decltype(cut_graph.m_pairings) pairings(num_segments, NIL_VALUE);
     decltype(cut_graph.m_node_info) branch_node_info(num_segments);
     std::vector<hpindex> last_branch_node(branch_nodes.size(), NIL_VALUE);
     LOG_DEBUG(5, "- found %d cut nodes: %s", branch_nodes.size(), utils::str(branch_nodes));
     LOG_DEBUG(5, "- computing positions of cut nodes in %u-gon", num_segments);
     unsigned int first_branch_offset = 0;
     const auto circuit_size = circuit.size();
     for (hpindex i = 0; i < circuit_size; i++) {
          edge_walker.e(circuit[i]);
          const auto segment_index = segments.size();
          const auto branch_node_it {std::lower_bound(begin(branch_nodes), end(branch_nodes), edge_walker.u())};
          if (branch_node_it != end(branch_nodes) && *branch_node_it == edge_walker.u()) {
               if (segment_index == 0)
                    first_branch_offset = i;
               assert(segments.size() < num_segments && "assertion segments.size() < num_segments failed");
               segments.emplace_back(i-first_branch_offset);
               auto& last_index = last_branch_node[std::distance(begin(branch_nodes), branch_node_it)];
               if (last_index != NIL_VALUE) {
                    branch_node_info[last_index].next = segment_index;
                    branch_node_info[segment_index].prev = last_index;
                    branch_node_info[segment_index].degree = 1 + branch_node_info[last_index].degree;
               } else {
                    branch_node_info[segment_index].prev = NIL_VALUE;
                    branch_node_info[segment_index].degree = 1;
               }
               last_index = segment_index;
               branch_node_info[segment_index].u = edge_walker.u();
               branch_node_info[segment_index].vo = edge_walker.v();
               branch_node_info[segment_index].vi = edge_walker.e(circuit[i > 0 ? i-1 : circuit_size-1]).u();
          }
     }
     segments.emplace_back(circuit.size());
     LOG_DEBUG(5, "  segment offsets: %s", segments);
// Rotate m_circuit such that first element is an edge emanating from the
// first branch node seen.
     std::rotate(circuit.begin(), circuit.begin()+first_branch_offset, circuit.end());
// Pre-compute side pairings; also fix up doubly-linked list.
// FIXME: maybe simplify loop at the expense of a few more iterations. Number
// of branch nodes is bounded by the topological genus within a constant
// factor anyway.
     hpindex branch_node_index = 0;
     for (auto last_index : last_branch_node) {
          branch_node_index++;
          const auto branch_degree = branch_node_info[last_index].degree;
          hpindex prev_index = last_index, segment_index;
          do {
               segment_index = prev_index;
               auto& info = branch_node_info[segment_index];
               info.degree = branch_degree;
// Look for paired segment among remaining copies of this branch node by
// walking doubly-linked list backwards from segment_index.
               for (auto scan_index = info.prev;
                    pairings[segment_index] == NIL_VALUE && scan_index != NIL_VALUE;
                    scan_index = branch_node_info[scan_index].prev
               ) {
                    if (branch_node_info[scan_index].vi != info.vo)
                         continue;
                    const auto paired_index = scan_index > 0 ? scan_index-1 : num_segments-1;
                    pairings[paired_index] = segment_index;
                    pairings[segment_index] = paired_index; // also forces loop condition false and breaks loop
               }
               prev_index = branch_node_info[segment_index].prev;
          } while (prev_index != NIL_VALUE);
          branch_node_info[segment_index].prev = last_index;
          branch_node_info[last_index].next = segment_index;
     }
     circuit.swap(cut_graph.m_circuit);
     segments.swap(cut_graph.m_segments);
     pairings.swap(cut_graph.m_pairings);
     branch_node_info.swap(cut_graph.m_node_info);
     return cut_graph;
}    // }}}1 cut_graph_from_edges()

/// Builds an array of indices mapping the index of a cut path segment to the
/// index of its paired counterpart.
///
/// \note The current implementation builds an array from pre-computed values
/// stored internally by CutGraph.
/// - space and time: linear in the number of segments (precomputed)
inline auto
segment_pairings(const CutGraph& cut_graph) {  // {{{1
     const auto num_segments = segment_count(cut_graph);
     std::vector<hpindex> pairings(num_segments);
     for (unsigned i = 0; i < num_segments; i++)
          pairings[i] = paired_side(cut_graph, i);
     return pairings;
};   // }}}1 segment_pairings()

/// Build array of neighbors of triangle fan representing our fundamental
/// region circumscribed by \p cut_graph.
///
/// \param[in] cut_graph
/// \param[in] source_mesh
///
/// \returns An array of length <tt>3*cut_graph.segment_count()</tt> holding
/// the indices of faces adjacent to face \a j at positions
/// <tt>3*j..3*j+2</tt>. We follow the convention that adjacent triangles are
/// given in counter-clockwise order. Assuming that each face is defined as
/// the vertex triple <tt>(v0, v1, v2)</tt> and \p v0 is the common vertex of
/// the macro fan, neighbors of triangle \a i in a fan of size \a n are given
/// in the order <tt>((i-1) mod n, opposite(i), (i+1) mod n))</tt>.
inline std::vector<hpindex>
make_neighbors(const CutGraph& cut_graph) { // {{{1
     const auto num_segments = segment_count(cut_graph);
     std::vector<hpindex> neighbors(3*num_segments);

     for (auto i = 0u; i < num_segments; i++) {
          neighbors[3*i + 0] = (i+num_segments-1) % num_segments;
          neighbors[3*i + 1] = (paired_side(cut_graph, i)+num_segments) % num_segments;
          neighbors[3*i + 2] = (i+1) % num_segments;
     }
     return neighbors;
}    // }}}1 make_neighbors()

/// Builds a mesh representation of a fundamental region corresponding to the
/// cut graph in the projective disk model.
///
/// The construction is based on a proof from [1; Thm. 7.16.2]
///
/// For a cut graph with \a n branch nodes and associated degrees <tt>d(i)
/// (i=1,...,n)</tt>, the result is a fan-shaped, triangulated <tt>n</tt>-gon
/// with center vertex at the origin and interior angles <tt>2*M_PI /
/// d(i)</tt> at corner \a i. Note that angles are to be measured in the
/// projective disk model of hyperbolic space: they are chosen such that
/// exactly <tt>d(i)</tt> copies of the fundamental region meet at corner \a i
/// without gaps or overlaps in the hyperbolic plane. Note also that because
/// of the properties of hyperbolic space, scaling the mesh will also change
/// these angles!
///
/// \note [1] Alan Beardon. The Geometry of Discrete Groups. (1983)
inline TriangleMesh<VertexP2>
make_fundamental_domain(const CutGraph& cut_graph) { // {{{1
// Collect information about number and valences of branch nodes and build a
// suitable polygon in the projective disk model.
//
// segment_count(CutGraph)
// - number of corners
// branch_node_degree(CutGraph, hpindex)
// - number of polygons meeting at a corner
// - quotient of 2*M_PI and degree gives the desired interior angle
     const auto num_segments = segment_count(cut_graph);
     std::vector<hpreal> thetas;
     thetas.reserve(num_segments);
     for (unsigned k = 0; k < num_segments; k++)
          thetas.emplace_back((2*M_PI) / branch_node_degree(cut_graph, k));
     //const auto corners {hyp_polygon_from_angles_P(thetas)};
     const auto corners {make_convex_polygon(thetas)};
     std::vector<VertexP2> vertices;
     vertices.reserve(num_segments+1);
     std::vector<hpindex> indices;
     indices.reserve(3*num_segments);
     vertices.emplace_back(VertexP2({0.0, 0.0}));
     for (unsigned k = 0; k < num_segments; k++) {
          vertices.emplace_back(hyp_CtoP(corners[k]));
          indices.emplace_back(0);
          indices.emplace_back(1 + k);
          indices.emplace_back(1 + ((k+1) % num_segments));
     }
     return make_triangle_mesh(std::move(vertices), std::move(indices));
} // }}}

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
remove_chords(CutGraph& cut_graph, const Mesh& mesh) { // {{{1
     using std::begin;
     using std::end;
     const auto num_segments = segment_count(cut_graph);
     auto& circuit {cut_edges(cut_graph)};
     std::vector<bool> processed(num_segments);
     LOG_DEBUG(5, "remove_chords(): processing %d sides", num_segments);
     auto edge_walker {make_edge_walker(mesh)};

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
     const auto get_range_containing = [num_segments, &edge_walker, &cut_graph, &sorted_insert](auto v, auto ei) {
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
          const auto other_index = paired_side(cut_graph, s_index);
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
               LOG_DEBUG(6, " -relaxing edges at current_node=%d", current_node);
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
                    const auto d = path_nodes[current_node].distance;
                    assert(pos != current_node);
                    if (pos < current_node)
                         continue;
#ifndef GREEDY_REMOVE_CHORDS
                    if (path_nodes[pos].visited && path_nodes[pos].distance <= d+1) {
                         next_node = current_node+1;        // indicator, see comment below
                         continue;
                    }
#endif
                    LOG_DEBUG(7, "  updating pos=%d: distance=%d, rev_edge=%d", pos, d+1, edge_walker.opp());
                    path_nodes[pos].visited = 1;
                    path_nodes[pos].distance = d+1;
                    path_nodes[pos].prev = current_node;
                    path_nodes[pos].rev_edge = edge_walker.opp();
#ifdef GREEDY_REMOVE_CHORDS
                    if (next_node < pos)
                         next_node = pos;
#else
// Could move this out of the do-while loop, but we keep it here as an
// indicator that at least one edge has been relaxed, which is verified
// below (next_node == current_node).
                    next_node = current_node+1;
#endif
               } while (edge_walker.flip().next().e() != relax_begin);
               if (next_node == current_node)
                    throw std::runtime_error("remove_chords(): could not find next edge in cut segment");
          }
          processed[other_index] = true;
          const int delta = s_length - path_nodes[next_node].distance;
// Sanity check (shouldn't happen)
          assert(delta >= 0 && "Internal error: path cannot grow longer");
          if (delta == 0)
               continue;      // no improvement
          LOG_DEBUG(3, "  segment %d shortened by %d edges", s_index, delta);
          did_remove = true;
// FIXME untested!
// Trace back over shortest path, updating circuit.
          auto fp = s_begin + path_nodes[s_length].distance;
          auto rp = begin(circuit) + segment_index(cut_graph, other_index);
          for (; next_node != 0; next_node = path_nodes[next_node].prev) {
               const auto rev_edge = path_nodes[next_node].rev_edge;
               const auto fwd_edge = edge_walker.e(rev_edge).opp();
               *(--fp) = fwd_edge;
               *(rp++) = rev_edge;
          }
// m_circuit layout: [...i..j..i+1..//..o..k..o+1...]
//   [i,i+1): i-th segment before; [i..j): after taking shortcuts
//   [o,o+1): o-th segment before; [o..k): after taking shortcuts
// Ranges [j, i+1) and [k, o+1) are to be erased; erase range starting at k
// first, keeping indices for first intact.
          assert(s_index < other_index);
          circuit.erase(rp, begin(circuit)+segment_index(cut_graph, other_index)+segment_length(cut_graph, other_index));
          circuit.erase(s_end-delta, s_end);
// Adjust segment indices: every range up to and including "other_index" moved
// by "delta"; everything following "other_index" by twice that, including the
// extra element at the end of m_segments that marks the length of the
// complete circuit (m_circuit.size()).
          for (unsigned j = s_index+1; j <= other_index; j++)
               cut_graph.m_segments[j] -= delta;
          for (unsigned j = other_index+1; j < cut_graph.m_segments.size(); j++)
               cut_graph.m_segments[j] -= 2*delta;
     }
     return did_remove;
}    // }}}1 remove_chords()

/// Compute transition matrices for macro tetrahedra corresponding to a
/// tesselation {3; p,6,6}, where p is the valence of the center vertex.
std::vector<hpreal>
compute_transitions_p_3(std::vector<hpindex> neighbors) {   // {{{1
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
/// - different versions for both mesh formats (TriangleGraph, TriangleMesh)
/// - complexity analysis
template<class Mesh>
std::tuple<DiskMesh, typename CutGraph::BoundaryInfo>
cut_to_disk(   // {{{1
     const CutGraph& cut_graph,
     const Mesh& source_mesh,
     const CutGraph::VertexMapping mapping_mode = CutGraph::VertexMapping::CONTIGUOUS_BOUNDARY
) {
     using namespace std::string_literals;
     using std::to_string;
     const bool contiguous_mode = mapping_mode == CutGraph::VertexMapping::CONTIGUOUS_BOUNDARY;
     const auto& circuit = cut_edges(cut_graph);
     const auto circuit_length = circuit.size();
     CutGraph::BoundaryInfo boundary_info(circuit_length);
     const auto num_triangles {source_mesh.getNumberOfTriangles()};
// TODO: note that we need the number of vertices in the graph, whereas the
// triangle graph instance may store more than are actually "used". See note
// in TriangleGraph.h.
     const auto num_vertices {source_mesh.getNumberOfVertices()};
// Start with copy of source_mesh face indices; its size does not change, even
// though technically we will have one additional face (the generally
// non-triangular "outer" face), which is not represented explicitly.
     auto disk_indices {detail::make_indices(source_mesh)};
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
//   (num_triangles/2) - num_vertices
// is equal to 2g-2, where g is the genus of source_mesh. This also follows
// from Euler's formula, taking into account the number of faces of the
// triangulated closed source_mesh (num_triangles) and the number of
// (undirected) edges (num_triangles*3 / 2).
     assert((circuit_length % 2) == 0);
     assert((num_triangles % 2) == 0);
     const int extra_vertices = 1 + (circuit_length/2) + (num_triangles/2) - num_vertices;
     assert(extra_vertices > 0);
     LOG_DEBUG(3, "  mesh stats: %d faces, %d vertices ==> genus %d",
          num_triangles, num_vertices, (2 + num_triangles/2 - num_vertices)/2);
     LOG_DEBUG(3, "  %d source vertices, %d extra vertices ==> %d disk vertices total",
          num_vertices, extra_vertices, num_vertices+extra_vertices);
// Initialize each vertex with a magic value for org_id in order to be able to
// detect uninitialized inner vertices later (Note: marked is_inner by default)
     const auto num_disk_vertices = num_vertices + extra_vertices;
     std::vector<DiskVertex> disk_vertices(num_disk_vertices, {{}, std::numeric_limits<hpindex>::max()});
// v_mapped keeps track of whether a particular source vertex index has
// already been seen while walking along the cut.
     std::unordered_set<hpindex> v_mapped(circuit_length/2);
     const hpindex first_v = num_vertices + extra_vertices - 1;
     hpindex new_id = first_v;
     LOG_DEBUG(5, "  first new ID: %d", first_v);
// prev_v holds the previous allocated vertex ID for updating its next_v
// member; it is not used in the first iteration, because we do not know the
// new ID of the last vertex along the cut yet.
     hpindex prev_v = 0;
     unsigned current_side = 0;
     hpindex current_segment_index = segment_index(cut_graph, current_side);
     hpindex current_segment_length = segment_length(cut_graph, current_side);
     auto prev = make_edge_walker(source_mesh);
     for (hpindex bi = 0; bi < circuit_length; bi++) {
          if (bi == current_segment_index+current_segment_length) {
               ++current_side;
               current_segment_index = bi;
               current_segment_length = segment_length(cut_graph, current_side);
          }
          auto current_ei = circuit[bi];
// Initialize edge walker "prev" with predecessor of current_ei in circuit
// and store the vertex ID it is pointing at in current_v.
          const auto prev_ei = circuit[bi > 0 ? bi-1 : circuit_length-1];
          auto current_v  = prev.e(prev_ei).v();
// Consistency check: preceding edge must point to the tail vertex of the
// current edge current_ei == circuit[bi].
          if (current_v != prev().e(current_ei).u()) {
               LOG_DEBUG(3, "Error: inconsistency in cut circuit at position %d: previous edge %d points to node #%d, expected #%d", bi, prev_ei, current_v, prev().e(current_ei).u());
               LOG_DEBUG(3, "Error: prev_ei = (%d -> %d), current_ei = (%d -> %d)", prev.u(), prev.v(), prev().e(current_ei).u(), prev().e(current_ei).v());
          }
          assert(current_v == prev().e(current_ei).u());
          const bool mapped = v_mapped.count(current_v) > 0;
          hpindex new_v = 0;
          if (!mapped)
               v_mapped.emplace(current_v);
          LOG_DEBUG(5, "  current_ei <- %d <- circuit[%d], side=%d, current_v=%d (mapped=%s)", current_ei, bi, current_side, current_v, mapped);
// Assign new ID if current_v has been seen before, or if a contiguous mapping
// of boundary vertices is requested.
          if (mapped || contiguous_mode)
               new_v = new_id--;
          else
               new_v = current_v;
          LOG_DEBUG(5, "    new_v <- %d", new_v);
          disk_vertices[new_v].org_id = current_v;
          if (bi != 0)        // last vertex is updated after the loop
               disk_vertices[prev_v].next_v = new_v;
          prev_v = new_v;
// Enumerate incoming edges of current_v in clockwise order between prev
// (inclusive) and the opposite of current_ei (exclusive) and replace
// current_v with its new index in each corresponding face, which is stored at
// position of the edge offset of the *following* edge with respect to that
// face.  This is only required if new_v actually differs from the original
// current_v. However, in CONTIGUOUS_BOUNDARY mode, we need to distinguish
// between old and new indices, and we mark new ones by adding
// disk_vertices.size() to prevent ambiguities. Otherwise, there may be some
// new boundary vertex ID matching an old inner vertex ID, making it
// impossible to decide further down which inner indices need to be remapped.
          if (new_v != current_v || contiguous_mode) {
               const auto updated_index = contiguous_mode ? new_v + num_disk_vertices : new_v;
               LOG_DEBUG(7, "  - updating index new_v=%d -> updated_index=%d", new_v, updated_index);
               do {
                    ++prev;   // [new interface?] prev.move_next()
                    const auto t = make_triangle_index(prev.e()), i = make_edge_offset(prev.e());
                    LOG_DEBUG(9, "    edge %d -> triangle %d, offset %d (recorded index: %d)", prev.e(), t, i, disk_indices[3*t+i]);
                    disk_indices[3*t + i] = updated_index;
// Bail out if we make a round-trip around current_v without finding
// current_ei. This should not happen after we ensured that prev_ei and
// current_ei meet at the same vertex at the beginning of the for loop.
                    if (prev.opp() == prev_ei)
                         throw std::runtime_error("Error: Inconsistent cut circuit at position "s+to_string(bi)+": no outgoing edge of vertex #"s+to_string(current_v)+" matched expected edge #"s+to_string(current_ei)+" at this position"s);
               } while (prev.flip().opp() != current_ei);   // [new interface?] prev.move_opp().opp()
          } else
               prev.e(current_ei).flip();    // match while condition above; [new interface?] prev.e(current_ei).move_opp()
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
                         assert(new_id != hpindex(-1));
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
          make_triangle_graph<DiskVertex>(disk_vertices, disk_indices),
          std::move(boundary_info)
     );
}    // }}}1 cut_to_disk()

/// Fix positions of boundary if \p disk_mesh according to \p coord_builder.
///
/// Additional information about boundary edges is provided in
/// \p boundary_info, which provides cached information about where on which
/// side of the fundamental region a particular edge is located.
template<class DiskMesh, class Coords>
void
set_boundary_constraints(     // {{{1
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
/// in the SIGGRAPH Asia 2008 Course Notes [1; Sec. 2.3, p. 10] or in the
/// survey on surface parameterization by Floater and Hormann [2].
///
/// \note [1] Kai Hormann, Konrad Polthier, Alla Sheffer. Mesh
/// Parameterization: Theory and Practice. (2008)
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
class MeanValueCoordinateGenerator {    // {{{1
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
};   // }}}1 class MeanValueCoordinateGenerator

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
compute_embedding_coeff( // {{{1
     const DiskMesh& disk_mesh,
     std::unordered_map<hpindex, hpindex> v_inner,
     std::unordered_map<hpindex, hpindex> v_boundary,
     Coords&& coord_builder
) {
     using InnerMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
     using BoundaryMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
     const auto num_inner = v_inner.size(), num_boundary = v_boundary.size();
     LOG_DEBUG(4, "compute_embedding_coeff(): %d inner and %d boundary nodes", num_inner, num_boundary);
     LOG_DEBUG(5, "compute_embedding_coeff(): %d vertices, %d faces", disk_mesh.getVertices().size(), disk_mesh.getNumberOfTriangles());
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
          visit_spokes(make_spokes_enumerator(disk_mesh.getEdges(), disk_mesh.getOutgoing(v)),
               [&] (auto ei) {
                    const auto w = disk_mesh.getEdge(ei).vertex;
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
// }}}1 compute_embedding_coeff()

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
/// - reconsider how we tell compute_disk_embedding which (barycentric)
///   coordinates to use, and where we fix the boundary
/// - what about methods with "free-floating" boundaries, how/when do they
///   work? There are some pointers in the related work section of e.g. Choi
///   et al. [1], who reference a method developed by Hormann et al. [2],
///   which according to table 1, p. 5 in [1] is bijective with free boundary.
///
/// \note [1] Gary Pui-Tung Choi, Lok Ming Lui. "A linear formulation for disk
/// conformal parameterization of simply-connected open surfaces." (2017)
/// \note [2] Kai Hormann, Günther Greiner. "MIPS: An efficient global
///   parametrization method." Curve and Surface Design (2000), pp. 153–162.
template<class DiskMesh, class Coords>
void
compute_disk_embedding(DiskMesh& disk_mesh, Coords&& coord_builder) {   // {{{1
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
     LOG_DEBUG(5, "  building coordinate columns for boundary vertices");
     Eigen::MatrixX2d boundary_coords(num_boundary, 2);
     for (const auto& bv : v_boundary) {
          const auto& vpos {disk_vertices[bv.first].position};
          boundary_coords(bv.second, 0) = vpos.x;
          boundary_coords(bv.second, 1) = vpos.y;
     }
     LOG_DEBUG(5, "  computing right-hand sides for x and y coordinates");
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
// }}}1 compute_disk_embedding()

/// Computes the projective structure corresponding to an embedding of a mesh
/// into the Klein disk model from a cut graph \p cut_graph, a corresponding
/// disk embedding \p disk_mesh with boundary information and a mesh representing the
/// fundamental region used for mapping boundary vertices.
///
/// \param[in] cut_graph the cut graph that was used to cut a closed manifold
///   mesh into the disk represented by \p disk_mesh
/// \param[in] boundary_info pairing information for boundary edges, generated by
///   cut_to_disk()
/// \param[in] disk_mesh 2D mesh with disk topology, representing the original
///   mesh after cutting it open; this mesh is expected to have exactly one
///   boundary, the cut circuit, and extra vertex data as provided by DiskVertex
/// \param[in] fp_mesh representation of a fundamental region as triangle mesh
///   in the projective disk model; it is used to compute deck transformations
///   and it is thus required that the boundary mapping used in the embedding
///   procedure placed boundary vertices onto matching boundary edges of this
///   mesh representation
///
/// Notes:
/// - We require the disk_mesh to be an embedding without any degenerate
///   triangles to ensure that transition maps are well-defined.
/// - The fundamental polygon mesh \p fp_mesh is expected to be a triangle fan
///   with vertex index 0 being the center vertex, and each face being
///   specified with indices [0, u, v] in counter-clockwise order.
///
/// TODO
/// - save some computations and visit every (undirected) edge only once
/// - declare fp_mesh with template parameter type?
inline
ProjectiveStructure
make_projective_structure( // {{{1
     const CutGraph& cut_graph,
     const CutGraph::BoundaryInfo& boundary_info,
     const DiskMesh& disk_mesh,
     const TriangleMesh<VertexP2>& fp_mesh
) {
     constexpr bool mps_debug_output = _MPS_DUMP_TRANSITIONS;
     const auto num_segments = segment_count(cut_graph);
     const auto num_faces = disk_mesh.getNumberOfTriangles();
// We rebuild ``neighbors'' from the disk_mesh triangle graph and side pairing
// information from ``boundary_info''.
     std::vector<hpindex> neighbors;
     neighbors.reserve(3*num_faces);
     std::vector<hpreal> transitions;
     transitions.reserve(9*num_faces);
// Compute deck transformations based on the given fundamental polygon and the
// side pairing obtained from the cut graph. The matrix at fuchsians[i] is the
// transformation that takes the paired side s_i' to s_i.
     using mat3 = glm::dmat3;
     std::vector<mat3> fuchsians;
     fuchsians.reserve(num_segments);
     for (unsigned k = 0; k < num_segments; k++) {
// In order for the transformation to be orientation-preserving, we have to
// flip the order of one pair of vertices. Note that in the corresponding
// polygonal schema, paired edges represent loops with opposite direction.
          const auto kp = paired_side(cut_graph, k);
          const auto& p1 = fp_mesh.getVertex(k, 2).position;
          const auto& q1 = fp_mesh.getVertex(k, 1).position;
          const auto& p2 = fp_mesh.getVertex(kp, 1).position;
          const auto& q2 = fp_mesh.getVertex(kp, 2).position;
// Compute map that takes edge (p2,q2) to (p1,q2)
          fuchsians.emplace_back(hyp_decktrans_P<mat3>(p2, q2, p1, q1));
     }
/* Consider consecutive triangles S, T, U in a triangle fan, where the
   double-stroked edge (v1,v2) is the boundary the fundamental region.
               c'
          .....o.....
              / \
       .     /   \     .
        .   /  T' \   .
         . /       \ .
       v2 o=========o v1
         / \       / \
        /   \  T  /   \
       /     \   /     \
      /   U   \ /   S   \
     o_________o_________o
   v3 .       .c.       . v0
       .     .   .     .

   Goal: Compute transitions phiTX for (half-)edges of each triangle T,
   transforming coordinates in the chart U(T) into coordinates in U(X). The
   matrix representation of phiTX depends on ordered bases for each chart,
   which we define by the following convention:

   For transition map phiTX, we use basis {w, v, u} for   |       x
   U(T) and {w, x, v} for U(X). Here, u, v, w, and x are  |       o
   the coordinate vectors of the corresponding points in  |      / \
   real projective 2-space (RP^2=R^3). Conversely, for    |     / X \
   phiXT we use the basis {v, w, w} for U(X) and          |  w o-----o v
   {v, u, w} for U(T).                                    |     \ T /
                                                          |      \ /
   If all four points are known, transition maps can be   |       o
   derived as follows:                                    |       u

   phiTS: U(T) -> U(S),
     phiTS  = inv([v1 | v0 | c ])*[v1 | c  | v2]
   phiTU: U(T) -> U(U),
     phiTU  = inv([c  | v3 | v2])*[c  | v2 | v1]
   phiTT': U(T) -> U(T'),
     phiTT' = inv([v2 | c' | v1])*[v2 | v1 | c ]

   Coordinates for the first two transition maps are readily available. For
   phiTT', we have to find the image c' of c under the deck transformation
   that takes the paired edge of S((v1, v2)) = (v1', v2') to (v1, v2). This
   transformation is a Fuchsian group element and can be computed by finding
   the Moebius transformation that takes v1' to v1 and v2' to v2. The paired
   side S((v1, v2)) is known from the computation of the cut graph.
*/
     using FrameMatrix = glm::dmat3;
     using Vec3 = typename FrameMatrix::col_type;
     //using Vec2 = glm::tvec2<typename FrameMatrix::value_type>;
     //using Vec2 = glm::tvec2<typename FrameMatrix::value_type,
     //     detail::glm_vec_traits<FrameMatrix::col_type>::precision>;
     if (mps_debug_output)
          std::cout << "{\n  \"transitions\": {\n";
     auto push_tr = [&transitions] (auto frame, auto w) {
          const auto coord {glm::inverse(frame) * w};
          transitions.emplace_back(coord.x);
          transitions.emplace_back(coord.y);
          transitions.emplace_back(coord.z);
          if (mps_debug_output) {
               std::cout << coord.x << "," << coord.y << "," << coord.z << "]";
               const auto sum = coord.x + coord.y + coord.z;
               if (std::abs(std::abs(sum) - 1.0) > 1e-4)
                    std::cout << " /* |sum| = |" << sum << "| > 1.0 */";
          }
     };
// compute_map(): (v0,v1) is the common edge; computes column of transition
// map that sends v2 onto vertex w of the adjacent triangle (v1,v0,w).
// Note: vertices v0, v1, v2 are given in counter-clockwise order.
     auto edge_walker {make_edge_walker(disk_mesh)};
     auto compute_map = [&](auto common_ei, auto v0, auto v1, auto v2) {
// assertion: common_ei == *make_edge_index(v0, v1)
          FrameMatrix frame;
          //LOG_DEBUG(3, "  compute_map: frame (%d, %d, %d) -> (%d, %d, %d)", v1, v0, v2, v1, w, v0);
          const auto& ref_v0 = disk_mesh.getVertex(v0);
          const auto& ref_v1 = disk_mesh.getVertex(v1);
// Simple case: Directly use coordinates if we stay within the same "copy" of
// our fundamental region, which is the case if the edge (v0,v1) is not on the
// boundary (i.e. not a cut edge). Note that it is not enough to just test if
// the common edge has both endpoints on the boundary: The edge might connect
// vertices on two different sides of the fundamental polygon.
          if (!boundary_info.count(common_ei)) {
               const auto walker = edge_walker().e(common_ei).flip();
               if (mps_debug_output) {
                    if (neighbors.size() > 0)
                         std::cout << ",\n";
                    std::cout << "    \"" << common_ei << "->" << walker.e() << "\": [";
               }
               neighbors.emplace_back(make_triangle_index(walker.e()));
               const auto w = walker().next().v();
               const auto& w_ref = disk_mesh.getVertex(w);
               //const Vec2 v1_pos(ref_v1.position.x, ref_v1.position.y);
               //const Vec2 w_pos(w_ref.position.x, w_ref.position.y);
               //const Vec2 v0_pos(ref_v0.position.x, ref_v0.position.y);
               const auto v1_pos = Vec3(ref_v1.position.x, ref_v1.position.y, 1.0);
               const auto w_pos = Vec3(w_ref.position.x, w_ref.position.y, 1.0);
               const auto v0_pos = Vec3(ref_v0.position.x, ref_v0.position.y, 1.0);
               // TODO: does it make sense to use higher precision here?
#ifdef TRANSITION_MAPS_USE_TETRAHEDRAL_PROJECTION
               //frame = FrameMatrix(hyp_PtoH(v1_pos), hyp_PtoH(w_pos), hyp_PtoH(v0_pos));
               static_assert(false, "not implemented");
#else
               //frame = FrameMatrix(Vec3(v1_pos.x, v1_pos.y, 1.0), Vec3(w_pos.x, w_pos.y, 1.0), Vec3(v0_pos.x, v0_pos.y, 1.0));
               frame = FrameMatrix(v1_pos, w_pos, v0_pos);
#endif
          } else {
// Boundary case: We have to switch to a translate of our fundamental
// region to find the proper preimage coordinates of w.
// Finding vertex w: The index of the "opposite" e' of edge (v0,v1) (on the
// paired side) is conveniently stored in the vertex data as rev_v0.opposite.
// Then w is obtained as the head of the edge following e'.
               const auto& edge_info = boundary_info.at(common_ei);
               assert(boundary_info.count(edge_info.opposite) > 0);
               if (mps_debug_output) {
                    if (neighbors.size() > 0)
                         std::cout << ",\n";
                    std::cout << "    \"*" << common_ei << "->" << edge_info.opposite << "\": [";
               }
               neighbors.emplace_back(make_triangle_index(edge_info.opposite));
               const auto w = edge_walker().e(edge_info.opposite).next().v();
               const auto& ref_w {disk_mesh.getVertex(w)};
// Walk along the orbit of w using the deck transformation, translating it
// into the "copy" adjacent to edge common_ei in the universal cover. As the
// result is only determined up to a scalar, we normalize by projecting it
// back into projective disk coordinates and then back onto the hyperboloid.
               const auto wprime_tmp = fuchsians[edge_info.side]*FrameMatrix::col_type(ref_w.position.x, ref_w.position.y, 1.0);
               //const Vec2 v1_pos(ref_v1.position.x, ref_v1.position.y);
               //const Vec2 wprime_pos(wprime_tmp.x/wprime_tmp.z, wprime_tmp.y/wprime_tmp.z);
               //const Vec2 v0_pos(ref_v0.position.x, ref_v0.position.y);
               const auto v1_pos = Vec3(ref_v1.position.x, ref_v1.position.y, 1.0);
               const auto wprime_pos = Vec3(wprime_tmp);
               const auto v0_pos = Vec3(ref_v0.position.x, ref_v0.position.y, 1.0);
#ifdef TRANSITION_MAPS_USE_TETRAHEDRAL_PROJECTION
               //frame = FrameMatrix(hyp_PtoH(v1_pos), hyp_PtoH(wprime_pos), hyp_PtoH(v0_pos));
               static_assert(false, "not implemented");
#else
               //frame = FrameMatrix(Vec3(v1_pos.x, v1_pos.y, 1.0), Vec3(wprime_pos.x, wprime_pos.y, 1.0), Vec3(v0_pos.x, v0_pos.y, 1.0));
               frame = FrameMatrix(v1_pos, wprime_pos, v0_pos);
#endif
               //const auto s_here = edge_info.side;
               //const auto s_paired = boundary_info.at(edge_info.opposite).side;
               //LOG_DEBUG(3, "  - mapping over boundary edge: s_here=%d, s_neigh=%d", s_here, s_neigh);
          }
          const auto& ref_v2 = disk_mesh.getVertex(v2);

#ifdef TRANSITION_MAPS_USE_TETRAHEDRAL_PROJECTION
          static_assert(false, "not implemented");
          //push_tr(frame, hyp_PtoH(Vec2(ref_v2.position.x, ref_v2.position.y)));
#else
          push_tr(frame, Vec3(ref_v2.position.x, ref_v2.position.y, 1.0));
#endif
          //LOG_DEBUG(1, "  - %s [sum=%.4f] coords of %s in frame %s",
          //     strmatrix(tr_v2), (tr_v2(0)+tr_v2(1)+tr_v2(2)),
          //     strmatrix(vector_from_vertex(disk_vertices, v2)), strmatrix(frame));
     };
// Iterate over triplets of edges, each representing a single face. Only visit
// num_faces triplets and ignore additional edges that represent the
// "outer" half-edges of the boundary.
     visit_triplets(std::cbegin(disk_mesh.getEdges()), num_faces, 3,
          [&](const auto& e0, const auto& e1, const auto& e2) {
               const hpindex v0 = e2.vertex, v1 = e0.vertex, v2 = e1.vertex;
               //LOG_DEBUG(1, "- [#%4d] triplet (%d, %d, %d)", make_triangle_index(e0), v0, v1, v2);
// Note: the order in which vertex pairs (v0, v1) etc. are processed is
// important and has to match the edge_index passed in the first argument.
// The use of make_edge_index() may not be strictly necessary (we could simply
// keep a counter), but it seems cleaner not to rely on implementation
// details; the interface of TriangleGraphh does not necessarily guarantee
// that getEdges()[i] corresponds to the edge with index i.
// FIXME: We do rely on such details in other places.
               compute_map(make_edge_index(e0), v0, v1, v2);
               compute_map(make_edge_index(e1), v1, v2, v0);
               compute_map(make_edge_index(e2), v2, v0, v1);
          });
     if (mps_debug_output)
          std::cout << "\n  }\n}\n";
     return make_projective_structure(std::move(neighbors), std::move(transitions));
}    // }}}

}    // namespace happah
#ifdef DISKEMBEDDING_HPP__LOG_DEBUG
#undef DISKEMBEDDING_HPP__LOG_DEBUG
#undef LOG_DEBUG
#endif
