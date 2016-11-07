// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <stack>
#include <utility>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/geometries/Vertex.h"
#include "happah/math/Space.h"
#include "happah/math/TriangleRefinementScheme.h"
#include "happah/utils/DeindexedArray.h"
#include "happah/utils/ProjectiveStructureUtils.h"
#include "happah/geometries/TriangleMeshUtils.h"
#include "happah/utils/VertexFactory.h"

//TODO: projective structure refinement validation
//   1. refine projective structure A to get projective structure B
//   2. convert A into triangle mesh a
//   3. convert B into triangle mesh b
//   4. refine a using same refinement scheme to get triangle mesh a'
//   5. check if b and a' are the same mesh

//TODO: reimplement projective structure using a different strategy; instead of specifying neighbors of tetrahedron specify vertex types, their neighbors and how to go about a vertex

namespace happah {

//NOTE: A projective structure is a set of tuples of n reals encoding transitions between simplices in n-space.  These simplices have one corner on the origin.  The n vectors formed by the edges adjacent to the origin form a basis for the n-space, which implies that the simplices are non-degenerate.
template<class Space>
class ProjectiveStructureBase {
public:
     using Border = std::vector<hpuint>;//every two integers specify the two neighbors that share a border edge
     using Indices = std::vector<hpuint>;
     using Neighbors = std::vector<hpuint>;
     using Transition = typename Space::POINT;
     using Transitions = std::vector<Transition>;

     Border getBorder() const { return toIntegers(m_border); }

     /**
      * Returns an index to a transition which belongs to the edge at the given simplex and index.
      *
      * @param[out] The transition index or -1 if no transition exists.
      */
     hpuint getIndex(hpuint simplex, hpuint face) const {
          assert(simplex < getNumberOfSimplices() && face < Space::DIMENSION);
          return m_indices[simplex * Space::DIMENSION + face];
     }

     const Indices& getIndices() const { return m_indices; }

     const Neighbors& getNeighbors() const { return m_neighbors; }

     hpuint getNumberOfSimplices() const { return m_indices.size() / Space::DIMENSION; }

     /**
      * Returns the transition at the specified index.
      */
     const Transition& getTransition(hpuint index) const {
          assert(index < m_transitions.size());
          return m_transitions[index];
     }

     const Transitions& getTransitions() const { return m_transitions; }

protected:
     Indices m_indices;
     Neighbors m_neighbors;
     Transitions m_transitions;
     boost::dynamic_bitset<> m_border;
     
     ProjectiveStructureBase() {}

     //NOTE: The matrices specify the transition functions between simplices.  The transitions encode which matrix to use when moving from a given simplex to one of its neighbors.  There are n transitions per simplex, where n is the dimension of the space in which the projective structure is.  Similarly, there are n neighbors per simplex.  An index in the neighbors array is the simplex' relative index; that is, the transitions of a neighbor with index i are found between ni and (n+1)i-1 inclusive.
     ProjectiveStructureBase(Transitions transitions, Indices indices, Neighbors neighbors, const Border& border) 
          : m_indices(std::move(indices)), m_neighbors(std::move(neighbors)), m_transitions(std::move(transitions)), m_border(toBitset(border)) {}

     ProjectiveStructureBase(Transitions transitions, Indices indices, Neighbors neighbors, boost::dynamic_bitset<> border) 
          : m_indices(std::move(indices)), m_neighbors(std::move(neighbors)), m_transitions(std::move(transitions)), m_border(std::move(border)) {}

     /**
      * @param[in]  mask  Mask specifying the refinement scheme. 
      * @param[out]       Returns the constructed projective structure in which every macro-simplex of this projective structure is subdivided into micro-simplices as specified by the refinement scheme.
      */
     //ProjectiveStructure<Space>* refine(const Mask& mask) { return nullptr; }

private:
     boost::dynamic_bitset<> toBitset(const std::vector<hpuint>& border) const {
          boost::dynamic_bitset<> temp(m_neighbors.size());
          auto i = border.cbegin();
          auto end = border.cend();
          while(i != end) {
               hpuint n0 = *i;
               hpuint n1 = *(++i);
               hpuint o0 = 3 * n0;
               hpuint o1 = 3 * n1;
               auto i0 = m_neighbors.cbegin() + o0;
               auto i1 = m_neighbors.cbegin() + o1;
               if(*i0 == n1) temp.set(o0);
               else if(*(++i0) == n1) temp.set(o0 + 1);
               else {
                    assert(*(++i0) == n1);
                    temp.set(o0 + 2);
               }
               if(*i1 == n0) temp.set(o1);
               else if(*(++i1) == n0) temp.set(o1 + 1);
               else {
                    assert(*(++i1) == n0);
                    temp.set(o1 + 2);
               }
               ++i;
          }
          return std::move(temp);
     }

     std::vector<hpuint> toIntegers(const boost::dynamic_bitset<>& border) const {
          std::vector<hpuint> temp;
          auto index = border.find_first();
          while(index != boost::dynamic_bitset<>::npos) {
               hpuint n0 = index / 3;
               hpuint n1 = *(m_neighbors.cbegin() + index);
               if(n0 < n1) {//avoid duplicates
                    temp.push_back(n0);
                    temp.push_back(n1);
               }
               index = border.find_next(index);
          }
          return std::move(temp);
     }

};//ProjectiveStructureBase

template<class Space>
class ProjectiveStructure;

namespace mode {

     //NOTE: Fans are formed by the simplices adjacent to a vertex.  Neighbors are the n neighbors with one vertex on the origin of a given n-simplex on the origin.  Transitions are the maps that map one simplex to its neighbor.
     enum class ProjectiveStructure3D { BORDER, TETRAHEDRA, TRANSITIONS, VERTICES };

}//namespace mode

namespace view {

     enum class ProjectiveStructure3D { TETRAHEDRA, VERTICES };

}

using ProjectiveStructure3D = ProjectiveStructure<Space3D>;

template<>
class ProjectiveStructure<Space3D> : public ProjectiveStructureBase<Space3D> {
     using Mode = mode::ProjectiveStructure3D;
     using Vertices = std::vector<hpuint>;
     using View = view::ProjectiveStructure3D;

     template<int t_dummy, View t_view, Mode... t_modes>
     class Iterator {
     public:
          using difference_type = hpuint;
          using iterator_category = std::random_access_iterator_tag;
          using value_type = std::tuple<typename Iterator<t_dummy, t_view, t_modes>::value_type...>;
          using pointer = value_type*;
          using reference = value_type&;

          Iterator(const ProjectiveStructure& projectiveStructure, hpuint tetrahedron)
               : m_i(Iterator<t_dummy, t_view, t_modes>(projectiveStructure, tetrahedron)...) {}

          difference_type operator-(const Iterator& iterator) { return std::get<0>(m_i) - std::get<0>(iterator.m_i); }

          Iterator operator+(hpuint offset) const {
               Iterator iterator(*this);
               iterator += offset;
               return iterator;
          }

          Iterator& operator++() {
               preincrement();
               return *this;
          }

          Iterator& operator--() {
               predecrement();
               return *this;
          }

          Iterator& operator+=(hpuint offset) {
               add(offset);
               return *this;
          }

          Iterator operator++(int) {
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          Iterator operator--(int) {
               Iterator iterator(*this);
               --(*this);
               return iterator;
          }

          bool operator!=(const Iterator& iterator) const { return std::get<0>(iterator.m_i) != std::get<0>(m_i); }

          value_type operator*() const { return getValue(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

     private:
          std::tuple<Iterator<t_dummy, t_view, t_modes>...> m_i;

          void add(hpuint offset) { add(offset, std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<unsigned long I0, unsigned long... Is>
          void add(hpuint offset, std::index_sequence<I0, Is...>) { doAdd(offset, std::get<I0>(m_i)); add(offset, std::index_sequence<Is...>()); }

          template<unsigned long I0>
          void add(hpuint offset, std::index_sequence<I0>) { doAdd(offset, std::get<I0>(m_i)); }

          template<Mode mode>
          void doAdd(hpuint offset, Iterator<t_dummy, t_view, mode>& i) { i += offset; }

          template<unsigned long... Is>
          value_type getValue(std::index_sequence<Is...>) const { return std::make_tuple(*(std::get<Is>(m_i))...); }

          BUILD_TUPLE_HANDLER_METHODS(predecrement, doPredecrement)

          void predecrement() { predecrement(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<Mode mode>
          void doPredecrement(Iterator<t_dummy, t_view, mode>& i) { --i; }

          BUILD_TUPLE_HANDLER_METHODS(preincrement, doPreincrement)

          void preincrement() { preincrement(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<Mode mode>
          void doPreincrement(Iterator<t_dummy, t_view, mode>& i) { ++i; }

     };//Iterator

     //NOTE: Tetrahedra are ordered counterclockwise.
     template<int t_dummy>
     class Iterator<t_dummy, View::VERTICES, Mode::TETRAHEDRA> {
     public:
          using value_type = const std::vector<hpuint>&;//TODO: maybe output two iterators (begin, end) Iterator<View::VERTEX, Mode::TETRAHEDRA>?

          Iterator(const ProjectiveStructure& projectiveStructure) 
               : m_last(projectiveStructure.m_neighbors[0]), m_neighbors(projectiveStructure.m_neighbors.cbegin()), m_tetrahedron(0) {
               m_todo.resize(projectiveStructure.m_neighbors.size());
               m_todo.set();
               m_fan = fan();
          }

          Iterator(const ProjectiveStructure& projectiveStructure, hpuint vertex) 
               : m_last(vertex == projectiveStructure.m_nVertices ? UNULL : projectiveStructure.m_neighbors[0]), m_neighbors(projectiveStructure.m_neighbors.cbegin()), m_tetrahedron(vertex == projectiveStructure.m_nVertices ? UNULL : 0) {
               m_todo.resize(projectiveStructure.m_neighbors.size());
               if(vertex < projectiveStructure.m_nVertices) {
                    m_todo.set();
                    m_fan = fan();
                    for(auto i = 0u; i < vertex; ++i) ++(*this);//TODO: create += operator and avoid generating fan in each increment
               }
          }

          bool end() const { return m_tetrahedron == UNULL; }

          Iterator& operator++() {
               if(m_todo.none()) {
                    m_tetrahedron = m_last = UNULL;
                    return *this;
               }
               do {
                    hpuint i = 3 * m_tetrahedron;
                    if(m_todo[i]) {
                         m_last = *(m_neighbors + i);
                         break;
                    } else if(m_todo[i + 1]) {
                         m_last = *(m_neighbors + (i + 1));
                         break;
                    } else if(m_todo[i + 2]) {
                         m_last = *(m_neighbors + (i + 2));
                         break;
                    }
                    ++m_tetrahedron;
               } while(true);
               m_fan = fan();
               return *this;
          }

          Iterator operator++(int) { 
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          bool operator==(const Iterator& iterator) const { return iterator.m_tetrahedron == m_tetrahedron && (m_tetrahedron == UNULL || (iterator.m_last == m_last && iterator.m_todo == m_todo)); }

          bool operator!=(const Iterator& iterator) const { return !(*this == iterator); }

          value_type operator*() const { return m_fan; }

     private:
          std::vector<hpuint> m_fan;
          hpuint m_last;//last tetrahedron we visit in current fan
          const Neighbors::const_iterator m_neighbors;
          hpuint m_tetrahedron;//first tetrahedron we visit in current fan 
          boost::dynamic_bitset<> m_todo;//flags for determining which vertices have been processed

          std::vector<hpuint> fan() {
               std::vector<hpuint> fan;
               hpuint current = m_tetrahedron;
               hpuint previous = m_last;
               do {
                    fan.push_back(current);
                    hpuint i = 3 * current;
                    auto n = m_neighbors + i;
                    hpuint n0 = *n;
                    hpuint n1 = *(++n);
                    hpuint n2 = *(++n);
                    if(previous == n0) {
                         m_todo[i] = false;
                         previous = current;
                         current = n2;
                    } else if(previous == n1) {
                         m_todo[i + 1] = false;
                         previous = current;
                         current = n0;
                    } else {
                         assert(previous == n2);
                         m_todo[i + 2] = false;
                         previous = current;
                         current = n1;
                    }
               } while(!(current == m_tetrahedron && previous == m_last));
               return fan;
          }

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, View::VERTICES, Mode::VERTICES> {
     public:
          using value_type = Vertices;//TODO: maybe output two iterators (begin, end) Iterator<View::VERTEX, Mode::VERTICES>?

          Iterator(const ProjectiveStructure& projectiveStructure, hpuint vertex)
               : m_i(projectiveStructure, vertex), m_vertex(vertex), m_vertices(projectiveStructure.m_vertices.cbegin()) {}

          Iterator& operator++() {
               ++m_i;
               ++m_vertex;
               return *this;
          }

          Iterator operator++(int) { 
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          bool operator==(const Iterator& iterator) const { return iterator.m_vertex == m_vertex; }

          bool operator!=(const Iterator& iterator) const { return !(*this == iterator); }

          value_type operator*() const {
               Vertices ring;
               auto& fan = *m_i;
               ring.reserve(fan.size());
               for(auto f = fan.cbegin(), end = fan.cend(); f != end; ++f) {
                    auto v = m_vertices + 3 * (*f);
                    auto v0 = *v;
                    auto v1 = *(++v);
                    auto v2 = *(++v);
                    if(v0 == m_vertex) ring.push_back(v1);
                    else if(v1 == m_vertex) ring.push_back(v2);
                    else {
                         assert(v2 == m_vertex);
                         ring.push_back(v0);
                    }
               }
               return ring;
          }

     private:
          Iterator<t_dummy, View::VERTICES, Mode::TETRAHEDRA> m_i;
          hpuint m_vertex;
          const Vertices::const_iterator m_vertices;

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, View::TETRAHEDRA, Mode::BORDER> {
     public:
          using difference_type = hpuint;
          using value_type = std::tuple<bool, bool, bool>;

          Iterator(const ProjectiveStructure& projectiveStructure, hpuint tetrahedron) 
               : m_border(projectiveStructure.m_border), m_i(3 * tetrahedron) {}

          difference_type operator-(const Iterator& iterator) const { return (m_i - iterator.m_i) / 3; }

          Iterator operator+(hpuint offset) const {
               Iterator iterator(*this);
               iterator += offset;
               return iterator;
          }

          Iterator& operator++() {
               m_i += 3;
               return *this;
          }

          Iterator& operator--() {
               m_i -= 3;
               return *this;
          }

          Iterator& operator+=(hpuint offset) {
               m_i += 3 * offset;
               return *this;
          }

          Iterator operator++(int) { 
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          Iterator operator--(int) {
               Iterator iterator(*this);
               --(*this);
               return iterator;
          }

          value_type operator[](hpuint offset) {
               auto i = m_i + 3 * offset;
               return std::make_tuple(m_border[i], m_border[i + 1], m_border[i + 2]);
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_i != m_i; }

          value_type operator*() const { return std::make_tuple(m_border[m_i], m_border[m_i + 1], m_border[m_i + 2]); }

     private:
          const boost::dynamic_bitset<>& m_border;
          hpuint m_i;

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, View::TETRAHEDRA, Mode::TETRAHEDRA> {
     public:
          using difference_type = hpuint;
          using value_type = std::tuple<hpuint, hpuint, hpuint>;

          Iterator(const ProjectiveStructure& projectiveStructure, hpuint tetrahedron) 
               : m_neighbor(projectiveStructure.m_neighbors.cbegin() + 3 * tetrahedron) {}

          difference_type operator-(const Iterator& iterator) const { return (m_neighbor - iterator.m_neighbor) / 3; }

          Iterator operator+(hpuint offset) const {
               Iterator iterator(*this);
               iterator += offset;
               return iterator;
          }

          Iterator& operator++() {
               m_neighbor += 3;
               return *this;
          }

          Iterator& operator--() {
               m_neighbor -= 3;
               return *this;
          }

          Iterator& operator+=(hpuint offset) {
               m_neighbor += 3 * offset;
               return *this;
          }

          Iterator operator++(int) { 
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          Iterator operator--(int) {
               Iterator iterator(*this);
               --(*this);
               return iterator;
          }

          value_type operator[](hpuint offset) {
               auto n = m_neighbor + 3 * offset;
               return std::make_tuple(*n, *(n + 1), *(n + 2));
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_neighbor != m_neighbor; }

          value_type operator*() const { return std::make_tuple(*m_neighbor, *(m_neighbor + 1), *(m_neighbor + 2)); }

     private:
          Neighbors::const_iterator m_neighbor;

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, View::TETRAHEDRA, Mode::TRANSITIONS> {
          using ProxyIterator = typename DeindexedArray<Transitions>::const_iterator;

     public:
          using difference_type = typename ProxyIterator::difference_type;
          using value_type = std::tuple<const Transition&, const Transition&, const Transition&>;

          Iterator(const ProjectiveStructure& projectiveStructure, hpuint tetrahedron) 
               : m_i(i(projectiveStructure, tetrahedron)) {}

          difference_type operator-(const Iterator& iterator) const { return (m_i - iterator.m_i) / 3; }

          Iterator operator+(hpuint offset) const {
               Iterator iterator(*this);
               iterator += offset;
               return iterator;
          }

          Iterator& operator++() { 
               m_i += 3;
               return *this;
          }

          Iterator& operator--() { 
               m_i -= 3;
               return *this;
          }

          Iterator& operator+=(hpuint offset) {
               m_i += 3 * offset;
               return *this;
          }

          Iterator operator++(int) { 
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          Iterator operator--(int) {
               Iterator iterator(*this);
               --(*this);
               return iterator;
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_i != m_i; }

          value_type operator*() const { return std::tie(*m_i, *(m_i + 1), *(m_i + 2)); }

     private:
          ProxyIterator m_i;

          static ProxyIterator i(const ProjectiveStructure& projectiveStructure, hpuint tetrahedron) { return deindex(projectiveStructure.m_transitions, projectiveStructure.m_indices).cbegin() + 3 * tetrahedron; }

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, View::TETRAHEDRA, Mode::VERTICES> {
     public:
          using difference_type = hpuint;
          using value_type = std::tuple<hpuint, hpuint, hpuint>;

          Iterator(const ProjectiveStructure& projectiveStructure, hpuint tetrahedron) 
               : m_vertex(projectiveStructure.m_vertices.cbegin() + 3 * tetrahedron) {}

          difference_type operator-(const Iterator& iterator) const { return (m_vertex - iterator.m_vertex) / 3; }

          Iterator operator+(hpuint offset) const {
               Iterator iterator(*this);
               iterator += offset;
               return iterator;
          }

          Iterator& operator++() {
               m_vertex += 3;
               return *this;
          }

          Iterator& operator--() {
               m_vertex -= 3;
               return *this;
          }

          Iterator& operator+=(hpuint offset) {
               m_vertex += 3 * offset;
               return *this;
          }

          Iterator operator++(int) { 
               Iterator iterator(*this);
               ++(*this);
               return iterator;
          }

          Iterator operator--(int) {
               Iterator iterator(*this);
               --(*this);
               return iterator;
          }

          value_type operator[](hpuint offset) {
               auto v = m_vertex + 3 * offset;
               return std::make_tuple(*v, *(v + 1), *(v + 2));
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_vertex != m_vertex; }

          value_type operator*() const { return std::make_tuple(*m_vertex, *(m_vertex + 1), *(m_vertex + 2)); }

     private:
          Vertices::const_iterator m_vertex;

     };//Iterator

     //NOTE: The refiner guarantees that all refined simplices are next to each other in the resulting projective structure but does not guarantee the order in which they appear.
     class Refiner {
     public:
          Refiner(const ProjectiveStructure3D& projectiveStructure, const TriangleRefinementScheme& scheme, hpreal epsilon = EPSILON)
               : m_containsCorners(false), m_scheme(scheme), m_oldBorder(projectiveStructure.m_border), m_oldIndices(projectiveStructure.m_indices), m_oldNeighbors(projectiveStructure.m_neighbors), m_oldTransitions(projectiveStructure.m_transitions), m_nMacroTriangles(m_oldIndices.size() / 3), m_nMicroTriangles(scheme.indices.size() / 3) {
               using hpuint7 = std::array<hpuint, 7>;

               //NOTE: There can be refinements without corner triangles, i.e. Powell-Sabin 6-split.  But if there is at least one corner triangle, then there are three corner triangles because of symmetry.

               auto nNewNeighbors = 3 * m_nMicroTriangles * m_nMacroTriangles;
               m_newBorder.resize(nNewNeighbors);
               m_newNeighbors.reserve(nNewNeighbors);
               m_newIndices.reserve(nNewNeighbors);

               hpuint7 corner01;
               hpuint7 corner12;
               hpuint7 corner20;

               std::vector<hpuint7> edge0;
               std::vector<hpuint7> edge1;
               std::vector<hpuint7> edge2;
               std::vector<hpuint7> inner;

               {//sort triangles
                    auto p0 = 0u, t = 0u;
                    for(auto i = scheme.indices.begin(), end = scheme.indices.end(); i != end; ++i, p0 += 3, ++t) {
                         auto i0 = *i;
                         auto i1 = *(++i);
                         auto i2 = *(++i);

                         auto p1 = p0 + 1;
                         auto p2 = p1 + 1;

                         auto& v0 = m_scheme.points[i0];
                         auto& v1 = m_scheme.points[i1];
                         auto& v2 = m_scheme.points[i2];

                         auto e01 = getEdgeType(v0, v1, epsilon);
                         auto e12 = getEdgeType(v1, v2, epsilon);
                         auto e20 = getEdgeType(v2, v0, epsilon);

                         //NOTE: There are only three ways that the three vertices of a triangle in the mask can be associated with the mask's defining triangle and these correspond to the three rotations of the triangle about its center.
                         if(e01 == 0) {
                              if(e12 == 1) {
                                   assert(e20 == UNULL);
                                   corner01 = { i0, i1, i2, p0, p1, p2, t };
                                   m_containsCorners = true;
                              } else if(e20 == 2) {
                                   assert(e12 == UNULL);
                                   corner20 = { i2, i0, i1, p2, p0, p1, t };
                                   m_containsCorners = true;
                              } else {
                                   assert(e12 == UNULL && e20 == UNULL);
                                   edge0.push_back({ i0, i1, i2, p0, p1, p2, t });
                              }
                         } else if(e01 == 1) {
                              if(e12 == 2) {
                                   assert(e20 == UNULL);
                                   corner12 = { i0, i1, i2, p0, p1, p2, t };
                                   m_containsCorners = true;
                              } else if(e20 == 0) {
                                   assert(e12 == UNULL);
                                   corner01 = { i2, i0, i1, p2, p0, p1, t };
                                   m_containsCorners = true;
                              } else {
                                   assert(e12 == UNULL && e20 == UNULL);
                                   edge1.push_back({ i0, i1, i2, p0, p1, p2, t });
                              }
                         } else if(e01 == 2) {
                              if(e12 == 0) {
                                   assert(e20 == UNULL);
                                   corner20 = { i0, i1, i2, p0, p1, p2, t };
                                   m_containsCorners = true;
                              } else if(e20 == 1) {
                                   assert(e12 == UNULL);
                                   corner12 = { i2, i0, i1, p2, p0, p1, t };
                                   m_containsCorners = true;
                              } else {
                                   assert(e12 == UNULL && e20 == UNULL);
                                   edge2.push_back({ i0, i1, i2, p0, p1, p2, t });
                              }
                         } else {
                              assert(e01 == UNULL);
                              if(e12 == 0) {
                                   if(e20 == 1) {
                                        corner01 = { i1, i2, i0, p1, p2, p0, t };
                                        m_containsCorners = true;
                                   } else {
                                        assert(e20 == UNULL);
                                        edge0.push_back({ i1, i2, i0, p1, p2, p0, t });
                                   }
                              } else if(e12 == 1) {
                                   if(e20 == 2) {
                                        corner12 = { i1, i2, i0, p1, p2, p0, t };
                                        m_containsCorners = true;
                                   } else {
                                        assert(e20 == UNULL);
                                        edge1.push_back({ i1, i2, i0, p1, p2, p0, t });
                                   }
                              } else if(e12 == 2) {
                                   if(e20 == 0) {
                                        corner20 = { i1, i2, i0, p1, p2, p0, t };
                                        m_containsCorners = true;
                                   } else {
                                        assert(e20 == UNULL);
                                        edge2.push_back({ i1, i2, i0, p1, p2, p0, t });
                                   }
                              } else {
                                   assert(e12 == UNULL);
                                   if(e20 == 0) edge0.push_back({ i2, i0, i1, p2, p0, p1, t });
                                   else if(e20 == 1) edge1.push_back({ i2, i0, i1, p2, p0, p1, t });
                                   else if(e20 == 2) edge2.push_back({ i2, i0, i1, p2, p0, p1, t });
                                   else {
                                        assert(e20 == UNULL);
                                        inner.push_back({ i0, i1, i2, p0, p1, p2, t });
                                   }
                              }
                         }
                    }
               }

               m_nEdgeTriangles = edge0.size();
               assert(edge1.size() == m_nEdgeTriangles);
               assert(edge2.size() == m_nEdgeTriangles);

               std::sort(edge0.begin(), edge0.end(), [&] (const hpuint7& a, const hpuint7& b) -> bool { return m_scheme.points[a[0]].y < m_scheme.points[b[0]].y; });
               std::sort(edge1.begin(), edge1.end(), [&] (const hpuint7& a, const hpuint7& b) -> bool { return m_scheme.points[a[0]].z < m_scheme.points[b[0]].z; });
               std::sort(edge2.begin(), edge2.end(), [&] (const hpuint7& a, const hpuint7& b) -> bool { return m_scheme.points[a[0]].x < m_scheme.points[b[0]].x; });

               m_schemeIndices.reserve(scheme.indices.size());
               m_vertexPermutation.reserve(3 * m_nMicroTriangles);
               m_trianglePermutation.reserve(edge0.size() + edge1.size() + edge2.size() + ((m_containsCorners) ? 3 : 0));
               if(m_containsCorners) {
                    m_schemeIndices.insert(m_schemeIndices.end(), corner01.begin(), corner01.begin() + 3);
                    m_schemeIndices.insert(m_schemeIndices.end(), corner12.begin(), corner12.begin() + 3);
                    m_schemeIndices.insert(m_schemeIndices.end(), corner20.begin(), corner20.begin() + 3);
                    m_vertexPermutation.insert(m_vertexPermutation.end(), corner01.begin() + 3, corner01.begin() + 6);
                    m_vertexPermutation.insert(m_vertexPermutation.end(), corner12.begin() + 3, corner12.begin() + 6);
                    m_vertexPermutation.insert(m_vertexPermutation.end(), corner20.begin() + 3, corner20.begin() + 6);
                    m_trianglePermutation.push_back(corner01[6]);
                    m_trianglePermutation.push_back(corner12[6]);
                    m_trianglePermutation.push_back(corner20[6]);
               }
               for(auto& i : edge0) {
                    m_schemeIndices.insert(m_schemeIndices.end(), i.begin(), i.begin() + 3);
                    m_vertexPermutation.insert(m_vertexPermutation.end(), i.begin() + 3, i.begin() + 6);
                    m_trianglePermutation.push_back(i[6]);
               }
               for(auto& i : edge1) {
                    m_schemeIndices.insert(m_schemeIndices.end(), i.begin(), i.begin() + 3);
                    m_vertexPermutation.insert(m_vertexPermutation.end(), i.begin() + 3, i.begin() + 6);
                    m_trianglePermutation.push_back(i[6]);
               }
               for(auto& i : edge2) {
                    m_schemeIndices.insert(m_schemeIndices.end(), i.begin(), i.begin() + 3);
                    m_vertexPermutation.insert(m_vertexPermutation.end(), i.begin() + 3, i.begin() + 6);
                    m_trianglePermutation.push_back(i[6]);
               }
               for(auto& i : inner) {
                    m_schemeIndices.insert(m_schemeIndices.end(), i.begin(), i.begin() + 3);
                    m_vertexPermutation.insert(m_vertexPermutation.end(), i.begin() + 3, i.begin() + 6);
               }
          }

          ProjectiveStructure3D refine(hpreal epsilon = EPSILON) {
               {//inner macro triangle neighbors and transitions
                    auto newNeighbors = make_neighbors(m_schemeIndices);

                    Indices newIndices;
                    newIndices.resize(newNeighbors.size());

                    auto i = m_schemeIndices.begin();
                    auto p = m_vertexPermutation.begin();
                    auto n = newNeighbors.begin();

                    if(m_containsCorners) {
                         handleInnerCornerTriangle(newIndices, i, p, n);
                         handleInnerCornerTriangle(newIndices, i, p, n);
                         handleInnerCornerTriangle(newIndices, i, p, n);
                    }


                    handleInnerEdgeTriangles(newIndices, i, p, n);
                    handleInnerEdgeTriangles(newIndices, i, p, n);
                    handleInnerEdgeTriangles(newIndices, i, p, n);

                    for(auto end = m_schemeIndices.end(); i != end; ++i, ++n, ++p) {
                         auto i0 = *i;
                         auto i1 = *(++i);
                         auto i2 = *(++i);

                         auto n0 = *n;
                         auto n1 = *(++n);
                         auto n2 = *(++n);

                         auto& v0 = m_scheme.points[i0];
                         auto& v1 = m_scheme.points[i1];
                         auto& v2 = m_scheme.points[i2];

                         newIndices[*p] = handleInnerEdge(n0, i0, i1, v0, v1, v2);
                         newIndices[*(++p)] = handleInnerEdge(n1, i1, i2, v1, v2, v0);
                         newIndices[*(++p)] = handleInnerEdge(n2, i2, i0, v2, v0, v1);
                    }

                    for(auto i = 0u; i < m_nMacroTriangles; ++i) m_newIndices.insert(m_newIndices.end(), newIndices.begin(), newIndices.end());
               }

               {//inter macro triangle neighbors and transitions
                    auto newNeighbors = make_neighbors(m_scheme.indices);
                    for(auto i = 0u, end = m_nMicroTriangles * m_nMacroTriangles; i < end; i += m_nMicroTriangles) {
                         for(auto n : newNeighbors) {
                              if(n == UNULL) m_newNeighbors.push_back(UNULL);
                              else m_newNeighbors.push_back(n + i);
                         }
                    }

                    auto b = 0u;
                    auto d = m_newIndices.begin();
                    auto o = m_oldNeighbors.begin();
                    auto r = m_newNeighbors.begin();
                    auto s = 0u;

                    for(auto macroTriangle = 0u; macroTriangle < m_nMacroTriangles; ++macroTriangle, ++b, ++o) {
                         auto b0 = m_oldBorder[b];
                         auto b1 = m_oldBorder[++b];
                         auto b2 = m_oldBorder[++b];

                         auto o0 = *o;
                         auto o1 = *(++o);
                         auto o2 = *(++o);

                         auto f0 = getNeighborOffset(o0, macroTriangle);
                         auto f1 = getNeighborOffset(o1, macroTriangle);
                         auto f2 = getNeighborOffset(o2, macroTriangle);

                         auto i = m_schemeIndices.begin();
                         auto p = m_vertexPermutation.begin();
                         auto t = m_trianglePermutation.begin();

                         if(m_containsCorners) {
                              m_newBorder[s + p[0]] = b0;
                              m_newBorder[s + p[7]] = b0;
                              m_newBorder[s + p[1]] = b1;
                              m_newBorder[s + p[3]] = b1;
                              m_newBorder[s + p[4]] = b2;
                              m_newBorder[s + p[6]] = b2;
                              handleInterCornerTriangle<0, 1>(macroTriangle, i, p, t, d, r, o0, o1, f0, f1);
                              handleInterCornerTriangle<1, 2>(macroTriangle, i, p, t, d, r, o1, o2, f1, f2);
                              handleInterCornerTriangle<2, 0>(macroTriangle, i, p, t, d, r, o2, o0, f2, f0);
                         }

                         handleInterEdgeTriangles<0>(macroTriangle, i, p, t, d, r, s, o0, f0, b0);
                         handleInterEdgeTriangles<1>(macroTriangle, i, p, t, d, r, s, o1, f1, b1);
                         handleInterEdgeTriangles<2>(macroTriangle, i, p, t, d, r, s, o2, f2, b2);

                         d += 3 * m_nMicroTriangles;
                         r += 3 * m_nMicroTriangles;
                         s += 3 * m_nMicroTriangles;
                    }
               }

               return ProjectiveStructure3D(std::move(m_newTransitions), std::move(m_newIndices), std::move(m_newNeighbors), std::move(m_newBorder));
          }

     private:
          bool m_containsCorners;//true if the refinement has triangles at the corners
          Indices m_schemeIndices;
          const TriangleRefinementScheme& m_scheme;
          boost::dynamic_bitset<> m_newBorder;
          Indices m_newIndices;
          Neighbors m_newNeighbors;
          Transitions m_newTransitions;
          const boost::dynamic_bitset<>& m_oldBorder;
          const Indices& m_oldIndices;
          const Neighbors& m_oldNeighbors;
          const Transitions& m_oldTransitions;
          hpuint m_nEdgeTriangles;
          const hpuint m_nMacroTriangles;
          const hpuint m_nMicroTriangles;
          Indices m_trianglePermutation;
          Indices m_vertexPermutation;

          static hpuint getEdgeType(const Point3D& u, const Point3D& v, hpreal epsilon) {
               if(std::abs(u.x) < epsilon && std::abs(v.x) < epsilon) return 1;
               else if(std::abs(u.y) < epsilon && std::abs(v.y) < epsilon) return 2;
               else if(std::abs(u.z) < epsilon && std::abs(v.z) < epsilon) return 0;
               else return UNULL;
          }
     
          static inline void swapRows(hpmat3x3& matrix, hpuint r0, hpuint r1) {
               assert(r0 < 3 && r1 < 3);
     
               auto t0 = matrix[0][r0];
               auto t1 = matrix[1][r0];
               auto t2 = matrix[2][r0];
     
               matrix[0][r0] = matrix[0][r1];
               matrix[1][r0] = matrix[1][r1];
               matrix[2][r0] = matrix[2][r1];
     
               matrix[0][r1] = t0;
               matrix[1][r1] = t1;
               matrix[2][r1] = t2;
          }

          hpuint getNeighborOffset(hpuint triangle, hpuint neighbor) {
               if(triangle == UNULL) return -1;
               auto n = m_oldNeighbors.cbegin() + 3 * triangle;
               if (*n == neighbor) return 0;
               else if (*(++n) == neighbor) return 1;
               else {
                    assert(*(++n) == neighbor);
                    return 2;
               }
          }

          template<class Iterator>
          void handleInnerCornerTriangle(Indices& newIndices, Iterator& i, Iterator& p, Iterator& n) {
               auto i0 = *i;
               auto i1 = *(++i);
               auto i2 = *(++i);
               ++i;

               ++n;
               const hpuint n2 = *(++n);
               ++n;

               auto& v0 = m_scheme.points[i0];
               auto& v1 = m_scheme.points[i1];
               auto& v2 = m_scheme.points[i2];

               newIndices[*p] = UNULL;
               newIndices[*(++p)] = UNULL;
               newIndices[*(++p)] = handleInnerEdge(n2, i2, i0, v2, v0, v1);
               ++p;
          }

          hpuint handleInnerEdge(hpuint neighbor, hpuint i0, hpuint i1, const Point3D& v0, const Point3D& v1, const Point3D& v2) {
               auto j = m_schemeIndices.cbegin() + 3 * neighbor;
               auto j0 = *j;
               auto j1 = *(++j);
               auto j2 = *(++j);
               hpuint i3;
               if(j0 == i0) {
                    if(j1 == i1) i3 = j2;
                    else if(j2 == i1) i3 = j1;
                    else assert(false);
               } else if(j0 == i1) {
                    if(j1 == i0) i3 = j2;
                    else if(j2 == i0) i3 = j1;
                    else assert(false);
               } else if((j1 == i0 && j2 == i1) || (j1 == i1 && j2 == i0)) i3 = j0;
               else assert(false);
               auto& v3 = m_scheme.points[i3];
               hpmat3x3 beta(v1, v0, v2);
               hpmat3x3 gamma(v1, v3, v0);
               auto betaInverse = glm::inverse(beta);
               auto result = betaInverse * gamma;
               return insert({ result[1][0], result[1][1], result[1][2] }); 
          }

          template<class Iterator>
          void handleInnerEdgeTriangles(Indices& newIndices, Iterator& i, Iterator& p, Iterator& n) {
               for(auto end = i + 3 * m_nEdgeTriangles; i != end; ++i, ++n, ++p) {
                    auto i0 = *i;
                    auto i1 = *(++i);
                    auto i2 = *(++i);

                    auto n1 = *(++n);
                    auto n2 = *(++n);

                    auto& v0 = m_scheme.points[i0];
                    auto& v1 = m_scheme.points[i1];
                    auto& v2 = m_scheme.points[i2];

                    newIndices[*p] = UNULL;
                    newIndices[*(++p)] = handleInnerEdge(n1, i1, i2, v1, v2, v0);
                    newIndices[*(++p)] = handleInnerEdge(n2, i2, i0, v2, v0, v1);
               }
          }

          template<hpuint edge0, hpuint edge1, class Iterator>
          void handleInterCornerTriangle(const hpuint macroTriangle, Iterator& i, Iterator& p, Iterator& t, Iterator d, Iterator r, hpuint o0, hpuint o1, hpuint f0, hpuint f1) {
               auto& v0 = m_scheme.points[*i];
               auto& v1 = m_scheme.points[*(++i)];
               auto& v2 = m_scheme.points[*(++i)];
               ++i;

               if(o0 != UNULL) {
                    hpuint neighbor;
                    if(f0 == 0) neighbor = 2;
                    else if(f0 == 1) neighbor = 0;
                    else neighbor = 1;
                    auto j = m_schemeIndices.begin() + 3 * neighbor;
                    auto& w0 = m_scheme.points[*j];
                    auto& w1 = m_scheme.points[*(++j)];
                    auto& w2 = m_scheme.points[*(++j)];
                    d[*p] = handleInterEdge<edge0>(macroTriangle, v1, v0, v2, w1, w0, w2, f0);
                    r[*p] = t[neighbor] + o0 * m_nMicroTriangles;
               }
               ++p;

               if(o1 != UNULL) {
                    auto neighbor = f1;
                    auto j = m_schemeIndices.begin() + 3 * neighbor;
                    auto& w0 = m_scheme.points[*j];
                    auto& w1 = m_scheme.points[*(++j)];
                    auto& w2 = m_scheme.points[*(++j)];
                    d[*p] = handleInterEdge<edge1>(macroTriangle, v2, v1, v0, w0, w2, w1, f1);
                    r[*p] = t[neighbor] + o1 * m_nMicroTriangles;
               }
               p += 2;
          }

          template<hpuint edge> 
          hpuint handleInterEdge(const hpuint macroTriangle, const Point3D& beta1, const Point3D& beta2, const Point3D& beta3, const Point3D& gamma1, const Point3D& gamma2, const Point3D& gamma3, const hpuint f) {
               hpmat3x3 beta(beta1, beta2, beta3);
               hpmat3x3 gamma(gamma1, gamma2, gamma3);
               if(edge == 0) swapRows(beta, 0, 1);
               else if(edge == 1) swapRows(beta, 0, 2);
               else swapRows(beta, 1, 2);
               if(f == 0) swapRows(gamma, 1, 2);
               else if(f == 1) swapRows(gamma, 0, 1);
               else swapRows(gamma, 0, 2);
               auto betaInverse = glm::inverse(beta);
               auto& macroTransition = m_oldTransitions[m_oldIndices[3 * macroTriangle + edge]];
               hpmat3x3 macro({ 1, 0, 0 }, macroTransition, { 0, 1, 0 });
               auto temp = betaInverse * macro * gamma;
               return insert({ temp[1][0], temp[1][1], temp[1][2] }); 
          }
     
          template<hpuint edge, class Iterator>
          void handleInterEdgeTriangles(const hpuint macroTriangle, Iterator& i, Iterator& p, Iterator& t, Iterator d, Iterator r, hpuint s, hpuint o, hpuint f, bool b) {
               if(o == UNULL) {
                    p += 3 * m_nEdgeTriangles;
                    i += 3 * m_nEdgeTriangles;
               } else {
                    auto index = 0u;
                    for(auto end = i + 3 * m_nEdgeTriangles; i != end; ++i, p += 3) {
                         auto& v0 = m_scheme.points[*i];
                         auto& v1 = m_scheme.points[*(++i)];
                         auto& v2 = m_scheme.points[*(++i)];

                         auto neighbor = (f + 1) * m_nEdgeTriangles - ++index;
                         if(m_containsCorners) neighbor += 3;

                         auto j = m_schemeIndices.begin() + 3 * neighbor;
                         auto& w0 = m_scheme.points[*j];
                         auto& w1 = m_scheme.points[*(++j)];
                         auto& w2 = m_scheme.points[*(++j)];

                         d[*p] = handleInterEdge<edge>(macroTriangle, v1, v0, v2, w0, w2, w1, f);
                         r[*p] = t[neighbor] + o * m_nMicroTriangles;
                         m_newBorder[s + *p] = b;
                    }
               }
          }

          /**
           * The transition is inserted into the vector if it is not already there.
           *
           * @return Returns the index of the transition in the vector.
           */
          hpuint insert(Transition transition) {
               auto result = std::find_if(m_newTransitions.begin(), m_newTransitions.end(), [&] (const Transition& other) -> bool { return (std::abs(transition.x - other.x) < EPSILON) && (std::abs(transition.y - other.y) < EPSILON) && (std::abs(transition.z - other.z) < EPSILON); });
               if(result != m_newTransitions.end()) return std::distance(m_newTransitions.begin(), result);
               m_newTransitions.push_back(std::move(transition));
               return m_newTransitions.size() - 1;
          }

     };//Refiner

public:
     static const ProjectiveStructure DOUBLE_TORUS_PROJECTIVE_STRUCTURE;
     template<View view, Mode... modes>
     using const_iterator = Iterator<0, view, modes...>;
     using default_const_iterator = const_iterator<View::TETRAHEDRA, Mode::TETRAHEDRA, Mode::TRANSITIONS>;

     //NOTE: Transitions are calculated based on the given points and their neighborhood relationships.
     ProjectiveStructure toProjectiveStructure(const std::vector<Point3D>& points, const std::vector<hpuint>& indices, Neighbors neighbors, Border border) {
          auto transitions = ProjectiveStructureUtils::toTransitions(points, indices, neighbors);
          return {std::move(transitions.first), std::move(transitions.second), std::move(neighbors), std::move(border)};
     }

     /** 
      * Constructs a projective structure.
      *
      * @param[in] transitions  Set of transitions, i.e. middle columns of the transitions from one tetrahedron to another.
      * @param[in] indices      Three indices per tetrahedron specifying which transition to use when moving to the neighbor.
      * @param[in] neighbors    Three indices per tetrahedron specifying its three neighbors.
      * @param[in] border       Two indices per border edge specifying the two incident tetrahedra.
      */
     ProjectiveStructure(Transitions transitions, Indices indices, Neighbors neighbors, Border border) 
          : ProjectiveStructureBase<Space3D>(std::move(transitions), std::move(indices), std::move(neighbors), std::move(border)), m_nVertices(0), m_vertices(vertices()) {}

     default_const_iterator cbegin() const { return default_const_iterator(*this, 0); };

     template<View view, Mode... modes>
     Iterator<0, view, modes...> cbegin() const { return Iterator<0, view, modes...>(*this, 0); };

     default_const_iterator cend() const { return default_const_iterator(*this, getNumberOfSimplices()); };

     template<View view, Mode... modes>
     Iterator<0, view, modes...> cend() const { return Iterator<0, view, modes...>(*this, (view == View::VERTICES) ? m_nVertices : getNumberOfSimplices()); };

     using ProjectiveStructureBase<Space3D>::getNeighbors;

     std::tuple<hpuint, hpuint, hpuint> getNeighbors(hpuint simplex) const {
          auto n = m_neighbors.begin() + simplex * 3;
          auto n0 = *n;
          auto n1 = *(++n);
          auto n2 = *(++n);
          return std::make_tuple(n0, n1, n2);
     }

     hpuint getNumberOfVertices() const { return m_nVertices; }

     const Vertices& getVertices() const { return m_vertices; } //TODO: iterator?

     std::tuple<hpuint, hpuint, hpuint> getVertices(hpuint simplex) const {
          auto v = m_vertices.begin() + simplex * 3;
          auto v0 = *v;
          auto v1 = *(++v);
          auto v2 = *(++v);
          return std::make_tuple(v0, v1, v2);
     }

     /**
      * The vertices and the neighbors are implicitly ordered according to the following scheme:
      *
      *                           (0,0,1)
      *                         /         \
      *            neighbor 2 /             \ neighbor 1
      *                     /                 \
      *                   /                     \
      *              (1,0,0) - - - - - - - - - (0,1,0)
      *                         neighbor 0
      *
      * How the scheme is applied to a particular tetrahedron depends on the order of the vertices of the macro tetrahedron.
      *
      * @param[in]  mask  Mask specifying the refinement scheme.  The only requirements on the refinement scheme is that all edges are refined in the exactly same way and that the indices of the triangles are arranged in counter-clockwise order.
      * @param[out]       Returns the constructed projective structure in which every macro-simplex of this projective structure is subdivided into micro-simplices as specified by the refinement scheme.
      */
     ProjectiveStructure3D refine(const TriangleRefinementScheme& scheme, hpreal epsilon = EPSILON) const { 
          Refiner refiner(*this, scheme, epsilon);
          return std::move(refiner.refine(epsilon));
     }

     /**
      * Constructs a triangle mesh represented by this projective structure where the indexed tetrahedron has the corners given by p0, p1, and p2 in counter-clockwise order respecting the implicit ordering in the neighbors array.
      *
      * @param[in] tetrahedron     Index of tetrahedron whose corners are given.
      * @param[in] p0              First corner of tetrahedron.
      * @param[in] p1              Second corner of tetrahedron.
      * @param[in] p2              Third corner of tetrahedron.
      */
     template<class Vertex = VertexP3, class VertexFactory = VertexFactory<Vertex> >
     TriangleMesh<Vertex> toTriangleMesh(hpuint tetrahedron, const Point3D& p0, const Point3D& p1, const Point3D& p2, VertexFactory factory = VertexFactory()) const {
          std::vector<Vertex> vertices;
          std::vector<hpuint> indices;
          indices.resize(m_indices.size());

          vertices.push_back(factory(p0));
          vertices.push_back(factory(p1));
          vertices.push_back(factory(p2));

          hpuint o = 3 * tetrahedron;

          auto i = indices.begin() + o;
          *i = 0;
          *(++i) = 1;
          *(++i) = 2;

          std::stack<hpuint> tetrahedra;
          auto n = m_neighbors.cbegin() + o;
          if(!m_border[o]) tetrahedra.push(*n);
          ++n;
          if(!m_border[o + 1]) tetrahedra.push(*n);
          ++n;
          if(!m_border[o + 2]) tetrahedra.push(*n);

          hpuint nTriangles = m_neighbors.size() / 3;
          boost::dynamic_bitset<> done(nTriangles);
          done[tetrahedron] = true;

          while(!tetrahedra.empty()) {
               hpuint tetrahedron = tetrahedra.top();
               tetrahedra.pop();
               if(done[tetrahedron]) continue;

               hpuint o = 3 * tetrahedron;

               auto n = m_neighbors.cbegin() + o;
               hpuint n0 = *n;
               hpuint n1 = *(++n);
               hpuint n2 = *(++n);

               bool b0 = m_border[o];
               bool b1 = m_border[o + 1];
               bool b2 = m_border[o + 2];

               hpuint v0 = -1;
               hpuint v1 = -1;
               hpuint v2 = -1;

               auto get = [&](hpuint neighbor) -> std::pair<hpuint, hpuint> {
                    hpuint o = 3 * neighbor;
                    auto i = indices.cbegin() + o;
                    hpuint i0 = *i;
                    hpuint i1 = *(++i);
                    hpuint i2 = *(++i);
                    auto n = m_neighbors.cbegin() + o;
                    if(*n == tetrahedron) return std::make_pair(i0, i1);
                    else if(*(++n) == tetrahedron) return std::make_pair(i1, i2);
                    else {
                         assert(*(++n) == tetrahedron);
                         return std::make_pair(i2, i0);
                    }
               };//get

               auto find = [&](hpuint neighborCW, hpuint neighborCCW, hpuint neighborOpposite, bool borderCW, bool borderCCW) -> hpuint {
                    if(!borderCCW) {
                         hpuint previous = tetrahedron;
                         hpuint current = neighborCCW;
                         while(!(current == tetrahedron && previous == neighborCW)) {//visiting counterclockwise
                              hpuint o = 3 * current;
                              auto n = m_neighbors.cbegin() + o;
                              hpuint n0 = *n;
                              hpuint n1 = *(++n);
                              hpuint n2 = *(++n);
                              if(n0 == previous) {
                                   if(m_border[o + 2]) goto l_visitClockwise;
                              } else if(n1 == previous) {
                                   if(m_border[o]) goto l_visitClockwise;
                              } else {
                                   assert(n2 == previous);
                                   if(m_border[o + 1]) goto l_visitClockwise;
                              }
                              if(done[current]) {
                                   if(n0 == previous) return indices[o];
                                   else if(n1 == previous) return indices[o + 1];
                                   else return indices[o + 2];
                              } else {
                                   if(n0 == previous) {
                                        previous = current;
                                        current = n2;
                                   } else if(n1 == previous) {
                                        previous = current;
                                        current = n0;
                                   } else {
                                        previous = current;
                                        current = n1;
                                   }
                              }
                         }
                         goto l_generateVertex;
                    }
                    l_visitClockwise:
                    if(!borderCW) {
                         hpuint previous = tetrahedron;
                         hpuint current = neighborCW;
                         while(!(current == tetrahedron && previous == neighborCCW)) {//visiting clockwise
                              hpuint o = 3 * current;
                              auto n = m_neighbors.cbegin() + o;
                              hpuint n0 = *n;
                              hpuint n1 = *(++n);
                              hpuint n2 = *(++n);
                              if(n0 == previous) {
                                   if(m_border[o + 1]) goto l_generateVertex;
                              } else if(n1 == previous) {
                                   if(m_border[o + 2]) goto l_generateVertex;
                              } else {
                                   assert(n2 == previous);
                                   if(m_border[o]) goto l_generateVertex;
                              }
                              if(done[current]) {
                                   if(n0 == previous) return indices[o + 1];
                                   else if(n1 == previous) return indices[o + 2];
                                   else return indices[o];
                              } else {
                                   if(n0 == previous) {
                                        previous = current;
                                        current = n1;
                                   } else if(n1 == previous) {
                                        previous = current;
                                        current = n2;
                                   } else {
                                        previous = current;
                                        current = n0;
                                   }
                              }
                         }
                    }
                    l_generateVertex:
                    hpuint v = vertices.size();
                    hpuint o = 3 * neighborOpposite;
                    auto i = indices.cbegin() + o;
                    auto& p0 = vertices[*i].position;
                    auto& p1 = vertices[*(++i)].position;
                    auto& p2 = vertices[*(++i)].position;
                    auto n = m_neighbors.cbegin() + o;
                    if(*n == tetrahedron) {
                         const Transition& t = m_transitions[m_indices[o]];
                         vertices.push_back(factory(t.x * p1 + t.y * p0 + t.z * p2));
                    } else if(*(++n) == tetrahedron) {
                         const Transition& t = m_transitions[m_indices[o + 1]];
                         vertices.push_back(factory(t.x * p2 + t.y * p1 + t.z * p0));
                    } else {
                         assert(*(++n) == tetrahedron);
                         const Transition& t = m_transitions[m_indices[o + 2]];
                         vertices.push_back(factory(t.x * p0 + t.y * p2 + t.z * p1));
                    }
                    return v;
               };//find

               hpuint v;
               if(done[n0] && !b0) {
                    std::tie(v1, v0) = get(n0);
                    if(done[n1] && !b1) std::tie(v2, v) = get(n1);
                    else if(done[n2] && !b2) std::tie(v, v2) = get(n2);
                    else v2 = find(n2, n1, n0, b2, b1);
               } else if(done[n2] && !b2) {
                    std::tie(v0, v2) = get(n2);
                    if(done[n1] && !b1) std::tie(v, v1) = get(n1);
                    else v1 = find(n1, n0, n2, b1, b0);
               } else {
                    std::tie(v2, v1) = get(n1);
                    v0 = find(n0, n2, n1, b0, b2);
               }

               done[tetrahedron] = true;

               auto i = indices.begin() + o;
               *i = v0;
               *(++i) = v1;
               *(++i) = v2;

               if(!done[n0] && !b0) tetrahedra.push(n0);
               if(!done[n1] && !b1) tetrahedra.push(n1);
               if(!done[n2] && !b2) tetrahedra.push(n2);
          }

          return make_triangle_mesh(std::move(vertices), std::move(indices));
     }

private:
     hpuint m_nVertices;
     Vertices m_vertices;

     ProjectiveStructure(Transitions transitions, Indices indices, Neighbors neighbors, boost::dynamic_bitset<> border)
          : ProjectiveStructureBase<Space3D>(std::move(transitions), std::move(indices), std::move(neighbors), std::move(border)), m_nVertices(0), m_vertices(vertices()) {}

     Vertices vertices() {
          Vertices vertices;
          vertices.resize(m_neighbors.size());
          for(auto i = Iterator<0, View::VERTICES, Mode::TETRAHEDRA>(*this); !i.end(); ++i) {
               auto& fan = *i;
               auto f = fan.cbegin();
               auto first = *f;
               auto e = fan.cend();
               auto current = first;
               ++f;
               while(f != e) {
                    auto next = *f;
                    auto n = m_neighbors.cbegin() + 3 * current;
                    if(next == *n) vertices[3 * current + 1] = m_nVertices;
                    else if(next == *(++n)) vertices[3 * current + 2] = m_nVertices;
                    else {
                         assert(next == *(++n));
                         vertices[3 * current] = m_nVertices;
                    }
                    ++f;
                    current = next;
               }
               auto next = first;
               auto n = m_neighbors.cbegin() + 3 * current;
               if(next == *n) vertices[3 * current + 1] = m_nVertices; 
               else if(next == *(++n)) vertices[3 * current + 2] = m_nVertices;
               else {
                    assert(next == *(++n));
                    vertices[3 * current] = m_nVertices;
               }
               ++m_nVertices;
          }
          return vertices;
     }

};//ProjectiveStructure

}//namespace happah

