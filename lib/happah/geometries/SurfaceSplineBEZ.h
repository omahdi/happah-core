// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/dynamic_bitset.hpp>
#include <type_traits>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/Surface.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/utils/DeindexedArray.h"
#include "happah/utils/SurfaceSplineSubdividerBEZ.h"
#include "happah/utils/SurfaceSplineUtilsBEZ.h"
#include "happah/utils/SurfaceUtilsBEZ.h"
#include "happah/utils/VertexFactory.h"

namespace happah {

template<class Space, hpuint t_degree>
class SurfaceSplineBEZ : public Surface<Space> {
     using Point = typename Space::POINT;

public:
     using ControlPoints = std::vector<Point>;
     using Indices = std::vector<hpuint>;
     using Mode = SurfaceSplineUtilsBEZ::Mode;
     using ParameterPoints = std::vector<Point2D>;

private:
     template<int t_dummy = 0, Mode... t_modes>
     class Iterator {
     public:
          using difference_type = hpuint;
          using value_type = std::tuple<typename Iterator<t_dummy, t_modes>::value_type...>;

          Iterator(const SurfaceSplineBEZ& surface, hpuint patch)
               : m_i(Iterator<t_dummy, t_modes>(surface, patch)...) {}

          difference_type operator-(const Iterator& iterator) { return std::get<0>(m_i) - std::get<0>(iterator.m_i); }

          Iterator& operator++() {
               preincrement();
               return *this;
          }

          Iterator& operator--() {
               predecrement();
               return *this;
          }

          Iterator operator++(int) const {
               Iterator iterator(*this);
               return ++iterator;
          }

          Iterator operator--(int) const {
               Iterator iterator(*this);
               return --iterator;
          }

          bool operator!=(const Iterator& iterator) const { return std::get<0>(iterator.m_i) != std::get<0>(m_i); }

          value_type operator*() const { return getValue(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

     private:
          std::tuple<Iterator<t_dummy, t_modes>...> m_i;

          template<unsigned long... Is>
          value_type getValue(std::index_sequence<Is...>) const { return std::make_tuple(*(std::get<Is>(m_i))...); }

          BUILD_TUPLE_HANDLER_METHODS(predecrement, doPredecrement)

          void predecrement() { predecrement(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<Mode t_mode>
          void doPredecrement(Iterator<t_dummy, t_mode>& i) { --i; }

          BUILD_TUPLE_HANDLER_METHODS(preincrement, doPreincrement)

          void preincrement() { preincrement(std::make_integer_sequence<std::size_t, std::tuple_size<decltype(m_i)>::value>()); }

          template<Mode t_mode>
          void doPreincrement(Iterator<t_dummy, t_mode>& i) { ++i; }

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, Mode::CONTROL_POINTS> {
          using ProxyIterator = typename DeindexedArray<ControlPoints, Indices>::const_iterator;

     public:
          using difference_type = typename ProxyIterator::difference_type;
          using value_type = ProxyIterator;

          Iterator(const SurfaceSplineBEZ& surface, hpuint patch) 
               : m_i(i(surface, patch)) {}

          difference_type operator-(const Iterator& iterator) const { return (m_i - iterator.m_i) / m_nControlPoints; }

          Iterator& operator++() { 
               m_i += m_nControlPoints;
               return *this;
          }

          Iterator& operator--() { 
               m_i -= m_nControlPoints;
               return *this;
          }

          Iterator operator++(int) const { 
               Iterator iterator(*this);
               return ++iterator;
          }

          Iterator operator--(int) const {
               Iterator iterator(*this);
               return --iterator;
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_i != m_i; }

          value_type operator*() const { return m_i; }

     private:
          ProxyIterator m_i;
          static constexpr hpuint m_nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;

          static ProxyIterator i(const SurfaceSplineBEZ& surface, hpuint patch) { return deindex(surface.m_controlPoints, surface.m_controlPointIndices).cbegin() + m_nControlPoints * patch; }

     };//Iterator

     template<int t_dummy>
     class Iterator<t_dummy, Mode::NEIGHBORS> {
     public:
          using difference_type = hpuint;
          using value_type = std::tuple<hpuint, hpuint, hpuint>;

          Iterator(const SurfaceSplineBEZ& surface, hpuint patch) 
               : m_patch(patch), m_surface(surface) {}

          difference_type operator-(const Iterator& iterator) const { return m_patch - iterator.m_patch; }

          Iterator& operator++() {
               ++m_patch;
               return *this;
          }

          Iterator& operator--() {
               --m_patch;
               return *this;
          }

          Iterator operator++(int) const { 
               Iterator iterator(*this);
               return ++iterator;
          }

          Iterator operator--(int) const {
               Iterator iterator(*this);
               return --iterator;
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_patch != m_patch; }

          value_type operator*() const { return m_surface.getPatchNeighbors(m_patch); }

     private:
          hpuint m_patch;
          const SurfaceSplineBEZ& m_surface;

     };//Iterator

public:
     //NOTE: There is no automatic way of figuring out the neighborhood of patches given only the control points.
     SurfaceSplineBEZ(ControlPoints controlPoints, Indices controlPointIndices)
          : m_controlPointIndices{std::move(controlPointIndices)}, m_controlPoints{std::move(controlPoints)} {}

     SurfaceSplineBEZ(ControlPoints controlPoints, Indices controlPointIndices, ParameterPoints parameterPoints, Indices parameterPointIndices)
          : m_controlPointIndices{std::move(controlPointIndices)}, m_controlPoints{std::move(controlPoints)}, m_parameterPointIndices{std::move(parameterPointIndices)}, m_parameterPoints{std::move(parameterPoints)} {}

     template<Mode... t_modes>
     Iterator<0, t_modes...> cbegin() const { return Iterator<0, t_modes...>(*this, 0); };

     template<Mode... t_modes>
     Iterator<0, t_modes...> cend() const { return Iterator<0, t_modes...>(*this, getNumberOfPatches()); };

     const ControlPoints& getControlPoints() const { return m_controlPoints; }

     template<class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, bool t_includeParameterPoints = is_relative_vertex<Vertex>::value>
     TriangleMesh<Vertex> getControlPolygon(VertexFactory&& factory = VertexFactory()) const { return toTriangleMesh<Vertex, VertexFactory, t_includeParameterPoints>(0, std::forward<VertexFactory>(factory)); }

     std::tuple<hpuint, hpuint, hpuint> getPatchNeighbors(hpuint patch) const { return TriangleMeshUtils::getNeighbors(m_parameterPointIndices, patch); }

     const ParameterPoints& getParameterPoints() const { return m_parameterPoints; }

     hpuint getNumberOfPatches() const { return m_controlPointIndices.size() / SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value; }

     template<class Visitor>
     void sample(hpuint nSamples, Visitor visit) const {
          auto matrix = SurfaceUtilsBEZ::getEvaluationMatrix<t_degree>(nSamples);
          //TODO; skip multiple computations on common edges and eliminate common points in array
          //TODO: t_degree == 1,2,4, general
          //TODO: optimize calculations below; (u,v,w) do not have to be recalculated each time
          switch(t_degree) {
          case 1:

               break;
          case 2:
 
               break;
          case 3: {
               auto c = m_controlPointIndices.cbegin();
               for(auto i = m_parameterPointIndices.cbegin(), end = m_parameterPointIndices.cend(); i != end; ++i) {
                    const Point2D& p0 = m_parameterPoints[*i];
                    const Point2D& p1 = m_parameterPoints[*(++i)];
                    const Point2D& p2 = m_parameterPoints[*(++i)];

                    const Point& c0 = m_controlPoints[*c];
                    const Point& c1 = m_controlPoints[*(++c)];
                    const Point& c2 = m_controlPoints[*(++c)];
                    const Point& c3 = m_controlPoints[*(++c)];
                    const Point& c4 = m_controlPoints[*(++c)];
                    const Point& c5 = m_controlPoints[*(++c)];
                    const Point& c6 = m_controlPoints[*(++c)];
                    const Point& c7 = m_controlPoints[*(++c)];
                    const Point& c8 = m_controlPoints[*(++c)];
                    const Point& c9 = m_controlPoints[*(++c)];
                    ++c;

                    auto m = matrix.cbegin();
                    SurfaceUtilsBEZ::sample(nSamples, [&] (hpreal u, hpreal v, hpreal w) {
                         Point temp = *m * c0;
                         temp += *(++m) * c1;
                         temp += *(++m) * c2;
                         temp += *(++m) * c3;
                         temp += *(++m) * c4;
                         temp += *(++m) * c5;
                         temp += *(++m) * c6;
                         temp += *(++m) * c7;
                         temp += *(++m) * c8;
                         temp += *(++m) * c9;
                         ++m;
                         visit(u * p0 + v * p1 + w * p2, temp);
                    });
               }
               break;
          }
          case 4: {//TODO: implementation in hez is exactly the same; refactor to one location
               auto c = m_controlPointIndices.cbegin();
               for(auto i = m_parameterPointIndices.cbegin(), end = m_parameterPointIndices.cend(); i != end; ++i) {
                    auto& p0 = m_parameterPoints[*i];
                    auto& p1 = m_parameterPoints[*(++i)];
                    auto& p2 = m_parameterPoints[*(++i)];

                    auto& c0 = m_controlPoints[*c];
                    auto& c1 = m_controlPoints[*(++c)];
                    auto& c2 = m_controlPoints[*(++c)];
                    auto& c3 = m_controlPoints[*(++c)];
                    auto& c4 = m_controlPoints[*(++c)];
                    auto& c5 = m_controlPoints[*(++c)];
                    auto& c6 = m_controlPoints[*(++c)];
                    auto& c7 = m_controlPoints[*(++c)];
                    auto& c8 = m_controlPoints[*(++c)];
                    auto& c9 = m_controlPoints[*(++c)];
                    auto& c10 = m_controlPoints[*(++c)];
                    auto& c11 = m_controlPoints[*(++c)];
                    auto& c12 = m_controlPoints[*(++c)];
                    auto& c13 = m_controlPoints[*(++c)];
                    auto& c14 = m_controlPoints[*(++c)];
                    ++c;

                    auto m = matrix.cbegin();
                    SurfaceUtilsBEZ::sample(nSamples, [&] (hpreal u, hpreal v, hpreal w) {
                         auto sample = *m * c0;
                         sample += *(++m) * c1;
                         sample += *(++m) * c2;
                         sample += *(++m) * c3;
                         sample += *(++m) * c4;
                         sample += *(++m) * c5;
                         sample += *(++m) * c6;
                         sample += *(++m) * c7;
                         sample += *(++m) * c8;
                         sample += *(++m) * c9;
                         sample += *(++m) * c10;
                         sample += *(++m) * c11;
                         sample += *(++m) * c12;
                         sample += *(++m) * c13;
                         sample += *(++m) * c14;
                         ++m;
                         visit(u * p0 + v * p1 + w * p2, sample);
                    });
               }
               break;
          }
          default:

               break;
          }
     }

     template<bool t_includeParameterPoints = false>
     SurfaceSplineBEZ subdivide(hpuint nSubdivisions) const {
          assert(!t_includeParameterPoints || (t_includeParameterPoints && m_parameterPoints.size() > 0));
          using Subdivider = SurfaceSplineSubdividerBEZ<Space, t_degree>;
          Subdivider subdivider(cbegin<Mode::CONTROL_POINTS>(), getNumberOfPatches());
          auto subdivided = subdivider.subdivide(nSubdivisions);
          if(t_includeParameterPoints) {
               auto temp = Subdivider::getParameterPoints(m_parameterPoints, m_parameterPointIndices, nSubdivisions);
               return SurfaceSplineBEZ(subdivided.first, subdivided.second, temp.first, temp.second);
          } else return SurfaceSplineBEZ(subdivided.first, subdivided.second);
     }

     template<class Vertex = VertexP<Space>, class VertexFactory = happah::VertexFactory<Vertex>, bool t_includeParameterPoints = is_relative_vertex<Vertex>::value, typename = typename std::enable_if<(t_degree > 0)>::type>
     TriangleMesh<Vertex> toTriangleMesh(hpuint nSubdivisions = 4, VertexFactory&& factory = VertexFactory()) const {
          if(nSubdivisions > 0) {
               auto subdivided = subdivide<t_includeParameterPoints>(nSubdivisions);
               return subdivided.getControlPolygon<Vertex, VertexFactory, t_includeParameterPoints>(std::forward<VertexFactory>(factory));
          }
          TriangleMeshBuilder<Vertex, t_includeParameterPoints> builder(*this);
          return builder.build(std::forward<VertexFactory>(factory));
     }

private:
     Indices m_controlPointIndices;
     ControlPoints m_controlPoints;
     static constexpr hpuint m_nControlPoints = SurfaceUtilsBEZ::get_number_of_control_points<t_degree>::value;
     Indices m_parameterPointIndices;
     ParameterPoints m_parameterPoints;

     template<class Vertex, bool t_includeParameterPoints>
     struct TriangleMeshBuilder;

     template<class Vertex>
     struct TriangleMeshBuilder<Vertex, false> {
          
          TriangleMeshBuilder(const SurfaceSplineBEZ& surface) : m_surface(surface) {}

          template<class VertexFactory>
          TriangleMesh<Vertex> build(VertexFactory&& factory) const {
               static_assert(std::is_base_of<Vertex, decltype(factory(Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");

               std::vector<Vertex> vertices;
               vertices.reserve(m_surface.m_controlPoints.size());

               for(const Point& controlPoint : m_surface.m_controlPoints)
                    vertices.push_back(factory(controlPoint));

               std::vector<hpuint> indices = SurfaceSplineUtilsBEZ::template buildTriangleMeshIndices<t_degree>(m_surface.m_controlPointIndices);

               return make_triangle_mesh(std::move(vertices), std::move(indices));
          }

     private:
          const SurfaceSplineBEZ& m_surface;

     };//TriangleMeshBuilder

     template<class Vertex>
     struct TriangleMeshBuilder<Vertex, true> {

          TriangleMeshBuilder(const SurfaceSplineBEZ& surface) : m_surface(surface) {}

          template<class VertexFactory>
          TriangleMesh<Vertex> build(VertexFactory&& factory) const {
               static_assert(std::is_base_of<Vertex, decltype(factory(Point2D(0.0), Point(0.0)))>::value, "The vertex generated by the factory must be a subclass of the vertex with which the triangle mesh is parameterized.");
               
               assert(m_surface.m_parameterPointIndices.size() > 0);
               assert(m_surface.m_parameterPoints.size() > 0);

               //TODO: Similar control points can be shared between different edges and/or patches.  Therefore, when including the parameter points, we cannot reuse the original control point indices but generate new ones.  Here, we do not optimize and assume each resulting vertex is different.  However, we can optimize slightly; any control points shared on an edge between two adjacent patches can be shared in the triangle mesh.  In fact, for usability in the graphical interface, these should be the same if they are the same originally.  Thus, we need a data structure that tells us which edges are shared and update the generated control point indices accordingly.

               std::vector<Vertex> vertices;
               vertices.reserve(m_surface.m_controlPointIndices.size());

               auto c = m_surface.m_controlPointIndices.cbegin();
               for(auto i = m_surface.m_parameterPointIndices.cbegin(), end = m_surface.m_parameterPointIndices.cend(); i != end; ++i) {
                    const Point2D& p0 = m_surface.m_parameterPoints[*i];
                    const Point2D& p1 = m_surface.m_parameterPoints[*(++i)];
                    const Point2D& p2 = m_surface.m_parameterPoints[*(++i)];

                    constexpr hpreal nth = 1.0 / (hpreal) t_degree;
                    for(hpint i2 = 0; i2 <= t_degree; ++i2) {
                         hpint i1 = 0;
                         hpint i0 = t_degree - i2;
                         while(i0 >= 0) {
                              vertices.push_back(factory(((hpreal)i0 * p0 + (hpreal)i1 * p1 + (hpreal)i2 * p2) * nth, m_surface.m_controlPoints[*c]));
                              ++i1;
                              --i0;
                              ++c;
                         }
                    }
               }

               std::vector<hpuint> controlPointIndices;
               controlPointIndices.reserve(m_surface.m_controlPointIndices.size());

               for(hpint i = 0, end = m_surface.m_controlPointIndices.size(); i < end; ++i) controlPointIndices.push_back(i);
               std::vector<hpuint> indices = SurfaceSplineUtilsBEZ::template buildTriangleMeshIndices<t_degree>(controlPointIndices);

               return make_triangle_mesh(std::move(vertices), std::move(indices));
          }

     private:
          const SurfaceSplineBEZ& m_surface;

     };//TriangleMeshBuilder

};//SurfaceSplineBEZ
template<class Space>
using ConstantSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 0>;
template<class Space>
using CubicSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 3>;
template<class Space>
using LinearSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 1>;
template<class Space>
using QuadraticSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 2>;
template<class Space>
using QuarticSurfaceSplineBEZ = SurfaceSplineBEZ<Space, 4>;

}//namespace happah

