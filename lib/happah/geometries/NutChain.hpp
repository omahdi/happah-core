// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/geometries/TriangleMesh.hpp"
#include "happah/utils/VertexFactory.hpp"

namespace happah {

//DEFINITIONS

class NutChain;

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const NutChain& chain, VertexFactory&& build = VertexFactory());

inline hpuint size(const NutChain& chain);

//DECLARATIONS

class NutChain {
public:
     NutChain(hpuint nNuts, hpreal outerLength, hpreal innerLength, hpreal thickness, hpreal padding)
          : m_nNuts(nNuts), m_innerLength(innterLength), m_outerLength(outerLength), m_padding(padding), m_thickness(thickness) {}

     hpreal getInnerLength() const { return m_innerLength; }

     hpuint getNumberOfNuts() const { return m_nNuts; }

     hpreal getOuterLength() const { return m_outerLength; }

     hpreal getPadding() const { return m_padding; }

     hpreal getThickness() const { return m_thickness; }

private:
     hpreal m_innerLength;
     hpuint m_nNuts;
     hpreal m_outerLength;
     hpreal m_padding;
     hpreal m_thickness;

};//NutChain

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const NutChain& chain, VertexFactory&& build) {
     auto indices = Indices();
     auto vertices = std::vector<Vertex>();

     //TODO: comment graphic how the vertices are ordered
     //TODO
     //auto nTriangles = ??;
     //indices.reserve(3 * nTriangles);
     //vertices.reserve(??);

     //TODO: generate vertices, then indices (be aware of counterclockwise order)

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const NutChain& chain) { return chain.getNumberOfNuts(); }

}//namespace happah

