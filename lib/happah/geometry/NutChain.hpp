// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.08 - Hedwig Amberg    - Implemented nutchain triangle mesh.

#pragma once

#include "happah/Happah.hpp"
#include "happah/geometry/TriangleMesh.hpp"
#include "happah/util/VertexFactory.hpp"

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
          : m_nNuts(nNuts), m_innerLength(innerLength), m_outerLength(outerLength), m_padding(padding), m_thickness(thickness) {}

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
     
     /*
     IN EVERY NUT
     order of vertices on the underside      oder of vertices between top
     (the uneven indices are the             and bottom side
     corresponding ones on the top)
       ______________________                 _______________________
     0                      20               |                       |
     |           10          |               |                       |
     |     6__________14     |               |      ____27_____      |
     |     |           |     |               |     |           |     |
     |     |           |     |               |     |           |     |
     |  4  |           | 18  |               24   26           29    25
     |     |           |     |               |     |           |     |
     |     |___________|     |               |     |____28_____|     |
     |     8          16     |               |                       |
     |           12          |               |                       |
     2 _____________________22               |_______________________|
     
     OPTIONAL EXTRA VERTICES
     if starting nut (between                if ending nut (between top
     top and bottom)                         and bottom)
      ___________30__________                
     |                       |               |     |___________|     |
     |                       |               |                       |
     |      ___________      |               |                       |
     |     |           |     |               |__________31___________|
     
     if connected nut (the middle vertex is doubled similar to above, 
     the even one is on the underside, the uneven one on top; the 
     others are between)
     
     |     |___________|     |
     |                       |
     |_______________________|
     |                       |
     34        (32/33)       35
     |                       |
     */
     
     int nNuts = chain.getNumberOfNuts();
     float outerL = chain.getOuterLength();
     float innerL = chain.getInnerLength();
     float padding = chain.getPadding();
     float thickness = chain.getThickness();
     hpuint nTriangles = 8 * ( nNuts * 7 + 2 * std::max(nNuts-1, 0) + 1 );
     indices.reserve(3 * nTriangles);
     vertices.reserve( 2 * (nNuts * 15 + 2 * std::max(nNuts-1, 0) + 1) );

     auto toTop = Point3D(0.0f, 0.0f, thickness);
     auto width = (outerL - innerL) / 2.0;
     
     hpuint startIndex = 0;
     hpuint sharedLOffset = 0;
     hpuint sharedROffset = 0;
     
     for (hpuint n = 0; n < nNuts; n++){
          auto start = Point3D(n * outerL + n * padding, 0, 0);
          
          hpuint count = 0;
          
          auto add_single_vertex = [&](float v0, float v1, float v2){
               auto offset = Point3D(v0, v1, v2);
               vertices.push_back(build( start + offset ));
               count++;
          };
          
          auto add_double_vertex = [&](float v0, float v1, float v2){
               auto offset = Point3D(v0, v1, v2);
               vertices.push_back(build( start + offset ));
               vertices.push_back(build( start + offset + toTop ));
               count += 2;
          };
          
          //top/bottom side vertices
          add_double_vertex(0, 0, 0);
          add_double_vertex(outerL, 0, 0);
          add_double_vertex(outerL / 2.0, width / 2.0, 0);
          add_double_vertex(width, width, 0);
          add_double_vertex(outerL - width, width, 0);
          add_double_vertex(width / 2.0, outerL / 2.0, 0);
          add_double_vertex(outerL - (width / 2.0), outerL / 2.0, 0);
          add_double_vertex(width, outerL - width, 0);
          add_double_vertex(outerL - width, outerL - width, 0);
          add_double_vertex(outerL / 2.0, outerL - (width / 2.0), 0);
          add_double_vertex(0, outerL, 0);
          add_double_vertex(outerL, outerL, 0);
          //left+right
          add_single_vertex(outerL / 2.0, 0, thickness / 2.0);
          add_single_vertex(outerL / 2.0, outerL, thickness / 2.0);
          //walls of central hole
          add_single_vertex(outerL / 2.0, width, thickness / 2.0);
          add_single_vertex(width, outerL / 2.0, thickness / 2.0);
          add_single_vertex(outerL - width, outerL / 2.0, thickness / 2.0);
          add_single_vertex(outerL / 2.0, outerL - width, thickness / 2.0);
          
          auto add_single_indices = [&](int v0, int v1, int v2){
               indices.push_back(startIndex + v0);
               indices.push_back(startIndex + v1);
               indices.push_back(startIndex + v2);
          };
          
          auto add_double_indices = [&](int v0, int v1, int v2){
               indices.push_back(startIndex + v0);
               indices.push_back(startIndex + v1);
               indices.push_back(startIndex + v2);
               indices.push_back(startIndex + v2 + 1);
               indices.push_back(startIndex + v1 + 1);
               indices.push_back(startIndex + v0 + 1);
          };
          
          //connect triangles of top/bottom side
          add_double_indices(4, 2, 0);
          add_double_indices(6, 4, 0);
          add_double_indices(8, 2, 4);
          add_double_indices(8, 4, 6);
          add_double_indices(10, 6, 0);
          add_double_indices(12, 2, 8);
          add_double_indices(20, 10, 0);
          add_double_indices(14, 6, 10);
          add_double_indices(16, 12, 8);
          add_double_indices(22, 2, 12);
          add_double_indices(14, 10, 20);
          add_double_indices(22, 12, 16);
          add_double_indices(18, 16, 14);
          add_double_indices(18, 14, 20);
          add_double_indices(22, 16, 18);
          add_double_indices(22, 18, 20);
          //connect left side
          add_single_indices(1, 24, 3);
          add_single_indices(1, 0, 24);
          add_single_indices(3, 24, 2);
          add_single_indices(24, 0, 2);
          //connect right side
          add_single_indices(23, 25, 21);
          add_single_indices(23, 22, 25);
          add_single_indices(21, 25, 20);
          add_single_indices(25, 22, 20);
          //connect walls of hole
          add_single_indices(9, 26, 7);
          add_single_indices(9, 8, 26);
          add_single_indices(7, 26, 6);
          add_single_indices(26, 8, 6);
          
          add_single_indices(7, 27, 15);
          add_single_indices(7, 6, 27);
          add_single_indices(15, 27, 14);
          add_single_indices(27, 6, 14);
          
          add_single_indices(17, 28, 9);
          add_single_indices(17, 16, 28);
          add_single_indices(9, 28, 8);
          add_single_indices(28, 16, 8);
          
          add_single_indices(15, 29, 17);
          add_single_indices(15, 14, 29);
          add_single_indices(17, 29, 16);
          add_single_indices(29, 14, 16);
          
          
          if(n == 0){
          //extra vertex and traingles for starting nut
               add_single_vertex(0, outerL / 2.0, thickness / 2.0);
               add_single_indices(21, count-1, 1);
               add_single_indices(21, 20, count-1);
               add_single_indices(1, count-1, 0);
               add_single_indices(count-1, 20, 0);
          } else {
          //connect vertices of previous padding to own
               add_double_indices(-sharedROffset, -4, -sharedLOffset);
               add_double_indices(-4, 20, 0);
               add_double_indices(-sharedLOffset, -4, 0);
               add_double_indices(-sharedROffset, 20, -4);
               
               add_single_indices(-sharedLOffset + 1, -2, 1);
               add_single_indices(-sharedLOffset + 1, -sharedLOffset, -2);
               add_single_indices(1, -2, 0);
               add_single_indices(-2, -sharedLOffset, 0);
               
               add_single_indices(21, -1, -sharedROffset + 1);
               add_single_indices(21, 20, -1);
               add_single_indices(-sharedROffset + 1, -1, -sharedROffset);
               add_single_indices(-1, 20, -sharedROffset);
          } 
          if ( n == (nNuts-1) ){
          //extra vertex and triangles for ending nut
               add_single_vertex(outerL, outerL / 2.0, thickness / 2.0);
               add_single_indices(3, count-1, 23);
               add_single_indices(3, 2, count-1);
               add_single_indices(23, count-1, 22);
               add_single_indices(count-1, 2, 22);
          } else {
          //create padding vertices, which will get connected by next nut
               add_double_vertex(outerL + (padding / 2.0), outerL / 2.0, 0);
               add_single_vertex(outerL + (padding / 2.0), 0, thickness / 2.0);
               add_single_vertex(outerL + (padding / 2.0), outerL, thickness / 2.0);
               sharedLOffset = count - 2;
               sharedROffset = count - 22;
          }
          
          startIndex += count;
     }
     
     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const NutChain& chain) { return chain.getNumberOfNuts(); }

}//namespace happah

