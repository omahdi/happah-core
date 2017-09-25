// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)


// 2017.09 - Hedwig Amberg    - Implemented nutring triangle mesh.

//TODO: fix placement of padding wall so it's no longer "indented"

#pragma once

#include "NutChain.hpp"

namespace happah {

//DEFINITIONS

class NutRing;

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const NutRing& ring, VertexFactory&& build = VertexFactory());

inline hpuint size(const NutRing& ring);

//DECLARATIONS

class NutRing {
     
public:
     NutRing(hpuint nNuts, hpreal outerLength, hpreal innerLength, hpreal thickness, hpreal padding)
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

};//NutRing

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const NutRing& ring, VertexFactory&& build) {

     /* order of vertices per segment
     
     doubled on top and bottom,              vertices between layers to form
     missing indices are on upper side       vertical walls
       ______________________                 __________26___________
     2                      22               |                       |
     |           10          |               |                       |
     |     6__________14     |               |      ____27_____      |
     |     |           |     |               |     |           |     |
     |     |           |     |               |     |           |     |
     |  4  |           | 18  |               24   25           29    30
     |     |           |     |               |     |           |     |
     |     |___________|     |               |     |____28_____|     |
     |     8          16     |               |                       |
     |           12          |               |                       |
     0 _____________________20               |_______________________|
      \                      /\               \                      /\
       \                    /  \               \                    /  31
        \                  /    \               \                  /    \
     */

     auto indices = Indices();
     auto vertices = std::vector<Vertex>();
     
     auto nNuts = ring.getNumberOfNuts();
     auto outerL = ring.getOuterLength();
     auto innerL = ring.getInnerLength();
     auto padding = ring.getPadding();
     auto thickness = ring.getThickness();
     
     if(nNuts == 1){
          auto chain = NutChain(nNuts, outerL, innerL, thickness, padding);
          return make_triangle_mesh<Vertex>(chain);
     }else if(nNuts == 2){
          auto chain = NutChain(nNuts, outerL, innerL, thickness, padding);
          return make_triangle_mesh<Vertex>(chain);
     }
     
     hpuint nTriangles = 68 * nNuts;
     indices.reserve(3 * nTriangles);
     vertices.reserve(32 * nNuts + 2);
     
     auto toTop = Point3D(0.0f, 0.0f, thickness);
     
     auto add_single_vertex = [&](float v0, float v1, float v2){
          auto offset = Point3D(v0, v1, v2);
          vertices.push_back(build( offset ));
     };
          
     auto add_double_vertex = [&](float v0, float v1, float v2){
          auto offset = Point3D(v0, v1, v2);
          vertices.push_back(build( offset ));
          vertices.push_back(build( offset + toTop ));
     };
     
     auto width = (outerL - innerL) / 2.0;
     auto alpha = (2.0 * M_PI) / nNuts;
     auto beta = (alpha * outerL) / (outerL + padding);
     auto gamma = alpha - beta;
     
     //FIRST SEGMENT VERTICES
     auto seg_height = (outerL/2.0) / (tan(beta/2.0)); //outerL / beta;
     auto seg_width = outerL / 2.0;
     auto inner_width = innerL / 2.0;
     
     auto h_thickness = thickness / 2.0;
     
     add_double_vertex(-seg_width, seg_height, 0);
     add_double_vertex(-seg_width, seg_height + outerL, 0);
     add_double_vertex(-seg_width + (width / 2.0), seg_height + seg_width, 0);
     add_double_vertex(-inner_width, seg_height + outerL - width, 0);
     add_double_vertex(-inner_width, seg_height + width, 0);
     add_double_vertex(0, seg_height + outerL - (width / 2.0), 0);
     add_double_vertex(0, seg_height + (width / 2.0), 0);
     add_double_vertex(inner_width, seg_height + outerL - width, 0);
     add_double_vertex(inner_width, seg_height + width, 0);
     add_double_vertex(seg_width - (width / 2.0), seg_height + seg_width, 0);
     add_double_vertex(seg_width, seg_height, 0);
     add_double_vertex(seg_width, seg_height + outerL, 0);
     
     add_single_vertex(-seg_width, seg_height + seg_width, h_thickness);
     add_single_vertex(-inner_width, seg_height + seg_width, h_thickness);
     add_single_vertex(0, seg_height + outerL, h_thickness);
     add_single_vertex(0, seg_height + outerL - width, h_thickness);
     add_single_vertex(0, seg_height + width, h_thickness);
     add_single_vertex(inner_width, seg_height + seg_width, h_thickness);
     add_single_vertex(seg_width, seg_height + seg_width, h_thickness);
     
     //TODO: calculation incorrect?
     auto delta = (beta / 2.0) + (gamma / 2.0);
     auto p_height = ((padding / 2.0) / tan(gamma / 2.0)) * cos(delta);
     auto p_width = p_height * tan(delta);
     
     add_single_vertex(p_width, p_height, h_thickness);
     
     //ROTATED VERTICES
     int verticesPerSegment = 32;
     
     auto sporeRotationMat = hpmat3x3(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, 0, 0, 1);
     auto rotMat = sporeRotationMat;
     
     for(int i = 1; i < nNuts; i++){
          
          for(int j = 0; j < verticesPerSegment; j++){
               auto old_vertex = Point3D(vertices[j].position);
               auto new_vertex = rotMat * old_vertex;
               add_single_vertex(new_vertex.x, new_vertex.y, new_vertex.z);
          }
          
          rotMat *= sporeRotationMat;
     }
     
     //MIDDLE VERTEX
     add_double_vertex(0, 0, 0);
     
     auto add_single_indices = [&](int v0, int v1, int v2){
          indices.push_back(v0);
          indices.push_back(v1);
          indices.push_back(v2);
     };
          
     auto add_double_indices = [&](int v0, int v1, int v2){
          indices.push_back(v0);
          indices.push_back(v1);
          indices.push_back(v2);
          indices.push_back(v2 + 1);
          indices.push_back(v1 + 1);
          indices.push_back(v0 + 1);
     };
     
     //FIRST SEGMENT TRIANGLES
     add_double_indices(2, 4, 0);
     add_double_indices(2, 6, 4);
     add_double_indices(4, 8, 0);
     add_double_indices(6, 8, 4);
     
     add_double_indices(2, 22, 10);
     add_double_indices(2, 10, 6);
     add_double_indices(22, 14, 10);
     add_double_indices(10, 14, 6);
     
     add_double_indices(8, 16, 12);
     add_double_indices(8, 12, 0);
     add_double_indices(16, 20, 12);
     add_double_indices(12, 20, 0);
     
     add_double_indices(14, 18, 16);
     add_double_indices(22, 18, 14);
     add_double_indices(18, 20, 16);
     add_double_indices(22, 20, 18);
     
     //outer
     add_single_indices(3, 24, 1);
     add_single_indices(3, 2, 24);
     add_single_indices(1, 24, 0);
     add_single_indices(24, 2, 0);
     
     add_single_indices(23, 26, 3);
     add_single_indices(23, 22, 26);
     add_single_indices(3, 26, 2);
     add_single_indices(26, 22, 2);
     
     add_single_indices(21, 30, 23);
     add_single_indices(21, 20, 30);
     add_single_indices(23, 30, 22);
     add_single_indices(30, 20, 22);
     
     //inner
     add_single_indices(9, 25, 7);
     add_single_indices(9, 8, 25);
     add_single_indices(7, 25, 6);
     add_single_indices(25, 8, 6);
     
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
     
     add_single_indices(33, 31, 21);
     add_single_indices(33, 32, 31);
     add_single_indices(21, 31, 20);
     add_single_indices(31, 32, 20);
     
     //SHIFTED INDICES
     int offset = verticesPerSegment;
     int total = verticesPerSegment * nNuts;
     int trianglesPerSegment = 64;
     
     for(int i = 1; i < nNuts; i++){
          
          for(int j = 0; j < trianglesPerSegment * 3; j++){
               indices.push_back( (indices[j] + offset) % total );
          }
          offset += verticesPerSegment;
     }
     
     //MIDDLE AND PADDING TRIANGLES
     offset = 0;
     int middleIndex = nNuts * verticesPerSegment;
     
     for(int i = 0; i < nNuts; i++){
          add_double_indices(offset, 20 + offset, middleIndex);
          add_double_indices(20 + offset, (32 + offset) % total, middleIndex);
          offset += verticesPerSegment;
     }
     
     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const NutRing& ring) { return ring.getNumberOfNuts(); }

}//namespace happah

