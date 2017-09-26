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
     
     auto width = (outerL - innerL) / 2.0;
     auto alpha = (2.0 * M_PI) / nNuts;
     auto beta = (alpha * outerL) / (outerL + padding);
     auto gamma = alpha - beta;
     
     //FIRST SEGMENT VERTICES
     auto seg_height = (outerL/2.0) / (tan(beta/2.0)); //outerL / beta;
     auto seg_width = outerL / 2.0;
     auto inner_width = innerL / 2.0;
     auto h_thickness = thickness / 2.0;
     auto delta = (beta / 2.0) + (gamma / 2.0);      //TODO: calculation incorrect?
     auto p_height = ((padding / 2.0) / tan(gamma / 2.0)) * cos(delta); //TODO: calculation incorrect?
     auto p_width = p_height * tan(delta); //TODO: calculation incorrect?
     
     vertices.assign({
     
     build(Point3D(-seg_width,                    seg_height,                             0)),
     build(Point3D(-seg_width,                    seg_height,                             thickness)),
     build(Point3D(-seg_width,                    seg_height + outerL,                    0)),
     build(Point3D(-seg_width,                    seg_height + outerL,                    thickness)),
     build(Point3D(-seg_width + (width / 2.0),    seg_height + seg_width,                 0)),
     build(Point3D(-seg_width + (width / 2.0),    seg_height + seg_width,                 thickness)),
     build(Point3D(-inner_width,                  seg_height + outerL - width,            0)),
     build(Point3D(-inner_width,                  seg_height + outerL - width,            thickness)),
     build(Point3D(-inner_width,                  seg_height + width,                     0)),
     build(Point3D(-inner_width,                  seg_height + width,                     thickness)),
     build(Point3D(0,                             seg_height + outerL - (width / 2.0),    0)),
     build(Point3D(0,                             seg_height + outerL - (width / 2.0),    thickness)),
     build(Point3D(0,                             seg_height + (width / 2.0),             0)),
     build(Point3D(0,                             seg_height + (width / 2.0),             thickness)),
     build(Point3D(inner_width,                   seg_height + outerL - width,            0)),
     build(Point3D(inner_width,                   seg_height + outerL - width,            thickness)),
     build(Point3D(inner_width,                   seg_height + width,                     0)),
     build(Point3D(inner_width,                   seg_height + width,                     thickness)),
     build(Point3D(seg_width - (width / 2.0),     seg_height + seg_width,                 0)),
     build(Point3D(seg_width - (width / 2.0),     seg_height + seg_width,                 thickness)),
     build(Point3D(seg_width,                     seg_height,                             0)),
     build(Point3D(seg_width,                     seg_height,                             thickness)),
     build(Point3D(seg_width,                     seg_height + outerL,                    0)),
     build(Point3D(seg_width,                     seg_height + outerL,                    thickness)),
     
     build(Point3D(-seg_width,     seg_height + seg_width,       h_thickness)),
     build(Point3D(-inner_width,   seg_height + seg_width,       h_thickness)),
     build(Point3D(0,              seg_height + outerL,          h_thickness)),
     build(Point3D(0,              seg_height + outerL - width,  h_thickness)),
     build(Point3D(0,              seg_height + width,           h_thickness)),
     build(Point3D(inner_width,    seg_height + seg_width,       h_thickness)),
     build(Point3D(seg_width,      seg_height + seg_width,       h_thickness)),
     build(Point3D(p_width,        p_height,                     h_thickness))
     });
     
     //ROTATED VERTICES
     hpuint verticesPerSegment = 32;
     auto sporeRotationMat = hpmat3x3(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, 0, 0, 1);
     auto rotMat = sporeRotationMat;
     
     for(hpuint i = 1; i < nNuts; i++){
          for(hpuint j = 0; j < verticesPerSegment; j++){
               auto old_vertex = Point3D(vertices[j].position);
               auto new_vertex = rotMat * old_vertex;
               vertices.push_back(build(new_vertex));
          }
          rotMat *= sporeRotationMat;
     }
     
     //MIDDLE VERTEX
     vertices.push_back(build(Point3D(0, 0, 0)));
     vertices.push_back(build(Point3D(0, 0, thickness)));
     
     //FIRST SEGMENT TRIANGLES
     indices.assign({
     
     hpuint(2), hpuint(4), hpuint(0),
     hpuint(1), hpuint(5), hpuint(3),
     hpuint(2), hpuint(6), hpuint(4),
     hpuint(5), hpuint(7), hpuint(3),
     hpuint(4), hpuint(8), hpuint(0),
     hpuint(1), hpuint(9), hpuint(5),
     hpuint(6), hpuint(8), hpuint(4),
     hpuint(5), hpuint(9), hpuint(7),
     
     hpuint(2), hpuint(22), hpuint(10),
     hpuint(11), hpuint(23), hpuint(3),
     hpuint(2), hpuint(10), hpuint(6),
     hpuint(7), hpuint(11), hpuint(3),
     hpuint(22), hpuint(14), hpuint(10),
     hpuint(11), hpuint(15), hpuint(23),
     hpuint(10), hpuint(14), hpuint(6),
     hpuint(7), hpuint(15), hpuint(11),
     
     hpuint(8), hpuint(16), hpuint(12),
     hpuint(13), hpuint(17), hpuint(9),
     hpuint(8), hpuint(12), hpuint(0),
     hpuint(1), hpuint(13), hpuint(9),
     hpuint(16), hpuint(20), hpuint(12),
     hpuint(13), hpuint(21), hpuint(17),
     hpuint(12), hpuint(20), hpuint(0),
     hpuint(1), hpuint(21), hpuint(13),
     
     hpuint(14), hpuint(18), hpuint(16),
     hpuint(17), hpuint(19), hpuint(15),
     hpuint(22), hpuint(18), hpuint(14),
     hpuint(15), hpuint(19), hpuint(23),
     hpuint(18), hpuint(20), hpuint(16),
     hpuint(17), hpuint(21), hpuint(19),
     hpuint(22), hpuint(20), hpuint(18),
     hpuint(19), hpuint(21), hpuint(23),
     //outer walls
     hpuint(3), hpuint(24), hpuint(1),
     hpuint(3), hpuint(2), hpuint(24),
     hpuint(1), hpuint(24), hpuint(0),
     hpuint(24), hpuint(2), hpuint(0),
     
     hpuint(23), hpuint(26), hpuint(3),
     hpuint(23), hpuint(22), hpuint(26),
     hpuint(3), hpuint(26), hpuint(2),
     hpuint(26), hpuint(22), hpuint(2),
     
     hpuint(21), hpuint(30), hpuint(23),
     hpuint(21), hpuint(20), hpuint(30),
     hpuint(23), hpuint(30), hpuint(22),
     hpuint(30), hpuint(20), hpuint(22),
     //inner walls
     hpuint(9), hpuint(25), hpuint(7),
     hpuint(9), hpuint(8), hpuint(25),
     hpuint(7), hpuint(25), hpuint(6),
     hpuint(25), hpuint(8), hpuint(6),
     
     hpuint(7), hpuint(27), hpuint(15),
     hpuint(7), hpuint(6), hpuint(27),
     hpuint(15), hpuint(27), hpuint(14),
     hpuint(27), hpuint(6), hpuint(14),
     
     hpuint(17), hpuint(28), hpuint(9),
     hpuint(17), hpuint(16), hpuint(28),
     hpuint(9), hpuint(28), hpuint(8),
     hpuint(28), hpuint(16), hpuint(8),
     
     hpuint(15), hpuint(29), hpuint(17),
     hpuint(15), hpuint(14), hpuint(29),
     hpuint(17), hpuint(29), hpuint(16),
     hpuint(29), hpuint(14), hpuint(16),
     
     hpuint(33), hpuint(31), hpuint(21),
     hpuint(33), hpuint(32), hpuint(31),
     hpuint(21), hpuint(31), hpuint(20),
     hpuint(31), hpuint(32), hpuint(20)
     });
     
     //SHIFTED INDICES
     hpuint offset = verticesPerSegment;
     hpuint total = verticesPerSegment * nNuts;
     hpuint trianglesPerSegment = 64;
     
     for(hpuint i = 1; i < nNuts; i++){
          
          for(hpuint j = 0; j < trianglesPerSegment * 3; j++){
               indices.push_back( (indices[j] + offset) % total );
          }
          offset += verticesPerSegment;
     }
     
     //MIDDLE AND PADDING TRIANGLES
     offset = 0;
     hpuint middleIndex = nNuts * verticesPerSegment;
     
     for(hpuint i = 0; i < nNuts; i++){
          auto temp = {
          
          hpuint(offset), hpuint(20 + offset), hpuint(middleIndex),
          hpuint(middleIndex + 1), hpuint(21 + offset), hpuint(offset + 1),
          
          hpuint(20 + offset), hpuint((32 + offset) % total), hpuint(middleIndex),
          hpuint(middleIndex + 1), hpuint((33 + offset) % total), hpuint(21 + offset),
          };
          
          indices.insert(std::end(indices), std::begin(temp), std::end(temp));
          
          offset += verticesPerSegment;
     }
     
     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const NutRing& ring) { return ring.getNumberOfNuts(); }

}//namespace happah

