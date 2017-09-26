// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.09 - Hedwig Amberg    - Implemented general nut triangle mesh.

//TODO: what should inner and outer length represent

#pragma once

namespace happah {

//DEFINITIONS

class Nut;

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const Nut& nut, VertexFactory&& build = VertexFactory());

inline hpuint size(const Nut& nut);

//DECLARATIONS

class Nut {

public:
     Nut(hpuint nSides, hpreal outerLength, hpreal innerLength, hpreal thickness)
          : m_nSides(nSides), m_innerLength(innerLength), m_outerLength(outerLength), m_thickness(thickness) {}

     hpreal getInnerLength() const { return m_innerLength; }

     hpuint getNumberOfSides() const { return m_nSides; }

     hpreal getOuterLength() const { return m_outerLength; }

     hpreal getThickness() const { return m_thickness; }

private:
     hpreal m_innerLength;
     hpuint m_nSides;
     hpreal m_outerLength;
     hpreal m_thickness;

};//Nut

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const Nut& nut, VertexFactory&& build) {

     auto indices = Indices();
     auto vertices = std::vector<Vertex>();
     
     auto nSides = nut.getNumberOfSides();
     auto outerL = nut.getOuterLength();
     auto innerL = nut.getInnerLength();
     auto thickness = nut.getThickness();
     
     assert(nSides >= 3);
     
     hpuint nTriangles = 8 * nSides;
     indices.reserve(3 * nTriangles);
     vertices.reserve(4 * nSides);
     
     auto alpha = (M_PI * (nSides - 2)) / nSides;
     auto inner_height = innerL / 2.0;
     auto outer_height = outerL / 2.0;
     auto inner_width = inner_height / tan(alpha / 2.0);
     auto outer_width = outer_height / tan(alpha / 2.0);
     
     //FIRST SEGMENT
     vertices.assign({
     build(Point3D(-inner_width,     inner_height,    0)),
     build(Point3D(-inner_width,     inner_height,    thickness)),
     build(Point3D(-outer_width,     outer_height,    0)),
     build(Point3D(-outer_width,     outer_height,    thickness))
     });
     
     auto angle = (2.0 * M_PI) / nSides;
     auto sporeRotationMat = hpmat3x3(cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1);
     auto rotMat = sporeRotationMat;
     
     //ROTATED SEGMENTS
     for(hpuint i = 1; i < nSides; i++){
          for(hpuint j = 0; j < 4; j++){
               auto old_vertex = Point3D(vertices[j].position);
               auto new_vertex = rotMat * old_vertex;
               vertices.push_back(build(new_vertex));
          }
          rotMat *= sporeRotationMat;
     }
     
     //FIRST SEGMENT TRIANGLES
     indices.assign({
     hpuint(2), hpuint(4), hpuint(0),
     hpuint(2), hpuint(6), hpuint(4),
     hpuint(3), hpuint(1), hpuint(5),
     hpuint(3), hpuint(5), hpuint(7),
     
     hpuint(1), hpuint(0), hpuint(4),
     hpuint(1), hpuint(4), hpuint(5),
     hpuint(7), hpuint(6), hpuint(3),
     hpuint(3), hpuint(6), hpuint(2)
     });
     
     //SHIFTED TRIAGLES
     hpuint offset = 4;
     hpuint total = 4 * nSides;
     
     for(hpuint i = 1; i < nSides; i++){
          for(hpuint j = 0; j < 24; j++){
               indices.push_back( (indices[j] + offset) % total );
          }
          offset += 4;
     }
     
     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const Nut& nut) { return nut.getNumberOfSides(); }

}//namespace happah

     