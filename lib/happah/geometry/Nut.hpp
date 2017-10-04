// Copyright 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.09 - Hedwig Amberg    - Implemented make_triangle_mesh.

#pragma once

#include <glm/gtc/constants.hpp>

namespace happah {

//DEFINITIONS

class Nut;

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const Nut& nut, VertexFactory&& build = VertexFactory());

inline hpuint size(const Nut& nut);

//DECLARATIONS

class Nut {
public:
     Nut(hpuint nSides, hpreal outerRadius, hpreal innerRadius, hpreal thickness)
          : m_nSides(nSides), m_innerRadius(innerRadius), m_outerRadius(outerRadius), m_thickness(thickness) {}

     auto getInnerRadius() const { return m_innerRadius; }

     auto getNumberOfSides() const { return m_nSides; }

     auto getOuterRadius() const { return m_outerRadius; }

     auto getThickness() const { return m_thickness; }

private:
     hpreal m_innerRadius;
     hpuint m_nSides;
     hpreal m_outerRadius;
     hpreal m_thickness;

};//Nut

template<class Vertex, class VertexFactory>
TriangleMesh<Vertex> make_triangle_mesh(const Nut& nut, VertexFactory&& build) {
     auto indices = Indices();
     auto vertices = std::vector<Vertex>();
     auto nSides = nut.getNumberOfSides();
     auto innerRadius = nut.getInnerRadius();
     auto outerRadius = nut.getOuterRadius();
     auto thickness = nut.getThickness();
     auto rotation = glm::angleAxis(glm::two_pi<hpreal>() / nSides, Vector3D(0, 0, 1));
     
     assert(nSides > 2);
     
     vertices.reserve(nSides << 2);
     indices.reserve(3 * (nSides << 3));
     
     vertices.assign({
          build(Point3D(innerRadius, 0, 0)),
          build(Point3D(innerRadius, 0, thickness)),
          build(Point3D(outerRadius, 0, 0)),
          build(Point3D(outerRadius, 0, thickness))
     });
     
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
     
     while(--nSides) {
          auto v = std::end(vertices) - 4;
          vertices.push_back(build(glm::rotate(rotation, v[0].position)));
          vertices.push_back(build(glm::rotate(rotation, v[1].position)));
          vertices.push_back(build(glm::rotate(rotation, v[2].position)));
          vertices.push_back(build(glm::rotate(rotation, v[3].position)));

          for(auto i = std::end(indices) - 24, end = std::end(indices); i != end; ++i) indices.push_back(*i + 4);
     }

     auto i = std::end(indices) - 24;
     i[ 1] = 0;
     i[ 4] = 2;
     i[ 5] = 0;
     i[ 8] = 1;
     i[10] = 1;
     i[11] = 3;
     i[14] = 0;
     i[16] = 0;
     i[17] = 1;
     i[18] = 3;
     i[19] = 2;
     i[22] = 2;
     
     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const Nut& nut) { return nut.getNumberOfSides(); }

}//namespace happah

     
