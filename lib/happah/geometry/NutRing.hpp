// Copyright 2017
//   Pawel Herman   - Karlsruhe Institute of Technology - pherman@ira.uka.de
//   Hedwig Amberg  - Karlsruhe Institute of Technology - hedwigdorothea@gmail.com
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// 2017.09 - Hedwig Amberg    - Implemented make_triangle_mesh.
// 2017.10 - Pawel Herman     - Fix bug with positioning of padding vertex.

#pragma once

#include <boost/range.hpp>
#include <glm/gtc/constants.hpp>

namespace happah {

//DEFINITIONS

class NutRing;

template<class Vertex, class VertexFactory = VertexFactory<Vertex> >
TriangleMesh<Vertex> make_triangle_mesh(const NutRing& ring, VertexFactory&& build = VertexFactory());

inline hpuint size(const NutRing& ring);

//DECLARATIONS

class NutRing {
     
public:
     NutRing(hpuint nNuts, hpreal innerRadius, hpreal outerRadius, hpreal thickness, hpreal padding)
          : m_nNuts(nNuts), m_innerRadius(innerRadius), m_outerRadius(outerRadius), m_padding(padding), m_thickness(thickness) {}

     auto& getInnerRadius() const { return m_innerRadius; }

     auto& getNumberOfNuts() const { return m_nNuts; }

     auto& getOuterRadius() const { return m_outerRadius; }

     auto& getPadding() const { return m_padding; }

     auto& getThickness() const { return m_thickness; }

private:
     hpreal m_innerRadius;
     hpuint m_nNuts;
     hpreal m_outerRadius;
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
     auto innerLength = 1.414213562 * ring.getInnerRadius();
     auto outerLength = 1.414213562 * ring.getOuterRadius();
     auto padding = ring.getPadding();
     auto thickness = ring.getThickness();
     auto width = (outerLength - innerLength) / 2.0;
     auto alpha = glm::pi<hpreal>() / nNuts;
     auto height0 = (padding + outerLength * glm::cos(alpha)) / (hpreal(2.0) * glm::sin(alpha));
     auto height1 = (outerLength + padding * glm::cos(alpha)) / (hpreal(2.0) * glm::sin(alpha));
     auto rotation = glm::angleAxis(-hpreal(2.0) * alpha, Vector3D(0, 0, 1));

     assert(nNuts > 2);
     
     indices.reserve(204 * nNuts);
     vertices.reserve(32 * nNuts + 2);
     
     vertices.assign({
          build(Point3D(0,                                  0,                                      0)),
          build(Point3D(0,                                  0,                                      thickness)),
          build(Point3D(-outerLength / 2.0,                 height0,                                0)),
          build(Point3D(-outerLength / 2.0,                 height0,                                thickness)),
          build(Point3D(-outerLength / 2.0,                 height0 + outerLength,                  0)),
          build(Point3D(-outerLength / 2.0,                 height0 + outerLength,                  thickness)),
          build(Point3D(-outerLength / 2.0 + (width / 2.0), height0 + outerLength / 2.0,            0)),
          build(Point3D(-outerLength / 2.0 + (width / 2.0), height0 + outerLength / 2.0,            thickness)),
          build(Point3D(-innerLength / 2.0,                 height0 + outerLength - width,          0)),
          build(Point3D(-innerLength / 2.0,                 height0 + outerLength - width,          thickness)),
          build(Point3D(-innerLength / 2.0,                 height0 + width,                        0)),
          build(Point3D(-innerLength / 2.0,                 height0 + width,                        thickness)),
          build(Point3D(0,                                  height0 + outerLength - (width / 2.0),  0)),
          build(Point3D(0,                                  height0 + outerLength - (width / 2.0),  thickness)),
          build(Point3D(0,                                  height0 + (width / 2.0),                0)),
          build(Point3D(0,                                  height0 + (width / 2.0),                thickness)),
          build(Point3D(innerLength / 2.0,                  height0 + outerLength - width,          0)),
          build(Point3D(innerLength / 2.0,                  height0 + outerLength - width,          thickness)),
          build(Point3D(innerLength / 2.0,                  height0 + width,                        0)),
          build(Point3D(innerLength / 2.0,                  height0 + width,                        thickness)),
          build(Point3D(outerLength / 2.0 - (width / 2.0),  height0 + outerLength / 2.0,            0)),
          build(Point3D(outerLength / 2.0 - (width / 2.0),  height0 + outerLength / 2.0,            thickness)),
          build(Point3D(outerLength / 2.0,                  height0,                                0)),
          build(Point3D(outerLength / 2.0,                  height0,                                thickness)),
          build(Point3D(outerLength / 2.0,                  height0 + outerLength,                  0)),
          build(Point3D(outerLength / 2.0,                  height0 + outerLength,                  thickness)),
          build(Point3D(-outerLength / 2.0,                 height0 + outerLength / 2.0,            thickness / 2.0)),
          build(Point3D(-innerLength / 2.0,                 height0 + outerLength / 2.0,            thickness / 2.0)),
          build(Point3D(0,                                  height0 + outerLength,                  thickness / 2.0)),
          build(Point3D(0,                                  height0 + outerLength - width,          thickness / 2.0)),
          build(Point3D(0,                                  height0 + width,                        thickness / 2.0)),
          build(Point3D(innerLength / 2.0,                  height0 + outerLength / 2.0,            thickness / 2.0)),
          build(Point3D(outerLength / 2.0,                  height0 + outerLength / 2.0,            thickness / 2.0)),
          build(Point3D(glm::sin(alpha) * height1,          glm::cos(alpha) * height1,              thickness / 2.0))
     });
     
     indices.assign({
          hpuint( 4), hpuint( 6), hpuint( 2),
          hpuint( 3), hpuint( 7), hpuint( 5),
          hpuint( 4), hpuint( 8), hpuint( 6),
          hpuint( 7), hpuint( 9), hpuint( 5),
          hpuint( 6), hpuint(10), hpuint( 2),
          hpuint( 3), hpuint(11), hpuint( 7),
          hpuint( 8), hpuint(10), hpuint( 6),
          hpuint( 7), hpuint(11), hpuint( 9),
          hpuint( 4), hpuint(24), hpuint(12),
          hpuint(13), hpuint(25), hpuint( 5),
          hpuint( 4), hpuint(12), hpuint( 8),
          hpuint( 9), hpuint(13), hpuint( 5),
          hpuint(24), hpuint(16), hpuint(12),
          hpuint(13), hpuint(17), hpuint(25),
          hpuint(12), hpuint(16), hpuint( 8),
          hpuint( 9), hpuint(17), hpuint(13),
          hpuint(10), hpuint(18), hpuint(14),
          hpuint(15), hpuint(19), hpuint(11),
          hpuint(10), hpuint(14), hpuint( 2),
          hpuint( 3), hpuint(15), hpuint(11),
          hpuint(18), hpuint(22), hpuint(14),
          hpuint(15), hpuint(23), hpuint(19),
          hpuint(14), hpuint(22), hpuint( 2),
          hpuint( 3), hpuint(23), hpuint(15),
          hpuint(16), hpuint(20), hpuint(18),
          hpuint(19), hpuint(21), hpuint(17),
          hpuint(24), hpuint(20), hpuint(16),
          hpuint(17), hpuint(21), hpuint(25),
          hpuint(20), hpuint(22), hpuint(18),
          hpuint(19), hpuint(23), hpuint(21),
          hpuint(24), hpuint(22), hpuint(20),
          hpuint(21), hpuint(23), hpuint(25),

          //outer walls
          hpuint( 5), hpuint(26), hpuint( 3),
          hpuint( 5), hpuint( 4), hpuint(26),
          hpuint( 3), hpuint(26), hpuint( 2),
          hpuint(26), hpuint( 4), hpuint( 2),
          hpuint(25), hpuint(28), hpuint( 5),
          hpuint(25), hpuint(24), hpuint(28),
          hpuint( 5), hpuint(28), hpuint( 4),
          hpuint(28), hpuint(24), hpuint( 4),
          hpuint(23), hpuint(32), hpuint(25),
          hpuint(23), hpuint(22), hpuint(32),
          hpuint(25), hpuint(32), hpuint(24),
          hpuint(32), hpuint(22), hpuint(24),

          //inner walls
          hpuint(11), hpuint(27), hpuint( 9),
          hpuint(11), hpuint(10), hpuint(27),
          hpuint( 9), hpuint(27), hpuint( 8),
          hpuint(27), hpuint(10), hpuint( 8),
          hpuint( 9), hpuint(29), hpuint(17),
          hpuint( 9), hpuint( 8), hpuint(29),
          hpuint(17), hpuint(29), hpuint(16),
          hpuint(29), hpuint( 8), hpuint(16),
          hpuint(19), hpuint(30), hpuint(11),
          hpuint(19), hpuint(18), hpuint(30),
          hpuint(11), hpuint(30), hpuint(10),
          hpuint(30), hpuint(18), hpuint(10),
          hpuint(17), hpuint(31), hpuint(19),
          hpuint(17), hpuint(16), hpuint(31),
          hpuint(19), hpuint(31), hpuint(18),
          hpuint(31), hpuint(16), hpuint(18),
          hpuint(35), hpuint(33), hpuint(23),
          hpuint(35), hpuint(34), hpuint(33),
          hpuint(23), hpuint(33), hpuint(22),
          hpuint(33), hpuint(34), hpuint(22),

          //center triangles
          hpuint( 2), hpuint(22), hpuint( 0),
          hpuint(22), hpuint(34), hpuint( 0),
          hpuint(23), hpuint( 3), hpuint( 1),
          hpuint(35), hpuint(23), hpuint( 1)
     });
    
     for(auto i = hpuint(1); i < nNuts; ++i) {
          for(auto& vertex : boost::make_iterator_range(std::end(vertices) - 32, std::end(vertices))) vertices.push_back(build(glm::rotate(rotation, vertex.position)));
          for(auto& i : boost::make_iterator_range(std::end(indices) - 204, std::end(indices))) indices.push_back(i + 32);
          auto temp = std::end(indices);
          temp[ -7] = hpuint(0);
          temp[-10] = hpuint(0);
          temp[ -1] = hpuint(1);
          temp[ -4] = hpuint(1);
     }
     for(auto& i : boost::make_iterator_range(std::end(indices) - 204, std::end(indices))) if(i >= 32 * nNuts + 2) i -= 32 * nNuts;
     
     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

inline hpuint size(const NutRing& ring) { return ring.getNumberOfNuts(); }

}//namespace happah

