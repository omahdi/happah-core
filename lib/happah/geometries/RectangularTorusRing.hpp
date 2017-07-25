// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/geometries/Circle.hpp"
#include "happah/geometries/TriangleMesh.hpp"
#include "happah/geometries/TriangleMeshUtils.hpp"
#include "happah/geometries/Vertex.hpp"
#include "happah/utils/VertexFactory.hpp"

namespace happah {

class RectangularTorusRing {
public:
     RectangularTorusRing(hpuint nHandles, hpreal side, hpreal spacing, hpreal hole);

     template<class Vertex = VertexP3, class VertexFactory = VertexFactory<Vertex> >
     TriangleMesh<Vertex> toTriangleMesh(VertexFactory&& factory = VertexFactory()) const {
          //TODO: m_space < epsilon or m_hole < epsilon
          std::vector<Vertex> vertices = getVertices(factory);
          std::vector<hpuint> indices(m_nHandles * 108);

          hpuint nSingleFaceIndices = m_nHandles * 30;
          BaseTrianglesIndicesSetter<Vertex, VertexFactory> visitor(vertices, indices, m_r0, m_r1, m_hole, m_nFaceVertices, (nSingleFaceIndices << 1), factory);
          vertices.reserve(vertices.size() + 6 * m_nHandles);
          Circle::visitNeighborhoodIndices(vertices, m_nHandles, 1, visitor);
          TriangleMeshUtils::extrude(vertices, indices.begin(), nSingleFaceIndices, indices.begin() + nSingleFaceIndices, m_side, factory);

          return TriangleMesh<Vertex>(std::move(vertices), std::move(indices));
     } 

private:
     hpreal m_baseRadius;
     hpreal m_hole;
     hpuint m_nFaceVertices;
     hpuint m_nHandles;
     hpreal m_r0, m_r1;
     hpreal m_side;
     hpreal m_spacing;
     hpreal m_theta;

     template<class Vertex, class VertexFactory>
     class BaseTrianglesIndicesSetter {
     public:
          BaseTrianglesIndicesSetter(std::vector<Vertex>& vertices, std::vector<hpuint>& indices, hpreal r0, hpreal r1, hpreal hole, hpuint deltav, hpuint deltai, VertexFactory& factory)
               : m_deltav(deltav), m_faceIndices(indices.begin()), m_factory(factory), m_i(vertices.size()), m_r0(r0), m_r1(r1), m_s0((1.0 - hole) / 2.0), m_s1(1.0 - m_s0), m_sideIndices(m_faceIndices + deltai), m_vertices(vertices) {}

          void visit(const Vertex& n, hpuint i, const Vertex& n1, hpuint i1, bool o1) {}
          //NOTE: ocw is true if the subwedge between n and ncw is the second subwedge in a pair.
          void visit(const Vertex& n, hpuint i, const Vertex& nccw, hpuint iccw, const Vertex& ncw, hpuint icw, bool ocw) {
               //base indices

               (*m_faceIndices) = 0; 
               ++m_faceIndices;
               (*m_faceIndices) = i; 
               ++m_faceIndices;
               (*m_faceIndices) = iccw; 
               ++m_faceIndices;

               if(ocw) {
                    Point3D p0 = Point3D(m_r0 * n.position.x - m_r1 * n.position.y, m_r1 * n.position.x + m_r0 * n.position.y, 0.0);
                    Point3D p1 = Point3D(m_r0 * nccw.position.x + m_r1 * nccw.position.y, m_r0 * nccw.position.y - m_r1 * nccw.position.x, 0.0);
                    m_vertices.push_back(m_factory(p0));
                    m_vertices.push_back(m_factory(p1));

                    Point3D q0, q1;
                    
                    q0 = Point3D(m_s0 * n.position.x + m_s1 * nccw.position.x, m_s0 * n.position.y + m_s1 * nccw.position.y, 0.0);
                    q1 = Point3D(m_s0 * p0.x + m_s1 * p1.x, m_s0 * p0.y + m_s1 * p1.y, 0.0);
                    m_vertices.push_back(m_factory(Point3D(m_s0 * q1.x + m_s1 * q0.x, m_s0 * q1.y + m_s1 * q0.y, 0.0)));
                    m_vertices.push_back(m_factory(Point3D(m_s1 * q1.x + m_s0 * q0.x, m_s1 * q1.y + m_s0 * q0.y, 0.0)));

                    q0 = Point3D(m_s1 * n.position.x + m_s0 * nccw.position.x, m_s1 * n.position.y + m_s0 * nccw.position.y, 0.0);
                    q1 = Point3D(m_s1 * p0.x + m_s0 * p1.x, m_s1 * p0.y + m_s0 * p1.y, 0.0);
                    m_vertices.push_back(m_factory(Point3D(m_s0 * q1.x + m_s1 * q0.x, m_s0 * q1.y + m_s1 * q0.y, 0.0)));
                    m_vertices.push_back(m_factory(Point3D(m_s1 * q1.x + m_s0 * q0.x, m_s1 * q1.y + m_s0 * q0.y, 0.0)));

                    hpuint iv = i + m_deltav, iccwv = iccw + m_deltav, icwv = icw + m_deltav;
                    hpuint mi1 = m_i + 1, mi2 = m_i + 2, mi3 = m_i + 3, mi4 = m_i + 4, mi5 = m_i + 5;
                    hpuint miv = m_i + m_deltav, mi1v = mi1 + m_deltav, mi2v = mi2 + m_deltav, mi3v = mi3 + m_deltav, mi4v = mi4 + m_deltav, mi5v = mi5 + m_deltav;

                    //square with hole indices

                    (*m_faceIndices) = iccw; 
                    ++m_faceIndices;
                    (*m_faceIndices) = i; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi2; 
                    ++m_faceIndices;

                    (*m_faceIndices) = mi2; 
                    ++m_faceIndices;
                    (*m_faceIndices) = i; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi4; 
                    ++m_faceIndices;

                    (*m_faceIndices) = mi4; 
                    ++m_faceIndices;
                    (*m_faceIndices) = i; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi5; 
                    ++m_faceIndices;

                    (*m_faceIndices) = mi5; 
                    ++m_faceIndices;
                    (*m_faceIndices) = i; 
                    ++m_faceIndices;
                    (*m_faceIndices) = m_i; 
                    ++m_faceIndices;

                    (*m_faceIndices) = mi5; 
                    ++m_faceIndices;
                    (*m_faceIndices) = m_i; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi1; 
                    ++m_faceIndices;

                    (*m_faceIndices) = mi3; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi5; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi1; 
                    ++m_faceIndices;

                    (*m_faceIndices) = mi2; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi3; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi1; 
                    ++m_faceIndices;

                    (*m_faceIndices) = iccw; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi2; 
                    ++m_faceIndices;
                    (*m_faceIndices) = mi1; 
                    ++m_faceIndices;

                    //square sides

                    (*m_sideIndices) = m_i;
                    ++m_sideIndices;
                    (*m_sideIndices) = i;
                    ++m_sideIndices;
                    (*m_sideIndices) = miv;
                    ++m_sideIndices;

                    (*m_sideIndices) = i;
                    ++m_sideIndices;
                    (*m_sideIndices) = iv;
                    ++m_sideIndices;
                    (*m_sideIndices) = miv;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi1;
                    ++m_sideIndices;
                    (*m_sideIndices) = m_i;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi1v;
                    ++m_sideIndices;

                    (*m_sideIndices) = m_i;
                    ++m_sideIndices;
                    (*m_sideIndices) = miv;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi1v;
                    ++m_sideIndices;

                    (*m_sideIndices) = iccw;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi1;
                    ++m_sideIndices;
                    (*m_sideIndices) = iccwv;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi1;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi1v;
                    ++m_sideIndices;
                    (*m_sideIndices) = iccwv;
                    ++m_sideIndices;

                    //hole sides

                    (*m_sideIndices) = mi3;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi2;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi3v;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi2;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi2v;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi3v;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi5;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi3;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi5v;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi3;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi3v;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi5v;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi4;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi5;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi4v;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi5;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi5v;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi4v;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi2;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi4;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi2v;
                    ++m_sideIndices;

                    (*m_sideIndices) = mi4;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi4v;
                    ++m_sideIndices;
                    (*m_sideIndices) = mi2v;
                    ++m_sideIndices;

                    //spacing side

                    (*m_sideIndices) = i;
                    ++m_sideIndices;
                    (*m_sideIndices) = icw;
                    ++m_sideIndices;
                    (*m_sideIndices) = iv;
                    ++m_sideIndices;

                    (*m_sideIndices) = icw;
                    ++m_sideIndices;
                    (*m_sideIndices) = icwv;
                    ++m_sideIndices;
                    (*m_sideIndices) = iv;
                    ++m_sideIndices;

                    m_i += 6;
               }
          }
     
     private:
          hpuint m_deltav;
          std::vector<hpuint>::iterator m_faceIndices;
          VertexFactory& m_factory;
          hpuint m_i;
          hpreal m_r0, m_r1;
          hpreal m_s0, m_s1;
          std::vector<hpuint>::iterator m_sideIndices;
          std::vector<Vertex>& m_vertices;

     };

     template<class VertexFactory>
     std::vector<typename VertexFactory::PRODUCT> getVertices(VertexFactory& factory) const {
          std::vector<typename VertexFactory::PRODUCT> vertices;
          vertices.reserve((m_nFaceVertices << 1) + (m_nHandles << 1));

          auto visit = [&](hpreal x, hpreal y) { return vertices.push_back(factory(Point3D(x, y, 0.0))); };
          visit(0.0, 0.0);
          Circle::sample(m_baseRadius, m_nHandles, m_theta, visit);

          return std::move(vertices);
     }

};//RectangularTorusRing

}//namespace happah

