// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <memory>

#include "happah/geometries/Geometry.h"
#include "happah/geometries/Mesh.h"
#include "happah/math/Space.h"

//NOTE: The technically most accurate name of this class would be LineSegmentMesh but this is shortened to SegmentMesh.
template<class Vertex>
class SegmentMesh : public Geometry1D<typename Vertex::SPACE>, public Mesh<Vertex> {
     using Space = typename Vertex::SPACE;

public:
     using Indices = typename Mesh<Vertex>::Indices;
     using Vertices = typename Mesh<Vertex>::Vertices;

     SegmentMesh(Vertices vertices, Indices indices)
          : Geometry1D<Space>(), Mesh<Vertex>(std::move(vertices), std::move(indices)) {}

     virtual ~SegmentMesh() {}

     class SegmentIterator {
     public:
          SegmentIterator(SegmentMesh<Vertex>& segmentMesh, std::vector<hpuint>::iterator i)
               : m_i(i), m_segmentMesh(segmentMesh) {}
          SegmentIterator(const SegmentIterator& segmentIterator)
               : m_i(segmentIterator.m_i), m_segmentMesh(segmentIterator.m_segmentMesh) {}
         
          SegmentIterator& operator++() { 
               m_i += 2;
               return *this;
          }
          SegmentIterator& operator--() { 
               m_i -= 2;
               return *this;
          }
          SegmentIterator operator++(int) { 
               SegmentIterator segmentIterator(*this);
               ++segmentIterator;
               return segmentIterator;
          }
          SegmentIterator operator--(int) { 
               SegmentIterator segmentIterator(*this);
               --segmentIterator;
               return segmentIterator;
          }
          bool operator!=(const SegmentIterator& segmentIterator) { return segmentIterator.m_i != m_i; }
          std::tuple<Vertex&, Vertex&> operator*() {
               std::vector<hpuint>::iterator i = m_i;
               Vertex& v0 = m_segmentMesh.m_vertices[*i];
               Vertex& v1 = m_segmentMesh.m_vertices[*(++i)];
               return std::make_tuple(std::ref(v0), std::ref(v1)); 
          }

     private:
          std::vector<hpuint>::iterator m_i;
          SegmentMesh<Vertex>& m_segmentMesh;

     };
     typedef SegmentIterator iterator;

     iterator pbegin() { return iterator(*this, this->m_indices.begin()); }
     iterator pend() { return iterator(*this, this->m_indices.end()); }

};
typedef SegmentMesh<VertexP<Space2D> > SegmentMesh2D;
typedef std::shared_ptr<SegmentMesh2D> SegmentMesh2D_ptr;
typedef SegmentMesh<VertexPN<Space3D> > SegmentMesh3D;
typedef std::shared_ptr<SegmentMesh3D> SegmentMesh3D_ptr;

template<class M, class Space = typename M::SPACE, class Vertex = typename M::VERTEX>
struct is_segment_mesh {
     static const bool value = std::is_base_of<SegmentMesh<Vertex>, M>::value && std::is_base_of<typename M::SPACE, Space>::value;
};

//TODO: can we combine SegmentIterator and TriangleIterator into MeshIterator?

//TODO: change name IteratorP
//TODO: generalize to return position or abscissa ordinate and not just for space2d
template<class Vertex, typename = void>
class SegmentIterator2D;

template<class Vertex>
class SegmentIterator2D<Vertex, typename enable_if_absolute_vertex<Vertex, Space2D>::type> {
public:
     typedef hpuint difference_type;

     SegmentIterator2D(std::shared_ptr<SegmentMesh<Vertex> > segmentMesh)
          : m_i(segmentMesh->getIndices().begin()), m_segmentMesh(segmentMesh) {}
     SegmentIterator2D(std::shared_ptr<SegmentMesh<Vertex> > segmentMesh, std::vector<hpuint>::iterator i)
          : m_i(i), m_segmentMesh(segmentMesh) {}
     SegmentIterator2D(const SegmentIterator2D<Vertex>& segmentIterator)
          : m_i(segmentIterator.m_i), m_segmentMesh(segmentIterator.m_segmentMesh) {}

     hpint operator-(const SegmentIterator2D<Vertex>& segmentIterator) const { return (this->m_i-segmentIterator.m_i) >> 1; } 
     SegmentIterator2D<Vertex>& operator++() { 
          m_i += 2;
          return *this;
     }
     SegmentIterator2D<Vertex>& operator--() { 
          m_i -= 2;
          return *this;
     }
     SegmentIterator2D<Vertex> operator++(int) { 
          SegmentIterator2D<Vertex> segmentIterator(*this);
          ++segmentIterator;
          return segmentIterator;
     }
     SegmentIterator2D<Vertex> operator--(int) { 
          SegmentIterator2D<Vertex> segmentIterator(*this);
          --segmentIterator;
          return segmentIterator;
     }
SegmentIterator2D<Vertex> operator[](hpuint index) { 
          SegmentIterator2D<Vertex> segmentIterator(*this);
          segmentIterator.m_i += (index << 1);
          return segmentIterator;
     }

     bool operator!=(const SegmentIterator2D<Vertex>& segmentIterator) { return segmentIterator.m_i != m_i; }
     std::tuple<Point2D&, Point2D&> operator*() {
          std::vector<hpuint>::iterator i = m_i;
          Point2D& p0 = m_segmentMesh->getVertices()[*i].position;
          Point2D& p1 = m_segmentMesh->getVertices()[*(++i)].position;
          return std::make_tuple(std::ref(p0), std::ref(p1)); 
     }

private:
     std::vector<hpuint>::iterator m_i;
     std::shared_ptr<SegmentMesh<Vertex> > m_segmentMesh;

};

