// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/geometries/Model.hpp"
#include "happah/geometries/SegmentMesh.hpp"
#include "happah/utils/ArraysIterator.hpp"//NOTE arraysiterator has been trashed; reimplement this class

template<class Vertex>
class Strips : public Model<Vertex> {
public:
     typedef std::shared_ptr<std::vector<hpuint> > Indices;
     typedef std::vector<hpuint>::iterator IteratorI;
     typedef ArraysIterator<Indices, IteratorI, Indices, IteratorI> IteratorS;
     typedef typename Model<Vertex>::Vertices Vertices;

     Strips(Vertices vertices, Indices indices, Indices ends)
          : Model<Vertex>(vertices), m_ends(ends), m_indices(indices) {}

     virtual ~Strips() {}

     IteratorS beginS() { return IteratorS(m_indices, m_ends, m_ends->begin()); }

     IteratorS endS() { return IteratorS(m_indices, m_ends, m_ends->end()); }

     Indices getEnds() const { return m_ends; }

     Indices getIndices() const { return m_indices; }

     SegmentMesh<Vertex>* toSegmentMesh() const {
          typedef typename SegmentMesh<Vertex>::Indices Indices;
          typedef typename SegmentMesh<Vertex>::Vertices Vertices;

          std::vector<Vertex>* vertices = new std::vector<Vertex>();
          vertices->reserve(m_indices->size());
          std::vector<hpuint>* indices = new std::vector<hpuint>();
          indices->reserve((m_indices->size() - m_ends->size()) << 1);

          std::vector<hpuint>::iterator e = m_ends->begin();
          hpuint index = 0, index1 = *e - 1;
          for(hpuint i : *m_indices) {
               vertices->push_back((*(this->m_vertices))[i]);
               if(index == index1) {
                    ++e;
                    index1 = *e - 1;
                    ++index;
               } else {
                    indices->push_back(index);
                    indices->push_back(++index);
               }
          }

          return new SegmentMesh<Vertex>(Vertices(vertices), Indices(indices));
     } 

private:
     Indices m_ends;
     Indices m_indices;
     
};

