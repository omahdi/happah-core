// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/geometries/Model.hpp"
#include "happah/geometries/SegmentMesh.hpp"
#include "happah/utils/Arrays.hpp"

namespace happah {

template<class Vertex>
class Loops : public Model<Vertex> {
public:
     using IndicesArrays = Arrays<hpuint>;
     using Vertices = typename Model<Vertex>::Vertices;

     Loops(Vertices vertices, IndicesArrays loops)
          : Model<Vertex>(std::move(vertices)), m_loops(std::move(loops)) {}

     virtual ~Loops() {}

     hpuint getNumberOfLoops() const { return m_loops.size(); }

     SegmentMesh<Vertex> toSegmentMesh() const {
          std::vector<Vertex> vertices;
          vertices.reserve(m_loops.data().size());
          std::vector<hpuint> indices;
          indices.reserve(m_loops.data().size() << 1);

          hpuint index = 0, index0 = 0;
          for(auto loop : m_loops) {
               auto i = loop.first;
               auto end = loop.second;
               if(i == end) continue;
               --end;
               while(i != end) {
                    vertices.push_back(this->m_vertices[*i]);
                    indices.push_back(index);
                    indices.push_back(++index);
                    ++i;
               }
               vertices.push_back(this->m_vertices[*i]);
               indices.push_back(index);
               indices.push_back(index0);
               ++index;
               index0 = index;
          }

          return { std::move(vertices), std::move(indices) };
     } 

private:
     IndicesArrays m_loops;

};

}//namespace happah

