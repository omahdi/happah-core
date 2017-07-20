// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/optional.hpp>
#include <vector>

#include "happah/Happah.hpp"
#include "happah/geometries/Model.hpp"
#include "happah/utils/Arrays.hpp"

namespace happah { 

/*
 * @section DESCRIPTION
 *
 * A mesh is a set of patches.  All patches have the same number of vertices. For
 * example, a triangle mesh consists of patches that are triangles.
 */

template<class Vertex>
class Mesh : public Model<Vertex> {
public:
     using IndicesArrays = Arrays<hpuint>;
     using Vertices = typename Model<Vertex>::Vertices;

     Mesh() {}

     Mesh(Vertices vertices, Indices indices)
          : Model<Vertex>(std::move(vertices)), m_indices(std::move(indices)), m_loops(boost::none) {}

     virtual ~Mesh() {}

     boost::optional<hpuint> getGenus() const {
          if(m_loops) return m_loops->size() >> 1;
          else return boost::none;
     }

     const boost::optional<IndicesArrays>& getHandleTunnelLoops() const { return m_loops; }

     const Indices& getIndices() const { return m_indices; }

     Indices& getIndices() { return m_indices; }

     void setHandleTunnelLoops(IndicesArrays loops) { m_loops = std::move(loops); }

protected:
     Indices m_indices;
     boost::optional<IndicesArrays> m_loops;

};

template<class M, class Space = typename M::SPACE, class Vertex = typename M::VERTEX>
struct is_mesh : std::integral_constant<bool, std::is_base_of<Mesh<Vertex>, M>::value && std::is_base_of<typename M::SPACE, Space>::value> {};

template<class Mesh, class Space = typename Mesh::SPACE, class Vertex = typename Mesh::VERTEX>
struct is_absolute_mesh : std::integral_constant<bool, is_mesh<Mesh, Space, Vertex>::value && is_absolute_vertex<Vertex>::value> {};

template<class Mesh, class Space = typename Mesh::SPACE, class Vertex = typename Mesh::VERTEX>
struct is_relative_mesh : std::integral_constant<bool, is_mesh<Mesh, Space, Vertex>::value && is_relative_vertex<Vertex>::value> {};

template<class Mesh, class Space = typename Mesh::SPACE, class Vertex = typename Mesh::VERTEX>
struct enable_if_mesh : std::enable_if<is_mesh<Mesh, Space, Vertex>::value> {};

template<class Mesh, class Space = typename Mesh::SPACE, class Vertex = typename Mesh::VERTEX>
struct enable_if_absolute_mesh : std::enable_if<is_absolute_mesh<Mesh, Space, Vertex>::value> {};

template<class Mesh, class Space = typename Mesh::SPACE, class Vertex = typename Mesh::VERTEX>
struct enable_if_relative_mesh : std::enable_if<is_relative_mesh<Mesh, Space, Vertex>::value> {};

}//namespace happah

