#ifndef TRIANGLESUBDIVISION_H
#define TRIANGLESUBDIVISION_H

#include "happah/geometries/TriangleMesh.h"

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <tuple>

namespace happah {

template<class Vertex>
class LoopSubdivision
{
private:
    using OutputMesh = TriangleMesh<Vertex, Format::SIMPLE>;
    using InputMesh = TriangleMesh<Vertex, Format::DIRECTED_EDGE>;
    using Edge = std::pair<hpuint, hpuint>;

    const InputMesh& m_mesh;
    const std::vector<Vertex>& m_vertices;
    const std::vector<hpuint>& m_indices;
    std::vector<Vertex> m_new_vertices;
    std::vector<hpuint> m_new_indices;
    std::unordered_map<Edge, hpuint, boost::hash<Edge>> m_edge_index;

    Vertex vertex_rule(hpuint v) const {
        auto neighbors_begin = m_mesh.template cbegin<View::VERTEX, Mode::VERTICES>(v);
        Vertex old_vertex = m_vertices[v];
        // compute mean vertex
        Vertex mean_vertex;
        auto it = neighbors_begin;
        hpuint valence = 0;
        do {
            // TODO: generalize for vertices with multiple attributes
            mean_vertex.position += m_vertices[*it].position;
            ++valence;
            ++it;
        } while (neighbors_begin != it);
        mean_vertex.position *= 1.f / valence;

        // compute weighted average
        float alpha = 3.f / 8.f + 2.f / 8.f * cos(2 * glm::pi<double>() / valence);
        alpha = 5.f / 8.f - alpha * alpha;

        // TODO: generalize for vertices with multiple attributes
        return Vertex(alpha * mean_vertex.position + (1 - alpha) * old_vertex.position);
    }

    Vertex edge_rule(hpuint v, hpuint w, hpuint x, hpuint y) const {
        // TODO: generalize for vertices with multiple attributes
        return (3.f * (m_vertices[v].position + m_vertices[w].position) + (m_vertices[x].position + m_vertices[y].position)) / 8.f;
    }

    hpuint edge_index(hpuint v, hpuint w) {
        if (v >= w) std::swap(v, w);
        return m_edge_index[Edge(v, w)];
    }

public:
    LoopSubdivision(const InputMesh& mesh) : m_mesh(mesh), m_vertices(mesh.getVertices()), m_indices(mesh.getIndices()) {}

    OutputMesh subdivide() {
        hpuint n = m_vertices.size();
        hpuint f = m_mesh.getNumberOfTriangles();
        hpuint e = m_mesh.getNumberOfEdges();

        m_new_vertices.clear();
        m_new_indices.clear();
        m_edge_index.clear();

        m_new_vertices.reserve(n + e / 2);
        m_new_indices.reserve(f * 3 * 4);
        m_edge_index.reserve(e / 2);

        // compute new vertex points
        for (hpuint ivertex = 0; ivertex < n; ++ivertex) {
            m_new_vertices.emplace_back(vertex_rule(ivertex));
        }

        // compute new edge points
        hpuint edge_vertex = m_new_vertices.size();
        for (hpuint iedge = 0; iedge < e; ++iedge) {
            auto edge = m_mesh.getEdge(iedge);
            auto edge_rev = m_mesh.getEdge(edge.opposite);
            hpuint v, w, x, y;
            v = edge.vertex;
            w = edge_rev.vertex;
            // only process each full edge once
            if (v < w) {
                // find neighboring vertices
                x = m_mesh.getEdge(edge.next).vertex;
                y = m_mesh.getEdge(m_mesh.getEdge(edge.opposite).next).vertex;
                m_edge_index.insert(std::make_pair(Edge(v, w), edge_vertex));
                m_new_vertices.emplace_back(edge_rule(v, w, x, y));
                ++edge_vertex;
            }
        }

        // build new index buffer
        for (hpuint iface = 0; iface < f; ++iface) {
            hpuint u = m_indices[3 * iface];     // vertex 0
            hpuint v = m_indices[3 * iface + 1]; // vertex 2
            hpuint w = m_indices[3 * iface + 2]; // vertex 5
            hpuint uv = edge_index(u, v);        // vertex 1
            hpuint vw = edge_index(v, w);        // vertex 4
            hpuint wu = edge_index(w, u);        // vertex 3

            m_new_indices.insert(m_new_indices.end(), {
                                     // see TriangleRefinementScheme.cpp
                                     u,  uv, wu, // 0 1 3
                                     uv,  v, vw, // 1 2 4
                                     uv, vw, wu, // 1 4 3
                                     vw,  w, wu, // 4 5 3
                                 });
        }

        return OutputMesh(m_new_vertices, m_new_indices);
    }
};

}

#endif // TRIANGLESUBDIVISION_H
