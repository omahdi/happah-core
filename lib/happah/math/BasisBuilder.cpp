// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/math/BasisBuilder.h"

namespace happah {

BasisBuilder::BasisBuilder(const ProjectiveStructure3D& projectiveStructure, hpuint tetrahedron, Point3D p0, Point3D p1, Point3D p2, const TriangleRefinementScheme& scheme, hpreal epsilon)
     : m_p0(std::move(p0)), m_p1(std::move(p1)), m_p2(std::move(p2)), m_projectiveStructure(projectiveStructure), m_scheme(scheme), m_tetrahedron(tetrahedron), m_refinedProjectiveStructure(refinedProjectiveStructure(epsilon)) {}

const ProjectiveStructure3D& BasisBuilder::getProjectiveStructure() const { return m_projectiveStructure; }

const ProjectiveStructure3D& BasisBuilder::getRefinedProjectiveStructure() const { return m_refinedProjectiveStructure; }

const TriangleRefinementScheme& BasisBuilder::getTriangleRefinementScheme() const { return m_scheme; }

ProjectiveStructure3D BasisBuilder::refinedProjectiveStructure(hpreal epsilon) {
     const TriangleRefinementScheme& scheme = m_scheme;
     auto refinedProjectiveStructure = m_projectiveStructure.refine(scheme, epsilon);
     auto& q0 = scheme.points[scheme.indices[0]];
     auto& q1 = scheme.points[scheme.indices[1]];
     auto& q2 = scheme.points[scheme.indices[2]];
     auto r0 = q0.x * m_p0 + q0.y * m_p1 + q0.z * m_p2;
     auto r1 = q1.x * m_p0 + q1.y * m_p1 + q1.z * m_p2;
     auto r2 = q2.x * m_p0 + q2.y * m_p1 + q2.z * m_p2;
     auto mesh = refinedProjectiveStructure.toTriangleMesh(m_tetrahedron * scheme.getNumberOfTriangles(), r0, r1, r2);
     m_parameterPoints.clear();
     m_parameterPoints.reserve(mesh.getVertices().size());
     for(auto& vertex : mesh.getVertices()) m_parameterPoints.push_back(vertex.position);
     m_parameterPointIndices = mesh.getIndices();
     return refinedProjectiveStructure;
}

void BasisBuilder::setTriangleRefinementScheme(const TriangleRefinementScheme& scheme, hpreal epsilon) {
     m_scheme = std::cref(scheme);
     m_refinedProjectiveStructure = refinedProjectiveStructure();
}

}//namespace happah

