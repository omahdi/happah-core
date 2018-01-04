# Copyright 2017
#   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

AM_CPPFLAGS = -I$(top_srcdir)/lib -I/usr/include/eigen3 -std=c++1z
LDADD = $(top_srcdir)/lib/libhappah.la -lboost_iostreams

check_PROGRAMS = \
     happah_geometry_BezierTriangleMesh_DeltasEnumerator \
     happah_geometry_BezierTriangleMesh_NablasEnumerator \
     happah_geometry_BezierTriangleMesh_RowEnumerator \
     happah_geometry_BezierTriangleMesh_is_g1 \
     happah_geometry_Circle_intersect \
     happah_geometry_QuadMesh_SpokesEnumerator \
     happah_geometry_TriangleGraph_undegenerate \
     happah_geometry_TriangleGraph_analyze
happah_geometry_BezierTriangleMesh_DeltasEnumerator_SOURCES = \
     happah/geometry/BezierTriangleMesh/DeltasEnumerator.cpp
happah_geometry_BezierTriangleMesh_NablasEnumerator_SOURCES = \
     happah/geometry/BezierTriangleMesh/NablasEnumerator.cpp
happah_geometry_BezierTriangleMesh_RowEnumerator_SOURCES = \
     happah/geometry/BezierTriangleMesh/RowEnumerator.cpp
happah_geometry_BezierTriangleMesh_is_g1_SOURCES = \
     happah/geometry/BezierTriangleMesh/is_g1.cpp
happah_geometry_Circle_intersect_SOURCES = \
     happah/geometry/Circle/intersect.cpp
happah_geometry_QuadMesh_SpokesEnumerator_SOURCES = \
     happah/geometry/QuadMesh/SpokesEnumerator.cpp
happah_geometry_TriangleGraph_undegenerate_SOURCES = \
     happah/geometry/TriangleGraph/undegenerate.cpp
happah_geometry_TriangleGraph_analyze_SOURCES = \
     happah/geometry/TriangleGraph/analyze.cpp

TESTS = $(check_PROGRAMS)

