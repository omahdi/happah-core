# Copyright 2015 - 2017
#   Pawel Herman - Karlsruhe Institute of Technology
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

lib_LTLIBRARIES = libhappah.la
libhappah_la_SOURCES = \
     happah/Happah.cpp \
     happah/format/obj.cpp \
     happah/format/off.cpp \
     happah/geometries/BezierTriangleMesh.cpp \
     happah/geometries/Circle.cpp \
     happah/geometries/Plane.cpp \
     happah/geometries/Sphere.cpp \
     happah/geometries/TriangleGraph.cpp \
     happah/geometries/TriangleMesh.cpp \
     happah/math/ProjectiveStructure.cpp \
     happah/math/TriangleRefinementScheme.cpp
nobase_include_HEADERS = \
     happah/Eigen.hpp \
     happah/Happah.hpp \
     happah/format/default.hpp \
     happah/format/hph.hpp \
     happah/format/obj.hpp \
     happah/format/off.hpp \
     happah/format/xyz.hpp \
     happah/format.hpp \
     happah/geometries/BezierTriangle.hpp \
     happah/geometries/BezierTriangleMesh.hpp \
     happah/geometries/Circle.hpp \
     happah/geometries/LoopBoxSpline.hpp \
     happah/geometries/LoopBoxSplineMesh.hpp \
     happah/geometries/NutChain.hpp \
     happah/geometries/NutRing.hpp \
     happah/geometries/Plane.hpp \
     happah/geometries/Ray.hpp \
     happah/geometries/Sphere.hpp \
     happah/geometries/SplineCurve.hpp \
     happah/geometries/Triangle.hpp \
     happah/geometries/TriangleArray.hpp \
     happah/geometries/TriangleGraph.hpp \
     happah/geometries/TriangleMesh.hpp \
     happah/geometries/VertexCloud.hpp \
     happah/geometries/Vertex.hpp \
     happah/math/CirclePackingMetric.hpp \
     happah/math/functions.hpp \
     happah/math/HyperbolicSpace.hpp \
     happah/math/ProjectiveStructure.hpp \
     happah/math/Space.hpp \
     happah/math/TriangleRefinementScheme.hpp \
     happah/utils/BezierTriangleSubdivider.hpp \
     happah/utils/DeindexedArray.hpp \
     happah/utils/VertexFactory.hpp \
     happah/utils/visitors.hpp
libhappah_la_CPPFLAGS = -I/usr/include/eigen3 -std=c++1z -Wno-unused-label -Wno-unused-parameter -Wno-unused-variable
libhappah_la_LDFLAGS = -version-info 0:0:0

