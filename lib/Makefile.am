# Copyright 2015 - 2017
#   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

lib_LTLIBRARIES = libhappah.la
libhappah_la_SOURCES = \
     happah/Happah.cpp \
     happah/format/obj.cpp \
     happah/format/off.cpp \
     happah/geometry/BezierTriangleMesh.cpp \
     happah/geometry/Circle.cpp \
     happah/geometry/Plane.cpp \
     happah/geometry/Sphere.cpp \
     happah/geometry/TriangleGraph.cpp \
     happah/geometry/TriangleMesh.cpp \
     happah/math/CirclePacking.cpp \
     happah/math/HyperbolicSpace.cpp \
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
     happah/geometry/BezierTriangle.hpp \
     happah/geometry/BezierTriangleMesh.hpp \
     happah/geometry/Circle.hpp \
     happah/geometry/DiskEmbedding.hpp \
     happah/geometry/LoopBoxSpline.hpp \
     happah/geometry/LoopBoxSplineMesh.hpp \
     happah/geometry/Nut.hpp \
     happah/geometry/NutChain.hpp \
     happah/geometry/NutRing.hpp \
     happah/geometry/Plane.hpp \
     happah/geometry/Ray.hpp \
     happah/geometry/Sphere.hpp \
     happah/geometry/SplineCurve.hpp \
     happah/geometry/Triangle.hpp \
     happah/geometry/TriangleArray.hpp \
     happah/geometry/TriangleGraph.hpp \
     happah/geometry/TriangleMesh.hpp \
     happah/geometry/VertexCloud.hpp \
     happah/geometry/Vertex.hpp \
     happah/math/CirclePacking.hpp \
     happah/math/functions.hpp \
     happah/math/HyperbolicSpace.hpp \
     happah/math/ProjectiveStructure.hpp \
     happah/math/Space.hpp \
     happah/math/TriangleRefinementScheme.hpp \
     happah/util/BezierTriangleSubdivider.hpp \
     happah/util/ProxyArray.hpp \
     happah/util/RectangularCuboidLayout.hpp \
     happah/util/VertexFactory.hpp \
     happah/util/visitors.hpp
libhappah_la_CPPFLAGS = -I/usr/include/eigen3 -std=c++1z -Wno-unused-label -Wno-unused-parameter -Wno-unused-variable
libhappah_la_LDFLAGS = -version-info 0:0:0

