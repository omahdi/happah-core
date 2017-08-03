# Copyright 2015 - 2017
#   Pawel Herman - Karlsruhe Institute of Technology
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

lib_LTLIBRARIES = libhappah.la
libhappah_la_SOURCES = \
     happah/Happah.cpp \
     happah/format/obj.cpp \
     happah/format/off.cpp \
     happah/geometries/Circle.cpp \
     happah/geometries/Plane.cpp \
     happah/geometries/Ray.cpp \
     happah/geometries/RectangularTorusChain.cpp \
     happah/geometries/RectangularTorusRing.cpp \
     happah/geometries/Sphere.cpp \
     happah/geometries/SurfaceSplineBEZ.cpp \
     happah/geometries/TriangleGraph.cpp \
     happah/geometries/TriangleMesh.cpp \
     happah/geometries/TriangleMeshUtils.cpp \
     happah/math/ProjectiveStructure.cpp \
     happah/math/Space.cpp \
     happah/math/TriangleRefinementScheme.cpp \
     happah/utils/ControlPointIndexer.cpp \
     happah/utils/SurfaceUtilsBEZ.cpp \
     happah/weighers/ConstantWeigher.cpp
nobase_include_HEADERS = \
     happah/Eigen.hpp \
     happah/Happah.hpp \
     happah/format/default.hpp \
     happah/format/hph.hpp \
     happah/format/obj.hpp \
     happah/format/off.hpp \
     happah/format/xyz.hpp \
     happah/format.hpp \
     happah/geometries/Circle.hpp \
     happah/geometries/Curve.hpp \
     happah/geometries/DiskEmbedding.hpp \
     happah/geometries/LoopBoxSpline.hpp \
     happah/geometries/LoopBoxSplineMesh.hpp \
     happah/geometries/Mesh.hpp \
     happah/geometries/Model.hpp \
     happah/geometries/Plane.hpp \
     happah/geometries/Ray.hpp \
     happah/geometries/RectangularTorusChain.hpp \
     happah/geometries/RectangularTorusRing.hpp \
     happah/geometries/Sphere.hpp \
     happah/geometries/SplineCurve.hpp \
     happah/geometries/SurfaceBEZ.hpp \
     happah/geometries/SurfaceHEZ.hpp \
     happah/geometries/SurfaceSplineBEZ.hpp \
     happah/geometries/SurfaceSplineHEZ.hpp \
     happah/geometries/Triangle.hpp \
     happah/geometries/TriangleArray.hpp \
     happah/geometries/TriangleGraph.hpp \
     happah/geometries/TriangleMesh.hpp \
     happah/geometries/TriangleMeshUtils.hpp \
     happah/geometries/VertexCloud.hpp \
     happah/geometries/Vertex.hpp \
     happah/math/CirclePackingMetric.hpp \
     happah/math/HexagonDecomposition.hpp \
     happah/math/HyperbolicSpace.hpp \
     happah/math/MathUtils.hpp \
     happah/math/PantsDecomposition.hpp \
     happah/math/ProjectiveStructure.hpp \
     happah/math/Space.hpp \
     happah/math/TriangleDecomposition.hpp \
     happah/math/TriangleRefinementScheme.hpp \
     happah/utils/Arrays.hpp \
     happah/utils/ControlPointIndexer.hpp \
     happah/utils/CurveUtilsBEZ.hpp \
     happah/utils/DeindexedArray.hpp \
     happah/utils/IteratorJoiner.hpp \
     happah/utils/MeshUtils.hpp \
     happah/utils/PantsDecomposer.hpp \
     happah/utils/SegmentUtils.hpp \
     happah/utils/ShortestPathFinder.hpp \
     happah/utils/subdividers.hpp \
     happah/utils/SurfaceSplineConstrainerBEZ.hpp \
     happah/utils/SurfaceSubdividerBEZ.hpp \
     happah/utils/SurfaceUtilsBEZ.hpp \
     happah/utils/VertexFactory.hpp \
     happah/utils/visitors.hpp \
     happah/weighers/ConstantWeigher.hpp \
     happah/weighers/EdgeLengthWeigher.hpp \
     happah/weighers/HoleyWallsWeigher.hpp \
     happah/weighers/HoleyWallWeigher.hpp \
     happah/weighers/TraversableEdgeLengthWeigher.hpp
libhappah_la_CPPFLAGS = -I/usr/include/eigen3 -std=c++1z -Wno-unused-label -Wno-unused-parameter -Wno-unused-variable
libhappah_la_LDFLAGS = -version-info 0:0:0

