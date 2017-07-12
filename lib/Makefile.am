# Copyright 2015 - 2017
#   Pawel Herman - Karlsruhe Institute of Technology
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

lib_LTLIBRARIES = libhappah.la
libhappah_la_SOURCES = \
     happah/Happah.cpp \
     happah/format/obj.cpp \
     happah/format/off.cpp \
     happah/geometries/Plane.cpp \
     happah/geometries/Ray.cpp \
     happah/geometries/RectangularTorusChain.cpp \
     happah/geometries/RectangularTorusRing.cpp \
     happah/geometries/Sphere.cpp \
     happah/geometries/SurfaceSplineBEZ.cpp \
     happah/geometries/TriangleMesh.cpp \
     happah/geometries/TriangleMeshUtils.cpp \
     happah/math/CirclePackingMetric.cpp \
     happah/math/BasisBuilder.cpp \
     happah/math/ProjectiveStructure.cpp \
     happah/math/Space.cpp \
     happah/math/TriangleRefinementScheme.cpp \
     happah/utils/ControlPointIndexer.cpp \
     happah/utils/GeometryUtils.cpp \
     happah/utils/SurfaceUtilsBEZ.cpp \
     happah/weighers/ConstantWeigher.cpp
nobase_include_HEADERS = \
     happah/Eigen.h \
     happah/Happah.h \
     happah/format/default.h \
     happah/format/hph.h \
     happah/format/obj.h \
     happah/format/off.h \
     happah/format/xyz.h \
     happah/format.h \
     happah/geometries/converters.h \
     happah/geometries/Circle.h \
     happah/geometries/Curve.h \
     happah/geometries/DiskEmbedding.h \
     happah/geometries/Geometry.h \
     happah/geometries/LoopBoxSpline.h \
     happah/geometries/LoopBoxSplineMesh.h \
     happah/geometries/Loops.h \
     happah/geometries/Mesh.h \
     happah/geometries/Model.h \
     happah/geometries/Plane.h \
     happah/geometries/QuadMesh.h \
     happah/geometries/Ray.h \
     happah/geometries/RectangularTorusChain.h \
     happah/geometries/RectangularTorusRing.h \
     happah/geometries/SegmentLoops.h \
     happah/geometries/SegmentMesh.h \
     happah/geometries/SegmentStrips.h \
     happah/geometries/Sphere.h \
     happah/geometries/SplineCurve.h \
     happah/geometries/Strips.h \
     happah/geometries/SurfaceBEZ.h \
     happah/geometries/Surface.h \
     happah/geometries/SurfaceHEZ.h \
     happah/geometries/SurfaceSplineBEZ.h \
     happah/geometries/SurfaceSplineHEZ.h \
     happah/geometries/Trapezoidation.h \
     happah/geometries/Triangle.h \
     happah/geometries/TriangleArray.h \
     happah/geometries/TriangleMesh.h \
     happah/geometries/TriangleMeshUtils.h \
     happah/geometries/VertexCloud.h \
     happah/geometries/Vertex.h \
     happah/math/BasisBuilder.h \
     happah/math/CirclePackingMetric.h \
     happah/math/HexagonDecomposition.h \
     happah/math/HyperbolicSpace.h \
     happah/math/MathUtils.h \
     happah/math/PantsDecomposition.h \
     happah/math/ProjectiveStructure.h \
     happah/math/RigidAffineTransformation.h \
     happah/math/Space.h \
     happah/math/TriangleDecomposition.h \
     happah/math/TriangleRefinementScheme.h \
     happah/utils/Arrays.h \
     happah/utils/ControlPointIndexer.h \
     happah/utils/CurveUtilsBEZ.h \
     happah/utils/DeindexedArray.h \
     happah/utils/GeometryUtils.h \
     happah/utils/InterpolatorPCT.h \
     happah/utils/InterpolatorSCT.h \
     happah/utils/IteratorJoiner.h \
     happah/utils/MeshUtils.h \
     happah/utils/PantsDecomposer.h \
     happah/utils/ProjectiveStructureUtils.h \
     happah/utils/SegmentUtils.h \
     happah/utils/ShortestPathFinder.h \
     happah/utils/subdividers.h \
     happah/utils/SurfaceSplineConstrainerBEZ.h \
     happah/utils/SurfaceSubdividerBEZ.h \
     happah/utils/SurfaceUtilsBEZ.h \
     happah/utils/VertexFactory.h \
     happah/utils/Visitor.h \
     happah/utils/visitors.h \
     happah/weighers/ConstantWeigher.h \
     happah/weighers/EdgeLengthWeigher.h \
     happah/weighers/HoleyWallsWeigher.h \
     happah/weighers/HoleyWallWeigher.h \
     happah/weighers/TraversableEdgeLengthWeigher.h
libhappah_la_CPPFLAGS = -I/usr/include/eigen3 -std=c++1z -Wno-unused-label -Wno-unused-parameter -Wno-unused-variable
libhappah_la_LDFLAGS = -version-info 0:0:0

