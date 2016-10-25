# Copyright 2015
#   Pawel Herman - Karlsruhe Institute of Technology
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

lib_LTLIBRARIES = libhappah.la
libhappah_la_SOURCES = \
     happah/geometries/Plane.cpp \
     happah/geometries/RectangularTorusChain.cpp \
     happah/geometries/RectangularTorusRing.cpp \
     happah/geometries/Sphere.cpp \
     happah/geometries/TriangleMesh.cpp \
     happah/geometries/TriangleMeshUtils.cpp \
     happah/io/readers/ReaderHPH.cpp \
     happah/io/readers/ReaderLOP.cpp \
     happah/io/readers/ReaderOFF.cpp \
     happah/io/writers/WriterHPH.cpp \
     happah/math/BasisBuilder.cpp \
     happah/math/ProjectiveStructure.cpp \
     happah/math/TriangleRefinementScheme.cpp \
     happah/utils/ControlPointIndexer.cpp \
     happah/utils/GeometryUtils.cpp \
     happah/utils/SurfaceExamplesBEZ.cpp \
     happah/weighers/ConstantWeigher.cpp
nobase_include_HEADERS = \
     happah/geometries/Circle.h \
     happah/geometries/Geometry.h \
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
     happah/geometries/Strips.h \
     happah/geometries/SurfaceBEZ.h \
     happah/geometries/Surface.h \
     happah/geometries/SurfaceHEZ.h \
     happah/geometries/SurfaceSBB.h \
     happah/geometries/SurfaceSplineBEZ.h \
     happah/geometries/SurfaceSplineHEZ.h \
     happah/geometries/Trapezoidation.h \
     happah/geometries/TriangleMesh.h \
     happah/geometries/TriangleMeshUtils.h \
     happah/geometries/VertexCloud.h \
     happah/geometries/Vertex.h \
     happah/io/readers/ReaderHPH.h \
     happah/io/readers/ReaderLOP.h \
     happah/io/readers/ReaderOFF.h \
     happah/io/readers/ReaderUtils.h \
     happah/io/writers/WriterHPH.h \
     happah/io/writers/WriterOFF.h \
     happah/math/BasisBuilder.h \
     happah/math/HexagonDecomposition.h \
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
     happah/utils/SurfaceExamplesBEZ.h \
     happah/utils/SurfaceSplineConstrainerBEZ.h \
     happah/utils/SurfaceSplineSubdividerBEZ.h \
     happah/utils/SurfaceSplineUtilsBEZ.h \
     happah/utils/SurfaceSubdividerBEZ.h \
     happah/utils/SurfaceUtilsBEZ.h \
     happah/utils/VertexFactory.h \
     happah/utils/Visitor.h \
     happah/weighers/ConstantWeigher.h \
     happah/weighers/EdgeLengthWeigher.h \
     happah/weighers/HoleyWallsWeigher.h \
     happah/weighers/HoleyWallWeigher.h \
     happah/weighers/TraversableEdgeLengthWeigher.h \
     happah/Eigen.h \
     happah/Happah.h
libhappah_la_CPPFLAGS = -I/usr/include/eigen3 -std=c++1y -Wno-unused-label -Wno-unused-parameter -Wno-unused-variable -fcilkplus
libhappah_la_LDFLAGS = -version-info 0:0:0

