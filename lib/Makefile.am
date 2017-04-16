# Copyright 2015 - 2017
#   Pawel Herman - Karlsruhe Institute of Technology
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

lib_LTLIBRARIES = libhappah.la
libhappah_la_SOURCES = \
     happah/Happah.cpp \
     happah/geometries/Plane.cpp \
     happah/geometries/Ray.cpp \
     happah/geometries/RectangularTorusChain.cpp \
     happah/geometries/RectangularTorusRing.cpp \
     happah/geometries/Sphere.cpp \
     happah/geometries/SurfaceSplineBEZ.cpp \
     happah/geometries/TriangleMesh.cpp \
     happah/geometries/TriangleMeshUtils.cpp \
     happah/math/BasisBuilder.cpp \
     happah/math/nlopt.cpp \
     happah/math/ProjectiveStructure.cpp \
     happah/math/Space.cpp \
     happah/math/TriangleRefinementScheme.cpp \
     happah/readers/ReaderHPH.cpp \
     happah/readers/ReaderOFF.cpp \
     happah/utils/ControlPointIndexer.cpp \
     happah/utils/GeometryUtils.cpp \
     happah/utils/SurfaceUtilsBEZ.cpp \
     happah/weighers/ConstantWeigher.cpp \
     happah/writers/WriterHPH.cpp \
     happah/writers/WriterOFF.cpp
nobase_include_HEADERS = \
     happah/Eigen.h \
     happah/Happah.h \
     happah/geometries/converters.h \
     happah/geometries/Circle.h \
     happah/geometries/Curve.h \
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
     happah/geometries/SplineCurve.h \
     happah/geometries/Strips.h \
     happah/geometries/SurfaceBEZ.h \
     happah/geometries/Surface.h \
     happah/geometries/SurfaceHEZ.h \
     happah/geometries/SurfaceSplineBEZ.h \
     happah/geometries/SurfaceSplineHEZ.h \
     happah/geometries/Trapezoidation.h \
     happah/geometries/Triangle.h \
     happah/geometries/TriangleMesh.h \
     happah/geometries/TriangleMeshUtils.h \
     happah/geometries/VertexCloud.h \
     happah/geometries/Vertex.h \
     happah/math/BasisBuilder.h \
     happah/math/HexagonDecomposition.h \
     happah/math/MathUtils.h \
     happah/math/nlopt.cpp \
     happah/math/PantsDecomposition.h \
     happah/math/ProjectiveStructure.h \
     happah/math/RigidAffineTransformation.h \
     happah/math/Space.h \
     happah/math/TriangleDecomposition.h \
     happah/math/TriangleRefinementScheme.h \
     happah/readers/ReaderHPH.h \
     happah/readers/ReaderOFF.h \
     happah/readers/ReaderUtils.h \
     happah/readers.h \
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
     happah/weighers/TraversableEdgeLengthWeigher.h \
     happah/writers/Writer.h \
     happah/writers/WriterHPH.h \
     happah/writers/WriterOFF.h \
     happah/writers.h
libhappah_la_CPPFLAGS = $(ipopt_CFLAGS) $(eigen3_CFLAGS) $(BOOST_CPPFLAGS)
libhappah_la_CXXFLAGS = -std=c++1y -Wno-unused-label -Wno-unused-parameter -Wno-unused-variable
libhappah_la_LDFLAGS = -version-info 0:0:0 $(ipopt_LDFLAGS) $(eigen3_LDFLAGS) $(BOOST_LDFLAGS)
libhappah_la_LIBADD = $(ipopt_LIBS) $(eigen3_LIBS) -llpsolve55 $(BOOST_IOSTREAMS_LIB)

