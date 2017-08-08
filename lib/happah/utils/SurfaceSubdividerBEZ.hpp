// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/geometries/Triangle.hpp"
#include "happah/utils/SurfaceUtilsBEZ.hpp"

namespace happah {

//TODO: micro-optimizations
//NOTE: Here we subdivide a single surface piece.
template<class Space, hpuint t_degree>
class SurfaceSubdividerBEZ {
     static_assert(t_degree > 0, "Surface subdivision only makes sense for degree greater than zero.");
     using Point = typename Space::POINT;
     using ControlPoints = std::vector<Point>;

public:

     template<class Point>
     static std::pair<std::vector<Point>, std::vector<hpuint> > getParameterPoints(const Point& p0, const Point& p1, const Point& p2, hpuint nSubdivisions) {
          assert(nSubdivisions > 0);

          std::vector<Point> points;
          std::vector<hpuint> indices;

          const hpuint degree = 1 << nSubdivisions;

          points.reserve(make_patch_size(degree));

          sample(degree + 1, [&] (hpreal u, hpreal v, hpreal w) {
               points.push_back(u * p0 + v * p1 + w * p2);
          });

          //NOTE: Make sure when arranging the vertices that their order respects the orientation of the patch in the spline.  In other words, there is only one correct counterclockwise arrangement possible for each patch's parameter triangle.
          //NOTE: The diagram in the notes has p0 in the upper corner and not p2 as usual.  This rotation of the parameter vertices is a side effect of the subdivision algorithm because we have to start the recursion in the case EEB.
          switch(degree) {
          case 2: 
               indices = {
                    0, 1, 3,
                    1, 2, 3,
                    2, 4, 3,
                    4, 5, 3
               };
               break;
          case 4:
               indices = {
                    0, 1, 5,
                    1, 2, 5,
                    2, 6, 5,
                    6, 9, 5,
                    2, 3, 7,
                    6, 2, 7,
                    9, 6, 7,
                    10, 9, 7,
                    9, 10, 12,
                    3, 4, 7,
                    4, 8, 7,
                    8, 11, 7,
                    11, 10, 7,
                    10, 11, 12,
                    11, 13, 12,
                    13, 14, 12
               };
               break;
          case 8:
               indices = {
                    0, 1, 9,
                    1, 2, 9,
                    2, 10, 9,
                    10, 17, 9,
                    2, 3, 11,
                    10, 2, 11,
                    17, 10, 11,
                    18, 17, 11,
                    17, 18, 24,
                    3, 4, 11,
                    4, 12, 11,
                    12, 19, 11,
                    19, 18, 11,
                    18, 19, 24,
                    19, 25, 24,
                    25, 30, 24,
                    4, 5, 13,
                    12, 4, 13,
                    19, 12, 13,
                    20, 19, 13,
                    19, 20, 26,
                    25, 19, 26,
                    30, 25, 26,
                    31, 30, 26,
                    30, 31, 35,
                    5, 6, 13,
                    6, 14, 13,
                    14, 21, 13,
                    21, 20, 13,
                    20, 21, 26,
                    21, 27, 26,
                    27, 32, 26,
                    32, 31, 26,
                    31, 32, 35,
                    32, 36, 35,
                    36, 39, 35,
                    6, 7, 15,
                    14, 6, 15,
                    21, 14, 15,
                    22, 21, 15,
                    21, 22, 28,
                    27, 21, 28,
                    32, 27, 28,
                    33, 32, 28,
                    32, 33, 37,
                    36, 32, 37,
                    39, 36, 37,
                    40, 39, 37,
                    39, 40, 42,
                    7, 8, 15,
                    8, 16, 15,
                    16, 23, 15,
                    23, 22, 15,
                    22, 23, 28,
                    23, 29, 28,
                    29, 34, 28,
                    34, 33, 28,
                    33, 34, 37,
                    34, 38, 37,
                    38, 41, 37,
                    41, 40, 37,
                    40, 41, 42,
                    41, 43, 42,
                    43, 44, 42
               };
               break;
          default: {
                    indices.reserve(3 * make_control_polygon_size(degree));

                    hpuint temp[] = { 0, 1, degree + 1, 1, 2, degree + 1, 2, degree + 2, degree + 1, degree + 2, (degree << 1) + 1, degree + 1 };
                    indices.insert(indices.end(), temp, temp + 12);

                    hpuint r = 2;//row
                    while(r < degree) {
                         hpuint c, vt, vb, o;

                         //even row
                         c = 0;
                         vt = r;
                         vb = vt + 1;
                         o = degree;
                         while(c < r) {
                              hpuint vnt = vt + o;
                              hpuint vnb = vb + o;
                              --o;

                              hpuint vs1[] = { vt, vb, vnb, vnt, vt, vnb };
                              indices.insert(indices.end(), vs1, vs1 + 6);

                              vt = vnt;
                              vb = vnb;

                              vnt += o;
                              vnb += o;
                              --o;

                              hpuint vs2[] = { vnt, vt, vb, vnb, vnt, vb };
                              indices.insert(indices.end(), vs2, vs2 + 6);

                              vt = vnt;
                              vb = vnb;

                              c += 2;
                         }

                         indices.push_back(vt);
                         indices.push_back(vb);
                         indices.push_back(vb + o);

                         //odd row
                         c = 0;
                         vt = r + 1;
                         vb = vt + 1;
                         o = degree;
                         while(c < r) {
                              hpuint vnt = vt + o;
                              hpuint vnb = vb + o;
                              --o;

                              hpuint vs1[] = { vt, vb, vnt, vb, vnb, vnt };
                              indices.insert(indices.end(), vs1, vs1 + 6);

                              vt = vnt;
                              vb = vnb;

                              vnt += o;
                              vnb += o;
                              --o;

                              hpuint vs2[] = { vb, vnb, vt, vnb, vnt, vt };
                              indices.insert(indices.end(), vs2, vs2 + 6);

                              vt = vnt;
                              vb = vnb;

                              c += 2;
                         }

                         hpuint vnt = vt + o;
                         hpuint vnb = vb + o;
                         --o;

                         hpuint vs[] = { vt, vb, vnt, vb, vnb, vnt, vnb, vnb + o, vnt };
                         indices.insert(indices.end(), vs, vs + 9);

                         r += 2;
                    }
               }
          };

          return std::make_pair(std::move(points), std::move(indices));
     }

     //NOTE: The control points are ordered left to right starting at the bottom row and ending at the top row, which contains one point.
     template<class Iterator>
     SurfaceSubdividerBEZ(Iterator controlPoints) {
          //TODO: specialize for degree 2, 3
          //organize the input control points into the arrays as specified by the algorithm
          m_oFacePoints.resize(FACE_STRIDE);
          m_oEdgePoints.resize(3 * EDGE_STRIDE);
          m_oVertexPoints.resize(3);

          auto v = m_oVertexPoints.begin();
          auto i = controlPoints;
          auto f = m_oFacePoints.begin();
          auto e0 = m_oEdgePoints.begin();
          auto e1 = e0 + EDGE_STRIDE;
          auto e2 = e1 + EDGE_STRIDE;

          *v = *i;
          ++i;
          for(auto end = i + EDGE_STRIDE; i != end; ++i, ++e0) *e0 = *i;
          *(++v) = *i;

          auto nMiddlePoints = t_degree - 1;
          while(nMiddlePoints > 0) {
               --nMiddlePoints;
               *e1 = *(++i);
               ++e1;
               auto end = f;
               f += nMiddlePoints;
               auto f0 = f;
               while(f0 != end) *(--f0) = *(++i);
               *e2 = *(++i);
               ++e2;
          }

          *(++v) = *(++i);
     }

     std::tuple<ControlPoints, Indices> subdivide(hpuint nSubdivisions) {
          assert(nSubdivisions > 0);
          //TODO: specialize for degree 2, 3

          //calculate number of points of each type
          auto nRows = 1u << nSubdivisions;//NOTE: A row consists of zero or more rectangles followed by a triangle.
          auto nRectangles = (nRows * (nRows - 1)) >> 1;//NOTE: The number of rectangles is the result of the arithmetic series for (nRows-1).  The first row has no rectangles and each following row has one rectangle more than the previous row.
          auto nTriangles = (nRectangles << 1) + nRows;
          auto nEdges = 3 * (nRectangles + nRows);//NOTE: We count the three bottom left edges in each triangle so as not to count the top edge and the right edge are not counted twice.  Then, three edges per triangle at the end of each row remain to be counted.
          auto nFacePoints = nTriangles * FACE_STRIDE;//NOTE: This is the number of points not on the border.
          auto nEdgePoints = nEdges * EDGE_STRIDE;
          auto nVertexPoints = nRectangles + (nRows << 1) + 1;

          auto facePointsA = ControlPoints(nFacePoints);
          auto facePointsB = ControlPoints(nFacePoints);
          auto facePointsC = ControlPoints(FACE_STRIDE << 1);//stores two temporary faces generated by first binary subdivision
          auto edgePointsA = ControlPoints(nEdgePoints);
          auto edgePointsB = ControlPoints(nEdgePoints);
          auto vertexPointsA = ControlPoints(nVertexPoints);
          auto vertexPointsB = ControlPoints(nVertexPoints);

          auto uf0 = facePointsC.begin();
          auto uf1 = uf0 + FACE_STRIDE;

          {
               auto svt = m_oVertexPoints.begin();//top source vertices
               auto svb = svt + 1;//bottom source vertices
               auto sev = m_oEdgePoints.begin();//non-horizontal source edges
               auto seh = sev + (EDGE_STRIDE << 1);//horizontal source edges
               auto sf = m_oFacePoints.begin();//source face

               auto tvt = vertexPointsA.begin();//top target vertices
               auto tvm = tvt + 1;//middle target vertices
               auto tvb = tvm + 2;//bottom target vertices
               auto tevt = edgePointsA.begin();//non-horizontal top target edges
               auto tevb = tevt + (EDGE_STRIDE << 1);//non-horizontal bottom target edges
               auto teht = tevt + (6 * EDGE_STRIDE);//horizontal top target edges
               auto tehb = teht + EDGE_STRIDE;//horizontal bottom target edges
               auto tft = facePointsA.begin();//top target faces
               auto tfb = tft + FACE_STRIDE;//bottom target faces

               subdivideEEB(svt, svb, sev, seh, sf, tvt, tvm, tvb, tevt, tevb, teht, tehb, tft, tfb, uf0, uf1);

               auto iv = m_oVertexPoints.begin();
               auto ov = vertexPointsA.begin();
               *ov = *iv;
               ++iv;
               *(ov + 3) = *iv;
               ++iv;
               *(ov + 5) = *iv;
          }

          auto sFacePoints = &facePointsA;//source points
          auto sEdgePoints = &edgePointsA;
          auto sVertexPoints = &vertexPointsA;
          auto tFacePoints = &facePointsB;//target points
          auto tEdgePoints = &edgePointsB;
          auto tVertexPoints = &vertexPointsB;
          for(auto i = 1u; i < nSubdivisions; ++i) {
               auto snRows = 1u << i;
               auto snNonHorizontalEdges = snRows * (snRows + 1);
               auto tnRows = 1u << (i + 1);
               auto tnNonHorizontalEdges = tnRows * (tnRows + 1);

               auto svt = sVertexPoints->begin();//same meaning as above
               auto svb = svt + 1;
               auto sev = sEdgePoints->begin();
               auto sehb = sev + (snNonHorizontalEdges * EDGE_STRIDE);
               auto seht = sehb;
               auto sf = sFacePoints->begin();

               auto tvt = tVertexPoints->begin();//same meaning as above
               auto tvm = tvt + 1;
               auto tvb = tvm + 2;
               auto tevt = tEdgePoints->begin();
               auto tevb = tevt + (EDGE_STRIDE << 1);
               auto tehm = tEdgePoints->begin() + (tnNonHorizontalEdges * EDGE_STRIDE);
               auto tehb = tehm + EDGE_STRIDE;
               auto teht = tehb;
               auto tft = tFacePoints->begin();
               auto tfb = tft + FACE_STRIDE;

               auto r = 0u;
               while(true) {
                    hpuint c;

                    c = 0;
                    while(c < r) {//subdivide triangles in even row
                         subdivideEEB<true>(svt, svb, sev, sehb, sf, tvt, tvm, tvb, tevt, tevb, tehm, tehb, tft, tfb, uf0, uf1, c > 0);
                         *tvt = *svt;

                         ++svb;
                         sev += EDGE_STRIDE;
                         sehb += EDGE_STRIDE;
                         sf += FACE_STRIDE;

                         ++tvm;
                         tvb += 2;
                         tevt += EDGE_STRIDE;
                         tevb += (3 * EDGE_STRIDE);
                         tehm += EDGE_STRIDE;
                         tehb += (EDGE_STRIDE << 1);
                         tft += FACE_STRIDE;
                         tfb += (3 * FACE_STRIDE);

                         subdivideEET(svt, svb, sev, seht, sf, tvt, tvm, tvb, tevt, tevb, teht, tehm, tft, tfb, uf0, uf1);

                         ++svt;
                         sev += EDGE_STRIDE;
                         seht += EDGE_STRIDE;
                         sf += FACE_STRIDE;

                         tvt += 2;
                         ++tvm;
                         tevt += (3 * EDGE_STRIDE);
                         tevb += EDGE_STRIDE;
                         teht += (EDGE_STRIDE << 1);
                         tehm += EDGE_STRIDE;
                         tft += (3 * FACE_STRIDE);
                         tfb += FACE_STRIDE;

                         ++c;

                         subdivideEOT(svt, svb, sev, seht, sf, tvt, tvm, tvb, tevt, tevb, teht, tehm, tft, tfb, uf0, uf1);
                         *tvt = *svt;

                         ++svt;
                         sev += EDGE_STRIDE;
                         seht += EDGE_STRIDE;
                         sf += FACE_STRIDE;

                         tvt += 2;
                         ++tvm;
                         tevt += (3 * EDGE_STRIDE);
                         tevb += EDGE_STRIDE;
                         teht += (EDGE_STRIDE << 1);
                         tehm += EDGE_STRIDE;
                         tft += (3 * FACE_STRIDE);
                         tfb += FACE_STRIDE;

                         subdivideEOB(svt, svb, sev, sehb, sf, tvt, tvm, tvb, tevt, tevb, tehm, tehb, tft, tfb, uf0, uf1);

                         ++svb;
                         sev += EDGE_STRIDE;
                         sehb += EDGE_STRIDE;
                         sf += FACE_STRIDE;

                         ++tvm;
                         tvb += 2;
                         tevt += EDGE_STRIDE;
                         tevb += (3 * EDGE_STRIDE);
                         tehm += EDGE_STRIDE;
                         tehb += (EDGE_STRIDE << 1);
                         tft += FACE_STRIDE;
                         tfb += (3 * FACE_STRIDE);

                         ++c;
                    }

                    //subdivide last triangle in row
                    subdivideEEB<true>(svt, svb, sev, sehb, sf, tvt, tvm, tvb, tevt, tevb, tehm, tehb, tft, tfb, uf0, uf1, c > 0);
                    *tvt = *svt;

                    ++svt;
                    svb += 2;
                    sev += (EDGE_STRIDE << 1);
                    sehb += EDGE_STRIDE;
                    sf += FACE_STRIDE;

                    tvt = tvm + 2;
                    tvm = tvb + 3;
                    tevt = tevb + (EDGE_STRIDE << 2);
                    tevb = tevt + (((r << 2) + 6) * EDGE_STRIDE);
                    teht = tehm + EDGE_STRIDE;
                    tehm = tehb + (EDGE_STRIDE << 1);
                    ++r;
                    tvb = tvm + ((r << 1) + 2);
                    tehb = tehm + (((r << 1) + 1) * EDGE_STRIDE);
                    tft = tfb + (3 * FACE_STRIDE);
                    tfb = tft + (((r << 2) + 1) * FACE_STRIDE);

                    c = 0;
                    bool notLastRow = (r+1) != snRows;
                    while(c < r) {//subdivide triangles in odd row; subdivide bottom horizontal edge if this is the last row
                         subdivideOET(svt, svb, sev, seht, sf, tvt, tvm, tvb, tevt, tevb, teht, tehm, tft, tfb, uf0, uf1, c > 0);
                         *tvt = *svt;

                         ++svt;
                         sev += EDGE_STRIDE;
                         seht += EDGE_STRIDE;
                         sf += FACE_STRIDE;

                         tvt += 2;
                         ++tvm;
                         tevt += (3 * EDGE_STRIDE);
                         tevb += EDGE_STRIDE;
                         teht += (EDGE_STRIDE << 1);
                         tehm += EDGE_STRIDE;
                         tft += (3 * FACE_STRIDE);
                         tfb += FACE_STRIDE;

                         subdivideOEB(svt, svb, sev, sehb, sf, tvt, tvm, tvb, tevt, tevb, tehm, tehb, tft, tfb, uf0, uf1, notLastRow);

                         ++svb;
                         sev += EDGE_STRIDE;
                         sehb += EDGE_STRIDE;
                         sf += FACE_STRIDE;

                         ++tvm;
                         tvb += 2;
                         tevt += EDGE_STRIDE;
                         tevb += (3 * EDGE_STRIDE);
                         tehm += EDGE_STRIDE;
                         tehb += (EDGE_STRIDE << 1);
                         tft += FACE_STRIDE;
                         tfb += (3 * FACE_STRIDE);

                         ++c;
                         if(c == r) break;

                         subdivideOOB(svt, svb, sev, sehb, sf, tvt, tvm, tvb, tevt, tevb, tehm, tehb, tft, tfb, uf0, uf1, notLastRow);

                         ++svb;
                         sev += EDGE_STRIDE;
                         sehb += EDGE_STRIDE;
                         sf += FACE_STRIDE;

                         ++tvm;
                         tvb += 2;
                         tevt += EDGE_STRIDE;
                         tevb += (3 * EDGE_STRIDE);
                         tehm += EDGE_STRIDE;
                         tehb += (EDGE_STRIDE << 1);
                         tft += FACE_STRIDE;
                         tfb += (3 * FACE_STRIDE);

                         subdivideOOT(svt, svb, sev, seht, sf, tvt, tvm, tvb, tevt, tevb, teht, tehm, tft, tfb, uf0, uf1);
                         *tvt = *svt;

                         ++svt;
                         sev += EDGE_STRIDE;
                         seht += EDGE_STRIDE;
                         sf += FACE_STRIDE;

                         tvt += 2;
                         ++tvm;
                         tevt += (3 * EDGE_STRIDE);
                         tevb += EDGE_STRIDE;
                         teht += (EDGE_STRIDE << 1);
                         tehm += EDGE_STRIDE;
                         tft += (3 * FACE_STRIDE);
                         tfb += FACE_STRIDE;

                         ++c;
                    }

                    //subdivide last triangle in row
                    subdivideOOB(svt, svb, sev, sehb, sf, tvt, tvm, tvb, tevt, tevb, tehm, tehb, tft, tfb, uf0, uf1, notLastRow);
                    *tvt = *svt;

                    ++svt;
                    svb += 2;
                    sev += (EDGE_STRIDE << 1);
                    sehb += EDGE_STRIDE;
                    sf += FACE_STRIDE;

                    tvt = tvm + 2;
                    tvm = tvb + 3;
                    tevt = tevb + (EDGE_STRIDE << 2);
                    tevb = tevt + (((r << 2) + 6) * EDGE_STRIDE);
                    teht = tehm + EDGE_STRIDE;
                    tehm = tehb + (EDGE_STRIDE << 1);
                    ++r;
                    if(r == snRows) break;
                    tvb = tvm + ((r << 1) + 2);
                    tehb = tehm + (((r << 1) + 1) * EDGE_STRIDE);
                    tft = tfb + (3 * FACE_STRIDE);
                    tfb = tft + (((r << 2) + 1) * FACE_STRIDE);
               }

               //add common control points at triangle vertices in bottom row to target
               for(hpuint j = 0, end = snRows + 1; j < end; ++j) {
                    *tvt = *svt;
                    ++svt;
                    tvt += 2;
               }
               
               if((i & 1) == 0) {
                    sFacePoints = &facePointsA;//source points
                    sEdgePoints = &edgePointsA;
                    sVertexPoints = &vertexPointsA;
                    tFacePoints = &facePointsB;//target points
                    tEdgePoints = &edgePointsB;
                    tVertexPoints = &vertexPointsB;
               } else {
                    sFacePoints = &facePointsB;//source points
                    sEdgePoints = &edgePointsB;
                    sVertexPoints = &vertexPointsB;
                    tFacePoints = &facePointsA;//target points
                    tEdgePoints = &edgePointsA;
                    tVertexPoints = &vertexPointsA;
               } 
          }

          ControlPoints points;

          if((nSubdivisions & 1) == 0) {
               points = std::move(facePointsB);
               points.reserve(nFacePoints + nEdgePoints + nVertexPoints);
               points.insert(points.end(), edgePointsB.begin(), edgePointsB.end());
               points.insert(points.end(), vertexPointsB.begin(), vertexPointsB.end());
          } else {
               points = std::move(facePointsA);
               points.reserve(nFacePoints + nEdgePoints + nVertexPoints);
               points.insert(points.end(), edgePointsA.begin(), edgePointsA.end());
               points.insert(points.end(), vertexPointsA.begin(), vertexPointsA.end());
          }

          Indices indices;
          indices.reserve(nTriangles * make_patch_size(t_degree));

          auto vt = nFacePoints + nEdgePoints;
          auto vb = vt + 1;
          auto ev = nFacePoints;
          auto ehb = ev + nRows * (nRows + 1) * EDGE_STRIDE;//top horizontal edge
          auto eht = ehb;//bottom horizontal edge
          auto f = 0u;
          hpuint limit, nMiddlePoints;
          hpuint tev;//holds index of ev temporarily
          auto r = 0u;
          while(r < nRows) {
               hpuint c;

               c = 0;
               while(c < r) {
                    //eeb
                    indices.push_back(vt);
                    limit = ev + EDGE_STRIDE;
                    while(ev < limit) indices.push_back(ev++);
                    indices.push_back(vb++);
                    tev = ev;
                    nMiddlePoints = t_degree - 2;
                    while(nMiddlePoints > 0) {
                         indices.push_back(tev++);
                         limit = f;
                         f += nMiddlePoints;
                         while(f > limit) indices.push_back(--f);
                         f += nMiddlePoints;
                         indices.push_back(ehb++);
                         --nMiddlePoints;
                    }
                    indices.push_back(tev++);
                    indices.push_back(ehb++);
                    indices.push_back(vb);

                    //eet
                    indices.push_back(vt + 1);
                    limit = eht;
                    eht += EDGE_STRIDE;
                    while(eht > limit) indices.push_back(--eht);
                    eht += EDGE_STRIDE;
                    indices.push_back(vt++);
                    tev = ev + EDGE_STRIDE;
                    nMiddlePoints = t_degree - 2;
                    while(nMiddlePoints > 0) {
                    indices.push_back(tev++);
                         limit = f + nMiddlePoints;
                         while(f < limit) indices.push_back(f++);
                         indices.push_back(ev++);
                         --nMiddlePoints;
                    }
                    indices.push_back(tev++);
                    indices.push_back(ev++);
                    indices.push_back(vb);

                    ++c;

                    //eot
                    indices.push_back(vt + 1);
                    limit = eht;
                    eht += EDGE_STRIDE;
                    while(eht > limit) indices.push_back(--eht);
                    eht += EDGE_STRIDE;
                    indices.push_back(vt++);
                    tev = ev + (EDGE_STRIDE << 1);
                    nMiddlePoints = t_degree - 2;
                    while(nMiddlePoints > 0) {
                         indices.push_back(--tev);
                         limit = f;
                         f += nMiddlePoints;
                         while(f > limit) indices.push_back(--f);
                         f += nMiddlePoints;
                         indices.push_back(ev++);
                         --nMiddlePoints;
                    }
                    indices.push_back(--tev);
                    indices.push_back(ev++);
                    indices.push_back(vb);

                    //eob
                    indices.push_back(vb + 1);
                    tev = ev + (EDGE_STRIDE << 1);
                    ev += EDGE_STRIDE;
                    while(tev > ev) indices.push_back(--tev);
                    indices.push_back(vt);
                    ehb += EDGE_STRIDE;
                    nMiddlePoints = t_degree - 2;
                    while(nMiddlePoints > 0) {
                         indices.push_back(--ehb);
                         limit = f + nMiddlePoints;
                         while(f < limit) indices.push_back(f++);
                         indices.push_back(--tev);
                         --nMiddlePoints;
                    }
                    indices.push_back(--ehb);
                    indices.push_back(--tev);
                    indices.push_back(vb++);
                    ehb += EDGE_STRIDE;

                    ++c;
               }

               //eeb
               indices.push_back(vt++);//lower left point
               limit = ev + EDGE_STRIDE;//first row
               while(ev < limit) indices.push_back(ev++);
               indices.push_back(vb++);
               nMiddlePoints = t_degree - 2;//middle rows
               while(nMiddlePoints > 0) {
                    indices.push_back(ev++);
                    limit = f;
                    f += nMiddlePoints;
                    while(f > limit) indices.push_back(--f);
                    f += nMiddlePoints;
                    indices.push_back(ehb++);
                    --nMiddlePoints;
               }
               indices.push_back(ev++);//last row
               indices.push_back(ehb++);
               indices.push_back(vb++);//top point
               
               ++r;

               c = 0;
               while(c < r) {
                    //oet
                    indices.push_back(vt++);
                    limit = ev + EDGE_STRIDE;
                    while(ev < limit) indices.push_back(ev++);
                    indices.push_back(vb);
                    tev = ev;
                    nMiddlePoints = t_degree - 2;
                    while(nMiddlePoints > 0) {
                         indices.push_back(eht++);
                         limit = f + nMiddlePoints;
                         while(f < limit) indices.push_back(f++);
                         indices.push_back(tev++);
                         --nMiddlePoints;
                    }
                    indices.push_back(eht++);
                    indices.push_back(tev++);
                    indices.push_back(vt);

                    //oeb
                    indices.push_back(vb++);
                    limit = ehb + EDGE_STRIDE;
                    while(ehb < limit) indices.push_back(ehb++);
                    indices.push_back(vb);
                    tev = ev + (EDGE_STRIDE << 1);
                    nMiddlePoints = t_degree - 2;
                    while(nMiddlePoints > 0) {
                         indices.push_back(ev++);
                         limit = f;
                         f += nMiddlePoints;
                         while(f > limit) indices.push_back(--f);
                         f += nMiddlePoints;
                         indices.push_back(--tev);
                         --nMiddlePoints;
                    }
                    indices.push_back(ev++);
                    indices.push_back(--tev);
                    indices.push_back(vt);

                    ++c;

                    if(c == r) break;

                    //oob
                    indices.push_back(vb++);
                    limit = ehb + EDGE_STRIDE;
                    while(ehb < limit) indices.push_back(ehb++);
                    indices.push_back(vb);
                    tev = ev + (EDGE_STRIDE << 1);
                    ev += EDGE_STRIDE;
                    nMiddlePoints = t_degree - 2;
                    while(nMiddlePoints > 0) {
                         indices.push_back(--ev);
                         limit = f + nMiddlePoints;
                         while(f < limit) indices.push_back(f++);
                         indices.push_back(--tev);
                         --nMiddlePoints;
                    }
                    indices.push_back(--ev);
                    indices.push_back(--tev);
                    indices.push_back(vt);
                    ev += EDGE_STRIDE;

                    //oot
                    indices.push_back(vb);
                    tev = ev + (EDGE_STRIDE << 1);
                    ev += EDGE_STRIDE;
                    while(tev > ev) indices.push_back(--tev);
                    indices.push_back(vt + 1);
                    nMiddlePoints = t_degree - 2;
                    eht += EDGE_STRIDE;
                    while(nMiddlePoints > 0) {
                         indices.push_back(--ev);
                         limit = f;
                         f += nMiddlePoints;
                         while(f > limit) indices.push_back(--f);
                         f += nMiddlePoints;
                         indices.push_back(--eht);
                         --nMiddlePoints;
                    }
                    indices.push_back(--ev);
                    indices.push_back(--eht);
                    indices.push_back(vt++);
                    ev += EDGE_STRIDE;
                    eht += EDGE_STRIDE;

                    ++c;
               }

               //oob
               indices.push_back(vb++);
               limit = ehb + EDGE_STRIDE;
               while(ehb < limit) indices.push_back(ehb++);
               indices.push_back(vb++);
               tev = ev + (EDGE_STRIDE << 1);
               ev += EDGE_STRIDE;
               nMiddlePoints = t_degree - 2;
               while(nMiddlePoints > 0) {
                    indices.push_back(--ev);
                    limit = f + nMiddlePoints;
                    while(f < limit) indices.push_back(f++);
                    indices.push_back(--tev);
                    --nMiddlePoints;
               }
               indices.push_back(--ev);
               indices.push_back(--tev);
               indices.push_back(vt++);
               ev = tev + EDGE_STRIDE;
               
               ++r;
          }

          return std::make_tuple(std::move(points), std::move(indices));
     }

private:
     using PointIterator = typename ControlPoints::iterator;
     using r = std::reverse_iterator<PointIterator>;

     static constexpr hpuint EDGE_STRIDE = t_degree - 1;
     static constexpr hpuint FACE_STRIDE = make_patch_size(t_degree) - 3 * t_degree;

     ControlPoints m_oEdgePoints;
     ControlPoints m_oFacePoints;
     ControlPoints m_oVertexPoints;

     //NOTE: Splits triangle into two pieces.
     template<
          bool t_skipEdge = false,
          class PointIteratorIV1,
          class PointIteratorIV2,
          class PointIteratorOV01,
          class PointIteratorIE0,
          class PointIteratorIE1,
          class PointIteratorIE2,
          class PointIteratorOE00,
          class PointIteratorOE01,
          class PointIteratorOE10,
          class PointIteratorIF0,
          class PointIteratorOF0,
          class PointIteratorOF1
     >
     void doBinaryTriangleSubdivisionStep(
          PointIteratorIV1 iv1,
          PointIteratorIV2 iv2,
          PointIteratorOV01 ov01,
          PointIteratorIE0 ie0,
          PointIteratorIE1 ie1,
          PointIteratorIE2 ie2,
          PointIteratorOE00 oe00,
          PointIteratorOE01 oe01,
          PointIteratorOE10 oe10,
          PointIteratorIF0 if0,
          PointIteratorOF0 of0,
          PointIteratorOF1 of1,
          bool skipEdge = false
     ) {
          if(!skipEdge || !t_skipEdge) doBinaryEdgeSubdivisionStep(iv1, iv2, ov01, ie1, oe00, oe10);

          hpuint nRows = t_degree - 2;
          if(nRows > 0) {
               std::vector<Point> middlePoints;
               middlePoints.reserve(t_degree - 3);
               while(nRows > 0) {
                    middlePoints.clear();
                    --nRows;
                    of0 += nRows;
                    of1 += nRows;
                    auto tf0 = of0;
                    auto tf1 = of1;
                    auto g = if0 + nRows;
                    *tf0 = (*ie0 + *g) * hpreal(0.5);
                    ++ie0;
                    hpuint limit = nRows;
                    auto rowLength = t_degree - 2;
                    auto f = if0 + rowLength;
                    while(limit > 0) {
                         --rowLength;//TODO: ewwww
                         --limit;
                         auto h = f + limit;
                         middlePoints.push_back((*h + *g) * hpreal(0.5));
                         f += rowLength;
                         g = h;
                    }
                    *tf1 = (*g + *ie2) * hpreal(0.5);
                    ++ie2;

                    hpuint nMiddlePoints = nRows;
                    while(nMiddlePoints > 0) {
                         auto m = middlePoints.begin();
                         auto t0 = tf0 - 1;
                         *t0 = (*tf0 + *m) * hpreal(0.5);
                         tf0 = t0;
                         hpuint limit = --nMiddlePoints;
                         while(limit > 0) {
                              auto a = m + 1;
                              *m = (*a + *m) * hpreal(0.5);
                              m = a;
                              --limit;
                         }
                         auto t1 = tf1 - 1;
                         *t1 = (*m + *tf1) * hpreal(0.5);
                         tf1 = t1;
                    }
                    *oe01 = (*tf0 + *tf1) * hpreal(0.5);
                    ++oe01;
                    
                    ++of0;
                    ++of1;
               }
          }
          *oe01 = (*ie0 + *ie2) * hpreal(0.5);
     }

     //Splits edge into two pieces.
     template<
          class PointIteratorIV0,
          class PointIteratorIV1,
          class PointIteratorOV1,
          class PointIteratorIE0,
          class PointIteratorOE0,
          class PointIteratorOE1
     >
     void doBinaryEdgeSubdivisionStep(
          PointIteratorIV0 iv0,
          PointIteratorIV1 iv1,
          PointIteratorOV1 ov1,
          PointIteratorIE0 ie0,
          PointIteratorOE0 oe0,
          PointIteratorOE1 oe1
     ) {
          hpuint nMiddlePoints = t_degree - 2;
          ControlPoints middlePoints;
          middlePoints.reserve(nMiddlePoints);

          *oe0 = (*iv0 + *ie0) * hpreal(0.5);//NOTE: This assumes that ie contains at least one point, which means degree must be greater than 1.
          hpuint limit = nMiddlePoints;//calculate first column of the de castlejau algorithm
          while(limit > 0) {
               auto a = ie0;
               middlePoints.push_back((*a + *(++ie0)) * hpreal(0.5));
               --limit;
          }
          *oe1 = (*ie0 + *iv1) * hpreal(0.5);

          while(nMiddlePoints > 0) {//calculate the rest of the columns of the de casteljau algorithm
               auto m = middlePoints.begin();
               auto t0 = oe0 + 1;
               *t0 = (*oe0 + *m) * hpreal(0.5);
               oe0 = t0;
               hpuint limit = --nMiddlePoints;
               while(limit > 0) {
                    auto a = m + 1;
                    *m = (*a + *m) * hpreal(0.5);
                    m = a;
                    --limit;
               }
               auto t1 = oe1 + 1;
               *t1 = (*m + *oe1) * hpreal(0.5);
               oe1 = t1;
          }
          *ov1 = (*oe0 + *oe1) * hpreal(0.5);
     }

     //even row, even column, bottom triangle
     template<bool t_skip1 = false>
     void subdivideEEB(
          PointIterator svt, 
          PointIterator svb, 
          PointIterator sev, 
          PointIterator seh, 
          PointIterator sf, 

          PointIterator tvt, 
          PointIterator tvm, 
          PointIterator tvb, 
          PointIterator tevt, 
          PointIterator tevb, 
          PointIterator teht, 
          PointIterator tehb, 
          PointIterator tft, 
          PointIterator tfb,

          PointIterator uf0,
          PointIterator uf1,

          bool skip2 = false
     ) {
          hpuint e = EDGE_STRIDE;

          auto iv2 = svt;
          auto iv1 = svb;
          auto iv0 = iv1 + 1;

          auto ie2 = sev;
          auto ie0 = ie2 + e;
          auto ie1 = seh;

          auto if0 = sf;

          auto ov2 = tvm;
          auto ov0 = ov2 + 1;
          auto ov1 = tvb + 1;

          auto oe7 = tevt;
          auto oe2 = oe7 + e;
          auto oe8 = tevb;
          auto oe4 = oe8 + e;
          auto oe6 = oe4 + e;
          auto oe0 = oe6 + e;
          auto oe3 = teht;
          auto oe5 = tehb;
          auto oe1 = oe5 + e;

          auto of1 = tft;
          auto of2 = tfb;
          auto of3 = of2 + FACE_STRIDE;
          auto of0 = of3 + FACE_STRIDE;

          doBinaryTriangleSubdivisionStep(iv2, iv0, ov0, ie2, ie0, r(ie1 + e), oe2, r(oe4 + e), r(oe0 + e), if0, uf0, uf1);
          doBinaryTriangleSubdivisionStep(iv2, iv1, ov2, oe2, ie2, oe4, oe7, oe3, r(oe8 + e), uf0, of1, of2, skip2);
          doBinaryTriangleSubdivisionStep<t_skip1>(iv0, iv1, ov1, r(oe0 + e), r(ie1 + e), oe4, r(oe1 + e), r(oe6 + e), oe5, uf1, of0, of3);
     }

     //even row, even column, top triangle
     void subdivideEET(
          PointIterator svt, 
          PointIterator svb, 
          PointIterator sev, 
          PointIterator seh, 
          PointIterator sf, 

          PointIterator tvt, 
          PointIterator tvm, 
          PointIterator tvb, 
          PointIterator tevt, 
          PointIterator tevb, 
          PointIterator teht, 
          PointIterator tehb, 
          PointIterator tft, 
          PointIterator tfb,

          PointIterator uf0,
          PointIterator uf1
     ) {
          hpuint e = EDGE_STRIDE;

          auto iv0 = svb;
          auto iv2 = svt;
          auto iv1 = iv2 + 1;

          auto ie0 = sev;
          auto ie1 = ie0 + e;
          auto ie2 = seh;

          auto if0 = sf;

          auto ov0 = tvm;
          auto ov1 = ov0 + 1;
          auto ov2 = tvt + 1;

          auto oe2 = tevt;
          auto oe3 = oe2 + e;
          auto oe4 = oe3 + e;
          auto oe5 = oe4 + e;
          auto oe0 = tevb;
          auto oe1 = oe0 + e;
          auto oe7 = teht;
          auto oe8 = oe7 + e;
          auto oe6 = tehb;

          auto of1 = tft;
          auto of2 = of1 + FACE_STRIDE;
          auto of3 = of2 + FACE_STRIDE;
          auto of0 = tfb;

          doBinaryTriangleSubdivisionStep<true>(iv2, iv0, ov0, ie2, ie0, r(ie1 + e), oe2, oe4, r(oe0 + e), if0, uf0, uf1);
          doBinaryTriangleSubdivisionStep(iv2, iv1, ov2, oe2, ie2, r(oe4 + e), oe7, oe3, r(oe8 + e), uf0, of1, of2);
          doBinaryTriangleSubdivisionStep(iv0, iv1, ov1, r(oe0 + e), r(ie1 + e), r(oe4 + e), r(oe1 + e), r(oe6 + e), oe5, uf1, of0, of3);
     }

     //even row, odd column, bottom triangle
     void subdivideEOB(
          PointIterator svt, 
          PointIterator svb, 
          PointIterator sev, 
          PointIterator seh, 
          PointIterator sf, 

          PointIterator tvt, 
          PointIterator tvm, 
          PointIterator tvb, 
          PointIterator tevt, 
          PointIterator tevb, 
          PointIterator teht, 
          PointIterator tehb, 
          PointIterator tft, 
          PointIterator tfb,

          PointIterator uf0,
          PointIterator uf1
     ) {
          hpuint e = EDGE_STRIDE;

          auto iv0 = svb;
          auto iv1 = iv0 + 1;
          auto iv2 = svt;

          auto ie0 = sev;
          auto ie2 = ie0 + e;
          auto ie1 = seh;

          auto if0 = sf;

          auto ov0 = tvm;
          auto ov2 = ov0 + 1;
          auto ov1 = tvb + 1;

          auto oe2 = tevt;
          auto oe7 = oe2 + e;
          auto oe0 = tevb;
          auto oe6 = oe0 + e;
          auto oe4 = oe6 + e;
          auto oe8 = oe4 + e;
          auto oe3 = teht;
          auto oe1 = tehb;
          auto oe5 = oe1 + e;

          auto of1 = tft;
          auto of0 = tfb;
          auto of3 = of0 + FACE_STRIDE;
          auto of2 = of3 + FACE_STRIDE;

          doBinaryTriangleSubdivisionStep<true>(iv2, iv0, ov0, ie2, r(ie0 + e), ie1, r(oe2 + e), oe4, oe0, if0, uf0, uf1);
          doBinaryTriangleSubdivisionStep(iv2, iv1, ov2, r(oe2 + e), ie2, r(oe4 + e), oe7, r(oe3 + e), r(oe8 + e), uf0, of1, of2);
          doBinaryTriangleSubdivisionStep<true>(iv0, iv1, ov1, oe0, ie1, r(oe4 + e), oe1, r(oe6 + e), r(oe5 + e), uf1, of0, of3);
     }

     //even row, odd column, top triangle
     void subdivideEOT(
          PointIterator svt, 
          PointIterator svb, 
          PointIterator sev, 
          PointIterator seh, 
          PointIterator sf, 

          PointIterator tvt, 
          PointIterator tvm, 
          PointIterator tvb, 
          PointIterator tevt, 
          PointIterator tevb, 
          PointIterator teht, 
          PointIterator tehb, 
          PointIterator tft, 
          PointIterator tfb,

          PointIterator uf0,
          PointIterator uf1
     ) {
          hpuint e = EDGE_STRIDE;

          auto iv0 = svb;
          auto iv1 = svt;
          auto iv2 = iv1 + 1;

          auto ie1 = sev;
          auto ie0 = ie1 + e;
          auto ie2 = seh;

          auto if0 = sf;

          auto ov1 = tvm;
          auto ov0 = ov1 + 1;
          auto ov2 = tvt + 1;

          auto oe5 = tevt;
          auto oe4 = oe5 + e;
          auto oe3 = oe4 + e;
          auto oe2 = oe3 + e;
          auto oe1 = tevb;
          auto oe0 = oe1 + e;
          auto oe8 = teht;
          auto oe7 = oe8 + e;
          auto oe6 = tehb;

          auto of3 = tft;
          auto of2 = of3 + FACE_STRIDE;
          auto of1 = of2 + FACE_STRIDE;
          auto of0 = tfb;

          doBinaryTriangleSubdivisionStep(iv2, iv0, ov0, r(ie2 + e), r(ie0 + e), r(ie1 + e), r(oe2 + e), r(oe4 + e), oe0, if0, uf0, uf1);
          doBinaryTriangleSubdivisionStep(iv2, iv1, ov2, r(oe2 + e), r(ie2 + e), oe4, r(oe7 + e), oe3, oe8, uf0, of1, of2);
          doBinaryTriangleSubdivisionStep<true>(iv0, iv1, ov1, oe0, r(ie1 + e), oe4, r(oe1 + e), oe6, oe5, uf1, of0, of3);
     }

     //odd row, even column, bottom triangle
     void subdivideOEB(
          PointIterator svt, 
          PointIterator svb, 
          PointIterator sev, 
          PointIterator seh, 
          PointIterator sf, 

          PointIterator tvt, 
          PointIterator tvm, 
          PointIterator tvb, 
          PointIterator tevt, 
          PointIterator tevb, 
          PointIterator teht, 
          PointIterator tehb, 
          PointIterator tft, 
          PointIterator tfb,

          PointIterator uf0,
          PointIterator uf1,

          bool skip2 = false
     ) {
          hpuint e = EDGE_STRIDE;

          auto iv0 = svt;
          auto iv2 = svb;
          auto iv1 = iv2 + 1;

          auto ie0 = sev;
          auto ie1 = ie0 + e;
          auto ie2 = seh;

          auto if0 = sf;

          auto ov0 = tvm;
          auto ov1 = ov0 + 1;
          auto ov2 = tvb + 1;

          auto oe0 = tevt;
          auto oe1 = oe0 + e;
          auto oe2 = tevb;
          auto oe3 = oe2 + e;
          auto oe4 = oe3 + e;
          auto oe5 = oe4 + e;
          auto oe6 = teht;
          auto oe7 = tehb;
          auto oe8 = oe7 + e;

          auto of0 = tft;
          auto of1 = tfb;
          auto of2 = of1 + FACE_STRIDE;
          auto of3 = of2 + FACE_STRIDE;

          doBinaryTriangleSubdivisionStep<true>(iv2, iv0, ov0, ie2, ie0, ie1, oe2, oe4, r(oe0 + e), if0, uf0, uf1);
          doBinaryTriangleSubdivisionStep(iv2, iv1, ov2, oe2, ie2, r(oe4 + e), oe7, r(oe3 + e), r(oe8 + e), uf0, of1, of2, skip2);
          doBinaryTriangleSubdivisionStep(iv0, iv1, ov1, r(oe0 + e), ie1, r(oe4 + e), oe1, r(oe6 + e), r(oe5 + e), uf1, of0, of3);
     }

     //odd row, even column, top triangle
     void subdivideOET(
          PointIterator svt, 
          PointIterator svb, 
          PointIterator sev, 
          PointIterator seh, 
          PointIterator sf, 

          PointIterator tvt, 
          PointIterator tvm, 
          PointIterator tvb, 
          PointIterator tevt, 
          PointIterator tevb, 
          PointIterator teht, 
          PointIterator tehb, 
          PointIterator tft, 
          PointIterator tfb,

          PointIterator uf0,
          PointIterator uf1,

          bool skip2 = false
     ) {
          hpuint e = EDGE_STRIDE;

          auto iv2 = svb;
          auto iv1 = svt;
          auto iv0 = iv1 + 1;

          auto ie2 = sev;
          auto ie0 = ie2 + e;
          auto ie1 = seh;

          auto if0 = sf;

          auto ov2 = tvm;
          auto ov0 = ov2 + 1;
          auto ov1 = tvt + 1;

          auto oe8 = tevt;
          auto oe4 = oe8 + e;
          auto oe6 = oe4 + e;
          auto oe0 = oe6 + e;
          auto oe7 = tevb;
          auto oe2 = oe7 + e;
          auto oe5 = teht;
          auto oe1 = oe5 + e;
          auto oe3 = tehb;

          auto of2 = tft;
          auto of3 = of2 + FACE_STRIDE;
          auto of0 = of3 + FACE_STRIDE;
          auto of1 = tfb;

          doBinaryTriangleSubdivisionStep(iv2, iv0, ov0, r(ie2 + e), ie0, r(ie1 + e), oe2, r(oe4 + e), r(oe0 + e), if0, uf0, uf1);
          doBinaryTriangleSubdivisionStep(iv2, iv1, ov2, oe2, r(ie2 + e), oe4, r(oe7 + e), oe3, oe8, uf0, of1, of2, skip2);
          doBinaryTriangleSubdivisionStep(iv0, iv1, ov1, r(oe0 + e), r(ie1 + e), oe4, r(oe1 + e), oe6, oe5, uf1, of0, of3);
     }

     //odd row, odd column, bottom triangle
     void subdivideOOB(
          PointIterator svt, 
          PointIterator svb, 
          PointIterator sev, 
          PointIterator seh, 
          PointIterator sf, 

          PointIterator tvt, 
          PointIterator tvm, 
          PointIterator tvb, 
          PointIterator tevt, 
          PointIterator tevb, 
          PointIterator teht, 
          PointIterator tehb, 
          PointIterator tft, 
          PointIterator tfb,

          PointIterator uf0,
          PointIterator uf1,

          bool skip2 = false
     ) {
          hpuint e = EDGE_STRIDE;

          auto iv0 = svt;
          auto iv1 = svb;
          auto iv2 = iv1 + 1;

          auto ie1 = sev;
          auto ie0 = ie1 + e;
          auto ie2 = seh;

          auto if0 = sf;

          auto ov1 = tvm;
          auto ov0 = ov1 + 1;
          auto ov2 = tvb + 1;

          auto oe1 = tevt;
          auto oe0 = oe1 + e;
          auto oe5 = tevb;
          auto oe4 = oe5 + e;
          auto oe3 = oe4 + e;
          auto oe2 = oe3 + e;
          auto oe6 = teht;
          auto oe8 = tehb;
          auto oe7 = oe8 + e;

          auto of0 = tft;
          auto of3 = tfb;
          auto of2 = of3 + FACE_STRIDE;
          auto of1 = of2 + FACE_STRIDE;

          doBinaryTriangleSubdivisionStep(iv2, iv0, ov0, r(ie2 + e), r(ie0 + e), ie1, r(oe2 + e), r(oe4 + e), oe0, if0, uf0, uf1);
          doBinaryTriangleSubdivisionStep(iv2, iv1, ov2, r(oe2 + e), r(ie2 + e), oe4, r(oe7 + e), r(oe3 + e), oe8, uf0, of1, of2, skip2);
          doBinaryTriangleSubdivisionStep<true>(iv0, iv1, ov1, oe0, ie1, oe4, oe1, oe6, r(oe5 + e), uf1, of0, of3);
     }

     //odd row, odd column, top triangle
     void subdivideOOT(
          PointIterator svt, 
          PointIterator svb, 
          PointIterator sev, 
          PointIterator seh, 
          PointIterator sf, 

          PointIterator tvt, 
          PointIterator tvm, 
          PointIterator tvb, 
          PointIterator tevt, 
          PointIterator tevb, 
          PointIterator teht, 
          PointIterator tehb, 
          PointIterator tft, 
          PointIterator tfb,

          PointIterator uf0,
          PointIterator uf1
     ) {
          hpuint e = EDGE_STRIDE;

          auto iv0 = svt;
          auto iv1 = iv0 + 1;
          auto iv2 = svb;

          auto ie1 = seh;
          auto ie0 = sev;
          auto ie2 = ie0 + e;

          auto if0 = sf;

          auto ov1 = tvt + 1;
          auto ov0 = tvm;
          auto ov2 = ov0 + 1;

          auto oe0 = tevt;
          auto oe6 = oe0 + e;
          auto oe4 = oe6 + e;
          auto oe8 = oe4 + e;
          auto oe2 = tevb;
          auto oe7 = oe2 + e;
          auto oe1 = teht;
          auto oe5 = oe1 + e;
          auto oe3 = tehb;
     
          auto of0 = tft;
          auto of3 = of0 + FACE_STRIDE;
          auto of2 = of3 + FACE_STRIDE;
          auto of1 = tfb;

          doBinaryTriangleSubdivisionStep<true>(iv2, iv0, ov0, r(ie2 + e), r(ie0 + e), ie1, r(oe2 + e), oe4, oe0, if0, uf0, uf1);
          doBinaryTriangleSubdivisionStep(iv2, iv1, ov2, r(oe2 + e), r(ie2 + e), r(oe4 + e), r(oe7 + e), r(oe3 + e), oe8, uf0, of1, of2);
          doBinaryTriangleSubdivisionStep(iv0, iv1, ov1, oe0, ie1, r(oe4 + e), oe1, oe6, r(oe5 + e), uf1, of0, of3);
     }

};//SurfaceSubdividerBEZ

}//namespace happah

