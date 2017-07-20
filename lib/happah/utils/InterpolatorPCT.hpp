// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "happah/geometries/SurfaceBEZ.hpp"
#include "happah/geometries/TriangleMesh.hpp"
#include <iostream>

//TODO: move sinWave/rand interpolation example into static test class
//NOTE: This is the planar Clough-Tocher interpolator.
//TODO: HermiteDataIterator
//TODO: clean this class up
//TODO: clean surfacesbb, manifold3hez, interpolatersct
//TODO: tensor product surface taking two functions
class InterpolatorPCT {
public:
     template<class Iterator>
     using Vertex = typename std::remove_reference<typename std::tuple_element<0, typename Iterator::value_type>::type>::type;

     template<class Iterator>
     static std::vector<CubicSurfaceBEZ<typename Vertex<Iterator>::SPACEO> >* interpolate(Iterator begin, Iterator end) {
          static_assert(std::is_same<typename Vertex<Iterator>::SPACEA, Space2D>::value, "The Clough-Tocher interpolator can only be used to construct surfaces.");
          return do_interpolate<Iterator>::exec(begin, end);
     }

private:

     //TODO: move to BezierUtils?
     static Point2D getAbscissa(const CubicSurfaceBEZ<Space1D>& surface, const hpuint& i0, const hpuint& i1, const  hpuint& i2) {
          //TODO: error case if sum_i > surface degree
          auto triangle = surface.getParameterTriangle();
          const Point2D& t0 = std::get<0>(triangle);
          const Point2D& t1 = std::get<1>(triangle);
          const Point2D& t2 = std::get<2>(triangle);
          return (1/((hpreal)i0+(hpreal)i1+(hpreal)i2)) * ((hpreal)i0 * t0 + (hpreal)i1 * t1 + (hpreal)i2 * t2);
     }

     static Point1D planeEval( const Vector3D& normal, const Point2D& originXY, const  Point1D& originZ, const Point2D& eval) {
          if(normal.z == 0.0) { std::cout<<"bad input : normal.z == 0" <<"\n"; return Point1D(originZ); }
          hpreal temp = ((normal.x * (eval.x - originXY.x)) + (normal.y * (eval.y - originXY.y))) / normal.z;
          return Point1D(-1*temp + originZ.x);
     }

          //TODO: maybe generalize this for n-splits
     static bool findNeighborMicroTriangle(const std::vector<CubicSurfaceBEZ<Space1D> >* surfacePatches, CubicSurfaceBEZ<Space1D>& neighbor, const Point2D& a0, const Point2D& a1, const hpuint& neighborIndex) {
          if(neighborIndex == happah::UNULL) return false;
          hpuint nIndex = 3 * neighborIndex;
          if(nIndex >= surfacePatches->size()) return false;
          CubicSurfaceBEZ<Space1D> n1 = (*surfacePatches)[nIndex];
          CubicSurfaceBEZ<Space1D> n2 = (*surfacePatches)[nIndex+1];
          neighbor = (*surfacePatches)[nIndex+2];
          std::tuple<const Point2D&, const Point2D&, const Point2D&> pA = n1.getParameterTriangle();
          std::tuple<const Point2D&, const Point2D&, const Point2D&> pB = n2.getParameterTriangle();
          Point2D neighborAFirst = std::get<0>(pA);
          Point2D neighborASecond = std::get<1>(pA);
          Point2D neighborBFirst = std::get<0>(pB);
          Point2D neighborBSecond = std::get<1>(pB);
          if((neighborAFirst == a0 && neighborASecond == a1) || (neighborAFirst == a1 && neighborASecond == a0)) neighbor = n1;
          if((neighborBFirst == a0 && neighborBSecond == a1) || (neighborBFirst == a1 && neighborBSecond == a0)) neighbor = n2;
          return true;
     }

     static void setControlPoint111(const std::vector<CubicSurfaceBEZ<Space1D> >* surfacePatches, CubicSurfaceBEZ<Space1D>& current,  Point2D& a0,  Point2D& a1, const hpuint& neighborIndex) {
          CubicSurfaceBEZ<Space1D> neighbor(Point2D(0), Point2D(0), Point2D(0));
          if(findNeighborMicroTriangle(surfacePatches, neighbor, a0, a1, neighborIndex)) {
               Point2D abs1 = getAbscissa(current,2,1,0);
               Point2D abs2 = getAbscissa(current,1,2,0);
               Point2D abs3 = getAbscissa(neighbor,1,1,1);
               Point1D ord1 = current.getControlPoint<2,1,0>();
               Point3D p1(abs1.x, abs1.y, ord1.x);
               Point3D p2(abs2.x, abs2.y, current.getControlPoint<1,2,0>().x);
               Point3D p3(abs3.x, abs3.y, neighbor.getControlPoint<1,1,1>().x);
               current.setControlPoint<1,1,1>(planeEval(glm::cross((p1 - p3), (p2 - p3)), abs1, ord1, getAbscissa(current,1,1,1)));
          } else {
               current.setControlPoint<1,1,1>(1/2.f * (current.getControlPoint<2,0,1>() + current.getControlPoint<0,2,1>()));//TODO:find good choice
          }
     }

     static void setCornerControlPoints(CubicSurfaceBEZ<Space1D>& current, const VertexA2O1N& vA, const VertexA2O1N& vB) {
          current.setControlPoint<3,0,0>(vA.ordinate);
          current.setControlPoint<0,3,0>(vB.ordinate);
     }

     static void setTangentControlPoints(CubicSurfaceBEZ<Space1D>& current,  const VertexA2O1N& vA,  const VertexA2O1N& vB) {
          current.setControlPoint<2,1,0>(planeEval(vA.normal, vA.abscissa, vA.ordinate, getAbscissa(current,2,1,0)));
          current.setControlPoint<2,0,1>(planeEval(vA.normal, vA.abscissa, vA.ordinate, getAbscissa(current,2,0,1)));
          current.setControlPoint<1,2,0>(planeEval(vB.normal, vB.abscissa, vB.ordinate, getAbscissa(current,1,2,0)));
     }

     static void setInnerSplitControlPoint(CubicSurfaceBEZ<Space1D>& current, const CubicSurfaceBEZ<Space1D>& prev) {
          Point1D splitOrdinate = 1/3.f * (current.getControlPoint<2,0,1>() + current.getControlPoint<1,1,1>() + prev.getControlPoint<1,1,1>());
          current.setControlPoint<1,0,2>(splitOrdinate);
     }

     static void joinFirstOpposite(CubicSurfaceBEZ<Space1D>& current, const CubicSurfaceBEZ<Space1D>& next) {
          current.setControlPoint<0,2,1>(next.getControlPoint<2,0,1>());
     }
     static void joinSecondOpposite(CubicSurfaceBEZ<Space1D>& current, const CubicSurfaceBEZ<Space1D>& next) {
          current.setControlPoint<0,1,2>(next.getControlPoint<1,0,2>());
     }

     //NOTE: This is the parametric Clough-Tocher implementation.
     template<class Iterator, typename = void>
     struct do_interpolate {
          static CubicSurfaceBEZ<typename Vertex<Iterator>::SPACEO> exec(Iterator begin, Iterator end) {}//TODO
     };

     //NOTE: This is the functional Clough-Tocher implementation.
     template<class Iterator>
     struct do_interpolate<Iterator, typename std::enable_if<std::is_same<typename Vertex<Iterator>::SPACEO, Space1D>::value>::type> {
          //TODO: change interface to std::tuple<Point&, Point&, Point&, Vector&, Vector&, Vector&, hpuint, hpuint, hpuint>
          static std::vector<CubicSurfaceBEZ<Space1D> >* exec(Iterator begin, Iterator end) {
               typedef Vertex<Iterator> Vertex;
               typedef typename Vertex::SPACEA::POINT Abscissa;
               typedef typename Vertex::SPACEO::POINT Ordinate;
               hpuint nTriangles = 3 * (end - begin);
               std::vector<CubicSurfaceBEZ<Space1D> >* surfacePCT = new std::vector<CubicSurfaceBEZ<Space1D> >();
               surfacePCT->reserve(nTriangles);
               for(Iterator i = begin; i != end; ++i) {
                    std::tuple<Vertex&, Vertex&, Vertex&, hpuint, hpuint, hpuint> temp = *i;
                    Vertex& v0 = std::get<0>(temp);
                    Vertex& v1 = std::get<1>(temp);
                    Vertex& v2 = std::get<2>(temp);
                    hpuint n0 = std::get<3>(temp);
                    hpuint n1 = std::get<4>(temp);
                    hpuint n2 = std::get<5>(temp);
                   
                    Abscissa split = 1/3.f * (v0.abscissa + v1.abscissa + v2.abscissa);
                    CubicSurfaceBEZ<Space1D> surfaceA(v0.abscissa, v1.abscissa, split);
                    CubicSurfaceBEZ<Space1D> surfaceB(v1.abscissa, v2.abscissa, split);
                    CubicSurfaceBEZ<Space1D> surfaceC(v2.abscissa, v0.abscissa, split);

                    setCornerControlPoints(surfaceA, v0, v1);
                    setCornerControlPoints(surfaceB, v1, v2);
                    setCornerControlPoints(surfaceC, v2, v0);                    

                    setTangentControlPoints(surfaceA, v0, v1);
                    setTangentControlPoints(surfaceB, v1, v2);
                    setTangentControlPoints(surfaceC, v2, v0);

                    joinFirstOpposite(surfaceA, surfaceB);
                    joinFirstOpposite(surfaceB, surfaceC);
                    joinFirstOpposite(surfaceC, surfaceA);

                    setControlPoint111(surfacePCT, surfaceA, v0.abscissa, v1.abscissa, n0);
                    setControlPoint111(surfacePCT, surfaceB, v1.abscissa, v2.abscissa, n1);
                    setControlPoint111(surfacePCT, surfaceC, v2.abscissa, v0.abscissa, n2);

                    setInnerSplitControlPoint(surfaceA, surfaceC);
                    setInnerSplitControlPoint(surfaceB, surfaceA);
                    setInnerSplitControlPoint(surfaceC, surfaceB);

                    Ordinate splitOrdinate = 1/3.f * (surfaceA.getControlPoint<1,0,2>() + surfaceB.getControlPoint<1,0,2>() + surfaceC.getControlPoint<1,0,2>());
                    surfaceA.setControlPoint<0,0,3>(splitOrdinate);
                    surfaceB.setControlPoint<0,0,3>(splitOrdinate);
                    surfaceC.setControlPoint<0,0,3>(splitOrdinate);

                    joinSecondOpposite(surfaceA, surfaceB);
                    joinSecondOpposite(surfaceB, surfaceC);
                    joinSecondOpposite(surfaceC, surfaceA);
                   
                    surfacePCT->push_back(surfaceA);
                    surfacePCT->push_back(surfaceB);
                    surfacePCT->push_back(surfaceC);
               }
               return surfacePCT;
          }
     };

};

