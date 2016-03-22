// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <vector>

#include "happah/geometries/Sphere.h"
#include "happah/geometries/SurfaceSBB.h"
#include "happah/geometries/TriangleMesh.h"

//NOTE: This is the spherical Clough-Tocher interpolator.
class InterpolatorSCT {
     template<class Iterator>
     using Space = typename std::remove_reference<typename std::tuple_element<3, typename Iterator::value_type>::type>::type;

public:
     //TODO: should actually return Manifold3HEZ?
     template<class Iterator>
     static std::vector<CubicSurfaceSBB<Space<Iterator> > >* interpolate(Iterator begin, Iterator end) { return do_interpolate<Iterator>::exec(begin, end); }

private:

     static Point1D planeEval( const Vector3D& normal, const Point2D& originXY, const Point1D& originZ, const Point2D& eval) {
          if(normal.z == 0.0) { std::cout<<"bad input : normal.z == 0" <<"\n"; return Point1D(originZ); }
          hpreal temp = ((normal.x * (eval.x - originXY.x)) + (normal.y * (eval.y - originXY.y))) / normal.z;
          return Point1D(-1*temp + originZ.x);
     }

     

     static bool findNeighborMicroTriangle(const std::vector<CubicSurfaceSBB<Space1D> >* surfacePatches, CubicSurfaceSBB<Space1D>& neighbor, const Point2D& a0, const Point2D& a1, const hpuint& neighborIndex) {
          if(neighborIndex == happah::UNULL) return false;
          hpuint nIndex = 3 * neighborIndex;
          if(nIndex >= surfacePatches->size()) return false;
          CubicSurfaceSBB<Space1D> n1 = (*surfacePatches)[nIndex];
          CubicSurfaceSBB<Space1D> n2 = (*surfacePatches)[nIndex+1];
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

     static Point1D calcSphericalTangentControlPoint(const Point3D& pA, const Point1D& ordinateA, const Point3D& pB, const Vector3D normal) {
          return Point1D(ordinateA/3.0 * (4 * glm::dot(pA, pB) - (glm::dot(normal, pB) / glm::dot(normal, pA))));
     }

     static void setTangentControlPoints(CubicSurfaceSBB<Space1D>& current, const VertexA2O1N& vA, const Point3D& pA, const VertexA2O1N& vB, const Point3D& pB,  const Point3D& split) {
          current.setControlPoint<2,1,0>(calcSphericalTangentControlPoint(pA, vA.ordinate, pB, vA.normal));
          current.setControlPoint<2,0,1>(calcSphericalTangentControlPoint(pA, vA.ordinate, split, vA.normal));
          current.setControlPoint<1,2,0>(calcSphericalTangentControlPoint(pB, vB.ordinate, pA, vB.normal));
     }

     static void setCornerControlPoints(CubicSurfaceSBB<Space1D>& current, const VertexA2O1N& vA, const VertexA2O1N& vB) {
          current.setControlPoint<3,0,0>(vA.ordinate);
          current.setControlPoint<0,3,0>(vB.ordinate);
     }

     static void joinFirstOpposite(CubicSurfaceSBB<Space1D>& current, const CubicSurfaceSBB<Space1D>& next) {
          current.setControlPoint<0,2,1>(next.getControlPoint<2,0,1>());
     }
     static void joinSecondOpposite(CubicSurfaceSBB<Space1D>& current, const CubicSurfaceSBB<Space1D>& next) {
          current.setControlPoint<0,1,2>(next.getControlPoint<1,0,2>());
     }


     //NOTE: This is the parametric Clough-Tocher implementation.
     template<class Iterator, typename = void>
     struct do_interpolate {
          static std::vector<CubicSurfaceSBB<Space<Iterator> > >* exec(Iterator begin, Iterator end) {}//TODO
     };

     //NOTE: This is the functional Clough-Tocher implementation.
     template<class Iterator>
     struct do_interpolate<Iterator, typename std::enable_if<std::is_same<Space<Iterator>, Space1D>::value>::type> {
          //TODO: change interface to std::tuple<Point2D&, Point2D&, Point2D&, Point&, Point&, Point&, Vector&, Vector&, Vector&, hpuint, hpuint, hpuint>
          static std::vector<CubicSurfaceSBB<Space1D> >* exec(Iterator begin, Iterator end) {
               typedef Vertex<Iterator> Vertex;
               typedef typename Vertex::SPACEA::POINT Abscissa;
               typedef typename Vertex::SPACEO::POINT Ordinate;
              
               std::vector<CubicSurfaceSBB<Space1D> >* spline = new std::vector<CubicSurfaceSBB<Space1D> >();
               spline->reserve(3 * (end - begin));

               for(Iterator i = begin; i != end; ++i) {
                    std::tuple<Vertex&, Vertex&, Vertex&, hpuint, hpuint, hpuint> temp = *i;
                    Vertex& v0 = std::get<0>(temp);
                    Vertex& v1 = std::get<1>(temp);
                    Vertex& v2 = std::get<2>(temp);
                    hpuint n0 = std::get<3>(temp);
                    hpuint n1 = std::get<4>(temp);
                    hpuint n2 = std::get<5>(temp); //TODO

                    //NOTE: Use 3D Coordinates to avoid errors from periodic Abscissa values
                    Point3D p0 = Sphere::Utils::getPoint(v0.abscissa);
                    Point3D p1 = Sphere::Utils::getPoint(v1.abscissa);
                    Point3D p2 = Sphere::Utils::getPoint(v2.abscissa);
                    Point3D splitPoint = 1/3.0f * (p0 + p1 + p2);

                    Abscissa split = Sphere::Utils::getAbscissa(splitPoint);
                    CubicSurfaceSBB<Space1D> surfaceA(v0.abscissa, v1.abscissa, split);
                    CubicSurfaceSBB<Space1D> surfaceB(v1.abscissa, v2.abscissa, split);
                    CubicSurfaceSBB<Space1D> surfaceC(v2.abscissa, v0.abscissa, split);

                    setCornerControlPoints(surfaceA, v0, v1);
                    setCornerControlPoints(surfaceB, v1, v2);
                    setCornerControlPoints(surfaceC, v2, v0); 

                    setTangentControlPoints(surfaceA, v0, p0, v1, p1, splitPoint);
                    setTangentControlPoints(surfaceA, v1, p1, v2, p2, splitPoint);
                    setTangentControlPoints(surfaceA, v2, p2, v0, p0, splitPoint);

                    //TODO: set other control points by using homogenity
               }
               return spline;
          }
     };

};

