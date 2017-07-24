// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/Happah.hpp"
#include "happah/geometries/Geometry.hpp"
#include "happah/geometries/Ray.hpp"
#include "happah/geometries/TriangleMesh.hpp"
#include "happah/geometries/VertexCloud.hpp"
#include "happah/utils/VertexFactory.hpp"

namespace happah {

//TODO: render spheres as points with radius using render vertex cloud program
class Sphere : public Geometry2D<Space3D> {
public:
     struct Utils {
          static Point2D getAbscissa(const Point3D& position);
          static Point2D getAbscissa(const Point3D& center, const Point3D& position);
          static Vector3D getNormal(hpreal latitude, hpreal longitude);
          static Vector3D getNormal(const Point2D& abscissa);
          static Point3D getPoint(hpreal latitude, hpreal longitude);
          static Point3D getPoint(const Point2D& abscissa);
          static Point3D getPoint(const Point2D& abscissa, const Point1D& ordinate);
          static hpuint intersect(const Point3D& center, hpreal radius, const Ray3D& ray, hpreal& foot, hpreal& delta, hpreal epsilon = happah::EPSILON);
          //NOTE: Computes first intersection, if it exists.
          static bool intersect(const Point3D& center, hpreal radius, const Ray3D& ray, Point3D& intersection, hpreal epsilon = happah::EPSILON);
          static hpuint intersect(const Point3D& center, hpreal radius, const Ray3D& ray, Point3D& intersection1, Point3D& intersection2, hpreal epsilon = happah::EPSILON);
          /** The method determines whether there exists a sphere with the given radius that touches the three points and if such a sphere exists, stores in the referenced locations the two possible centers for this sphere.
           *  @param p1 first point
           *  @param p2 second point
           *  @param p3 third point
           *  @param radius radius of the desired sphere
           *  @param center1 reference to where the first center should be stored
           *  @param center2 reference to where the second center should be stored
           *  @return false if such a sphere does not exist and true if there exists such a sphere
           */
          static bool isSphere(const Point3D& p1, const Point3D& p2, const Point3D& p3, hpreal radius, Point3D& center1, Point3D& center2);

     private:
          template<class T>
          static T getCoordinates(hpreal latitude, hpreal longitude);

     };

     Sphere(const Point3D& center, hpreal radius);
     ~Sphere();

     Point2D getAbscissa(const Point3D& position) const;
     const Point3D& getCenter() const;
     hpreal getRadius() const;
     Point3D getPoint(const Point2D& abscissa) const;
     Point3D getPoint(const Point2D& abscissa, const Point1D& ordinate) const;
     hpuint intersect(const Ray3D& ray, hpreal& foot, hpreal& delta, hpreal epsilon = happah::EPSILON) const;
     //NOTE: Computes first intersection, if it exists.
     bool intersect(const Ray3D& ray, Point3D& intersection, hpreal epsilon = happah::EPSILON) const;
     hpuint intersect(const Ray3D& ray, Point3D& intersection1, Point3D& intersection2, hpreal epsilon = happah::EPSILON) const;
     //TODO: take advantage of sphere symmetry (8 quadrants) as in Circle::Utils::sample
     template<class Vertex>
     std::vector<Vertex>* sample(hpuint nLatitudes, hpuint nLongitudes) const {
          static_assert(is_vertex<Vertex, Space3D>::value, "The sample method can only be parameterized by a vertex in 3D space.");

          hpuint nLatitudesPlus1 = nLatitudes+1;
          hpreal phi = M_PI/nLatitudesPlus1;
          hpreal theta = 2.0*M_PI/nLongitudes;

          std::vector<Vertex>* vertices = new std::vector<Vertex>();
          vertices->reserve(nLatitudes*nLongitudes+2);
          vertices->push_back(happah::VertexFactory<Vertex>()(Point3D(m_center.x, m_center.y, m_center.z+m_radius), Vector3D(0.0,0.0,1.0)));
          for(hpuint i = 1; i < nLatitudesPlus1; ++i) {
               hpreal iPhi = M_PI/2.0-i*phi;
               hpreal cosIPhi = cos(iPhi);
               hpreal sinIPhi = sin(iPhi);
               hpreal z = m_radius*sinIPhi;
               for(hpuint j = 0; j < nLongitudes; ++j) {
                    hpreal jTheta = j*theta;
                    hpreal cosJTheta = cos(jTheta);
                    hpreal sinJTheta = sin(jTheta);
                    hpreal cosIPhicosJTheta = cosIPhi*cosJTheta;
                    hpreal cosIPhisinJTheta = cosIPhi*sinJTheta;
                    hpreal x = m_radius*cosIPhicosJTheta;
                    hpreal y = m_radius*cosIPhisinJTheta;
                    vertices->push_back(happah::VertexFactory<Vertex>()(Point3D(x+m_center.x,y+m_center.y,z+m_center.z), Vector3D(cosIPhicosJTheta,cosIPhisinJTheta,sinIPhi)));
               }
          }
          vertices->push_back(happah::VertexFactory<Vertex>()(Point3D(m_center.x, m_center.y, m_center.z-m_radius), Vector3D(0.0,0.0,-1.0)));
          return vertices;
     }
     PointCloud3D* toPointCloud(hpuint nLatitudes = 50, hpuint nLongitudes = 50) const;//TODO: toVertexCloud?
     //TODO: is it possible to make these to...Mesh methods more efficient?
     /*template<class Vertex = typename SegmentMesh3D::VERTEX>
     SegmentMesh<Vertex>* toSegmentMesh(hpuint nLatitudes = 50, hpuint nLongitudes = 50) const {
          static_assert(is_vertex<Vertex, Space3D>::value, "The toSegmentMesh method can only be parameterized by a vertex in 3D space.");

          std::vector<Vertex>* vertices = sample<Vertex>(nLatitudes, nLongitudes);

          std::vector<hpuint>* indices = new std::vector<hpuint>();
          indices->reserve(2*nLatitudes*nLongitudes+nLongitudes);
          hpuint nLongitudesPlus1 = nLongitudes+1;
          for(hpuint j = 1; j < nLongitudesPlus1; ++j) {
               indices->push_back(0);
               indices->push_back(j);
          }
          hpuint indexLastPointMinus1 = nLatitudes*nLongitudes;
          hpuint indexLastPoint = indexLastPointMinus1+1;
          hpuint j0 = indexLastPoint-nLongitudes;
          for(hpuint j = j0; j < indexLastPointMinus1;) {
               indices->push_back(indexLastPoint);
               indices->push_back(j);
               indices->push_back(j);
               indices->push_back(++j);
          }
          indices->push_back(indexLastPoint);
          indices->push_back(indexLastPointMinus1);
          indices->push_back(indexLastPointMinus1);
          indices->push_back(j0);
          hpuint j = 0;
          hpuint nLatitudesMinus1 = nLatitudes-1;
          for(hpuint i = 0; i < nLatitudesMinus1; ++i) {
               hpuint jLimit = j+nLongitudes;
               ++j;
               hpuint j0 = j;
               while(j < jLimit) {
                    indices->push_back(j);
                    indices->push_back(j+nLongitudes);
                    indices->push_back(j);
                    indices->push_back(++j);
               }
               indices->push_back(j);
               indices->push_back(j+nLongitudes);
               indices->push_back(j);
               indices->push_back(j0);
          }
          return new SegmentMesh<Vertex>(vertices, indices);
     }*/
     template<class Vertex>
     TriangleMesh<Vertex>* toTriangleMesh(hpuint nLatitudes = 50, hpuint nLongitudes = 50) const {
          static_assert(is_vertex<Vertex, Space3D>::value, "The toTriangleMesh method can only be parameterized by a vertex in 3D space.");

          std::vector<Vertex>* vertices = sample<Vertex>(nLatitudes, nLongitudes);

          std::vector<hpuint>* indices = new std::vector<hpuint>();
          indices->reserve(6*nLongitudes*nLatitudes);
          for(hpuint j = 1; j < nLongitudes;) {
               indices->push_back(j);
               indices->push_back(0);
               indices->push_back(++j);
          }
          indices->push_back(nLongitudes);
          indices->push_back(0);
          indices->push_back(1);
          hpuint indexLastPointMinus1 = nLatitudes*nLongitudes;
          hpuint indexLastPoint = indexLastPointMinus1+1;
          hpuint j0 = indexLastPoint-nLongitudes;
          for(hpuint j = j0; j < indexLastPointMinus1;) {
               indices->push_back(indexLastPoint);
               indices->push_back(j);
               indices->push_back(++j);
          }
          indices->push_back(indexLastPoint);
          indices->push_back(indexLastPointMinus1);
          indices->push_back(j0);
          for(hpuint i = 1; i < indexLastPoint-nLongitudes;) {
               hpuint iPlusNLongitudes = i+nLongitudes;
               indices->push_back(i);
               indices->push_back(iPlusNLongitudes);
               indices->push_back(iPlusNLongitudes-1);
               indices->push_back(i);
               indices->push_back(++i);
               indices->push_back(iPlusNLongitudes);
          }

          return new TriangleMesh<Vertex>(vertices, indices);
     }

private:
     Point3D m_center;
     hpreal m_radius;

};
typedef std::shared_ptr<Sphere> Sphere_ptr;

}//namespace happah

