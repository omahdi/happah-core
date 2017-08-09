// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <type_traits>
#include <vector>

#include "happah/Happah.hpp"

namespace happah {

template<class T>
struct get_dimension;

template<>
struct get_dimension<hpvec1> : std::integral_constant<hpuint, 1> {};

template<>
struct get_dimension<hpvec2> : std::integral_constant<hpuint, 2> {};

template<>
struct get_dimension<hpvec3> : std::integral_constant<hpuint, 3> {};

template<>
struct get_dimension<hpvec4> : std::integral_constant<hpuint, 4> {};

template<class Point, class Vector = Point>
class SpaceBase {
public:
     using POINT = Point;
     using VECTOR = Vector;

     static constexpr hpuint DIMENSION = get_dimension<Point>::value;
     static_assert(get_dimension<Vector>::value == DIMENSION, "The point and vector parameterizing a space must be of the same dimension.");

};

template<class Point, class Vector = Point>
class Space : public SpaceBase<Point, Vector> {};

template<class Space>
using Point = typename Space::POINT;

template<class Space>
using Vector = typename Space::VECTOR;

template<>
class Space<Point2D, Vector2D> : public SpaceBase<Point2D, Vector2D> {
public:
     static std::vector<hpuint> getQuadIndices(hpuint nx, hpuint ny) {
          assert(!(nx < 2 || ny < 2));//TODO: exception

          std::vector<hpuint> indices;

          hpuint limit = nx * ny - ny;
          indices.reserve((limit - nx + 1) << 2);
          hpuint i = 1;
          hpuint index = 0;
          --limit;
          while(index < limit) {
               if(i == ny) {
                    ++index;
                    i = 1;    
               }

               hpuint opposite = index+ny;

               indices.push_back(index);
               indices.push_back(opposite);
               indices.push_back(++index);
               indices.push_back(++opposite);

               ++i;           
          }

          return std::move(indices);
     }

     static std::vector<hpuint> getSegmentIndices(hpuint nx, hpuint ny, Diagonals diagonals) {
          assert(!(nx == 0 || ny == 0 || (nx == 1 && ny == 1)));//TODO: exception

          std::vector<hpuint> indices;

          if(nx == 1 || ny == 1) {
               hpuint limit;
               if(nx == 1) limit = ny - 1;
               else limit = nx - 1;
               hpuint index = 0;
               indices.reserve(limit << 1);
               while(index < limit) {
                    indices.push_back(index);
                    indices.push_back(++index);
               }
               return std::move(indices);
          }

          hpuint nVertices = nx * ny;
          hpuint nSegments;
          if(diagonals & A || diagonals & B)
               nSegments = nVertices * 3 - ((nx + ny) << 1) + 1;
          else if(diagonals & AB)
               nSegments = (nVertices << 2) - ((nx + ny) * 3) + 2;
          else
               nSegments = (nVertices << 1) - (nx + ny);
          indices.reserve(nSegments << 1);
          hpuint i = 1;
          hpuint index = 0;
          hpuint limit2 = nVertices - 1;
          hpuint limit1 = limit2 - ny;
          if(diagonals & A) {
               while(index < limit1) {
                    if(i == ny) {
                         indices.push_back(index);
                         indices.push_back(index+ny);
                         ++index;
                         i = 1;
                    }

                    hpuint opposite = index+ny;

                    indices.push_back(index);
                    indices.push_back(opposite);

                    indices.push_back(index);
                    indices.push_back(++opposite);
               
                    indices.push_back(index);
                    indices.push_back(++index);
                    
                    ++i;           
               }
          } else if(diagonals & B) {
               while(index < limit1) {
                    if(i == ny) {
                         indices.push_back(index);
                         indices.push_back(index+ny);
                         ++index;
                         i = 1;
                    }

                    hpuint opposite = index+ny;

                    indices.push_back(index);
                    indices.push_back(opposite);
               
                    indices.push_back(index);
                    indices.push_back(++index);

                    indices.push_back(index);
                    indices.push_back(opposite);
                    
                    ++i;           
               }
          } else if(diagonals & AB) {
               while(index < limit1) {
                    if(i == ny) {
                         indices.push_back(index);
                         indices.push_back(index+ny);
                         ++index;
                         i = 1;
                    }

                    hpuint opposite = index+ny;

                    indices.push_back(index);
                    indices.push_back(opposite);

                    indices.push_back(index);
                    indices.push_back(opposite+1);
               
                    indices.push_back(index);
                    indices.push_back(++index);

                    indices.push_back(index);
                    indices.push_back(opposite);
                    
                    ++i;           
               }
          } else {
               while(index < limit1) {
                    if(i == ny) {
                         indices.push_back(index);
                         indices.push_back(index+ny);
                         ++index;
                         i = 1;
                    }

                    indices.push_back(index);
                    indices.push_back(index+ny);
               
                    indices.push_back(index);
                    indices.push_back(++index);
                    
                    ++i;           
               }
          }

          indices.push_back(index);
          indices.push_back(index+ny);
          ++index;
          while(index < limit2) {
               indices.push_back(index);
               indices.push_back(++index);
          }

          return std::move(indices);
     }
          
     static std::vector<hpuint> getTriangleIndices(hpuint nx, hpuint ny) {
          assert(!(nx < 2 || ny < 2));

          std::vector<hpuint> indices;
          
          hpuint limit = nx * ny - ny;
          indices.reserve((limit - nx + 1) * 6);
          hpuint i = 1;
          hpuint index = 0;
          --limit;
          while(index < limit) {
               if(i == ny) {
                    ++index;
                    i = 1;    
               }

               hpuint opposite = index+ny;

               indices.push_back(index);
               indices.push_back(opposite);
               indices.push_back(++opposite);

               indices.push_back(index);
               indices.push_back(opposite);
               indices.push_back(++index);

               ++i;           
          }

          return std::move(indices);
     }

     template<class Visitor>
     static void sample(hpreal xEdgeLength, hpreal yEdgeLength, hpuint nx, hpuint ny, Visitor& visit) {
          hpreal xDelta = xEdgeLength / hpreal(nx - 1);
          hpreal yDelta = yEdgeLength / hpreal(ny - 1);
          hpreal x = -0.5f * xEdgeLength;
          hpreal yMin = -0.5f * yEdgeLength;
          for(hpuint i = 0; i < nx; ++i) {
               hpreal y = yMin;
               for(hpuint j = 0; j < ny; ++j) {
                    visit(x, y);
                    y += yDelta;
               }
               x += xDelta;
          }
     }

     template<class Space, class Visitor>
     static void sample(const Point<Space>& origin, const Vector<Space>& xAxis, const Vector<Space>& yAxis, hpreal xEdgeLength, hpreal yEdgeLength, hpuint nx, hpuint ny, Visitor& visit) {
          Vector<Space> xDelta = (xEdgeLength / hpreal(nx - 1)) * xAxis;
          Vector<Space> yDelta = (yEdgeLength / hpreal(ny - 1)) * yAxis;
          Point<Space> position = -0.5f * (xEdgeLength * xAxis + yEdgeLength * yAxis) + origin;
          for(hpuint i = 0; i < nx; ++i) {
               Vector<Space> temp = position;
               for(hpuint j = 0; j < ny; ++j) {
                    visit(position);
                    position += yDelta;
               }
               position = temp + xDelta;
          }
     }

     template<class T, class Visitor>
     static void visitNeighborhoods(const std::vector<T>& ts, hpuint nx, hpuint ny, Visitor& visit) {
          assert(nx > 0 && ny > 0 && (nx > 1 || ny > 1));

          if(nx == 1 || ny == 1) {
               auto c = ts.cbegin(), l = c, e = ts.cend() - 1;
               visit(*c, *(++c));
               while(c != e) {
                    visit(*c, *l, *(++c));
                    ++l;
               }
               visit(*c, *l);
               return;
          }

          hpuint column = 1;
          hpuint columne = ny - 1;
          hpuint row = 1;
          hpuint rowe = nx - 1;

          auto l = ts.cbegin(), bl = l, tl = l + ny, temp = l;//left, bottom left, top left
          
          //first point first row
          visit(*l, *(l+1), *tl, *(tl+1));//c, c+1, c+ny, c+ny+1

          //middle points first row
          while(column < columne) {
               temp = l;
               visit(*(++l), *temp, *(l+1), *tl, *(++tl), *(tl+1));//c, c-1, c+1, c+ny-1, c+ny, c+ny+1
               ++column;
          }

          //last point first row
          temp = l;
          visit(*(++l), *temp, *tl, *(++tl));//c, c-1, c+ny-1, c+ny

          ++l;
          ++tl;

          column = 1;
          
          //middle rows
          while(row < rowe) {
               //first point middle row
               visit(*l, *bl, *(bl+1), *(l+1), *tl, *(tl+1));//c, c-ny, c-ny+1, c+1, c+ny, c+ny+1

               //middle points middle row
               while(column < columne) {
                    temp = l;
                    visit(*(++l), *bl, *(++bl), *(bl+1), *temp, *(l+1), *tl, *(++tl), *(tl+1));//c, c-ny-1, c-ny, c-ny+1, c-1, c+1, c+ny-1, c+ny, c+ny+1
                    ++column;
               }

               //last point middle row
               temp = l;
               visit(*(++l), *bl, *(++bl), *temp, *tl, *(++tl));//c, c-ny-1, c-ny, c-1, c+ny-1, c+ny

               ++bl;
               ++l;
               ++tl;

               column = 1;
               ++row;
          }

          //first point last row
          visit(*l, *bl, *(bl+1), *(l+1));//c, c-ny, c-ny+1, c+1

          //middle points last row
          while(column < columne) {
               temp = l;
               visit(*(++l), *bl, *(++bl), *(bl+1), *temp, *(l+1));//c, c-ny-1, c-ny, c-ny+1, c-1, c+1
               ++column;
          }

          //last point last row
          temp = l;
          visit(*(++l), *bl, *(++bl), *temp);//c, c-ny-1, c-ny, c-1
     }

};//Space

template<>
class Space<Point3D, Vector3D> : public SpaceBase<Point3D, Vector3D> {
public:
     static Point3D toPoint(const Point1D& a, const Point1D& b, const Point1D& c) { return {a.x, b.x, c.x}; }

     static Point3D toPoint(const Point2D& a, const Point1D& b) { return {a.x, a.y, b.x}; }

     static Point3D toPoint(const Point1D& a, const Point2D& b) { return {a.x, b.x, b.y}; }

};//Space

using Space1D = Space<Point1D, Vector1D>;
using Space2D = Space<Point2D, Vector2D>;
using Space3D = Space<Point3D, Vector3D>;
using Space4D = Space<Point4D, Vector4D>;

}//namespace happah

