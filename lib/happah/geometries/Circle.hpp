// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/optional.hpp>
#include <cmath>
#include <vector>

#include "happah/Happah.hpp"
#include "happah/math/Space.hpp"

namespace happah {

//DECLARATIONS

class Circle;

Circle make_circle(Point2D center, hpreal radius);

boost::optional<std::tuple<Point2D, Point2D> > intersect(const Circle& circle0, const Circle& circle1);

//DEFINITIONS

class Circle {
public:
     Circle(Point2D center, hpreal radius);

     const Point2D& getCenter() const;

     hpreal getRadius() const;

private:
     Point2D m_center;
     hpreal m_radius;

};//Circle

template<class Visitor>
void sample(hpreal radius, hpuint nWedges, hpreal offset, Visitor&& visit) {
     hpreal delta = 2.0 * M_PI / (hpreal) nWedges;
     offset = std::fmod(offset, delta);
     hpreal phi = offset, theta = 0, x, y;

     if(!(nWedges & 3)) {
          visit(radius, 0.0);
          visit(0.0, radius);
          visit(-radius, 0.0);
          visit(0.0, -radius);

          x = radius * std::cos(phi);
          y = radius * std::sin(phi);
          visit(x, y);
          visit(-y, x);
          visit(-x, -y);
          visit(y, -x);

          hpuint limit = (nWedges >> 2) - 1;
          for(hpuint i = 0; i < limit; ++i) {
               theta += delta;
               phi += delta;

               x = radius * std::cos(theta);
               y = radius * std::sin(theta);
               visit(x, y);
               visit(-y, x);
               visit(-x, -y);
               visit(y, -x);

               x = radius * std::cos(phi);
               y = radius * std::sin(phi);
               visit(x, y);
               visit(-y, x);
               visit(-x, -y);
               visit(y, -x);
          }
     } else if(!(nWedges & 1)) {
          visit(radius, 0.0);
          visit(-radius, 0.0);

          x = radius * std::cos(phi);
          y = radius * std::sin(phi);
          visit(x, y);
          visit(-x, -y);

          hpuint limit = (nWedges >> 1) - 1;
          for(hpuint i = 0; i < limit; ++i) {
               theta += delta;
               phi += delta;

               x = radius * std::cos(theta);
               y = radius * std::sin(theta);
               visit(x, y);
               visit(-x, -y);

               x = radius * std::cos(phi);
               y = radius * std::sin(phi);
               visit(x, y);
               visit(-x, -y);
          }
     } else {
          visit(radius, 0.0);

          x = radius * std::cos(phi);
          y = radius * std::sin(phi);
          visit(x, y);

          hpuint limit = nWedges - 1;
          for(hpuint i = 0; i < limit; ++i) {
               theta += delta;
               phi += delta;

               x = radius * std::cos(theta);
               y = radius * std::sin(theta);
               visit(x, y);

               x = radius * std::cos(phi);
               y = radius * std::sin(phi);
               visit(x, y);
          }
     }
}

//NOTE: This method needs to be adapted if more sample methods are introduced.
template<class T, class NeighborhoodIndicesVisitor>
void visitNeighborhoodIndices(std::vector<T>& ts, hpuint nWedges, hpuint offset, NeighborhoodIndicesVisitor& visitor) {
     typename std::vector<T>::iterator n = ts.begin() + offset, nccw, ncw;
     hpuint i = offset, iccw, icw;
     bool ocw;
     hpuint nSamples = nWedges << 1;
     if(!(nWedges & 3)) {
          nccw = n + 4;
          iccw = i + 4;
          ncw = n + nSamples - 1;     
          icw = i + nSamples - 1;
          visitor.visit(*n, i, *nccw, iccw, *ncw, icw, true);
          visitor.visit(*ncw, icw, *n, i, *(ncw - 4), icw - 4, false);
        
          ncw -= 3; 
          icw -= 3; 
          visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *ncw, icw, true);
          visitor.visit(*ncw, icw, *n, i, *(ncw - 4), icw - 4, false);

          visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, true);
          visitor.visit(*ncw, icw, *n, i, *(ncw - 4), icw - 4, false);

          visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, true);
          visitor.visit(*ncw, icw, *n, i, *(ncw - 4), icw - 4, false);

          nccw = n + 4;
          iccw = i + 4;
          ncw = n - 4;
          icw = i - 4;
          ocw = false;
          hpuint limit = (nSamples - 8) >> 2;
          for(hpuint j = 0; j < limit; ++j) {
               visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, ocw);
               visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, ocw);
               visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, ocw);
               visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, ocw);
               ocw = !ocw;
          }          
     } else if (!(nWedges & 1)) {
          nccw = n + 2;
          iccw = i + 2;
          ncw = n + nSamples - 1;
          icw = i + nSamples - 1;
          visitor.visit(*n, i, *nccw, iccw, *ncw, icw, true);
          visitor.visit(*ncw, icw, *n, i, *(ncw - 2), icw - 2, false);

          --ncw;
          --icw;
          visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *ncw, icw, true);
          visitor.visit(*ncw, icw, *n, i, *(ncw - 2), icw - 2, false);

          nccw = n + 2;
          iccw = i + 2;
          ncw = n - 2;
          icw = i - 2;
          ocw = false;
          hpuint limit = (nSamples - 4) >> 1;
          for(hpuint j = 0; j < limit; ++j) {
               visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, ocw);
               visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, ocw);
               ocw = !ocw;
          }        
     } else if (nWedges == 1) {
          nccw = n + 1;
          iccw = i + 1;
          visitor.visit(*n, i, *nccw, iccw, false);
          visitor.visit(*nccw, iccw, *n, i, true);
     } else {
          nccw = n + 1;
          iccw = i + 1;
          ncw = n + nSamples - 1;
          icw = i + nSamples - 1;
          visitor.visit(*n, i, *nccw, iccw, *ncw, icw, true);
          visitor.visit(*ncw, icw, *n, i, *(ncw - 1), icw - 1, false);

          ncw = n - 1; 
          icw = i - 1;
          ocw = false;
          hpuint limit = ((nWedges - 1) << 1);
          for(hpuint j = 0; j < limit; ++j) {
               visitor.visit(*(++n), ++i, *(++nccw), ++iccw, *(++ncw), ++icw, ocw);
               ocw = !ocw;
          }
     }
}

}//namespace happah

