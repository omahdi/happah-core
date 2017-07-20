// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

//#include <boost/dynamic_bitset>
#include <vector>
#include <stack>
#include <glm/glm.hpp>
#include "happah/geometries/Vertex.hpp"
#include "happah/math/Space.hpp"
#include "happah/utils/SegmentUtils.hpp"

//TODO: reference Hess thesis

//NOTE: Trapezoidation implementation is based on Seidel. 
//TODO: reference Seidel paper
template<class Iterator, class Points>
class Trapezoidation {
     static bool isRelativCounterClockwise(const std::tuple<Point2D&, Point2D&> segment, const Point2D& point, const hpreal& epsilon = happah::EPSILON) { return doIsLeftOf(std::get<0>(segment), std::get<1>(segment), point, epsilon); }
     static bool isRelativCounterClockwise(const Point2D& p1, const Point2D& p2, const Point2D& point, const hpreal& epsilon = happah::EPSILON) { return doIsLeftOf(p1,p2,point,epsilon); }

     static bool doIsLeftOf(const Point2D& p1, const Point2D& p2, const Point2D& point, const hpreal& epsilon = happah::EPSILON) {
          Vector2D v = point - p1;
          hpreal temp = glm::dot(p2 - p1, v);
          if(glm::abs(temp) < epsilon) return v.x < 0.0;
          return temp < 0.0;
     }

public:
     struct Trapezoid { //TODO: choose better names
          hpuint topPoint;
          hpuint bottomPoint;
          hpuint leftSegment;
          hpuint rightSegment;
          hpuint topTrapezoid1;
          hpuint topTrapezoid2;
          hpuint bottomTrapezoid1;
          hpuint bottomTrapezoid2;
          hpuint leaf;

          Trapezoid() :
                    topPoint(-1), bottomPoint(-1), leftSegment(-1), rightSegment(-1), topTrapezoid1(-1), topTrapezoid2(-1), bottomTrapezoid1(-1), bottomTrapezoid2(-1), leaf(0) {
          }
          Trapezoid(Trapezoid& t) :
                    topPoint(t.topPoint), bottomPoint(t.bottomPoint), leftSegment(t.leftSegment), rightSegment(t.rightSegment), topTrapezoid1(t.topTrapezoid1), topTrapezoid2(t.topTrapezoid2), bottomTrapezoid1(t.bottomTrapezoid1), bottomTrapezoid2(
                              t.bottomTrapezoid2), leaf(t.topPoint) {
          }

     };

     //template<class Iterator, class Points>
     Trapezoidation(Iterator start, Iterator end, const Points& points) {
          /* Iterator get every time an hpuint. Every second is a segment.*/
     }

     //Only for testing - do not use this construktor otherwise
     Trapezoidation(Iterator begin, Iterator end, const Points& points, std::vector<Trapezoid>& trapezoids, std::vector<hpuint>& mapping) {
          m_begin = begin;
          m_end = end;
          m_points = points;
          m_trapezoids = trapezoids;
          m_point_trapezoid = mapping;

     }

     const Trapezoid& find(const Point2D& point) const {
          return doFind(point);
     }

     template<class Vertex>
     void getHorizontalExtensions(std::vector<Vertex>& vertices, std::vector<hpuint>& indices) {
          //to create a vertex:    VertexFactory<Vertex>::build(point);
          for (int i = 0; i < m_trapezoids.size(); ++i) {
               int t = m_trapezoids[i].topPoint, b = m_trapezoids[i].bottomPoint, l = m_trapezoids[i].leftSegment, r = m_trapezoids[i].rightSegment;
               int last = vertices.size();
               if (l == -1 || r == -1) {
                    // no output
               } else {
                    if (t == -1) {
                         // get second topmost point
                         std::tuple<Point2D, hpuint> tmp = getCut(l, r, true);
                         vertices.push_pack(VertexFactory < Vertex > ::build(std::get < 0 > (tmp)));
                         indices.push_back(last);
                         last++;
                         indices.push_back(std::get < 1 > (tmp));
                    } else {
                         // intersect l and r segment with topHE
                         vertices.push_back(VertexFactory < Vertex > ::build(SegmentUtils::intersect(getSegment(l), getPoint(t))));
                         vertices.push_back(VertexFactory < Vertex > ::build(SegmentUtils::intersect(getSegment(r), getPoint(t))));
                         indices.push_back(last);
                         last++;
                         indices.push_back(last);
                         last++;
                    }
                    if (b == -1) {
                         // cut with opposite segment
                         std::tuple<Point2D, hpuint> tmp = getCut(l, r, false);
                         vertices.push_pack(VertexFactory < Vertex > ::build(std::get < 0 > (tmp)));
                         indices.push_back(last);
                         last++;
                         indices.push_back(std::get < 1 > (tmp));
                    } else {
                         // intersect l and r segment with bottomHE
                         vertices.push_back(VertexFactory < Vertex > ::build(SegmentUtils::intersect(getSegment(l), getPoint(b))));
                         vertices.push_back(VertexFactory < Vertex > ::build(SegmentUtils::intersect(getSegment(r), getPoint(b))));
                         indices.push_back(last);
                         last++;
                         indices.push_back(last);
                         last++;
                    }
               }
          }
     }

     void getUnimonotonePolygons(std::vector<hpuint>& chains, std::vector<hpuint>& lengths) {
          /* The input should include i_1.. i_n n indicies being point (i1, ... in is diagonal of monotone Polygon) and length[j] = n, being n the length of jth chain*/
          std::stack<hpuint> stack1, stack2;
          std::vector<hpuint> indices(m_points.size());
          for (int i = 0; i < indices.size(); i++) {
               indices[i] = i;
          }
          for (int i = 0; i < m_points.size(); ++i) {
               if (stack2.top() == i) {
                    int last = i - 1;
                    while (stack2.top() == i) {
                         int size = 0;
                         hpuint walker = stack2.top();
                         chains.push_back(walker);
                         walker--;
                         while (walker != stack1.top()) {
                              if (indices[walker] == -1) {
                                   walker = indices[walker - 1];
                              }
                              chains.push_back(walker);
                              if (walker != stack1.top()) {
                                   indices[walker] = -1;
                                   walker--;
                              }
                         }
                         indices[stack2.top() - 1] = -1;
                         indices[stack2.top() - 2] = stack1.top();

                         last = stack1.top();
                         stack1.pop();
                         stack2.pop();
                    }
               }
               hpuint t = hasTwoOppositePoints(m_trapezoids[m_point_trapezoid[i]]);

               if (t != -1 && isInside(m_trapezoids[m_point_trapezoid[i]])) {
                    stack1.push(i);
                    stack2.push(t);
               }
          }
     }
     std::vector<std::tuple<hpuint, hpuint>> getDiagonalsOfTriangulatedUnimonotonePolygon(Points& points, std::vector<hpuint>& chain) {

          bool bottomMostIsZero = (m_points[chain[0]].y < m_points[chain[chain.size() - 1]].y);
          // bMiz == 0: 0 is minimal point, size-1 is topmost point
          int curr = 0;
          if (SegmentUtils::isLeftOf(points[0], points[chain.size() - 1], points[1])) {
               // kette ist links
               curr = bottomMostIsZero ? 0 : chain.size() - 1;
          } else {
               // kette ist rechts
               curr = bottomMostIsZero ? chain.size() - 1 : 0;
          }
          std::vector<std::tuple<hpuint, hpuint>> result(chain.size() - 2); // lemma ...
          std::vector<hpuint> copy(chain.size());
          for (int i = 0; i < copy.size(); ++i) {
               copy[i] = i;
          }
          int size = chain.size();
          //int add = bottomMostIsZero ? -1 : 1; // TODO Check
          hpuint next = 1, prev = 0;
          while (size > 3) {
               hpuint next = doNext(copy, curr); // TODO
               if (isConvex(copy[prev], copy[curr], copy[next]), bottomMostIsZero) {
                    // füge diagonale ein
                    result.push_back(std::tuple<hpuint, hpuint>(prev, next));
                    // delete current
                    chain[curr] = -1;
               } else {
                    prev = curr;
               }
               curr = doNext(copy, curr); // ??
          }
          return result;
     }

     void buildUpViaSeidel() {

          // zufällige Kette erzeugen
          std::vector<hpuint> randomSet;
          hpuint n = randomSet.size();
          hpuint logstarn = Trapezoidation::log_star(n);
          hpuint lastN = 0;
          for (hpuint i = 0; i < logstarn; i++) {
               hpuint nextN = N(i,n);
               for (int j = lastN; j < nextN; j++)
                    insert(randomSet[j]);
               lastN = nextN;
               trace();
          }
          for (hpuint i = logstarn; i < n; i++)
               insert(randomSet[i]);
     }


private:
     struct Node {
          // TODO are there some issus in ordering?

          unsigned data :32;
          bool leaf :1;
          unsigned leftChild :31;
          unsigned rightChild :31;
          bool x :1;

          Node(bool isLeaf, bool isXNode, hpuint data, hpuint leftChild, hpuint rightChild) :
                    leaf(isLeaf), x(isXNode), data(data), leftChild(leftChild), rightChild(rightChild) {
          }

     };
     Points m_points;

     // I assume a vector<tuple<hpuint, hpuint>>
     Iterator m_begin;
     Iterator m_end;
     hpuint m_last_node;
     std::vector<Node> m_queryStructure;
     hpuint m_last_trapezoid = 0;
     std::vector<Trapezoid> m_trapezoids;
     // not a very good way - a better way would be m_point_node
     std::vector<hpuint> m_point_trapezoid; // table representing mapping from point-index to trapez-index
     std::vector<hpuint> m_point_node;



     void trace() {
          for (int i = 0; i < m_points.size(); ++i) {
               doFind(i);
          }
     }

     void insert(hpuint segment) {
          hpuint p1 = segment; // TODO
          hpuint p2 = segment + 1; // TODO
          hpuint currentTrapez = -1;

          hpuint startPoint = p1, endPoint = p2;
          std::vector<hpuint> toSplit;

          bool directionDown = false;
          if (m_point_trapezoid[p1] != -1) {
               currentTrapez = m_point_trapezoid[p1];
          } else if (m_point_trapezoid[p2] != -1) {
               currentTrapez = m_trapezoids[p2];
               startPoint = p2;
               endPoint = p1;
          } else {
               currentTrapez = doFind(p1);
               startPoint = p1;
               endPoint = p2;

               directionDown = m_points[startPoint] > m_points[endPoint];
               hpuint newTrapezoid = splitHorizontal(currentTrapez, p1, directionDown);
               // set last trapezoid
               // get trapez where to split first

               // determine cutted trapezoids by segment
               toSplit.push_back(directionDown ? newTrapezoid : currentTrapez);
          }
          directionDown = m_points[startPoint] > m_points[endPoint];



          // can be done otherwise: if segment topdown: while (current.bottomHE > endPoint.y)
          //                            else while (current.topHE < endPoint.y
          while (directionDown ? (m_points[m_trapezoids[currentTrapez].topPoint].y < m_points[endPoint].y)
                    : (m_points[m_trapezoids[currentTrapez].bottomPoint].y > m_points[endPoint].y)) {

               currentTrapez = getNextTrapezoid(currentTrapez, segment, directionDown);

               // TODO is inside test
               //if (!isInside(currentTrapez, whereToGo))
               //   toSplit.push_back(currentTrapez);
          }

          if (m_point_trapezoid[endPoint] != -1) {
               // split also horizontal
               hpuint newTrapezoid = splitHorizontal(currentTrapez, endPoint, directionDown);
               toSplit.push_back(directionDown ?currentTrapez : newTrapezoid);
          }

          // was muss ich speichern?
          hpuint last_xnode = -1, last_leftTrapezoid = -1, last_rightTrapezoid = -1; // have to be uninited
          for (size_t i = 0; i < toSplit.size(); ++i) {
               // erster fall und letzter fall sind spezical cases.
               // split trapezoids and fix tree

               if (i == 0) {

               } else if (i == toSplit.size()-1){

               } else {
                    //Trapezoid& leftTrapez(m_queryStructure[m_last_trapezoid]);
                    Trapezoid leftTrapez(m_trapezoids[toSplit[i]]);
                    Trapezoid rightTrapez(leftTrapez);

                    leftTrapez.rightSegment = segment;
                    rightTrapez.leftSegment = segment;

                    //m_queryStructure[m_last_node] = Node(true, false, m_last_trapezoid, -1,-1);

                    // now it is cutted

                    Node& xnode = m_queryStructure[leftTrapez.leaf];

                    xnode.leaf = false;
                    xnode.x = true;

                    bool fusioned = false;;
                    if (canBeFusioned(leftTrapez, m_trapezoids[last_leftTrapezoid])) {
                         // then i got a great problem
                         // fall: linke müssen verschmolzen werden, rechts ist noch nicht in der liste
                         // => neues trapez geht auf die stelle vom alten
                         // hat das folgen für andere Punkte?
                         // -> nein sofern sie nur einen zeiger auf den node haben, dieser verweist ja immer noch richtig

                         // fusion left and last_left

                         Trapezoid& lastLeft = m_trapezoids[last_leftTrapezoid];
                         if (directionDown) {
                              // the last is upperside
                              lastLeft.bottomPoint = leftTrapez.bottomPoint;
                              lastLeft.bottomTrapezoid1 = leftTrapez.bottomTrapezoid1;
                              lastLeft.bottomTrapezoid2 = leftTrapez.bottomTrapezoid2;

                         } else {
                              lastLeft.topPoint = leftTrapez.topPoint;
                              lastLeft.topTrapezoid1 = leftTrapez.topTrapezoid1;
                              lastLeft.topTrapezoid2 = leftTrapez.topTrapezoid2;
                         }

                         rightTrapez.leaf = m_last_node;
                         m_trapezoids[toSplit[i]] = rightTrapez;
                         xnode.leftChild = lastLeft.leaf;
                         xnode.rightChild = m_last_node;
                         m_queryStructure[m_last_node] = Node(true, false, toSplit[i], -1,-1);
                         m_last_node++;


                         fusioned = true;
                    }
                    if (canBeFusioned(rightTrapez, m_trapezoids[last_rightTrapezoid])) {
                         // fall: right side will be fusioned. no trapezoid will be inserted in structure
                         // => no set of trapez.leaf
                         // => neues trapez muss nicht mehr in die liste eingefügt werden. last right wird nur aktualisiert
                         // =>

                         Trapezoid& lastRight = m_trapezoids[last_rightTrapezoid];

                         if (directionDown) {
                              lastRight.bottomPoint = rightTrapez.bottomPoint;
                              lastRight.bottomTrapezoid1 = rightTrapez.bottomTrapezoid1;
                              lastRight.bottomTrapezoid2 = rightTrapez.bottomTrapezoid2;
                         } else {
                              lastRight.topPoint = rightTrapez.topPoint;
                              lastRight.topTrapezoid1 = rightTrapez.topTrapezoid1;
                              lastRight.topTrapezoid2 = rightTrapez.topTrapezoid2;
                         }

                         m_trapezoids[toSplit[i]] = leftTrapez;
                         // no new trapezoid

                         xnode.leftChild = m_trapezoids[toSplit[i]].leaf;
                         xnode.rightChild = lastRight.leaf;

                         fusioned = true;
                    }

                    if (!fusioned) {
                         // TODO
                         // nothing to fusion, both are new
                         // case shouldn't be very often
                         Trapezoid& left = m_trapezoids[toSplit[i]];
                         Trapezoid right(left);
                         hpuint index_new;

                         left.rightSegment = segment;
                         right.leftSegment = segment;
                         Trapezoid& last_left = m_trapezoids[last_leftTrapezoid];
                         Trapezoid& last_right = m_trapezoids[last_rightTrapezoid];


                         if (directionDown) {
                              // TODO
                              // fix neighbors
                              last_left.bottomTrapezoid1 = toSplit[i];
                              last_left.bottomTrapezoid2 = -1;

                              last_right.bottomTrapezoid1 = index_new;
                              last_right.bottomTrapezoid2 = -1;

                              left.topTrapezoid1 = last_leftTrapezoid;
                              left.topTrapezoid2 = -1;

                              right.topTrapezoid1 = last_rightTrapezoid;
                              right.topTrapezoid2 = -1;
                         } else {
                              //TODO
                              // fix neighbors
                              last_left.topTrapezoid1 = toSplit[i];
                              last_left.topTrapezoid2 = -1;

                              last_right.topTrapezoid1 = index_new;
                              last_right.topTrapezoid2 = -1;
                         }




                    }


               }

               /*  what do i need for splitting vertical ? segment and trapezoid
                What do i need for fusion? trapezoids 1,2,3,4 x-node 1,2, or last x-node
                */


               // fix also horizontal endpoints
          }

     }

     bool canBeFusioned(Trapezoid& trapez1, Trapezoid& trapez2) {
          return trapez1.leftSegment == trapez2.leftSegment && trapez1.rightSegment == trapez2.rightSegment
                    && (trapez1.bottomPoint == trapez2.topPoint || trapez1.topPoint == trapez2.bottomPoint);
     }

     void fixNeighboursVertical(hpuint neighbour, hpuint oldT, hpuint newT) {
          Trapezoid& t = m_trapezoids[neighbour];
          if (t.bottomTrapezoid1 == oldT || t.bottomTrapezoid2 == oldT) {
               t.bottomTrapezoid1 = oldT;
               t.bottomTrapezoid2 = newT;
          } else {
               t.topTrapezoid1 == oldT;
               t.bottomTrapezoid2 = newT;
          }

     }

     void fixNeighboursHorizontal(hpuint neighbour, hpuint oldT, hpuint newT) {
          Trapezoid& tau = m_trapezoids[neighbour];
          if (tau.topTrapezoid1 == oldT)
               tau.topTrapezoid1 = newT;
          else
               tau.topTrapezoid2 = newT;
     }

     /* Returns the new trapezoid */
     hpuint splitHorizontal(hpuint trapez, hpuint point, bool lineDirectionDown) {
          Trapezoid& t = m_trapezoids[trapez];
          Node& node = m_queryStructure[t.leaf];

          Trapezoid& newTrapez(t);
          hpuint indexNewTrapez = m_last_trapezoid;
          // these lines of code will fix the problem of isExtensionDefinedBySegment
          if (lineDirectionDown)
               t.bottomPoint = point;
          else
               t.topPoint = point;

          t.bottomTrapezoid1 = indexNewTrapez;
          t.bottomTrapezoid2 = -1;

          newTrapez.topTrapezoid1 = trapez;
          newTrapez.topTrapezoid1 = -1;

          // fix neighbour
          fixNeighboursHorizontal(newTrapez.bottomTrapezoid1, trapez, indexNewTrapez);
          fixNeighboursHorizontal(newTrapez.bottomTrapezoid2, trapez, indexNewTrapez);

          m_trapezoids[indexNewTrapez] = newTrapez; // oder m_last_trapezoid++
          m_last_trapezoid++;

          // build up query-structure
          node.x = false;
          node.leaf = false;
          // push into list
          m_queryStructure[m_last_node] = Node(true, false, trapez, -1, -1);
          m_queryStructure[m_last_node]= Node(true, false, newTrapez, -1, -1);

          // fix childs
          node.leftChild = m_last_node;
          node.rightChild = m_last_node+1;
          m_last_node += 2;
          return indexNewTrapez;
     }

     /**
      * Return -1 if false else the point
      *
      */
     hpuint isExtensionDefinedBySegment(hpuint trapezoid, bool top) {
          Trapezoid& t = m_trapezoids[trapezoid];
          // if we assume there is no bad case
          return top ? t.topPoint : t.bottomPoint;
          // if we assume there can be bad cases - check - there will be no :) - if splitHorizontal & fusion is correct implemented
     }

     hpuint getNextTrapezoid(hpuint trapez, hpuint segment, bool directionDown) {
          // other implementierung:
          // only if there is an HE defined by a horizontal extension of third segment there is a branch, otherwise there is only one trapezoid
          hpuint pointOfHE = isExtensionDefinedBySegment(trapez, directionDown);

          if (pointOfHE != -1) {
               if (directionDown)
                    return SegmentUtils::isLeftOf(getSegment(segment), getPoint(pointOfHE)) ? m_trapezoids[trapez].bottomTrapezoid1 : m_trapezoids[trapez].bottomTrapezoid2;
               else
                    return SegmentUtils::isLeftOf(getSegment(segment), getPoint(pointOfHE)) ? m_trapezoids[trapez].topTrapezoid1 : m_trapezoids[trapez].topTrapezoid2;

          } else {
               if (directionDown)
                    return m_trapezoids[trapez].bottomTrapezoid1 != -1 ? m_trapezoids[trapez].bottomTrapezoid1 : m_trapezoids[trapez].bottomTrapezoid2;
               else
                    return m_trapezoids[trapez].topTrapezoid1 != -1 ? m_trapezoids[trapez].topTrapezoid1 : m_trapezoids[trapez].topTrapezoid2;
          }
          /* OLD CODE
           Trapezoid& t = m_trapezoids[trapez];
           if (t.bottomTrapezoid1 != lastTrapez && cross(t.bottomTrapezoid1, segment)) return t.bottomTrapezoid1;
           else if (t.bottomTrapezoid2 != lastTrapez && cross(t.bottomTrapezoid2, segment)) return t.bottomTrapezoid2;
           else if (t.topTrapezoid1 != lastTrapez && cross(t.topTrapezoid1, segment)) return t.topTrapezoid1;
           else if (t.topTrapezoid2 != lastTrapez && cross(t.topTrapezoid2, segment)) return t.topTrapezoid2;
           return -1; // error
           */
     }

     hpuint getMin(hpuint p1, hpuint p2) {
          return getPoint(p1).y < getPoint(p2).y ? p1 : p2;
     }
     hpuint getMax(hpuint p1, hpuint p2) {
          return getPoint(p1).y < getPoint(p2).y ? p1 : p2;
     }

     std::tuple<Point2D&, hpuint> getCut(hpuint i_seg1, hpuint i_seg2, bool top) {
          hpuint p = top ? getUpperHorizontalExtension(i_seg1, i_seg2) : getBottomHorizontalExtension(i_seg1, i_seg2);
          // TODO is wrong
          if (p & (i_seg1 >> 1) == i_seg1)
               return std::tuple<Point2D&, hpuint>(SegmentUtils::intersect(getSegment(i_seg2), getPoint(p)), p); // p is in seg1
          else
               return std::tuple<Point2D&, hpuint>(SegmentUtils::intersect(getSegment(i_seg1), getPoint(p)), p); // p is in seg2
     }
     bool isConvex(hpuint p1, hpuint p2, hpuint p3, bool oriented) {
          return oriented ? glm::dot(glm::normalize(m_points[p1] - m_points[p2]), glm::normalize(m_points[p3] - m_points[p2])) > 0 : isConvex(p3, p2, p1, true);
     }

     Point2D& getPoint(hpuint index) {
          return m_points[index];
          /* OLD VERSION
          std::tuple<Point2D&, Point2D&> ps = *(m_begin[(index >> 1)]);
          if (index & 1)
               return std::get < 1 > (ps);
          else
               return std::get < 0 > (ps);
               */
     }

     std::tuple<Point2D&, Point2D&> getSegment(hpuint index) {
          return tuple<Point2D, Point2D&>(m_points[std::get<0>(m_begin[index])], m_points[std::get<1>(m_begin[index])]);
          //return (m_begin[(index >> 1)]); // OLD VERSION
     }

     hpuint& doFind(hpuint point) {
          hpuint node = 0;
          if (0 <= point && point < m_point_trapezoid.size() && m_point_trapezoid[point] != -1) {
               node = m_trapezoids[m_point_trapezoid[point]]->leaf;
          }
          while (!m_queryStructure[node].leaf) {
               if (!m_queryStructure[node].x) {
                    if (SegmentUtils::isLeftOf(getSegment(m_queryStructure[node].data), m_points[point]))
                         node = m_queryStructure[node].leftChild;
                    else
                         node = m_queryStructure[node].rightChild;
               } else if (getPoint(m_queryStructure[node].data).y < m_points[point].y)
                    node = m_queryStructure[node].leftChild;
               else
                    node = m_queryStructure[node].rightChild;
          }
          m_point_trapezoid[point] = m_queryStructure[node].data;
          return m_point_trapezoid[point];
     }

     Trapezoid& doFind(const Point2D& point) {
          Node& node = m_queryStructure[0];
          while (!node.leaf) {
               if (node.x) {
                    std::tuple<Point2D&, Point2D&> segment = getSegment(node.data);
                    if (SegmentUtils::isLeftOf(segment, point))
                         node = m_queryStructure[node.leftChild];
                    else
                         node = m_queryStructure[node.rightChild];
               } else {
                    Point2D& p = getPoint(node.data);
                    if (p.y < point.y)
                         node = m_queryStructure[node.leftChild];
                    else
                         node = m_queryStructure[node.rightChild];
               }
          }
          return m_trapezoids[node.data];
     }

     hpuint doNext(std::vector<hpuint> chain, hpuint curr) {
          // TODO wrong indicies
          int size = chain.size();
          return chain[(curr + 1 + size) % size] == -1 ? chain[(curr + 2 + size) % size] : (curr + 1 + size) % size;
     }
     hpuint doPrev(std::vector<hpuint> chain, hpuint curr) {
          // TODO wrong indicies
          int size = chain.size();
          int prev = (curr - 1 + size) % size;
          return chain[prev] == -1 ? chain[(prev - 2 + size) % size] : chain[prev];
     }

     bool isUpsideDown(hpuint segment) {
          std::tuple<Point2D&, Point2D&> seg = getSegment(segment);
          return (std::get < 0 > (seg).y > std::get < 1 > (seg).y);
     }

     bool isInsideOfPolygon(hpuint trapezoid) { // in polygon
          // idea the right segment must be upsidedown and the left segment must be bottomup
          hpuint rightSeg = m_trapezoids[trapezoid].rightSegment, leftSeg = m_trapezoids[trapezoid].leftSegment;
          return (rightSeg == -1 || isUpsideDown(rightSeg)) && (leftSeg == -1 || !isUpsideDown(leftSeg));
     }

// Is point in trapezoid
     bool isInsideOfTrapezoid(hpuint trapezoid, hpuint point) {
          return (isRelativCounterClockwise(getSegment(m_trapezoids[trapezoid].leftSegment)) && isRelativCounterClockwise(getSegment(m_trapezoids[trapezoid].rightSegment)));
     }

// Only use if u know what u do
     hpuint getUpperHorizontalExtension(hpuint i_seg1, hpuint i_seg2) {
          return getMin(getMax(std::get < 0 > (getSegment(i_seg1)), std::get < 0 > (getSegment(i_seg1))), getMax(std::get < 0 > (getSegment(i_seg2)), std::get < 1 > (getSegment(i_seg2))));
     }
     hpuint getBottomHorizontalExtension(hpuint i_seg1, hpuint i_seg2) {
          return getMax(getMix(std::get < 0 > (getSegment(i_seg1)), std::get < 0 > (getSegment(i_seg1))), getMix(std::get < 0 > (getSegment(i_seg2)), std::get < 1 > (getSegment(i_seg2))));
     }
     hpuint hasTwoOppositePoints(Trapezoid& t) {
          hpuint u = getUpperHorizontalExtension(t.leftSegment, t.rightSegment), b = getBottomHorizontalExtension(t.leftSegment, t.rightSegment);
          if (t.leftSegment == -1 || t.rightSegment != -1)
               return -1;
          if (t.topPoint != -1) {
               return max(t.topPoint, b);
          }
          if (t.bottomPoint != -1) {
               return max(t.bottomPoint, u);
          }
          return (b >> 1) != (u >> 1) ? max(b, u) : -1;
     }
     static hpuint log_star(hpuint n) {
          return ackermannInverse(2, n);
     }
     static hpuint ackermannInverse(hpuint k, hpuint n) {
          // view: take position of first digit of hpuint that is not null, log_2 position of that digit
          if (k <= 1) {
               return n >> 1;
          } else if (k >= 2) {
               if (n == 1)
                    return 0;
               else
                    return 1 + ackermannInverse(k, ackermannInverse(k - 1, n));
          } else {
               return 0;
          }
     }
     static hpuint logh(hpuint h, hpuint n) {
          return h <= 1 ? std::log2(n) : std::log2(logh(h-1,n));
     }
     static hpuint N(hpuint h, hpuint n){
          return ceil(n/logh(h,n));
     }
};
