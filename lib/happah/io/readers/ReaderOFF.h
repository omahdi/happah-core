/**********************************************************************************
 * Copyright (c) 2012-2015  See the COPYRIGHT-HOLDERS file in the top-level
 * directory of this distribution and at http://github.com/happah-graphics/happah.
 *
 * This file is part of Happah. It is subject to the license terms in the LICENSE
 * file found in the top-level directory of this distribution and at
 * http://github.com/happah-graphics/happah. No part of Happah, including this
 * file, may be copied, modified, propagated, or distributed except according to
 * the terms contained in the LICENSE file.
 **********************************************************************************/

#pragma once

#include <boost/optional.hpp>
#include <vector>

#include "happah/Happah.h"
#include "happah/geometries/TriangleMesh.h"
#include "happah/io/readers/ReaderUtils.h"
#include "happah/utils/VertexFactory.h"

namespace happah {

//TODO: nice error messages/codes to explain why parsing failed; have ReadException inherit from runtime_error and pass message
//TODO: accomodate normal, color, texture
class ReaderOFF {
     struct Header {
          bool color;
          hpuint dimension;
          hpuint nEdges;
          hpuint nFaces;
          bool normal;
          hpuint nVertices;
          bool texture;

     };//Header

public:
     using Indices = std::vector<hpuint>;
     using ReadException = ReaderUtils::ReadException;
     using Reals = std::vector<hpreal>;

     struct Output {
          hpuint dimension;
          Indices faces;
          boost::optional<Indices> lengths;
          hpuint length;
          Reals vertices;

     };//Output

     template<class Iterator>
     static Output read(Iterator begin, Iterator end) {
          Output output;
          auto header = getHeader(begin, end);
          if(header.dimension == 0) throw ReadException();
          output.dimension = header.dimension;
          output.vertices = getVertices(begin, end, header);
          std::tie(output.faces, output.lengths, output.length) = getFaces(begin, end, header);
          return std::move(output);
     }

     static Output read(const char* path);

     template<class Vertex = VertexP3, class VertexFactory = VertexFactory<Vertex> >
     static TriangleMesh<Vertex> toTriangleMesh(const Output& output, VertexFactory&& factory = VertexFactory()) {
          if(output.length != 3) throw MissingImplementationException();
          std::vector<Vertex> vertices;
          for(auto i = output.vertices.cbegin(), end = output.vertices.cend(); i != end; ++i) {
               auto i0 = *i;
               auto i1 = *(++i);
               auto i2 = *(++i);
               vertices.push_back(factory({i0, i1, i2}));
          }
          return TriangleMesh<Vertex>(std::move(vertices), output.faces);
     }

private:
     template<class Iterator>
     static Header getHeader(Iterator& begin, const Iterator& end) {
          Header header;
          ReaderUtils::consumeWhitespaceStar(begin, end);
          if(begin != end && *begin == 'S') {
               ++begin;
               if(begin != end && *begin == 'T') {
                    ++begin;
                    header.texture = true;
               } else throw ReadException();
          }
          if(begin != end && *begin == 'C') {
               ++begin;
               header.color = true;
          }
          if(begin != end && *begin == 'N') {
               ++begin;
               header.normal = true;
          }
          bool n = false;
          if(begin != end) { 
               if(*begin == '4') {
                    ++begin;//TODO: 4n
                    header.dimension = 4;
               } else if(*begin == 'n') {
                    ++begin;
                    n = true;
               } else header.dimension = 3;
          } else throw ReadException();
          ReaderUtils::consume(begin, end, "OFF");
          ReaderUtils::consumeLine(begin, end);
          if(n) {
               ReaderUtils::consumeWhitespacePlus(begin, end);
               header.dimension = ReaderUtils::readUnsignedInt(begin, end);
               ReaderUtils::consumeLine(begin, end);
          }
          ReaderUtils::consumeWhitespaceStar(begin, end);
          header.nVertices = ReaderUtils::readUnsignedInt(begin, end);
          ReaderUtils::consumeBlanksPlus(begin, end);
          header.nFaces = ReaderUtils::readUnsignedInt(begin, end);
          ReaderUtils::consumeBlanksPlus(begin, end); 
          header.nEdges = ReaderUtils::readUnsignedInt(begin, end);
          ReaderUtils::consumeLine(begin, end);
          return std::move(header);
     }

     template<class Iterator>
     static Reals getVertices(Iterator& begin, const Iterator& end, const Header& header) {
          Reals vertices;
          vertices.reserve(header.nVertices * header.dimension);
          for(hpuint i = 0; i < header.nVertices; ++i) {
               ReaderUtils::consumeWhitespaceStar(begin, end);
               for(hpuint j = 0; j < header.dimension - 1; ++j) {
                    vertices.push_back(ReaderUtils::readDecimal<hpreal>(begin, end));
                    ReaderUtils::consumeBlanksPlus(begin, end);
               }
               vertices.push_back(ReaderUtils::readDecimal<hpreal>(begin, end));
               ReaderUtils::consumeLine(begin, end);
          }
          return std::move(vertices);
     }

     template<class Iterator>
     static std::tuple<Indices, boost::optional<Indices>, hpuint> getFaces(Iterator& begin, const Iterator& end, const Header& header) {
          if(header.nFaces == 0) throw ReadException();
          ReaderUtils::consumeWhitespaceStar(begin, end);
          hpuint nFirstFaceIndices = ReaderUtils::readUnsignedInt(begin, end);
          Indices faces;
          for(hpuint i = 0; i < nFirstFaceIndices; ++i) {
               ReaderUtils::consumeBlanksPlus(begin, end);
               faces.push_back(ReaderUtils::readUnsignedInt(begin, end));
          }
          Indices lengths;
          hpuint i = 0;
          bool same = true;
          while((++i) < header.nFaces) {
               ReaderUtils::consumeWhitespaceStar(begin, end);
               hpuint nIndices = ReaderUtils::readUnsignedInt(begin, end);
               for(hpuint j = 0; j < nIndices; ++j) {
                    ReaderUtils::consumeBlanksPlus(begin, end);
                    faces.push_back(ReaderUtils::readUnsignedInt(begin, end));
               }
               ReaderUtils::consumeLine(begin, end);
               if(nIndices != nFirstFaceIndices) {
                    lengths.reserve(header.nFaces);
                    lengths.insert(lengths.end(), i, nFirstFaceIndices);
                    lengths.push_back(nIndices);
                    same = false;
                    break;
               }
          }
          while((++i) < header.nFaces) {
               ReaderUtils::consumeWhitespaceStar(begin, end);
               hpuint nIndices = ReaderUtils::readUnsignedInt(begin, end);
               for(hpuint j = 0; j < nIndices; ++j) {
                    ReaderUtils::consumeBlanksPlus(begin, end);
                    faces.push_back(ReaderUtils::readUnsignedInt(begin, end));
               }
               ReaderUtils::consumeLine(begin, end);
               lengths.push_back(nIndices);
          }
          if(same) {
               boost::optional<Indices> temp = boost::none;
               return std::make_tuple(std::move(faces), temp, nFirstFaceIndices);
          } else {
               boost::optional<Indices> temp = lengths;
               return std::make_tuple(std::move(faces), temp, UNULL);
          }
     }

};//ReaderOFF

}//namespace happah

