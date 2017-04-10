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

#include "happah/Happah.h"
#include "happah/readers/ReaderUtils.h"
#include "happah/utils/Arrays.h"

//TODO: for flat triangle mesh shading compute normals in geometry shader
//TODO: nice error messages/codes to explain why parsing failed
class ReaderLOP {
public:
     using IndicesArrays = Arrays<hpuint>;
     using ReadException = ReaderUtils::ReadException;

     template<class Iterator>
     static IndicesArrays read(Iterator begin, Iterator end) {
          IndicesArrays arrays;
          ReaderUtils::consumeWhitespaceStar(begin, end);
          ReaderUtils::consume(begin, end, "LOP");
          ReaderUtils::consumeWhitespacePlus(begin, end);
          hpuint genus = ReaderUtils::readUnsignedInt(begin, end);
          arrays.reserve(genus << 1);
          while(begin != end && genus > 0) {
               ReaderUtils::consumeWhitespacePlus(begin, end);
               getArray(begin, end, IndicesArrays::ArrayAppender(arrays));
               ReaderUtils::consumeWhitespacePlus(begin, end);
               getArray(begin, end, IndicesArrays::ArrayAppender(arrays));
               --genus;     
          }
          if(genus != 0) throw ReadException();
          return std::move(arrays);
     }

     static IndicesArrays read(const char* path);

private:
     template<class Iterator>
     static void getArray(Iterator& begin, const Iterator& end, IndicesArrays::ArrayAppender appender) {
          hpuint nVertices = ReaderUtils::readUnsignedInt(begin, end);
          appender.reserve(nVertices);
          while(begin != end && nVertices > 0) {
               ReaderUtils::consumeWhitespacePlus(begin, end);
               hpuint index = ReaderUtils::readUnsignedInt(begin, end);
               appender.push_back(index);
               --nVertices;
          }
          if(nVertices != 0) throw ReadException();
     }

};

