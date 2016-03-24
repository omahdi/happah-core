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

class ReaderUtils {
public:
     struct ReadException : public std::exception {};

     template<class Iterator>
     static void consume(Iterator& begin, const Iterator& end, const char* literal) {
          if(begin == end) throw ReadException();
          do {
               if(*begin != *literal) throw ReadException();
               ++begin;
               ++literal;
          } while(begin != end && *literal != '\0');
          if(*literal != '\0') throw ReadException();
     }

     template<class Iterator>
     static void consumeBlanksPlus(Iterator& begin, const Iterator& end) {
          if(begin == end || !isblank(*begin)) throw ReadException();
          do ++begin;
          while(begin != end && isblank(*begin));
     }

     template<class Iterator>
     static void consumeBlanksStar(Iterator& begin, const Iterator& end) { while(begin != end && isblank(*begin)) ++begin; }

     template<class Iterator>
     static void consumeLine(Iterator& begin, const Iterator& end) {
          if(begin == end) return;
          while(*begin != '\n') {
               ++begin;
               if(begin == end) return;
          }
          ++begin;
     }

     template<class Iterator>
     static void consumeMacro(Iterator& begin, const Iterator& end) {
          while(begin < end) {//TODO: fix
               if(*begin == '\\') {//TODO: if inside string literal ignore
                    ++begin;
                    if(begin < end && *begin == '\n') ++begin;
                    else throw ReadException();
               } else if(*begin == '\n') {
                    ++begin;
                    return;
               } else ++begin;
          }
     };

     template<class Iterator>
     static void consumeWhitespacePlus(Iterator& begin, const Iterator& end) {
          if(begin == end || !isspace(*begin)) throw ReadException();
          do ++begin;
          while(begin != end && isspace(*begin));
     }

     template<class Iterator>
     static void consumeWhitespaceStar(Iterator& begin, const Iterator& end) { while(begin != end && isspace(*begin)) ++begin; }

     template<typename Decimal, class Iterator>
     static Decimal readDecimal(Iterator& begin, const Iterator& end) {
          if(begin == end) throw ReadException();
          Decimal r = 0.0;
          bool negate = false;
          if(*begin == '-') {
               negate = true;
               ++begin;
               if(begin == end) throw ReadException();
          }
          if(*begin == '.') {
               ++begin;
               readMantissa(begin, end, r);
          } else if(*begin >= '0' && *begin <= '9') {
               do {
                    r = (r * 10.0) + (*begin - '0');
                    ++begin;
               } while(begin != end && *begin >= '0' && *begin <= '9');
               if(*begin == '.') {
                    ++begin;
                    readMantissa(begin, end, r);
               } else if(*begin == 'e' || *begin == 'E') {
                    ++begin;
                    r *= std::pow(10.0, readInt(begin, end));
               }
          } else throw ReadException();
          return (negate) ? -r : r;
     }

     template<class Iterator>
     static double readDouble(Iterator& begin, const Iterator& end) { return readDecimal<double>(begin, end); }

     template<class Iterator>
     static float readFloat(Iterator& begin, const Iterator& end) { return readDecimal<float>(begin, end); }

     template<class Iterator>
     static int readInt(Iterator& begin, const Iterator& end) {
          if(begin == end) throw ReadException();
          bool negate = false;
          if(*begin == '-') {
               negate = true;
               ++begin;
               if(begin == end) throw ReadException();
          } else if(*begin == '+') {
               ++begin;
               if(begin == end) throw ReadException();
          }
          if(*begin >= '0' && *begin <= '9') {
               int n = 0;
               do {
                    n = (n * 10) + (*begin - '0');
                    ++begin;
               } while(begin != end && *begin >= '0' && *begin <= '9');
               return (negate) ? -n : n;
          } else throw ReadException();
     }

     template<class Iterator>
     static unsigned int readUnsignedInt(Iterator& begin, const Iterator& end) {
          if(begin == end) throw ReadException();
          if(*begin == '+') {
               ++begin;
               if(begin == end) throw ReadException();
          }
          if(*begin >= '0' && *begin <= '9') {
               unsigned int n = 0;
               do {
                    n = (n * 10) + (*begin - '0');
                    ++begin;
               } while(begin != end && *begin >= '0' && *begin <= '9');
               return n;
          } else throw ReadException();
     }

private:
     template<typename Decimal, class Iterator>
     static void readMantissa(Iterator& begin, const Iterator& end, Decimal& r) {
          if(begin == end) throw ReadException();
          Decimal f = 0.0;
          int n = 0;
          while (*begin >= '0' && *begin <= '9') {
               f = (f * 10.0) + (*begin - '0');
               ++begin;
               ++n;
          }
          r += f / std::pow(10.0, n);
          if(*begin == 'e' || *begin == 'E') {
               ++begin;
               r *= std::pow(10.0, readInt(begin, end));
          }
     }

};

