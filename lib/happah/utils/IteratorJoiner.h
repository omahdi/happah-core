// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <vector>

template<class Iterator>
class IteratorJoiner {
     class ProxyIterator {
     public:
          typedef typename Iterator::difference_type difference_type;
          typedef std::input_iterator_tag iterator_category;
          typedef typename Iterator::pointer pointer;
          typedef typename Iterator::reference reference;
          typedef typename Iterator::value_type value_type;

          static difference_type distance(ProxyIterator& begin, ProxyIterator& end) {
               if(begin.m_iterator == end.m_iterator) return std::distance(begin.m_i, end.m_i);
               difference_type result = std::distance(begin.m_i, begin.m_end);
               auto i = begin.m_iterator + 2;
               while(i != end.m_iterator) {
                    result += std::distance(begin.m_i, begin.m_end);
                    i += 2;
               }
               result += std::distance(i.m_begin, end.m_i);
               return result;
          }

          ProxyIterator(std::vector<Iterator>& iterators, typename std::vector<Iterator>::iterator iterator, Iterator i)
               : m_begin(*iterator), m_end(*(++iterator)), m_i(i), m_iterator(++iterator), m_iterators(iterators) {}

          ProxyIterator(const ProxyIterator& iterator)
               : m_begin(iterator.m_begin), m_end(iterator.m_end), m_i(iterator.m_i), m_iterator(iterator.m_iterator), m_iterators(iterator.m_iterators) {}

          ProxyIterator& operator++() {
               ++m_i;
               if(m_i == m_end && m_iterator != m_iterators.end()) {
                    m_begin = *m_iterator;
                    m_i = m_begin;
                    m_end = *(++m_iterator);
                    ++m_iterator;
               }
               return *this;
          }

          ProxyIterator& operator--() {
               if(m_i == m_begin) {
                    m_iterator -= 2;
                    if(m_iterator != m_iterators.begin()) {
                         m_end = *(m_iterator-1);
                         m_begin = *(m_end-1);
                         m_i = m_end;
                    }
               }
               --m_i;
               return *this;
          }

          bool operator!=(const ProxyIterator& iterator) { return iterator.m_i != m_i; }

          bool operator==(const ProxyIterator& iterator) { return iterator.m_i == m_i; }

          value_type operator*() { return *m_i; }

     private:
          Iterator m_begin;
          Iterator m_end;
          Iterator m_i;
          typename std::vector<Iterator>::iterator m_iterator;
          std::vector<Iterator>& m_iterators;
          
     };

public:
     typedef ProxyIterator iterator;

     iterator begin() { return ProxyIterator(m_iterators, m_iterators.begin(), *(m_iterators.begin())); }

     iterator end() { return ProxyIterator(m_iterators, m_iterators.end()-2, m_iterators.back()); }

     void pop_back() {
          m_iterators.pop_back();
          m_iterators.pop_back();
     }

     void push_back(Iterator begin, Iterator end) {
          m_iterators.push_back(begin);
          m_iterators.push_back(end);
     }

private:
     std::vector<Iterator> m_iterators;

};

namespace std {

template<class Iterator>
typename IteratorJoiner<Iterator>::iterator::difference_type distance(typename IteratorJoiner<Iterator>::iterator& begin, typename IteratorJoiner<Iterator>::iterator& end) { return typename IteratorJoiner<Iterator>::iterator::distance(begin, end); }

}
