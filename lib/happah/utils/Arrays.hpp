// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <type_traits>
#include <vector>

#include "happah/Happah.hpp"

template<class T>
class Arrays {
     using Data = std::vector<T>;
     using Indices = std::vector<hpuint>;

     template<bool t_const = false>
     class Iterator {
     public:
          using DataIterator = typename std::conditional<t_const, typename Data::const_iterator, typename Data::iterator>::type;
          using IndicesIterator = typename std::conditional<t_const, typename Indices::const_iterator, typename Indices::iterator>::type;
          using iterator = DataIterator;
          using value_type = std::pair<DataIterator, DataIterator>;

          Iterator(DataIterator array, IndicesIterator length)
               : m_array(array), m_length(length) {}

          iterator begin() const { return m_array; }

          iterator end() const { return m_array + *m_length; }

          Iterator operator+(hpuint n) const {
               IndicesIterator end = m_length + n;
               DataIterator array = m_array;
               for(IndicesIterator i = m_length; i != end; ++i)
                    array += *i;
               return Iterator(array, end); 
          }

          Iterator operator-(hpuint n) const {
               IndicesIterator end = m_length - n;
               DataIterator array = m_array;
               IndicesIterator i = m_length;
               do {
                    --i;
                    array -= *i;
               } while(i != end);
               return Iterator(array, end); 
          }

          Iterator& operator++() {
               m_array += *m_length;
               ++m_length;
               return *this;
          }

          Iterator& operator--() {
               --m_length;
               m_array -= *m_length;
               return *this;
          }

          Iterator operator++(int) const { 
               Iterator iterator(*this);
               return ++iterator;
          }

          Iterator operator--(int) const { 
               Iterator iterator(*this);
               return --iterator;
          }

          bool operator!=(const Iterator& iterator) const { return iterator.m_length != m_length; }

          bool operator==(const Iterator& iterator) const { return iterator.m_length == m_length; }

          value_type operator*() const { return std::make_pair(begin(), end()); }

     private:
          DataIterator m_array;
          IndicesIterator m_length;

          friend class Arrays<T>;

     };

public:
     class ArrayAppender {
     public:
          using value_type = typename Data::value_type;

          ArrayAppender(Arrays& arrays)
               : m_arrays(arrays.m_arrays), m_lengths(arrays.m_lengths), m_size(m_arrays.size()) {}

          ~ArrayAppender() {
               hpuint length = m_arrays.size() - m_size;
               if(length > 0) m_lengths.push_back(length);
          }

          void push_back(const T& t) { m_arrays.push_back(t); }

          void push_back(T&& t) { m_arrays.push_back(t); }

          void reserve(hpuint n) { m_arrays.reserve(m_arrays.size() + n); }

          void resize(hpuint n, const T& t) { m_arrays.resize(n, t); }//TODO: m_arrays.size() + n!

          T& operator[](hpuint n) { return m_arrays[n]; }

     private:
          Data& m_arrays;
          Indices& m_lengths;
          hpuint m_size;

     };//ArrayAppender

     using const_iterator = Iterator<true>;
     using iterator = Iterator<false>;

     iterator begin() { return { m_arrays.begin(), m_lengths.begin() }; }

     const_iterator begin() const { return { m_arrays.cbegin(), m_lengths.cbegin() }; }

     const_iterator cbegin() const { return { m_arrays.cbegin(), m_lengths.cbegin() }; }

     const_iterator cend() const { return { m_arrays.cend(), m_lengths.cend() }; }

     const Data& data() const { return m_arrays; }
     
     Data& data() { return m_arrays; }

     iterator end() { return { m_arrays.end(), m_lengths.end() }; }

     const_iterator end() const { return { m_arrays.cend(), m_lengths.cend() }; }

     void erase(iterator i) {
          m_arrays.erase(i.begin(), i.end());
          m_lengths.erase(i.m_length);
     }

     void erase(iterator first, iterator last) {
          m_arrays.erase(first.begin(), last.begin());
          m_lengths.erase(first.m_length, last.m_length);
     }

     template<class InputIterator>
     void insert(iterator i, InputIterator first, InputIterator last) {
          hpuint length = std::distance(first, last);
          m_arrays.insert(i.begin(), first, last);
          m_lengths.insert(i.m_length, length);
     }

     template<class InputIterator>
     void push_back(InputIterator first, InputIterator last) { insert(end(), first, last); }

     void reserve(hpuint nArrays) { m_lengths.reserve(nArrays); }

     void reserve(hpuint nArrays, hpuint nData) {
          m_arrays.reserve(nData);
          m_lengths.reserve(nArrays); 
     }

     hpuint size() const { return m_lengths.size(); }

     typename Iterator<false>::value_type operator[](hpuint i) { return *(begin() + i); }

     typename Iterator<true>::value_type operator[](hpuint i) const { return *(cbegin() + i); }

     template<class Stream>
     friend Stream& operator<<(Stream& stream, const Arrays& arrays) {
          stream << arrays.m_arrays;
          stream << arrays.m_lengths;
          return stream;
     }

     template<class Stream>
     friend Stream& operator>>(Stream& stream, Arrays& arrays) {
          stream >> arrays.m_arrays;
          stream >> arrays.m_lengths;
          return stream;
     }

private:
     Data m_arrays;
     Indices m_lengths;

};//Arrays

