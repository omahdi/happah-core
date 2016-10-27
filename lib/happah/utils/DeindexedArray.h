// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <iterator>

#include "happah/Happah.h"

namespace happah {

template<class Data>
class DeindexedArray {
     class Iterator : std::iterator<std::bidirectional_iterator_tag, typename Data::const_iterator::value_type, typename Indices::const_iterator::difference_type> {
     public:
          using difference_type = typename Indices::const_iterator::difference_type;
          using value_type = typename Data::const_iterator::value_type;

          Iterator(const DeindexedArray& array, hpuint offset) 
               : m_data(array.m_data), m_i(array.m_indices.cbegin() + offset) {}

          difference_type operator-(const Iterator& iterator) const { return m_i - iterator.m_i; }

          Iterator operator+(hpuint offset) const {
               Iterator iterator(*this);
               iterator += offset;
               return iterator;
          }

          Iterator& operator++() { 
               ++m_i;
               return *this;
          }

          Iterator& operator--() { 
               --m_i;
               return *this;
          }

          Iterator& operator+=(hpuint offset) {
               m_i += offset;
               return *this;
          }

          Iterator& operator-=(hpuint offset) {
               m_i -= offset;
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

          bool operator!=(const Iterator& iterator) const { return iterator.m_i != m_i; }

          const value_type& operator*() const { return m_data[*m_i]; }

     private:
          const Data& m_data;
          typename Indices::const_iterator m_i;

     };//Iterator

public:
     using const_iterator = Iterator;

     DeindexedArray(const Data& data, const Indices& indices)
          : m_data(data), m_indices(indices) {}

     const_iterator begin() const { return Iterator(*this, 0); }

     const_iterator cbegin() const { return begin(); }

     const_iterator cend() const { return end(); }

     const_iterator end() const { return Iterator(*this, m_indices.size()); }

private:
     const Data& m_data;
     const Indices& m_indices;

};//DeindexedArray

template<class Data>
DeindexedArray<Data> deindex(const Data& data, const Indices& indices) { return { data, indices }; }

}//namespace happah

//**********************************************************************************************************************************
//TODO: cleanup below
//**********************************************************************************************************************************

#include "happah/math/Space.h"

namespace std {

template<>
struct iterator_traits<typename happah::DeindexedArray<std::vector<Point3D> >::const_iterator> {
     using iterator_category = std::bidirectional_iterator_tag;
     using value_type = Point3D;
     using difference_type = typename happah::Indices::const_iterator::difference_type;
     using pointer = Point3D*;
     using reference = Point3D&;

};

template<>
struct iterator_traits<typename happah::DeindexedArray<std::vector<Point1D> >::const_iterator> {
     using iterator_category = std::bidirectional_iterator_tag;
     using value_type = Point1D;
     using difference_type = typename happah::Indices::const_iterator::difference_type;
     using pointer = Point1D*;
     using reference = Point1D&;

};

}//namespace std

