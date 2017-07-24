// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <iterator>

#include "happah/Happah.hpp"

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

          const value_type& operator[](hpuint offset) const { return m_data[m_i[offset]]; }

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

     auto& operator[](hpuint offset) const { return m_data[m_indices[offset]]; }

private:
     const Data& m_data;
     const Indices& m_indices;

};//DeindexedArray

template<class Data>
DeindexedArray<Data> deindex(const Data& data, const Indices& indices) { return { data, indices }; }

template<class Data>
DeindexedArray<Data> deindex(const std::tuple<const Data&, const Indices&>& array) { return { std::get<0>(array), std::get<1>(array) }; }

}//namespace happah

//**********************************************************************************************************************************
//TODO: cleanup below
//**********************************************************************************************************************************

#include "happah/math/Space.hpp"

namespace std {

template<>
struct iterator_traits<typename happah::DeindexedArray<std::vector<happah::Point3D> >::const_iterator> {
     using iterator_category = std::bidirectional_iterator_tag;
     using value_type = happah::Point3D;
     using difference_type = typename happah::Indices::const_iterator::difference_type;
     using pointer = happah::Point3D*;
     using reference = happah::Point3D&;

};

template<>
struct iterator_traits<typename happah::DeindexedArray<std::vector<happah::Point2D> >::const_iterator> {
     using iterator_category = std::bidirectional_iterator_tag;
     using value_type = happah::Point2D;
     using difference_type = typename happah::Indices::const_iterator::difference_type;
     using pointer = happah::Point2D*;
     using reference = happah::Point2D&;

};

template<>
struct iterator_traits<typename happah::DeindexedArray<std::vector<happah::Point1D> >::const_iterator> {
     using iterator_category = std::bidirectional_iterator_tag;
     using value_type = happah::Point1D;
     using difference_type = typename happah::Indices::const_iterator::difference_type;
     using pointer = happah::Point1D*;
     using reference = happah::Point1D&;

};

template<>
struct iterator_traits<typename happah::DeindexedArray<std::vector<happah::VertexP3> >::const_iterator> {
     using iterator_category = std::bidirectional_iterator_tag;
     using value_type = happah::VertexP3;
     using difference_type = typename happah::Indices::const_iterator::difference_type;
     using pointer = happah::VertexP3*;
     using reference = happah::VertexP3&;

};

}//namespace std

