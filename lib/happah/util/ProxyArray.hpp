// Copyright 2015 - 2016
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <iterator>

#include "happah/Happah.hpp"

namespace happah {

//DECLARATIONS

template<class Data>
class Deindexer;
template<class Data, class Indices>
class ProxyArray;

//DEFINITIONS

template<class Data>
class Deindexer {
     using Iterator = typename Indices::const_iterator;

public:
     Deindexer(Data& data, Iterator i)
          : m_data(data), m_i(i) {}

     auto operator-(const Deindexer& deindexer) const { return m_i - deindexer.m_i; }

     auto operator+(hpuint n) const {
          auto copy = *this;

          return copy += n;
     }

     auto& operator++() { 
          ++m_i;
          return *this;
     }

     auto& operator--() { 
          --m_i;
          return *this;
     }

     auto& operator+=(hpuint n) {
          m_i += n;
          return *this;
     }

     auto& operator-=(hpuint n) {
          m_i -= n;
          return *this;
     }

     auto operator++(int) const { 
          auto copy = *this;

          return ++copy;
     }

     auto operator--(int) const {
          auto copy = *this;

          return --copy;
     }

     auto& operator[](hpindex n) const { return m_data[m_i[n]]; }

     bool operator!=(const Deindexer& deindexer) const { return deindexer.m_i != m_i; }

     auto& operator*() const { return m_data[*m_i]; }

private:
     Data& m_data;
     Iterator m_i;

};//Deindexer

template<class Data>
class ProxyArray<Data, Indices> {
public:
     ProxyArray(Data& data, const Indices& indices)
          : m_data(data), m_indices(indices) {}

     Deindexer<Data> begin() const { return { m_data, std::begin(m_indices) }; }

     Deindexer<Data> end() const { return { m_data, std::end(m_indices) }; }

     auto& operator[](hpindex n) const { return m_data[m_indices[n]]; }

private:
     Data& m_data;
     const Indices& m_indices;

};//ProxyArray

template<class Data, class Indices>
class ProxyArray {
public:
     ProxyArray(Data& data, const Indices& indices)
          : m_data(data), m_indices(indices) {}

     Deindexer<Data> begin() const { return { m_data, std::begin(m_indices) }; }

     Deindexer<Data> end() const { return { m_data, std::end(m_indices) }; }

     auto& operator[](hpindex n) const { return m_data[m_indices[n]]; }

     Deindexer<Data> operator()(hpindex n) { return { m_data, m_indices(n) }; }

private:
     Data& m_data;
     const Indices& m_indices;

};//ProxyArray

template<class Data, class Indices>
ProxyArray<Data, Indices> deindex(Data& data, const Indices& indices) { return { data, indices }; }

template<class Data, class Indices>
ProxyArray<Data, Indices> deindex(const std::tuple<Data&, Indices&>& array) { return { std::get<0>(array), std::get<1>(array) }; }

}//namespace happah

