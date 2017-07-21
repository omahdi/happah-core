// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define GLM_FORCE_RADIANS

#include <iostream>//TODO: remove
#include <string>
#include <tuple>
#include <experimental/tuple>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/vec1.hpp>

template<class... T>
auto make_vector(T&&... args) {
    using U = typename std::tuple_element<0, std::tuple<T...>>::type;
    return std::vector<U>{std::forward<T>(args)...};
}

typedef glm::vec4 hpcolor;
typedef int hpint;
typedef glm::mat2x2 hpmat2x2;
typedef glm::mat3x3 hpmat3x3;
typedef glm::mat4x4 hpmat4x4;
typedef glm::mediump_float hpreal;
typedef glm::uvec4 hpucolor;
typedef unsigned int hpuint;
using hpindex = unsigned int;
using hpsize = std::size_t;//TODO: index and size have to be comparable: index < size_t; make everything int and -Wsign-compare to avoid compiler warnings when using std containers
typedef glm::ivec3 hpint3;
typedef glm::uvec3 hpuint3;
typedef glm::vec1 hpvec1;
typedef glm::vec2 hpvec2;
typedef glm::vec3 hpvec3;
typedef glm::vec4 hpvec4;

enum Diagonals {//TODO: QuadDiagonals?
     A = 0x1,
     AB = 0x4,
     B = 0x2,
     NONE = 0x0
};

namespace happah {

namespace Color {
     static const hpcolor BLUE(0.0,0.0,1.0,1.0);
     static const hpcolor COPPER(0.7038,0.27048,0.0828,1.0);
     static const hpcolor GREEN(0.0,1.0,0.0,1.0);
     static const hpcolor MAGENTA(1.0,0.0,1.0,1.0);
     static const hpcolor RED(1.0,0.0,0.0,1.0);
     static const hpcolor WHITE(1.0);
}

template<typename>
struct is_tuple : std::false_type {};

template<typename... T>
struct is_tuple<std::tuple<T...> > : std::true_type {};

namespace detail {

template<class T, typename = void>
struct apply {
     template<class Function>
     static auto call(Function&& function, T&& t) { return function(std::forward<T>(t)); }
};

template<class T>
struct apply<T, typename std::enable_if<is_tuple<T>::value>::type> {
     template<class Function>
     static auto call(Function&& function, T&& t) { return std::experimental::fundamentals_v1::apply(std::forward<Function>(function), std::forward<T>(t)); }
};

}//namespace detail

template<class Function, class T>
auto apply(Function&& function, T&& t) { return detail::apply<T>::call(std::forward<Function>(function), std::forward<T>(t)); }

template<class Enumerator, class Transformer>
class EnumeratorTransformer {
public:
     EnumeratorTransformer(Enumerator&& e, Transformer&& transform)
          : m_e(e), m_transform(transform) {}

     explicit operator bool() const { return bool(m_e); }

     auto operator*() const { return apply(m_transform, *m_e); }

     auto& operator++() {
          ++m_e;
          return *this;
     }

private:
     Enumerator m_e;
     Transformer m_transform;

};//EnumeratorTransformer

constexpr hpreal EPSILON = 1e-5;
constexpr hpuint UNULL = std::numeric_limits<hpuint>::max();

using Indices = std::vector<hpindex>;

Indices::iterator defrag(Indices::iterator begin, Indices::iterator end);

Indices::iterator defrag(Indices& indices);

Indices make_indices(const std::string& path);

std::vector<hpreal> make_reals(const std::string& path);

template<class T>
hpuint size(const std::vector<T>& ts) { return ts.size(); }

std::string slurp(const std::string& path);

struct hpir {
     hpuint i;
     hpreal r;

     hpir(hpuint i, hpreal r)
          : i(i), r(r) {}

};
struct hpijr {
     hpuint i;
     hpuint j;
     hpreal r;

     hpijr(hpuint i, hpuint j, hpreal r)
          : i(i), j(j), r(r) {}

};

struct hpijkr {
     hpuint i;
     hpuint j;
     hpuint k;
     hpreal r;

     hpijkr(hpuint i, hpuint j, hpuint k, hpreal r)
          : i(i), j(j), k(k), r(r) {}

};

struct hpijklr {
     hpuint i;
     hpuint j;
     hpuint k;
     hpuint l;
     hpreal r;

     hpijklr(hpuint i, hpuint j, hpuint k, hpuint l, hpreal r)
          : i(i), j(j), k(k), l(l), r(r) {}

};

struct MissingImplementationException : std::exception {};

template<class Container>
class back_inserter {
public:
     back_inserter(Container& container)
          : m_container(container) {}

     template<typename T>
     void operator()(const T& value) { m_container.push_back(value); }

     template<typename T, typename... Ts>
     void operator()(const T& value, const Ts&... values) {
          m_container.push_back(value);
          this->operator()(values...);
     }

private:
     Container& m_container;

};

template<class Container>
back_inserter<Container> make_back_inserter(Container& container) { return back_inserter<Container>(container); }

template<typename F>
void repeat(unsigned n, F f) { while(n--) f(); }

template<class Enumerator, class Visitor>
void visit(Enumerator e, Visitor&& visit) { while(e) { apply(visit, *e); ++e; } }

template<class Enumerator, class Transformer, class Visitor>
void visit(EnumeratorTransformer<Enumerator, Transformer> e, Visitor&& visit) { while(e) { apply(visit, *e); ++e; } }

}

#define BUILD_TUPLE_HANDLER_METHODS(NAME, HANDLER) \
template<unsigned long I0, unsigned long... Is>\
void NAME(std::index_sequence<I0, Is...>) { HANDLER(std::get<I0>(m_i)); NAME(std::index_sequence<Is...>()); }\
template<unsigned long I0>\
void NAME(std::index_sequence<I0>) { HANDLER(std::get<I0>(m_i)); }

//TODO: automatic code styler
//TODO: deb/ppa package for happah install, -dev, -dev-eclipse, -dev-vim, plus documentation about use
//TODO: doxygen docs generation
//TODO: document and organize eclipse use
//TODO: delete default copy constructor and assignment operators
//TODO: example segment mesh, triangle mesh, point cloud code snippets
//TODO: marching cubes or similar: www.imm.dtu.dk/~janba/gallery/polygonization.html
//TODO: isosurface blending: http://www.povray.org/documentation/view/3.6.1/73/
//TODO: TriangleStrips
//TODO: update sphere, plane, cylinder, etc. with toTriangleStrips methods for memory efficient representations
//TODO: ArraysIterator cbegin
//TODO: figure out const and non-const iterator design

