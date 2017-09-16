// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define GLM_FORCE_RADIANS

#include <experimental/tuple>
#include <experimental/filesystem>
#include <glm/glm.hpp>
#include <glm/gtc/vec1.hpp>
#include <glm/gtx/norm.hpp>
#include <iostream>//TODO: remove
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace happah {

//DECLARATIONS

template<class Container>
class back_inserter;

template<class Enumerator, class Transformer>
class EnumeratorTransformer;

using hpcolor = glm::vec4;
using hpindex = unsigned int;
using hpint = int;

struct hpir;

struct hpijr;

struct hpijkr;

struct hpijklr;

using hpmat2x2 = glm::mat2x2;
using hpmat3x3 = glm::mat3x3;
using hpmat4x4 = glm::mat4x4;
using hprandom = std::mt19937;
using hpreal = glm::mediump_float;
using hpucolor = glm::uvec4;
using hpuint = unsigned int;
using hpvec1 = glm::vec1;
using hpvec2 = glm::vec2;
using hpvec3 = glm::vec3;
using hpvec4 = glm::vec4;

using Indices = std::vector<hpindex>;

template<typename>
struct is_tuple;

using Point1D = hpvec1;
using Point2D = hpvec2;
using Point3D = hpvec3;
using Point4D = hpvec4;

template<class T, typename = void>
struct remove_tuple;

using Vector1D = hpvec1;
using Vector2D = hpvec2;
using Vector3D = hpvec3;
using Vector4D = hpvec4;

namespace detail {

template<class T, typename = void>
struct apply;

}//namespace detail

template<class Function, class T>
auto apply(Function&& function, T&& t);

Indices::iterator defrag(Indices::iterator begin, Indices::iterator end);

inline Indices::iterator defrag(Indices& indices);

template<class Enumerator>
auto expand(Enumerator e);

inline hpreal length2(const Point3D& point);

template<class Enumerator>
auto make(Enumerator e);

template<class Container>
back_inserter<Container> make_back_inserter(Container& container);

//Convert a string representation in HPH format.
Indices make_indices(const std::string& indices);

//Import data stored in the given file in HPH format.
Indices make_indices(const std::experimental::filesystem::path& indices);

//Convert a string representation in HPH format.
std::vector<hpreal> make_reals(const std::string& reals);

//Import data stored in the given file in HPH format.
std::vector<hpreal> make_reals(const std::experimental::filesystem::path& reals);

inline Point3D mix(const Point3D& point, hpreal lambda);

inline Point2D mix(const Point2D& p0, hpreal u, const Point2D& p1, hpreal v, const Point2D& p2, hpreal w);

inline std::experimental::filesystem::path p(const std::string& path) { return { path }; }

template<typename F>
void repeat(unsigned n, F f);

template<class T>
hpuint size(const std::vector<T>& ts);

std::string slurp(const std::string& path);

template<class Enumerator, class Visitor>
void visit(Enumerator e, Visitor&& visit);

//DEFINITIONS

extern hprandom happah_random_engine;

constexpr hpreal EPSILON = 1e-5;
constexpr hpuint UNULL = std::numeric_limits<hpuint>::max();

template<class Container>
class back_inserter {
public:
     back_inserter(Container& container)
          : m_container(container) {}

     template<typename T>
     void operator()(const T& value) const { m_container.push_back(value); }

     template<typename T, typename... Ts>
     void operator()(const T& value, const Ts&... values) const {
          m_container.push_back(value);
          this->operator()(values...);
     }

private:
     Container& m_container;

};//back_inserter

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

struct hpir {
     hpuint i;
     hpreal r;

     hpir(hpuint i, hpreal r)
          : i(i), r(r) {}

};//hpir

struct hpijr {
     hpuint i;
     hpuint j;
     hpreal r;

     hpijr(hpuint i, hpuint j, hpreal r)
          : i(i), j(j), r(r) {}

};//hpijr

struct hpijkr {
     hpuint i;
     hpuint j;
     hpuint k;
     hpreal r;

     hpijkr(hpuint i, hpuint j, hpuint k, hpreal r)
          : i(i), j(j), k(k), r(r) {}

};//hpijkr

struct hpijklr {
     hpuint i;
     hpuint j;
     hpuint k;
     hpuint l;
     hpreal r;

     hpijklr(hpuint i, hpuint j, hpuint k, hpuint l, hpreal r)
          : i(i), j(j), k(k), l(l), r(r) {}

};//hpijklr

template<typename>
struct is_tuple : std::false_type {};

template<typename... T>
struct is_tuple<std::tuple<T...> > : std::true_type {};

template<class T, typename>
struct remove_tuple { using type = T; };

template<class T>
struct remove_tuple<T, typename std::enable_if<is_tuple<T>::value>::type> { using type = typename std::tuple_element<0, T>::type; };

namespace detail {

template<class T, typename>
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

inline Indices::iterator defrag(Indices& indices) { return defrag(std::begin(indices), std::end(indices)); }

template<class Enumerator>
auto expand(Enumerator e) {
     using T = typename std::remove_const<typename std::remove_reference<typename remove_tuple<decltype(*e)>::type>::type>::type;

     auto ts = std::vector<T>();
     auto push_back = make_back_inserter(ts);
     do apply(push_back, *e); while(++e);
     return ts;
}

inline hpreal length2(const Point3D& point) { return glm::length2(point); }

template<class Enumerator>
auto make(Enumerator e) {
     using T = typename std::remove_const<typename std::remove_reference<decltype(*e)>::type>::type;

     auto ts = std::vector<T>();
     do ts.push_back(*e); while(++e);
     return ts;
}

template<class Container>
back_inserter<Container> make_back_inserter(Container& container) { return back_inserter<Container>(container); }

inline Point3D mix(const Point3D& point, hpreal lambda) { return point * lambda; }

inline Point2D mix(const Point2D& p0, hpreal u, const Point2D& p1, hpreal v, const Point2D& p2, hpreal w) { return p0 * u + p1 * v + p2 * w; }

template<typename F>
void repeat(unsigned n, F f) { while(n--) f(); }

template<class T>
hpuint size(const std::vector<T>& ts) { return ts.size(); }

template<class Enumerator, class Visitor>
void visit(Enumerator e, Visitor&& visit) { do apply(visit, *e); while(++e); }

namespace Color {
     static const hpcolor BLUE(0.0,0.0,1.0,1.0);
     static const hpcolor COPPER(0.7038,0.27048,0.0828,1.0);
     static const hpcolor GREEN(0.0,1.0,0.0,1.0);
     static const hpcolor MAGENTA(1.0,0.0,1.0,1.0);
     static const hpcolor RED(1.0,0.0,0.0,1.0);
     static const hpcolor WHITE(1.0);
}

}

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
//TODO: figure out const and non-const iterator design

