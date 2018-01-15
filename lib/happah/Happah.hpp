// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define GLM_FORCE_RADIANS

#include <algorithm>
#include <boost/serialization/strong_typedef.hpp>
#include <boost/variant.hpp>
#include <experimental/tuple>
#include <experimental/filesystem>
#include <glm/glm.hpp>
#include <glm/gtc/vec1.hpp>
#include <glm/gtx/norm.hpp>
#include <iostream>//TODO: remove
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace happah {

//DECLARATIONS

template<class Container>
class back_inserter;

template<class Enumerator, class Transformer>
class EnumeratorTransformer;
class quax;//quaternary index
template<typename T>
class Triples;
class trix;//ternary index
template<typename T>
class Tuples;
template<typename T>
class Quadruples;

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
BOOST_STRONG_TYPEDEF(unsigned int, quat);//quaternary bit
using hpreal = glm::mediump_float;
using hpucolor = glm::uvec4;
using hpuint = unsigned int;
BOOST_STRONG_TYPEDEF(unsigned int, trit);//ternary bit
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

namespace alg {

struct loop {};

}//namespace alg

namespace detail {

template<class T, typename = void>
struct apply;

template<class R, class F>
class LambdaVisitor;

}//namespace detail

template<class Function, class T>
auto apply(Function&& function, T&& t);

Indices::iterator defrag(Indices::iterator begin, Indices::iterator end);

inline Indices::iterator defrag(Indices& indices);

template<class Enumerator>
auto expand(Enumerator e);

template<typename T>
trix find(const Triples<T>& triples, const T& value);

template<typename T>
quax find(const Quadruples<T>& quadruples, const T& value);

template<typename T>
trit find(const Triples<T>& triples, hpindex t, const T& value);

template<typename T>
quat find(const Quadruples<T>& quadruples, hpindex q, const T& value);

inline hpreal length2(const Point3D& point);

template<class Enumerator>
auto make(Enumerator e);

template<typename... T>
std::array<typename std::common_type<T...>::type, sizeof...(T)> make_array(T&&... t);

template<class Container>
back_inserter<Container> make_back_inserter(Container& container);

//Convert a string representation in HPH format.
Indices make_indices(const std::string& indices);

//Import data stored in the given file in HPH format.
Indices make_indices(const std::experimental::filesystem::path& indices);

template<class Value>
auto make_map(hpuint n);

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
hpuint size(const std::vector<T>& ts);//TODO: move to std?

template<class Enumerator>
hpuint size(Enumerator e);

std::string slurp(const std::string& path);

template<class Enumerator, class Transformer>
EnumeratorTransformer<Enumerator, Transformer> transform(Enumerator&& e, Transformer&& transform);

template<class Enumerator, class Visitor>
void visit(Enumerator e, Visitor&& visit);

template<typename T, class Visitor>
void visit(const Triples<T>& triples, Visitor&& visit);

template<typename T, class Visitor>
void visit(const Triples<T>& triples, hpuint n, Visitor&& visit);

template<typename T, class Visitor>
void visit(const Quadruples<T>& quadruples, Visitor&& visit);

template<class R, class Visitor, class... T>
R visit(boost::variant<T...>& variant, Visitor&& visit);

//DEFINITIONS

constexpr hpreal EPSILON = 1e-5;
const quat QUAT0 = quat(0);
const quat QUAT1 = quat(1);
const quat QUAT2 = quat(2);
const quat QUAT3 = quat(3);
const trit TRIT0 = trit(0);
const trit TRIT1 = trit(1);
const trit TRIT2 = trit(2);

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

     auto operator*() const { return happah::apply(m_transform, *m_e); }

     auto& operator++() {
          ++m_e;
          return *this;
     }

     auto& operator+=(hpuint n) {
          while(n--) ++(*this);
          return *this;
     }

     auto operator+(hpuint n) const {
          auto copy = *this;
          return copy += n;
     }

private:
     Enumerator m_e;
     Transformer m_transform;

};//EnumeratorTransformer

class quax {
public:
     quax()
          : m_n(std::numeric_limits<hpuint>::max()) {}

     quax(hpindex q, quat i)
          : m_n((q << 2) | i) {}

     auto getNext() const { return quax(getQuadruple(), quat((++getOffset()) & 3)); }//TODO: rename to operator++?

     quat getOffset() const { return quat(m_n & 3); }

     auto getPrevious() const { return quax(getQuadruple(), quat((--getOffset()) & 3)); }//TODO: rename to operator--?

     hpindex getQuadruple() const { return (m_n & m_mask) >> 2; }

     auto operator==(const quax& x) const { return m_n == x.m_n; }

     operator hpindex() const { return m_n; }

private:
     hpuint m_n;
     static constexpr hpuint m_mask = hpuint(-1) - 3;

};//quax

class trix {
public:
     trix()
          : m_n(std::numeric_limits<hpuint>::max()) {}

     trix(hpindex t, trit i)
          : m_n(t | (i << m_shift)) {}

     auto getNext() const {
          static const trit o[3] = { TRIT1, TRIT2, TRIT0 };

          return trix(getTriple(), o[getOffset()]);
     }

     trit getOffset() const { return trit((m_n & m_mask0) >> m_shift); }

     auto getPrevious() const {
          static const trit o[3] = { TRIT2, TRIT0, TRIT1 };

          return trix(getTriple(), o[getOffset()]);
     }

     hpindex getTriple() const { return m_n & m_mask1; }

     auto operator==(const trix& x) const { return m_n == x.m_n; }

     operator hpindex() const { return 3 * getTriple() + getOffset(); }

private:
     hpuint m_n;
     static constexpr hpuint m_shift = std::numeric_limits<hpuint>::digits - 2;
     static constexpr hpuint m_mask0 = 3 << m_shift;
     static constexpr hpuint m_mask1 = hpuint(-1) - m_mask0;

};//trix

//Triples is a vector whose size is a multiple of three.
template<typename T>
class Triples : public std::vector<T> {
public:
     using std::vector<T>::vector;

     auto getLength() const { return hpuint(3); }

     auto operator()(hpindex t) const { return std::begin(*this) + 3 * t; }

     auto operator()(hpindex t) { return std::begin(*this) + 3 * t; }

     auto& operator()(hpindex t, trit i) const { return (*this)[3 * t + i]; }

     auto& operator()(hpindex t, trit i) { return (*this)[3 * t + i]; }

};//Triples

//Tuples is a vector whose size is a multiple of a given number.
template<typename T>
class Tuples : public std::vector<T> {
public:
     Tuples(hpuint length)
          : std::vector<T>(), m_length(length) {}

     Tuples(hpuint length, hpuint n)
          : std::vector<T>(n), m_length(length) {}

     Tuples(hpuint length, hpuint n, const T& value)
          : std::vector<T>(n, value), m_length(length) {}

     Tuples(hpuint length, std::initializer_list<T> list)
          : std::vector<T>(std::move(list)), m_length(length) {}

     Tuples(const Tuples& other)
          : std::vector<T>(other), m_length(other.m_length) {}

     Tuples(Tuples&& other)
          : std::vector<T>(std::move(other)), m_length(other.m_length) {}

     auto getLength() const { return m_length; }

     auto operator()(hpindex t) const { return std::begin(*this) + m_length * t; }

     auto operator()(hpindex t) { return std::begin(*this) + m_length * t; }

     auto& operator()(hpindex t, hpindex i) const { return (*this)[m_length * t + i]; }

     auto& operator()(hpindex t, hpindex i) { return (*this)[m_length * t + i]; }

private:
     hpuint m_length;

};//Tuples

//Quadruples is a vector whose size is a multiple of four.
template<typename T>
class Quadruples : public std::vector<T> {
public:
     using std::vector<T>::vector;

     auto getLength() const { return hpuint(4); }

     auto operator()(hpindex q) const { return std::begin(*this) + (q << 2); }

     auto operator()(hpindex q) { return std::begin(*this) + (q << 2); }

     auto& operator()(hpindex q, quat i) const { return (*this)[(q << 2) + i]; }

     auto& operator()(hpindex q, quat i) { return (*this)[(q << 2) + i]; }

};//Quadruples

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

template<class R, class F>
class LambdaVisitor : boost::static_visitor<R> {
public:
     LambdaVisitor(F f) : m_f(std::move(f)) {}

     template<class... T>
     R operator()(T&... ts) const { return m_f(ts...); }

private:
     F m_f;

};//LambdaVisitor

}//namespace detail

template<class Function, class T>
auto apply(Function&& function, T&& t) { return detail::apply<T>::call(std::forward<Function>(function), std::forward<T>(t)); }

inline Indices::iterator defrag(Indices& indices) { return defrag(std::begin(indices), std::end(indices)); }

template<class Enumerator>
auto expand(Enumerator e) {
     using T = typename std::remove_const<typename std::remove_reference<typename remove_tuple<decltype(*e)>::type>::type>::type;

     auto ts = std::vector<T>();
     auto push_back = make_back_inserter(ts);
     do happah::apply(push_back, *e); while(++e);
     return ts;
}

template<typename T>
trix find(const Triples<T>& triples, const T& value) {
     auto n = std::distance(std::begin(triples), std::find(std::begin(triples), std::end(triples), value));
     auto t = hpindex(n / 3);
     auto i = trit(n - 3 * t);

     return { t, i };
}

template<typename T>
quax find(const Quadruples<T>& quadruples, const T& value) {
     auto n = std::distance(std::begin(quadruples), std::find(std::begin(quadruples), std::end(quadruples), value));
     auto q = hpindex(n >> 2);
     auto i = quat(n - (q << 2));

     return { q, i };
}

template<typename T>
trit find(const Triples<T>& triples, hpindex t, const T& value) {
     auto i = std::begin(triples) + 3 * t;

     if(value == i[0]) return TRIT0;
     if(value == i[1]) return TRIT1;
     assert(value == i[2]);
     return TRIT2;
}

template<typename T>
quat find(const Quadruples<T>& quadruples, hpindex q, const T& value) {
     auto i = std::begin(quadruples) + (q << 2);

     if(value == i[0]) return QUAT0;
     if(value == i[1]) return QUAT1;
     if(value == i[2]) return QUAT2;
     assert(value == i[3]);
     return QUAT3;
}

inline hpreal length2(const Point3D& point) { return glm::length2(point); }

template<class Enumerator>
auto make(Enumerator e) {
     using T = typename std::remove_const<typename std::remove_reference<decltype(*e)>::type>::type;

     auto ts = std::vector<T>();
     do ts.push_back(*e); while(++e);
     return ts;
}

template<typename... T>
std::array<typename std::common_type<T...>::type, sizeof...(T)> make_array(T&&... t) { return { std::forward<T>(t)... }; }

template<class Container>
back_inserter<Container> make_back_inserter(Container& container) { return back_inserter<Container>(container); }

template<class Value>
auto make_map(hpuint n) {
     using Key = std::pair<hpuint, hpuint>;

     auto getHash = [](const Key& k) -> uint64_t {
          int32_t d = k.first - k.second;
          int32_t min = k.second + (d & d >> 31);
          int32_t max = k.first - (d & d >> 31);
          return ((uint64_t)max << 32 | min);
     };

     auto isKeysEqual = [](const Key& k1, const Key& k2) { return (k1.first == k2.first && k1.second == k2.second) || (k1.first == k2.second && k1.second == k2.first); };

     using Map = std::unordered_map<Key, Value, decltype(getHash), decltype(isKeysEqual)>;

     return Map(n, getHash, isKeysEqual);
}

inline Point3D mix(const Point3D& point, hpreal lambda) { return point * lambda; }

inline Point2D mix(const Point2D& p0, hpreal u, const Point2D& p1, hpreal v, const Point2D& p2, hpreal w) { return p0 * u + p1 * v + p2 * w; }

template<typename F>
void repeat(unsigned n, F f) { while(n--) f(); }

template<class T>
hpuint size(const std::vector<T>& ts) { return ts.size(); }

template<class Enumerator, class Transformer>
EnumeratorTransformer<Enumerator, Transformer> transform(Enumerator&& e, Transformer&& transform) { return { std::forward<Enumerator>(e), std::forward<Transformer>(transform) }; }

template<class Enumerator, class Visitor>
void visit(Enumerator e, Visitor&& visit) { do ::happah::apply(visit, *e); while(++e); }

template<typename T, class Visitor>
void visit(const Triples<T>& triples, Visitor&& visit) { happah::visit(triples, hpuint(triples.size() / 3), std::forward<Visitor>(visit)); }

template<typename T, class Visitor>
void visit(const Triples<T>& triples, hpuint n, Visitor&& visit) { for(auto i = std::begin(triples), end = std::begin(triples) + 3 * n; i != end; i += 3) visit(i[0], i[1], i[2]); }

template<typename T, class Visitor>
void visit(const Tuples<T>& tuples, Visitor&& visit) {
     auto n = tuples.getLength();

     for(auto i = std::begin(tuples), end = std::end(tuples); i != end; i += n) visit(i);
}

template<typename T, class Visitor>
void visit(const Quadruples<T>& quadruples, Visitor&& visit) { for(auto i = std::begin(quadruples), end = std::end(quadruples); i != end; i += 4) visit(i[0], i[1], i[2], i[3]); }

template<class R, class Visitor, class... T>
R visit(boost::variant<T...>& variant, Visitor&& visit) { return boost::apply_visitor(detail::LambdaVisitor<R, Visitor>(std::forward<Visitor>(visit)), variant); }

namespace Color {
     static const hpcolor BLUE(0.0,0.0,1.0,1.0);
     static const hpcolor COPPER(0.7038,0.27048,0.0828,1.0);
     static const hpcolor GREEN(0.0,1.0,0.0,1.0);
     static const hpcolor MAGENTA(1.0,0.0,1.0,1.0);
     static const hpcolor RED(1.0,0.0,0.0,1.0);
     static const hpcolor WHITE(1.0);
}//namespace Color

}//namespace happah

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

