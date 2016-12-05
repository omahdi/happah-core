// Copyright 2015
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define GLM_FORCE_RADIANS

#include <iostream>//TODO: remove
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/vec1.hpp>

typedef glm::vec4 hpcolor;
typedef int hpint;
typedef glm::mat2x2 hpmat2x2;
typedef glm::mat3x3 hpmat3x3;
typedef glm::mat4x4 hpmat4x4;
typedef glm::mediump_float hpreal;
typedef glm::uvec4 hpucolor;
typedef unsigned int hpuint;
using hpindex = int;
using hpsize = std::size_t;//TODO: index and size have to be comparable: index < size_t; make everything int and -Wsign-compare to avoid compiler warnings when using std containers
typedef glm::ivec3 hpint3;
typedef glm::uvec3 hpuint3;
typedef glm::vec1 hpvec1;
typedef glm::vec2 hpvec2;
typedef glm::vec3 hpvec3;
typedef glm::vec4 hpvec4;

namespace Color {
     static const hpcolor BLUE(0.0,0.0,1.0,1.0);
     static const hpcolor COPPER(0.7038,0.27048,0.0828,1.0);
     static const hpcolor GREEN(0.0,1.0,0.0,1.0);
     static const hpcolor MAGENTA(1.0,0.0,1.0,1.0);
     static const hpcolor RED(1.0,0.0,0.0,1.0);
     static const hpcolor WHITE(1.0);
}

enum Diagonals {//TODO: QuadDiagonals?
     A = 0x1,
     AB = 0x4,
     B = 0x2,
     NONE = 0x0
};

namespace happah {

constexpr hpreal EPSILON = 1e-5;
constexpr hpuint UNULL = -1;

using Indices = std::vector<hpuint>;

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

template<class Stream>
Stream& operator<<(Stream& out, const hpvec2& v) { 
     out << v.x << ' ' << v.y;
     return out; 
}

template<class Stream>
Stream& operator<<(Stream& out, const hpvec3& v) { 
     out << v.x << ' ' << v.y << ' ' << v.z;
     return out; 
}

template<class Stream>
Stream& operator<<(Stream& out, const hpvec4& v) { 
     out << v.x << ' ' << v.y << ' ' << v.z << ' ' << v.w;
     return out; 
}

struct MissingImplementationException : std::exception {};

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

