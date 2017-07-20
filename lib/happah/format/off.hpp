// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_tuple.hpp>
#include <boost/spirit/home/x3.hpp>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "happah/Happah.hpp"
#include "happah/format/default.hpp"
#include "happah/geometries/TriangleMesh.hpp"
#include "happah/geometries/TriangleGraph.hpp"
#include "happah/geometries/Vertex.hpp"

namespace happah {

namespace format {

namespace off {

struct Header {
     bool color;
     hpuint dimension;
     hpuint nEdges;
     hpuint nFaces;
     bool normal;
     hpuint nVertices;
     bool texture;

};//Header

struct Content {
     Header header;
     std::vector<hpreal> vertices;
     struct {
          Indices vertices;
          std::vector<hpreal> colors;
          Indices indices;
     } faces;

};//Content

template<class Vertex>
Header make_header(hpuint nFaces, hpuint nVertices);

hpuint make_vertex_size(const Header& header);

template<class Iterator>
Content read(Iterator begin, Iterator end);

Content read(const std::string& path);

template<class Vertex>
void write(const TriangleMesh<Vertex>& mesh, const std::string& path);

namespace detail {

namespace x3 = boost::spirit::x3;

auto increment = [](auto& context) { ++x3::_attr(context); };

const x3::rule<class _color, std::tuple<hpreal, hpreal, hpreal, hpreal> > color = "color";
const auto color_def = x3::double_ >> x3::double_ >> x3::double_ >> (x3::double_ | x3::attr(std::numeric_limits<hpreal>::max()));
BOOST_SPIRIT_DEFINE(color);

const x3::rule<class _header, Header, true> header = "header";
const auto header_def = x3::matches[x3::lit("ST")] >> x3::matches[x3::lit('C')] >> x3::matches[x3::lit('N')] >> ((x3::lit("4nOFF") >> x3::int_[increment]) | (x3::lit("4OFF") >> x3::attr(4)) | (x3::lit("nOFF") >> x3::int_) | (x3::lit("OFF") >> x3::attr(3)) | x3::attr(3)) >> x3::int_ >> x3::int_ >> x3::int_;
BOOST_SPIRIT_DEFINE(header);

struct face_parser : x3::parser<face_parser> {

     using attribute_type = Indices;
     static const bool has_attribute = true;

     template<typename Iterator, typename Context, typename Attribute>
     bool parse(Iterator& first, const Iterator& last, const Context& context, x3::unused_type, Attribute& attr) const {
          auto n = 0u;

          if(!x3::uint_.parse(first, last, context, x3::unused, n)) throw std::runtime_error("Failed to parse face.");
          x3::traits::push_back(attr, n);
          while(n--) {
               auto temp = 0u;
               if(!x3::uint_.parse(first, last, context, x3::unused, temp)) throw std::runtime_error("Failed to parse face.");
               x3::traits::push_back(attr, temp);
          }

          return true;
     }

};
const face_parser face_ = {};

}//namespace detail

using happah::format::operator<<;

template<class Stream>
Stream& operator<<(Stream& stream, const Header& header) {
     if(header.color) stream << 'C';
     if(header.normal) stream << 'N';
     if(header.dimension == 3) stream << "OFF\n";
     else if(header.dimension == 4) stream << "4OFF\n";
     else {
          stream << "nOFF\n";
          stream << header.dimension << '\n';
     }
     stream << header.nVertices << ' ' << header.nFaces << " 0";
     return stream;
}

template<class Stream, class T>
Stream& operator<<(Stream& stream, const std::vector<T>& ts) {
     for(auto t = std::begin(ts), end = std::end(ts) - 1; t != end; ++t) {
          stream << *t;
          stream << '\n';
     }
     stream << ts.back();
     return stream;
}

template<class Vertex>
Header make_header(hpuint nFaces, hpuint nVertices) {
     auto header = Header();
     header.color = contains_color<Vertex>::value;
     header.dimension = Vertex::SPACE::DIMENSION;
     header.nFaces = nFaces;
     header.normal = contains_normal<Vertex>::value;
     header.nVertices = nVertices;
     return header;
}

//TODO: more precise parsing error messages
//TODO: compare this parsing with a hand-written recursive descent parser
//TODO: change struct Faces to separate vertices into vertices and offsets; the most likely case is that all polygons have the same valence (3 for triangle meshes, 4 for quadrilateral meshes); then we only have to move the vertices into TriangleMesh and do not need to post-process the offsets
//TODO: normalize(Content): add missing alpha channels to colors, convert integer colors to floating point representation
//TODO: analyze(Iterator,Iterator) count how many face indices, face colors, detect ambiguous case; call this inside read but before parsing into content
template<class Iterator>
Content read(Iterator begin, Iterator end) {
     namespace x3 = boost::spirit::x3;

     auto content = Content();
     auto& header = content.header;
     auto& vertices = content.vertices;
     auto& faces = content.faces;

     if(!phrase_parse(begin, end, detail::header, x3::ascii::space, header)) throw std::runtime_error("Failed to parse header.");

     if(header.nVertices > 0) {
          auto n = make_vertex_size(header);
          vertices.reserve(n * header.nVertices);
          if(header.color) {
               if(!phrase_parse(begin, end, x3::repeat(header.nVertices)[x3::repeat(n - 1)[x3::double_] >> (x3::double_ | x3::attr(std::numeric_limits<hpreal>::max())) >> (+x3::eol | x3::eoi)], x3::ascii::blank, vertices)) throw std::runtime_error("Failed to parse vertices.");
          } else if(!phrase_parse(begin, end, x3::repeat(n * header.nVertices)[x3::double_], x3::ascii::space, vertices)) throw std::runtime_error("Failed to parse vertices.");
     }

     if(header.nFaces > 0) {
          auto face = 0u;

          auto increment = [&](auto& context) { ++face; };
          auto push_back = [&](auto& context) {
               auto& color = x3::_attr(context);
               faces.colors.push_back(std::get<0>(color));
               faces.colors.push_back(std::get<1>(color));
               faces.colors.push_back(std::get<2>(color));
               faces.colors.push_back(std::get<3>(color));
               faces.indices.push_back(face);
          };

          faces.vertices.reserve(header.nFaces << 2);
          if(!phrase_parse(begin, end, x3::repeat(header.nFaces)[(detail::face_ >> x3::omit[-detail::color[push_back]] >> (+x3::eol | x3::eoi))[increment]], x3::ascii::blank, faces.vertices)) throw std::runtime_error("Failed to parse faces.");
     }

     return content;
}

template<class Vertex>
void write(const TriangleGraph<Vertex>& graph, const std::string& path) {
     const auto& vertices = graph.getVertices();
     auto stream = std::ofstream();

     stream.exceptions(std::ofstream::failbit);
     stream.open(path);
     stream << make_header<Vertex>(size(graph), size(vertices)) << "\n\n";
     stream << std::fixed;
     stream << vertices << "\n\n";
     visit_triplets(std::begin(graph.getEdges()), size(graph), 3, [&](const auto& edge0, const auto& edge1, const auto& edge2) { stream << "3 " << edge2.vertex << ' ' << edge0.vertex << ' ' << edge1.vertex << '\n'; });
}

template<class Vertex>
void write(const TriangleMesh<Vertex>& mesh, const std::string& path) {
     const auto& vertices = mesh.getVertices();
     const auto& indices = mesh.getIndices();
     auto stream = std::ofstream();

     stream.exceptions(std::ofstream::failbit);
     stream.open(path);
     stream << make_header<Vertex>(size(mesh), size(vertices)) << "\n\n";
     stream << std::fixed;
     stream << vertices << "\n\n";
     visit_triplets(indices, [&](auto i0, auto i1, auto i2) { stream << "3 " << i0 << ' ' << i1 << ' ' << i2 << '\n'; });
}

}//namespace off

}//namespace format

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(const format::off::Content& content) {
     auto vertices = std::vector<Vertex>();
     auto indices = Indices();
     auto& header = content.header;
     auto n = format::off::make_vertex_size(header);
     vertices.reserve(header.nVertices);
     indices.reserve(3 * header.nFaces);

     for(auto i = std::begin(content.vertices), end = std::end(content.vertices); i != end; i += n) vertices.push_back(make_vertex<Vertex>(i));
     for(auto i = std::begin(content.faces.vertices), end = std::end(content.faces.vertices); i != end; ++i) {
          indices.push_back(*(++i));
          indices.push_back(*(++i));
          indices.push_back(*(++i));
     }

     return make_triangle_mesh(std::move(vertices), std::move(indices));
}

}//namespace happah

BOOST_FUSION_ADAPT_STRUCT(
     happah::format::off::Header,
     (bool, texture)
     (bool, color)
     (bool, normal)
     (hpuint, dimension)
     (hpuint, nVertices)
     (hpuint, nFaces)
     (hpuint, nEdges)
)

