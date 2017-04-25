// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_tuple.hpp>
#include <boost/spirit/home/x3.hpp>
#include <stdexcept>
#include <tuple>

#include "happah/Happah.h"
#include "happah/formats/off.h"

namespace happah {

namespace off {

namespace x3 = boost::spirit::x3;

namespace detail {

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

auto make_vertex_size(const Header& header) {
     auto n = header.dimension;
     if(header.normal) n <<= 1;
     if(header.texture) n += 2;
     if(header.color) n += 4;
     return n;
}

}//namespace detail

//NOTE: All colors must be represented using floating point numbers.
//NOTE: A vertex color must be in RGBA format.  https://sourceforge.net/p/geomview/mailman/message/2004718/
//NOTE: Officially, integer colors are only allowed for faces and is a backwards-compatibility hack.  https://sourceforge.net/p/geomview/mailman/message/2004728/
//NOTE: If at least one face has a color, then all faces have a color.  A face color does not necessarily
template<class Iterator>
Content import(Iterator begin, Iterator end) {
     auto header = Header();
     auto vertices = Vertices();
     auto faces = Faces();

     if(!phrase_parse(begin, end, detail::header, x3::ascii::space, header)) throw std::runtime_error("Failed to parse header.");

     if(header.nVertices > 0) {
          auto n = detail::make_vertex_size(header);
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

          if(!phrase_parse(begin, end, x3::repeat(header.nFaces)[(detail::face_ >> x3::omit[-detail::color[push_back]] >> (+x3::eol | x3::eoi))[increment]], x3::ascii::blank, faces.vertices)) throw std::runtime_error("Failed to parse faces.");
     }

     return std::make_tuple(std::move(header), std::move(vertices), std::move(faces));
}

}//namespace off

}//namespace happah

BOOST_FUSION_ADAPT_STRUCT(
     happah::off::Header,
     (bool, texture)
     (bool, color)
     (bool, normal)
     (hpuint, dimension)
     (hpuint, nVertices)
     (hpuint, nFaces)
     (hpuint, nEdges)
)

