// Copyright 2015 - 2017
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_tuple.hpp>
#include <boost/spirit/home/x3.hpp>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "happah/Happah.h"

namespace happah {

namespace format {

namespace obj {

struct Content {
     struct  {
          std::vector<hpreal> coordinates;
          hpuint dimension;
     } vs;
     struct {
          std::vector<hpuint> indices;
          std::vector<hpuint> offsets;
          bool normal;//whether fs contain indices to normals
          bool texture;//whether fs contain indices to textures

     } fs;

};

template<class Iterator>
Content read(Iterator begin, Iterator end);

Content read(const std::string& path);

namespace detail {

namespace x3 = boost::spirit::x3;

const x3::rule<class _v, std::tuple<hpreal, hpreal, hpreal, hpreal> > v = "v";
const auto v_def = x3::double_ >> x3::double_ >> x3::double_ >> (x3::double_ | x3::attr(std::numeric_limits<hpreal>::max()));
BOOST_SPIRIT_DEFINE(v);

}//namespace detail

template<class Iterator>
Content read(Iterator begin, Iterator end) {
     namespace x3 = boost::spirit::x3;

     auto content = Content();

     if(!phrase_parse(begin, end, +detail::v, x3::ascii::blank)) throw std::runtime_error("Failed to parse content.");

     return content;
}

/*
template<>
class Reader<OBJ> {
public:
     template<typename Iterator>
     static std::tuple<Positions, boost::optional<Normals>, Faces> load(Iterator first, Iterator last) {// load obj
          namespace qi = boost::spirit::qi;
          namespace ascii = boost::spirit::ascii;
          namespace phoenix = boost::phoenix;

          using qi::_1;
          using qi::_a;
          using qi::_b;
          using qi::_c;
          using qi::_val;
          using qi::eol;
          using qi::char_;
          using qi::float_;
          using qi::uint_;
          using qi::eps;
          using boost::spirit::lexeme;

          Positions positions;
          boost::optional<Normals> normals = Normals();
          Faces faces;
          Face face;

          auto destroy_normals = [&]() { normals.reset(); };
          auto push_face = [&]() { faces.push_back(std::move(face)); };
          auto push_normal = [&](glm::vec3 normal) { normals->push_back(normal); };

          qi::rule<Iterator, glm::vec2(), ascii::blank_type> vec2_ = float_ >> float_;

          qi::rule<Iterator, glm::vec3(), ascii::blank_type> vec3_ = float_ >> float_ >> float_;

          qi::rule<Iterator, Corner(), qi::locals<unsigned int, boost::optional<unsigned int>, boost::optional<unsigned int> >, ascii::blank_type> corner_ = (uint_[_a = _1] >> -('/' >> ((-uint_[_b = _1] >> '/' >> uint_[_c = _1]) | uint_[_b = _1])))[_val = phoenix::construct<Corner>(_a, _b, _c)];

          qi::rule<Iterator> comment_line = '#' >> *(char_ - eol);
          qi::rule<Iterator, ascii::blank_type> f_line = 'f' >> +(corner_[phoenix::push_back(phoenix::ref(face), _1)]);
          qi::rule<Iterator, ascii::blank_type> g_line = 'g' >> *(char_ - eol);
          qi::rule<Iterator, ascii::blank_type> mtllib_line = "mtllib" >> +(char_ - eol);
          qi::rule<Iterator, ascii::blank_type> o_line = 'o' >> *(char_ - eol);
          qi::rule<Iterator, ascii::blank_type> s_line = 's' >> (uint_|lexeme["off"]);
          qi::rule<Iterator, ascii::blank_type> usemtl_line = "usemtl" >> +(char_ - eol);
          qi::rule<Iterator, ascii::blank_type> v_line = 'v' >> vec3_[phoenix::push_back(phoenix::ref(positions), _1)];
          qi::rule<Iterator, ascii::blank_type> vn_line = "vn" >> vec3_[push_normal];
          qi::rule<Iterator, ascii::blank_type> vt_line = "vt" >> vec2_ >> -float_;

          auto grammar = (
               -((comment_line | mtllib_line | o_line) % eol >> +eol)
               >> v_line % eol
               >> -(+eol >> comment_line % eol)
               >> -(+eol >> vt_line % eol)
               >> -(+eol >> comment_line % eol)
               >> ((+eol >> vn_line % eol) | eps[destroy_normals])
               >> -(+eol >> comment_line % eol)
               >> -(+eol >> (f_line[push_face] | s_line | g_line | usemtl_line) % eol)
               >> -(+eol >> (comment_line | s_line) % eol)
          );

          if(!qi::phrase_parse(first, last, grammar, ascii::blank)) throw std::domain_error("Content is not a valid OBJ.");

          return std::make_tuple(positions, normals, faces);
     }

};

*/

}//namespace obj

}//namespace format

}//namespace happah

