/**********************************************************************************
 * Copyright (c) 2012-2015  See the COPYRIGHT-HOLDERS file in the top-level
 * directory of this distribution and at http://github.com/happah-graphics/happah.
 *
 * This file is part of Happah. It is subject to the license terms in the LICENSE
 * file found in the top-level directory of this distribution and at
 * http://github.com/happah-graphics/happah. No part of Happah, including this
 * file, may be copied, modified, propagated, or distributed except according to
 * the terms contained in the LICENSE file.
 **********************************************************************************/

#include "happah/readers/ReaderLOP.h"

#include <boost/iostreams/device/mapped_file.hpp>

auto ReaderLOP::read(const char* path) -> IndicesArrays {
     using boost::iostreams::mapped_file;
     mapped_file file(path, mapped_file::readonly);
     auto begin = file.const_data();
     auto end = begin + file.size();
     return read(begin, end);
}

