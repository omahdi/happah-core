# Copyright 2015 - 2017
#   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

AC_INIT(libhappah, 0.0, pherman@ira.uka.de)
AC_CONFIG_AUX_DIR(build)
AC_CONFIG_MACRO_DIR(m4)
AC_ENABLE_SHARED
AC_DISABLE_STATIC

AM_INIT_AUTOMAKE(foreign subdir-objects)

LT_INIT

AC_PROG_CXX
AM_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG
# Instruct autoconf to use the c++ compiler for subsequent compile checks
AC_LANG_PUSH([C++])
AC_CHECK_HEADER([glm/glm.hpp], [], [AC_MSG_ERROR([Could not find glm/glm.hpp.])])
AC_CHECK_HEADER([lpsolve/lp_lib.h], [], [AC_MSG_ERROR([Could not find lpsolve/lp_lib.h.])])
AC_LANG_POP([C++])
AX_BOOST_BASE
AX_BOOST_IOSTREAMS
PKG_CHECK_MODULES([eigen3], [eigen3])

AC_CONFIG_FILES(Makefile lib/Makefile test/Makefile)
AC_OUTPUT

