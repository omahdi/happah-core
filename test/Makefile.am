# Copyright 2017
#   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

AM_CPPFLAGS = -I$(top_srcdir)/lib -I/usr/include/eigen3 -std=c++1z
LDADD = $(top_srcdir)/lib/libhappah.la -lboost_iostreams

check_PROGRAMS = \
     happah_geometries_SurfaceSplineBEZ_DeltasEnumerator
happah_geometries_SurfaceSplineBEZ_DeltasEnumerator_SOURCES = \
     happah/geometries/SurfaceSplineBEZ/DeltasEnumerator.cpp

TESTS = $(check_PROGRAMS)

