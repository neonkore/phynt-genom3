#
# Copyright (c) 2018 LAAS/CNRS
# All rights reserved.
#
# Redistribution and use  in source  and binary  forms,  with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   1. Redistributions of  source  code must retain the  above copyright
#      notice and this list of conditions.
#   2. Redistributions in binary form must reproduce the above copyright
#      notice and  this list of  conditions in the  documentation and/or
#      other materials provided with the distribution.
#
#                                      Anthony Mallet on Mon Jun 11 2018
#
lib_LTLIBRARIES = libphynt_codels.la

libphynt_codels_la_SOURCES=	phynt_c_types.h
libphynt_codels_la_SOURCES+=	phynt_codels.c
libphynt_codels_la_SOURCES+=	phynt_main_codels.c

libphynt_codels_la_SOURCES+=	observer.cc
libphynt_codels_la_SOURCES+=	admittance.cc
libphynt_codels_la_SOURCES+=	codels.h

libphynt_codels_la_CPPFLAGS=	$(requires_CFLAGS)
libphynt_codels_la_LDFLAGS=	$(requires_LIBS)
libphynt_codels_la_LDFLAGS+=	-release $(PACKAGE_VERSION)


# idl  mappings
BUILT_SOURCES=	phynt_c_types.h
CLEANFILES=	${BUILT_SOURCES}
phynt_c_types.h: ${top_srcdir}/phynt.gen
	${GENOM3}  mappings \
	  -MD -MF .deps/$@.d -MT $@ --signature -l c $< >$@

-include .deps/phynt_c_types.h.d
