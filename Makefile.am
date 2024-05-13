AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 --install

VERSION_MAJOR   = @VERSION_MAJOR@
VERSION_MINOR   = @VERSION_MINOR@
VERSION_PATCH   = @VERSION_PATCH@

src_main_c = src/ex1.c src/doubleexchange.c src/readgraphmllib.c  src/utils.c src/functions.c src/get_proj_general.c
src_main_h = src/doubleexchange.h src/readgraphmllib.h src/utils.h src/functions.h src/get_proj_general.h

bin_PROGRAMS = $(top_builddir)/bin/ex1
__top_builddir__bin_ex1_SOURCES = $(src_main_c)
__top_builddir__bin_ex1_CPPFLAGS = $(AM_CPPFLAGS)
__top_builddir__bin_ex1_CFLAGS = $(AM_CFLAGS)
__top_builddir__bin_ex1_CFLAGS += "-I${top_srcdir}/include"  "-I${top_srcdir}/src"
__top_builddir__bin_ex1_LDFLAGS = $(AM_LDFLAGS)
