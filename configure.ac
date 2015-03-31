# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.62])
DUNE_AC_INIT # gets module version from dune.module file
AC_CONFIG_SRCDIR([src/immoplate.cc])
AC_CONFIG_HEADERS([config.h])
#AC_CHECK_HEADERS([ getopt.h ])

# we need no more than the standard DE-stuff
# this module depends on dune-common dune-geometry dune-grid
# this implies checking for [dune-common], [dune-geometry], [dune-grid]
DUNE_CHECK_ALL


# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

#AC_ARG_WITH(gsl,
#            [  --with-gsl=DIR        Directory where the GSL is installed (optional)],
#            GSL_DIR="$withval")
#
#GSL_CONFIG=
#AC_PATH_PROG([GSL_CONFIG], [gsl-config],, [$GSL_DIR/bin:$PATH])
#                                                                                                                
#GSL_INCLUDES=
#GSL_LIBS=
                                                                                                                 
#if test -x "$GSL_CONFIG"; then
#  GSL_INCLUDES=`$GSL_CONFIG --cflags`
#  GSL_LIBS=`$GSL_CONFIG --libs`
#else
#  AC_ERROR([GSL not found])
#fi
#AC_SUBST(GSL_INCLUDES)
#AC_SUBST(GSL_LIBS)

#AC_CHECK_FUNC([getopt_long], [HAVE_GETOPT_LONG="yes"])
#if test "x$HAVE_GETOPT_LONG" = "xyes" ; then
#  AC_DEFINE(HAVE_GETOPT_LONG, [], [Define to 1 if you have the 'getopt_long' function])
#  GETOPT_LONG_LIBS=""
#  AC_SUBST(GETOPT_LONG_LIBS)
#else
dnl check for gnugetopt library
#  AC_CHECK_LIB([gnugetopt], [getopt_long], [HAVE_GETOPT_LONG="yes"])
#  if test "x$HAVE_GETOPT_LONG" = "xyes" ; then
#    AC_DEFINE(HAVE_GETOPT_LONG, [], [Define to 1 if you have the 'getopt_long' function])
#    GETOPT_LONG_LIBS="-lgnugetopt"
#    AC_SUBST(GETOPT_LONG_LIBS)
#  fi
#fi

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/tests/Makefile
  cmake/Makefile
  cmake/modules/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/immo/Makefile
  m4/Makefile
  dune-immo.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
