# dependencies.m4 generated by dunecontrol

m4_define([DUNE_AC_INIT],[
  AC_INIT([dune-immo], [0.3.0], [valiska@fz-juelich.de])
  AM_INIT_AUTOMAKE([foreign 1.9 tar-pax])
  m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([no])])
  AC_SUBST([DUNE_MOD_VERSION], [0.3.0])
  AC_SUBST([DUNE_MOD_NAME], [dune-immo])
  AC_SUBST([DUNE_MAINTAINER_NAME], ["valiska@fz-juelich.de"])
  DUNE_PARSE_MODULE_VERSION([dune-immo], [0.3.0])
  REQUIRES=" dune-common   dune-geometry   dune-grid  "
  AC_SUBST(REQUIRES, [])
  AC_CONFIG_MACRO_DIR([m4])
])

AC_DEFUN([DUNE_CHECK_MOD_DEPENDENCIES], [
  ### check dependency dune-common
  # invoke checks required by this module
  AC_REQUIRE([DUNE_COMMON_CHECKS])
  # invoke check for this module
  AC_REQUIRE([DUNE_COMMON_CHECK_MODULE])
  if test x$with_dune_common = xno; then
    AC_MSG_ERROR([could not find required module _dune_name])
  fi
  ### check dependency dune-geometry
  # invoke checks required by this module
  AC_REQUIRE([DUNE_GEOMETRY_CHECKS])
  # invoke check for this module
  AC_REQUIRE([DUNE_GEOMETRY_CHECK_MODULE])
  if test x$with_dune_geometry = xno; then
    AC_MSG_ERROR([could not find required module _dune_name])
  fi
  ### check dependency dune-grid
  # invoke checks required by this module
  AC_REQUIRE([DUNE_GRID_CHECKS])
  # invoke check for this module
  AC_REQUIRE([DUNE_GRID_CHECK_MODULE])
  if test x$with_dune_grid = xno; then
    AC_MSG_ERROR([could not find required module _dune_name])
  fi
  ### invoke checks for dune-immo
  AC_REQUIRE([DUNE_IMMO_CHECKS])
])
