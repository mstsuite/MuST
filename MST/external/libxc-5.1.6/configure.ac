## Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
##               2019 X. Andrade
##               2020 Susi Lehtola
##
## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

AC_INIT([libxc],[5.1.6],[libxc@tddft.org],[libxc],[http://www.tddft.org/programs/Libxc])
AC_CONFIG_SRCDIR([src/xc.h])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Check if libm needs to be explicitly linked
AC_SEARCH_LIBS([sinh], [m])
AC_SEARCH_LIBS([atan2], [m])
AC_SEARCH_LIBS([tanh], [m])
AC_SEARCH_LIBS([cosh], [m])
AC_SEARCH_LIBS([sin], [m])
AC_SEARCH_LIBS([atan], [m])
AC_SEARCH_LIBS([exp], [m])
AC_SEARCH_LIBS([cos], [m])
AC_SEARCH_LIBS([log], [m])
AC_SEARCH_LIBS([pow], [m])
AC_SEARCH_LIBS([round], [m])
AC_SEARCH_LIBS([sqrt], [m])
AC_SEARCH_LIBS([erf], [m])

# Library versioning (C:R:A == current:revision:age)
# See the libtool manual for an explanation of the numbers
#
# libxc-2.0.0   1:0:0
# libxc-2.0.1   1:1:0
# libxc-2.0.2   1:2:0
# libxc-2.0.3   1:3:0
# libxc-2.1.0   2:0:1
# libxc-2.1.1   2:1:1
# libxc-2.1.2   2:2:1
# libxc-2.1.3   2:3:1
# libxc-2.2.0   3:0:0
# libxc-2.2.1   3:1:0
# libxc-2.2.3   3:3:0
# libxc-3.0.0   4:0:0
# libxc-3.0.1   4:1:0
# libxc-4.0.0   5:0:0
# libxc-4.0.1   5:1:0
# libxc-4.0.2   5:2:0
# libxc-4.0.3   5:3:0
# libxc-4.0.4   5:4:0
# libxc-4.0.5   5:5:0
# libxc-4.1.0   6:0:1
# libxc-4.1.1   6:1:1
# libxc-4.2.0   7:0:2
# libxc-4.2.1   7:1:2
# libxc-4.2.2   7:2:2
# libxc-4.2.3   7:3:2
# libxc-4.3.0   8:0:3
# libxc-4.3.2   8:1:3
# libxc-4.3.2   8:2:3
# libxc-4.3.3   8:3:3
# libxc-4.3.4   8:4:3
# libxc-5.0.0   9:0:0
# libxc-5.1.0  10:0:1
# libxc-5.1.1  10:1:1
# libxc-5.1.2  10:2:1
# libxc-5.1.3  10:3:1
# libxc-5.1.4  10:4:1
# libxc-5.1.5  10:5:1
# libxc-5.1.6  10:6:1
# libxc-6.0.0  11:0:0
# 
# Note that libtool versioning was not used prior to libxc-2.0.
# Because of this, we will start counting at 1.
#
# How to update library version number
# ====================================
# 
# C: increment if the interface has additions, changes, removals.
#
# R: increment any time the source changes; set to 0 if you
# incremented CURRENT
#
# A: increment if any interfaces have been added; set to 0 if any
# interfaces have been removed. removal has precedence over adding,
# so set to 0 if both happened.

XC_CURRENT=10
XC_REVISION=6
XC_AGE=1

XC_LT_VERSION="${XC_CURRENT}:${XC_REVISION}:${XC_AGE}"
AC_SUBST(XC_LT_VERSION)

# Split VERSION into XC_VERSION_MAJOR and XC_VERSION_MINOR
# Follows AX_SPLIT_VERSION macro from AC-Archive
# Rhys Ulerich <rhys.ulerich@gmail.com>
AC_PROG_SED
XC_MAJOR_VERSION=`echo "$VERSION" | $SED 's/\([[^.]][[^.]]*\).*/\1/'`
XC_MINOR_VERSION=`echo "$VERSION" | $SED 's/[[^.]][[^.]]*.\([[^.]][[^.]]*\).*/\1/'`
XC_MICRO_VERSION=`echo "$VERSION" | $SED 's/[[^.]][[^.]]*.[[^.]][[^.]]*.\(.*\)/\1/'`
AC_SUBST(XC_MAJOR_VERSION)
AC_SUBST(XC_MINOR_VERSION)
AC_SUBST(XC_MICRO_VERSION)


# Installation prefix by default
AC_PREFIX_DEFAULT([/opt/etsf])

AC_DISABLE_SHARED
LT_PREREQ([2.2])
LT_INIT

# enable cuda first, since this is required for the compiler selection
AC_ARG_ENABLE(cuda, AS_HELP_STRING([--enable-cuda], [Cuda support (experimental)]))

if test x"${enable_cuda}" = x"yes" ; then
  AC_DEFINE(HAVE_CUDA, 1, [defined if cuda support is enabled])
fi

AM_CONDITIONAL([HAVE_CUDA], [test x"${enable_cuda}" = x"yes"])

# Checks for programs.
AC_PROG_CC

#we may need to set the linker for the cuda case, so pass that to automake
if test x"$CCLD" = x ; then
   CCLD=$CC
fi
AC_SUBST(CCLD)

AC_LANG_PREPROC
AM_PROG_CC_C_O
AC_C_INLINE

# Checks for header files.
AC_HEADER_STDC
AC_FUNC_ALLOCA

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_SIZEOF(void*)

# Need BASH to run tests
AC_CHECK_PROG([BASH_SHELL],[bash],[yes])

# Check for presence of some C99 libm functions
AC_CHECK_LIB([m],[cbrt],[AC_DEFINE([HAVE_CBRT],[1],[libm includes cbrt])])
AC_CHECK_LIB([m],[feenableexcept],[AC_DEFINE([HAVE_FEENABLEEXCEPT],[1],[libm includes feenableexcept])])

AC_ARG_WITH([custom-fmoddir],
	    [AS_HELP_STRING([--with-custom-fmoddir=DIR],
			    [Instalation directory for Fortran module files (includedir by default)])],
	    [AC_SUBST([FMODDIR],["$withval"])],
	    [AC_SUBST([FMODDIR],['${includedir}'])])

# Enable/disable derivatives
AC_ARG_ENABLE([vxc],
	      AS_HELP_STRING([--disable-vxc], [don't build first derivatives]),
	      [ac_cv_build_vxc=$enableval],
	      [ac_cv_build_vxc=yes])
if test x"${ac_cv_build_vxc}" = x"no"; then
  AC_DEFINE(XC_DONT_COMPILE_VXC, 1, [defined if first derivatives are disabled])
  AC_DEFINE(XC_DONT_COMPILE_FXC, 1, [defined if second derivatives are disabled])
  AC_DEFINE(XC_DONT_COMPILE_KXC, 1, [defined if third derivatives are disabled])
  AC_DEFINE(XC_DONT_COMPILE_LXC, 1, [defined if fourth derivatives are disabled])
fi

AC_ARG_ENABLE([fxc],
	      AS_HELP_STRING([--disable-fxc], [don't build second derivatives]),
	      [ac_cv_build_fxc=$enableval],
	      [ac_cv_build_fxc=yes])
if test x"${ac_cv_build_fxc}" = x"no"; then
  AC_DEFINE(XC_DONT_COMPILE_FXC, 1, [defined if second derivatives are disabled])
  AC_DEFINE(XC_DONT_COMPILE_KXC, 1, [defined if third derivatives are disabled])
  AC_DEFINE(XC_DONT_COMPILE_LXC, 1, [defined if fourth derivatives are disabled])
fi

AC_ARG_ENABLE([kxc],
	      AS_HELP_STRING([--disable-kxc], [don't build third derivatives]),
	      [ac_cv_build_kxc=$enableval],
	      [ac_cv_build_kxc=no])
if test x"${ac_cv_build_kxc}" = x"no"; then
  AC_DEFINE(XC_DONT_COMPILE_KXC, 1, [defined if third derivatives are disabled])
  AC_DEFINE(XC_DONT_COMPILE_LXC, 1, [defined if fourth derivatives are disabled])
fi

AC_ARG_ENABLE([lxc],
	      AS_HELP_STRING([--disable-lxc], [don't build fourth derivatives]),
	      [ac_cv_build_lxc=$enableval],
	      [ac_cv_build_lxc=no])
if test x"${ac_cv_build_lxc}" = x"no"; then
  AC_DEFINE(XC_DONT_COMPILE_LXC, 1, [defined if fourth derivatives are disabled])
fi

# build/install fortran bits unless --disable-fortran is used
AC_ARG_ENABLE([fortran],
	      AS_HELP_STRING([--disable-fortran], [don't build/install fortran components]),
	      [ac_cv_build_fortran=$enableval],
	      [ac_cv_build_fortran=yes])
AM_CONDITIONAL([ENABLE_FORTRAN], [test $ac_cv_build_fortran = yes])

if test $ac_cv_build_fortran = yes; then
  HAVE_FORTRAN=1
  AC_DEFINE(HAVE_FORTRAN, [1], [Defined if libxc is compiled with fortran support])

  dnl try to find out what is the default FORTRAN 90 compiler
  acx_save_fcflags="${FCFLAGS}"
  AC_PROG_FC([], [Fortran 90])
  if test x"$FC" = x; then
    AC_MSG_ERROR([could not find Fortran 90 compiler])
  fi

  AC_LANG_PUSH(Fortran)
  AC_FC_SRCEXT(f90)
  FCFLAGS="${acx_save_fcflags}"

  dnl Fortran default flags
  ACX_FCFLAGS

  dnl libxc needs the preprocessor. The result goes to FCCPP
  ACX_FCCPP

  dnl how Fortran mangles function names
  AC_FC_WRAPPERS

  dnl check whether we have or not a compiler that allows for very long lines...
  ACX_LONG_FORTRAN_LINES

  dnl check whether the Fortran 90 compiler accepts line markers cast by the preprocessor.
  ACX_F90_ACCEPTS_LINE_NUMBERS

  AX_F90_MODULE_EXTENSION

  dnl Need to know the size of a Fortran integer
  ACX_FC_INTEGER_SIZE
  ACX_CC_FORTRAN_INT

  dnl check whether the Fortran compiler supports Fortran 2003 iso_c_binding
  ACX_FC_ISO_C_BINDING([], [AC_MSG_ERROR([Could not find Fortran 2003 iso_c_binding. Fortran 2003 interface can not be compiled.])])

AC_LANG_POP([Fortran])
else
  HAVE_FORTRAN=0
fi
AC_SUBST(ax_cv_f90_modext)
AM_CONDITIONAL(F90_MOD_UPPERCASE, [test x$ax_f90_mod_uppercase = xyes])
AM_CONDITIONAL(COMPILE_FORTRAN, test x${HAVE_FORTRAN} = x1)
AM_CONDITIONAL([COMPILE_FORTRAN03], [test x$ac_cv_build_fortran03 = xyes])


AC_CONFIG_FILES([Makefile
  xc_version.h
  src/Makefile
  src/maple2c/Makefile
  src/maple2c/lda_exc/Makefile
  src/maple2c/lda_vxc/Makefile
  src/maple2c/gga_exc/Makefile
  src/maple2c/gga_vxc/Makefile
  src/maple2c/mgga_exc/Makefile
  src/maple2c/mgga_vxc/Makefile
  testsuite/Makefile
  testsuite/pytest/Makefile
  pkgconfig/Makefile
  pkgconfig/libxc.pc
  pkgconfig/libxcf90.pc
  pkgconfig/libxcf03.pc
  scripts/Makefile
  maple/Makefile
  maple/lda_exc/Makefile
  maple/lda_vxc/Makefile
  maple/gga_exc/Makefile
  maple/gga_vxc/Makefile
  maple/mgga_exc/Makefile
  maple/mgga_vxc/Makefile
  pylibxc/Makefile
  cmake/Makefile
])

AC_OUTPUT

