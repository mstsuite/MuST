## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
##

################################################
# Check whether the compiler accepts very long lines.
# ----------------------------------
AC_DEFUN([ACX_LONG_FORTRAN_LINES],
[AC_MSG_CHECKING([whether the compiler accepts very long lines])
AC_COMPILE_IFELSE( AC_LANG_PROGRAM( [], [
write(*, *) '456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678904567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789001234567890123456789045678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890'
 ]), 
 [acx_long_lines_ok=yes; AC_DEFINE(LONG_LINES, 1, [compiler supports long lines])], [acx_long_lines_ok=no])
AC_SUBST([LONG_LINES], [$acx_long_lines_ok])
AC_MSG_RESULT($acx_long_lines_ok)
])


################################################
# Check whether the compiler accepts preprocessor "# line-number" lines.
# ----------------------------------
AC_DEFUN([ACX_F90_ACCEPTS_LINE_NUMBERS],
[
AC_MSG_CHECKING([whether the compiler accepts "line-number" lines cast by the preprocessor])
AC_COMPILE_IFELSE(
    AC_LANG_PROGRAM( [], [# 1]),
    [acx_f90_accepts_line_numbers_ok=yes
    AC_DEFINE(F90_ACCEPTS_LINE_NUMBERS, 1, [compiler supports line-number lines])],
    [acx_f90_accepts_line_numbers_ok=no])
AC_SUBST(F90_ACCEPTS_LINE_NUMBERS, $acx_f90_accepts_line_numbers_ok)
AC_MSG_RESULT($acx_f90_accepts_line_numbers_ok)
]
)

################################################
# Check for the presence of a given function in Fortran.
# It substitutes AC_CHECK_FUNC, since the latter
# seems to fail with some autotools versions, due to a call to some broken
# version of AC_LANG_FUNC_LINK_TRY.
AC_DEFUN([ACX_FORTRAN_CHECK_FUNC],
[
AC_MSG_CHECKING([for $1])
AC_LANG_PUSH(Fortran)dnl
AC_LINK_IFELSE([AC_LANG_CALL([], [$1])], 
[
acx_fortran_check_func=yes
AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_$1]),1, [Define if the $1 function can be called from Fortran])], 
[
acx_fortran_check_func=no
])dnl
AC_LANG_POP(Fortran)dnl
AC_MSG_RESULT($acx_fortran_check_func)
])


################################################
# AC_LANG_FUNC_LINK_TRY(Fortran)(FUNCTION)
# ----------------------------------
m4_define([AC_LANG_FUNC_LINK_TRY(Fortran)],
[AC_LANG_PROGRAM([], [call [$1]])])

################################################ 
# Fortran preprocessing
# --------------------------- 

# like built-in AC_EGREP_CPP, only using FCCPP, .F90 file, and regular grep
# some cpp's will behave differently on .F90 and on .c files
# using $GREP or $EGREP would imply AC_REQUIRE(AC_PROG_GREP) which leads to this warning:
# configure.ac:189: *GNU* is m4_require'd but not m4_defun'd
AC_DEFUN([ACX_GREP_FCCPP],[
     AC_ARG_VAR(FCCPP, [Fortran preprocessor])

     echo "$2" > conftest.F90

     # see if the attempt to preprocess raises an error
     if ! $FCCPP conftest.F90 > /dev/null 2>&5 ; then
       $4
     fi

     if (eval "$FCCPP conftest.F90" 2>&5) |
       grep -q "$1" 2>&5; then :
       $3
     else
       $4
    fi
    rm -f conftest*
])

AC_DEFUN([ACX_FCCPP],[
     # "gcc -E -x c" means treat the file as if it were C. For some reason, when gcc identifies the source
     # as Fortran, it will not concatenate tokens in preprocessing, so we must trick it.
     for FCCPP_base in "$FCCPP" "/lib/cpp" "$CPP" "$CPP -x c" "`which cpp`"; do
         # cycle if blank
         if test -z "$FCCPP_base"; then
           continue
         fi

         for FCCPP in "$FCCPP_base" "$FCCPP_base -ansi" "$FCCPP_base -C -ffreestanding"; do
           AC_MSG_CHECKING([whether $FCCPP is usable for Fortran preprocessing])
	   acx_fpp_ok=yes

      	   ACX_GREP_FCCPP([anything], AC_LANG_PROGRAM([],[anything]),
	     [], [acx_fpp_ok=no; AC_MSG_RESULT([preprocessor cannot be run]); break])
	     # very unlikely that adding -ansi will allow it to be run at all

      	   ACX_GREP_FCCPP([hi], AC_LANG_PROGRAM([],[
#define ADD_I(x) x ## i
ADD_I(h)]),
	     [], [acx_fpp_ok=no; AC_MSG_RESULT([preprocessor does not concatenate tokens])])

           # in Fortran this is string concatenation, must not be stripped
	   # some cpp's (e.g. icc -E -ansi) might actually insert a space between // too which is not acceptable
           ACX_GREP_FCCPP([rout // ine], AC_LANG_PROGRAM([],[
#define PUSH_SUB(x) x // ine
PUSH_SUB(rout)]),
	     [], [acx_fpp_ok=no; AC_MSG_RESULT([preprocessor mangles C++ style comment])])

	  if test x"$acx_fpp_ok" = xyes; then
            AC_MSG_RESULT([yes])
	    break
	  fi
       done
       if test x"$acx_fpp_ok" = xyes; then
	  break
       fi
     done

     if test x"$acx_fpp_ok" = xno; then
     	AC_MSG_ERROR([Could not find preprocessor usable for Fortran.])
     fi

     AC_SUBST(FCCPP)
])
