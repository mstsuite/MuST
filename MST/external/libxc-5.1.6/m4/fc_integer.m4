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
# Check size of a fortran integer
# ----------------------------------
AC_DEFUN([ACX_FC_INTEGER_SIZE],[
  AC_REQUIRE([AC_PROG_FC])

  if test -z "$FC_INTEGER_SIZE"; then
    AC_MSG_CHECKING([for the size of a Fortran integer])
    AC_RUN_IFELSE([AC_LANG_PROGRAM([],[
  integer    :: i
  integer(8) :: i8

  i8 = huge(i)

  select case(i8)
  case(127_8);                 i = 1
  case(32767_8);               i = 2
  case(2147483647_8);          i = 4
  case(9223372036854775807_8); i = 8
  case default; write(*,'(a,i20)') "unrecognized size ", i8; stop 999
  end select

  open(1, file='conftest.out')
  write(1,'(i1)') i
  close(1)
])],
      [ac_fcintegersize=`cat conftest.out`],
      [AC_MSG_FAILURE(f90 program to find the size of a Fortran integer failed)],
      [ac_fcintegersize=4; echo -n "cross-compiling; assuming... "])
    AC_MSG_RESULT([${ac_fcintegersize} bytes])
  else
    ac_fcintegersize=$FC_INTEGER_SIZE
  fi
  AC_DEFINE_UNQUOTED(FC_INTEGER_SIZE, ${ac_fcintegersize}, [The size of a Fortran integer])
])

################################################
# Check which C type corresponds to Fortran int
# ----------------------------------
AC_DEFUN([ACX_CC_FORTRAN_INT],[
  AC_MSG_CHECKING([for which C type corresponds to Fortran integer])
  AC_REQUIRE([ACX_FC_INTEGER_SIZE])
  AC_REQUIRE([AC_PROG_CC])
  if test -z "$CC_FORTRAN_INT"; then
    AC_LANG_PUSH([C])
    AC_RUN_IFELSE([AC_LANG_PROGRAM([
#include <stdio.h>
],[
  FILE* fp;
  fp = fopen("conftest.out", "w");
  if(${ac_fcintegersize} == sizeof(char))
    fprintf(fp, "char");
  else if(${ac_fcintegersize} == sizeof(short))
    fprintf(fp, "short");
  else if(${ac_fcintegersize} == sizeof(int))
    fprintf(fp, "int");
  else if(${ac_fcintegersize} == sizeof(long))
    fprintf(fp, "long");
  else
    return 1;
])],
      [ac_ccfortranint=`cat conftest.out`],
      [AC_MSG_FAILURE(C program failed to find the C type of a Fortran integer)],
      [ac_ccfortranint="int"; echo -n "cross-compiling; assuming... "])
    AC_LANG_POP([C])
    AC_MSG_RESULT([${ac_ccfortranint}])
  else
    ac_ccfortranint=$CC_FORTRAN_INT
  fi
  AC_DEFINE_UNQUOTED(CC_FORTRAN_INT, ${ac_ccfortranint}, [The C type of a Fortran integer])
])
