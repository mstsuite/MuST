## Copyright (C) 2008 T. Burnus
## Copyright (C) 2015 M. Oliveira
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
##

#
# Check for Fortran 2003 iso_c_bindings support
# ------------------------------------

AC_DEFUN([ACX_FC_ISO_C_BINDING], [

AC_MSG_CHECKING([for Fortran 2003 iso_c_binding])

testprog="AC_LANG_PROGRAM([],[
  use iso_c_binding
  implicit none
  type(c_ptr) :: ptr
  ptr = c_null_ptr
  if (c_associated(ptr)) stop 3])"

acx_iso_c_binding_ok=no
AC_LINK_IFELSE($testprog, [acx_iso_c_binding_ok=yes], [])

AC_MSG_RESULT([$acx_iso_c_binding_ok])
if test x"$acx_iso_c_binding_ok" = xyes; then
  AC_DEFINE(ISO_C_BINDING, 1, [compiler supports Fortran 2003 iso_c_binding])
  $1
else
  $2
fi
])
