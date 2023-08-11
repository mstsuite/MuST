dnl @synopsis AX_F90_MODULE_EXTENSION
dnl
dnl Find Fortran 90 modules file extension. The module extension is
dnl stored in the cached variable ax_cv_f90_modext, or "unknown" if the
dnl extension cannot be found.
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2005-06-17
dnl @license AllPermissive

AC_DEFUN([AX_F90_MODULE_EXTENSION],[
AC_CACHE_CHECK([Fortran 90 modules extension],
ax_cv_f90_modext,
[AC_LANG_PUSH(Fortran)
ax_f90_mod_uppercase=no
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
AC_COMPILE_IFELSE([module conftest_module
   contains
   subroutine conftest_routine
   write(*,'(a)') 'gotcha!'
   end subroutine conftest_routine
   end module conftest_module
  ],
  [
dnl prune out 'submodule' file .smod from gcc version 6 and later
   ax_cv_f90_modext=`ls | grep -v smod | sed -n 's,conftest_module\.,,p'`
   if test x"$ax_cv_f90_modext" = x ; then
dnl Some F90 compilers put module filename in uppercase letters
     ax_cv_f90_modext=`ls | grep -v SMOD | sed -n 's,CONFTEST_MODULE\.,,p'`
     if test x"$ax_cv_f90_modext" = x ; then
       ax_cv_f90_modext=unknown
     else
       ax_f90_mod_uppercase=yes
     fi
   fi
  ],
  [ax_cv_f90_modext=unknown])
cd ..
rm -fr tmpdir_$i
AC_LANG_POP(Fortran)
])])
