subroutine print_version(myunit)
use KindParamModule, only : IntKind
implicit none
   integer(kind=IntKind), intent(in) :: myunit

   write(myunit,'(''#'',79(''-''))')
   write(myunit,'(''#'',/,''#'')')
   write(myunit,'(''#'',11x,a)')'********************************************************'
   write(myunit,'(''#'',11x,a)')'*                                                      *'
   write(myunit,'(''#'',11x,a)')'*    Full-Potential Multiple Scattering Theory Based   *'
   write(myunit,'(''#'',11x,a)')'*                                                      *'
   write(myunit,'(''#'',11x,a)')'*  Ab Initio Electronic Structure Calculation Package  *'
   write(myunit,'(''#'',11x,a)')'*                                                      *'
!   write(myunit,'(12x,a)')'*                    Version 2.1                       *'
!   write(myunit,'(12x,a)')'*                                                      *'
   write(myunit,'(''#'',11x,a)')'********************************************************'
   write(myunit,'(''#'',/,''#'')')
!
#include "print_version_include.h"
#include "git_version.h"
!
   write(myunit,'(''#'',79(''-''))')

!   write(myunit,'(/,80(''=''))')

end subroutine print_version
