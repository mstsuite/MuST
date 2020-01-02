module Matrix3dModule 
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : TEN2m8
!
public :: invm3, &
          determ3
!
   interface invm3
      module procedure invm3_real1, invm3_real2,                      &
                       invm3_cmplx1, invm3_cmplx2
   end interface
!
   interface determ3
      module procedure determ3_real1, determ3_real2,                  &
                       determ3_cmplx1, determ3_cmplx2
   end interface
!
private
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine invm3_real1(xm,xminv,ierr)
!  ===================================================================
!
   implicit   none
!
   integer (kind=IntKind), intent(out) :: ierr
!
   real (kind=RealKind), intent(in) :: xm(9)
   real (kind=RealKind), intent(out) :: xminv(9)
!
   real (kind=RealKind) :: det
!
   det=determ3(xm)
!
   if (abs(det) .lt. ten2m8) then
      ierr=1
      return
   else
      ierr=0
      xminv(1)=+(xm(5)*xm(9)-xm(8)*xm(6))/det
      xminv(2)=-(xm(2)*xm(9)-xm(8)*xm(3))/det
      xminv(3)=+(xm(2)*xm(6)-xm(3)*xm(5))/det
      xminv(4)=-(xm(4)*xm(9)-xm(7)*xm(6))/det
      xminv(5)=+(xm(1)*xm(9)-xm(7)*xm(3))/det
      xminv(6)=-(xm(1)*xm(6)-xm(3)*xm(4))/det
      xminv(7)=+(xm(4)*xm(8)-xm(7)*xm(5))/det
      xminv(8)=-(xm(1)*xm(8)-xm(2)*xm(7))/det
      xminv(9)=+(xm(1)*xm(5)-xm(2)*xm(4))/det
   end if
!
   end subroutine invm3_real1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine invm3_real2(xm,xminv,ierr)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(out) :: ierr
!
   real (kind=RealKind), intent(in) :: xm(3,3)
   real (kind=RealKind), intent(out) :: xminv(3,3)
!
   real (kind=RealKind) :: det
!
   det=determ3(xm)
!
   if (abs(det) .lt. ten2m8) then
      ierr=1
      return
   else
      ierr=0
      xminv(1,1)=+(xm(2,2)*xm(3,3)-xm(2,3)*xm(3,2))/det
      xminv(2,1)=-(xm(2,1)*xm(3,3)-xm(2,3)*xm(3,1))/det
      xminv(3,1)=+(xm(2,1)*xm(3,2)-xm(3,1)*xm(2,2))/det
      xminv(1,2)=-(xm(1,2)*xm(3,3)-xm(1,3)*xm(3,2))/det
      xminv(2,2)=+(xm(1,1)*xm(3,3)-xm(1,3)*xm(3,1))/det
      xminv(3,2)=-(xm(1,1)*xm(3,2)-xm(3,1)*xm(1,2))/det
      xminv(1,3)=+(xm(1,2)*xm(2,3)-xm(1,3)*xm(2,2))/det
      xminv(2,3)=-(xm(1,1)*xm(2,3)-xm(2,1)*xm(1,3))/det
      xminv(3,3)=+(xm(1,1)*xm(2,2)-xm(2,1)*xm(1,2))/det
   end if
!
   end subroutine invm3_real2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine invm3_cmplx1(xm,xminv,ierr)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(out) :: ierr
!
   Complex (kind=CmplxKind), intent(in) :: xm(9)
   Complex (kind=CmplxKind), intent(out) :: xminv(9)
!
   Complex (kind=CmplxKind) :: det
!
   det=determ3(xm)
!
   if (abs(det) .lt. ten2m8) then
      ierr=1
      return
   else
      ierr=0
      xminv(1)=+(xm(5)*xm(9)-xm(8)*xm(6))/det
      xminv(2)=-(xm(2)*xm(9)-xm(8)*xm(3))/det
      xminv(3)=+(xm(2)*xm(6)-xm(3)*xm(5))/det
      xminv(4)=-(xm(4)*xm(9)-xm(7)*xm(6))/det
      xminv(5)=+(xm(1)*xm(9)-xm(7)*xm(3))/det
      xminv(6)=-(xm(1)*xm(6)-xm(3)*xm(4))/det
      xminv(7)=+(xm(4)*xm(8)-xm(7)*xm(5))/det
      xminv(8)=-(xm(1)*xm(8)-xm(2)*xm(7))/det
      xminv(9)=+(xm(1)*xm(5)-xm(2)*xm(4))/det
   end if
!
   end subroutine invm3_cmplx1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine invm3_cmplx2(xm,xminv,ierr)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(out) :: ierr
!
   Complex (kind=CmplxKind), intent(in) :: xm(3,3)
   Complex (kind=CmplxKind), intent(out) :: xminv(3,3)
!
   Complex (kind=CmplxKind) :: det
!
   det=determ3(xm)
!
   if (abs(det) .lt. ten2m8) then
      ierr=1
      return
   else
      ierr=0
      xminv(1,1)=+(xm(2,2)*xm(3,3)-xm(2,3)*xm(3,2))/det
      xminv(2,1)=-(xm(2,1)*xm(3,3)-xm(2,3)*xm(3,1))/det
      xminv(3,1)=+(xm(2,1)*xm(3,2)-xm(3,1)*xm(2,2))/det
      xminv(1,2)=-(xm(1,2)*xm(3,3)-xm(1,3)*xm(3,2))/det
      xminv(2,2)=+(xm(1,1)*xm(3,3)-xm(1,3)*xm(3,1))/det
      xminv(3,2)=-(xm(1,1)*xm(3,2)-xm(3,1)*xm(1,2))/det
      xminv(1,3)=+(xm(1,2)*xm(2,3)-xm(1,3)*xm(2,2))/det
      xminv(2,3)=-(xm(1,1)*xm(2,3)-xm(2,1)*xm(1,3))/det
      xminv(3,3)=+(xm(1,1)*xm(2,2)-xm(2,1)*xm(1,2))/det
   end if
!
   end subroutine invm3_cmplx2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function determ3_real1(xm)
!  ===================================================================
   implicit  none
!
   real (kind=RealKind), intent(in) :: xm(9)
   real (kind=RealKind) :: determ3_real1
!
!  *******************************************************************
!
!  determ = the determinant of the 3x3 matrix xm.
!
!  *******************************************************************
   determ3_real1 = xm(1)*(xm(5)*xm(9)-xm(8)*xm(6))      &
                  -xm(4)*(xm(2)*xm(9)-xm(8)*xm(3))      &
                  +xm(7)*(xm(2)*xm(6)-xm(5)*xm(3))
!
   end function determ3_real1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function determ3_real2(xm)
!  ===================================================================
   implicit  none
!
   real (kind=RealKind), intent(in) :: xm(3,3)
   real (kind=RealKind) :: determ3_real2
!
!  *******************************************************************
!
!  determ = the determinant of the 3x3 matrix xm.
!
!  *******************************************************************
   determ3_real2 = xm(1,1)*(xm(2,2)*xm(3,3)-xm(2,3)*xm(3,2))      &
                  -xm(1,2)*(xm(2,1)*xm(3,3)-xm(2,3)*xm(3,1))      &
                  +xm(1,3)*(xm(2,1)*xm(3,2)-xm(2,2)*xm(3,1))
!
   end function determ3_real2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function determ3_cmplx1(xm)
!  ===================================================================
   implicit  none
!
   complex (kind=CmplxKind), intent(in) :: xm(9)
   complex (kind=CmplxKind) :: determ3_cmplx1
!
!  *******************************************************************
!
!  determ = the determinant of the 3x3 matrix xm.
!
!  *******************************************************************
   determ3_cmplx1 = xm(1)*(xm(5)*xm(9)-xm(8)*xm(6))      &
                   -xm(4)*(xm(2)*xm(9)-xm(8)*xm(3))      &
                   +xm(7)*(xm(2)*xm(6)-xm(5)*xm(3))
!
   end function determ3_cmplx1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function determ3_cmplx2(xm)
!  ===================================================================
   implicit  none
!
   complex (kind=CmplxKind), intent(in) :: xm(3,3)
   complex (kind=CmplxKind) :: determ3_cmplx2
!
!  *******************************************************************
!
!  determ = the determinant of the 3x3 matrix xm.
!
!  *******************************************************************
   determ3_cmplx2 = xm(1,1)*(xm(2,2)*xm(3,3)-xm(2,3)*xm(3,2))      &
                   -xm(1,2)*(xm(2,1)*xm(3,3)-xm(2,3)*xm(3,1))      &
                   +xm(1,3)*(xm(2,1)*xm(3,2)-xm(2,2)*xm(3,1))
!
   end function determ3_cmplx2
!  ===================================================================
!
end module Matrix3dModule
