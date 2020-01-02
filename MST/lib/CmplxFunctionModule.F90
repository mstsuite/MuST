module CmplxFunctionModule
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use KindParamModule, only : CmplxKind
!
public :: ComplexSign,    &
          ImaginaryPart,  &
          RealPart,       &
          ComplexConvert
!
   interface ComplexSign
      module procedure ComplexSign1, ComplexSign2
   end interface
!
private
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function ImaginaryPart(c) result (i)
!  ===================================================================
   use MathParamModule, only : sqrtm1
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: c
!
   real (kind=RealKind) :: i
!
   i=real(-sqrtm1*c, RealKind)
!
   end function ImaginaryPart
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function ComplexConvert(r,i) result (c)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind) :: c
!
   real (kind=RealKind), intent(in) :: r, i
!
   c = cmplx(r, i, CmplxKind)
!
   end function ComplexConvert
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function RealPart(c) result (r)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: c
!
   real (kind=RealKind) :: r
!
   r=real(c, RealKind)
!
   end function RealPart
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function ComplexSign1(c) result (s)
!  ===================================================================
   use MathParamModule, only : cone
   use MathParamModule, only : zero
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: c
!
   complex (kind=CmplxKind) :: s
!
   real (kind=RealKind) :: r
!
   r=abs(c)
   if (r.gt.zero) then
      s=c/r
   else
      s=cone
   endif
!
   end function ComplexSign1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function ComplexSign2(c1, c2) result (s)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: c1, c2
!
   complex (kind=CmplxKind) :: s
!
   s=abs(c1)*ComplexSign1(c2)
!
   end function ComplexSign2
!  ===================================================================
end module CmplxFunctionModule
