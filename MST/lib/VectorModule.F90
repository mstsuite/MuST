module VectorModule
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
   use MathParamModule, only : ZERO
!
   implicit none
!
public :: getVecLength,    &
          getDotProduct,   &
          getCrossProduct
private
   integer (kind=IntKind) :: nstore = 0
!
   real (kind=RealKind), allocatable, target :: vstore(:)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVecLength(n,vec) result(r)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: vec(n)
   real (kind=RealKind) :: r
!
   if (n < 1) then
      call ErrorHandler('getVecLength','invalid dimension',n)
   else if (n == 1) then
      r = abs(vec(1))
   else if (n == 2) then
      r = sqrt(vec(1)*vec(1) + vec(2)*vec(2))
   else if (n == 3) then
      r = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
   else
      r = vec(1)*vec(1)
      do i = 2, n
         r = r + vec(i)*vec(i)
      enddo
      r = sqrt(r)
   endif
!
   end function getVecLength
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDotProduct(n,vec1,vec2,absolute) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: vec1(n), vec2(n)
   real (kind=RealKind) :: p
!
   logical, optional, intent(in) :: absolute
!
   if (n < 1) then
      call ErrorHandler('getDotProduct','invalid dimension',n)
   else if (n == 1) then
      p = vec1(1)*vec2(1)
   else if (n == 2) then
      p = vec1(1)*vec2(1) + vec1(2)*vec2(2)
   else if (n == 3) then
      p = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
   else
      p = vec1(1)*vec2(1)
      do i = 2, n
         p = p + vec1(i)*vec2(i)
      enddo
   endif
!
   if (present(absolute)) then
      if (absolute .and. p < ZERO) then
         p = -p
      endif
   endif
!
   end function getDotProduct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCrossProduct(v1,v2) result(v)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: v1(3), v2(3)
   real (kind=RealKind) :: v(3)
!
!  *******************************************************************
!  returns:
!           ->  ->   ->
!           v = v1 x v2
!
!  *******************************************************************
!
   v(1)=v1(2)*v2(3)-v1(3)*v2(2)
   v(2)=v1(3)*v2(1)-v1(1)*v2(3)
   v(3)=v1(1)*v2(2)-v1(2)*v2(1)
!
   end function getCrossProduct
!  ===================================================================
end module VectorModule
