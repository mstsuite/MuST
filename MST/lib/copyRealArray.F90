!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine copyRealArray(a,b,n)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
!
   real (kind=RealKind), intent(in) :: a(n)
   real (kind=RealKind), intent(out) :: b(n)
!
   b(1:n) = a(1:n)
!
   end subroutine copyRealArray
!  ===================================================================