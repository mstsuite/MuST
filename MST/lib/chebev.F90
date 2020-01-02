!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function chebev(a,b,c,m,x)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : zero
   use MathParamModule, only : half
   use MathParamModule, only : two
   
   use ErrorHandlerModule 
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: m
   integer :: j
!
   real (kind=RealKind) :: chebev
   real (kind=RealKind), intent(in) :: a,b,x,c(m)
   real (kind=RealKind) :: d,dd,sv,y,y2
!
   if ((x-a)*(x-b).gt.zero) then
      call ErrorHandler('CHEBEV','x not in range in chebev',x)
   endif
!
   d=zero
   dd=zero
   y=(two*x-a-b)/(b-a)
   y2=two*y
   do j=m,2,-1
      sv=d
      d=y2*d-dd+c(j)
      dd=sv
   enddo
   chebev=y*d-dd+half*c(1)
   end function chebev
