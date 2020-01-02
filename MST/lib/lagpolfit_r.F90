!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine lagpolfit_r( step, nfit, xn, m, yn, x, yfit )
!  ==================================================================
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ONE
   implicit none
!
   integer (kind=IntKind), intent(in) :: step, nfit, m
!
   real (kind=RealKind), intent(in) :: x, xn(nfit), yn(m)
!
   real (kind=RealKind), intent(out) :: yfit(m)
!
   integer (kind=IntKind) :: i, l
!
   real (kind=RealKind) :: pm1
!
   do i = 1,m
      pm1 = ONE
      do l = 1,step-1
         pm1 = pm1*(x-xn(l))/(xn(step)-xn(l))
      enddo
      do l = step+1,nfit
         pm1 = pm1*(x-xn(l))/(xn(step)-xn(l))
      enddo
      yfit(i) = yfit(i) + yn(i)*pm1
   enddo
!
   end subroutine lagpolfit_r
!  ==================================================================