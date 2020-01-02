!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine lagpolfit_c( step, nfit, xn, m, yn, x, yfit )
!  ==================================================================
   use KindParamModule, only : IntKind, CmplxKind
   use MathParamModule, only : CONE
   implicit none
!
   integer (kind=IntKind), intent(in) :: step, nfit, m
!
   complex (kind=CmplxKind), intent(in) :: x, xn(nfit), yn(m)
!
   complex (kind=CmplxKind) :: yfit(m)
!
   integer (kind=IntKind) :: i, l
!
   complex (kind=CmplxKind) :: pm1
!
   do i = 1,m
      pm1 = CONE
      do l = 1,step-1
         pm1 = pm1*(x-xn(l))/(xn(step)-xn(l))
      enddo
      do l = step+1,nfit
         pm1 = pm1*(x-xn(l))/(xn(step)-xn(l))
      enddo
      yfit(i) = yfit(i) + yn(i)*pm1
   enddo
!
   end subroutine lagpolfit_c
!  ==================================================================
