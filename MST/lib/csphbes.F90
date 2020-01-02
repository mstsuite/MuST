!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine csphbes(n,x,sj,sy,sjp,syp)
!  ===================================================================
   use KindparamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : sqrt_pio2
   use MathParamModule, only : half
   use MathParamModule, only : two
   
   use ErrorHandlerModule
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   complex (kind=CmplxKind), intent(in) :: x
   complex (kind=CmplxKind), intent(out) :: sj,sjp,sy,syp
!
   real (kind=RealKind) :: order
   complex (kind=CmplxKind) :: factor,rj,rjp,ry,ryp
!
   if (n.lt.0) then
      call ErrorHandler('SPHBES','bad argument in sphbes',n)
   endif
!
   order=n+half
!  -------------------------------------------------------------------
   call cbessjy(x,order,rj,ry,rjp,ryp)
!  -------------------------------------------------------------------
   factor=sqrt_pio2/sqrt(x)
   sj=factor*rj
   sy=factor*ry
   sjp=factor*rjp-sj/(two*x)
   syp=factor*ryp-sy/(two*x)
!
   end subroutine csphbes
