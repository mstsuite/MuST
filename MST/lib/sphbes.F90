!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sphbes(n,x,sj,sy,sjp,syp)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : sqrt_pio2
   use MathParamModule, only : half
   use MathParamModule, only : two
  
   use ErrorHandlerModule
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   real (kind=RealKind), intent(in) :: x
   real (kind=RealKind), intent(out) :: sj,sjp,sy,syp
!
   real (kind=RealKind) :: factor,order,rj,rjp,ry,ryp
!
   if (n.lt.0) then
      call ErrorHandler('SPHBES','bad argument in sphbes',n)
   else if (x.le.0.) then
      call ErrorHandler('SPHBES','bad argument in sphbes',x)
   endif
!
   order=n+half
!  -------------------------------------------------------------------
   call bessjy(x,order,rj,ry,rjp,ryp)
!  -------------------------------------------------------------------
   factor=sqrt_pio2/sqrt(x)
   sj=factor*rj
   sy=factor*ry
   sjp=factor*rjp-sj/(two*x)
   syp=factor*ryp-sy/(two*x)
!
   end subroutine sphbes
