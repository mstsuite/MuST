module LegendreModule
!
public :: legendre
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine legendre(lmax,x,plm)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
   use MathParamModule, only : zero
   use MathParamModule, only : one
   use MathParamModule, only : two
   use MathParamModule, only : ten2m8
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(out) :: plm((lmax+1)*(lmax+2)/2)
   real (kind=RealKind), intent(in) :: x
!
   real (kind=RealKind) :: pmm
   real (kind=RealKind) :: somx2
   real (kind=RealKind) :: fact
!
!  ===================================================================
!  Calclates Associated Legendre function, p(l,m), up to lmax.
!  Based on the formulae given in "Numerical Recipes" pages 180-183
!  (Equations 6.6.7, 6.6.8 and 6.6.9)
!  W. H. Press, B. P. Flannery, S A Teukolsky and W. T. Vetterling.
!  Cambridge Univ Press 1986.
!
!  N.B. The definition of p(l,m) has been modified.
!  p(l,m) of this code = [(-1)**m]*p(l,m) of "Numerical Recipes".
!  ===================================================================
   if(lmax.lt.0) then
      call ErrorHandler('Legendre','bad lmax: lmax < 0',lmax)
   else if(abs(x).gt.one) then
      call ErrorHandler('Legendre','bad arguments: abs(x) > 1',x)
   endif
!
   plm = zero
   jmax=(lmax+1)*(lmax+2)/2
   if((one-abs(x)).le.ten2m8) then
      plm(1:jmax)=zero
      if(x.lt.zero) then
         do l=0,lmax
            i=(l+1)*(l+2)/2-l
            plm(i)=one-two*mod(l,2)
         enddo
      else
         do l=0,lmax
            i=(l+1)*(l+2)/2-l
            plm(i)=one
         enddo
      endif
      return
   endif
!
!  ===================================================================
!  begin calculation of p(l,m)'s...................................
!  ===================================================================
   if(lmax.eq.0) then
!     ================================================================
!     special case lmax=0..........................................
!     ================================================================
      plm(1)=one
   else if(lmax.eq.1) then
!     ================================================================
!     special case lmax=1..........................................
!     ================================================================
      plm(1)=one
      plm(2)=x
      plm(3)=sqrt((one-x)*(one+x))
   else
      plm(1)=one
      plm(2)=x
      somx2=sqrt((one-x)*(one+x))
      do m=1,lmax-1
!        =============================================================
!                                 m       m
!        calculate the first two P   and P
!                                 m       m+1
!        =============================================================
         pmm=one
         fact=one
         do i=1,m
            pmm=pmm*fact*somx2
            fact=fact+two
         enddo
         plm(m*(m+1)/2+m+1)=pmm
         plm((m+1)*(m+2)/2+m+1)=x*(2*m+1)*pmm
      enddo
      pmm=one
      fact=one
      do i=1,lmax
         pmm=pmm*fact*somx2
         fact=fact+two
      enddo
      plm(lmax*(lmax+1)/2+lmax+1)=pmm
!     ================================================================
!                         m        m
!     calculate the rest P     to P
!                         m+2      lmax
!     ================================================================
      do m=0,lmax
         do l=m+2,lmax
            plm(l*(l+1)/2+m+1)=( x*(2*l-1)*plm((l-1)*l/2+m+1) -        &
                                 (l+m-1)*plm((l-2)*(l-1)/2+m+1) )/dble(l-m)
         enddo
      enddo
   endif
!
   end subroutine legendre
!  ===================================================================
end module LegendreModule
