! ********************************************************************
! *  Given a function f(x) value at x0, its derivative at x0, and    *
! *  its second derivative at x0, find a cubic polynormial           *
! *       p(x) =  a +  0  + c*x^2 + d*x^3,         if l = 0;         *
! *            = (0 + b*x + c*x^2 + d*x^3)*x,      if l = 1;         *
! *            = (a + b*x + c*x^2 +   0  )*x^l,    otherwise.        *
! *  that smoothly extends f(x) to x < x0. A requirement is that:    *
! *       p'(x=0) = 0                                                *
! *  Given x0, f=f(x0), df=f'(x0), ddf=f"(x0), and l, this routine   *
! *  finds a, b, c, and d.                                           *
! ********************************************************************
subroutine findCubicFit(l,x0,f,df,ddf,a,b,c,d)
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, TWO, THREE, FOUR, SIX, CZERO
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use Matrix3dModule, only : invm3
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: l
   integer (kind=IntKind) :: ierr
!
   complex (kind=CmplxKind), intent(in) :: f, df, ddf
   complex (kind=CmplxKind), intent(out) :: a, b, c, d
!
   real (kind=RealKind), intent(in) :: x0
   real (kind=RealKind) :: x2l, x2lm1, x2lm2
   real (kind=RealKind) :: cmat(3,3),inv_cmat(3,3)
!
   if (l < 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('findCubicFit','l < 0',l)
!     ----------------------------------------------------------------
   else if (l == 0) then
      cmat(1,1)=ONE;      cmat(2,1)=ZERO;          cmat(3,1)=ZERO
      cmat(1,2)=x0*x0;    cmat(2,2)=TWO*x0;        cmat(3,2)=TWO
      cmat(1,3)=x0*x0*x0; cmat(2,3)=THREE*x0*x0;   cmat(3,3)=SIX*x0
   else if (l == 1) then
      cmat(1,1)=x0*x0;    cmat(2,1)=TWO*X0;        cmat(3,1)=TWO
      cmat(1,2)=x0*x0*x0; cmat(2,2)=THREE*x0*x0;   cmat(3,2)=SIX*x0
      cmat(1,3)=x0**4;    cmat(2,3)=FOUR*x0*x0*x0; cmat(3,3)=FOUR*THREE*x0*x0
   else
      if (l == 2) then
         x2l = x0*x0;     x2lm1 = x0;              x2lm2 = ONE
      else
         x2l = x0**l;     x2lm1 = x0**(l-1);       x2lm2 = x0**(l-2)
      endif
      cmat(1,1)=x2l;       cmat(2,1)=l*x2lm1;      cmat(3,1)=l*(l-1)*x2lm2
      cmat(1,2)=x2l*x0;    cmat(2,2)=(l+1)*x2l;    cmat(3,2)=(l+1)*l*x2lm1
      cmat(1,3)=x2l*x0*x0; cmat(2,3)=(l+2)*x2l*x0; cmat(3,3)=(l+2)*(l+1)*x2l
   endif
!
!  -------------------------------------------------------------------
   call invm3(cmat,inv_cmat,ierr)
!  -------------------------------------------------------------------
   if (ierr > 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('findCubicFit','can not find a suitable cubic fit')
!     ----------------------------------------------------------------
   endif
   if (l == 0) then
      a = inv_cmat(1,1)*f + inv_cmat(1,2)*df + inv_cmat(1,3)*ddf
      b = CZERO
      c = inv_cmat(2,1)*f + inv_cmat(2,2)*df + inv_cmat(2,3)*ddf
      d = inv_cmat(3,1)*f + inv_cmat(3,2)*df + inv_cmat(3,3)*ddf
   else if (l == 1) then
      a = CZERO
      b = inv_cmat(1,1)*f + inv_cmat(1,2)*df + inv_cmat(1,3)*ddf
      c = inv_cmat(2,1)*f + inv_cmat(2,2)*df + inv_cmat(2,3)*ddf
      d = inv_cmat(3,1)*f + inv_cmat(3,2)*df + inv_cmat(3,3)*ddf
   else
      a = inv_cmat(1,1)*f + inv_cmat(1,2)*df + inv_cmat(1,3)*ddf
      b = inv_cmat(2,1)*f + inv_cmat(2,2)*df + inv_cmat(2,3)*ddf
      c = inv_cmat(3,1)*f + inv_cmat(3,2)*df + inv_cmat(3,3)*ddf
      d = CZERO
   endif
!
!  *******************************************************************
!  check if there is a saddle point within [0,x0].
!  *******************************************************************
!  if (c > TEN2m8) then
!  endif
!
end subroutine findCubicFit
