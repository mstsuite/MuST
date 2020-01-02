! ********************************************************************
! *  Given a function f(x) value at x0, its derivative at x0         *
! *  find a quadratic polynormial:                                   *
! *       p(x) =  a +  0  + c*x^2,                 if l = 0;         *
! *            = (0 + b*x + c*x^2)*x,              if l = 1;         *
! *            = (a + b*x +   0  )*x^l,            otherwise.        *
! *  that smoothly extends f(x) to x < x0. A requirement is that:    *
! *       p'(x=0) = 0                                                *
! *  Given x0, f=f(x0), df=f'(x0), and l, this routine finds a, b,   *
! *  and c.                                                          *
! ********************************************************************
subroutine findQuadraticFit(l,x0,f,df,a,b,c)
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, TWO, THREE, FOUR, SIX, CZERO, TEN2m6
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: l
!
   complex (kind=CmplxKind), intent(in) :: f, df
   complex (kind=CmplxKind), intent(out) :: a, b, c
!
   real (kind=RealKind), intent(in) :: x0
   real (kind=RealKind) :: x2l, x2lm1, det
   real (kind=RealKind) :: cmat(2,2),inv_cmat(2,2)
!
   if (l < 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('findQuadraticFit','l < 0',l)
!     ----------------------------------------------------------------
   else if (l == 0) then
      cmat(1,1)=ONE;      cmat(2,1)=ZERO
      cmat(1,2)=x0*x0;    cmat(2,2)=TWO*x0
   else if (l == 1) then
      cmat(1,1)=x0*x0;    cmat(2,1)=TWO*X0
      cmat(1,2)=x0*x0*x0; cmat(2,2)=THREE*x0*x0
   else
      if (l == 2) then
         x2l = x0*x0;     x2lm1 = x0
      else
         x2l = x0**l;     x2lm1 = x0**(l-1)
      endif
      cmat(1,1)=x2l;      cmat(2,1)=l*x2lm1
      cmat(1,2)=x2l*x0;   cmat(2,2)=(l+1)*x2l
   endif
!
!  -------------------------------------------------------------------
   det = cmat(1,1)*cmat(2,2)-cmat(1,2)*cmat(2,1)
   if (abs(det) < TEN2m6 .and. abs(det) < x2l*x2l) then
!     ----------------------------------------------------------------
      call ErrorHandler('findQuadraticFit','det = 0',det)
!     ----------------------------------------------------------------
   endif
   inv_cmat(1,1) = cmat(2,2)/det; inv_cmat(2,1) =-cmat(2,1)/det
   inv_cmat(1,2) =-cmat(1,2)/det; inv_cmat(2,2) = cmat(1,1)/det
!  -------------------------------------------------------------------
   if (l == 0) then
      a = inv_cmat(1,1)*f + inv_cmat(1,2)*df
      b = CZERO
      c = inv_cmat(2,1)*f + inv_cmat(2,2)*df
   else if (l == 1) then
      a = CZERO
      b = inv_cmat(1,1)*f + inv_cmat(1,2)*df
      c = inv_cmat(2,1)*f + inv_cmat(2,2)*df
   else
      a = inv_cmat(1,1)*f + inv_cmat(1,2)*df
      b = inv_cmat(2,1)*f + inv_cmat(2,2)*df
      c = CZERO
   endif
!
end subroutine findQuadraticFit
