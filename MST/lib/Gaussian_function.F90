function Gaussian_function(sigma, mu, e) result(f)
   use KindParamModule, only : RealKind
   use MathParamModule, only : PI2, TWO, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   real (kind=RealKind), intent(in) :: sigma, mu, e
   real (kind=RealKind) :: f
   real (kind=RealKind), parameter :: delta = TEN2m8
!
   if (sigma < delta) then
      call ErrorHandler('Gaussian_function','sigma is too small',sigma)
   else
      f = exp(-(e-mu)*(e-mu)/(TWO*sigma*sigma))/(sigma*sqrt(PI2))
   endif
end function Gaussian_function
