function ErrorStep_function(width, R0, r) result(f)
   use KindParamModule, only : RealKind
   use MathParamModule, only : ZERO, ONE, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   real (kind=RealKind), intent(in) :: width, R0, r
   real (kind=RealKind) :: f, x, f0
!
   interface
      function erfc(x) result(f)
         real*8, intent(in) :: x
         real*8 :: f
      end function erfc
   end interface
!
   if (width < TEN2m8 .or. R0 < TEN2m8 .or. r < ZERO) then
      call ErrorHandler('ErrorStep_function','parameter is too small')
   else
      f0 = erfc(-R0/width)
      x = (r-R0)/width
      f = erfc(x)/f0
   endif
end function ErrorStep_function
