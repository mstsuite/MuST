function Fermi_function(kT,Ef,e) result(f)
   use KindParamModule, only : RealKind
   use MathParamModule, only : ZERO, ONE, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   real (kind=RealKind), intent(in) :: kT, Ef, e
   real (kind=RealKind) :: f
   real (kind=RealKind), parameter :: delta = TEN2m6
!
   if (kT < ZERO) then
      call ErrorHandler('Fermi_function','kT < 0',kT)
   else if (kT <= delta) then
      if (e <= ef) then
         f = ONE
      else
         f = ZERO
      endif
   else
      f = ONE/(ONE+exp((e-ef)/kT))
   endif
end function Fermi_function
