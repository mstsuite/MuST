!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFermiDiracFunc(z,mu,kBT) result(fd)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, ONE, Ten2m6
!
   implicit none
!
   real (kind=RealKind), intent(in) :: mu
   real (kind=RealKind), intent(in) :: kBT
   real (kind=RealKind), parameter :: temp_tol = TEN2m6
!
   real (kind=RealKind), intent(in) :: z
   real (kind=RealKind) :: fd
!
   if (kBT < temp_tol) then
      if (z <= mu) then
         fd = ONE
      else
         fd = ZERO
      endif
   else
      fd = ONE/(ONE+exp((z-mu)/kBT))
   endif
!
   end function getFermiDiracFunc
