!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function adjustEnergy_r(is,e) result(energy)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO
!
   use PotentialModule, only : getVdif
!
   use ScfDataModule, only : isInterstitialElectronPolarized
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind), pointer :: p_vdif(:)
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind) :: energy, vdif
!
   if (isInterstitialElectronPolarized()) then
      p_vdif => getVdif()
      vdif = p_vdif(1)
   else
      vdif = ZERO
   endif
!
   energy = e - (is-1)*vdif
!
   end function adjustEnergy_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function adjustEnergy_c(is,e) result(energy)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO
!
   use PotentialModule, only : getVdif
!
   use ScfDataModule, only : isInterstitialElectronPolarized
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind), pointer :: p_vdif(:)
   real (kind=RealKind) :: vdif
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: energy
!
   if (isInterstitialElectronPolarized()) then
      p_vdif => getVdif()
      vdif = p_vdif(1)
   else
      vdif = ZERO
   endif
!
   energy = e - (is-1)*vdif
!
   end function adjustEnergy_c
!  ===================================================================
