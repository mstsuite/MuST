module SpectralFunctionModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, HALF, ONE, TWO, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initSpectralFunction,     &
          endSpectralFunction,      &
          computeSpectralFunction,  &
          getSpectralFunction
!
private
   integer (kind=IntKind) :: LocalNumSites
!
   complex (kind=CmplxKind) :: energy
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSpectralFunction()
!  ===================================================================
   implicit none
!
!
   end subroutine initSpectralFunction
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSpectralFunction()
!  ===================================================================
   implicit none
!
!
   end subroutine endSpectralFunction
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeMedium(site,e)
!  ===================================================================
   use ScfDataModule, only : isKKRCPA
!
   use CPAMediumModule, only : computeCPAMedium
!
   use AccessoryMatrixModule, only : computeAccessoryMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
!
   complex (kind=CmplxKind), intent(in) :: e
!
   if (.not.isKKRCPA()) then
      call ErrorHandler('computeMedium',                              &
                        'This functionality is only implemented for the KKR-CPA method')
   else if (site < 1 .or. site > LocalNumSites) then
      call ErrorHandler('computeMedium','site index is out of range',site)
   endif
!
   energy = e  
!
!  -------------------------------------------------------------------
   call computeCPAMedium(energy)
!  -------------------------------------------------------------------
   call computeAccessoryMatrix(energy)
!  -------------------------------------------------------------------
!
   end subroutine computeMedium
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSpectralFunction(site,kvec) result (Ake)
!  ===================================================================
   use AccessoryMatrixModule, only : getFcMatrix, getFccMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
!
   real (kind=RealKind), intent(in) :: kvec(3)
   real (kind=RealKind) :: Ake
!
   complex (kind=CmplxKind), pointer :: Fc(:), Fcc(:)
!
   if (site < 1 .or. site > LocalNumSites) then
      call ErrorHandler('getSpectralFunction','site index is out of range',site)
   endif
!
   end function getSpectralFunction
!  ===================================================================
end module SpectralFunctionModule
