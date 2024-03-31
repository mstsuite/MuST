module DMFTModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, THREE, PI, PI2, PI4, &
                               CZERO, CONE, TWO, HALF, SQRTm1, Y0, TEN2m8
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
   use IntegerFactorsModule, only : lofk, lofj, mofj, m1m, mofk, jofk
   use PublicTypeDefinitionsModule, only : GridStruct, LloydStruct
   use TimerModule, only : getTime, getRoutineCallTiming
   use MPPModule, only : MyPE, syncAllPEs
!
public :: initDMFT,          &
          endDMFT,           &
          runDMFTsolver,     &
          getSelfEnergy
!
private
   integer (kind=IntKind) :: kmax_kkr
   integer (kind=IntKind) :: kmax_gf
   integer (kind=IntKind) :: NumLocalAtoms
   integer (kind=IntKind) :: NumEs
!
   complex (kind=CmplxKind), allocatable :: e_contour(:)
   complex (kind=CmplxKind), target, allocatable :: sigma(:,:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initDMFT(lmax_kkr,lmax_gf,na,ne,emesh)
!  ===================================================================
   use localGFModule, only : initLocalGF
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax_kkr, lmax_gf, na, ne
!
   complex (kind=CMplxKind) :: emesh(ne)
!
   kmax_kkr = (lmax_kkr+1)**2
   kmax_gf = (lmax_gf+1)*2
   NumLocalAtoms = na
   NumEs = ne
!
   allocate(sigma(kmax_gf,na), e_contour(ne))
!
   e_contour = emesh
!
   call initLocalGF(lmax_kkr,lmax_gf,na,ne)
!
   end subroutine initDMFT
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endDMFT()
!  ===================================================================
   use localGFModule, only : endLocalGF
!
   implicit none
!
   deallocate(sigma, e_contour)
!
   call endLocalGF()
!
   end subroutine endDMFT
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine runDMFTsolver()
!  ===================================================================
   use LocalGFModule, only : computeLocalGF, getLocalGF
!
   implicit none
!
   integer (kind=IntKind) :; ie, id
!
   complex (kind=CmplxKind), pointer :: gf(:)
!
   do ie = 1, NumEs
      call computeLocalGF(ie,e_contour(ie))
   enddo
!
   do id = 1, NumLocalAtoms
      gf => getLocalGF(id)
   enddo
!
   end subroutine runDMFTsolver
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSelfEnergy(id) result (s)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   complex (kind=CmplxKind), pointer :: s(:)
!
   s => sigma(:,id)
!
   end function getSelfEnergy
!  ===================================================================
end module DMFTModule
