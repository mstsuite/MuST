!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!!    This module contains calculates and stores the Fc and Fcc matrices,
!!!    to be used by the spectral function calculations.
!!!===================================================================
module AccessoryMatrixModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, HALF, ONE, TWO, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initAccessoryMatrix,     &
          endAccessoryMatrix,      &
          computeAccessoryMatrix,  &
          getFcMatrix,             &
          getFccMatrix
!
private
   integer (kind=IntKind) :: LocalNumSites
   integer (kind=IntKind) :: kmax_max
   integer (kind=IntKind) :: ncs_max
!
   type AccMatStruct
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind) :: kmax_kkr
      complex (kind=CmplxKind), allocatable :: Fc(:,:)
      complex (kind=CmplxKind), allocatable :: Fcc(:,:)
      complex (kind=CmplxKind), allocatable :: Fab(:,:,:,:)
      complex (kind=CmplxKind), allocatable :: Da(:,:,:)
      complex (kind=CmplxKind), allocatable :: DaT(:,:,:)
   end type AccMatStruct
!
   type (AccMatStruct), allocatable :: SiteAccMat(:)
!
   complex (kind=CmplxKind), allocatable, target :: PhiMat(:,:,:)
   complex (kind=CmplxKind), allocatable, target :: sine_inv(:,:)
   complex (kind=CmplxKind), allocatable, target :: wspace(:)
!
   logical :: Initialized = .false.
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initAccessoryMatrix()
!  ===================================================================
   use MediumHostModule, only : getLocalNumSites, getNumSpecies
   use MediumHostModule, only : getLmaxKKR, getGlobalSiteIndex
!
   implicit none
!
   integer (kind=IntKind) :: i, ig, n, kmax
!
   LocalNumSites = getLocalNumSites()
!
   if (LocalNumSites < 1) then
      call ErrorHandler('initAccessoryMatrix','LocalNumSites < 1',LocalNumSites)
   else if(kmax_kkr < 1) then
      call ErrorHandler('initAccessoryMatrix','kmax_kkr < 1',kmax_kkr)
   endif
!
   allocate(SiteAccMat(LocalNumSites))
!
   kmax_max = 0
   ncs_max = 0
   do i = 1, LocalNumSites
      ig = getGlobalSiteIndex(i)
      n = getNumSpecies(ig)
      ncs_max = max(n,ncs_max)
      kmax = (getLmaxKKR(ig)+1)**2
      kmax_max = max(kmax,kmax_max)
      SiteAccMat(i)%NumSpecies = n
      SiteAccMat(i)%kmax_kkr = kmax
      allocate(SiteAccMat(i)%Fc(kmax,kmax))
      allocate(SiteAccMat(i)%Fcc(kmax,kmax))
      allocate(SiteAccMat(i)%Fab(kmax,kmax,n,n))
      allocate(SiteAccMat(i)%Da(kmax,kmax,n))
      allocate(SiteAccMat(i)%DaT(kmax,kmax,n))
   enddo
!
   allocate(PhiMat(kmax_max*kmax_max,ncs_max,ncs_max))
   allocate(sine_inv(kmax_max*kmax_max,ncs_max))
   allocate(wspace(kmax_max*kmax_max))
!
   Initialized = .true.
!
   end subroutine initAccessoryMatrix
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endAccessoryMatrix()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i
!
   do i = 1, LocalNumSites
      deallocate(SiteAccMat(i)%Fc)
      deallocate(SiteAccMat(i)%Fcc)
      deallocate(SiteAccMat(i)%Fab)
      deallocate(SiteAccMat(i)%Da)
      deallocate(SiteAccMat(i)%DaT)
   enddo
!
   deallocate(SiteAccMat)
!
   deallocate(PhiMat)
   deallocate(sine_inv)
   deallocate(wspace)
!
   Initialized = .false.
!
   end subroutine endAccessoryMatrix
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeAccessoryMatrix(e)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IntegrationModule, only : calIntegration
!
   use SSSolverModule, only : getRegSolution, getSolutionRmeshSize
   use SSSolverModule, only : getSineMatrix
!
   use CPAMediumModule, only : getCPAMatrix
!
   use PotentialTypeModule, only : isASAPotential, isFullPotential
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use RadialGridModule, only : getGrid
!
   implicit none
!
   integer (kind=IntKind) :: ic, jc, n, nr
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: tcpa(:,:), tau_cpa(:,:), tau_a(:,:)
   complex (kind=CmplxKind), pointer :: Fc(:,:), Fcc(:,:), Da(:,:), DaT(:,:), Fab(:,:,:)
   complex (kind=CmplxKind), pointer :: phil_a(:,:,:), phil_b(:,:,:)
   complex (kind=CmplxKind), pointer :: phil_a_L(:,:), phil_b_L(:,:)
   complex (kind=CmplxKind), pointer :: sine_mat(:,:), sinv(:,:), Fhi_ab(:,:)
   complex (kind=CmplxKind), allocatable :: ffr(:)
!
   type (GridStruct), pointer :: Grid
!
   if (.not.Initialized) then
      call ErrorHandler('computeAccessoryMatrix','Need to call initAccessoryMatrix first')
   endif
!
   nr = 0
   do n = 1, LocalNumSites
      nr = max(nr,getSolutionRmeshSize(n))
   enddo
   allocate(ffr(0:nr), sqrt_r(0:nr))
!
   do n = 1, LocalNumSites
      Grid => getGrid(n)
      nr = getSolutionRmeshSize(n)
      sqrt_r(0) = ZERO
      do ir = 1, nr
         sqrt_r(ir) = sqrt(Grid%r_mesh(ir))
      enddo
!     ----------------------------------------------------------------
      tau_cpa => getCPAMatrix('Tau',site=n,dsize=kmax)
!     ----------------------------------------------------------------
      if (kmax /= SiteAccMat(n)%kmax_kkr) then
!        -------------------------------------------------------------
         call ErrorHandler('computeAccessoryMatrix','Returned kmax <> kmax_kkr', &
                           kmax,SiteAccMat(n)%kmax_kkr)
!        -------------------------------------------------------------
      endif
      Phi_ab => aliasArray2_c(wspace(kmax,kmax) 
      do ic = 1, SiteAccMat(n)%NumSpecies
         sine_mat => getSineMatrix(site=n,atom=ic)
         sinv => aliasArray2_c(sine_inv(:,ic),kmax,kmax)
         sinv = sine_mat
!        -------------------------------------------------------------
         call MtxInv_LU(sinv,kmax)
!        -------------------------------------------------------------
      enddo
      Fc => SiteAccMat(n)%Fc
      Fcc => SiteAccMat(n)%Fcc
!
      if (isFullPotential()) then
         nr = getSolutionRmeshSize(site=n)
      else if (isASAPotential()) then
         nr = getSolutionRmeshSize(site=n)
         if (nr < Grid%jws) then
!           ----------------------------------------------------------
            call ErrorHandler('computeAccessoryMatrix','nr < jws',    &
                              nr,Grid%jws)
!           ----------------------------------------------------------
         endif
         nr = Grid%jws
      else
         nr = getSolutionRmeshSize(site=n)
         if (nr < Grid%jmt) then
!           ----------------------------------------------------------
            call ErrorHandler('computeAccessoryMatrix','nr < jmt',    &
                              nr,Grid%jmt)
!           ----------------------------------------------------------
         endif
         nr = Grid%jmt
      endif
      do ic = 1, SiteAccMat(n)%NumSpecies
!        -------------------------------------------------------------
         tcpa => getCPAMatrix('Tcpa',site=n,atom=ic)
         tau_a => getCPAMatrix('Tau',site=n,atom=ic)
!        -------------------------------------------------------------
         phil_a => getRegSolution(site=n,atom=ic)
         Da => SiteAccMat(n)%Da(:,:,ic)
         DaT => SiteAccMat(n)%DaT(:,:,ic)
         Fab => SiteAccMat(n)%Fab(:,:,:,ic)
         do jc = 1, SiteAccMat(n)%NumSpecies
            phil_b => getRegSolution(site=n,atom=jc)
! Needs to take care compelx conjugate........................
            if (isFullPotential()) then
               do kl = 1, kmax
                  phil_a_L => phil_a(:,:,kl)
                  do klp = 1, kmax
                     phil_b_L => phil_b(:,:,klp)
                     Fhi_ab(klp,kl) = getVolumeIntegration(n,nr,Grid%r_mesh, &
                                                           kmax,1,phil_b_L,  &
                                                           kmax,1,phil_a_L,v_mr)
                  enddo
               enddo
            else
               Fhi_ab = CZERO
               do kl = 1, kmax
                  phil_a_L => phil_a(:,:,kl)
                  phil_b_L => phil_b(:,:,klp)
                  ffr(0) = CZERO
                  do ir = 1, nr
                     ffr(ir) = phil_a_L(ir,kl)*phil_b_L(ir,kl) ! which includes r^2
                  enddo
!                 ----------------------------------------------------
                  call calIntegration(nr+1,sqrt_r(0:nr),ffr(0:nr),Fhi_ab(kl,kl),1)
!                 ----------------------------------------------------
                  Fhi_ab(kl,kl) = TWO*Fhi_ab(kl,kl) ! Factor TWO is needed because of 
                                                    ! integration over sqrt(r)
               enddo
            endif
Da = 
Fab = 
         enddo
      enddo
   enddo
!
   deallocate(ffr, sqrt_r)
!
   end subroutine computeAccessoryMatrix
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFccMatrix(i,kmax) result (fm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind), intent(out), optional :: kmax
!
   complex (kind=CmplxKind), pointer :: fm(:,:)
!
   if (i < 1 .or. i > LocalNumSites) then
      call ErrorHandler('getFccMatrix','i is out of range',i)
   else if (.not.Initialized) then
      call ErrorHandler('getFccMatrix','Need to call initAccessoryMatrix first')
   endif
!
   fm => SiteAccMat(i)%Fcc(:,:)
!
   if (present(kmax)) then
      kmax = SiteAccMat(i)%kmax_kkr
   endif
!
   end function getFccMatrix
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFcMatrix(i,kmax) result (fm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind), intent(out), optional :: kmax
!
   complex (kind=CmplxKind), pointer :: fm(:,:)
!
   if (i < 1 .or. i > LocalNumSites) then
      call ErrorHandler('getFcMatrix','i is out of range',i)
   else if (.not.Initialized) then
      call ErrorHandler('getFcMatrix','Need to call initAccessoryMatrix first')
   endif
!
   fm => SiteAccMat(i)%Fc(:,:)
!
   if (present(kmax)) then
      kmax = SiteAccMat(i)%kmax_kkr
   endif
!
   end function getFcMatrix
!  ===================================================================
end module AccessoryMatrixModule
