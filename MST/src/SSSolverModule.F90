module SSSolverModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : TEN2m12, TEN2m10, TEN2m8, TEN2m6, TEN2m4, TEN2m11
   use MathParamModule, only : ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, NINE, TEN
   use MathParamModule, only : PI, PI2, SQRT_PI, HALF, TEN2m5, TEN2m16
   use MathParamModule, only : CZERO, CONE, SQRTm1, Y0
!
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
!
   use IntegerFactorsModule, only : lofk, mofk, jofk, lofj, mofj, kofj, m1m
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
public :: initSSSolver,          &
          endSSSolver,           &
          solveSingleScattering, &
          computePhaseShift,     &
          computeDOS,            &
          computePDOS,           &
          computeGreenFunction,  &
          getScatteringMatrix,   &
          getSineMatrix,         &
          getCosineMatrix,       &
          getJostMatrix,         &
          getJostInvMatrix,      &
          getOmegaHatMatrix,     &
          getOmegaHatInvMatrix,  &
          getTMatrix,            &
          getSMatrix,            &
          getRegSolution,        &
          getSolutionFlags,      &
          getSolutionRmeshSize,  & 
          getDOS,                &
          getPDOS,               &
          getOutsideDOS,         &
          getCellDOS,            &
          getMTSphereDOS,        &
          getCellPDOS,           &
          getMTSpherePDOS,       &
          getPhaseShift,         &
          getGreenFunction,      &
          getDOSDerivative,      &
          getRegSolutionDerivative, &
          getGreenFunctionDerivative
!
!tmat_global          convertTmatToGlobalFrame
!
          interface computeGreenFunction
             module procedure computeGF0, computeGF1
          end interface computeGreenFunction
!
private
!
   integer (kind=IntKind) :: LocalNumSites
   integer (kind=IntKind) :: NumSpins
!
   type SolutionStruct
      real (kind=RealKind), pointer :: phase_shift(:)
!
      complex (kind=CmplxKind), pointer :: sin_mat(:,:)
      complex (kind=CmplxKind), pointer :: cos_mat(:,:)
      complex (kind=CmplxKind), pointer :: Omega_mat(:,:)
      complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)
      complex (kind=CmplxKind), pointer :: OmegaHatInv_mat(:,:)
      complex (kind=CmplxKind), pointer :: S_mat(:,:)
      complex (kind=CmplxKind), pointer :: t_mat(:,:)
      complex (kind=CmplxKind), pointer :: t_mat_inv(:,:)
      complex (kind=CmplxKind), pointer :: jost_mat(:,:)
      complex (kind=CmplxKind), pointer :: jinv_mat(:,:)
      complex (kind=CmplxKind), pointer :: reg_sol(:,:,:)  ! Phi_{Lp,L}(r)*r
      complex (kind=CmplxKind), pointer :: irr_sol(:,:,:)  ! H_{Lp,L}(r)*r
      complex (kind=CmplxKind), pointer :: reg_dsol(:,:,:) ! r*d{Phi_{Lp,L}(r)}/dr
      complex (kind=CmplxKind), pointer :: irr_dsol(:,:,:)
      complex (kind=CmplxKind), pointer :: green(:,:)
      complex (kind=CmplxKind), pointer :: der_green(:,:)
      complex (kind=CmplxKind), pointer :: dos(:,:)
      complex (kind=CmplxKind), pointer :: der_dos(:,:)
      complex (kind=CmplxKind), pointer :: pdos(:,:,:)
   end type SolutionStruct
!
   type ScatterStruct
      integer (kind=IntKind) :: numrs
      integer (kind=IntKind) :: numrs_cs
      integer (kind=IntKind) :: numrs_mt
      integer (kind=IntKind) :: numrs_trunc
      integer (kind=IntKind) :: lmax_kkr
      integer (kind=IntKind) :: kmax_kkr
      integer (kind=IntKind) :: lmax_phi
      integer (kind=IntKind) :: kmax_phi
      integer (kind=IntKind) :: kmax_int
      integer (kind=IntKind) :: lmax_pot
      integer (kind=IntKind) :: jmax_pot
      integer (kind=IntKind) :: lmax_pot_trunc
      integer (kind=IntKind) :: jmax_pot_trunc
      integer (kind=IntKind) :: lmax_step
      integer (kind=IntKind) :: jmax_step
      integer (kind=IntKind) :: lmax_green
      integer (kind=IntKind) :: kmax_green
      integer (kind=IntKind) :: NumSpecies
!
      integer (kind=IntKind), pointer :: AtomicNumber(:)
      integer (kind=IntKind), pointer :: sol_flags(:,:)
      integer (kind=IntKind), pointer :: lm_index_sph(:,:)
!
      complex (kind=CmplxKind), pointer :: step_rad(:,:)
!tmat_global      complex (kind=CmplxKind), pointer :: tmat_global(:,:,:)
!
      type (SolutionStruct), pointer :: Solutions(:,:)
   end type ScatterStruct
!
   type (ScatterStruct), allocatable :: Scatter(:)
!
   integer (kind=IntKind), allocatable :: npi(:,:,:,:)
   real (kind=RealKind), allocatable :: phase_ref(:,:,:,:)
!
   type (GridStruct), pointer :: Grid
!
   character (len=50) :: stop_routine
   character (len=1) :: IrrSolType
!
   logical :: FullSolver
   logical :: isSphericalSolverOn
   logical :: PrintPotT_done
   logical :: isIrrSolOn
   logical :: isDiagSolv
   logical :: isScaleOn = .false.
   logical :: isS2C2invCS = .true.
   logical :: isIsoparam = .false.
   logical :: CheckWronskian = .false.
!
   integer :: SSSMethod
!
   integer (kind=IntKind) :: LocalIndex, LocalSpin
!
   integer (kind=IntKind) :: lmax_phi, kmax_phi, kmax_int
   integer (kind=IntKind) :: lmax_kkr, kmax_kkr
   integer (kind=IntKind) :: lmax_pot, jmax_pot
   integer (kind=IntKind) :: lmax_step, jmax_step
   integer (kind=IntKind) :: lmax_green, kmax_green
   integer (kind=IntKind) :: lmax_trunc, jmax_trunc
   integer (kind=IntKind) :: jmax_pot_solver
   integer (kind=IntKind) :: jmax_pot_solver_default
   integer (kind=IntKind) :: jmax_potmp
   integer (kind=IntKind) :: kmax_max_kkr, kmax_max_phi, lmax_max_kkr, kmax_max_int
   integer (kind=IntKind) :: MaxSpecies
!
   integer (kind=IntKind), allocatable :: print_instruction(:)
!
   integer (kind=IntKind) :: iend    ! no. of r-mesh needed
   integer (kind=IntKind) :: numrs_cs
   integer (kind=IntKind) :: numrs_mt
   integer (kind=IntKind) :: numrs_trunc
   integer (kind=IntKind) :: spin_index
!
   integer (kind=1) :: Relativity
   integer (kind=1), parameter :: NonRelativistic = 0_1
   integer (kind=1), parameter :: ScalarRelativistic = 1_1
   integer (kind=1), parameter :: FullyRelativistic = 2_1
!
   integer (kind=1), parameter :: RegularSolution = 0_1
   integer (kind=1), parameter :: IrregularSolution = 1_1
   integer (kind=1), parameter :: NonDerivative = 0_1
   integer (kind=1), parameter :: Derivative = 1_1
!
   real (kind=RealKind) :: c2inv
   real (kind=RealKind) :: etol = FIVE*TEN2m16
!
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind) :: PotShift
!
   complex (kind=CmplxKind), allocatable :: c_sph_mt(:)
   complex (kind=CmplxKind), allocatable :: s_sph_mt(:)
   complex (kind=CmplxKind), allocatable :: c_mtrx(:)
   complex (kind=CmplxKind), allocatable :: s_mtrx(:)
!
!  v0r = r*(spherical part of potential), i.e., v (r)*r
!                                                0
   complex (kind=CmplxKind), allocatable :: v0r(:)
   complex (kind=CmplxKind), allocatable :: v0r_save(:)
!
   complex (kind=CmplxKind), allocatable, target :: cm0(:)
   complex (kind=CmplxKind), allocatable, target :: bjl(:,:)
   complex (kind=CmplxKind), allocatable, target :: bnl(:,:)
   complex (kind=CmplxKind), allocatable, target :: bhl(:,:)
   complex (kind=CmplxKind), allocatable, target :: dbjl(:,:)
   complex (kind=CmplxKind), allocatable, target :: dbnl(:,:)
   complex (kind=CmplxKind), allocatable, target :: dbhl(:,:)
!
   complex (kind=CmplxKind), allocatable, target :: pl_regwrk(:,:)
   complex (kind=CmplxKind), allocatable, target :: ql_regwrk(:,:)
   complex (kind=CmplxKind), allocatable, target :: pl_irrwrk(:,:)
   complex (kind=CmplxKind), allocatable, target :: ql_irrwrk(:,:)
!
   complex (kind=CmplxKind), allocatable :: nlrmt(:), dnlrmt(:)
!
   logical :: Initialized = .false.
   logical :: IrrSpaceAllocated = .false.
   logical :: isSphPotZeroOutsideRmt
   logical :: isPotZeroOutsideRmt
!
   integer (kind=IntKind) :: iend_max
   integer (kind=IntKind) :: NumDiffSurfR
   integer (kind=IntKind) :: lmax_sol_cutoff
   integer (kind=IntKind), allocatable :: JHunt(:)
   integer (kind=IntKind), allocatable :: IPVT(:)
!
   integer (kind=IntKind) :: NumSurfPoints
!
   integer (kind=IntKind), pointer :: GaussR2SurfR(:)
   integer (kind=IntKind), pointer :: nj3(:,:)
   integer (kind=IntKind), pointer :: kj3(:,:,:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
   real (kind=RealKind), pointer :: DiffSurfR(:)
   real (kind=RealKind), pointer :: n_dot_r(:)
   real (kind=RealKind), pointer :: weight(:)
   real (kind=RealKind), pointer :: r_surf(:)
!
   complex (kind=CmplxKind), pointer :: smpicm(:,:)
   complex (kind=CmplxKind), pointer :: pl_reg(:)
   complex (kind=CmplxKind), pointer :: ql_reg(:)
   complex (kind=CmplxKind), pointer :: pl_irr(:)
   complex (kind=CmplxKind), pointer :: ql_irr(:)
   complex (kind=CmplxKind), pointer :: qlhat_reg(:,:)
   complex (kind=CmplxKind), pointer :: qlhat_irr(:,:)
   complex (kind=CmplxKind), pointer :: plhat_reg(:,:)
   complex (kind=CmplxKind), pointer :: plhat_irr(:,:)
   complex (kind=CmplxKind), pointer :: ylm(:,:)
   complex (kind=CmplxKind), pointer :: nDotGradY(:,:)
   complex (kind=CmplxKind), pointer :: pot_jl(:,:)
   complex (kind=CmplxKind), pointer :: pot_trunc_jl(:,:)
!
   complex (kind=CmplxKind), allocatable :: sjl(:,:)
   complex (kind=CmplxKind), allocatable :: snl(:,:)
   complex (kind=CmplxKind), allocatable :: dsjl(:,:)
   complex (kind=CmplxKind), allocatable :: dsnl(:,:)
!
   integer (kind=IntKind), allocatable, target   :: flags_jl(:)
   integer (kind=IntKind), allocatable, target   :: flags_trunc_jl(:)
!
   real (kind=RealKind), allocatable, target :: wks_PS(:)
!
   complex (kind=CmplxKind), allocatable, target :: wks_sinmat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_cosmat(:)
!
   complex (kind=CmplxKind), allocatable, target :: wks_tmat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_tmat_inv(:)
!tmat_global   complex (kind=CmplxKind), allocatable, target :: wks_tmatg(:)
   complex (kind=CmplxKind), allocatable, target :: wks_jostmat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_jinvmat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_Omega(:)
   complex (kind=CmplxKind), allocatable, target :: wks_OmegaHat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_OmegaHatInv(:)
   complex (kind=CmplxKind), allocatable, target :: wks_S(:)
!
   complex (kind=CmplxKind), allocatable, target :: wks_plhat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_qlhat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_regsol(:)
   complex (kind=CmplxKind), allocatable, target :: wks_irrsol(:)
   complex (kind=CmplxKind), allocatable, target :: wks_dregsol(:)
   complex (kind=CmplxKind), allocatable, target :: wks_dirrsol(:)
!
   complex (kind=CmplxKind), allocatable, target :: wks_green(:)
   complex (kind=CmplxKind), allocatable, target :: wks_dgreen(:)
   complex (kind=CmplxKind), allocatable, target :: wks_dos(:)
   complex (kind=CmplxKind), allocatable, target :: wks_ddos(:)
   complex (kind=CmplxKind), allocatable, target :: wks_pdos(:)
!
   complex (kind=CmplxKind), allocatable, target :: wks_mtx1(:)
   complex (kind=CmplxKind), allocatable, target :: wks_mtx2(:)
!
   complex (kind=CmplxKind), allocatable, target :: wksl_invscm(:)
   complex (kind=CmplxKind), allocatable, target :: wksl_smpicm(:)
   complex (kind=CmplxKind), allocatable, target :: wksl_invscmt(:)
   complex (kind=CmplxKind), allocatable, target :: wksl_smpicmt(:)
   complex (kind=CmplxKind), allocatable, target :: TmpSpace(:)
   complex (kind=CmplxKind), allocatable, target :: wksl_tmpirr(:)
   complex (kind=CmplxKind), allocatable, target :: wksl_plsurf(:)
   complex (kind=CmplxKind), allocatable, target :: wksl_qlsurf(:)
   complex (kind=CmplxKind), allocatable, target :: wks_step(:)
!
   complex (kind=CmplxKind), allocatable :: al(:),bl(:)
   complex (kind=CmplxKind), allocatable :: scp(:),hp(:)
   complex (kind=CmplxKind), allocatable :: bjtmp(:), bntmp(:)
!
   complex (kind=CmplxKind), allocatable :: factor_L(:)
!
   integer (kind=IntKind), allocatable :: WronskianIndx(:,:)
!
   real (kind=RealKind), allocatable, target :: wks_int_rtmp(:)
   real (kind=RealKind), allocatable :: space_integrated_dos_cell(:,:,:)
   real (kind=RealKind), allocatable :: space_integrated_dos_mt(:,:,:)
   real (kind=RealKind), allocatable, target :: space_integrated_pdos_cell(:,:,:,:)
   real (kind=RealKind), allocatable, target :: space_integrated_pdos_mt(:,:,:,:)
!
   complex (kind=CmplxKind), allocatable, target :: wks_int_lsph(:)
   complex (kind=CmplxKind), allocatable, target :: wks_int_rsph(:)
   complex (kind=CmplxKind), allocatable, target :: wks_int_res(:)
   complex (kind=CmplxKind), allocatable, target :: wks_int_fg(:)
   complex (kind=CmplxKind), allocatable, target :: wks_int_bjl(:)
   complex (kind=CmplxKind), allocatable, target :: wks_int_bnl(:)
   complex (kind=CmplxKind), allocatable, target :: wks_scvphi(:,:)
   complex (kind=CmplxKind), allocatable, target :: cpotstep(:,:)
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, kGID
!
   logical :: isDosSymmOn = .false.
   logical :: rad_deriv = .false.
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSSSolver( na, getNumSpecies, getAtomicNumber,       &
                            lkkr, lphi, lpot, lstep, lgreen, nsp, nsc, rel,  &
                            istop, iprint, derivative )
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
!
   use IntegerFactorsModule, only : initIntegerFactors, pushIntegerFactorsToAccel
!
   use PhysParamModule, only : LightSpeed
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   use RadialGridModule, only : getGrid, pushRadialGridToAccel
!
   use PotentialModule, only : getPotLmax, getTruncPotLmax
!
   use LdaCorrectionModule, only : checkLdaCorrection
!
   use ScfDataModule, only : getSingleSiteSolverType, getSingleSiteSolverMethod
   use ScfDataModule, only : getLmaxPotSolver
   use ScfDataModule, only : getLmaxSolution
   use ScfDataModule, only : isChargeSymm
!
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
   use GauntFactorsModule, only : pushGauntFactorsToAccel
!
   use StepFunctionModule, only : getNumGaussRs
   use StepFunctionModule, only : getRadialStepFunction
!
   use SurfElementsModule, only : getNumDiffSurfR
!
   implicit   none
!
   character (len=*), intent(in) :: istop
!
   logical, intent(in), optional :: derivative
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: lkkr(na)
   integer (kind=IntKind), intent(in) :: lphi(na)
   integer (kind=IntKind), intent(in) :: lpot(na)
   integer (kind=IntKind), intent(in) :: lstep(na)
   integer (kind=IntKind), intent(in) :: lgreen(na)
   integer (kind=IntKind), intent(in) :: nsp, nsc
   integer (kind=IntKind), intent(in) :: rel
   integer (kind=IntKind), intent(in) :: iprint(na)
   integer (kind=IntKind) :: is, ia, ic, kl, l, m, lmax, kmax, jmax, lmax_step, n, npout, kpm, kpsm, i
   integer (kind=IntKind) :: index_end, lmax_trunc, lmax_potmp
   integer (kind=IntKind) :: jmax_max, jmax_trunc_max, lmax_pot_in, jmax_pot_in, jmax_step, lmax_max_green
   integer (kind=IntKind) :: solv_type, ng, ng_max
   integer (kind=IntKind) :: max_NDiffSurfR, NDiffSurfR, INFO, js
!
   integer (kind=IntKind) :: sz_1, sz_2, skk_1, skk_2, spk_1, spk_2, spp_1, spp_2
!
   integer (kind=IntKind) :: sz_ind_intint
   integer (kind=IntKind) :: sz_ind_intkkr
   integer (kind=IntKind) :: sz_ind_kkrkkr
   integer (kind=IntKind) :: sz_ind_iendphikkr
   integer (kind=IntKind) :: sz_ind_kkrkkrns
   integer (kind=IntKind) :: sz_ind_step
!
   integer (kind=IntKind), allocatable :: ind_intkkr(:,:,:)
   integer (kind=IntKind), allocatable :: ind_kkrkkr(:,:,:)
   integer (kind=IntKind), allocatable :: ind_intint(:,:,:)
   integer (kind=IntKind), allocatable :: ind_iendphikkr(:,:,:)
   integer (kind=IntKind), allocatable :: ind_kkrkkrns(:,:,:)
   integer (kind=IntKind), allocatable :: ind_step(:)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: smdexp
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   interface 
      function getNumSpecies(id) result(n)
         use KindParamModule, only : IntKind
         implicit none
         integer (kind=IntKind), intent(in) :: id
         integer (kind=IntKind) :: n
      end function getNumSpecies
!
      function getAtomicNumber(id,ic) result(n)
         use KindParamModule, only : IntKind
         implicit none
         integer (kind=IntKind), intent(in) :: id
         integer (kind=IntKind), intent(in), optional :: ic
         integer (kind=IntKind) :: n
      end function getAtomicNumber
   end interface
!
   if (Initialized) then
      call ErrorHandler('initSSSolver','intend to initialize twice')
   else if (na < 1) then
      call ErrorHandler('initSSSolver','LocalNumSites < 1',na)
   else if (nsc < 1 .or. nsc > 2) then
      call ErrorHandler('initSSSolver','invalid spin canting index',nsc)
   else if (nsp < 1 .or. nsp > 2 .or. nsp < nsc) then
      call ErrorHandler('initSSSolver','invalid spin polarization index',nsp)
   endif
!
   LocalNumSites = na
   NumSpins = nsc
!
   stop_routine = istop
!
   if (rel == 0) then           ! non-relativistic
      c2inv = ZERO
      Relativity = NonRelativistic
   else if (rel == 1) then      ! semi-relativistic
      c2inv = ONE/(LightSpeed*LightSpeed)
      Relativity = ScalarRelativistic
   else if (rel == 2) then      ! fully-relativistic
      c2inv = ONE/(LightSpeed*LightSpeed)
      Relativity = FullyRelativistic
   else
      call ErrorHandler('initSSSolver',                       &
                        'invalid relativity parameter',rel)
   endif

   solv_type = getSingleSiteSolverType()
   SSSMethod = getSingleSiteSolverMethod()
   if ( SSSMethod < -1 .or. SSSMethod > 2) then
      SSSMethod = 2
   endif
!
   if ( solv_type==0 .or. solv_type==1 ) then
      FullSolver = .false.
   else if (solv_type==2 .or. solv_type==3 ) then
      FullSolver = .true.
   endif
   isSphericalSolverOn = .not.FullSolver
!
   if ( solv_type==0 .or. solv_type==2 ) then
      if (rel>1) then
         call ErrorHandler("initSSSolver",                    &
              "Inconsitency in use of solver type and relativity",    &
               solv_type,rel)
      endif
   else if (solv_type==1 .or. solv_type==3 ) then
      if (rel<2) then
         call ErrorHandler("initSSSolver",                    &
              "Inconsitency in use of solver type and relativity",    &
               solv_type,rel)
      endif
   else
     call ErrorHandler("initSSSolver",                        &
                       "Undefined Single Site Solver", solv_type,rel)
   endif
!
   if (present(derivative)) then
      rad_deriv = derivative
   else
      rad_deriv = .false.
   endif
!
   allocate( Scatter(LocalNumSites) )
   allocate( print_instruction(LocalNumSites) )
!
   MaxSpecies = 1
   do ia=1,LocalNumSites
      n = getNumSpecies(ia)
      Scatter(ia)%NumSpecies = n
      allocate(Scatter(ia)%AtomicNumber(n), Scatter(ia)%Solutions(n,NumSpins))
      do i = 1, n
         Scatter(ia)%AtomicNumber(i) = getAtomicNumber(ia,i)
         if (Scatter(ia)%AtomicNumber(i) < 0 .or. Scatter(ia)%AtomicNumber(i) > 108) then
            call ErrorHandler('initSSSolver','invalid atomic number',Scatter(ia)%AtomicNumber(i))
         endif
      enddo
      MaxSpecies = max(n, MaxSpecies)
   enddo
!
   allocate( ind_intkkr(MaxSpecies, NumSpins, LocalNumSites) )
   allocate( ind_kkrkkr(MaxSpecies, NumSpins, LocalNumSites) )
   allocate( ind_intint(MaxSpecies, NumSpins, LocalNumSites) )
   allocate( ind_iendphikkr(MaxSpecies, NumSpins, LocalNumSites) )
   allocate( ind_kkrkkrns(MaxSpecies, NumSpins, LocalNumSites) )
   allocate( ind_step(LocalNumSites) )
!
   sz_ind_intkkr     = 0
   sz_ind_kkrkkr      = 0
   sz_ind_intint      = 0
   sz_ind_iendphikkr = 0
   sz_ind_kkrkkrns    = 0
   sz_ind_step        = 0
!
   lmax = 0
   lmax_phi = 0
   jmax_max = 0
   iend_max = 0
   jmax_trunc_max = 0
   lmax_max_kkr = 0
   kmax_max_kkr = 0
   kmax_max_phi = 0
   kmax_max_int = 0
   lmax_max_green = 0
   max_NDiffSurfR = 0
   ng_max = 0
   lmax_potmp = 0
!
   do ia=1,LocalNumSites
      print_instruction(ia) = iprint(ia)
      if (lkkr(ia) < 0) then
         call ErrorHandler('initSSSolver','lkkr < 0',lkkr(ia))
      else if (lkkr(ia) > lphi(ia)) then
         call ErrorHandler('initSSSolver','lkkr > lphi',      &
                           lkkr(ia),lphi(ia))
      else if (lpot(ia) < 0) then
         call ErrorHandler('initSSSolver','lpot < 0',lpot(ia))
      else if (lpot(ia) > getPotLmax(ia)) then
         call ErrorHandler('initSSSolver','lpot > lmax_pot',lpot(ia))
      endif
      Grid => getGrid(ia)
!
      if ( FullSolver ) then
         iend = Grid%jend_plus_n
         lmax_pot = lpot(ia)
      else
         iend = min( max(Grid%jmt+11,Grid%jend),Grid%jend_plus_n )
         lmax_pot = 0
      endif
      iend_max = max(iend_max,iend)
      lmax = max(lphi(ia)+lkkr(ia), lpot(ia), lmax)
      lmax_trunc = getTruncPotLmax(ia)
      lmax_potmp = max(lmax,lmax_trunc)
      lmax_step = lstep(ia)
      jmax_step = ((lmax_step+1)*(lmax_step+2))/2
!
      jmax_pot = ((lmax_pot+1)*(lmax_pot+2))/2
      jmax_max = max(jmax_max,jmax_pot)
      if ( lmax_trunc>-1 ) then 
         jmax_trunc = ((lmax_trunc+1)*(lmax_trunc+2))/2
         jmax_trunc_max = max(jmax_trunc_max,jmax_trunc)
      else
         jmax_trunc = -1
         jmax_trunc_max = -1
      endif
!
      if (SSSMethod == -1) then
         NDiffSurfR = getNumDiffSurfR(ia)
         max_NDiffSurfR = max(max_NDiffSurfR,NDiffSurfR)
      endif
!
      ng = getNumGaussRs(ia)
      ng_max = max(ng_max,ng)
!
      lmax_kkr = lkkr(ia)
      lmax_phi = max(lmax_phi,lphi(ia))
      kmax_kkr = (lkkr(ia)+1)*(lkkr(ia)+1)
      kmax_phi = (lphi(ia)+1)*(lphi(ia)+1)
      kmax_int = getKmaxIntern(kmax_kkr,kmax_phi)
      lmax_max_kkr = max(lmax_max_kkr,lmax_kkr)
      kmax_max_kkr = max(kmax_max_kkr,kmax_kkr)
      kmax_max_phi = max(kmax_max_phi,kmax_phi)
      kmax_max_int = max(kmax_max_int,kmax_int)
      if ( .not. FullSolver ) then
         lmax_pot = 0
      else
         lmax_pot = lpot(ia)
      endif
!
      allocate( Scatter(ia)%lm_index_sph(kmax_kkr,Scatter(ia)%NumSpecies) )
      Scatter(ia)%lm_index_sph = 0
      do ic = 1, Scatter(ia)%NumSpecies
         index_end = 0
         kl = 0
         do l = 0, lmax_kkr
            if (checkLdaCorrection(ia,ic,l)) then
               do m = -l, l
                  kl = kl + 1
                  index_end = index_end + 1
                  Scatter(ia)%lm_index_sph(kl,ic) = index_end
               enddo
            else
               index_end = index_end + 1
               do m = -l, l
                  kl = kl + 1
                  Scatter(ia)%lm_index_sph(kl,ic) = index_end
               enddo
            endif
         enddo
      enddo
!
      allocate( Scatter(ia)%sol_flags(kmax_phi,kmax_kkr) )
      Scatter(ia)%sol_flags = 0
!
      Scatter(ia)%numrs    = iend
      Scatter(ia)%numrs_cs = Grid%jend
      Scatter(ia)%numrs_mt = Grid%jmt
      Scatter(ia)%numrs_trunc = Grid%jend-Grid%jmt+1 ! includes the muffin-tin point
      Scatter(ia)%lmax_phi = lphi(ia)
      Scatter(ia)%lmax_kkr = lkkr(ia)
      Scatter(ia)%lmax_pot = lmax_pot
      Scatter(ia)%kmax_kkr = kmax_kkr
      Scatter(ia)%kmax_phi = kmax_phi
      Scatter(ia)%kmax_int = kmax_int
      Scatter(ia)%jmax_pot = jmax_pot
      Scatter(ia)%lmax_pot_trunc = lmax_trunc
      Scatter(ia)%jmax_pot_trunc = jmax_trunc
      Scatter(ia)%lmax_step = lmax_step
      Scatter(ia)%jmax_step = jmax_step
      Scatter(ia)%lmax_green = lgreen(ia)
      Scatter(ia)%kmax_green = (lgreen(ia)+1)**2
!
      lmax_max_green = max(lmax_max_green,Scatter(ia)%lmax_green)
!
      do ic = 1, Scatter(ia)%NumSpecies
         do is=1,NumSpins
            ind_intkkr(ic,is,ia)     = sz_ind_intkkr + 1
            ind_kkrkkr(ic,is,ia)      = sz_ind_kkrkkr + 1
            ind_intint(ic,is,ia)      = sz_ind_intint + 1
            ind_iendphikkr(ic,is,ia) = sz_ind_iendphikkr + 1
            ind_kkrkkrns(ic,is,ia)    = sz_ind_kkrkkrns + 1
!
            sz_ind_intkkr     = sz_ind_intkkr + kmax_kkr*kmax_int
            sz_ind_kkrkkr      = sz_ind_kkrkkr  + kmax_kkr*kmax_kkr
            sz_ind_intint      = sz_ind_intint  + kmax_int*kmax_int
            sz_ind_iendphikkr = sz_ind_iendphikkr + iend*kmax_kkr*kmax_phi
            sz_ind_kkrkkrns    = sz_ind_kkrkkrns + kmax_kkr*kmax_kkr*NumSpins
         enddo
      enddo
      ind_step(ia) = sz_ind_step + 1
      sz_ind_step  = sz_ind_step + Scatter(ia)%numrs_trunc*jmax_step
   enddo
   jmax_pot_solver = max(jmax_max,jmax_trunc_max)
   lmax_pot_in = getLmaxPotSolver() ! = -1, if it is undefined from the input
   jmax_pot_in = ((lmax_pot_in+1)*(lmax_pot_in+2))/2
   jmax_potmp = (lmax_potmp+1)*(lmax_potmp+2)/2
   if ( lmax_pot_in >= 0 ) then
      jmax_pot_solver = min(jmax_pot_solver,jmax_pot_in)
   endif
   jmax_pot_solver_default = jmax_pot_solver
   lmax_sol_cutoff = getLmaxSolution()
!
   allocate( WronskianIndx(kmax_max_phi,kmax_max_kkr) )
   allocate( wks_sinmat(sz_ind_intkkr) )
   allocate( wks_cosmat(sz_ind_intkkr) )
!
   allocate( wks_tmat(sz_ind_intint) )
   allocate( wks_tmat_inv(sz_ind_intint) )
   allocate( wks_jostmat(sz_ind_intkkr) )
   allocate( wks_jinvmat(sz_ind_intkkr) )
   allocate( wks_Omega(sz_ind_kkrkkr) )
   allocate( wks_OmegaHat(sz_ind_kkrkkr) )
   allocate( wks_OmegaHatInv(sz_ind_kkrkkr) )
   allocate( wks_S(sz_ind_kkrkkr) )
   allocate( wks_PS(2*kmax_max_phi*LocalNumSites*NumSpins*MaxSpecies) )
!
!tmat_global   if (NumSpins == 2) then
!tmat_global      allocate(wks_tmatg(sz_ind_intint*4))
!tmat_global   endif
!
   allocate( wks_regsol(sz_ind_iendphikkr), stat=INFO )
   if (INFO > 0) then
      write(6,'(a,2i8)')'Error in allocating wks_regsol: ',MyPE,sz_ind_iendphikkr
   endif
   allocate( wks_step(sz_ind_step), stat=INFO )
   if (INFO > 0) then
      write(6,'(a,2i8)')'Error in allocating wks_step: ',MyPE,sz_ind_step
   endif
!
   allocate( wks_plhat(iend_max*kmax_max_phi), stat=INFO )
   if (INFO > 0) then
      write(6,'(a,2i8)')'Error in allocating wks_plhat: ',MyPE,iend_max*kmax_max_phi
   endif
   allocate( wks_qlhat(iend_max*kmax_max_phi), stat=INFO )
   if (INFO > 0) then
      write(6,'(a,2i8)')'Error in allocating wks_qlhat: ',MyPE,iend_max*kmax_max_phi
   endif
   allocate( wks_mtx1(kmax_max_phi*kmax_max_kkr+kmax_max_phi+kmax_max_kkr), stat=INFO )
   if (INFO > 0) then
      write(6,'(a,2i8)')'Error in allocating wks_mtx1: ',MyPE,kmax_max_phi*kmax_max_kkr+kmax_max_phi+kmax_max_kkr
   endif
   allocate( wks_mtx2(kmax_max_phi*kmax_max_kkr+kmax_max_phi+kmax_max_kkr), stat=INFO )
   if (INFO > 0) then
      write(6,'(a,2i8)')'Error in allocating wks_mtx2: ',MyPE,kmax_max_phi*kmax_max_kkr+kmax_max_phi+kmax_max_kkr
   endif
   allocate( wks_dregsol(sz_ind_iendphikkr), stat=INFO )
   if (INFO > 0) then
      write(6,'(a,2i8)')'Error in allocating wks_dregsol: ',MyPE,sz_ind_iendphikkr
   endif
!
   allocate( wks_int_lsph((iend_max+ng_max)*(lmax_max_kkr+1)) )
   allocate( wks_int_rsph((iend_max+ng_max)*(lmax_max_kkr+1)) )
   allocate( wks_int_rtmp(ng_max+iend_max) )
   allocate( wks_int_res(iend_max+ng_max), wks_int_fg(iend_max+ng_max) )
   allocate( wks_int_bjl(ng_max*(lmax_max_kkr+1)),  &
             wks_int_bnl(ng_max*(lmax_max_kkr+1)) )
   allocate( IPVT(kmax_max_phi) )
!
   n = 1
!tmat_global   sz_tmatg = 1
   do ia=1,LocalNumSites
      Grid => getGrid(ia)
      r_mesh => Grid%r_mesh
!
      iend     = Scatter(ia)%numrs
      kmax_kkr = Scatter(ia)%kmax_kkr
      kmax_phi = Scatter(ia)%kmax_phi
      kmax_int = Scatter(ia)%kmax_int
!
      do is=1,NumSpins
         do ic = 1, Scatter(ia)%NumSpecies
            skk_1 = ind_kkrkkr(ic,is,ia)
            skk_2 = ind_kkrkkr(ic,is,ia) + kmax_kkr*kmax_kkr-1
            spp_1 = ind_intint(ic,is,ia)
            spp_2 = ind_intint(ic,is,ia) + kmax_int*kmax_int-1
            spk_1 = ind_intkkr(ic,is,ia)
            spk_2 = ind_intkkr(ic,is,ia) + kmax_int*kmax_kkr-1
            p1c => wks_sinmat(spk_1:spk_2)
            Scatter(ia)%Solutions(ic,is)%sin_mat => aliasArray2_c( p1c, kmax_int, kmax_kkr )
            p1c => wks_cosmat(spk_1:spk_2)
            Scatter(ia)%Solutions(ic,is)%cos_mat => aliasArray2_c( p1c, kmax_int, kmax_kkr )
            p1c => wks_tmat(spp_1:spp_2)
            Scatter(ia)%Solutions(ic,is)%t_mat => aliasArray2_c( p1c, kmax_int, kmax_int)
            p1c => wks_tmat_inv(spp_1:spp_2)
            Scatter(ia)%Solutions(ic,is)%t_mat_inv => aliasArray2_c( p1c, kmax_int, kmax_int)
            p1c => wks_jostmat(spk_1:spk_2)
            Scatter(ia)%Solutions(ic,is)%jost_mat => aliasArray2_c( p1c, kmax_int, kmax_kkr)
            p1c => wks_jinvmat(spk_1:spk_2)
            Scatter(ia)%Solutions(ic,is)%jinv_mat => aliasArray2_c( p1c, kmax_int, kmax_kkr)
            p1c => wks_Omega(skk_1:skk_2)
            Scatter(ia)%Solutions(ic,is)%Omega_mat => aliasArray2_c(p1c, kmax_kkr, kmax_kkr)
            p1c => wks_OmegaHat(skk_1:skk_2)
            Scatter(ia)%Solutions(ic,is)%OmegaHat_mat => aliasArray2_c(p1c, kmax_kkr, kmax_kkr)
            p1c => wks_OmegaHatInv(skk_1:skk_2)
            Scatter(ia)%Solutions(ic,is)%OmegaHatInv_mat => aliasArray2_c(p1c, kmax_kkr, kmax_kkr)
            p1c => wks_S(skk_1:skk_2)
            Scatter(ia)%Solutions(ic,is)%S_mat => aliasArray2_c(p1c, kmax_kkr, kmax_kkr)
            Scatter(ia)%Solutions(ic,is)%phase_shift => wks_PS(n:n+kmax_phi-1)
            n = n + kmax_phi
!
            sz_1 = ind_iendphikkr(ic,is,ia)
            sz_2 = ind_iendphikkr(ic,is,ia) + iend*kmax_phi*kmax_kkr-1
            p1c => wks_regsol(sz_1:sz_2)
            Scatter(ia)%Solutions(ic,is)%reg_sol => aliasArray3_c( p1c, iend, kmax_phi, kmax_kkr)
            p1c => wks_dregsol(sz_1:sz_2)
            Scatter(ia)%Solutions(ic,is)%reg_dsol => aliasArray3_c( p1c, iend, kmax_phi, kmax_kkr)
!
            Scatter(ia)%Solutions(ic,is)%reg_sol = CZERO
            Scatter(ia)%Solutions(ic,is)%reg_dsol = CZERO
         enddo
      enddo
!
!tmat_global      if (NumSpins == 2) then
!tmat_global         nsz = kmax_int*kmax_int*Scatter(ia)%NumSpecies
!                    p1c => wks_tmatg(sz_tmatg:sz_tmatg+nsz)
!tmat_global         allocate(Scatter(ia)%tmat_global => aliasArray3_c( p1c,kmax_int,kmax_int, &
!tmat_global                                                            Scatter(ia)%NumSpecies)
!tmat_global         sz_tmatg = sz_tmatg + nsz
!tmat_global      endif
!
      sz_1 = ind_step(ia)
      sz_2 = ind_step(ia) -1 + Scatter(ia)%numrs_trunc*Scatter(ia)%jmax_step
      p1c => wks_step(sz_1:sz_2)
      Scatter(ia)%step_rad => aliasArray2_c( p1c, Scatter(ia)%numrs_trunc, Scatter(ia)%jmax_step )
      if ( FullSolver ) then
         lmax_step = Scatter(ia)%lmax_step
         call getRadialStepFunction( ia, Scatter(ia)%numrs_mt, Scatter(ia)%numrs_cs, r_mesh, &
                                     lmax_step, Scatter(ia)%step_rad )
      else
         Scatter(ia)%step_rad = CZERO
      endif
   enddo
!
   deallocate( ind_intkkr, ind_kkrkkr, ind_intint, ind_iendphikkr, ind_kkrkkrns, ind_step )
!
   nullify( Grid, r_mesh )
!
   allocate( v0r(iend_max), cm0(iend_max) )
   allocate( v0r_save(iend_max) )
   allocate( bjl(iend_max,0:lmax_phi+1), bnl(iend_max,0:lmax_phi+1),  &
             bhl(iend_max,0:lmax_phi+1) ) 
   allocate( dbjl(iend_max,0:lmax_phi), dbnl(iend_max,0:lmax_phi),    &
             dbhl(iend_max,0:lmax_phi) )
   allocate( bjtmp(0:lmax_phi+1), bntmp(0:lmax_phi+1) )
   allocate( TmpSpace(1:5*iend_max) )
   allocate( flags_jl(jmax_max) )
!
   if ( FullSolver ) then
      if ( SSSMethod == -1 ) then
         allocate( wksl_smpicm(1:kmax_max_kkr*kmax_max_kkr) )
         allocate( JHunt(1:max_NDiffSurfR) )
         allocate( sjl(0:lmax_phi,1:max_NDiffSurfR),                    &
                   snl(0:lmax_phi,1:max_NDiffSurfR),                    &
                   dsjl(0:lmax_phi,1:max_NDiffSurfR),                   &
                   dsnl(0:lmax_phi,1:max_NDiffSurfR) )
         allocate( nlrmt(0:lmax_phi), dnlrmt(0:lmax_phi))
         allocate( wksl_plsurf(kmax_max_phi*max_NDiffSurfR) )
         allocate( wksl_qlsurf(kmax_max_phi*max_NDiffSurfR) )
      else if ( SSSMethod == 0 ) then
         npout = 1; kpm = 1; kpsm = 1
         do ia=1,LocalNumSites
            npout = max(npout, Scatter(ia)%numrs_cs - Scatter(ia)%numrs_mt + 1)
            kpm = max( kpm, min((2*Scatter(ia)%lmax_phi+1)**2,(Scatter(ia)%lmax_pot+Scatter(ia)%lmax_step+1)**2), &
                           (Scatter(ia)%lmax_pot+1)**2, (Scatter(ia)%lmax_step+1)**2 )
            kpsm = max( kpsm,(Scatter(ia)%lmax_pot+Scatter(ia)%lmax_step+1)**2 )
         enddo
         allocate( wks_scvphi(npout*kpm,2) )
         allocate( cpotstep(npout,kpsm) )
      else if ( SSSMethod == 1 ) then
         npout = 1; kpm = 1
         do ia=1,LocalNumSites
            npout = max(npout, Scatter(ia)%numrs_cs - Scatter(ia)%numrs_mt + 1)
            kpm = max(kpm, (Scatter(ia)%lmax_phi+1)**2,(Scatter(ia)%lmax_pot+1)**2)
         enddo
         allocate( wks_scvphi(npout*kpm,2) )
      endif
!
      if ( SSSMethod > 0 ) then
         allocate( flags_trunc_jl(jmax_trunc_max) )
      endif
!
   endif
!  ===================================================================
!  -------------------------------------------------------------------
   call initIntegerFactors(3*lmax)
!  -------------------------------------------------------------------
!
#ifdef TMP_ACCEL
   call pushIntegerFactorsToAccel()
#endif
   lmax = lofk(kmax_max_kkr)
   allocate( pl_regwrk(1:iend_max,0:lmax),    &
             ql_regwrk(1:iend_max,0:lmax) )
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   kmax = max(kmax_max_kkr,kmax_max_phi)
   lmax = lofk(kmax)
   allocate( factor_L(0:lmax) )
!
   lmax = lofk(kmax_max_kkr)
   allocate( s_sph_mt(0:lmax), c_sph_mt(0:lmax), s_mtrx(0:lmax),  c_mtrx(0:lmax) )
!
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call initialize_calphilr(c2inv,iend_max,lmax_phi,                  &
                            lofj(jmax_max),lofj(jmax_trunc_max),      &
                            kmax_max_kkr,kmax_max_phi)
!  -------------------------------------------------------------------
   call pushGauntFactorsToAccel()
   call pushRadialGridToAccel()
#endif
!
!
   energy = CZERO
   PotShift = CZERO
   Initialized = .true.
   IrrSpaceAllocated = .false.
   isIrrSolOn = .false.
   PrintPotT_done = .false.
   spin_index = 0
!
   if ( maxval(print_instruction) >= 0 ) then
      write(6,'(/,80(''-''))')
      write(6,'(a)')'      ************************************'
      write(6,'(a)')'      *     Output from initSSSolver     *'
      write(6,'(a)')'      ************************************'
!
      write(6,'(a,i5)') "   Number of local Atoms:   ", LocalNumSites
      write(6,'(a,i5)') "   Number of Spins:         ", NumSpins
      write(6,*)"   Single Site Solver methods: "
      if ( solv_type==0 .or. solv_type==1 ) then
         write(6,'(a)')"      1. Spherical-Potential Solver"
      else if (solv_type==2 .or. solv_type==3 ) then
         write(6,'(a)')"      1. Full-Potential Solver"
         if ( SSSMethod==-1 ) then
            write(6,'(a)')"      2. Solver for untruncated potential, surface integration"
         else if ( SSSMethod==0 ) then
            write(6,'(a)')"      2. Solver for untruncated potential, volume integration for untruncated potential"
         else if ( SSSMethod==1 ) then
            write(6,'(a)')"      2. Solver for untruncated potential, volume integration for truncated potential"
         else
            write(6,'(a)')"      2. Solver for truncated potential"
         endif
      endif
!
      write(6,'(a,i5)')"         jmax potential solver:    ", jmax_pot_solver
      if (lmax_pot_in >= 0) then
         write(6,'(a,i5)')"         jmax potential solver in: ", jmax_pot_in
      endif
      write(6,'(a,i5)')"         jmax_trunc_max :          ", jmax_trunc_max
      write(6,'(a,i5)')"         lmax_max_kkr   :          ", lofk(kmax_max_kkr)
      write(6,'(a,i5)')"         lmax_max_phi   :          ", lofk(kmax_max_phi)
      write(6,'(a,i5)')"         lmax_max_pot   :          ", lofj(jmax_max)
      write(6,'(a,i5)')"         lmax_max_green :          ", lmax_max_green
      if ( SSSMethod > 0 ) then
         write(6,'(a,i5)')"         lmax_pot_trunc :          ", lofj(jmax_trunc_max)
      else
         write(6,'(a,i5)')"         lmax_pot_trunc :          ", lofj(jmax_max)
      endif
      write(6,'(40(''=''))')
      call FlushFile(6)
   endif
!
!  ===================================================================
!  Use the existing "K-Mesh" MPI group to create a parallelization
!  over the kl loop in the single site solver...
!  -------------------------------------------------------------------
   kGID = getGroupID('K-Mesh')
   NumPEsInGroup = getNumPEsInGroup(kGID)
   MyPEinGroup = getMyPEinGroup(kGID)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Determine the starting phase shift at e0 = (1.0d-6, 0.0d0)
!  ===================================================================
   allocate( phase_ref(kmax_max_phi,MaxSpecies,nsp,LocalNumSites),   &
             npi(kmax_max_phi,MaxSpecies,nsp,LocalNumSites) )
   do ia = 1, LocalNumSites
      do is = 1, nsp ! Loop over spin polarization index
         js = min(is, NumSpins)
!        -------------------------------------------------------------
         call solveSingleScattering(is, ia, CZERO+TEN2m6, CZERO)
!        -------------------------------------------------------------
         do ic = 1, Scatter(ia)%NumSpecies
!           ==========================================================
!           Triangularize the S_matrix
!           ----------------------------------------------------------
            call ZGETRF_nopivot(Scatter(ia)%kmax_phi, Scatter(ia)%kmax_phi, &
                                Scatter(ia)%Solutions(ic,js)%S_mat,         &
                                Scatter(ia)%kmax_phi, IPVT, INFO)
!           ----------------------------------------------------------
            do kl = 1, Scatter(ia)%kmax_phi
               smdexp = log(Scatter(ia)%Solutions(ic,js)%S_mat(kl,kl))
               phase_ref(kl,ic,is,ia) = aimag(smdexp)/TWO
               npi(kl,ic,is,ia) = 0
            enddo
         enddo
      enddo
   enddo
!
   isDosSymmOn = isChargeSymm()
!
   end subroutine initSSSolver
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSSSolver()
!  ===================================================================
use MPPModule, only : MyPE, syncAllPEs
   use IntegerFactorsModule, only : endIntegerFactors, deleteIntegerFactorsOnAccel
   use GauntFactorsModule, only: deleteGauntFactorsOnAccel
   use PotentialModule, only: deletePotentialOnAccel,deleteTruncatedPotOnAccel
   use RadialGridModule, only: deleteRadialGridOnAccel
   implicit  none
!
   integer (kind=IntKind) :: ia, is, ic
!
   if (.not.Initialized) then
      call ErrorHandler('endSSSolver','module not initialized')
   endif
!
   nullify( DiffSurfR, GaussR2SurfR, n_dot_r, weight, r_surf, nDotGradY )
   nullify( ylm, pot_jl, nj3, kj3, cgnt )
   nullify( smpicm, pl_reg, ql_reg, pl_irr, ql_irr, qlhat_reg, qlhat_irr, &
            plhat_reg, plhat_irr, pot_trunc_jl )
!
!  ===================================================================
!  free the allocated space to save memory.........................
!  -------------------------------------------------------------------
   do ia = 1,LocalNumSites
      do is = 1, NumSpins
         do ic = 1, Scatter(ia)%NumSpecies
            nullify( Scatter(ia)%Solutions(ic,is)%sin_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%cos_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%t_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%t_mat_inv )
            nullify( Scatter(ia)%Solutions(ic,is)%jost_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%jinv_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%Omega_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%OmegaHat_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%OmegaHatInv_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%S_mat )
            nullify( Scatter(ia)%Solutions(ic,is)%reg_sol )
            nullify( Scatter(ia)%Solutions(ic,is)%irr_sol )
            nullify( Scatter(ia)%Solutions(ic,is)%reg_dsol )
            nullify( Scatter(ia)%Solutions(ic,is)%irr_dsol )
            nullify( Scatter(ia)%Solutions(ic,is)%green )
            nullify( Scatter(ia)%Solutions(ic,is)%dos )
            nullify( Scatter(ia)%Solutions(ic,is)%pdos )
            nullify( Scatter(ia)%Solutions(ic,is)%phase_shift )
         enddo
      enddo
      nullify( Scatter(ia)%step_rad )
!tmat_global    nullify( Scatter(ia)%tmat_global )
      deallocate( Scatter(ia)%lm_index_sph, Scatter(ia)%sol_flags )
      deallocate( Scatter(ia)%AtomicNumber, Scatter(ia)%Solutions )
   enddo
   deallocate( Scatter, print_instruction )
   deallocate( phase_ref, npi )
!  ------------------------------------------------------------------
   call endIntegerFactors()
#ifdef TMP_ACCEL
   call deleteIntegerFactorsOnAccel()
   call delete_calphilr_data()
   call deleteGauntFactorsOnAccel()
   call deleteTruncatedPotOnAccel()
   call deletePotentialOnAccel()
   call deleteRadialGridOnAccel()
#endif
!  ------------------------------------------------------------------
!
   deallocate( wks_sinmat, wks_cosmat )
!
   deallocate( wks_tmat, wks_tmat_inv, wks_jostmat, wks_jinvmat )
   deallocate( wks_Omega, wks_OmegaHat, wks_OmegaHatInv, wks_S, wks_PS, wks_step )
!
!tmat_global   if (NumSpins == 2) then
!tmat_global      deallocate( wks_tmatg )
!tmat_global   endif
!
   deallocate( wks_plhat, wks_qlhat )
   deallocate( wks_regsol )
!
   deallocate( IPVT )
   deallocate( wks_mtx1, wks_mtx2 )
   deallocate( pl_regwrk, ql_regwrk )
   deallocate( wks_dregsol )
!
   if (allocated(wks_green)) then
      deallocate( wks_green )
   endif
   if (allocated(wks_dgreen)) then
      deallocate( wks_dgreen )
   endif
   if (allocated(wks_dos)) then
      deallocate( wks_dos )
   endif
   if (allocated(wks_ddos)) then
      deallocate( wks_ddos )
   endif
   if (allocated(wks_pdos)) then
      deallocate( wks_pdos )
   endif
!
   deallocate( bjl, bnl, bhl )
   deallocate( dbjl, dbnl, dbhl )
   deallocate( bjtmp, bntmp )
   deallocate( TmpSpace )
   deallocate(flags_jl)
!
   if ( FullSolver ) then
      if ( SSSMethod == -1 ) then
         deallocate( wksl_smpicm )
         deallocate( JHunt )
         deallocate( sjl, snl, dsjl, dsnl )
         deallocate( nlrmt, dnlrmt )
         deallocate( wksl_plsurf, wksl_qlsurf )
      else if ( SSSMethod == 0 ) then
         deallocate( cpotstep )
      else
         deallocate( flags_trunc_jl )
      endif
      if (SSSMethod == 0 .or. SSSMethod == 1) then
         deallocate( wks_scvphi )
      endif
   endif
   deallocate( v0r, cm0 )
   deallocate( v0r_save )
!
   deallocate( factor_L )
!
   deallocate( WronskianIndx )
   deallocate( s_sph_mt, c_sph_mt, s_mtrx,  c_mtrx )
!
   deallocate( wks_int_rtmp, wks_int_lsph, wks_int_rsph )
!
   deallocate( wks_int_res, wks_int_fg, wks_int_bjl, wks_int_bnl )
!
   if (allocated(space_integrated_dos_cell) ) then
      deallocate( space_integrated_dos_cell, space_integrated_dos_mt, &
                  space_integrated_pdos_cell, space_integrated_pdos_mt )
   endif
!
   if ( IrrSpaceAllocated ) then
      call deallocateIrrSpace()
   endif
!
   energy = CZERO
   PotShift = CZERO
   isDosSymmOn = .false.
   Initialized = .false.
!
   end subroutine endSSSolver
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initGlobalVariables(spin, site, e, vshift)
!  ===================================================================
   use RadialGridModule, only : getGrid
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin
   integer (kind=IntKind), intent(in) :: site
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in) :: vshift
!
   spin_index = spin
   LocalIndex = site
!
   energy = e
   kappa = sqrt(energy+energy*energy*c2inv)
   PotShift = vshift
!
!  -------------------------------------------------------------------
   Grid => getGrid(LocalIndex)
!  -------------------------------------------------------------------
!
   LocalSpin = min(spin, NumSpins)    ! if non-spin-canted, LocalSpin is always 1.
!
   lmax_kkr = Scatter(site)%lmax_kkr
   kmax_kkr = Scatter(site)%kmax_kkr
   lmax_phi = Scatter(site)%lmax_phi
   kmax_phi = Scatter(site)%kmax_phi
   kmax_int = Scatter(site)%kmax_int
   lmax_pot = Scatter(site)%lmax_pot
   jmax_pot = Scatter(site)%jmax_pot
   lmax_step = Scatter(site)%lmax_step
   jmax_step = Scatter(site)%jmax_step
   jmax_trunc = Scatter(site)%jmax_pot_trunc
   lmax_trunc = Scatter(site)%lmax_pot_trunc
   lmax_green = Scatter(site)%lmax_green
   kmax_green = Scatter(site)%kmax_green
   iend = Scatter(site)%numrs
   numrs_mt = Scatter(site)%numrs_mt
   numrs_cs = Scatter(site)%numrs_cs
   numrs_trunc = Scatter(site)%numrs_trunc
!
   TmpSpace = CZERO
!
   end subroutine initGlobalVariables
!  ===================================================================
!
!  *******************************************************************
!  
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKmaxIntern(kkr, phi) result(kl)
!  ===================================================================
   implicit none
!  
   integer (kind=IntKind), intent(in) :: kkr, phi
!
   integer (kind=IntKind) :: kl
!  
   kl = kkr     ! set internal kmax to the kmax for kkr
!
   end function getKmaxIntern
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocateIrrSpace()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: ia, is, ic, lmax_phi_max
   integer (kind=IntKind) :: iend, kmax_kkr, kmax_phi, sz, sz_1, sz_2
!
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   lmax_phi_max = 0
   sz = 0
   do ia=1,LocalNumSites
      iend     = Scatter(ia)%numrs
      kmax_kkr = Scatter(ia)%kmax_kkr
      kmax_phi = Scatter(ia)%kmax_phi
      lmax_phi_max = max(lmax_phi_max,Scatter(ia)%lmax_phi)
      sz = sz + iend*kmax_phi*kmax_kkr*NumSpins*Scatter(ia)%NumSpecies
   enddo
!
   allocate( wks_irrsol(sz), wks_dirrsol(sz) )
   allocate( pl_irrwrk(1:iend_max,0:lmax_phi_max), ql_irrwrk(1:iend_max,0:lmax_phi_max))
!
   sz_2 = 0
   do ia=1,LocalNumSites
      iend     = Scatter(ia)%numrs
      kmax_kkr = Scatter(ia)%kmax_kkr
      kmax_phi = Scatter(ia)%kmax_phi
      do is=1,NumSpins
         do ic = 1, Scatter(ia)%NumSpecies
            sz_1 = sz_2 + 1
            sz_2 = sz_1 - 1 + iend*kmax_phi*kmax_kkr
            p1c => wks_irrsol(sz_1:sz_2)
            Scatter(ia)%Solutions(ic,is)%irr_sol => aliasArray3_c( p1c, iend, kmax_phi, kmax_kkr)
            p1c => wks_dirrsol(sz_1:sz_2)
            Scatter(ia)%Solutions(ic,is)%irr_dsol => aliasArray3_c( p1c, iend, kmax_phi, kmax_kkr)
            Scatter(ia)%Solutions(ic,is)%irr_sol = CZERO
            Scatter(ia)%Solutions(ic,is)%irr_dsol= CZERO
         enddo
      enddo
   enddo
!
   if ( FullSolver .and. SSSMethod == -1 ) then
      allocate( wksl_invscm(1:kmax_max_kkr*kmax_max_kkr) )
      allocate( wksl_invscmt(1:kmax_max_kkr*kmax_max_kkr) )
      allocate( wksl_smpicmt(1:kmax_max_kkr*kmax_max_kkr) )
      allocate( wksl_tmpirr(1:iend_max*kmax_max_phi*kmax_max_kkr) )
      allocate( al(1:kmax_max_kkr), bl(1:kmax_max_kkr) )
      allocate( scp(1:kmax_max_kkr), hp(1:kmax_max_kkr) )
   endif
!
   IrrSpaceAllocated = .true.
!
   end subroutine allocateIrrSpace
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deallocateIrrSpace()
!  ===================================================================
   implicit none
!
   deallocate( wks_irrsol, wks_dirrsol )
   deallocate( pl_irrwrk, ql_irrwrk)
!
   if ( FullSolver .and. SSSMethod == -1 ) then
      deallocate( wksl_invscm )
      deallocate( wksl_invscmt )
      deallocate( wksl_smpicmt )
      deallocate( wksl_tmpirr )
      deallocate( al, bl )
      deallocate( scp, hp )
   endif
!
   IrrSpaceAllocated = .false.
!
   end subroutine deallocateIrrSpace
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isIrrSolTypeValid(useIrrSol) result(y)
!  ===================================================================
   implicit none
!
   character (len=1) :: useIrrSol
!
   logical :: y
!
   if (useIrrSol == 'J' .or. useIrrSol == 'j' .or. useIrrSol == 'H' .or. &
       useIrrSol == 'h') then
      y = .true.
   else
      y = .false.
   endif
!
   end function isIrrSolTypeValid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine solveSingleScattering(spin, site, e, vshift,            &
                                    atom,                             &
                                    isSphSolver, useIrrSol, isCheckWronsk)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
   use MPPModule, only : nbsendMessage, nbrecvMessage
   use MPPModule, only : sendMessage, recvMessage
   use MPPModule, only : waitMessage, setCommunicator, resetCommunicator
!
   use GroupCommModule, only : GlobalSumInGroup, getGroupCommunicator
!
   use MatrixInverseModule, only : MtxInv_GE, MtxInv_LU
!
   use BesselModule, only : SphericalBessel
   use BesselModule, only : SphericalNeumann
!
   use SurfElementsModule, only : getWronskianIntegration
   use SurfElementsModule, only : getNumDiffSurfR
!
   use WriteMatrixModule,  only : writeMatrix
!
   use LdaCorrectionModule, only : checkLdaCorrection, getPotentialCorrection
   use LdaCorrectionModule, only : transformScatteringMatrix
   use LdaCorrectionModule, only : transformWaveFunction
!
   use StepFunctionModule, only : truncate
!
   use PotentialModule, only : getPotential, isPotComponentZero, getPotComponentFlag
   use PotentialModule, only : getTruncatedPotential, getTruncatedPotComponentFlag
!
   use TimerModule, only : getTime
!
   implicit none
!
   character (len=22), parameter :: sname='solveSingleScattering'
   character (len=1), intent(in), optional :: useIrrSol
!
   logical, optional, intent(in) :: isSphSolver, isCheckWronsk
   logical :: isPotSpherical
!
   integer (kind=IntKind), intent(in) :: spin
   integer (kind=IntKind), intent(in) :: site
   integer (kind=IntKind), intent(in), optional :: atom
   integer (kind=IntKind) :: ia, l, lp, mp, mpp, m, np
   integer (kind=IntKind) :: kl, klp, jl, kl1, kl2
   integer (kind=IntKind) :: lpt, idx, idx_sav
   integer (kind=IntKind) :: ir, j, klr, next_pe, prev_pe, comm
   integer (kind=IntKind) :: send_msgid1, send_msgid2, recv_msgid1, recv_msgid2
   integer (kind=IntKind) :: send_msgid3, recv_msgid3
   integer (kind=IntKind) :: kls, klps, klpp, klpps, kl_save, kmax_step, kmax_pot
   integer (kind=IntKind) :: lmax
   integer (kind=IntKind) :: info, n_iter
   integer (kind=IntKind), pointer :: pflag(:)
!
   real (kind=RealKind) :: t0, t1
!
   complex (kind=CmplxKind) :: vcorr, detJ
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in) :: vshift
!
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), pointer :: sin_mat_inv(:,:)
   complex (kind=CmplxKind), pointer :: cos_mat(:,:)
   complex (kind=CmplxKind), pointer :: Omega_mat(:,:)
   complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)
   complex (kind=CmplxKind), pointer :: OmegaHatInv_mat(:,:)
   complex (kind=CmplxKind), pointer :: S_mat(:,:)
   complex (kind=CmplxKind), pointer :: sin_t(:,:), cos_t(:,:)
   complex (kind=CmplxKind), pointer :: dcs_mat(:,:)
   complex (kind=CmplxKind), pointer :: wfr_reg(:,:,:)
   complex (kind=CmplxKind), pointer :: wfr_irr(:,:,:)
   complex (kind=CmplxKind), pointer :: wfdr_reg(:,:,:)
   complex (kind=CmplxKind), pointer :: wfdr_irr(:,:,:)
!
   complex (kind=CmplxKind), pointer :: t_mat(:,:)
   complex (kind=CmplxKind), pointer :: t_mat_inv(:,:)
   complex (kind=CmplxKind), pointer :: jost_mat(:,:)
   complex (kind=CmplxKind), pointer :: jinv_mat(:,:)
   complex (kind=CmplxKind), pointer :: pm(:,:)
   complex (kind=CmplxKind), pointer :: fmem(:,:)
!
   complex (kind=CmplxKind), pointer :: pot_0(:)
   complex (kind=CmplxKind), pointer :: pot_jltr0(:,:)
   complex (kind=CmplxKind), pointer :: p_sml(:), p_cml(:), p_hml(:)
   complex (kind=CmplxKind), pointer :: p_wfr(:,:), p_wfdr(:,:)
!
   complex (kind=CmplxKind) :: cfac, x, rmt, sit
!
   real (kind=RealKind) :: vk_max, vkls, vect(3)
!
!* ===================================================================
!*    Program to solve Schrodinger equation for full potential and 
!*    complex energy using generalized Calogero phase method.
!*
!*    Atomic units employed:
!*        m_e = 1/2, hbar = 1, c = 2/alpha
!*        M(r)= 1 + (energy-v(r,0))*c2inv
!*
!*    Obtain wave functions using :               
!*        regular     :      Pl=r*W and     Ql=r*dW/dr/M
!*        irregular   :      Pl_irr=r*J     Ql_irr=r*dJ/dr/M
!*    where, W and J are regular and irregular solutions of the single
!*    site Schrodinger equation.
!*
!*    Note: assuming one-electron full potential:
!*
!*            ->                ^          ->
!*          v(r ) = v (r) + sum v (r) * Y (r ),   where L=(l,m), l > 0
!*                   0       L   L       L
!*
!*    In this routine :
!*
!*       pot_l   : v (r), l >= 0, m >= 0
!*                  L
!*    Calls:
!*       solveRadEqn
!*       getPhiLr 
!*    Input:
!*       energy  : energy mesh along a contour (e0) 
!*       kappa   : sqrt(e0+e0*e0*c2inv)
!*       pot_l   : potential L-expansion  
!*    Returns:
!*       cos_mat : full potential cosine matrix
!*       sin_mat : full potential cosine matrix
!*       wfr_reg : r*W(r), W(r) is the regular wave function
!*                    *
!*       wfr_irr : r*H (r), H(r) is the irregular wave function
!*
!*    Written by Y. Wang and G. Malcolm Stocks, June 1995
!*    Updated by Y. Wang, May 1996
!*    Updated by Y. Wang, Jun 1996
!*    Updated by Y. Wang, Jan 1998
!*    Updated by Y. Wang, Apr 2001
!* ===================================================================
!
!  -------------------------------------------------------------------
   call initGlobalVariables(spin, site, e, vshift)
!  -------------------------------------------------------------------
!
   if(abs(kappa) < TEN2m6) then
      call ErrorHandler('solveSingleScattering','|kappa| < 1.0d-06',kappa)
   endif
!
!  ===================================================================
!  Set up Riccatti bessel function for use by solveRadEqn
!      bjl(x) = x*jl(x),  bnl(x)=x*nl(x), bhl(x)=x*hl(x)
!      dbjl(x) = x*djl(x),  dbhl(x)=x*dhl(x)
!  ===================================================================
   lmax = lmax_phi+1
   do ir=1,iend
      x=kappa*Grid%r_mesh(ir)
!     ----------------------------------------------------------------
      call SphericalBessel(lmax,x,bjtmp)
      call SphericalNeumann(lmax,x,bntmp)
!     ----------------------------------------------------------------
      do lp=0,lmax
         bjl(ir,lp)=x*bjtmp(lp)
         bnl(ir,lp)=x*bntmp(lp)
         bhl(ir,lp)=bjl(ir,lp)+SQRTm1*bnl(ir,lp)
      enddo
   enddo
   do lp=0,lmax_phi
      cfac=lp/kappa
      do ir=1,iend
         dbjl(ir,lp)=cfac*bjl(ir,lp)/Grid%r_mesh(ir)-bjl(ir,lp+1)
      enddo
      do ir=1,iend
         dbnl(ir,lp)=cfac*bnl(ir,lp)/Grid%r_mesh(ir)-bnl(ir,lp+1)
      enddo
      do ir=1,iend
         dbhl(ir,lp)=cfac*bhl(ir,lp)/Grid%r_mesh(ir)-bhl(ir,lp+1)
      enddo 
   enddo
!
   if ( isPotComponentZero(LocalIndex,1) ) then
      call ErrorHandler('solveSingleScattering','spherical potential is 0')
   endif
!
   if (present(isSphSolver)) then
      isSphericalSolverOn = isSphSolver
   else
      isSphericalSolverOn = .not.FullSolver
   endif
!
   if (present(isCheckWronsk)) then
      CheckWronskian = isCheckWronsk
   else
      CheckWronskian = .false.
   endif
!
   if (present(useIrrSol)) then
      if (.not.isIrrSolTypeValid(useIrrSol) ) then
         call ErrorHandler('solveSingleScattering','Invalid irregular solution boundary condition type', &
                           useIrrSol)
      else if (.not.isIrrSolOn .and. maxval(print_instruction) >= 0 ) then
         write(6,'(a,/)')'In solveSingleScattering, solving irregular solution is enabled.'
      endif
      IrrSolType = useIrrSol
      isIrrSolOn = .true.
   else
      if (isIrrSolOn .and. maxval(print_instruction) >= 0 ) then
         write(6,'(a,/)')'In solveSingleScattering, solving irregular solution is disabled.'
      endif
      isIrrSolOn = .false.
   endif
!
   if (isIrrSolOn .and. .not.IrrSpaceAllocated) then
!     ----------------------------------------------------------------
      call allocateIrrSpace()
!     ----------------------------------------------------------------
   endif
!
!  if ( FullSolver .and.  SSSMethod==-1 ) then
   if (.not.isSphericalSolverOn .and. SSSMethod==-1 ) then
      NumDiffSurfR = getNumDiffSurfR(LocalIndex)
      smpicm => aliasArray2_c(wksl_smpicm,kmax_kkr,kmax_kkr )
!     ----------------------------------------------------------------
      call calSurfJNY()
!     ----------------------------------------------------------------
      rmt = Grid%r_mesh(Grid%jmt)
!     ----------------------------------------------------------------
      call SphericalNeumann(lmax_phi,kappa*rmt,nlrmt(0:lmax_phi),dnlrmt(0:lmax_phi))
!     ----------------------------------------------------------------
   endif
!
   LOOP_ia: do ia = 1, Scatter(site)%NumSpecies
      if (present(atom)) then
         if (ia /= atom .and. atom > 0) then
            cycle LOOP_ia
         endif
      endif
!     ================================================================
      sin_mat => Scatter(site)%Solutions(ia,LocalSpin)%sin_mat ! (1:kmax_int,1:kmax_kkr)
      cos_mat => Scatter(site)%Solutions(ia,LocalSpin)%cos_mat ! (1:kmax_int,1:kmax_kkr)
      jost_mat => Scatter(site)%Solutions(ia,LocalSpin)%jost_mat ! (1:kmax_int,1:kmax_kkr)
      jinv_mat => Scatter(site)%Solutions(ia,LocalSpin)%jinv_mat ! (1:kmax_int,1:kmax_kkr)
      Omega_mat => Scatter(site)%Solutions(ia,LocalSpin)%Omega_mat ! (1:kmax_kkr,1:kmax_kkr)
      OmegaHat_mat => Scatter(site)%Solutions(ia,LocalSpin)%OmegaHat_mat ! (1:kmax_kkr,1:kmax_kkr)
      OmegaHatInv_mat => Scatter(site)%Solutions(ia,LocalSpin)%OmegaHatInv_mat ! (1:kmax_kkr,1:kmax_kkr)
      S_mat => Scatter(site)%Solutions(ia,LocalSpin)%S_mat ! (1:kmax_kkr,1:kmax_kkr)
      wfr_reg => Scatter(site)%Solutions(ia,LocalSpin)%reg_sol ! (1:iend,1:kmax_phi,1:kmax_kkr)
      wfdr_reg => Scatter(site)%Solutions(ia,LocalSpin)%reg_dsol ! (1:iend,1:kmax_phi,1:kmax_kkr)
      t_mat => Scatter(site)%Solutions(ia,LocalSpin)%t_mat ! (1:kmax_kkr,1:kmax_kkr)
      t_mat_inv => Scatter(site)%Solutions(ia,LocalSpin)%t_mat_inv ! (1:kmax_kkr,1:kmax_kkr)
!     ================================================================
!     if ( FullSolver ) then
      if (.not.isSphericalSolverOn) then
         pot_jl => getPotential(LocalIndex,ia,spin_index)
         pflag => getPotComponentFlag(LocalIndex)
         flags_jl(1:jmax_pot) = pflag(1:jmax_pot)
!        =============================================================
!        calculate v0r, the spherical component of potential pot_l....
!        =============================================================
         pot_0 => pot_jl(1:numrs_cs,1)
         do ir=1,numrs_cs
            v0r(ir)=(Y0*pot_0(ir)+PotShift)*Grid%r_mesh(ir)
            cm0(ir)=CONE+c2inv*(energy-pot_0(ir)*Y0-PotShift)
         enddo
         do ir=numrs_cs+1,iend
            v0r(ir)=CZERO
            cm0(ir)=CONE
         enddo
         if ( SSSMethod==1 .or. SSSMethod==2 ) then
            pot_trunc_jl => getTruncatedPotential(LocalIndex,ia,spin_index)
            pflag => getTruncatedPotComponentFlag(LocalIndex)
            flags_trunc_jl(1:jmax_trunc) = pflag(1:jmax_trunc)
            pot_0 => pot_trunc_jl(1:numrs_trunc,1)
            do ir = 1,numrs_trunc
               v0r(numrs_mt+ir-1) = (Y0*pot_0(ir)+PotShift)*Grid%r_mesh(numrs_mt+ir-1)
               cm0(numrs_mt+ir-1) = CONE+c2inv*(energy-pot_0(ir)*Y0-PotShift)
            enddo
            nullify(pot_0)
         endif
      else
         flags_jl(1) = 1
         pot_0 => getPotential(LocalIndex,ia,spin_index,1)
!
!        =============================================================
!        calculate v0r, the spherical component of potential pot_l.......
!
!        calculate M (r) and store it in cm0, where
!                   0
!
!        M (r) = 1 + c2inv * ( e - v (r) )
!         0                         0
!        =============================================================
         do ir=1,numrs_cs
            v0r(ir)=(Y0*pot_0(ir)+PotShift)*Grid%r_mesh(ir)
            cm0(ir)=CONE+c2inv*(energy-pot_0(ir)*Y0-PotShift)
         enddo
         do ir=numrs_cs+1,iend
            v0r(ir)=CZERO
            cm0(ir)=CONE
         enddo
         nullify(pot_0)
      endif
!
!     ================================================================
!     Uncomment the following three lines allows to test the code for
!     the empty cell case
!     ================================================================
!     v0r = CZERO
!     cm0 = CONE+c2inv*(energy-PotShift)
!     Scatter(site)%AtomicNumber(ia) = 0
!     ================================================================
!
!     ================================================================
!     check if the spherical component of the potential is zero beyond the
!     muffin-tin radius
!     ================================================================
      isSphPotZeroOutsideRmt = .true.
      LOOP_ir: do ir=numrs_mt+1,iend
         if (abs(v0r(ir)) > TEN2m6) then
            isSphPotZeroOutsideRmt = .false.
            exit LOOP_ir
         endif
      enddo LOOP_ir
!
!     ================================================================
!     Turn on spherical solver if the potential is spherical
!     ================================================================
      if (.not.present(isSphSolver) .and. .not.isSphericalSolverOn) then
         isPotSpherical = .true.
         LOOP_jl0: do jl = 2, jmax_pot
            if (flags_jl(jl) /= 0) then
               do ir=1,numrs_cs
                  if (abs(pot_jl(ir,jl)) > TEN2m6) then
                     isPotSpherical = .false.
                     exit LOOP_jl0
                  endif
               enddo
            endif
         enddo LOOP_jl0
         isSphericalSolverOn = isPotSpherical
      endif
!
!     if ( FullSolver ) then
      if (.not.isSphericalSolverOn) then
         if (.not.isSphPotZeroOutsideRmt) then
            isPotZeroOutsideRmt = .false.
         else
            isPotZeroOutsideRmt = .true.
            LOOP_jl: do jl = 1, jmax_pot
               if (flags_jl(jl) /= 0) then
                  do ir=1,Grid%jmt
                     if (abs(pot_jl(ir,jl)) > TEN2m6) then
                        isSphericalSolverOn = .false.
                        exit LOOP_jl
                     endif
                  enddo
                  do ir=Grid%jmt+1,numrs_cs
                     if (abs(pot_jl(ir,jl)) > TEN2m6) then
                        isPotZeroOutsideRmt = .false.
                        exit LOOP_jl
                     endif
                  enddo
               endif
            enddo LOOP_jl
         endif
#ifdef TMP_ACCEL
!        -------------------------------------------------------------
         call push_calphilr_data(spin,ia,energy,kappa,iend_max,lmax_phi,bjl,dbjl,bnl,dbnl,cm0)
!        -------------------------------------------------------------
#endif
      else
         isPotZeroOutsideRmt = isSphPotZeroOutsideRmt
      endif
!
!     ================================================================
!     Check if the Orbital dependent LDA correction is needed
!     ================================================================
      if (checkLdaCorrection(LocalIndex,ia)) then
         v0r_save(1:Grid%jmt) = v0r(1:Grid%jmt)
      endif
!
!     ================================================================
!     obtain Sine and Cosine matrices (L,Lp) with one Lp at a time.
!         l  = 0, 1, ..., lmax_phi
!         lp = 0, 1, ..., lmax_kkr
!     ================================================================
      lpt=-1
      idx_sav = 0
      fmem => aliasArray2_c(TmpSpace,iend,4)
      do kl = 1,kmax_kkr
         idx = Scatter(site)%lm_index_sph(kl,ia)
         l=lofk(kl)
!
         if (checkLdaCorrection(site,ia,l)) then
            vcorr = getPotentialCorrection(LocalIndex,ia,spin,l,mofk(kl))
            v0r(1:Grid%jmt) = v0r(1:Grid%jmt) + vcorr*Grid%r_mesh(1:Grid%jmt)
         endif
!
         if ( idx > idx_sav ) then
            pl_regwrk(1:iend,l) = CZERO
            ql_regwrk(1:iend,l) = CZERO
            pl_reg => pl_regwrk(1:iend,l)
            ql_reg => ql_regwrk(1:iend,l)
!           ==========================================================
!           solving scalar-relativistic or non-relativistic equation for
!           spherical part of the potential.
!           ----------------------------------------------------------
            call solveRadEqn4Reg(Scatter(site)%AtomicNumber(ia),l,fmem)
!           ----------------------------------------------------------
!           write(6,'(a,i3,2x,2d15.7,2x,2d15.7)') 'sl, cl = ',l,s_sph_mt(l),c_sph_mt(l)
!           write(6,'(a,i3,2x,2d15.7,2x,2d15.7)') 'pl, ql = ',l,pl_reg(numrs_mt), ql_reg(numrs_mt)
            idx_sav = idx
         endif
!
         if (checkLdaCorrection(LocalIndex,ia,l)) then
            v0r(1:Grid%jmt) = v0r_save(1:Grid%jmt)
         endif
      enddo
!
      wfr_reg = CZERO
      wfdr_reg = CZERO
!
!     if ( SphericalSolver .or. isPotZeroOutsideRmt ) then
      if (isSphericalSolverOn) then
         do kl = 1,kmax_kkr
            l=lofk(kl)
!           ==========================================================
!           Solve coupled integral equations, and calculate sin_mat and 
!           cos_mat using surface integration technique.
!
!           Muffin-tin or ASA potential assign :
!           cos_mat and sin_mat contribution from spherical part potential
!           inside the muffin-tin sphere or Wigner-Seitz sphere
!           ==========================================================
            sin_mat(1:kmax_int,kl)=CZERO
            cos_mat(1:kmax_int,kl)=CZERO
            t_mat(1:kmax_int,kl)=CZERO
            t_mat_inv(1:kmax_int,kl)=CZERO
            jost_mat(1:kmax_int,kl)=CZERO
            jinv_mat(1:kmax_int,kl)=CZERO
            Omega_mat(1:kmax_kkr,kl)=CZERO
            OmegaHat_mat(1:kmax_kkr,kl)=CZERO
            OmegaHatInv_mat(1:kmax_kkr,kl)=CZERO
            S_mat(1:kmax_kkr,kl)=CZERO
            if ( isSphPotZeroOutsideRmt ) then !why do we have this if-statement???
               cos_mat(kl,kl)=c_sph_mt(l)   ! spherical part
               sin_mat(kl,kl)=s_sph_mt(l)   ! spherical part
            else
               cos_mat(kl,kl)=c_mtrx(l)   ! spherical part
               sin_mat(kl,kl)=s_mtrx(l)   ! spherical part
            endif
!
            jost_mat(kl,kl)=sqrtm1*sin_mat(kl,kl)-cos_mat(kl,kl)
            jinv_mat(kl,kl)=CONE/jost_mat(kl,kl)
            t_mat(kl,kl)=sin_mat(kl,kl)*jinv_mat(kl,kl)/kappa
            t_mat_inv(kl,kl)=kappa*(sqrtm1-cos_mat(kl,kl)/sin_mat(kl,kl))
            Omega_mat(kl,kl) = CONE/(sin_mat(kl,kl)**2+cos_mat(kl,kl)**2)
            OmegaHatInv_mat(kl,kl) = sin_mat(kl,kl)*jost_mat(kl,kl)
            OmegaHat_mat(kl,kl) = CONE/OmegaHatInv_mat(kl,kl)
            S_mat(kl,kl) = CONE-TWO*SQRTm1*sin_mat(kl,kl)*jinv_mat(kl,kl)
!           ----------------------------------------------------------
            call zcopy(iend,pl_regwrk(1:iend,l),1,wfr_reg(1:iend,kl,kl),1)
            call zcopy(iend,ql_regwrk(1:iend,l),1,wfdr_reg(1:iend,kl,kl),1)
!           ----------------------------------------------------------
            Scatter(site)%sol_flags(kl,kl) = 1
         enddo
      else
         WronskianIndx = 0
         TmpSpace = CZERO
         fmem => aliasArray2_c(TmpSpace,iend,5)
         t0 = getTime()
!
!        =============================================================
!        calculate sum [gaunt_factor * pot * step_function]
!        This routine call is independent of energy, so could be moved out
!        of the energy loop.
!        =============================================================
         if ( SSSMethod == 0 ) then
!           ----------------------------------------------------------
            call setupSSSM0(spin,LocalIndex,ia)
!           ----------------------------------------------------------
         endif
!
#ifdef TMP_ACCEL
!        -------------------------------------------------------------
         call accel_calphilr(numrs_mt,numrs_cs,iend,kmax_kkr,kmax_phi,   &
                             pl_regwrk,ql_regwrk,wfr_reg,wfdr_reg,       &
                             Scatter(site)%sol_flags)
!        -------------------------------------------------------------
#else
         if (NumPEsInGroup > 1) then
!           ----------------------------------------------------------
            comm = getGroupCommunicator(kGID)
            call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!           ----------------------------------------------------------
         endif
         np = mod(kmax_kkr,NumPEsInGroup)
         do kl = MyPEinGroup+1, kmax_kkr-np, NumPEsInGroup
            l = lofk(kl)
            pl_reg => pl_regwrk(1:iend,l)
            ql_reg => ql_regwrk(1:iend,l)
!
            plhat_reg=>wfr_reg(1:iend,1:kmax_phi,kl)
            qlhat_reg=>wfdr_reg(1:iend,1:kmax_phi,kl)
            do j = 1, NumPEsInGroup
               if (j < NumPEsInGroup) then
                  prev_pe = mod(MyPEInGroup-j+NumPEsInGroup,NumPEsInGroup)
                  klr = kl-MyPEinGroup + prev_pe
                  p_wfr => wfr_reg(1:iend,1:kmax_phi,klr)
                  p_wfdr => wfdr_reg(1:iend,1:kmax_phi,klr)
!                 ----------------------------------------------------
                  recv_msgid1 = nbrecvMessage(p_wfr,iend,kmax_phi,klr*1000+1,prev_pe)
                  recv_msgid2 = nbrecvMessage(p_wfdr,iend,kmax_phi,klr*1000+2,prev_pe)
!                 ----------------------------------------------------
               endif
!
               if (j == 1) then
!                 ====================================================
!                 solving coupled integral equation for plhat_reg
!                 and then calculate qlhat_reg
!                 ====================================================
                  t1 = getTime()
!                 ----------------------------------------------------
!                 call solveIntEqn(ia,kl,fmem,.true.,.false.)
                  call calPhiLr(kl,fmem,4,pl_reg,ql_reg)
!                 ----------------------------------------------------
#ifdef TIMING
                  write(6,'(a,f10.5)')'In solveSS, time for each calPhiLr = ',getTime()-t1
#endif
               endif
               if (j < NumPEsInGroup) then
                  next_pe = mod(MyPEInGroup+j,NumPEsInGroup)
!                 ----------------------------------------------------
                  send_msgid1 = nbsendMessage(plhat_reg,iend,kmax_phi,kl*1000+1,next_pe)
                  send_msgid2 = nbsendMessage(qlhat_reg,iend,kmax_phi,kl*1000+2,next_pe)
                  call waitMessage(recv_msgid1)
                  call waitMessage(recv_msgid2)
                  call waitMessage(send_msgid1)
                  call waitMessage(send_msgid2)
!                 ----------------------------------------------------
               endif
            enddo
         enddo
         if (NumPEsInGroup > 1) then
            nullify(p_wfr, p_wfdr)
!           ----------------------------------------------------------
            call resetCommunicator()
!           ----------------------------------------------------------
         endif
         do kl = kmax_kkr-np+1, kmax_kkr
            l = lofk(kl)
            pl_reg => pl_regwrk(1:iend,l)
            ql_reg => ql_regwrk(1:iend,l)
!
            plhat_reg=>wfr_reg(1:iend,1:kmax_phi,kl)
            qlhat_reg=>wfdr_reg(1:iend,1:kmax_phi,kl)
!
!           ==========================================================
!           solving coupled integral equation for plhat_reg
!           and then calculate qlhat_reg
!           ----------------------------------------------------------
!           call solveIntEqn(ia,kl,fmem,.true.,.false.)
            call calPhiLr(kl,fmem,4,pl_reg,ql_reg)
!           ----------------------------------------------------------
         enddo
#endif
#ifdef TIMING
         write(6,'(a,f10.5)')'In solveSS, cumulative time for calPhiLr = ',getTime()-t0
#endif
!
         do kl = 1,kmax_kkr
            plhat_reg=>wfr_reg(1:iend,1:kmax_phi,kl)
            qlhat_reg=>wfdr_reg(1:iend,1:kmax_phi,kl)
!
            p_sml => sin_mat(1:kmax_int,kl)
            p_cml => cos_mat(1:kmax_int,kl)
            if (SSSMethod == -1) then
!              -------------------------------------------------------
               call calSCMatrix( kl, p_sml, p_cml )
!              -------------------------------------------------------
            else if (SSSMethod == 0 .or. SSSMethod == 1) then
               TmpSpace = CZERO
!              -------------------------------------------------------
               call calSCMatrixVolInt( ia, kl, p_sml, p_cml, fmem )
!              -------------------------------------------------------
            else if (SSSMethod == 2) then
               do klp = 1, kmax_int
                  lp = lofk(klp)
                  sin_mat(klp,kl) = bjl(numrs_cs,lp)*qlhat_reg(numrs_cs,klp)*cm0(numrs_cs) &
                                   -kappa*plhat_reg(numrs_cs,klp)*dbjl(numrs_cs,lp)
                  cos_mat(klp,kl) = bnl(numrs_cs,lp)*qlhat_reg(numrs_cs,klp)*cm0(numrs_cs) &
                                   -kappa*plhat_reg(numrs_cs,klp)*dbnl(numrs_cs,lp)
               enddo
            else
!              -------------------------------------------------------
               call ErrorHandler(sname,'Not knowing how to calculate the sine/cosine matrices')
!              -------------------------------------------------------
            endif
            do klp = 1, kmax_int
               jost_mat(klp,kl) = SQRTm1*sin_mat(klp,kl) - cos_mat(klp,kl)
            enddo
!           do klp = 1, kmax_int
!              lp = lofk(klp)
!              write(6,'(a,2i4,2d16.8,2x,2d16.8)')'phi = ',klp,kl,              &
!                                             bnl(numrs_cs,lp)*sin_mat(klp,kl)- &
!                                             bjl(numrs_cs,lp)*cos_mat(klp,kl), &
!                                             kappa*plhat_reg(numrs_cs,klp)
!           enddo
         enddo
!        =============================================================
!        The following code is for testing purposes
!        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        do kl = 1,kmax_kkr
!           call checkWaveFuncAtR(kl,numrs_cs,ONE,ONE,ONE)
!        enddo
!
!        sin_t => aliasArray2_c(wks_mtx1,kmax_int,kmax_kkr)
!        do kl = 1,kmax_kkr
!           m = mofk(kl)
!           do klp = 1, kmax_int
!              mp = mofk(klp)
!              sin_t(klp,kl) = m1m(m+mp)*sin_mat(klp-2*mp,kl-2*m)
!           enddo
!        enddo
!        -------------------------------------------------------------
!        call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_int, CONE,    &
!                    sin_t, kmax_int, jost_mat, kmax_int, CZERO, OmegaHat_mat, kmax_kkr)
!        -------------------------------------------------------------
!        call writeMatrix('S^{T*}*Jost matrix',OmegaHat_mat,kmax_kkr,kmax_kkr,TEN2m6)
!        stop 'Ok'
!        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        The end of testing
!        =============================================================
!
!        =============================================================
!        calculate the OmegaHat = [iS^{*T}xS-S^{*T}*C]^{-1}
!                      OmegaHatInv = iS^{*T}xS-S^{*T}*C
!        =============================================================
         sin_t => aliasArray2_c(wks_mtx1,kmax_int,kmax_kkr)
         do kl = 1,kmax_kkr
            m = mofk(kl)
            do klp = 1,kmax_int
               mp = mofk(klp)
               sin_t(klp,kl) = m1m(m+mp)*sin_mat(klp-2*mp,kl-2*m)
            enddo
         enddo
!
         if (.false.) then
            cos_t  => aliasArray2_c(wks_mtx2,kmax_int,kmax_kkr)
            do kl = 1,kmax_kkr
               m = mofk(kl)
               do klp = 1,kmax_int
                  mp = mofk(klp)
                  cos_t(klp,kl) = m1m(m+mp)*cos_mat(klp-2*mp,kl-2*m)
               enddo
            enddo
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_int, CONE, &
                        sin_t, kmax_int, sin_mat, kmax_int, CZERO, OmegaHatInv_mat, kmax_kkr)
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_int, CONE, &
                        cos_t, kmax_int, cos_mat, kmax_int, CONE, OmegaHatInv_mat, kmax_kkr)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_int, CONE, &
                        sin_t, kmax_int, jost_mat, kmax_int, CZERO, OmegaHatInv_mat, kmax_kkr)
!           ----------------------------------------------------------
         endif
!        call writeMatrix('Sine Matrix',sin_mat,kmax_int,kmax_kkr,5.0d0*TEN2m8)
!        call writeMatrix('Cosine Matrix',cos_mat,kmax_int,kmax_kkr,5.0d0*TEN2m8)
!        call writeMatrix('OmegaHatInv Matrix',OmegaHatInv_mat,kmax_kkr,kmax_kkr,5.0d0*TEN2m8)
!        -------------------------------------------------------------
         OmegaHat_mat = OmegaHatInv_mat
         call MtxInv_LU(OmegaHat_mat,kmax_kkr)
!        -------------------------------------------------------------
!
!        =============================================================
!        It might be better to use jinv_mat = OmegaHat*S^{*T}
!        =============================================================
! !!     jinv_mat = jost_mat
! !!     call MtxInv_LU(jinv_mat,kmax_kkr)
!        -------------------------------------------------------------
         call zgemm( 'n', 't', kmax_kkr, kmax_int, kmax_kkr, CONE,    &
                     OmegaHat_mat, kmax_kkr, sin_t, kmax_int, CZERO, jinv_mat, kmax_kkr)
!        -------------------------------------------------------------
!
!        =============================================================
!        For testing purposes, calculate C^{T*}*S - S^{T*}*C
!        =============================================================
!        dcs_mat => aliasArray2_c( wks_jinv, kmax_kkr, kmax_kkr)
!        -------------------------------------------------------------
!        call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_int, CONE,    &
!                    cos_t, kmax_int, sin_mat, kmax_int, CZERO, dcs_mat, kmax_kkr)
!        -------------------------------------------------------------
!        call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_int, -CONE,   &
!                    sin_t, kmax_int, cos_mat, kmax_int, CONE, dcs_mat, kmax_kkr)
!        -------------------------------------------------------------
!        call writeMatrix('C^{T*}*S-S^{T*}*C Matrix',dcs_mat,kmax_kkr,kmax_kkr,TEN2m8)
!        -------------------------------------------------------------
!
!        =============================================================
!        calculate the t-matrix
!        =============================================================
         if (kmax_kkr == kmax_int) then
!           ==========================================================
!           Use t = S*[i*S - C]^{-1}/kappa
!           ----------------------------------------------------------
            call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, CONE, &
                       sin_mat, kmax_kkr, jinv_mat, kmax_kkr, CZERO, t_mat, kmax_kkr)
!           ----------------------------------------------------------
         else
!           ==========================================================
!           Use t = S*OmegaHat*S^{T*}/kappa. 
!           ==========================================================
            cos_t  => aliasArray2_c(wks_mtx2,kmax_int,kmax_kkr) ! use cos_t as temp space
!           ----------------------------------------------------------
            call zgemm( 'n', 'n', kmax_int, kmax_kkr, kmax_kkr, -CONE,&
                       sin_mat, kmax_int, OmegaHat_mat, kmax_kkr, CZERO, cos_t, kmax_int)
!           ----------------------------------------------------------
            call zgemm( 'n', 't', kmax_int, kmax_int, kmax_kkr, CONE, &
	               cos_t, kmax_int, sin_t, kmax_int, CZERO, t_mat, kmax_int)
!           ----------------------------------------------------------
         endif
!
!        =============================================================
!        Use t^{-1} = kappa*[i*S - C]*S^{-1} = kappa*Jost*S^{-1}
!        =============================================================
         t_mat_inv = jost_mat
         sin_mat_inv => S_mat  ! Use S_mat as a swapping space
         sin_mat_inv = sin_mat
         call MtxInv_LU(sin_mat_inv,kmax_kkr)
!        -------------------------------------------------------------
         call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, kappa,   &
                    jost_mat, kmax_kkr, sin_mat_inv, kmax_kkr, CZERO, &
                    t_mat_inv, kmax_kkr)
!        -------------------------------------------------------------
         nullify(sin_mat_inv)
!
!        =============================================================
!        calculate the S-matrix = 1 - 2*i*kappa*t-matrix
!        =============================================================
         S_mat = CZERO
         do kl = 1, kmax_kkr
            S_mat(kl,kl) = CONE
         enddo
         cfac = TWO*SQRTm1
         S_mat = S_mat - cfac*t_mat
!
         t_mat = t_mat/kappa
      endif
!
!     call writeMatrix('Sine Matrix',sin_mat,kmax_int,kmax_kkr,HALF*TEN2m6)
!     call writeMatrix('Cosine Matrix',cos_mat,kmax_int,kmax_kkr,HALF*TEN2m6)
!     call writeMatrix('Omega Matrix',OmegaHat_mat,kmax_kkr,kmax_kkr,HALF*TEN2m6)
!     call writeMatrix('T Matrix',t_mat,kmax_int,kmax_int,TEN2m6)
!     call writeMatrix('S Matrix',S_mat,kmax_int,kmax_int,TEN2m6)
!
      nullify(sin_t, cos_t, dcs_mat, t_mat_inv)
!
      if (checkLdaCorrection(site,ia)) then
         call transformScatteringMatrix(site,ia,spin,Scatter(LocalIndex)%Solutions(ia,LocalSpin)%sin_mat)
         call transformScatteringMatrix(site,ia,spin,Scatter(LocalIndex)%Solutions(ia,LocalSpin)%cos_mat)
         call transformScatteringMatrix(site,ia,spin,Scatter(LocalIndex)%Solutions(ia,LocalSpin)%t_mat)
         call transformScatteringMatrix(site,ia,spin,Scatter(LocalIndex)%Solutions(ia,LocalSpin)%t_mat_inv)
         call transformWaveFunction(site,ia,spin,iend,kmax_phi,wfr_reg)
      endif
!
!     ================================================================
!     calculate the irregular solutions.
!     ================================================================
      if ( isIrrSolOn ) then
!        pm => aliasArray2_c(wks_mtx1,kmax_int,kmax_kkr)
!        pm = SQRTm1*sin_mat+cos_mat
         pm => sin_mat
         fmem => aliasArray2_c(TmpSpace,iend,4)
!        -------------------------------------------------------------
!        call computeIrrSol('H',ia,fmem,pm)
         call computeIrrSol(useIrrSol,ia,fmem)
!        -------------------------------------------------------------
         if (CheckWronskian) then
            wfr_irr => Scatter(site)%Solutions(ia,LocalSpin)%irr_sol
            wfdr_irr => Scatter(site)%Solutions(ia,LocalSpin)%irr_dsol
            if (isSphericalSolverOn) then
               do l = 0, lmax_kkr
                  kl1 = (l+1)**2-l
!                 ----------------------------------------------------
                  call printWronskian( kl1, kl1, wfr_reg, wfdr_reg, wfr_irr, &
                                       wfdr_irr, 1, numrs_mt-20, 50 )
                  call printWronskian( kl1, kl1, wfr_reg, wfdr_reg, wfr_irr, &
                                       wfdr_irr, numrs_mt-250, numrs_mt-100, 10)
                  call printWronskian( kl1, kl1, wfr_reg, wfdr_reg, wfr_irr, &
                                       wfdr_irr, numrs_mt-20, numrs_mt+10, 5 )
!                 ----------------------------------------------------
               enddo
            else
               do kl2 = 1,kmax_kkr
                  do kl1 = 1,kmax_kkr
!                    -------------------------------------------------
                     call printWronskian( kl1, kl2, wfr_reg, wfdr_reg, wfr_irr, &
                                          wfdr_irr, 1, numrs_mt-20, 50 )
                     call printWronskian( kl1, kl2, wfr_reg, wfdr_reg, wfr_irr, &
                                          wfdr_irr, numrs_mt-20, numrs_cs+10, 5 )
!                    -------------------------------------------------
                  enddo
               enddo
            endif
            nullify(wfdr_irr,wfr_irr)
         endif
      endif
   enddo LOOP_ia
!
#ifdef TMP_ACCEL
!  if ( FullSolver .and. .not.isSphPotZeroOutsideRmt ) then
   if (.not.isSphericalSolverOn) then
!     ----------------------------------------------------------------
!     call delete_calphilr_data()
!     ----------------------------------------------------------------
   endif 
#endif
!
   nullify( plhat_reg, plhat_irr )
   nullify( sin_mat, cos_mat, wfr_reg, wfdr_reg )
!
   if(stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine solveSingleScattering
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupSSSM0(spin,site,atom)
!  ===================================================================
   use PotentialModule, only : getPotential, getPotComponentFlag
   use LdaCorrectionModule, only : checkLdaCorrection
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, site, atom
   integer (kind=IntKind) :: kmax_step, kmax_pot, npout, loc_spin
   integer (kind=IntKind) :: klv, mlv, jlv, kls, mls, jls, ir, nj3_i3, i3, kl3
   integer (kind=IntKind), pointer :: kj3_i3(:)
   integer (kind=IntKind), pointer :: pflag(:)
!
   real (kind=RealKind), pointer :: cgnt_i3(:)
!
   complex (kind=CmplxKind) :: vcorr
   complex (kind=CmplxKind), pointer :: tmp_pot(:,:), tmp_step(:,:), pot_l(:,:)
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   loc_spin = min(spin, NumSpins)
!
   kmax_step = (Scatter(site)%lmax_step+1)**2
   kmax_pot = (Scatter(site)%lmax_pot+1)**2
   npout = Scatter(site)%numrs_cs - Scatter(site)%numrs_mt + 1
   p1c => wks_scvphi(1:kmax_pot*npout,1)
   tmp_pot => aliasArray2_c( p1c,npout,kmax_pot )
   p1c => wks_scvphi(1:kmax_step*npout,2)
   tmp_step => aliasArray2_c( p1c,npout,kmax_step )
!
   if (checkLdaCorrection(site,atom)) then
      call ErrorHandler('setupSSSM0','This piece of code needs to be put within a Lp loop')
   else
      vcorr = CZERO
   endif
!
   pot_l => getPotential(site,atom,spin)
   pflag => getPotComponentFlag(site)
   tmp_pot = CZERO
   do klv = 1, kmax_pot
      mlv = mofk(klv)
      jlv = jofk(klv)
      if (pflag(jlv) /= 0) then
         if (klv==1) then
            do ir = 1,npout
               tmp_pot(ir,klv) = pot_l(numrs_mt+ir-1,jlv)+(vcorr+PotShift)/Y0
            enddo                     
         else if (mlv>=0) then
            do ir = 1,npout
               tmp_pot(ir,klv) = pot_l(numrs_mt+ir-1,jlv)
            enddo
         else
            do ir = 1,npout
               tmp_pot(ir,klv) = m1m(mlv)*conjg(pot_l(numrs_mt+ir-1,jlv))
            enddo
         endif
      endif
   enddo
   do kls = 1, kmax_step
      mls = mofk(kls)
      jls = jofk(kls)
      if (mls>=0) then
         do ir = 1,npout
            tmp_step(ir,kls) = Scatter(site)%step_rad(ir,jls)
         enddo
      else
         do ir = 1,npout
            tmp_step(ir,kls) = m1m(mls)*conjg(Scatter(site)%step_rad(ir,jls))
         enddo
      endif
   enddo
!
   cpotstep = CZERO
   do klv = 1, kmax_pot
      jlv = jofk(klv)
      if (pflag(jlv) /= 0) then
         do kls = 1, kmax_step
            nj3_i3 = nj3(kls,klv)
            kj3_i3 => kj3(1:nj3_i3,kls,klv)
            cgnt_i3 => cgnt(1:nj3_i3,kls,klv)
            do i3 = 1, nj3_i3
               kl3 = kj3_i3(i3)
               do ir = 1, npout
                  cpotstep(ir,kl3) = cpotstep(ir,kl3) + cgnt_i3(i3)*tmp_pot(ir,klv)*conjg(tmp_step(ir,kls))
               enddo
            enddo
         enddo
      endif
   enddo
!
   nullify(tmp_step, tmp_pot)
!
   end subroutine setupSSSM0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine solveRadEqn4Reg(atomic_number,l,fmem)
!  ===================================================================
   use InterpolationModule, only : PolyInterp, FitInterp
!
   use IntegrationModule, only : calIntegration
!
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
!
   implicit none
!
   character (len= 15), parameter :: sname='solveRadEqn4Reg'
!
   integer (kind=IntKind), intent(in) :: atomic_number, l
!
   integer (kind=IntKind) :: lp1
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: j_inter
   integer (kind=IntKind) :: npower
   integer (kind=IntKind) :: ibeg
!
   integer (kind=1) :: SolutionType
   integer (kind=IntKind), parameter :: n_inter=5
!
   real (kind=RealKind) :: Zi
   real (kind=RealKind) :: beta
   real (kind=RealKind) :: betap1
   real (kind=RealKind) :: fz2oc2
   real (kind=RealKind) :: tzoc2
   real (kind=RealKind) :: llp1
   real (kind=RealKind) :: rp
   real (kind=RealKind) :: xmt
   real (kind=RealKind) :: z_err_est

   real (kind=RealKind), pointer :: xp(:)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), intent(inout), target :: fmem(iend,4)
!
   complex (kind=CmplxKind), pointer :: yp(:)

   complex (kind=CmplxKind) :: v0, v1
   complex (kind=CmplxKind), pointer :: cl(:)
   complex (kind=CmplxKind), pointer :: sl(:)
   complex (kind=CmplxKind), pointer :: dcl(:)
   complex (kind=CmplxKind), pointer :: dsl(:)
!
   complex (kind=CmplxKind) :: norm
!
   complex (kind=CmplxKind) :: clp_reg(0:3)
   complex (kind=CmplxKind) :: slp_reg(0:3)
   complex (kind=CmplxKind) :: dclp_reg(0:3)
   complex (kind=CmplxKind) :: dslp_reg(0:3)
!
   complex (kind=CmplxKind) :: eoc2p1
   complex (kind=CmplxKind) :: e2oc2
   complex (kind=CmplxKind) :: em0
   complex (kind=CmplxKind) :: a0
   complex (kind=CmplxKind) :: a1p
   complex (kind=CmplxKind) :: a1
   complex (kind=CmplxKind) :: a2
   complex (kind=CmplxKind) :: a2p
   complex (kind=CmplxKind) :: a3
   complex (kind=CmplxKind) :: a3p
   complex (kind=CmplxKind) :: pr
   complex (kind=CmplxKind) :: ddcl, ddsl, dgpl, dfpl, b1, b2, em, emr
   complex (kind=CMplxKind) :: rvr, drvr
   complex (kind=CmplxKind) :: c_sphere_mt
   complex (kind=CmplxKind) :: s_sphere_mt
   complex (kind=CmplxKind) :: c_matrix
   complex (kind=CmplxKind) :: s_matrix
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
! *===================================================================
! *   this subroutine solves the scalar relativistic equation.........
! *   Modified by Y.W. in July 4th, 1994..............................
! *
! *   Input...........................................................
! *       l          : current l value................................
! *       energy     : electon energy.................................
! *       kappa      : Relativistic(non-relativistic) momentum........
! *
! *   Output..........................................................
! *       c_matrix   : Cosine matrix  @ r_mesh(jend) for sph. potential
! *       s_matrix   : Sine   matrix  @ r_mesh(jend) for sph. potential
! *       c_sphere_mt: Cosine matrix  @ r_mesh(jmt) used in full potential
! *       s_sphere_mt: Sine   matrix  @ r_mesh(jmt) used in full potential
! *
! *   Grid information................................................
! *       jmt        : index of muffin-tin radius.....................
! *       jend       : = jmt for muffin-tin potential.................
! *                  : = jws for ASA-potential........................
! *                  : = jcirc for full-potential.....................
! *       iend       : end of grid, a few points (> 5) beyond jend....
! *       r_mesh(i)  : radial mesh....................................
! *       x_mesh(i)  : Double uniform spacing log grid [x=ln(r)]......
! *       hin        : x_mesh spacing for x_mesh(1)  < x <x_mesh(jmt).
! *       hout       : x_mesh spacing for x_mesh(jmt)< x <x_mesh(iend)
! *
! *   Potential.......................................................
! *       v0r(i)     : r * v (r) * Y    on logrid.....................
! *                         0       0,0
! *   Solutions (in SSteringModule)...........................
! *       pl_reg  : regular solution   :: ql_reg = d(pl_reg)/dr-pl_reg/r.
! *       pl_irr  : irregular solution :: ql_irr = d(pl_irr)/dr-pl_irr/r.
! *
! *===================================================================
!
!  ===================================================================
!  cannot solve for energy equal to 0.0
!  ===================================================================
   if( abs(energy) < TEN2m10 ) then
      call ErrorHandler(sname,'energy is too small',energy)
   endif
!
!  ===================================================================
   cl  => fmem(1:iend,1)
   sl  => fmem(1:iend,2)
   dcl => fmem(1:iend,3)
   dsl => fmem(1:iend,4)
   xmt=Grid%x_mesh(Grid%jmt)
   r_mesh => Grid%r_mesh(1:iend)
   Zi = atomic_number
!  ===================================================================
!
!  ===================================================================
!  get parameters: 1/c, 1/c**2, e**2/c**2, 1+e/c**2, 1+(e-v0)/c**2, 
!  etc.............................................................
!  ===================================================================
!
   lp1=l+1
   llp1=l*lp1
   fz2oc2=FOUR*Zi*Zi*c2inv
   tzoc2=TWO*Zi*c2inv
!  ... BEGIN
!  ===================================================================
!  Note: the following line is the original code as of 09/28/03.......
!  ===================================================================
!   v0=(v0r(1)+TWO*Zi)/Grid%r_mesh(1)
!
!  ===================================================================
!  Note: the following lines are copied over from gms on 09/29/03......
!  ===================================================================
!   v1 = ( Grid%r_mesh(2)*( v0r(1)+two*AtomicNumber(LocalIndex) ) -     &
!          Grid%r_mesh(1)*( v0r(2)+two*AtomicNumber(LocalIndex) ) ) /   &
!        ( two*lp1+1 )
!   v0 = v1*( Grid%r_mesh(1)+Grid%r_mesh(2) ) + ( v0r(1)-v0r(2) )/      &
!        ( Grid%r_mesh(1)-Grid%r_mesh(2) )
!  ===================================================================
!  ... END
   if ( Zi >= 0 ) then
      v1=(r_mesh(2)*(v0r(1)+2*Zi)-r_mesh(1)*(v0r(2)+2*Zi))/ &
         (r_mesh(1)*r_mesh(2)*(r_mesh(1)-r_mesh(2)))
      v0=v1*(r_mesh(1)+r_mesh(2))+(v0r(1)-v0r(2))/(r_mesh(1)-r_mesh(2))
   else
      v1 = ( r_mesh(2)*(v0r(1)+two*Zi) - r_mesh(1)*(v0r(2)+two*Zi) )/ &
           ( two*lp1+1 )
      v0 = v1*( r_mesh(1)+r_mesh(2) ) + ( v0r(1)-v0r(2) )/   &
           ( r_mesh(1)-r_mesh(2) )
   endif
!
   eoc2p1=CONE+energy*c2inv
   e2oc2=energy*energy*c2inv
   if (Relativity == NonRelativistic) then
      betap1=lp1
   else
      betap1=sqrt( lp1*lp1-l-fz2oc2 )
   endif
   beta=betap1-ONE
!
!  ===================================================================
!  start off with small r expansion of wave-functions................. 
!  ===================================================================
   if (Relativity == ScalarRelativistic) then
!     ================================================================
!     scalar-relativistic case........................................ 
!     need to be checked out for AtomicNumber = 0 case................
!     ================================================================
      em0= CONE + (energy-v0)*c2inv
      a0 = CONE             ! In full-potential case, this will give wrong
                            ! sine and cosine matrices. We will worry
                            ! about this later
      norm = CONE
      if (atomic_number == 0) then
         a1=ZERO
         a2=-(energy-v0)*em0/(FOUR*l+SIX)
         a3=ZERO
      else
         a1p= ( fz2oc2 + em0*(beta-TWO*fz2oc2) )/(TWO*betap1+ONE)
         a1 = a0*a1p/tzoc2
         a2p= ( a1p*(fz2oc2 + em0*(THREE*betap1+ONE - TWO*fz2oc2)) +  &
                em0*em0*(fz2oc2 - TWO*beta) )/(FOUR*betap1+FOUR)
         a2 = a0*a2p/(tzoc2*tzoc2)
         a3p= ( a2p*(fz2oc2 + em0*(FIVE*betap1+FIVE - TWO*fz2oc2)) -  &
                a1p*(FOUR*betap1+ONE - fz2oc2)*em0*em0 +              &
                (THREE*beta-fz2oc2)*em0*em0*em0 )/(SIX*betap1+NINE)
         a3 = a0*a3p/(tzoc2*tzoc2*tzoc2)
      endif
      npower=1
   else if (Relativity == NonRelativistic) then
!     ================================================================
!     non-relativistic case........................................... 
!     a0 is arbitrary, but this make Rl a spherical Bessel at origin.. 
!     a0 = kappa**l/(2*l+1)!! 
!
!     Note: as of 09/29/03, gms code is written as:
!     I commented out the following line and introduced norm for
!     normalizing the sine and cosine matrices, and the regular solution later
!     ================================================================
!gms  a0 = kappa*lp1
!
      a0=CONE
      norm = CONE
      do j=1,l
         norm=norm*kappa/real(2*j+1,RealKind)
      enddo
!     norm = bjl(1,l)/Grid%r_mesh(1)**(l+1)/kappa
!
      if (atomic_number >= 1) then
         a1p= -Zi/real(lp1,RealKind)
         a1 = a0*a1p
         a2p= (v0-energy - a1p*TWO*Zi)/real(4*l+6,RealKind)
         a2 = a0*a2p
         a3p= (a1p*(v0-energy) - a2p*TWO*Zi)/real(6*l+12,RealKind)
!     ================================================================
!     Note: as of 09/29/03, gms code has the following extra lines for
!           a3p expression
!     ================================================================
!gms     a3p= ( a1p*( (two*AtomicNumber(LocalIndex))**2 +             &
!gms            dble(3*l+4)*(v0-kappa) )/(two*l+3) + v1 )/dble(6*l+12)
         a3 = a0*a3p
         npower=1
      else if ( atomic_number==0 ) then
         a1 = ZERO
         a2 = (v0-energy)/real(4*l+6,RealKind)
         a3 = ZERO
         npower=1
      else
         a1=a0*(v0-energy)/dble(4*l+6)
         a2=a1*(v0-energy)/dble(8*l+20)
         a3=a2*(v0-energy)/dble(12*l+42)
         npower=2
      endif
!     ================================================================
   else
      call ErrorHandler(sname,'Not setup for Dirac equation yet!')
   endif
!
!  ===================================================================
!  get first 4 points from power series expansion..................... 
!  ===================================================================
   do j=1,4
!     ================================================================
!     cm0,pl_reg and ql_reg correspond to M(r), Pl(r) and Ql(r) in 
!     the notes respectively.......................................... 
!     ================================================================
      rp=Grid%r_mesh(j)**npower
      pl_reg(j)=Grid%r_mesh(j)**betap1*(a0+(a1+a2*rp+a3*rp*rp)*rp)
      ql_reg(j)=Grid%r_mesh(j)**beta*                             &
                ( a0*beta+( a1*(beta+npower)+a2*(beta+2*npower)*rp+   &
                            a3*(beta+3*npower)*rp*rp )*rp )/cm0(j)
!     if ( atomic_nnumber==0 ) then
!        pl_reg(j)=bjl(j,l); ql_reg(j)=dbjl(j,l)/kappa
!     endif
!     ================================================================
!     get cl's and sl's
!     ================================================================
!ywg  cl(j)=-bnl(j,1)*pl_reg(j) - bnl(j,0)*ql_reg(j)/kappa
!ywg  sl(j)=-bjl(j,1)*pl_reg(j) - bjl(j,0)*ql_reg(j)/kappa
      cl(j) = ql_reg(j)*bnl(j,l)-pl_reg(j)*dbnl(j,l)*kappa
      sl(j) = ql_reg(j)*bjl(j,l)-pl_reg(j)*dbjl(j,l)*kappa
!     ================================================================
!     get derivatives of cl and sl :: for predictor-corrector use 
!     Note: dcl and dsl correspond to d{Cl(r)}/dx = r*d{Cl(r)}/dr
!           and d{Sl(r)}/dx = r*d{Sl(r)}/dr
!     ================================================================
!ywg  dcl(j)=-bnl(j,1)*(cm0(j)-CONE)*ql_reg(j)*Grid%r_mesh(j)-        &
!ywg         ( llp1/(cm0(j)*Grid%r_mesh(j))+v0r(j)+                   &
!ywg           e2oc2*Grid%r_mesh(j) )*pl_reg(j)*bnl(j,0)/kappa
!ywg  dsl(j)=-bjl(j,1)*(cm0(j)-CONE)*ql_reg(j)*Grid%r_mesh(j)-        &
!ywg         ( llp1/(cm0(j)*Grid%r_mesh(j))+v0r(j)+                   &
!ywg           e2oc2*Grid%r_mesh(j) )*pl_reg(j)*bjl(j,0)/kappa
      dcl(j)=bnl(j,l)*( llp1*(CONE-cm0(j))/Grid%r_mesh(j)+v0r(j)+     &
                        e2oc2*Grid%r_mesh(j) )*pl_reg(j)
      dsl(j)=bjl(j,l)*( llp1*(CONE-cm0(j))/Grid%r_mesh(j)+v0r(j)+     &
                        e2oc2*Grid%r_mesh(j) )*pl_reg(j)
!     ================================================================
   enddo
!
!  ===================================================================
!  get regular solution of scalar relativistic eqs up to the muffin
!  tin sphere......................................................... 
!  ===================================================================
   SolutionType=RegularSolution
!  if ( SphericalSolver .or. isSphPotZeroOutsideRmt ) then
   if (isSphericalSolverOn) then
!     ----------------------------------------------------------------
      call precor4(l,5,Grid%jmt,Grid%hin,SolutionType,                &
                   e2oc2,Grid%r_mesh,cl,dcl,sl,dsl)
!     ----------------------------------------------------------------
   else if (abs(Grid%hin-Grid%hout) < TEN2m8) then ! single logrithmic grid
!     ----------------------------------------------------------------
      call precor4(l,5,iend,Grid%hin,SolutionType,                    &
                   e2oc2,Grid%r_mesh,cl,dcl,sl,dsl)
!     ----------------------------------------------------------------
   else                                            ! double logrithmic grid
!     ----------------------------------------------------------------
      call precor4(l,5,Grid%jmt,Grid%hin,SolutionType,                &
                   e2oc2,Grid%r_mesh,cl,dcl,sl,dsl)
!     ----------------------------------------------------------------
!
!     =================================================================
!     generate pseudo-grid points inside mt sphere and pass equation
!     solution values..................................................
!     =================================================================
      j_inter=Grid%jmt
      clp_reg(0)=cl(j_inter)
      slp_reg(0)=sl(j_inter)
      dclp_reg(0)=dcl(j_inter)
      dslp_reg(0)=dsl(j_inter)
      do j=1,3
         rp=exp(xmt-Grid%hout*j)
!        -------------------------------------------------------------
         call hunt(iend,Grid%r_mesh,rp,j_inter)
!        -------------------------------------------------------------
         j_inter=min(j_inter-(n_inter-1)/2,Grid%jmt-n_inter+1)
         xp=>Grid%r_mesh(j_inter:j_inter+n_inter-1)
         yp=>cl(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,rp,clp_reg(j),z_err_est)
!        -------------------------------------------------------------
         yp=>sl(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,rp,slp_reg(j),z_err_est)
!        -------------------------------------------------------------
         yp=>dcl(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,rp,dclp_reg(j),z_err_est)
!        -------------------------------------------------------------
         yp=>dsl(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,rp,dslp_reg(j),z_err_est)
!        -------------------------------------------------------------
      enddo
!
!     ================================================================
!     get regular solution between jmt+1 and iend regions.............
!     ================================================================
      do j=0,3
         cl(Grid%jmt-j)=clp_reg(j)
         sl(Grid%jmt-j)=slp_reg(j)
         dcl(Grid%jmt-j)=dclp_reg(j)
         dsl(Grid%jmt-j)=dslp_reg(j)
      enddo
!     ----------------------------------------------------------------
      call precor4(l,Grid%jmt+1,iend,Grid%hout,SolutionType,          &
                   e2oc2,Grid%r_mesh,cl,dcl,sl,dsl)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  get normalization etc. at sphere boundary.......................
!  At R, need spherical Bessel functions for all l to get normalization and 
!  physical phase-shifts to determine cl and sl from 
!  Pl=r*g=r*( cl*jl - sl*nl ) and Ql=r*dg/dr/M for all l. 
!  Modify expression to account for difference between spherical and 
!  Ricatti-Bessel fcts. 
!
!  calculate sine and cosine matrices, s_matrix and c_matrix, at the 
!  sphere radius (MT, ASA, or Circumscribing) = r_mesh(jend).......
!  ===================================================================
!
!  ===================================================================
!  get cosine and sine matrices for r = Rmt and store them in
!  c_sphere_mt and s_sphere_mt, respectively.......................... 
!  ===================================================================
   c_sphere_mt = -( l*bnl(Grid%jmt,l)/Grid%r_mesh(Grid%jmt)-          &
                 kappa*bnl(Grid%jmt,lp1) )*pl_reg(Grid%jmt)+          &
                 bnl(Grid%jmt,l)*ql_reg(Grid%jmt)*cm0(Grid%jmt)
   s_sphere_mt = -( l*bjl(Grid%jmt,l)/Grid%r_mesh(Grid%jmt)-          &
                 kappa*bjl(Grid%jmt,lp1) )*pl_reg(Grid%jmt)+          &
                 bjl(Grid%jmt,l)*ql_reg(Grid%jmt)*cm0(Grid%jmt)
!  print *,'delta = ',(s_sphere_mt*bnl(Grid%jmt,l)-                   &
!                      c_sphere_mt*bjl(Grid%jmt,l))/kappa-pl_reg(Grid%jmt)
!
!  ===================================================================
!  In MuffinTin, ASA, and MuffinTin-ASA cases, calculates the regular
!  solution from rmt upto the bounding sphere.
!  ===================================================================
!  if ( SphericalSolver .or. isPotZeroOutsideRmt ) then
!  if (isSphericalSolverOn) then
   if ( isSphericalSolverOn .or. isPotZeroOutsideRmt ) then
      do j=Grid%jmt+1,iend
!        pl_reg(j) = ( s_sphere_mt*bnl(j,l) -c_sphere_mt*bjl(j,l)  )/kappa
         pl_reg(j) = ( sl(Grid%jmt)*bnl(j,l) -cl(Grid%jmt)*bjl(j,l)  )/kappa
      enddo
      do j=Grid%jmt+1,iend
!        ql_reg(j) = ( s_sphere_mt*dbnl(j,l)-c_sphere_mt*dbjl(j,l) )/cm0(j)
         ql_reg(j) = ( sl(Grid%jmt)*dbnl(j,l)-cl(Grid%jmt)*dbjl(j,l) )/cm0(j)
      enddo
   endif
!
!  ===================================================================
!  get cosine and sine matrices for r = bounding sphere and store them
!  in c_matrix and s_matrix, respectively.......................... 
!  ===================================================================
   s_matrix=-(l*bjl(Grid%jend,l)/Grid%r_mesh(Grid%jend)-              &
              kappa*bjl(Grid%jend,lp1) )*pl_reg(Grid%jend)+           &
            bjl(Grid%jend,l)*ql_reg(Grid%jend)*cm0(Grid%jend)
   c_matrix=-(l*bnl(Grid%jend,l)/Grid%r_mesh(Grid%jend)-              &
              kappa*bnl(Grid%jend,lp1) )*pl_reg(Grid%jend)+           &
            bnl(Grid%jend,l)*ql_reg(Grid%jend)*cm0(Grid%jend)
!
!  ===================================================================
!  normalize sine and cosine matrices, and the regular solution
!  ===================================================================
   s_sphere_mt = s_sphere_mt*norm
   c_sphere_mt = c_sphere_mt*norm
   s_matrix = s_matrix*norm
   c_matrix = c_matrix*norm
   s_sph_mt(l) = s_sphere_mt
   c_sph_mt(l) = c_sphere_mt
   s_mtrx(l) = s_matrix
   c_mtrx(l) = c_matrix
!  -------------------------------------------------------------------
   call zscal(iend,norm,pl_reg,1)
   call zscal(iend,norm,ql_reg,1)
!  -------------------------------------------------------------------
!
!  write(6,'(a,i4,4d15.8)')'l,s-mat = ',l,s_sphere_mt,s_matrix
!  write(6,'(a,i4,4d15.8)')'l,c-mat = ',l,c_sphere_mt,c_matrix
!
!  ===================================================================
   nullify( xp, yp, cl, sl, dsl, dcl )
!  ===================================================================
!
   if (stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine solveRadEqn4Reg
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine precor4(l,n1,n2,h0,SolutionType,e2oc2,r_mesh,cl,dcl,sl,dsl)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: l
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=1), intent(in) :: SolutionType
!
   integer (kind=IntKind) :: nstep
   integer (kind=IntKind) :: nstep1
   integer (kind=IntKind) :: nstep2
   integer (kind=IntKind) :: nstep3
   integer (kind=IntKind) :: nstep4
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: icor
   integer (kind=IntKind), parameter :: icmax=5
!
   real (kind=RealKind), intent(in) :: h0
   real (kind=RealKind), intent(in) :: r_mesh(iend)
!
   real (kind=RealKind) :: h
   real (kind=RealKind) :: llp1
!
   real (kind=RealKind), parameter :: nineteen = 19.0d0
   real (kind=RealKind), parameter :: twenty4 = 24.0d0
   real (kind=RealKind), parameter :: thirty7 = 37.0d0
   real (kind=RealKind), parameter :: fifty5 = 55.0d0
   real (kind=RealKind), parameter :: fifty9 = 59.0d0
!
   complex (kind=CmplxKind), intent(in) :: e2oc2
   complex (kind=CmplxKind), intent(inout) :: cl(iend)
   complex (kind=CmplxKind), intent(inout) :: dcl(iend)
   complex (kind=CmplxKind), intent(inout) :: sl(iend)
   complex (kind=CmplxKind), intent(inout) :: dsl(iend)
!
   complex (kind=CmplxKind), pointer :: pl(:)
   complex (kind=CmplxKind), pointer :: ql(:)
   complex (kind=CmplxKind) :: pold
   complex (kind=CmplxKind) :: qold
   complex (kind=CmplxKind) :: clold
   complex (kind=CmplxKind) :: slold
!
   if (SolutionType == RegularSolution) then
      pl=>pl_reg(:)
      ql=>ql_reg(:)
   else if (SolutionType == IrregularSolution) then
      pl=>pl_irr(:)
      ql=>ql_irr(:)
   else
      j=SolutionType
!     ----------------------------------------------------------------
      call ErrorHandler('PRECOR4','Incorrect SolutionType',j)
!     ----------------------------------------------------------------
   endif
!
   if(n2 >= n1) then
      nstep=1
   else
      nstep=-1
   endif
   h=h0*nstep
   nstep1=  nstep
   nstep2=2*nstep
   nstep3=3*nstep
   nstep4=4*nstep
!
   llp1=l*(l+1)
   do j=n1,n2,nstep
!     ================================================================
!     evaluate predictor using Adams-Bashforth formular of order 4
!     ================================================================
      cl(j)=cl(j-nstep1)+h*(fifty5*dcl(j-nstep1)-fifty9*dcl(j-nstep2)+  &
                            thirty7*dcl(j-nstep3)- NINE*dcl(j-nstep4))/twenty4
      sl(j)=sl(j-nstep1)+h*(fifty5*dsl(j-nstep1)-fifty9*dsl(j-nstep2)+  &
                            thirty7*dsl(j-nstep3)- NINE*dsl(j-nstep4))/twenty4
!     ================================================================
!     evaluate corrector using Adams-Moulton formular of order 4...
!     Note: dcl and dsl correspond to d{Cl(r)}/dx = r*d{Cl(r)}/dr
!           and d{Sl(r)}/dx = r*d{Sl(r)}/dr
!     The algebra here needs to be checked for the scalar-relativistic case
!     ================================================================
      do icor=1,icmax
         clold=cl(j)
         slold=sl(j)
!        =============================================================
!ywg     pold=clold*bjl(j,0)-slold*bnl(j,0)
!ywg     qold=-kappa*(clold*bjl(j,1)-slold*bnl(j,1))
!ywg     dcl(j)=-bnl(j,1)*(cm0(j)-CONE)*qold*r_mesh(j)                &
!ywg            -( llp1/(cm0(j)*r_mesh(j))+v0r(j)+e2oc2*r_mesh(j) )   &
!ywg            *pold*bnl(j,0)/kappa
!ywg     dsl(j)=-bjl(j,1)*(cm0(j)-CONE)*qold*r_mesh(j) &
!ywg            -( llp1/(cm0(j)*r_mesh(j))+v0r(j)+e2oc2*r_mesh(j) )   &
!ywg            *pold*bjl(j,0)/kappa
!
         pold = (slold*bnl(j,l)-clold*bjl(j,l))/kappa
         qold = (slold*dbnl(j,l)-clold*dbjl(j,l))/cm0(j)
!        dcl(j)=bnl(j,l)*( llp1/(cm0(j)*r_mesh(j))+v0r(j)+e2oc2*r_mesh(j) )*pold
!        dsl(j)=bjl(j,l)*( llp1/(cm0(j)*r_mesh(j))+v0r(j)+e2oc2*r_mesh(j) )*pold
         dcl(j)=bnl(j,l)*( llp1*(CONE-cm0(j))/Grid%r_mesh(j)+v0r(j)+e2oc2*Grid%r_mesh(j) )*pold
         dsl(j)=bjl(j,l)*( llp1*(CONE-cm0(j))/Grid%r_mesh(j)+v0r(j)+e2oc2*Grid%r_mesh(j) )*pold
!        =============================================================
         cl(j)=cl(j-nstep1)+h*( NINE*dcl(j)+nineteen*dcl(j-nstep1)    &
                               -FIVE*dcl(j-nstep2)+dcl(j-nstep3) )/twenty4
         sl(j)=sl(j-nstep1)+h*( NINE*dsl(j)+nineteen*dsl(j-nstep1)    &
                               -FIVE*dsl(j-nstep2)+dsl(j-nstep3) )/twenty4
      enddo
!     ================================================================
!ywg  pl(j)=cl(j)*bjl(j,0)-sl(j)*bnl(j,0)
!ywg  ql(j)=-kappa*(cl(j)*bjl(j,1)-sl(j)*bnl(j,1))
      pl(j)=(sl(j)*bnl(j,l)-cl(j)*bjl(j,l))/kappa
      ql(j)=(sl(j)*dbnl(j,l)-cl(j)*dbjl(j,l))/cm0(j)
!     ================================================================
   enddo
!
   nullify(pl,ql)
!
   end subroutine precor4
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine precor4p(l,n1,n2,h0,SolutionType,e2oc2,r_mesh,cl,dcl,sl,dsl)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: l
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=1), intent(in) :: SolutionType
!
   integer (kind=IntKind) :: nstep
   integer (kind=IntKind) :: nstep1
   integer (kind=IntKind) :: nstep2
   integer (kind=IntKind) :: nstep3
   integer (kind=IntKind) :: nstep4
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: icor
   integer (kind=IntKind), parameter :: icmax=5
!
   real (kind=RealKind), intent(in) :: h0
   real (kind=RealKind), intent(in) :: r_mesh(iend)
!
   real (kind=RealKind) :: h
   real (kind=RealKind) :: llp1
!
   real (kind=RealKind), parameter :: nineteen = 19.0d0
   real (kind=RealKind), parameter :: twenty4 = 24.0d0
   real (kind=RealKind), parameter :: thirty7 = 37.0d0
   real (kind=RealKind), parameter :: fifty5 = 55.0d0
   real (kind=RealKind), parameter :: fifty9 = 59.0d0
!
   complex (kind=CmplxKind), intent(in) :: e2oc2
   complex (kind=CmplxKind), intent(inout) :: cl(iend)
   complex (kind=CmplxKind), intent(inout) :: dcl(iend)
   complex (kind=CmplxKind), intent(inout) :: sl(iend)
   complex (kind=CmplxKind), intent(inout) :: dsl(iend)
!
   complex (kind=CmplxKind), pointer :: pl(:)
   complex (kind=CmplxKind), pointer :: ql(:)
   complex (kind=CmplxKind) :: pold
   complex (kind=CmplxKind) :: qold
   complex (kind=CmplxKind) :: clold
   complex (kind=CmplxKind) :: slold
!
   if (SolutionType == RegularSolution) then
      pl=>pl_reg(:)
      ql=>ql_reg(:)
   else if (SolutionType == IrregularSolution) then
      pl=>pl_irr(:)
      ql=>ql_irr(:)
   else
      j=SolutionType
!     ----------------------------------------------------------------
      call ErrorHandler('PRECOR4','Incorrect SolutionType',j)
!     ----------------------------------------------------------------
   endif
!
   if(n2 >= n1) then
      nstep=1
   else
      nstep=-1
   endif
   h=h0*nstep
   nstep1=  nstep
   nstep2=2*nstep
   nstep3=3*nstep
   nstep4=4*nstep
!
   llp1=l*(l+1)
   do j=n1,n2,nstep
!     ================================================================
!     evaluate predictor using Adams-Bashforth formular of order 4
!     ================================================================
      cl(j)=cl(j-nstep1)+h*(fifty5*dcl(j-nstep1)-fifty9*dcl(j-nstep2)+  &
                            thirty7*dcl(j-nstep3)- NINE*dcl(j-nstep4))/twenty4
      sl(j)=sl(j-nstep1)+h*(fifty5*dsl(j-nstep1)-fifty9*dsl(j-nstep2)+  &
                            thirty7*dsl(j-nstep3)- NINE*dsl(j-nstep4))/twenty4
!     ================================================================
!     evaluate corrector using Adams-Moulton formular of order 4...
!     Note: dcl and dsl correspond to d{Cl(r)}/dx = r*d{Cl(r)}/dr
!           and d{Sl(r)}/dx = r*d{Sl(r)}/dr
!     ================================================================
      do icor=1,icmax
         clold=cl(j)
         slold=sl(j)
!        =============================================================
         pold=clold*bjl(j,0)-slold*bnl(j,0)
         qold=-kappa*(clold*bjl(j,1)-slold*bnl(j,1))
         dcl(j)=-bnl(j,1)*(cm0(j)-CONE)*qold*r_mesh(j)                &
                -( llp1/(cm0(j)*r_mesh(j))+v0r(j)+e2oc2*r_mesh(j) )   &
                *pold*bnl(j,0)/kappa
         dsl(j)=-bjl(j,1)*(cm0(j)-CONE)*qold*r_mesh(j) &
                -( llp1/(cm0(j)*r_mesh(j))+v0r(j)+e2oc2*r_mesh(j) )   &
                *pold*bjl(j,0)/kappa
!
!new     pold = (slold*bnl(j,l)-clold*bjl(j,l))/kappa
!new     qold = slold*dbnl(j,l)-clold*dbjl(j,l)
!new     dcl(j)=bnl(j,l)*( llp1*(CONE-cm0(j))/r_mesh(j)+v0r(j)+       &
!new                       e2oc2*r_mesh(j) )*pold
!new     dsl(j)=bjl(j,l)*( llp1*(CONE-cm0(j))/r_mesh(j)+v0r(j)+       &
!new                       e2oc2*r_mesh(j) )*pold
!        =============================================================
         cl(j)=cl(j-nstep1)+h*( NINE*dcl(j)+nineteen*dcl(j-nstep1)    &
                               -FIVE*dcl(j-nstep2)+dcl(j-nstep3) )/twenty4
         sl(j)=sl(j-nstep1)+h*( NINE*dsl(j)+nineteen*dsl(j-nstep1)    &
                               -FIVE*dsl(j-nstep2)+dsl(j-nstep3) )/twenty4
      enddo
!     ================================================================
      pl(j)=cl(j)*bjl(j,0)-sl(j)*bnl(j,0)
      ql(j)=-kappa*(cl(j)*bjl(j,1)-sl(j)*bnl(j,1))
!new  pl(j)=(sl(j)*bnl(j,l)-cl(j)*bjl(j,l))/kappa
!new  ql(j)=sl(j)*dbnl(j,l)-cl(j)*dbjl(j,l)
!     ================================================================
   enddo
!
   nullify(pl,ql)
!
   end subroutine precor4p
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calPhiLr(kl,fmem,np,wfr,dwfr,irr_in)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
   use InterpolationModule, only : FitInterp
   use TimerModule, only : getTime
   implicit   none
!
   character (len=8), parameter :: sname='calPhiLr'
!
   integer (kind=IntKind), intent(in) :: kl, np
!
   logical, intent(in), optional :: irr_in
   logical :: irr
!
   integer (kind=IntKind) :: ir, jmt, jend
   integer (kind=IntKind) :: l, klp, lp, jmax_pot_loc, kmax_pot_loc
   integer (kind=IntKind) :: llp1, j, kl1, l1, m1, jl1, klpp, nj3_i3, lpp, i1
   integer (kind=IntKind), pointer :: flags(:), kj3_i3(:)
!
   real (kind=RealKind) :: h, hstep, t0, t1
   real (kind=RealKind), pointer :: r_mesh(:), cgnt_i3(:)
!
   complex (kind=CmplxKind), intent(inout), target :: fmem(iend,2)
   complex (kind=CmplxKind), intent(in) :: wfr(np), dwfr(np)
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
   complex (kind=CmplxKind), pointer :: sx(:,:)
   complex (kind=CmplxKind), pointer :: cx(:,:)
   complex (kind=CmplxKind), pointer :: sxklp(:)
   complex (kind=CmplxKind), pointer :: cxklp(:)
   complex (kind=CmplxKind) :: dcx(kmax_phi,4)
   complex (kind=CmplxKind) :: dsx(kmax_phi,4)
   complex (kind=CmplxKind) :: potmp(jmax_potmp)
   complex (kind=CmplxKind) :: cpl, dcpl, e2oc2
   complex (kind=CmplxKind) :: gaunt_pot(kmax_phi,kmax_phi)
   complex (kind=CmplxKind) :: vpsum(kmax_phi), pwave(kmax_phi)
!
   if (present(irr_in)) then
      irr = irr_in
   else
      irr = .false.
   endif
!
!  ===================================================================
!  calculate: plhat_reg, plhat_irr, qlhat_reg, and qlhat_irr
!
!  use plhat_reg and qlhat_reg as swapping space for sx and cx.
!  ===================================================================
!
   e2oc2=energy*energy*c2inv
   r_mesh=>Grid%r_mesh(1:iend)
   sxklp=>fmem(1:iend,1)
   cxklp=>fmem(1:iend,2)
   l=lofk(kl)
   llp1 = l*(l+1)
!
!  ===================================================================
!  Using Adams-Moulton predictor-corrector method to solve for S(r) and
!  C(r). Here, we are not truncating the scalar relativistic term. This
!  could potentially be a problem.
!  ===================================================================
   if (irr) then  ! The irregular solution needs to be checked...
      sx=>plhat_irr(1:iend,1:kmax_phi)
      cx=>qlhat_irr(1:iend,1:kmax_phi)
      sx = CZERO; dsx = CZERO
      cx = CZERO; dcx = CZERO
!     ================================================================
!     calculate sx, dsx, cx, dcx, at points: 1, 2, 3, 4 outside numrs_cs
!     ================================================================
      do j = 1, np
         ir = numrs_cs+np-j+1
         sx(ir,kl) = bjl(ir,l)*dwfr(j)*cm0(ir)-kappa*wfr(j)*dbjl(ir,l)
         cx(ir,kl) = bnl(ir,l)*dwfr(j)*cm0(ir)-kappa*wfr(j)*dbnl(ir,l)
      enddo
!
      h=Grid%hout
      pot_l => pot_trunc_jl
      flags => flags_trunc_jl
      jmax_pot_loc = jmax_trunc
      kmax_pot_loc = kofj(jmax_pot_loc)
!
!     ================================================================
!     When using predictor-corrector scheme, it is found necessary to
!     calcutate the initial dsx and dcx on a potential extrapolated
!     from the truncated potential beyond the bounding sphere.
!     ================================================================
      kmax_pot_loc = kofj(jmax_pot_loc)
      do j = 1, 4
         ir = numrs_cs+5-j
         potmp = CZERO
         do jl1 = 1, jmax_pot_loc
            if (flags(jl1) /= 0) then
!              -------------------------------------------------------
               call FitInterp(numrs_cs-numrs_mt,Grid%r_mesh(numrs_mt+1:numrs_cs), &
                                 pot_l(:,jl1),r_mesh(ir),cpl,dcpl)
!              -------------------------------------------------------
               if (jl1 == 1) then
                  potmp(jl1) = llp1*(CONE-cm0(ir))/(Y0*r_mesh(ir)**2)+cpl+(e2oc2+PotShift)/Y0
               else
                  potmp(jl1) = cpl
               endif
            endif
         enddo
         do klp=1,kmax_phi
            do klpp=1,kmax_phi
               gaunt_pot(klpp,klp) = CZERO
               nj3_i3 = nj3(klpp,klp)
               kj3_i3 => kj3(1:nj3_i3,klpp,klp)
               cgnt_i3 => cgnt(1:nj3_i3,klpp,klp)
               do i1=1,nj3_i3
                  kl1 = kj3_i3(i1)
                  m1 = mofk(kl1); jl1 = jofk(kl1)
                  if (m1 >= 0) then
                     gaunt_pot(klpp,klp) = gaunt_pot(klpp,klp) + cgnt_i3(i1)*potmp(jl1)
                  else
                     gaunt_pot(klpp,klp) = gaunt_pot(klpp,klp) + cgnt_i3(i1)*m1m(m1)*conjg(potmp(jl1))
                  endif
               enddo
            enddo
         enddo
!
         do klpp=1,kmax_phi
            lpp=lofk(klpp)
            pwave(klpp) = (sx(ir,klpp)*bnl(ir,lpp)-cx(ir,klpp)*bjl(ir,lpp))/kappa
         enddo
         call zgemv('t',kmax_phi,kmax_phi,CONE,gaunt_pot,kmax_phi,pwave,1,CZERO,vpsum,1)
         do klp=1,kmax_phi
            lp=lofk(klp)
!           ==========================================================
!           update dsx and dcx at j, and use them for the corrector
!           nstep here gives correct sign, depending outgoing or incoming solutions
!           ==========================================================
            dsx(klp,j)=r_mesh(ir)*vpsum(klp)*bjl(ir,lp)
            dcx(klp,j)=r_mesh(ir)*vpsum(klp)*bnl(ir,lp)
         enddo
      enddo
!     ----------------------------------------------------------------
      call solveSCr(kl,numrs_cs,numrs_mt+1,h,sx,cx,dsx,dcx,numrs_mt-1,jmax_pot_loc,pot_l,flags)
!     ----------------------------------------------------------------
!
      do j = 1, 4
         ir = numrs_mt+5-j
         potmp = CZERO
         do jl1 = 1, jmax_pot_loc
            if (flags(jl1) /= 0) then
               if (jl1 == 1) then
                  potmp(jl1) = llp1*(CONE-cm0(ir))/(Y0*r_mesh(ir)**2)+pot_l(ir-numrs_mt+1,jl1)+(e2oc2+PotShift)/Y0
               else
                  potmp(jl1) = pot_l(ir-numrs_mt+1,jl1)
               endif
            endif
         enddo
         do klp=1,kmax_phi
            do klpp=1,kmax_phi
               gaunt_pot(klpp,klp) = CZERO
               nj3_i3 = nj3(klpp,klp)
               kj3_i3 => kj3(1:nj3_i3,klpp,klp)
               cgnt_i3 => cgnt(1:nj3_i3,klpp,klp)
               do i1=1,nj3_i3
                  kl1 = kj3_i3(i1)
                  m1 = mofk(kl1); jl1 = jofk(kl1)
                  if (m1 >= 0) then
                     gaunt_pot(klpp,klp) = gaunt_pot(klpp,klp) + cgnt_i3(i1)*potmp(jl1)
                  else
                     gaunt_pot(klpp,klp) = gaunt_pot(klpp,klp) + cgnt_i3(i1)*m1m(m1)*conjg(potmp(jl1))
                  endif
               enddo
            enddo
         enddo
!
         do klpp=1,kmax_phi
            lpp=lofk(klpp)
            pwave(klpp) = (sx(ir,klpp)*bnl(ir,lpp)-cx(ir,klpp)*bjl(ir,lpp))/kappa
         enddo
         call zgemv('t',kmax_phi,kmax_phi,CONE,gaunt_pot,kmax_phi,pwave,1,CZERO,vpsum,1)
         do klp=1,kmax_phi
            lp=lofk(klp)
!           ==========================================================
!           update dsx and dcx at j, and use them for the corrector
!           nstep here gives correct sign, depending outgoing or incoming solutions
!           ==========================================================
            dsx(klp,j)=r_mesh(ir)*vpsum(klp)*bjl(ir,lp)
            dcx(klp,j)=r_mesh(ir)*vpsum(klp)*bnl(ir,lp)
         enddo
      enddo
!
      h=Grid%hin
      pot_l => pot_jl
      flags => flags_jl
      jmax_pot_loc = jmax_pot
!     ----------------------------------------------------------------
      call solveSCr(kl,numrs_mt,1,h,sx,cx,dsx,dcx,0,jmax_pot_loc,pot_l,flags)
!     ----------------------------------------------------------------
!
!     ================================================================
!     Calculate plhat_irr and qlhat_irr, which correspond to the deri-
!     vative of plhat_irr.
!     ================================================================
      do klp=1,kmax_phi
!        ==============================================================
!        Since plhat_reg and sx share the same space, use sxklp to swap
!        Since qlhat_reg and cx share the same space, use cxklp to swap
!        -------------------------------------------------------------
         call zcopy(iend,sx(1,klp),1,sxklp,1)
         call zcopy(iend,cx(1,klp),1,cxklp,1)
!        -------------------------------------------------------------
!
!        ==============================================================
!        calculate plhat_irr and qlhat_irr
!        ==============================================================
         lp=lofk(klp)
         do ir=1,iend
            plhat_irr(ir,klp)=(bnl(ir,lp)*sxklp(ir)-bjl(ir,lp)*cxklp(ir))/kappa
         enddo
         do ir=1,iend
            qlhat_irr(ir,klp)=dbnl(ir,lp)*sxklp(ir)-dbjl(ir,lp)*cxklp(ir)
         enddo
      enddo
   else
      h=Grid%hin
      pot_l => pot_jl
      flags => flags_jl
      jmax_pot_loc = jmax_pot
!
      sx=>plhat_reg(1:iend,1:kmax_phi)
      cx=>qlhat_reg(1:iend,1:kmax_phi)
      dsx = CZERO; dcx = CZERO
!     sx = CZERO; dsx = CZERO
!     cx = CZERO; dcx = CZERO
!     ================================================================
!     calculate sx, dsx, cx, dcx, at points: 1, 2, 3, 4 starting from the origin
!     ================================================================
      do j = 1, 4
         sx(j,kl) = bjl(j,l)*dwfr(j)*cm0(j)-kappa*wfr(j)*dbjl(j,l)
         cx(j,kl) = bnl(j,l)*dwfr(j)*cm0(j)-kappa*wfr(j)*dbnl(j,l)
         dsx(kl,j)=r_mesh(j)*bjl(j,l)*(llp1*(CONE-cm0(j))/(Y0*r_mesh(j)**2)+pot_l(j,1)+(e2oc2+PotShift)/Y0)*wfr(j)
         dcx(kl,j)=r_mesh(j)*bnl(j,l)*(llp1*(CONE-cm0(j))/(Y0*r_mesh(j)**2)+pot_l(j,1)+(e2oc2+PotShift)/Y0)*wfr(j)
!        if (abs(sx(j,kl)) > TEN2m8 .or. abs(cx(j,kl)) > TEN2m8) then
!           write(6,'(a,i4,8d15.7)')'sx,cx = ',kl,sx(j,kl),cx(j,kl),dsx(kl,j),dcx(kl,j)
!        endif
      enddo
      t0 = getTime()
!     ----------------------------------------------------------------
      call solveSCr(kl,5,numrs_mt,h,sx,cx,dsx,dcx,0,jmax_pot_loc,pot_l,flags)
!     ----------------------------------------------------------------
      t1 = getTime() - t0
!     do klp = 1, kmax_phi
!        if (abs(sx(numrs_mt,klp)) > TEN2m8 .or. abs(cx(numrs_mt,klp)) > TEN2m8) then
!           write(6,'(a,2i4,4d15.7)')'sx,cx at numrs_mt = ',kl,klp,sx(numrs_mt,klp),cx(numrs_mt,klp)
!        endif
!     enddo
!
      kmax_pot_loc = kofj(jmax_pot_loc)
      do j = 1, 4
         ir = numrs_mt-4+j
         potmp = CZERO
         do jl1 = 1, jmax_pot_loc
            if (flags(jl1) /= 0) then
               if (jl1 == 1) then
                  potmp(jl1) = llp1*(CONE-cm0(ir))/(Y0*r_mesh(ir)**2)+pot_l(ir,jl1)+(e2oc2+PotShift)/Y0
               else
                  potmp(jl1) = pot_l(ir,jl1)
               endif
            endif
         enddo
         do klp=1,kmax_phi
            do klpp=1,kmax_phi
               gaunt_pot(klpp,klp) = CZERO
               nj3_i3 = nj3(klpp,klp)
               kj3_i3 => kj3(1:nj3_i3,klpp,klp)
               cgnt_i3 => cgnt(1:nj3_i3,klpp,klp)
               do i1=1,nj3_i3
                  kl1 = kj3_i3(i1)
                  m1 = mofk(kl1); jl1 = jofk(kl1)
                  if (m1 >= 0) then
                     gaunt_pot(klpp,klp) = gaunt_pot(klpp,klp) + cgnt_i3(i1)*potmp(jl1)
                  else
                     gaunt_pot(klpp,klp) = gaunt_pot(klpp,klp) + cgnt_i3(i1)*m1m(m1)*conjg(potmp(jl1))
                  endif
               enddo
            enddo
         enddo
!
         do klpp=1,kmax_phi
            lpp=lofk(klpp)
            pwave(klpp) = (sx(ir,klpp)*bnl(ir,lpp)-cx(ir,klpp)*bjl(ir,lpp))/kappa
         enddo
         call zgemv('t',kmax_phi,kmax_phi,CONE,gaunt_pot,kmax_phi,pwave,1,CZERO,vpsum,1)
         do klp=1,kmax_phi
            lp=lofk(klp)
!           ==========================================================
!           update dsx and dcx at j, and use them for the corrector
!           nstep here gives correct sign, depending outgoing or incoming solutions
!           ==========================================================
            dsx(klp,j)=r_mesh(ir)*vpsum(klp)*bjl(ir,lp)
            dcx(klp,j)=r_mesh(ir)*vpsum(klp)*bnl(ir,lp)
!           if (abs(dsx(klp,j)) > TEN2m8 .or. abs(dcx(klp,j)) > TEN2m8) then
!              write(6,'(a,2i4,4d15.7)')'dsx,dcx at j,klp = ',j,klp,dsx(klp,j),dcx(klp,j)
!           endif
         enddo
      enddo
!
      h=Grid%hout
!     write(6,'(a,d15.7)')'hout = ',h
      t0 = getTime()
      if ( SSSMethod==2 ) then
         pot_l => pot_trunc_jl
         flags => flags_trunc_jl
         jmax_pot_loc = jmax_trunc
!        =============================================================
!        The reference point was set to be N0=numrs_mt-1. This should be
!        a bug that causes the FCC structure having perculiar E vs a0 curve.
!        It is now N0 = numrs_mt
!        -- Y.W., Aug. 10, 2015
!        As a matter of fact, pot_trunc_jl(1) is the truncated potential  
!        at  r = r_mesh(numrs_mt). The original code is right, and is now
!        put back.
!        -- Y.W., Aug. 16, 2015
!        -------------------------------------------------------------
         call solveSCr(kl,numrs_mt+1,numrs_cs,h,sx,cx,dsx,dcx,numrs_mt-1,jmax_pot_loc,pot_l,flags)
!        -------------------------------------------------------------
!        do jl1 = 1, jmax_pot_loc
!           if (abs(pot_l(1,jl1)) > TEN2m8) then
!              write(6,'(a,3i4,4d15.7)')'pot_l at rcs,jl1 = ',lofj(jl1),mofj(jl1),flags(jl1), &
!                                       pot_l(numrs_cs-numrs_mt+1,jl1), &
!                                       Scatter(LocalIndex)%step_rad(numrs_cs-numrs_mt+1,jl1)
!           endif
!        enddo
      else
!        -------------------------------------------------------------
         call solveSCr(kl,numrs_mt+1,numrs_cs,h,sx,cx,dsx,dcx,0,jmax_pot_loc,pot_l,flags)
!        -------------------------------------------------------------
      endif
#ifdef TIMING
      t1 = t1 + (getTime() - t0)
      write(6,'(a,f10.5)')'In calPhiLr, total timing for solveSCr = ',t1
#endif
      do klp = 1, kmax_phi
!        if (abs(sx(numrs_cs,klp)) > TEN2m8 .or. abs(cx(numrs_cs,klp)) > TEN2m8) then
!           write(6,'(a,2i4,4d15.7)')'sx,cx at kl,klp = ',kl,klp,sx(numrs_cs,klp),cx(numrs_cs,klp)
!        endif
         do ir = numrs_cs+1, iend
            sx(ir,klp) = sx(numrs_cs,klp)
            cx(ir,klp) = cx(numrs_cs,klp)
         enddo
         if (abs(sx(numrs_cs,klp)) > TEN2m8 .or. abs(cx(numrs_cs,klp)) > TEN2m8) then
            Scatter(LocalIndex)%sol_flags(klp,kl) = 1
         endif
      enddo
!
!     ================================================================
!     Calculate plhat_reg and qlhat_reg, which correspond to the deri-
!     vative of plhat_reg.
!     ================================================================
      do klp=1,kmax_phi
!        ==============================================================
!        Since plhat_reg and sx share the same space, use sxklp to swap
!        Since qlhat_reg and cx share the same space, use cxklp to swap
!        -------------------------------------------------------------
         call zcopy(iend,sx(1,klp),1,sxklp,1)
         call zcopy(iend,cx(1,klp),1,cxklp,1)
!        -------------------------------------------------------------
!
!        ==============================================================
!        calculate plhat_reg and qlhat_reg
!        ==============================================================
         lp=lofk(klp)
         do ir=1,iend
            plhat_reg(ir,klp)=(bnl(ir,lp)*sxklp(ir)-bjl(ir,lp)*cxklp(ir))/kappa
         enddo
         do ir=1,iend
            qlhat_reg(ir,klp)=dbnl(ir,lp)*sxklp(ir)-dbjl(ir,lp)*cxklp(ir)
         enddo
      enddo
   endif
!write(6,'(a,3i5)')'In the end calPhiLr ...',MyPE,MyPEInGroup,kl
!  -------------------------------------------------------------------
   nullify(sx, cx, sxklp, cxklp)
!  -------------------------------------------------------------------
!
   if(stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine calPhiLr
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine solveIntEqn(atom,kl,fmem,reg,irr)
!  ===================================================================
   implicit   none
!
   character (len=11), parameter :: sname='solveIntEqn'
!
   logical, intent(in) :: reg, irr
!
   integer (kind=IntKind), intent(in) :: atom, kl
!
   integer (kind=1) :: SolutionType
!
   integer (kind=IntKind) :: klp,lp
!
   complex (kind=CmplxKind), intent(inout), target :: fmem(iend,5)
!
   complex (kind=CmplxKind) :: plrs, qlrs
!#ifdef CHECK_WRONSKI
!   integer (kind=IntKind)   :: i
!   complex (kind=CmplxKind) :: norm
!   complex (kind=CmplxKind), pointer :: fint(:,:)
!#endif
!
!  ===================================================================
!  calculate: plhat_reg, plhat_irr, qlhat_reg, and qlhat_irr
!  ===================================================================
!
   if (reg) then
      qlhat_reg(1:iend,1:kmax_phi)=CZERO
   endif
!
   if ( isIrrSolOn .and. irr ) then
      qlhat_irr(1:iend,1:kmax_phi)=CZERO
   endif
!
!  if ( FullSolver .and. .not.isPotZeroOutsideRmt ) then
   if (.not.isSphericalSolverOn) then
      if ( reg ) then
         isDiagSolv = .true.
         SolutionType=RegularSolution
!        -------------------------------------------------------------
         call intkvpr(atom,SolutionType,kl,fmem)
!        -------------------------------------------------------------
         isDiagSolv = .false.
!        -------------------------------------------------------------
         call intkvpr(atom,SolutionType,kl,fmem)
!        -------------------------------------------------------------
      endif
      if ( isIrrSolOn .and. irr ) then
         isDiagSolv = .true.
         SolutionType = IrregularSolution
!        -------------------------------------------------------------
         call intkvpr(atom,SolutionType,kl,fmem)
!        -------------------------------------------------------------
         isDiagSolv = .false.
!        -------------------------------------------------------------
         call intkvpr(atom,SolutionType,kl,fmem)
!        -------------------------------------------------------------
      endif
   endif
!
!  ===================================================================
!  add pl_reg into plhat_reg, and ql_reg into qlhat_red
!  add pl_irr into plhat_irr, and ql_irr into qlhat_irr
!  ===================================================================
   if ( reg ) then
!     plhat_reg(1:iend,kl) = plhat_reg(1:iend,kl)+pl_reg(1:iend)
!     qlhat_reg(1:iend,kl) = qlhat_reg(1:iend,kl)+ql_reg(1:iend)
!     ----------------------------------------------------------------
      call zaxpy(iend,CONE,pl_reg,1,plhat_reg(1,kl),1)
      call zaxpy(iend,CONE,ql_reg,1,qlhat_reg(1,kl),1)
!     ----------------------------------------------------------------
   endif
   if ( isIrrSolOn .and. irr ) then
!     plhat_irr(1:iend,kl) = plhat_irr(1:iend,kl)+pl_irr(1:iend)
!     qlhat_irr(1:iend,kl) = qlhat_irr(1:iend,kl)+ql_irr(1:iend)
!     ----------------------------------------------------------------
      call zaxpy(iend,CONE,pl_irr,1,plhat_irr(1,kl),1)
      call zaxpy(iend,CONE,ql_irr,1,qlhat_irr(1,kl),1)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  calculate : smpicm = smat + i*cmat
!  Note: Here both smat and cmat are sine and cosine matrices for the
!        smooth potential
!  ===================================================================
!  if ( FullSolver .and. SSSMethod==-1 .and. reg ) then
   if (SSSMethod==-1 .and. reg .and. .not.isSphericalSolverOn) then
      do klp = 1,kmax_kkr
         plrs = plhat_reg(numrs_cs,klp)
         qlrs = qlhat_reg(numrs_cs,klp)
         lp   = lofk(klp)
         smpicm(klp,kl) = bhl(numrs_cs,lp)*qlrs -                     &
                          kappa*dbhl(numrs_cs,lp)*plrs
      enddo
   endif
!
!!#ifdef CHECK_WRONSKI
!  ===================================================================
!  check Wronskian relations.
!  ===================================================================
!   write(6,'(/,'' kl = '',1i5)')kl
!! do klp=1,kmax_phi
!!    write(6,'(1i5,5x,2d20.13)')klp,plhat_reg(Grid%jmt+10,klp)
!! enddo
!  -------------------------------------------------------------------
!   allocate(fint(1:iend,1))
!  -------------------------------------------------------------------
!   fint(1:iend,1)=czero
!   do klp=1,kmax_phi
!      fmem(1:iend,1)=plhat_reg(1:iend,klp)
!      fmem(1:iend,2)=plhat_irr(1:iend,klp)
!      fmem(1:iend,3)=qlhat_reg(1:iend,klp)
!      fmem(1:iend,4)=qlhat_irr(1:iend,klp)
!      fint(1:iend,1)=fint(1:iend,1)+fmem(1:iend,2)*fmem(1:iend,3)-    &
!                                    fmem(1:iend,1)*fmem(1:iend,4)
!   enddo
!   norm=fint(Grid%jend,1)
!   if(abs(norm) > TEN2m8) then
!      do i=1,iend,20
!         write(6,'('' i, Wronskian ='',1i5,2d15.8)')i,fint(i,1)/norm
!      enddo
!      write(6,'('' i, Wronskian ='',1i5,2d15.8)')iend,fint(iend,1)/norm
!      do i=Grid%jmt-5,Grid%jmt+5
!          write(6,'('' i, Wronskian ='',1i5,2d15.8)')i,fint(i,1)/norm
!      enddo
!      do i=Grid%jend-10,Grid%jend
!         write(6,'('' i, Wronskian ='',1i5,2d15.8)')i,fint(i,1)/norm
!      enddo
!   endif
!  -------------------------------------------------------------------
!   deallocate(fint)
!  -------------------------------------------------------------------
!#endif
!
   if(stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine solveIntEqn
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine solveSCr(kl,N1,N2,h,sx,cx,dsx,dcx,N0,jmax_pot_loc,pot_l,flags)
!  ===================================================================
!
! *===================================================================
! *                                                                  *
! * Using Runge-Kutta fourth-order method to solve for the ordinary  *
! * differential equations:                                          *
! *                                                                  *
! * ^p                                                               *
! * S  (r) = f  (r) +                                                *
! *  LpL      LpL                                                    *
! *                L1   ^             ^      ^      ^      ^         *
! *          sum  C   * j (kr)*v (r)*[n (kr)*S  (r)-j (kr)*C  (r)]/k *
! *         L'',L1 L''Lp  lp    L1     l''    L''L   l''    L''L     *
! *                                                                  *
! *                                                                  *
! * ^p                                                               *
! * C  (r) = g  (r) +                                                *
! *  LpL      LpL                                                    *
! *                L1   ^             ^      ^      ^      ^         *
! *          sum  C   * n (kr)*v (r)*[n (kr)*S  (r)-j (kr)*C  (r)]/k *
! *         L'',L1 L''Lp lp     L1     l''    L''L   l''    L''L     *
! *                                                                  *
! * where,                                                           *
! *                                                                  *
! * k = sqrt(e)                                                      *
! *                                                                  *
! *  L1                                                              *
! * C    : Gaunt factors                                             *
! *  L''Lp                                                           *
! *                                                                  *
! * ^                                                                *
! * j (kr) = k*r*j (kr), j (kr) is the spherical Bessel function     *
! *  l            l       l                                          *
! *                                                                  *
! * ^                                                                *
! * n (kr) = k*r*n (kr), n (kr) is the spherical Neumann function    *
! *  l            l       l                                          *
! *                                                                  *
! *                L1   ^                                            *
! * f  (r) = sum  C   * j (kr)*v (r)*p (r)                           *
! *  LpL     L1>0  LLp   lp     L1    L                              *
! *                                                                  *
! *                L1   ^                                            *
! * g  (r) = sum  C   * n (kr)*v (r)*p (r)                           *
! *  LpL     L1>0  LLp   lp     L1    L                              *
! *                                                                  *
! *                                                                  *
! * p (r) = r*phi (r)                                                *
! *  L           l                                                   *
! *                                                                  *
! * phi (r): the solution of Schrodinger equation corresponding to   *
! *    l     the spherical part of the potential.                    *
! *                                                                  *
! *===================================================================
!
   use LdaCorrectionModule, only : checkLdaCorrection, getPotentialCorrection
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: kl
   integer (kind=IntKind), intent(in) :: N1, N2, N0, jmax_pot_loc
   integer (kind=IntKind), intent(in) :: flags(jmax_pot_loc)
!
   integer (kind=IntKind) :: nstep, icor, kmax_pot_loc
   integer (kind=IntKind) :: ir, irm1, irmn0, j, jm1, jm2, jm3
   integer (kind=IntKind) :: i1, klp, klpp, kl1, jl1, m1, l1, lp, lpp
   integer (kind=IntKind) :: nj3_i3, llp1
!
   integer (kind=IntKind), parameter :: icmax = 5
!
   integer (kind=IntKind), pointer :: kj3_i3(:)
!
   real (kind=RealKind), intent(in) :: h
!
   real (kind=RealKind) :: hstep, hfac
   real (kind=RealKind), pointer :: cgnt_i3(:)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), parameter :: Nine    =  9.0d0
   real (kind=RealKind), parameter :: Twenty4 = 24.0d0
   real (kind=RealKind), parameter :: Thirty7 = 37.0d0
   real (kind=RealKind), parameter :: Fifty5  = 55.0d0
   real (kind=RealKind), parameter :: Fifty9  = 59.0d0
   real (kind=RealKind), parameter :: Nineteen = 19.0d0
!
   complex (kind=CmplxKind), intent(inout) :: cx(iend,kmax_phi)
   complex (kind=CmplxKind), intent(inout) :: sx(iend,kmax_phi)
   complex (kind=CmplxKind), intent(inout) :: dcx(kmax_phi,4)
   complex (kind=CmplxKind), intent(inout) :: dsx(kmax_phi,4)
!
   complex (kind=CmplxKind), intent(in) :: pot_l(:,:)
!
   complex (kind=CmplxKind) :: cfac, e2oc2
   complex (kind=CmplxKind) :: gaunt_pot(kmax_phi,kmax_phi)
   complex (kind=CmplxKind) :: vpsum(kmax_phi), pwave(kmax_phi)
   complex (kind=CmplxKind) :: sxir(kmax_phi), cxir(kmax_phi)
   complex (kind=CmplxKind) :: sxirm1(kmax_phi), cxirm1(kmax_phi)
   complex (kind=CmplxKind) :: potmp(jmax_potmp)
!
   e2oc2=energy*energy*c2inv
!
   if (N2 > N1) then
      nstep=1
   else if (N2 < N1) then
      nstep=-1
   else
      call ErrorHandler('solveSCr','improper N1 and N2 value',N1,N2)
   endif
!
   r_mesh=>Grid%r_mesh(1:iend)
!
   llp1 = lofk(kl)*(lofk(kl)+1)
   hstep=abs(h)*nstep
   hfac=hstep/Twenty4
   kmax_pot_loc = kofj(jmax_pot_loc)
!
!  ===================================================================
!  Solve the rest of the points using Adams-Bashforth 4-step method
!  ===================================================================
   j=4; jm1 = 3; jm2 = 2; jm3 = 1
   do klp=1,kmax_phi
      sxirm1(klp) = sx(N1-nstep,klp)
      cxirm1(klp) = cx(N1-nstep,klp)
   enddo
   do ir=N1, N2, nstep
      irm1=ir-nstep
      irmn0=ir-N0
!     ================================================================
!     calculate sx and cx using dsx(j,klp) and dcx(j,klp), j=1,2,3,4
!     ================================================================
      do klp=1,kmax_phi
         sxir(klp)=sxirm1(klp)+hfac*( Fifty5 *dsx(klp,j)              &
                                     -Fifty9 *dsx(klp,jm1)            &
                                     +Thirty7*dsx(klp,jm2)            &
                                     -Nine   *dsx(klp,jm3) )
         cxir(klp)=cxirm1(klp)+hfac*( Fifty5 *dcx(klp,j)              &
                                     -Fifty9 *dcx(klp,jm1)            &
                                     +Thirty7*dcx(klp,jm2)            &
                                     -Nine   *dcx(klp,jm3) )
      enddo
!     ================================================================
      potmp = CZERO
      do jl1 = 1, jmax_pot_loc
         if (flags(jl1) /= 0) then
            if (jl1 == 1) then
               potmp(jl1) = llp1*(CONE-cm0(ir))/(Y0*r_mesh(ir)**2)+pot_l(irmn0,jl1)+(e2oc2+PotShift)/Y0
            else
               potmp(jl1) = pot_l(irmn0,jl1)
            endif
         endif
      enddo
!
      do klp=1,kmax_phi
         do klpp=1,kmax_phi
            gaunt_pot(klpp,klp) = CZERO
            nj3_i3 = nj3(klpp,klp)
            kj3_i3 => kj3(1:nj3_i3,klpp,klp)
            cgnt_i3 => cgnt(1:nj3_i3,klpp,klp)
            do i1=1,nj3_i3
               kl1 = kj3_i3(i1)
               m1 = mofk(kl1); jl1 = jofk(kl1)
               if (m1 >= 0) then
                  gaunt_pot(klpp,klp) = gaunt_pot(klpp,klp) + cgnt_i3(i1)*potmp(jl1)
               else
                  gaunt_pot(klpp,klp) = gaunt_pot(klpp,klp) + cgnt_i3(i1)*m1m(m1)*conjg(potmp(jl1))
               endif
            enddo
         enddo
      enddo
!
!     ================================================================
!     evaluate corrector using Adams-Moulton formular of order 4...
!     Note: dcl and dsl correspond to d{Cl(r)}/dx = r*d{Cl(r)}/dr
!           and d{Sl(r)}/dx = r*d{Sl(r)}/dr
!     ================================================================
      j  =jm3
      jm1=mod(j+2,4)+1
      jm2=mod(j+1,4)+1
      jm3=mod(j,  4)+1
      do icor=1,icmax
         do klpp=1,kmax_phi
            lpp=lofk(klpp)
            pwave(klpp) = (sxir(klpp)*bnl(ir,lpp)-cxir(klpp)*bjl(ir,lpp))/kappa
         enddo
         call zgemv('t',kmax_phi,kmax_phi,CONE,gaunt_pot,kmax_phi,pwave,1,CZERO,vpsum,1)
         do klp=1,kmax_phi
            lp=lofk(klp)
!           ==========================================================
!           update dsx and dcx at j, and use them for the corrector
!           nstep here gives correct sign, depending outgoing or incoming solutions
!           ==========================================================
            dsx(klp,j)=nstep*r_mesh(ir)*vpsum(klp)*bjl(ir,lp)
            dcx(klp,j)=nstep*r_mesh(ir)*vpsum(klp)*bnl(ir,lp)
         enddo
!
!        =============================================================
!        update sx and cx at ir with the corrector
!        =============================================================
         do klp=1,kmax_phi
            sxir(klp) = sxirm1(klp) + hfac*( NINE    *dsx(klp,j)      &
                                            +Nineteen*dsx(klp,jm1)    &
                                            -FIVE    *dsx(klp,jm2)    &
                                                     +dsx(klp,jm3) )
            cxir(klp) = cxirm1(klp) + hfac*( NINE    *dcx(klp,j)      &
                                            +Nineteen*dcx(klp,jm1)    &
                                            -FIVE    *dcx(klp,jm2)    &
                                                     +dcx(klp,jm3) )
         enddo
      enddo
      do klp=1,kmax_phi
         sx(ir,klp) = sxir(klp)
         cx(ir,klp) = cxir(klp)
      enddo
      sxirm1 = sxir
      cxirm1 = cxir
   enddo
!
   nullify(kj3_i3, cgnt_i3, r_mesh)
!
   end subroutine solveSCr
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine intkvpr(atom,SolutionType,kl,fmem)
!  ===================================================================
!
!**===================================================================
!* calculate:                                                        *
!*                                                                   *
!*                    L1                                             *
!*  fint = r * sum   C     * int dr1*r1 * K (r,r1;k)*v (r1)*P (r1;k) *
!*             L1,L2  L2,L    r1           l          L1     L2      *
!*                                                                   *
!* where,                                                            *
!*                                                                   *
!*  summation over L1 is from 1 to jmax_pot;                         *
!*  summation over L2 is from 1 to klpp_max;                         *
!*  integration is performed on r1 from 0 to r, if regular solution  *
!*                              or from r to infinity, if irregular  *
!*                                                         solution; *
!*                                                                   *
!*  k = kappa;                                                       *
!*                                                                   *
!*  P = r*(wf_reg(r)*delta_L''L*theta_L1+wfhat_reg(r)), for regular  *
!*   L''                                                   solution  *
!*                                                                   *
!*    = r*(wf_irr(r)*delta_L''L*theta_L1+wfhat_irr(r)), for irregular*
!*                                                         solution  *
!*                                                                   *
!*  delta_L''L = 1,  if L''= L                                       *
!*             = 0,  else                                            *
!*                                                                   *
!*  theta_L1  = 1,  if L1 = 0                                        *
!*            = 0,  else                                             *
!*                                                                   *
!* and,                                                              *
!*                       ^         ^          ^         ^            *
!*    r*K (r,rp;k)*rp = [n (k*r) * j (k*rp) - j (k*r) * n (k*rp)]/k, *
!*       lp               lp        lp         lp        lp          *
!*                                                                   *
!*                                              for non-derivative;  *
!*                                                                   *
!*                        lp    ^         ^            ^             *
!*                    = [---- * n (k*r) - n   (k*r)] * j (k*rp)      *
!*                       r*k     lp        lp+1         lp           *
!*                                                                   *
!*                        lp    ^         ^            ^             *
!*                     -[---- * j (k*r) - j   (k*r)] * n (k*rp),     *
!*                       r*k     lp        lp+1         lp           *
!*                                                                   *
!*                                              for derivative       *
!*                                                                   *
!**===================================================================
   use IntegrationModule, only : calIntegration
   use DerivativeModule, only : derv5
!
   use LdaCorrectionModule, only : checkLdaCorrection, getPotentialCorrection
!
   implicit   none
!
   integer (kind=1), intent(in) :: SolutionType
!
   integer (kind=IntKind), intent(in) :: atom, kl
!
   logical :: isJl1One
   logical :: isLargeL
!
   integer (kind=IntKind) :: klp
   integer (kind=IntKind) :: klps
   integer (kind=IntKind) :: klpp
   integer (kind=IntKind) :: jl1
   integer (kind=IntKind) :: kl1
   integer (kind=IntKind) :: kls, jls, mls, mp
   integer (kind=IntKind) :: l1, m1, mlv, klv, jlv
   integer (kind=IntKind) :: i1, i2, i1_save
   integer (kind=IntKind) :: lp, ll, lv, ls, lpp
   integer (kind=IntKind) :: ir, ir0
   integer (kind=IntKind) :: it
   integer (kind=IntKind) :: mode
   integer (kind=IntKind) :: nit
   integer (kind=IntKind) :: kmax_pot
   integer (kind=IntKind) :: kmax_step
   integer (kind=IntKind) :: kmax_phi0
   integer (kind=IntKind) :: kmax_phi1
   integer (kind=IntKind) :: kstep, d_csmt
   integer (kind=IntKind) :: klp_0, klp_1
!
   integer (kind=IntKind), parameter :: nit_kk = 13
!
   integer (kind=IntKind), pointer :: flags(:)
   integer (kind=IntKind), pointer :: flags_trunc(:)
!
   integer (kind=IntKind) :: nj0_i3
   integer (kind=IntKind) :: nj3_i3
   integer (kind=IntKind), pointer :: kj0_i3(:)
   integer (kind=IntKind), pointer :: kj3_i3(:)
!
   real (kind=RealKind), pointer :: cgnt0_i3(:)
   real (kind=RealKind), pointer :: cgnt_i3(:)
!
!
   real (kind=RealKind) :: r2lp
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), intent(inout), target :: fmem(iend,5)
!
   complex (kind=CmplxKind), pointer :: plhat(:,:), plhat_old(:,:)
   complex (kind=CmplxKind), pointer :: qlhat(:,:)
   complex (kind=CmplxKind), pointer :: pl(:)
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
   complex (kind=CmplxKind), pointer :: pot_trunc_l(:,:)
   complex (kind=CmplxKind), pointer :: f1mem(:)
   complex (kind=CmplxKind), pointer :: f2mem(:)
   complex (kind=CmplxKind), pointer :: f3mem(:)
   complex (kind=CmplxKind), pointer :: f4mem(:)
   complex (kind=CmplxKind), pointer :: f5mem(:)
   complex (kind=CmplxKind), pointer :: pv1(:), pv2(:), pv3(:), pv4(:), pv5(:), pv0(:)
   complex (kind=CmplxKind), pointer :: pv11(:), pv12(:), pv31(:), pv32(:)
   complex (kind=CmplxKind), pointer :: pv41(:), pv42(:), pv51(:), pv52(:)
   complex (kind=CmplxKind), pointer :: pvr1(:), pvr2(:), pvr3(:), pvr4(:)
!
   complex (kind=CmplxKind) :: cfac, cfacp
   complex (kind=CmplxKind) :: vcorr
!
   f1mem=>fmem(1:iend,1)
   f2mem=>fmem(1:iend,2)
   f3mem=>fmem(1:iend,3)
   f4mem=>fmem(1:iend,4)
   f5mem=>fmem(1:iend,5)
!
   plhat_old => aliasArray2_c(wks_plhat,iend,kmax_phi)
!
   if (SolutionType.eq.RegularSolution) then
      mode=0
      pl=>pl_reg(1:iend)
      plhat=>plhat_reg(1:iend,1:kmax_phi)
      qlhat=>qlhat_reg(1:iend,1:kmax_phi)
      kmax_phi0 = 1
      kmax_phi1 = kmax_phi
      kstep = 1
      klp_0 = 1
      klp_1 = kmax_phi
!
      nit = nit_kk+7
!
   else if (SolutionType.eq.IrregularSolution) then
      mode=1
      pl=>pl_irr(1:iend)
      plhat=>plhat_irr(1:iend,1:kmax_phi)
      qlhat=>qlhat_irr(1:iend,1:kmax_phi)
      kmax_phi0 = 1
      kmax_phi1 = kmax_phi
      kstep = 1
      klp_0 = 1
      klp_1 = kmax_phi
!
      nit = nit_kk+13
!
   else
      i1=SolutionType
!     ----------------------------------------------------------------
      call ErrorHandler('INTKVPR','Unknown SolutionType',i1)
!     ----------------------------------------------------------------
   endif
!
   if ( isPotZeroOutsideRmt ) then
      kmax_phi0 = kl
      kmax_phi1 = kl
!      nit = 10
   endif
!
   if ( isDiagSolv ) then
      klp_0 = kl
      klp_1 = kl
!      nit = 10
   endif
!
   r_mesh=>Grid%r_mesh(1:iend)
!
   d_csmt = numrs_cs-numrs_mt+1
!
   ll = lofk(kl)
   if (checkLdaCorrection(LocalIndex,atom)) then
      m1 = mofk(kl)
      vcorr = getPotentialCorrection(LocalIndex,atom,spin_index,lp,m1)
   else
      vcorr = CZERO
   endif
!
!  ===================================================================
!  Start iterations to solve the integral equations for solution: plhat 
!  Note: nit is the total number of iterations.
!        It is found that for kl=1, nit=10 is necessary to get
!        converged irregular solutions. In other cases, nit=5 should
!        be enough to get converged regular and irregular solutions.
!  ===================================================================
   flags => flags_jl
   pot_l => pot_jl
   kmax_pot = (lmax_pot+1)*(lmax_pot+1)
   kmax_step = (lmax_step+1)*(lmax_step+1)
!
   if ( SSSMethod==2 ) then
      flags_trunc => flags_trunc_jl
      pot_trunc_l => pot_trunc_jl
   endif
!
   do it=1,nit
      plhat_old(1:numrs_cs,1:kmax_phi) = plhat(1:numrs_cs,1:kmax_phi)
      LoopKlp: do klp = klp_0,klp_1
         lp=lofk(klp)
         mp=mofk(klp)
         klps = klp-2*mp
         f1mem=CZERO
!
         LoopKlpp1: do klpp = kmax_phi0,kmax_phi1,kstep
            pv2=>f2mem(1:numrs_cs)
            pv5=>plhat_old(1:numrs_cs,klpp)
!
            pv2 = pv5
!
            if (klpp.eq.kl) then
               pv4 => pl(1:numrs_cs)
!
               pv2 => f4mem(1:numrs_cs)
               pv2 = pv5 + pv4
!
            endif
!
            nj3_i3 = nj3(klpp,klp)
            kj3_i3 => kj3(1:nj3_i3,klpp,klp)
            cgnt_i3 => cgnt(1:nj3_i3,klpp,klp)
!
            pv31 => f3mem(1:numrs_mt)
            pv32 => f3mem(numrs_mt+1:numrs_cs)
            pv3  => f3mem(1:numrs_cs)
            pv3  =  CZERO
!
            isJl1One = .false.
            LOOP_i2_1: do i1 = nj3_i3,1,-1
               kl1 = kj3_i3(i1)
               jl1 = jofk(kl1)
!
               if ( jl1 > jmax_pot_solver ) then
                  cycle LOOP_i2_1
               else if ( jl1>jmax_pot .and. jl1>jmax_trunc ) then
                  cycle LOOP_i2_1
               else if ( jl1<=jmax_pot ) then
                  if ( flags(jl1)==0 .and. flags_trunc(jl1)==0 ) then
                     cycle LOOP_i2_1
                  endif
               else if ( jl1<=jmax_trunc ) then
                  if ( flags_trunc(jl1)==0 ) then
                     cycle LOOP_i2_1
                  endif
               endif
!
               m1=mofk(kl1)
!
               cfac=cgnt_i3(i1)
               cfacp=cfac*m1m(m1)
!
               pv41 => pot_l(1:numrs_mt,jl1)
               if (jl1 == 1) then
!                  pv31 = pv31+cfac*(pv41+(vcorr+PotShift)/Y0)
!                  pv32 = pv32+cfac*(pv42+(vcorr+PotShift)/Y0)
                  isJl1One = .true.
                  i1_save = i1
               else if (m1 >= 0) then
                  if ( jl1<=jmax_pot ) then
                     if ( flags(jl1) /=0 ) then
                        pv31 = pv31 + cfac*pv41
                     endif
                  endif
                  if ( flags_trunc(jl1) /=0 ) then
                     pv42 => pot_trunc_l(2:d_csmt,jl1)
                     pv32 = pv32 + cfac*pv42
                  endif
               else
                  if ( jl1<=jmax_pot ) then
                     if ( flags(jl1) /=0 ) then
                        pv31 = pv31 + cfacp*conjg(pv41)
                     endif
                  endif
                  if ( flags_trunc(jl1) /=0 ) then
                     pv42 => pot_trunc_l(2:d_csmt,jl1)
                     pv32 = pv32 + cfacp*conjg(pv42)
                  endif
               endif
            enddo LOOP_i2_1
!
            if (klpp/=kl) then
               pv1 => f1mem(1:numrs_cs)
!
               pv1 = pv1+pv3*pv2
!
               if ( isJl1One ) then
                  cfac=cgnt_i3(i1_save)
                  pv41 => pot_l(1:numrs_mt,1)
                  pv42 => pot_trunc_l(2:d_csmt,1)
                  pv51 => plhat_old(1:numrs_mt,klpp)
                  pv52 => plhat_old(numrs_mt+1:numrs_cs,klpp)
                  pv11 => f1mem(1:numrs_mt)
                  pv12 => f1mem(numrs_mt+1:numrs_cs)
!
                  pv11 = pv11+cfac*pv51*(pv41+(vcorr+PotShift)/Y0)
                  pv12 = pv12+cfac*pv52*(pv42+(vcorr+PotShift)/Y0)
!
                  isJl1One = .false.
               endif
            else
               pv2 = pv3*pv2
               if ( isJl1One ) then
                  cfac=cgnt_i3(i1_save)
                  pv41 => pot_l(1:numrs_mt,1)
                  pv42 => pot_trunc_l(2:d_csmt,1)
                  pv51 => plhat_old(1:numrs_mt,klpp)
                  pv52 => plhat_old(numrs_mt+1:numrs_cs,klpp)
                  pv11 => f4mem(1:numrs_mt)
                  pv12 => f4mem(numrs_mt+1:numrs_cs)
!
                  pv11 = pv11+cfac*pv51*(pv41+(vcorr+PotShift)/Y0)
                  pv12 = pv12+cfac*pv52*(pv42+(vcorr+PotShift)/Y0)
!
                  isJl1One = .false.
               endif
            endif
         enddo LoopKlpp1
         pv1 = pv1+f4mem(1:numrs_cs)
!
         pv1 => f1mem(1:numrs_cs)
         pvr1 => bjl(1:numrs_cs,lp)
!
         pv2 => f2mem(1:numrs_cs)
         pv2 = pv1*pvr1
!
         pvr2 => bnl(1:numrs_cs,lp)
!
         pv3 => f3mem(1:numrs_cs)
         pv3 = pv1*pvr2
!
!        -------------------------------------------------------------
         call calIntegration( mode, numrs_cs, r_mesh, f2mem, f1mem )
!        -------------------------------------------------------------
         call calIntegration( mode, numrs_cs, r_mesh, f3mem, f4mem )
!        -------------------------------------------------------------
!
         pv4 => f4mem(1:numrs_cs)
!
         if ( mode == 0 ) then
            pv1 = pv1 + HALF*f2mem(1)*r_mesh(1)
            pv4 = pv4 + HALF*f3mem(1)*r_mesh(1)
            f1mem(numrs_cs+1:iend) = CZERO
            f4mem(numrs_cs+1:iend) = CZERO
         else
            f1mem(numrs_cs+1:iend) = CZERO
            f4mem(numrs_cs+1:iend) = CZERO
         endif
!
         pv1 => f1mem(1:iend)
         pv4 => f4mem(1:iend)
         pvr1 => bnl(1:iend,lp)
         pvr2 => bjl(1:iend,lp)
         pv5 => plhat(1:iend,klp)
!
         do ir = 1,numrs_cs
            pv5(ir) = (pvr1(ir)*pv1(ir)-pvr2(ir)*pv4(ir))/kappa
         enddo
         pv5(numrs_cs+1:iend) = CZERO
!
         if ( it.eq.nit_kk ) then
!           pv1 => qlhat(1:numrs_cs,klp)
!           call derv5(pv5,pv1,r_mesh,numrs_cs)
!           pv1(numrs_cs+1:iend) = CZERO
            cfac = lp/kappa
            do ir = 1,numrs_cs
               qlhat(ir,klp)=                                          &
                  ((cfac*bnl(ir,lp)/r_mesh(ir)-bnl(ir,lp+1))*f1mem(ir)-&
                   (cfac*bjl(ir,lp)/r_mesh(ir)-bjl(ir,lp+1))*f4mem(ir))
            enddo
         endif
      enddo LoopKlp
   enddo
!
   nullify( kj3_i3, cgnt_i3, plhat_old )
!
   end subroutine intkvpr
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSCMatrixVolInt(atom,kl,sin_mat,cos_mat,fmem)
!  ===================================================================
   use IntegrationModule, only : calIntegration
!
   use LdaCorrectionModule, only : checkLdaCorrection, getPotentialCorrection
!
   use PotentialModule, only : getTruncatedPotential, getTruncatedPotComponentFlag
!
   use IsoparametricIntegrationModule, only : getVolumeIntegral
!
   implicit none
!
   character (len=16), parameter :: sname='calSCMatrixVolInt'
!
   integer (kind=IntKind), intent(in) :: atom, kl
!
   logical :: isJl1One
!
   integer (kind=IntKind) :: klp, klpp, klppp, klppp_max, lppp_max
   integer (kind=IntKind) :: jl1
   integer (kind=IntKind) :: kl1
   integer (kind=IntKind) :: m1, l1
   integer (kind=IntKind) :: i1, i2, i1_save, kmax_pot, kmax_step
   integer (kind=IntKind) :: lp, lpp, lpot, kpot
   integer (kind=IntKind) :: klv, lv, mlv, jlv, kls, ls, mls, jls
   integer (kind=IntKind) :: ir, npout, nr, nr0
   integer (kind=IntKind), pointer :: potflag(:)
   integer (kind=IntKind), allocatable :: WFflag(:,:)
   integer (kind=IntKind), allocatable :: tmp_potflag(:)
!
   integer (kind=IntKind) :: nj3_i3, nj0_i3, jmax_pot_internal
   integer (kind=IntKind), pointer :: kj0_i3(:)
   integer (kind=IntKind), pointer :: kj3_i3(:)
!
!   real (kind=RealKind) :: sl_abs, cl_abs
   real (kind=RealKind), pointer :: cgnt0_i3(:)
   real (kind=RealKind), pointer :: cgnt_i3(:)
   real (kind=RealKind), pointer :: pr_mesh(:)
!
   complex (kind=CmplxKind) :: vcorr, s_tmp, c_tmp
   complex (kind=CmplxKind) :: cfac, cfacp
!
   complex (kind=CmplxKind), intent(out) :: sin_mat(kmax_int)
   complex (kind=CmplxKind), intent(out) :: cos_mat(kmax_int)
!
   complex (kind=CmplxKind), target :: fmem(iend,4)
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
   complex (kind=CmplxKind), pointer :: tmp_pot(:,:)
   complex (kind=CmplxKind), pointer :: tmp_phi(:,:)
   complex (kind=CmplxKind), pointer :: tmp_cf(:,:)
   complex (kind=CmplxKind), pointer :: tmp_cvs(:,:)
!
   complex (kind=CmplxKind), pointer :: f1mem(:)
   complex (kind=CmplxKind), pointer :: f2mem(:)
   complex (kind=CmplxKind), pointer :: f3mem(:)
   complex (kind=CmplxKind), pointer :: f4mem(:)
   complex (kind=CmplxKind), pointer :: pv1(:), pv2(:), pv3(:), pv4(:)
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   pr_mesh => Grid%r_mesh
!
!  ===================================================================
!  calculate the volume integration over the inscribing sphere which is
!  converted into the surface integration.
!  ===================================================================
   sin_mat(1:kmax_int)=CZERO
   cos_mat(1:kmax_int)=CZERO
!
   nr = numrs_mt
   do klp = 1,kmax_int
      lp=lofk(klp)
      sin_mat(klp) = ( bjl(nr,lp)*qlhat_reg(nr,klp) -                 &
                       kappa*dbjl(nr,lp)*plhat_reg(nr,klp) )
      cos_mat(klp) = ( bnl(nr,lp)*qlhat_reg(nr,klp) -                 &
                       kappa*dbnl(nr,lp)*plhat_reg(nr,klp) )
   enddo
!
   if ( isPotZeroOutsideRmt) then
      return
   endif
!
   npout = numrs_cs-numrs_mt+1
!
   if ( isIsoparam ) then
      vcorr = CZERO
      pot_l   => pot_jl
      potflag => flags_jl
      jmax_pot_internal = jmax_pot
      lpot = lofj(jmax_pot)
      kpot = (lpot+1)**2
      p1c => wks_scvphi(1:kpot*npout,1)
      tmp_pot => aliasArray2_c( p1c,npout,kpot)
      p1c => wks_scvphi(1:kmax_phi*npout,2)
      tmp_phi => aliasArray2_c( p1c,npout,kmax_phi)
      allocate(tmp_potflag(kpot), WFflag(kmax_phi,kmax_int))
      do klp=1,kpot
         m1=mofk(klp)
         jl1=jofk(klp)
         tmp_potflag(klp) = potflag(jl1)
         do ir =1,npout
            if(klp==1) then
               tmp_pot(ir,klp) = (pot_l(numrs_mt+ir-1,jl1)+(vcorr+PotShift)/Y0)
            else
               if (m1>=0) then
                  tmp_pot(ir,klp) = pot_l(numrs_mt+ir-1,jl1)
               else
                  tmp_pot(ir,klp) = m1m(m1)*conjg(pot_l(numrs_mt+ir-1,jl1))
               endif
            endif
         enddo
      enddo
      do klp=1,kmax_phi
         do ir =1,npout
            tmp_phi(ir,klp) = plhat_reg(ir+numrs_mt-1,klp)/(pr_mesh(ir+numrs_mt-1)**2)
         enddo
      enddo
      nr0 = numrs_mt
      LoopKlp1: do klp = 1,kmax_int
         lp=lofk(klp)
         sin_mat(klp) = sin_mat(klp) + getVolumeIntegral( LocalIndex, npout,                 &
                                                          pr_mesh(numrs_mt:numrs_cs),        &
                                                          klp,bjl(numrs_mt:numrs_cs,lp),     &
                                                          jmax_pot,kpot,tmp_pot,tmp_potflag, &
                                                          jofk(kmax_phi),kmax_phi,tmp_phi,   &
                                                          WFflag(1:kmax_phi,kl))
         cos_mat(klp) = cos_mat(klp) + getVolumeIntegral( LocalIndex, npout,                 &
                                                          pr_mesh(numrs_mt:numrs_cs),        &
                                                          klp,bnl(numrs_mt:numrs_cs,lp),     &
                                                          jmax_pot,kpot,tmp_pot,tmp_potflag, &
                                                          jofk(kmax_phi),kmax_phi,tmp_phi,   &
                                                          WFflag(1:kmax_phi,kl))
      enddo LoopKlp1
      deallocate(tmp_potflag, WFflag)
   else
      f1mem=>fmem(1:npout,1)
      f2mem=>fmem(1:npout,2)
      f3mem=>fmem(1:npout,3)
      f4mem=>fmem(1:npout,4)
!
      if (SSSMethod == 1) then
         jmax_pot_internal = jmax_trunc
         pot_l   => getTruncatedPotential(LocalIndex,atom,spin_index)
         potflag => getTruncatedPotComponentFlag(LocalIndex)
      else
         jmax_pot_internal = jmax_pot
         potflag => flags_jl
         pot_l => pot_jl
      endif
      nr0 = 1
!
      if (checkLdaCorrection(LocalIndex,atom)) then
         vcorr = getPotentialCorrection(LocalIndex,atom,spin_index,lofk(kl),mofk(kl))
      else
         vcorr = CZERO
      endif
!
      lv = lofj(jmax_pot)
      kmax_pot = (lv+1)**2
      kmax_step = (lmax_step+1)*(lmax_step+1)
      lppp_max = min(2*lmax_phi,lmax_step+lv) ! It is the internal L's upper limit for double Gaunt factors case
      klppp_max = (lppp_max+1)**2
!
      LoopKlp0: do klp = 1,kmax_int
         lp=lofk(klp)
         f1mem=CZERO
         f3mem=CZERO
         if (SSSMethod == 1) then
            do klpp = kmax_phi,1,-1
               lpp = lofk(klpp)
               f2mem=CZERO
               nj3_i3 = nj3(klpp,klp)
               kj3_i3 => kj3(1:nj3_i3,klpp,klp)
               cgnt_i3 => cgnt(1:nj3_i3,klpp,klp)
!
               isJl1One = .false.
               LOOP_i1: do i1 = nj3_i3,1,-1
                  kl1=kj3_i3(i1)
                  jl1=jofk(kl1)
!
                  if ( jl1>jmax_pot_internal .or. (jl1>jmax_pot_solver   &
                       .and. jmax_pot_solver>0) ) then
                     cycle LOOP_i1
                  else if (potflag(jl1) == 0) then
                     cycle LOOP_i1
                  endif
!
                  m1=mofk(kl1)
                  cfac=cgnt_i3(i1)
                  cfacp=cfac*m1m(m1)
!
                  if (jl1 == 1) then
                     isJl1One = .true.
                     i1_save = i1
                  else if (m1 >= 0) then
                     do ir=1,npout
                        f2mem(ir) = f2mem(ir) + cfac*pot_l(ir+nr0-1,jl1)
                     enddo
                  else
                     do ir=1,npout
                        f2mem(ir) = f2mem(ir) +                          &
                                    cfacp*conjg(pot_l(ir+nr0-1,jl1))
                     enddo
                  endif
               enddo LOOP_i1
               if ( isJl1One ) then
                  cfac = cgnt_i3(i1_save)
                  do ir=1,npout
                     f2mem(ir) = f2mem(ir) + &
                                 cfac*(pot_l(ir+nr0-1,1)+(vcorr+PotShift)/Y0)
                  enddo
                  isJl1One = .false.
               endif
               do ir=1,npout
                  f1mem(ir) = f1mem(ir) +                                &
                              f2mem(ir)*plhat_reg(ir+numrs_mt-1,klpp)
               enddo
            enddo
         else ! This method involves the products of double Gaunt factors
!           ==========================================================
!           Note: tmp_cf = sum [gaunt_factor * reg_sol]
!           ==========================================================
            p1c => wks_scvphi(1:klppp_max*npout,1)
            tmp_cf => aliasArray2_c(p1c,npout,klppp_max)
            tmp_cf = CZERO
            do klpp = kmax_phi,1,-1
               nj3_i3 = nj3(klpp,klp)
               kj3_i3 => kj3(1:nj3_i3,klpp,klp)
               cgnt_i3 => cgnt(1:nj3_i3,klpp,klp)
               do i1 = 1, nj3_i3
                  klppp = kj3_i3(i1)
                  if (klppp <= klppp_max) then
                     do ir=1,npout
                        tmp_cf(ir,klppp) = tmp_cf(ir,klppp) + cgnt_i3(i1)*plhat_reg(ir+numrs_mt-1,klpp)
                     enddo
                  endif
               enddo
            enddo
!           ==========================================================
!           Note: cpotstep = sum [gaunt_factor * pot * step]
!                 f1mem = sum [gaunt_factor * gaunt_factor * pot * step * reg_sol]
!           ==========================================================
            do klppp = 1, klppp_max
               do ir = 1, npout
                  f1mem(ir) = f1mem(ir) + tmp_cf(ir,klppp) * cpotstep(ir,klppp)
               enddo
            enddo
         endif
!
         do ir=1,npout
            f3mem(ir) = bjl(ir+numrs_mt-1,lp)*f1mem(ir)
            f4mem(ir) = bnl(ir+numrs_mt-1,lp)*f1mem(ir)
         enddo
!
         f1mem = CZERO
!        -------------------------------------------------------------
!         call calIntegration( npout, r_mesh, f3mem, f1mem, 0 )
         call calIntegration( 0, npout, pr_mesh(numrs_mt:numrs_cs), f3mem, f1mem )
!        -------------------------------------------------------------
         sin_mat(klp) = sin_mat(klp) + f1mem(npout)
!
         f1mem = CZERO
!        -------------------------------------------------------------
!         call calIntegration( npout, r_mesh, f4mem, f1mem, 0 )
         call calIntegration( 0, npout, pr_mesh(numrs_mt:numrs_cs), f4mem, f1mem )
!        -------------------------------------------------------------
         cos_mat(klp) = cos_mat(klp) + f1mem(npout)
!
      enddo LoopKlp0
      nullify( kj3_i3, cgnt_i3, pot_l)
      nullify(f1mem, f2mem, f3mem, f4mem)
   endif
!
   end subroutine calSCMatrixVolInt
!  ===================================================================
!
!  *******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSCMatrix(kl,sin_mat,cos_mat)
!  ====================================================================
   use InterpolationModule, only : PolyInterp
!
   use BesselModule, only : SphericalBessel
   use BesselModule, only : SphericalNeumann
!
   use SurfElementsModule, only : getNumGaussPoints
   use SurfElementsModule, only : getNumDiffSurfR
   use SurfElementsModule, only : getGaussYlm, getGaussGradYlmDotNvec
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kl
!
   complex (kind=CmplxKind), intent(out) :: sin_mat(kmax_int)
   complex (kind=CmplxKind), intent(out) :: cos_mat(kmax_int)
!
   character (len=12), parameter :: sname='calSCMatrix'
!
   integer (kind=IntKind) :: klp, klc, klpi, klpf
   integer (kind=IntKind) :: l,lp, mp
   integer (kind=IntKind) :: i, mfac, jsurf, j_inter
   integer (kind=IntKind), parameter :: n_inter=7
!
   real (kind=RealKind), pointer :: xp(:)
   real (kind=RealKind) :: err_est
!   real (kind=RealKind) :: sl_abs, cl_abs
!
   complex (kind=CmplxKind) :: philm, nDotGradPhi
   complex (kind=CmplxKind) :: JYcc, NYcc, nDotGradJYcc, nDotGradNYcc
!   complex (kind=CmplxKind) :: smix,cmix,dsmix,dcmix
   complex (kind=CmplxKind) :: nfac, pcm0
   complex (kind=CmplxKind) :: s_tmp, c_tmp
   complex (kind=CmplxKind) :: c_sphere_mt
   complex (kind=CmplxKind) :: s_sphere_mt
!
   complex (kind=CmplxKind), pointer :: pl_surf(:,:)
   complex (kind=CmplxKind), pointer :: ql_surf(:,:)
   complex (kind=CmplxKind), pointer :: yp(:)
   complex (kind=CmplxKind), pointer :: pv1(:),pv2(:),pv3(:),pv4(:)
!
!  ===================================================================
!  interpolate pl and ql on the surface.
!
!  Note that,
!
!                                   r   d
!     p (r) = r*phi (r),    q (r) = - * -- phi (r)
!      l           l         l      M   dr    l
!                                    0
!
!  ===================================================================
   NumSurfPoints = getNumGaussPoints(LocalIndex)
   NumDiffSurfR = getNumDiffSurfR(LocalIndex)
   ylm => getGaussYlm(LocalIndex)
   nDotGradY => getGaussGradYlmDotNvec(LocalIndex)
!
   pl_surf=>aliasArray2_c(wksl_plsurf,kmax_phi,NumDiffSurfR)
   ql_surf=>aliasArray2_c(wksl_qlsurf,kmax_phi,NumDiffSurfR)
!
   pl_surf(1:kmax_phi,1:NumDiffSurfR)=CZERO
   ql_surf(1:kmax_phi,1:NumDiffSurfR)=CZERO
!
!  ===================================================================
   if (jmax_pot == 1 .and. isSphPotZeroOutsideRmt) then
      klpi=kl
      klpf=kl
   else
      klpi=1
      klpf=kmax_phi
   endif
!
   do jsurf=1,NumDiffSurfR
      j_inter=JHunt(jsurf)
      j_inter=min(max(j_inter - (n_inter-1)/2,1), iend+1-n_inter)
      xp=>Grid%r_mesh(j_inter:j_inter+n_inter-1)
      if (jmax_pot == 1 .and. isSphPotZeroOutsideRmt) then
         yp=>cm0(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,DiffSurfR(jsurf),pcm0,err_est)
!        -------------------------------------------------------------
         l = lofk(kl)
         s_sphere_mt = s_sph_mt(l)
         c_sphere_mt = c_sph_mt(l)
         pl_surf(kl,jsurf) = ( s_sphere_mt*snl(l,jsurf)               &
                          -c_sphere_mt*sjl(l,jsurf) )*DiffSurfR(jsurf)
         ql_surf(kl,jsurf) = ( s_sphere_mt*dsnl(l,jsurf)              &
                          -c_sphere_mt*dsjl(l,jsurf) )*DiffSurfR(jsurf)/pcm0
      else
         LoopKlp2: do klp=klpi,klpf
            yp=>plhat_reg(j_inter:j_inter+n_inter-1,klp)
!           ----------------------------------------------------------
            call PolyInterp(n_inter,xp,yp,DiffSurfR(jsurf),           &
                            pl_surf(klp,jsurf),err_est)
!           ----------------------------------------------------------
            yp=>qlhat_reg(j_inter:j_inter+n_inter-1,klp)
!           ----------------------------------------------------------
            call PolyInterp(n_inter,xp,yp,DiffSurfR(jsurf),           &
                            ql_surf(klp,jsurf),err_est)
!           ----------------------------------------------------------
         enddo LoopKlp2
      endif
!
!!    l = lofk(kl)
!!    if (abs(pl_surf(kl,jsurf)-(s_sphere_mt*snl(l,jsurf)               &
!!           -c_sphere_mt*sjl(l,jsurf) )*DiffSurfR(jsurf)) > TEN2m6) then
!!       write(6,'(a,2i5,2d15.6,2x,2d15.6)')'kl,jsurf,pl = ',kl,jsurf,  &
!!            pl_surf(kl,jsurf),(s_sphere_mt*snl(l,jsurf)               &
!!                        -c_sphere_mt*sjl(l,jsurf))*DiffSurfR(jsurf)
!!    endif
!!    if (abs(ql_surf(kl,jsurf)-(s_sphere_mt*dsnl(l,jsurf)              &
!!           -c_sphere_mt*dsjl(l,jsurf) )*DiffSurfR(jsurf)) > TEN2m6) then
!!       write(6,'(a,2i5,2d15.6,2x,2d15.6)')'kl,jsurf,ql = ',kl,jsurf,  &
!!            ql_surf(kl,jsurf), ( s_sphere_mt*dsnl(l,jsurf)            &
!!                        -c_sphere_mt*dsjl(l,jsurf) )*DiffSurfR(jsurf)
!!    endif
!!    pl_surf(kl,jsurf) = ( s_sphere_mt*snl(l,jsurf)                    &
!!                         -c_sphere_mt*sjl(l,jsurf) )*DiffSurfR(jsurf)
!!    ql_surf(kl,jsurf) = ( s_sphere_mt*dsnl(l,jsurf)                   &
!!                         -c_sphere_mt*dsjl(l,jsurf))*DiffSurfR(jsurf)
!
   enddo
   nullify(xp, yp)
!
!  ===================================================================
!  integration of function over the triangle.
!  loop over the points on the boundary planes that defines the 
!  Voronoi Polyhedron.
!  ===================================================================
   sin_mat(1:kmax_int)=CZERO
   cos_mat(1:kmax_int)=CZERO
   l = lofk(kl)
   do i=1,NumSurfPoints
!     ================================================================
!     calculate Phi, Grad{Phi}, for (l, m), and store them in 
!     philm, and nDotGradPhi.
!     ================================================================
      jsurf=GaussR2SurfR(i)
      philm=CZERO
      nDotGradPhi=CZERO
      pv1=>pl_surf(klpi:klpf,jsurf)
      pv2=>ql_surf(klpi:klpf,jsurf)
      pv3=>ylm(klpi:klpf,i)
      pv4=>nDotGradY(klpi:klpf,i)
      philm=sum(pv1*pv3)
      nDotGradPhi=sum(pv2*pv3*n_dot_r(i)+pv4*pv1)
!      LoopKlp1: do klp = klpf, klpi, -1
!         philm = philm+pl_surf(klp,jsurf)*ylm(klp,i)
!         nDotGradPhi = nDotGradPhi +                                  &
!                            ql_surf(klp,jsurf)*n_dot_r(i)*ylm(klp,i)+ &
!                            nDotGradY(klp,i)*pl_surf(klp,jsurf))
!      enddo LoopKlp1
      philm=philm/r_surf(i)/dnlrmt(l)
      nDotGradPhi=nDotGradPhi/r_surf(i)/dnlrmt(l)
! 
!     ================================================================
!     calculate sin_mat and cos_mat matrices
!     ================================================================
      LoopKlp: do klp = 1,kmax_int
         lp=lofk(klp)
         mp=mofk(klp)
         klc=klp-2*mp
         mfac=m1m(mp)
         JYcc=mfac*sjl(lp,jsurf)*ylm(klc,i)
         NYcc=mfac*snl(lp,jsurf)*ylm(klc,i)/dnlrmt(lp)
         nDotGradJYcc=mfac*(n_dot_r(i)*dsjl(lp,jsurf)*ylm(klc,i)+     &
                            sjl(lp,jsurf)*nDotGradY(klc,i))
         nDotGradNYcc=mfac*(n_dot_r(i)*dsnl(lp,jsurf)*ylm(klc,i)+     &
                            snl(lp,jsurf)*nDotGradY(klc,i))/dnlrmt(lp)
         sin_mat(klp)=sin_mat(klp)+                                   &
                      weight(i)*(JYcc*nDotGradPhi-nDotGradJYcc*philm)
         cos_mat(klp)=cos_mat(klp)+                                   &
                      weight(i)*(NYcc*nDotGradPhi-nDotGradNYcc*philm)
      enddo LoopKlp
   enddo
   nfac = kappa*dnlrmt(l)
   sin_mat=nfac*sin_mat
   do klp=1,kmax_int
      lp=lofk(klp)
      cos_mat(klp)=nfac*cos_mat(klp)*dnlrmt(lp)
   enddo
!
   if ( real(-SQRTm1*energy,kind=RealKind) == ZERO ) then
      if ( real(energy,kind=RealKind)<0 ) then
         do klp = 1,kmax_int
            s_tmp = sin_mat(klp)
            c_tmp = cos_mat(klp)
            sin_mat(klp) = cmplx(ZERO,real(-SQRTm1*s_tmp,kind=RealKind),kind=CmplxKind)
            cos_mat(klp) = cmplx(real(c_tmp,kind=RealKind),ZERO,kind=CmplxKind)
         enddo
      else
         do klp = 1,kmax_int
            s_tmp = sin_mat(klp)
            c_tmp = cos_mat(klp)
            sin_mat(klp) = cmplx(real(s_tmp,kind=RealKind),ZERO,kind=CmplxKind)
            cos_mat(klp) = cmplx(real(c_tmp,kind=RealKind),ZERO,kind=CmplxKind)
         enddo
      endif
   endif
!
!  ===================================================================
!  make corrections to sin_mat and cos_mat matrices
!  ===================================================================
!! allocate(sin_cor(1:kmax_int),cos_cor(1:kmax_int))
!! sin_cor(1:kmax_int)=CZERO
!! cos_cor(1:kmax_int)=CZERO
!! do i=1,NumSurfPoints
!     ================================================================
!     get corrections to sin_mat and cos_mat matrices
!     ================================================================
!!    smix=CZERO; cmix=CZERO
!!    dsmix=CZERO; dcmix=CZERO
!!    jsurf=GaussR2SurfR(i)
!!    do klp=1,kmax_int
!!       lp=lofk(klp)
!!       mp=mofk(klp)
!!       smix=smix+sin_mat(klp)*snl(lp,jsurf)*ylm(klp,i)
!!       cmix=cmix+cos_mat(klp)*sjl(lp,jsurf)*ylm(klp,i)
!!       dsmix=dsmix+sin_mat(klp)*(n_dot_r(i)*dsnl(lp,jsurf)*ylm(klp,i)+  &
!!                                 snl(lp,jsurf)*nDotGradY(klp,i))
!!       dcmix=dcmix+cos_mat(klp)*(n_dot_r(i)*dsjl(lp,jsurf)*ylm(klp,i)+  &
!!                                 sjl(lp,jsurf)*nDotGradY(klp,i))
!!    enddo
!!    philm=CZERO
!!    nDotGradPhi=CZERO
!!    do klp=klpi,klpf
!!       philm=philm+pl_surf(klp,jsurf)*ylm(klp,i)
!!       nDotGradPhi=nDotGradPhi+ql_surf(klp,jsurf)*n_dot_r(i)*ylm(klp,i)+ &
!!                               nDotGradY(klp,i)*pl_surf(klp,jsurf)
!!    enddo
!!    philm=philm/r_surf(i)
!!    nDotGradPhi=nDotGradPhi/r_surf(i)
!     ================================================================
!     re-calculate sin_mat and cos_mat matrices
!     ================================================================
!!    do klp=1,kmax_int
!!       lp=lofk(klp)
!!       mp=mofk(klp)
!!       klc=klp-2*mp
!!       mfac=m1m(mp)
!!       JYcc=mfac*sjl(lp,jsurf)*ylm(klc,i)
!!       NYcc=mfac*snl(lp,jsurf)*ylm(klc,i)
!!       nDotGradJYcc=mfac*(n_dot_r(i)*dsjl(lp,jsurf)*ylm(klc,i)+     &
!!                          sjl(lp,jsurf)*nDotGradY(klc,i))
!!       nDotGradNYcc=mfac*(n_dot_r(i)*dsnl(lp,jsurf)*ylm(klc,i)+     &
!!                          snl(lp,jsurf)*nDotGradY(klc,i))
!!       sin_cor(klp)=sin_cor(klp)+                                   &
!!                    weight(i)*(JYcc*(nDotGradPhi-dsmix)-            &
!!                                     nDotGradJYcc*(philm-smix))
!!       cos_cor(klp)=cos_cor(klp)+                                   &
!!                    weight(i)*(NYcc*(nDotGradPhi+dcmix)-            &
!!                                     nDotGradNYcc*(philm+cmix))
!!    enddo
!! enddo
!! sin_mat(1:kmax_int)=sin_mat(1:kmax_int)-kappa*sin_cor(1:kmax_int)
!! cos_mat(1:kmax_int)=cos_mat(1:kmax_int)-kappa*cos_cor(1:kmax_int)
!! deallocate(sin_cor,cos_cor)
!
   if (stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine calSCMatrix
!  =================================================================== 
!
!  ******************************************************************* 
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
   subroutine calJostMatrix(kl,jost_mat)
!  =================================================================== 
   use InterpolationModule, only : PolyInterp
!
   use SurfElementsModule, only : getNumGaussPoints
   use SurfElementsModule, only : getNumDiffSurfR
   use SurfElementsModule, only : getGaussYlm, getGaussGradYlmDotNvec
!
   implicit none
!
   character (len=12), parameter :: sname='calJostMatrix'
!
   integer (kind=IntKind), intent(in) :: kl
!
   integer (kind=IntKind) :: klp, klc, klpi, klpf
   integer (kind=IntKind) :: l,lp, mp
   integer (kind=IntKind) :: i, mfac, jsurf, j_inter
   integer (kind=IntKind), parameter :: n_inter=5
!
   real (kind=RealKind), pointer :: xp(:)
   real (kind=RealKind) :: err_est
!
   complex (kind=CmplxKind), intent(out) :: jost_mat(kmax_int)
!
   complex (kind=CmplxKind) :: philm, nDotGradPhi
   complex (kind=CmplxKind) :: HYcc, nDotGradHYcc
   complex (kind=CmplxKind) :: shl, dshl, nfac, eoc2p1
!
   complex (kind=CmplxKind), pointer :: pl_surf(:,:)
   complex (kind=CmplxKind), pointer :: ql_surf(:,:)
   complex (kind=CmplxKind), pointer :: yp(:)
!
!  ===================================================================
!  interpolate pl and ql on the surface.
!
!  Note that,
!
!                                 r   d
!     p (r) = r*H (r),    q (r) = - * -- H (r)
!      l         l         l      M   dr  l
!                                  0
!
!  ===================================================================
   NumSurfPoints = getNumGaussPoints(LocalIndex)
   NumDiffSurfR = getNumDiffSurfR(LocalIndex)
   ylm => getGaussYlm(LocalIndex)
   nDotGradY => getGaussGradYlmDotNvec(LocalIndex)
!
   pl_surf=>aliasArray2_c(wksl_plsurf,kmax_phi,NumDiffSurfR)
   ql_surf=>aliasArray2_c(wksl_qlsurf,kmax_phi,NumDiffSurfR)
!
   pl_surf(1:kmax_kkr,1:NumDiffSurfR)=CZERO
   ql_surf(1:kmax_kkr,1:NumDiffSurfR)=CZERO
!
!  ===================================================================
   if (jmax_pot == 1) then
      klpi=kl
      klpf=kl
   else
      klpi=1
      klpf=kmax_kkr
   endif
!
   eoc2p1=CONE+energy*c2inv
   do jsurf=1,NumDiffSurfR
      j_inter=JHunt(jsurf)
      j_inter=min(max(j_inter - (n_inter-1)/2,1), iend+1-n_inter)
      xp=>Grid%r_mesh(j_inter:j_inter+n_inter-1)
      if (jmax_pot == 1 .and. isSphPotZeroOutsideRmt) then
         l = lofk(kl)
         pl_surf(kl,jsurf) = ( sjl(l,jsurf)+SQRTm1*snl(l,jsurf))      &
                            *DiffSurfR(jsurf)
         ql_surf(kl,jsurf) = ( dsjl(l,jsurf)+SQRTm1*dsnl(l,jsurf) )   &
                            *DiffSurfR(jsurf)/eoc2p1
      else
         do klp = klpf, klpi, -1
            yp=>plhat_irr(j_inter:j_inter+n_inter-1,klp)
!           ----------------------------------------------------------
            call PolyInterp(n_inter,xp,yp,DiffSurfR(jsurf),           &
                            pl_surf(klp,jsurf),err_est)
!           ----------------------------------------------------------
            yp=>qlhat_irr(j_inter:j_inter+n_inter-1,klp)
!           ----------------------------------------------------------
            call PolyInterp(n_inter,xp,yp,DiffSurfR(jsurf),           &
                            ql_surf(klp,jsurf),err_est)
!           ----------------------------------------------------------
         enddo
      endif
   enddo
   nullify(xp, yp)
!
!  ===================================================================
!  integration of function over the triangle.
!  loop over the points on the boundary planes that defines the 
!  Voronoi Polyhedron.
!  ===================================================================
   jost_mat(1:kmax_int)=CZERO
   l = lofk(kl)
   do i=1,NumSurfPoints
!     ================================================================
!     calculate Phi, Grad{Phi}, for (l, m), and store them in 
!     philm, and nDotGradPhi.
!     ================================================================
      jsurf=GaussR2SurfR(i)
      philm=CZERO
      nDotGradPhi=CZERO
      do klp = klpf, klpi, -1
         philm=philm+pl_surf(klp,jsurf)*ylm(klp,i)
         nDotGradPhi=nDotGradPhi+ql_surf(klp,jsurf)*n_dot_r(i)*ylm(klp,i)  &
                                +nDotGradY(klp,i)*pl_surf(klp,jsurf)
      enddo
      philm=philm/r_surf(i)/dnlrmt(l)
      nDotGradPhi=nDotGradPhi/r_surf(i)/dnlrmt(l)
! 
!     ================================================================
!     calculate jost_mat matrix
!     ================================================================
      do klp=1,kmax_int
         lp=lofk(klp)
         mp=mofk(klp)
         klc=klp-2*mp
         mfac=m1m(mp)
         shl=(sjl(lp,jsurf)+SQRTm1*snl(lp,jsurf))/dnlrmt(lp)
         dshl=(dsjl(lp,jsurf)+SQRTm1*dsnl(lp,jsurf))/dnlrmt(lp)
         HYcc=mfac*shl*ylm(klc,i)
         nDotGradHYcc=mfac*(n_dot_r(i)*dshl*ylm(klc,i)+shl*nDotGradY(klc,i))
         jost_mat(klp)=jost_mat(klp)+weight(i)*(HYcc*nDotGradPhi-nDotGradHYcc*philm)
      enddo
   enddo
   nfac = SQRTm1*kappa*dnlrmt(l)
   do klp = 1,kmax_int
      lp = lofk(klp)
      jost_mat(klp)=nfac*dnlrmt(lp)*jost_mat(klp)
   enddo
!
   if (stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine calJostMatrix
!  =================================================================== 
!
!  *******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calTMatrix(kl,t_mat)          ! This routine needs more work
!  ====================================================================
   use InterpolationModule, only : PolyInterp
!
   use SurfElementsModule, only : getNumGaussPoints
   use SurfElementsModule, only : getNumDiffSurfR
   use SurfElementsModule, only : getGaussYlm, getGaussGradYlmDotNvec
!
   implicit none
!
   character (len=12), parameter :: sname='calTMatrix'
!
   integer (kind=IntKind), intent(in) :: kl
!
   integer (kind=IntKind) :: klp, klc, klpi, klpf
   integer (kind=IntKind) :: l,lp, mp
   integer (kind=IntKind) :: i, mfac, jsurf, j_inter
   integer (kind=IntKind), parameter :: n_inter=5
!
   real (kind=RealKind), pointer :: xp(:)
   real (kind=RealKind) :: err_est
!
   complex (kind=CmplxKind), intent(out) :: t_mat(kmax_kkr)
!
   complex (kind=CmplxKind) :: philm, nDotGradPhi
   complex (kind=CmplxKind) :: JYcc, nDotGradJYcc, pcm0
   complex (kind=CmplxKind) :: c_sphere_mt
   complex (kind=CmplxKind) :: s_sphere_mt
!
   complex (kind=CmplxKind), pointer :: pl_surf(:,:)
   complex (kind=CmplxKind), pointer :: ql_surf(:,:)
   complex (kind=CmplxKind), pointer :: yp(:)
!
!  ===================================================================
!  interpolate pl and ql on the surface.
!
!  Note that,
!
!                                    r   d
!      p (r) = r*phi (r),    q (r) = - * -- phi (r)
!       l           l         l      M   dr    l
!                                     0
!
!  ===================================================================
   NumSurfPoints = getNumGaussPoints(LocalIndex)
   NumDiffSurfR = getNumDiffSurfR(LocalIndex)
   ylm => getGaussYlm(LocalIndex)
   nDotGradY => getGaussGradYlmDotNvec(LocalIndex)
!
   pl_surf=>aliasArray2_c(wksl_plsurf,kmax_phi,NumDiffSurfR)
   ql_surf=>aliasArray2_c(wksl_qlsurf,kmax_phi,NumDiffSurfR)
!
   pl_surf(1:kmax_phi,1:NumDiffSurfR)=CZERO
   ql_surf(1:kmax_phi,1:NumDiffSurfR)=CZERO
   do jsurf=1,NumDiffSurfR
      j_inter=JHunt(jsurf)
      j_inter=min(max(j_inter - (n_inter-1)/2,1), iend+1-n_inter)
      xp=>Grid%r_mesh(j_inter:j_inter+n_inter-1)
      if (jmax_pot == 1 .and. isSphPotZeroOutsideRmt) then
         yp=>cm0(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,DiffSurfR(jsurf),pcm0,err_est)
!        -------------------------------------------------------------
         l = lofk(kl)
         s_sphere_mt = s_sph_mt(l)
         c_sphere_mt = c_sph_mt(l)
         pl_surf(kl,jsurf) = ( s_sphere_mt*snl(l,jsurf)               &
                          -c_sphere_mt*sjl(l,jsurf) )*DiffSurfR(jsurf)
         ql_surf(kl,jsurf) = ( s_sphere_mt*dsnl(l,jsurf)              &
                          -c_sphere_mt*dsjl(l,jsurf) )*DiffSurfR(jsurf)/pcm0
      else
         do klp=1,kmax_phi
            yp=>plhat_reg(j_inter:j_inter+n_inter-1,klp)
!           ----------------------------------------------------------
            call PolyInterp(n_inter,xp,yp,DiffSurfR(jsurf),           &
                            pl_surf(klp,jsurf),err_est)
!           ----------------------------------------------------------
            yp=>qlhat_reg(j_inter:j_inter+n_inter-1,klp)
!           ----------------------------------------------------------
            call PolyInterp(n_inter,xp,yp,DiffSurfR(jsurf),           &
                            ql_surf(klp,jsurf),err_est)
!           ----------------------------------------------------------
         enddo
      endif
   enddo
   nullify(xp, yp)
!
   if (jmax_pot == 1) then
      klpi=kl
      klpf=kl
   else
      klpi=1
      klpf=kmax_phi
   endif
!
!  ===================================================================
!  integration of function over the triangle.
!  loop over the points on the boundary planes that defines the 
!  Voronoi Polyhedron.
!  ===================================================================
   t_mat(1:kmax_kkr)=CZERO
   l = lofk(kl)
   do i=1,NumSurfPoints
!     ================================================================
!     calculate Phi, Grad{Phi}, for (l, m), and store them in 
!     philm, and nDotGradPhi.
!     ================================================================
      jsurf=GaussR2SurfR(i)
      philm=CZERO
      nDotGradPhi=CZERO
      do klp=klpi,klpf
         philm=philm+pl_surf(klp,jsurf)*ylm(klp,i)
         nDotGradPhi=nDotGradPhi+ql_surf(klp,jsurf)*n_dot_r(i)*ylm(klp,i)+ &
                                 nDotGradY(klp,i)*pl_surf(klp,jsurf)
      enddo
      philm=philm/r_surf(i)/dnlrmt(l)
      nDotGradPhi=nDotGradPhi/r_surf(i)/dnlrmt(l)
! 
!     ================================================================
!     calculate t_mat
!     ================================================================
      do klp=1,kmax_kkr
         lp=lofk(klp)
         mp=mofk(klp)
         klc=klp-2*mp
         mfac=m1m(mp)
         JYcc=mfac*sjl(lp,jsurf)*ylm(klc,i)
         nDotGradJYcc=mfac*(n_dot_r(i)*dsjl(lp,jsurf)*ylm(klc,i)+     &
                            sjl(lp,jsurf)*nDotGradY(klc,i))
         t_mat(klp)=t_mat(klp)+weight(i)*(JYcc*nDotGradPhi-nDotGradJYcc*philm)
      enddo
   enddo
   do klp=1,kmax_kkr
      t_mat(klp)=t_mat(klp)*dnlrmt(l)
   enddo
!
   if (stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine calTMatrix
!  =================================================================== 
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSurfJNY()
!  ===================================================================
   use SurfElementsModule, only : getNumGaussPoints
   use SurfElementsModule, only : getGaussRvecDotNvec, getGaussWght, getGaussR
   use SurfElementsModule, only : getGaussYlm, getGaussGradYlmDotNvec
   use SurfElementsModule, only : getNumDiffSurfR, getDiffSurfR
   use SurfElementsModule, only : getGaussR2SurfR
!
   use BesselModule, only : SphericalBessel
   use BesselModule, only : SphericalNeumann
!
   implicit none
!
   integer (kind=IntKind) :: j
!
   complex (kind=CmplxKind) :: kr
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   NumSurfPoints = getNumGaussPoints(LocalIndex)
   n_dot_r => getGaussRvecDotNvec(LocalIndex)
   weight => getGaussWght(LocalIndex)
   r_surf => getGaussR(LocalIndex)
   ylm => getGaussYlm(LocalIndex)
   nDotGradY => getGaussGradYlmDotNvec(LocalIndex)
!
   NumDiffSurfR = getNumDiffSurfR(LocalIndex)
   DiffSurfR => getDiffSurfR(LocalIndex)
   GaussR2SurfR => getGaussR2SurfR(LocalIndex)
   do j=1,NumDiffSurfR
      JHunt(j)=Grid%jmt
!     ----------------------------------------------------------------
      call hunt(Grid%jend_plus_n, Grid%r_mesh, DiffSurfR(j), JHunt(j))
!     ----------------------------------------------------------------
   enddo
!
   do j=1,NumDiffSurfR
      kr=kappa*DiffSurfR(j)
!     ----------------------------------------------------------------
      call SphericalBessel(lmax_phi,kr,sjl(0:lmax_phi,j),dsjl(0:lmax_phi,j))
!     ----------------------------------------------------------------
      dsjl(0:lmax_phi,j)=kappa*dsjl(0:lmax_phi,j)
!     ----------------------------------------------------------------
      call SphericalNeumann(lmax_phi,kr,snl(0:lmax_phi,j),dsnl(0:lmax_phi,j))
!     ----------------------------------------------------------------
      dsnl(0:lmax_phi,j)=kappa*dsnl(0:lmax_phi,j)
   enddo
!
   end subroutine calSurfJNY
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSineMatrix(spin, site, atom) result(sin_mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getSineMatrix','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   sin_mat => Scatter(id)%Solutions(ic,is)%sin_mat
!
   end function getSineMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSMatrix(spin, site, atom) result(S_mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: S_mat(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getS2C2MtxInv','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   S_mat => Scatter(id)%Solutions(ic,is)%S_mat
!
   end function getSMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGreenFunction(spin, site, atom) result(green)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: green(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getGreenFunction','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   green => Scatter(id)%Solutions(ic,is)%green
!
   end function getGreenFunction
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRegSolutionDerivative(spin, site, atom) result(drs)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: drs(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getRegSolutionDerivative','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   drs => Scatter(id)%Solutions(ic,is)%reg_dsol
!
   end function getRegSolutionDerivative
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGreenFunctionDerivative(spin, site, atom) result(dg)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: dg(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getGreenFunctionDerivative','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   dg => Scatter(id)%Solutions(ic,is)%der_green
!
   end function getGreenFunctionDerivative
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPhaseShift(spin, site, atom) result(ps)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   real (kind=RealKind), pointer :: ps(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getPhaseShift','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   ps => Scatter(id)%Solutions(ic,is)%phase_shift
!
   end function getPhaseShift
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDOS(spin, site, atom) result(dos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: dos(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getDOS','module not initialized')
   else if (.not.allocated(wks_dos)) then
!     ----------------------------------------------------------------
      call ErrorHandler('getDOS','Need to call computeDOS first')
!     ----------------------------------------------------------------
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!  
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   dos => Scatter(id)%Solutions(ic,is)%dos
!
   end function getDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDOSDerivative(spin, site, atom) result(der_dos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: der_dos(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getDOS','module not initialized')
   else if (.not.allocated(wks_ddos)) then
!     ----------------------------------------------------------------
      call ErrorHandler('getDOSDerivative','Need to call computeDOS first')
!     ----------------------------------------------------------------
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!  
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   der_dos => Scatter(id)%Solutions(ic,is)%der_dos
!
   end function getDOSDerivative
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPDOS(spin, site, atom) result(pdos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: pdos(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getPDOS','module not initialized')
   else if (.not.allocated(wks_pdos)) then
!     ----------------------------------------------------------------
      call ErrorHandler('getPDOS','Need to call computePDOS first')
!     ----------------------------------------------------------------
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   pdos => Scatter(id)%Solutions(ic,is)%pdos
!
   end function getPDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCellDOS(spin, site, atom) result(dos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   real (kind=RealKind) :: dos
!
   if (.not.Initialized) then
      call ErrorHandler('getCellDOS','module not initialized')
   else if (.not.allocated(space_integrated_dos_cell)) then
!     ----------------------------------------------------------------
      call ErrorHandler('getCellDOS','Need to call computePDOS first')
!     ----------------------------------------------------------------
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else 
      id = LocalIndex
   endif
!  
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   dos = space_integrated_dos_cell(ic,is,id)
!
   end function getCellDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMTSphereDOS(spin, site, atom) result(dos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   real (kind=RealKind) :: dos
!
   complex (kind=CmplxKind), pointer :: pdos(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getMTSphereDOS','module not initialized')
   else if (.not.allocated(space_integrated_dos_mt)) then
!     ----------------------------------------------------------------
      call ErrorHandler('getMTSphereDOS','Need to call computePDOS first')
!     ----------------------------------------------------------------
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   dos = space_integrated_dos_mt(ic,is,id)
!
   end function getMTSphereDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCellPDOS(spin, site, atom) result(pdos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   real (kind=RealKind), pointer :: pdos(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getCellPDOS','module not initialized')
   else if (.not.allocated(space_integrated_pdos_cell)) then
!     ----------------------------------------------------------------
      call ErrorHandler('getCellPDOS','Need to call computePDOS first')
!     ----------------------------------------------------------------
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   pdos => space_integrated_pdos_cell(:,ic,is,id)
!
   end function getCellPDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMTSpherePDOS(spin, site, atom) result(pdos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   real (kind=RealKind), pointer :: pdos(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getMTSpherePDOS','module not initialized')
   else if (.not.allocated(space_integrated_pdos_mt)) then
!     ----------------------------------------------------------------
      call ErrorHandler('getMTSpherePDOS','Need to call computePDOS first')
!     ----------------------------------------------------------------
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   pdos => space_integrated_pdos_mt(:,ic,is,id)
!
   end function getMTSpherePDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOutsideDOS(spin, site, atom, rs) result(dos_out)
!  ===================================================================
   use PolyhedraModule, only : getVolume, getInscrSphRadius, getOutscrSphRadius
   use BesselModule, only : IntegrateSphHankelSq
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   real (kind=RealKind), intent(in), optional :: rs
!
   logical :: outside_rs = .false.
!
   integer (kind=IntKind) :: is, id, ic, l, kl, klp
!
   real (kind=RealKind) :: dos_out, rmt, rc
   complex (kind=CmplxKind) :: cdos_out
   complex (kind=CmplxKind) :: fint(0:lmax_phi), gint(kmax_phi,kmax_phi)
!
   if (.not.Initialized) then
      call ErrorHandler('getDOSoutside','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   if (present(rs)) then
      if (rs < getInscrSphRadius(id)) then
         call ErrorHandler('getDOSoutside','rs < Rmt',rs,getInscrSphRadius(id))
      else if (rs < getOutscrSphRadius(id) .and. .not.isPotZeroOutsideRmt) then
         call ErrorHandler('getDOSoutside','rs < Rc',rs,getOutscrSphRadius(id))
      else
         outside_rs = .true.
         rc = rs
      endif
!  else if (isPotZeroOutsideRmt .or. isSphericalSolverOn) then ! "Outside" means out of muffin-tin sphere
!     outside_rs = .true.
!     rc = getInscrSphRadius(id)
   endif
!
   cdos_out = CZERO
   if (outside_rs) then
!     ----------------------------------------------------------------
!     call calIntSphHankelSq0(lmax_phi,rc,energy,fint)
      call IntegrateSphHankelSq(lmax_phi,rc,energy,fint)
!     ----------------------------------------------------------------
      do kl = kmax_phi, 1, -1
         l = lofk(kl)
!!       cdos_out = cdos_out + Scatter(id)%Solutions(ic,is)%t_mat(kl,kl)*fint(l)
         cdos_out = cdos_out + (Scatter(id)%Solutions(ic,is)%S_mat(kl,kl)-CONE)*fint(l)
      enddo
   else ! calculating the DOS outside the atomic cell
      rmt = getInscrSphRadius(id)
!     ----------------------------------------------------------------
!     call calIntSphHankelSq0(lmax_phi,rmt,energy,fint)
      call IntegrateSphHankelSq(lmax_phi,rmt,energy,fint)
!     ----------------------------------------------------------------
      do kl = kmax_phi, 1, -1
         l = lofk(kl)
!!       cdos_out = cdos_out + Scatter(id)%Solutions(ic,is)%t_mat(kl,kl)*fint(l)
         cdos_out = cdos_out + (Scatter(id)%Solutions(ic,is)%S_mat(kl,kl)-CONE)*fint(l)
      enddo
!
      rc = getOutscrSphRadius(id)
!     ----------------------------------------------------------------
      call calIntSphHankelSq1(id,lmax_phi,rmt,rc,energy,gint)
!     ----------------------------------------------------------------
      do kl = 1, kmax_phi
         do klp = 1, kmax_phi
!!          cdos_out = cdos_out - Scatter(id)%Solutions(ic,is)%t_mat(klp,kl)*gint(klp,kl)
            if (klp == kl) then
               cdos_out = cdos_out - (Scatter(id)%Solutions(ic,is)%S_mat(kl,kl)-CONE)*gint(kl,kl)
            else
               cdos_out = cdos_out - Scatter(id)%Solutions(ic,is)%S_mat(klp,kl)*gint(klp,kl)
            endif
         enddo
      enddo
   endif
!  ===================================================================
!  Note: dos_out is the total density od states outside the cell or rs
!        per spin.....
!  ===================================================================
!! dos_out = aimag(energy*cdos_out)/PI
   dos_out = real(kappa*cdos_out,kind=RealKind)/PI2
!
   end function getOutsideDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCosineMatrix(spin, site, atom) result(cos_mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: cos_mat(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getCosineMatrix','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   cos_mat => Scatter(id)%Solutions(ic,is)%cos_mat
!
   end function getCosineMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getJostMatrix(spin, site, atom) result(jost_mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: jost_mat(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getJostMatrix','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   jost_mat => Scatter(id)%Solutions(ic,is)%jost_mat
!
   end function getJostMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getJostInvMatrix(spin, site, atom) result(jinv_mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: jinv_mat(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getJostInvMatrix','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   jinv_mat => Scatter(id)%Solutions(ic,is)%jinv_mat
!
   end function getJostInvMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOmegaHatMatrix(spin, site, atom) result(OmegaHat_mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getJostMatrix','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   OmegaHat_mat => Scatter(id)%Solutions(ic,is)%OmegaHat_mat
!
   end function getOmegaHatMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOmegaHatInvMatrix(spin, site, atom) result(OmegaHatInv_mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: OmegaHatInv_mat(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getJostMatrix','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   OmegaHatInv_mat => Scatter(id)%Solutions(ic,is)%OmegaHatInv_mat
!
   end function getOmegaHatInvMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!tmat_global   function getTMatrix(spin, site, atom, global_frame) result(t_mat)
   function getTMatrix(spin, site, atom, dsize) result(t_mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind) :: is, id, ic
!tmat_global   logical, intent(in), optional :: global_frame
!tmat_global   logical :: gframe
!
   complex (kind=CmplxKind), pointer :: t_mat(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getTMatrix','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
!tmat_global   if (present(global_frame)) then
!tmat_global      gframe = global_frame
!tmat_global   else
!tmat_global      gframe = .false.
!tmat_global   endif
!
!tmat_global   if (gframe) then
!tmat_global      t_mat => Scatter(id)%tmat_global(:,:,ic)
!tmat_global   else
                  t_mat => Scatter(id)%Solutions(ic,is)%t_mat
!tmat_global   endif
   if (present(dsize)) then
      dsize = Scatter(id)%kmax_int
   endif
!
   end function getTMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getScatteringMatrix(sm_type, spin, site, atom, dsize) result(sm)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: sm_type
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: sm(:,:)
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('getScatteringMatrix','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   if (nocaseCompare(sm_type,'T-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%t_mat
   else if (nocaseCompare(sm_type,'TInv-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%t_mat_inv
   else if (nocaseCompare(sm_type,'S-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%S_mat
   else if (nocaseCompare(sm_type,'Sine-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%sin_mat
   else if (nocaseCompare(sm_type,'Cosine-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%cos_mat
   else if (nocaseCompare(sm_type,'Jost-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%jost_mat
   else if (nocaseCompare(sm_type,'JostInv-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%jinv_mat
   else if (nocaseCompare(sm_type,'Omega-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%Omega_mat
   else if (nocaseCompare(sm_type,'OmegaHat-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%OmegaHat_mat
   else if (nocaseCompare(sm_type,'OmegaHatInv-Matrix')) then
      sm => Scatter(id)%Solutions(ic,is)%OmegaHatInv_mat
   else
      call ErrorHandler('getScatteringMatrix','Scattering matrix type is invalid', &
                        sm_type(1:15))
   endif
!
   if (present(dsize)) then
      dsize = Scatter(id)%kmax_int
   endif
!
   end function getScatteringMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRegSolution(spin, site, atom) result(reg_sol)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: is, id, ic
!
   complex (kind=CmplxKind), pointer :: reg_sol(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getRegSolution','module not initialized')
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (present(atom)) then
      ic = atom
   else
      ic = 1
   endif
!
   reg_sol => Scatter(id)%Solutions(ic,is)%reg_sol
!
   end function getRegSolution
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSolutionFlags(site) result(sol_flags)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: site
   integer (kind=IntKind) :: id
!
   integer (kind=IntKind), pointer :: sol_flags(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getRegSolution','module not initialized')
   else if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   sol_flags => Scatter(id)%sol_flags
!
   end function getSolutionFlags
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSolutionRmeshSize(site) result(nr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: site
   integer (kind=IntKind) :: id
   integer (kind=IntKind) :: nr
!
   if (.not.Initialized) then
      call ErrorHandler('getSolutionRmeshSize','module not initialized')
   else if (present(site)) then
      if (site < 1 .or. site > LocalNumSites) then
         call ErrorHandler('getSolutionRmeshSize','invalid number of local atoms', &
                           LocalNumSites)
      endif
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   nr = Scatter(id)%numrs
!
   end function getSolutionRmeshSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printWronskian(kl1,kl2,wfr1,dwfr1,wfr2,dwfr2,n0,n1,ns)
!  ===================================================================
   use RadialGridModule, only : getGrid
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kl1, kl2, n0, n1, ns
   integer (kind=IntKind) :: klp, ir, kl1s, klps
!
   complex (kind=CmplxKind) :: norm, cfac
   complex (kind=CmplxKind), intent(in) :: wfr1(:,:,:)
   complex (kind=CmplxKind), intent(in) :: dwfr1(:,:,:)
   complex (kind=CmplxKind), intent(in) :: wfr2(:,:,:)
   complex (kind=CmplxKind), intent(in) :: dwfr2(:,:,:)
   complex (kind=CmplxKind) :: wronskian
!
   kl1s = kl1 - 2*mofk(kl1)
!
!  ===================================================================
!  check and print the Wronskian relation between two solutions.
!  ===================================================================
   norm=CZERO
   do klp=1,kmax_phi
      klps = klp - 2*mofk(klp)
      cfac = m1m(mofk(klp)+mofk(kl1))
      norm = norm + cfac*( wfr1(iend,klps,kl1s)*dwfr2(iend,klp,kl2) -  &
                           dwfr1(iend,klps,kl1s)*wfr2(iend,klp,kl2) )
   enddo
   write(6,'(/,a,1i5,5x,a,i5,10x,a,2d16.8)')                        &
            'Wronskian: kl1 = ',kl1,'kl2 = ',kl2,'Norm = ',norm
   if (abs(norm) < TEN2m8) then
      return
   endif
!
   do ir = n0, n1, ns
      wronskian=CZERO
      do klp=1,kmax_phi
         klps = klp - 2*mofk(klp)
         cfac = m1m(mofk(klp)+mofk(kl1))
         wronskian = wronskian + cfac*( wfr1(ir,klps,kl1s)*dwfr2(ir,klp,kl2) - &
                                        dwfr1(ir,klps,kl1s)*wfr2(ir,klp,kl2) )
      enddo
      if(abs(norm) > TEN2m8) then
!        write(6,'(1i5,8d15.7)')kl1,wfr1(ir,kl1,kl1),dwfr1(ir,kl1,kl1),wfr2(ir,kl2,kl2),dwfr2(ir,kl2,kl2)
         write(6,'(1i5,4d16.8)')ir,wronskian,CONE-wronskian/norm
      else
         write(6,'(1i5,2d16.8)')ir,wronskian
      endif
   enddo
!
   call FlushFile(6)
!
   end subroutine printWronskian
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkWaveFuncAtR(kl,ir,x,y,z)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kl, ir
!
   integer (kind=IntKind) :: l, lp, klp, ic
!
   real (kind=RealKind), intent(in) :: x,y,z
!
   complex (kind=CmplxKind) :: ylmr(kmax_phi), phi_L
!
!  -------------------------------------------------------------------
   call calYlm(x,y,z,lmax_phi,ylmr)
!  -------------------------------------------------------------------
!
   do ic =1, Scatter(LocalIndex)%NumSpecies
      phi_L = CZERO
      do klp = 1, kmax_phi
         phi_L = phi_L + ylmr(klp)*Scatter(LocalIndex)%Solutions(ic,LocalSpin)%reg_sol(ir,klp,kl)
      enddo
      write(6,'(a,i2,2x,f12.7,2x,2d15.7)')'Phi_L at r = ',kl,Grid%r_mesh(ir),phi_L
   enddo
!
   end subroutine checkWaveFuncAtR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeIrrSol(bc,atom,fmem,pm)
!  ===================================================================
   use LdaCorrectionModule, only : checkLdaCorrection, transformWaveFunction
!
!  *******************************************************************
!  calculate the irregular solution with the following boundary condition
!     sum_L' pm(L',L)*j_l'(kappa*r)*Y_L' r >= Rc
!  where max(L') = kmax_phi, max(L) = kmax_kkr. !
!  pm is a matrix to be multiplied to the sperical Bessel (or Hankel) function
!  at the boundary. If pm is not present, no such transformation is needed.
!  *******************************************************************
   implicit none
!
   character (len=1), intent(in) :: bc
!
   integer (kind=IntKind), intent(in) :: atom
   integer (kind=IntKind) :: l, l_save, kl, klp, ir, lp, i
   integer (kind=IntKind) :: m, mp, kls, klps, klpp
!
   complex (kind=CmplxKind), optional, intent(in) :: pm(:,:)
   complex (kind=CmplxKind), pointer :: bsf(:,:)
   complex (kind=CmplxKind), pointer :: dbsf(:,:)
   complex (kind=CmplxKind), intent(inout), target :: fmem(iend,4)
!
   complex (kind=CmplxKind), pointer :: invscm(:,:), invscmt(:,:)
   complex (kind=CmplxKind), pointer :: smpicmt(:,:)
   complex (kind=CmplxKind), pointer :: tmp_irr(:,:,:)
   complex (kind=CmplxKind), pointer :: wfr_irr(:,:,:)
   complex (kind=CmplxKind), pointer :: wfdr_irr(:,:,:)
   complex (kind=CmplxKind), pointer :: pw(:), wfr(:,:), dwfr(:,:)
   complex (kind=CmplxKind) :: eoc2p1
!
   if (bc == 'H' .or. bc == 'h') then
      bsf => bhl
      dbsf => dbhl
   else if (bc == 'J' .or. bc == 'j') then
      bsf => bjl
      dbsf => dbjl
   else if (bc == 'N' .or. bc == 'n') then
      bsf => bnl
      dbsf => dbnl
   else
!     ----------------------------------------------------------------
      call ErrorHandler('computeIrrSol','Invalid boundary condition type',bc)
!     ----------------------------------------------------------------
   endif
!
   wfr_irr => Scatter(LocalIndex)%Solutions(atom,LocalSpin)%irr_sol(1:iend,1:kmax_phi,1:kmax_kkr)
   wfdr_irr => Scatter(LocalIndex)%Solutions(atom,LocalSpin)%irr_dsol(1:iend,1:kmax_phi,1:kmax_kkr)
   wfr_irr = CZERO
   wfdr_irr = CZERO
   eoc2p1=CONE+energy*c2inv
!
!  if ( SphericalSolver .or. isSphPotZeroOutsideRmt ) then
   if (isSphericalSolverOn) then
      m = (iend-numrs_mt+1)*kmax_kkr
      pw => wks_plhat(1:m)
      wfr => aliasArray2_c(pw,iend-numrs_mt+1,kmax_kkr)
      pw => wks_plhat(m+1:2*m)
      dwfr => aliasArray2_c(pw,iend-numrs_mt+1,kmax_kkr)
      if ( present(pm) ) then
         do kl = 1, kmax_kkr
            do i = 1, iend-numrs_mt+1
               wfr(i,kl) = CZERO; dwfr(i,kl) = CZERO
               ir = iend+1-i
               do klp = 1, kmax_phi
                  lp = lofk(klp)
                  wfr(i,kl) = wfr(i,kl) + pm(klp,kl)*bsf(ir,lp)/kappa
                  dwfr(i,kl) = dwfr(i,kl) + pm(klp,kl)*dbsf(ir,lp)/eoc2p1
               enddo
            enddo
         enddo
      else
         do kl = 1, kmax_kkr
            l = lofk(kl)
            do i = 1, iend-numrs_mt+1
               ir = iend+1-i
               wfr(i,kl) = bsf(ir,l)/kappa
               dwfr(i,kl) = dbsf(ir,l)/eoc2p1
            enddo
         enddo
      endif
!
      l_save = -1
      do kl = 1,kmax_kkr
         l=lofk(kl)
         if (l > l_save) then
            pl_irr => pl_irrwrk(1:iend,l)
            ql_irr => ql_irrwrk(1:iend,l)
            pl_irrwrk(1:iend,l) = CZERO
            ql_irrwrk(1:iend,l) = CZERO
!           ----------------------------------------------------------
            call  solveRadEqn4Irr(kl,iend-numrs_mt+1,wfr(:,kl),dwfr(:,kl),fmem)
!           ----------------------------------------------------------
            l_save = l
         endif
!        -------------------------------------------------------------
         call zcopy(iend,pl_irrwrk(1:iend,l),1,wfr_irr(1:iend,kl,kl),1)
         call zcopy(iend,ql_irrwrk(1:iend,l),1,wfdr_irr(1:iend,kl,kl),1)
!        -------------------------------------------------------------
      enddo
   else
      m = (iend-numrs_cs)*kmax_kkr
      pw => wks_plhat(1:m)
      wfr => aliasArray2_c(pw,iend-numrs_cs,kmax_kkr)
      pw => wks_plhat(m+1:2*m)
      dwfr => aliasArray2_c(pw,iend-numrs_cs,kmax_kkr)
      if ( present(pm) ) then
         do kl = 1, kmax_kkr
            do i = 1, iend-numrs_cs
               wfr(i,kl) = CZERO; dwfr(i,kl) = CZERO
               ir = iend+1-i
               do klp = 1, kmax_phi
                  lp = lofk(klp)
                  wfr(i,kl) = wfr(i,kl) + pm(klp,kl)*bsf(ir,lp)/kappa
                  dwfr(i,kl) = dwfr(i,kl) + pm(klp,kl)*dbsf(ir,lp)/eoc2p1
               enddo
            enddo
         enddo
      else
         do kl = 1, kmax_kkr
            l = lofk(kl)
            do i = 1, iend-numrs_cs
               ir = iend+1-i
               wfr(i,kl) = bsf(ir,l)/kappa
               dwfr(i,kl) = dbsf(ir,l)/eoc2p1
            enddo
         enddo
      endif
!
      do kl = 1,kmax_kkr
         plhat_irr=>wfr_irr(1:iend,1:kmax_phi,kl)
         qlhat_irr=>wfdr_irr(1:iend,1:kmax_phi,kl)
!        -------------------------------------------------------------
!        call solveIntEqn(kl,fmem,.false.,.true.)
!        plhat_irr(1:iend,kl)=plhat_irr(1:iend,kl)+pl_irr(1:iend)
!        -------------------------------------------------------------
         call calPhiLr(kl,fmem,iend-numrs_cs,wfr(:,kl),dwfr(:,kl),.true.)
!        -------------------------------------------------------------
      enddo
   endif
!
   if (checkLdaCorrection(LocalIndex,atom)) then
      call transformWaveFunction(LocalIndex,atom,spin_index,iend,kmax_phi,wfr_irr)
   endif
!
   nullify( wfr_irr, wfdr_irr, pw, dwfr, wfr )
!
   end subroutine computeIrrSol
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine solveRadEqn4Irr(kl,np,wfr,dwfr,fmem)
!  ===================================================================
   use InterpolationModule, only : PolyInterp, FitInterp
!
!  *******************************************************************
!  calculate irregular solution for the radial schrodinger (or scalar-
!  relativistic) eqn for the spherical component of the potential
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in) :: kl, np
   integer (kind=IntKind) :: klp, ir, l, lp, ibeg, j, j_inter, llp1
   integer (kind=IntKind), parameter :: n_inter=5
!
   real (kind=RealKind), pointer :: xp(:)
   real (kind=RealKind) :: rp
   real (kind=RealKind) :: xmt
   real (kind=RealKind) :: z_err_est
!
   complex (kind=CmplxKind), intent(in) :: wfr(np), dwfr(np)
   complex (kind=CmplxKind), intent(inout), target :: fmem(iend,4)
   complex (kind=CmplxKind), pointer :: yp(:)
   complex (kind=CmplxKind), pointer :: cl(:)
   complex (kind=CmplxKind), pointer :: sl(:)
   complex (kind=CmplxKind), pointer :: dcl(:)
   complex (kind=CmplxKind), pointer :: dsl(:)
   complex (kind=CmplxKind) :: clp_irr(0:3)
   complex (kind=CmplxKind) :: slp_irr(0:3)
   complex (kind=CmplxKind) :: dclp_irr(0:3)
   complex (kind=CmplxKind) :: dslp_irr(0:3)
   complex (kind=CmplxKind) :: eoc2p1
   complex (kind=CmplxKind) :: e2oc2
   complex (kind=CmplxKind) :: ddcl, ddsl, dgpl, dfpl, b1, b2, em, emr
   complex (kind=CMplxKind) :: rvr, drvr
!
   e2oc2=energy*energy*c2inv
   eoc2p1=CONE+energy*c2inv
   l = lofk(kl)
   llp1 = l*(l+1)
!
!  if ( SphericalSolver .or. isSphPotZeroOutsideRmt ) then
   if (isSphericalSolverOn) then
      ibeg = numrs_mt
   else
      ibeg = numrs_cs
   endif
   if ( ibeg+2 >  iend ) then
!     ----------------------------------------------------------------
      call ErrorHandler('solveRadEq4Irr','ibeg+2 > iend',ibeg+2,iend)
!     ----------------------------------------------------------------
   endif
!
   cl  => fmem(1:iend,1)
   sl  => fmem(1:iend,2)
   dcl => fmem(1:iend,3)
   dsl => fmem(1:iend,4)
   xmt=Grid%x_mesh(Grid%jmt)
!
!  ===================================================================
!  starting with free solution ( v(r)=0 r > R) to get Pl_irr and 
!  Ql_irr, but using extrapolated potentials at and beyond R to 
!  get derivative terms of cl and sl.
!  Note: effect of potential on ql_irr here ignored: it is small 
!  ===================================================================
   do j = 1, np
      ir = ibeg + np - j
      pl_irr(ir) = wfr(j)
      ql_irr(ir) = dwfr(j)
   enddo
!
   do j=ibeg,ibeg+3
      if (.true.) then
!        =============================================================
!        the following code is the original lines, which were commented
!        out on 11/19/03, and revised by YW.
!        The idea is to extrapolate the muffin-tin potential so that the
!        quantity dcl and dsl has a smooth extension beyond the muffin-tin
!        so to allow for better numerics when using the predictor-corrector algorithm
!        for the grid points near the muffin-tin radius. This helps to
!        improve the Wronskian relations near the muffin-tin radius
!        =============================================================
!ywg     cl(j)=-bnl(j,1)*pl_irr(j) - bnl(j,0)*ql_irr(j)/kappa
!ywg     sl(j)=-bjl(j,1)*pl_irr(j) - bjl(j,0)*ql_irr(j)/kappa
!ywg     dcl(j)=-bnl(j,1)*(cm0(j)-CONE)*ql_irr(j)*Grid%r_mesh(j)         &
!ywg            -( llp1/(cm0(j)*Grid%r_mesh(j))+v0r(j)+                  &
!ywg               e2oc2*Grid%r_mesh(j) )*pl_irr(j)*bnl(j,0)/kappa
!ywg     dsl(j)=-bjl(j,1)*(cm0(j)-CONE)*ql_irr(j)*Grid%r_mesh(j)         &
!ywg            -( llp1/(cm0(j)*Grid%r_mesh(j))+v0r(j)+                  &
!ywg               e2oc2*Grid%r_mesh(j) )*pl_irr(j)*bjl(j,0)/kappa
         cl(j) = ql_irr(j)*bnl(j,l)-pl_irr(j)*dbnl(j,l)*kappa
         sl(j) = ql_irr(j)*bjl(j,l)-pl_irr(j)*dbjl(j,l)*kappa
         if (j == ibeg) then
            dcl(j)=bnl(j,l)*( llp1*(CONE-cm0(j))/Grid%r_mesh(j)+v0r(j)+e2oc2*Grid%r_mesh(j) )*pl_irr(j)
            dsl(j)=bjl(j,l)*( llp1*(CONE-cm0(j))/Grid%r_mesh(j)+v0r(j)+e2oc2*Grid%r_mesh(j) )*pl_irr(j)
         else
!           ----------------------------------------------------------
            call FitInterp(ibeg,Grid%r_mesh,v0r,Grid%r_mesh(j),rvr,drvr)
!           ----------------------------------------------------------
            dcl(j)= bnl(j,l)*( llp1*(CONE-cm0(j))/Grid%r_mesh(j)+rvr+e2oc2*Grid%r_mesh(j) )*pl_irr(j)
            dsl(j)= bjl(j,l)*( llp1*(CONE-cm0(j))/Grid%r_mesh(j)+rvr+e2oc2*Grid%r_mesh(j) )*pl_irr(j)
         endif
      else
!        =============================================================
!        in the mean time the following code is copied-pasted over from
!        version 1.7 (on 11/19/03).
!        =============================================================
         if(j == ibeg) then
!ywg        cl(j)=-bnl(j,1)*pl_irr(j) - bnl(j,0)*ql_irr(j)/kappa
!ywg        sl(j)=-bjl(j,1)*pl_irr(j) - bjl(j,0)*ql_irr(j)/kappa
            cl(j) = ql_irr(j)*bnl(j,l)-pl_irr(j)*dbnl(j,l)*kappa
            sl(j) = ql_irr(j)*bjl(j,l)-pl_irr(j)*dbjl(j,l)*kappa
            ddcl=CZERO
            ddsl=CZERO
            dgpl=CZERO
            dfpl=CZERO
         endif
!ywg     dgpl=dgpl+(ddcl*bjl(j,0)-ddsl*bnl(j,0))*Grid%hout
!ywg     dfpl=dfpl-kappa*(ddcl*bjl(j,1)-ddsl*bnl(j,1))*Grid%hout
         dgpl=dgpl+(ddcl*bjl(j,0)-ddsl*bnl(j,0))*Grid%hout
         dfpl=dfpl-kappa*(ddcl*bjl(j,1)-ddsl*bnl(j,1))*Grid%hout
!        -------------------------------------------------------------
         call FitInterp(ibeg,Grid%r_mesh,v0r,Grid%r_mesh(j),rvr,drvr)
!        -------------------------------------------------------------
         em = eoc2p1 - c2inv*rvr/Grid%r_mesh(j)
         emr= em*Grid%r_mesh(j)
!        =============================================================
!        Note: effect of potential on ql_irr here ignored: it is small.
!        =============================================================
         b1 = (em-CONE)*(ql_irr(j)+dfpl)*Grid%r_mesh(j)
         b2 = ( l*(l+1)/emr + rvr +                                   &
                (eoc2p1-CONE)*energy*Grid%r_mesh(j) )*(pl_irr(j)+dgpl)/kappa
         dcl(j) = -bnl(j,1)*b1 - bnl(j,0)*b2
         dsl(j) = -bjl(j,1)*b1 - bjl(j,0)*b2
         b1 = (em-CONE)*dfpl*Grid%r_mesh(j)
         b2 = b2 - ( l*(l+1)/emr +                                    &
                     (eoc2p1-CONE)*energy*Grid%r_mesh(j) )*pl_irr(j)/kappa
         ddcl = -bnl(j,1)*b1-bnl(j,0)*b2
         ddsl = -bjl(j,1)*b1-bjl(j,0)*b2
      endif
   enddo
!
!  if ( .not. FullSolver .or. abs(Grid%hin-Grid%hout) < TEN2m8 .or.    & 
!       isSphPotZeroOutsideRmt ) then
   if (abs(Grid%hin-Grid%hout) < TEN2m8 .or. isSphericalSolverOn) then
!     ----------------------------------------------------------------
      call precor4(l,ibeg-1,1,Grid%hin,IrregularSolution,             &
                    e2oc2,Grid%r_mesh,cl,dcl,sl,dsl)
!     ----------------------------------------------------------------
   else                                            ! double logrithmic grid
      write(6,'(/,a,/)')'It needs to be checked! *********************'
!     ================================================================
!     get irregular solution inside rcirc to rmt......................
!     ================================================================
      if(ibeg-1 > Grid%jmt) then
!        -------------------------------------------------------------
         call precor4(l,ibeg-1,Grid%jmt,Grid%hout,IrregularSolution,  &
                       e2oc2,Grid%r_mesh,cl,dcl,sl,dsl)
!        -------------------------------------------------------------
      endif
!
!     ================================================================
!     generate pseudo-grid points outside mt sphere and pass equation
!     solution values.................................................
!     ================================================================
      j_inter=Grid%jmt
      clp_irr(0)=cl(j_inter)
      slp_irr(0)=sl(j_inter)
      dclp_irr(0)=dcl(j_inter)
      dslp_irr(0)=dsl(j_inter)
      do j=1,3
         rp=exp(xmt+Grid%hin*j)
!        -------------------------------------------------------------
         call hunt(iend,Grid%r_mesh,rp,j_inter)
!        -------------------------------------------------------------
         j_inter=min(max(j_inter-(n_inter-1)/2,Grid%jmt),iend+1-n_inter)
         xp=>Grid%r_mesh(j_inter:j_inter+n_inter-1)
         yp=>cl(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,rp,clp_irr(j),z_err_est)
!        -------------------------------------------------------------
         yp=>sl(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,rp,slp_irr(j),z_err_est)
!        -------------------------------------------------------------
         yp=>dcl(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,rp,dclp_irr(j),z_err_est)
!        -------------------------------------------------------------
         yp=>dsl(j_inter:j_inter+n_inter-1)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,xp,yp,rp,dslp_irr(j),z_err_est)
!        -------------------------------------------------------------
      enddo
!
!     ================================================================
!     get irregular solution inside Muffin-tin sphere.................
!     ================================================================
      do j=0,3
         cl(Grid%jmt+j)=clp_irr(j)
         sl(Grid%jmt+j)=slp_irr(j)
         dcl(Grid%jmt+j)=dclp_irr(j)
         dsl(Grid%jmt+j)=dslp_irr(j)
      enddo
!     ----------------------------------------------------------------
      call precor4(l,Grid%jmt-1,1,Grid%hin,IrregularSolution,        &
                    e2oc2,Grid%r_mesh,cl,dcl,sl,dsl)
!     ----------------------------------------------------------------
   endif
!
   nullify( xp, yp, cl, sl, dsl, dcl )
!
   end subroutine solveRadEqn4Irr
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!tmat_global   subroutine convertTmatToGlobalFrame(site)
!  ===================================================================
!tmat_global   use SpinRotationModule, only : rotateLtoG
!tmat_global   implicit none
!
!tmat_global   integer (kind=IntKind), intent(in), optional :: site
!tmat_global   integer (kind=IntKind) :: id, ic, dsize
!
!tmat_global   complex (kind=CmplxKind), pointer :: tm1(:,:), tm2(:,:), tm0(:,:)
!
!tmat_global   if (present(site)) then
!tmat_global      id = site
!tmat_global   else
!tmat_global      id = LocalIndex
!tmat_global   endif
!
!tmat_global   if (NumSpins == 2) then
!tmat_global      do ic = 1, Scatter(id)%NumSpecies
!tmat_global         tm1 => Scatter(id)%Solutions(ic,1)%tmat
!tmat_global         tm2 => Scatter(id)%Solutions(ic,2)%tmat
!tmat_global         tm0 => Scatter(id)%tmat_global(:,:,ic)
!tmat_global         dsize = Scatter(id)%kmax_int
!        =============================================================
!        NOTE:
!        id is the local index for the medium host, therefore it needs to
!        find appropriate index for the following call. This will be looked
!        into further in the future (12/31/2016 - YW)
!        -------------------------------------------------------------
!tmat_global         call rotateLtoG(id, dsize, dsize, tm1, tm2, tm0)
!        -------------------------------------------------------------
!tmat_global      enddo
!tmat_global   endif
!
!tmat_global   end subroutine convertTmatToGlobalFrame
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computePhaseShift(spin,site,atom)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: id, is, ic, kmax, INFO, kl
!
   real (kind=RealKind) :: si, sr
!
   complex (kind=CmplxKind) :: smdexp
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   kmax = Scatter(id)%kmax_phi
!
   LOOP_ic: do ic = 1, Scatter(id)%NumSpecies
      if (present(atom)) then
         if (ic /= atom .and. atom > 0) then
            cycle LOOP_ic
         endif
      endif
!     ================================================================
!     Triangularize the S_matrix
!     ----------------------------------------------------------------
      call ZGETRF_nopivot(kmax, kmax, Scatter(id)%Solutions(ic,is)%S_mat, &
                          kmax, IPVT, INFO)
!     ----------------------------------------------------------------
      if (INFO /= 0) then
         call ErrorHandler('computePhaseShift','Unsuccessful S-matrix diagonalization',INFO)
      endif
!
!     ================================================================
!     Determine the generialized partial phase shift based on
!     triangularized S_matrix (Due to Balazs Gyorffy)
!           The partial phase shifts are stored in diags array
!           Notes: adding or subtracting a PI to the partial phase
!                  shift is necessary when a diagnonal element of
!                  the S-matrix passes the real energy axis from the
!                  2nd quadrant to the 3rd quadrant, or from from the
!                  3rd quadrant to the 2nd quadrant.
!           For a given (l,m), S_matrix(kl,kl) = (sr,si),
!                              2*phase*i = log(S_matrix(kl,kl)).
!           If phase stays positive, which is likely the case for l > 0(?),
!           as energy increases, the phase may jumps from postive to
!           negative if S_matrix crosses from the 2nd quadrant into
!           the 3rd quadrant. In this happes, a PI needs to be added to phase.
!           In the case of l = 0, as energy increases, the phase changes clockwisely
!           and may enter from the 3rd quadrant into the 2nd quadrant. If this happens
!           the phase needs to be subtracted by a PI.
!     The following scheme will work if the phase is calculated along
!     the real energy axis from left to right.
!     ================================================================
      do kl = 1, kmax
         smdexp = log(Scatter(id)%Solutions(ic,is)%S_mat(kl,kl))
         si = aimag(Scatter(id)%Solutions(ic,is)%S_mat(kl,kl))
         sr = real(Scatter(id)%Solutions(ic,is)%S_mat(kl,kl), kind=RealKind)
!        =============================================================
!        The following scheme for adding or subtracting PI is commented out
!        It will be done outside the module..............................
         if (npi(kl,ic,spin_index,id) == 0 .and. Scatter(id)%AtomicNumber(ic) > 0) then
            if (lofk(kl) == 0 .and. (sr < ZERO .and. si > ZERO)) then
               npi(kl,ic,spin_index,id) = -1
            else if (lofk(kl) == 2 .and. si < ZERO) then ! This happens to Zn: the phase of
               npi(kl,ic,spin_index,id) = 1              ! S matrix jumps from 2nd quadrant
            endif                                        ! to 4th quadrant
         endif
         Scatter(id)%Solutions(ic,is)%phase_shift(kl) = aimag(smdexp)/TWO &
                        + npi(kl,ic,spin_index,id)*PI - phase_ref(kl,ic,spin_index,id)
!        =============================================================
      enddo
   enddo LOOP_ic
!
   end subroutine computePhaseShift
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeGF0(spin,site,atom)
!  ===================================================================
   use SystemSymmetryModule, only : getSymmetryFlags
!
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: id, is, ic
   integer (kind=IntKind) :: kl1, kl2, kl1p, kl2p, klg, ir, i
   integer (kind=IntKind) :: kl2pc, kl2c, ma, m2, m2p, kl, klp, m, mp
   integer (kind=IntKind) :: sz, loc1, loc2, nr, km, ia, js
   integer (kind=IntKind) :: kmax_green_loc, kmax_kkr_loc, kmax_phi_loc
   integer (kind=IntKind), pointer :: green_flags(:)
!
   complex (kind=CmplxKind) :: cmat, cfac
   complex (kind=CmplxKind), pointer :: green(:,:)
   complex (kind=CmplxKind), pointer :: der_green(:,:)
   complex (kind=CmplxKind), pointer :: p_tmp(:,:), q_tmp(:,:)
   complex (kind=CmplxKind), pointer :: wfr_reg(:,:,:)
   complex (kind=CmplxKind), pointer :: wfr_irr(:,:,:)
   complex (kind=CmplxKind), pointer :: der_wfr_reg(:,:,:)
   complex (kind=CmplxKind), pointer :: der_wfr_irr(:,:,:)
   complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)
   complex (kind=CmplxKind), pointer :: mat(:,:), sin_t(:,:), sin_mat(:,:)
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   if (IrrSolType /= 'H' .and. IrrSolType /= 'h') then
      call ErrorHandler('computeGF0',                  &
              'Green function is not yet implemented for this bounday condition',IrrSolType)
   endif
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (.not. allocated(wks_green)) then
      sz = 0
      do ia = 1, LocalNumSites
         nr = Scatter(ia)%numrs_cs
         km = Scatter(ia)%kmax_green
         sz = sz + nr*km*NumSpins*Scatter(ia)%NumSpecies
      enddo
      allocate( wks_green(sz) )
      loc1 = 0
      do ia = 1, LocalNumSites
         nr = Scatter(ia)%numrs_cs
         km = Scatter(ia)%kmax_green
         do js = 1, NumSpins
            do ic = 1, Scatter(ia)%NumSpecies
               loc2 = loc1 + nr*km
               p1c => wks_green(loc1+1:loc2)
               Scatter(ia)%Solutions(ic,js)%green => aliasArray2_c(p1c,nr,km)
               loc1 = loc2
            enddo
         enddo
      enddo
      if (rad_deriv) then
         allocate( wks_dgreen(sz) )
         loc1 = 0
         do ia = 1, LocalNumSites
            nr = Scatter(ia)%numrs_cs
            km = Scatter(ia)%kmax_green
            do js = 1, NumSpins
               do ic = 1, Scatter(ia)%NumSpecies
                  loc2 = loc1 + nr*km
                  p1c => wks_dgreen(loc1+1:loc2)
                  Scatter(ia)%Solutions(ic,js)%der_green => aliasArray2_c(p1c,nr,km)
                  loc1 = loc2
               enddo
            enddo
         enddo
      endif
   endif
!
   kmax_green_loc = Scatter(id)%kmax_green
   kmax_kkr_loc = Scatter(id)%kmax_kkr
   kmax_phi_loc = Scatter(id)%kmax_phi
   sin_t => aliasArray2_c(wks_mtx1,kmax_phi,kmax_kkr)
   mat => aliasArray2_c(wks_mtx2,kmax_kkr,kmax_phi)
   nr = Scatter(id)%numrs_cs
!
   LOOP_ic: do ic = 1, Scatter(id)%NumSpecies
      if (present(atom)) then
         if (atom > 0 .and. ic /= atom) then
            cycle LOOP_ic
         endif
      endif
      wfr_reg => Scatter(id)%Solutions(ic,is)%reg_sol
      wfr_irr => Scatter(id)%Solutions(ic,is)%irr_sol
      green => Scatter(id)%Solutions(ic,is)%green
      green = CZERO
      if (rad_deriv) then
         der_wfr_reg => Scatter(id)%Solutions(ic,is)%reg_dsol
         der_wfr_irr => Scatter(id)%Solutions(ic,is)%irr_dsol
         der_green => Scatter(id)%Solutions(ic,is)%der_green
         der_green = CZERO
      endif
!
!     ================================================================
!     compute [i*S - C]^{-1} 
!     ================================================================
      sin_mat => Scatter(id)%Solutions(ic,is)%sin_mat
      do kl = 1,kmax_kkr
         m = mofk(kl)
         do klp = 1,kmax_phi
            mp = mofk(klp)
            sin_t(klp,kl) = m1m(m+mp)*sin_mat(klp-2*mp,kl-2*m)
         enddo
      enddo
      OmegaHat_mat => Scatter(id)%Solutions(ic,is)%OmegaHat_mat
!     ----------------------------------------------------------------
      call zgemm( 'n', 't', kmax_kkr, kmax_phi, kmax_kkr, CONE,       &
                  OmegaHat_mat, kmax_kkr, sin_t, kmax_phi, CZERO, mat, kmax_kkr)
!     ----------------------------------------------------------------
!
!     if ( SphericalSolver .or. isSphPotZeroOutsideRmt ) then
      if (isSphericalSolverOn) then
         do kl1 = 1, kmax_kkr_loc
            do i = 1, nj3(kl1,kl1)
               klg = kj3(i,kl1,kl1)
               if (klg <= kmax_green_loc) then
                  cfac = cgnt(i,kl1,kl1)*mat(kl1,kl1)
                  do ir = 1, nr
                     green(ir,klg) = green(ir,klg) + cfac*wfr_reg(ir,kl1,kl1)*wfr_irr(ir,kl1,kl1)
                  enddo
                  if (rad_deriv) then
                     do ir = 1, nr
                        der_green(ir,klg) = der_green(ir,klg)                  &
                            + cfac*der_wfr_reg(ir,kl1,kl1)*wfr_irr(ir,kl1,kl1) &
                            + cfac*wfr_reg(ir,kl1,kl1)*der_wfr_irr(ir,kl1,kl1)
                     enddo
                  endif
               endif
            enddo
         enddo
      else
         p_tmp => aliasArray2_c(wks_plhat,Scatter(id)%numrs,kmax_phi_loc)
         if (rad_deriv) then
            q_tmp => aliasArray2_c(wks_qlhat,Scatter(id)%numrs,kmax_phi_loc)
         endif
         do kl2 = 1, kmax_kkr_loc
!!          p_tmp = CZERO
!!          do kl1 = 1, kmax_kkr_loc
!!             do kl1p = 1, kmax_phi_loc
!!                do ir = 1, nr
!!                   p_tmp(ir,kl1p) = p_tmp(ir,kl1p) + mat(kl1,kl2)*wfr_reg(ir,kl1p,kl1)
!!                enddo
!!             enddo
!!          enddo
!           ----------------------------------------------------------
            call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                       CONE,wfr_reg,Scatter(id)%numrs*kmax_phi_loc,     &
                       mat(1:kmax_kkr_loc,kl2),1,CZERO,p_tmp,1)
!           ----------------------------------------------------------
            if (rad_deriv) then
!              -------------------------------------------------------
               call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                          CONE,der_wfr_reg,Scatter(id)%numrs*kmax_phi_loc, &
                          mat(1:kmax_kkr_loc,kl2),1,CZERO,q_tmp,1)
!              -------------------------------------------------------
            endif
            m2 = mofk(kl2); kl2c = kl2-2*m2
            do kl2p = 1, kmax_phi_loc
               m2p = mofk(kl2p); kl2pc = kl2p-2*m2p
               ma = abs(m2+m2p); ma = 1-2*mod(ma,2)
               do kl1p = 1, kmax_phi_loc
                  do i = 1, nj3(kl2p,kl1p)
                     klg = kj3(i,kl2p,kl1p)
                     if (klg <= kmax_green_loc) then
                        cfac = ma*cgnt(i,kl2p,kl1p)
                        do ir = 1, nr
                           green(ir,klg) = green(ir,klg) + cfac*p_tmp(ir,kl1p)*wfr_irr(ir,kl2pc,kl2c)
                        enddo
                        if (rad_deriv) then
                           do ir = 1, nr
                              der_green(ir,klg) = der_green(ir,klg) +          &
                                  cfac*q_tmp(ir,kl1p)*wfr_irr(ir,kl2pc,kl2c) + &
                                  cfac*p_tmp(ir,kl1p)*der_wfr_irr(ir,kl2pc,kl2c)
                           enddo
                        endif
                     endif
                  enddo
               enddo
            enddo
         enddo
         nullify( p_tmp, q_tmp )
      endif
      cfac = -SQRTm1*kappa
      green = cfac*green
      if (rad_deriv) then
         der_green = cfac*der_green
      endif
!
!     ================================================================
!     Symmetrizing the green function, if needed.  Added by Yang on 09-21-2018
!     ================================================================
      if (isDosSymmOn) then
         green_flags => getSymmetryFlags(id)
         do klg = 1, kmax_green_loc
            if (green_flags(jofk(klg)) == 0) then
               green(:,klg) = CZERO
               if (rad_deriv) then
                  der_green(:,klg) = CZERO
               endif
            endif
         enddo
      endif
!     ================================================================
   enddo LOOP_ic
!
   nullify( wfr_reg, wfr_irr, mat, OmegaHat_mat,                      &
            green, der_green, der_wfr_reg, der_wfr_irr )
!
   end subroutine computeGF0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeGF1(rvec,gf,spin,site,atom,dgf)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
   use InterpolationModule, only : PolyInterp
   use PublicTypeDefinitionsModule, only : GridStruct
   use RadialGridModule, only : getGrid
!
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: id, is, ic
   integer (kind=IntKind) :: kl1, kl2, klp1, klp2, m, mp, kl1s, klp1s, kl, klp
   integer (kind=IntKind) :: nr, lmax_phi_loc, kmax_kkr_loc, kmax_phi_loc
   integer (kind=IntKind), parameter :: n_inter=5
   integer (kind=IntKind) :: j_inter
!
   real (kind=RealKind), intent(in) :: rvec(3)
   real (kind=RealKind) :: r
   real (kind=RealKind) :: z_err_est
   real (kind=RealKind), pointer :: xp(:)
!
   complex (kind=CmplxKind), intent(out) :: gf(:)
   complex (kind=CmplxKind), intent(out), optional :: dgf(:)
   complex (kind=CmplxKind) :: cmat, yic
   complex (kind=CmplxKind), pointer :: yp(:)
   complex (kind=CmplxKind), pointer :: wfr_reg(:,:,:)
   complex (kind=CmplxKind), pointer :: wfr_irr(:,:,:)
   complex (kind=CmplxKind), pointer :: mat(:,:), sin_t(:,:), sin_mat(:,:)
   complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)
   complex (kind=CmplxKind), pointer :: regf(:)
   complex (kind=CmplxKind), pointer :: irrf(:)
   complex (kind=CmplxKind), pointer :: ylmv(:)
   complex (kind=CmplxKind), pointer :: ylmvc(:)
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   type (GridStruct), pointer :: Grid_loc
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
!  if (present(atom)) then
!     ic = atom
!  else
!     ic = 1
!  endif
!
   lmax_phi_loc = Scatter(id)%lmax_phi
   kmax_kkr_loc = Scatter(id)%kmax_kkr
   kmax_phi_loc = Scatter(id)%kmax_phi
   nr = Scatter(id)%numrs
   Grid_loc => getGrid(id)
   r = sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
   j_inter=Grid_loc%jmt
!  -------------------------------------------------------------------
   call hunt(nr,Grid_loc%r_mesh,r,j_inter)
!  -------------------------------------------------------------------
   if (j_inter > nr-(n_inter-1)/2) then
      j_inter=nr-n_inter+1
   else if (2*j_inter+1 > n_inter) then
      j_inter=j_inter-(n_inter-1)/2
   else
      j_inter=1
   endif
   xp=>Grid_loc%r_mesh(j_inter:j_inter+n_inter-1)
!
   ylmv  => wks_mtx1(1:kmax_phi_loc)
   ylmvc => wks_mtx2(1:kmax_phi_loc)
   regf => wks_mtx1(kmax_phi_loc+1:kmax_phi_loc+kmax_kkr_loc)
   irrf => wks_mtx2(kmax_phi_loc+1:kmax_phi_loc+kmax_kkr_loc)
   p1c => wks_mtx1(kmax_phi_loc+kmax_kkr_loc+1:)
   sin_t => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
   p1c => wks_mtx2(kmax_phi_loc+kmax_kkr_loc+1:)
   mat => aliasArray2_c(p1c,kmax_kkr,kmax_phi)
!  -------------------------------------------------------------------
   call calYlm(rvec,lmax_phi_loc,ylmv)
!  -------------------------------------------------------------------
   do klp1 = 1, kmax_phi_loc
      ylmvc(klp1) = conjg(ylmv(klp1))
   enddo
!
   LOOP_ic: do ic = 1, Scatter(id)%NumSpecies
      if (present(atom)) then
         if (atom > 0 .and. ic /= atom) then
            cycle LOOP_ic
         endif
      endif
      wfr_reg => Scatter(id)%Solutions(ic,is)%reg_sol
      wfr_irr => Scatter(id)%Solutions(ic,is)%irr_sol
!     ================================================================
!     compute [i*S - C]^{-1} 
!     ================================================================
      sin_mat => Scatter(id)%Solutions(ic,is)%sin_mat
      do kl = 1,kmax_kkr
         m = mofk(kl)
         do klp = 1,kmax_phi
            mp = mofk(klp)
            sin_t(klp,kl) = m1m(m+mp)*sin_mat(klp-2*mp,kl-2*m)
         enddo
      enddo
      OmegaHat_mat => Scatter(id)%Solutions(ic,is)%OmegaHat_mat
!     ----------------------------------------------------------------
      call zgemm( 'n', 't', kmax_kkr, kmax_phi, kmax_kkr, CONE,          &
                  OmegaHat_mat, kmax_kkr, sin_t, kmax_phi, CZERO, mat, kmax_kkr)
!     ----------------------------------------------------------------
!
      do kl1 = 1,kmax_kkr_loc
         regf(kl1) = CZERO
         irrf(kl1) = CZERO
         m = mofk(kl1); kl1s = kl1-2*m
         do klp1 = 1,kmax_phi_loc
!           ----------------------------------------------------------
            yp=>wfr_reg(j_inter:j_inter+n_inter-1,klp1,kl1)
            call PolyInterp(n_inter,xp,yp,r,yic,z_err_est)
            regf(kl1) = regf(kl1)+yic*ylmv(klp1)
!           ----------------------------------------------------------
            mp = mofk(klp1); klp1s = klp1-2*mp
            yp=>wfr_irr(j_inter:j_inter+n_inter-1,klp1s,kl1s)
            call PolyInterp(n_inter,xp,yp,r,yic,z_err_est)
            irrf(kl1) = irrf(kl1)+yic*ylmvc(klp1)
!           ----------------------------------------------------------
         enddo
      enddo
!
      gf(ic) = CZERO
      do kl2 = 1,kmax_kkr_loc
         do kl1 = 1,kmax_kkr_loc
            gf(ic) = gf(ic) + mat(kl1,kl2)*regf(kl1)*irrf(kl2)
         enddo
      enddo
      gf(ic) = -SQRTm1*kappa*gf(ic)
!
      if (present(dgf)) then
         call ErrorHandler('computeGF1','For calculating the Green function derivative, it is not implemented yet')
      endif
   enddo LOOP_ic
!
   nullify( wfr_reg, wfr_irr, mat, OmegaHat_mat, ylmv, ylmvc, yp, xp, regf, irrf )
!
   end subroutine computeGF1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeDOS(add_highl_fec,spin,site,atom)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use RadialGridModule, only : getGrid
!
   use SystemSymmetryModule, only : getSymmetryFlags
!
   implicit none
!
   logical, intent(in), optional :: add_highl_fec   ! Add high L contribution by free electron
!
   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: id, is, ic
   integer (kind=IntKind) :: jl, kl, klc, l, m, kl1, kl1p, kl2, kl2c, kl2p, kl2pc
   integer (kind=IntKind) :: ir, m1, m2, m2p, ma, i, l2
   integer (kind=IntKind) :: js, ia, nr, lm, jm, sz, loc1, loc2
   integer (kind=IntKind) :: lmax_dos, jmax_dos, kmax_kkr_loc, kmax_phi_loc, np
   integer (kind=IntKind), pointer :: dos_flags(:)
!
   complex (kind=CmplxKind), pointer :: p_tmp(:,:), q_tmp(:,:), sjm(:,:)
   complex (kind=CmplxKind), pointer :: dos(:,:), der_dos(:,:)
!
   complex (kind=CmplxKind) :: cfac, cmat
   complex (kind=CmplxKind), pointer :: green(:,:), der_green(:,:)
   complex (kind=CmplxKind), pointer :: wfr_reg(:,:,:)
   complex (kind=CmplxKind), pointer :: der_wfr_reg(:,:,:)
   complex (kind=CmplxKind), pointer :: Omega_mat(:,:)
   complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (.not. allocated(wks_dos)) then
      sz = 0
      do ia = 1, LocalNumSites
         nr = Scatter(ia)%numrs_cs
         lm = Scatter(ia)%lmax_green
         jm = (lm+1)*(lm+2)/2
         sz = sz + nr*jm*NumSpins*Scatter(ia)%NumSpecies
      enddo
      allocate( wks_dos(sz) )
      loc1 = 0
      do ia = 1, LocalNumSites
         nr = Scatter(ia)%numrs_cs
         lm = Scatter(ia)%lmax_green
         jm = (lm+1)*(lm+2)/2
         do js = 1, NumSpins
            do ic = 1, Scatter(ia)%NumSpecies
               loc2 = loc1 + nr*jm
               p1c => wks_dos(loc1+1:loc2)
               Scatter(ia)%Solutions(ic,js)%dos => aliasArray2_c(p1c,nr,jm)
               loc1 = loc2
            enddo
         enddo
      enddo
      if (rad_deriv) then
         allocate( wks_ddos(sz) )
         loc1 = 0
         do ia = 1, LocalNumSites
            nr = Scatter(ia)%numrs_cs
            lm = Scatter(ia)%lmax_green
            jm = (lm+1)*(lm+2)/2
            do js = 1, NumSpins
               do ic = 1, Scatter(ia)%NumSpecies
                  loc2 = loc1 + nr*jm
                  p1c => wks_ddos(loc1+1:loc2)
                  Scatter(ia)%Solutions(ic,js)%der_dos => aliasArray2_c(p1c,nr,jm)
                  loc1 = loc2
               enddo
            enddo
         enddo
      endif
   endif
!
   lmax_dos = Scatter(id)%lmax_green
   jmax_dos = (lmax_dos+1)*(lmax_dos+2)/2
   kmax_kkr_loc = Scatter(id)%kmax_kkr
   kmax_phi_loc = Scatter(id)%kmax_phi
   nr = Scatter(id)%numrs_cs
   Grid => getGrid(id)
!
   LOOP_ic: do ic = 1, Scatter(id)%NumSpecies
      if (present(atom)) then
         if (atom > 0 .and. ic /= atom) then
            cycle LOOP_ic
         endif
      endif
      dos => Scatter(id)%Solutions(ic,is)%dos
      if (rad_deriv) then
         der_dos => Scatter(id)%Solutions(ic,is)%der_dos
      endif
      Omega_mat => Scatter(id)%Solutions(ic,is)%Omega_mat
      OmegaHat_mat => Scatter(id)%Solutions(ic,is)%OmegaHat_mat
      if (abs(aimag(energy)) < TEN2m8) then
         wfr_reg => Scatter(id)%Solutions(ic,is)%reg_sol
         dos = CZERO
         if (rad_deriv) then
            der_wfr_reg => Scatter(id)%Solutions(ic,is)%reg_dsol
            der_dos = CZERO
         endif
!        if ( SphericalSolver .or. isSphPotZeroOutsideRmt ) then
         if (isSphericalSolverOn) then
            do kl1 = 1, kmax_kkr_loc
               do i = 1, nj3(kl1,kl1)
                  kl = kj3(i,kl1,kl1); m = mofk(kl); l = lofk(kl)
                  if (m >= 0 .and. l<= lmax_dos) then
                     jl = (l+1)*(l+2)/2 - l + m
                     cfac = cgnt(i,kl1,kl1)*Omega_mat(kl1,kl1) ! Here uses Omega, not OmegaHat!!!
                     do ir = 1, nr
                        dos(ir,jl) = dos(ir,jl) + cfac*wfr_reg(ir,kl1,kl1)**2
                     enddo
                     if (rad_deriv) then
                        do ir = 1, nr
                           der_dos(ir,jl) = der_dos(ir,jl) +                &
                               TWO*cfac*wfr_reg(ir,kl1,kl1)*der_wfr_reg(ir,kl1,kl1)
                        enddo
                     endif
                  endif
               enddo
            enddo
            cfac = kappa/PI
         else if (.false.) then
            p_tmp => aliasArray2_c(wks_plhat,Scatter(id)%numrs,kmax_phi_loc)
            if (rad_deriv) then
               q_tmp => aliasArray2_c(wks_qlhat,Scatter(id)%numrs,kmax_phi_loc)
            endif
            do kl2 = 1, kmax_kkr_loc
!              =======================================================
!!             p_tmp = CZERO
!!             do kl1 = 1, kmax_kkr_loc
!!                do kl1p = 1, kmax_phi_loc
!!                   do ir = 1, Scatter(id)%numrs
!!                      p_tmp(ir,kl1p) = p_tmp(ir,kl1p) + OmegaHat_mat(kl1,kl2)*wfr_reg(ir,kl1p,kl1)
!!                   enddo
!!                enddo
!!             enddo
!              -------------------------------------------------------
               call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                          CONE,wfr_reg,Scatter(id)%numrs*kmax_phi_loc,     &
                          OmegaHat_mat(1:kmax_kkr_loc,kl2),1,CZERO,p_tmp,1)
!              -------------------------------------------------------
               if (rad_deriv) then
!                 ----------------------------------------------------
                  call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                             CONE,der_wfr_reg,Scatter(id)%numrs*kmax_phi_loc, &
                             OmegaHat_mat(1:kmax_kkr_loc,kl2),1,CZERO,q_tmp,1)
!                 ----------------------------------------------------
               endif
               m2 = mofk(kl2); kl2c = kl2-2*m2
               do kl2p = 1, kmax_phi_loc
                  m2p = mofk(kl2p); kl2pc = kl2p-2*m2p
                  ma = abs(m2+m2p); ma = 1-2*mod(ma,2)
                  do kl1p = 1, kmax_phi_loc
                     do i = 1, nj3(kl2p,kl1p)
                        kl = kj3(i,kl2p,kl1p); m = mofk(kl); l = lofk(kl)
                        if (m >= 0 .and. l<= lmax_dos) then
                           cfac = ma*cgnt(i,kl2p,kl1p)
                           jl = (l+1)*(l+2)/2 - l + m
                           do ir = 1, nr
                              dos(ir,jl) = dos(ir,jl) + cfac*p_tmp(ir,kl1p)*wfr_reg(ir,kl2pc,kl2c)
                           enddo
                           if (rad_deriv) then
                              do ir = 1, nr
                                 der_dos(ir,jl) = der_dos(ir,jl) +                &
                                     cfac*q_tmp(ir,kl1p)*wfr_reg(ir,kl2pc,kl2c) + &
                                     cfac*p_tmp(ir,kl1p)*der_wfr_reg(ir,kl2pc,kl2c)
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
            cfac = kappa/PI
            nullify( p_tmp, q_tmp )
         else
            sjm => aliasArray2_c(wks_mtx1,kmax_kkr_loc,kmax_kkr_loc)
            do kl2 = 1, kmax_kkr_loc
               do kl1 = 1, kmax_kkr_loc
                  sjm(kl1,kl2)=kappa*OmegaHat_mat(kl1,kl2)-conjg(kappa*OmegaHat_mat(kl2,kl1))
               enddo
            enddo
            p_tmp => aliasArray2_c(wks_plhat,Scatter(id)%numrs,kmax_phi_loc)
            if (rad_deriv) then
               q_tmp => aliasArray2_c(wks_qlhat,Scatter(id)%numrs,kmax_phi_loc)
            endif
            np = mod(kmax_kkr_loc,NumPEsInGroup)
            do kl2 = MyPEinGroup+1, kmax_kkr_loc-np, NumPEsInGroup
!           do kl2 = 1, kmax_kkr_loc
!              -------------------------------------------------------
               call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                          CONE,wfr_reg,Scatter(id)%numrs*kmax_phi_loc,     &
                          sjm(1,kl2),1,CZERO,p_tmp,1)
!              -------------------------------------------------------
               if (rad_deriv) then
!                 ----------------------------------------------------
                  call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                             CONE,der_wfr_reg,Scatter(id)%numrs*kmax_phi_loc, &
                             sjm(1,kl2),1,CZERO,q_tmp,1)
!                 ----------------------------------------------------
               endif
               m2 = mofk(kl2); kl2c = kl2-2*m2
               do kl2p = 1, kmax_phi_loc
                  m2p = mofk(kl2p); kl2pc = kl2p-2*m2p
                  ma = abs(m2+m2p); ma = 1-2*mod(ma,2)
                  do kl1p = 1, kmax_phi_loc
                     do i = 1, nj3(kl2p,kl1p)
                        kl = kj3(i,kl2p,kl1p); m = mofk(kl); l = lofk(kl)
                        if (m >= 0 .and. l<= lmax_dos) then
                           cfac = ma*cgnt(i,kl2p,kl1p)
                           jl = (l+1)*(l+2)/2 - l + m
                           do ir = 1, nr
                              dos(ir,jl) = dos(ir,jl) + cfac*p_tmp(ir,kl1p)*wfr_reg(ir,kl2pc,kl2c)
                           enddo
                           if (rad_deriv) then
                              do ir = 1, nr
                                 der_dos(ir,jl) = der_dos(ir,jl) +                &
                                     cfac*q_tmp(ir,kl1p)*wfr_reg(ir,kl2pc,kl2c) + &
                                     cfac*p_tmp(ir,kl1p)*der_wfr_reg(ir,kl2pc,kl2c)
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
            if (NumPEsInGroup > 1) then
!              -------------------------------------------------------
               call GlobalSumInGroup(kGID,dos,nr,jmax_dos)
!              -------------------------------------------------------
               if (rad_deriv) then
!                 ----------------------------------------------------
                  call GlobalSumInGroup(kGID,der_dos,nr,jmax_dos)
!                 ----------------------------------------------------
               endif
            endif
            do kl2 = kmax_kkr_loc-np+1, kmax_kkr_loc
!              -------------------------------------------------------
               call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                          CONE,wfr_reg,Scatter(id)%numrs*kmax_phi_loc,     &
                          sjm(1,kl2),1,CZERO,p_tmp,1)
!              -------------------------------------------------------
               if (rad_deriv) then
!                 ----------------------------------------------------
                  call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                             CONE,der_wfr_reg,Scatter(id)%numrs*kmax_phi_loc, &
                             sjm(1,kl2),1,CZERO,q_tmp,1)
!                 ----------------------------------------------------
               endif
               m2 = mofk(kl2); kl2c = kl2-2*m2
               do kl2p = 1, kmax_phi_loc
                  m2p = mofk(kl2p); kl2pc = kl2p-2*m2p
                  ma = abs(m2+m2p); ma = 1-2*mod(ma,2)
                  do kl1p = 1, kmax_phi_loc
                     do i = 1, nj3(kl2p,kl1p)
                        kl = kj3(i,kl2p,kl1p); m = mofk(kl); l = lofk(kl)
                        if (m >= 0 .and. l<= lmax_dos) then
                           cfac = ma*cgnt(i,kl2p,kl1p)
                           jl = (l+1)*(l+2)/2 - l + m
                           do ir = 1, nr
                              dos(ir,jl) = dos(ir,jl) + cfac*p_tmp(ir,kl1p)*wfr_reg(ir,kl2pc,kl2c)
                           enddo
                           if (rad_deriv) then
                              do ir = 1, nr
                                 der_dos(ir,jl) = der_dos(ir,jl) +                &
                                     cfac*q_tmp(ir,kl1p)*wfr_reg(ir,kl2pc,kl2c) + &
                                     cfac*p_tmp(ir,kl1p)*der_wfr_reg(ir,kl2pc,kl2c)
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
            cfac = sqrtm1/PI2
            nullify( p_tmp, q_tmp, sjm )
         endif
      else if ( isIrrSolOn ) then
!        -------------------------------------------------------------
         call computeGF0()
!        -------------------------------------------------------------
         green => Scatter(id)%Solutions(ic,is)%green
         if (rad_deriv) then
            der_green => Scatter(id)%Solutions(ic,is)%der_green
         endif
         do jl = 1, jmax_dos
            kl = kofj(jl)
            m = mofj(jl)
            klc = kl-2*m
            do ir = 1, nr
               dos(ir,jl) = green(ir,kl) - m1m(m)*conjg(green(ir,klc))
            enddo
            if (rad_deriv) then
               do ir = 1, nr
                  der_dos(ir,jl) = der_green(ir,kl) - m1m(m)*conjg(der_green(ir,klc))
               enddo
            endif
         enddo
         cfac = sqrtm1/PI2
      else
         call ErrorHandler('computeDOS','For energy in the complex plane, the irregular solution is required')
      endif
!
      dos = cfac*dos
      if (rad_deriv) then
         der_dos = cfac*der_dos
      endif
!
      if ( present(add_highl_fec) ) then
         if (add_highl_fec) then
            TmpSpace = CZERO
            do l = lofk(kmax_kkr_loc), 0, -1
               do ir = 1, nr
                  TmpSpace(ir) = TmpSpace(ir) + (2*l+1)*bjl(ir,l)**2/kappa**2
               enddo
            enddo
            cfac = kappa*sqrt(PI)/(2.0d0*PI**2) ! = kappa/(4.0d0*PI**2)/Y00
            do ir = 1, nr
               dos(ir,1) = dos(ir,1) + cfac*(Grid%r_mesh(ir)**2-TmpSpace(ir))
            enddo
            if (rad_deriv) then
               TmpSpace = CZERO
               do l = lofk(kmax_kkr_loc), 0, -1
                  do ir = 1, nr
                     TmpSpace(ir) = TmpSpace(ir) + (2*l+1)*dbjl(ir,l)**2/kappa**2
                  enddo
               enddo
               do ir = 1, nr
                  der_dos(ir,1) = der_dos(ir,1) - cfac*TmpSpace(ir)
               enddo
            endif
         endif
      endif
!
!     ================================================================
!     Symmetrizing the DOS, if needed.  Added by Yang on 09-21-2018
!     ================================================================
      if (isDosSymmOn) then
         dos_flags => getSymmetryFlags(id)
         do jl = 1, jmax_dos
            if (dos_flags(jl) == 0) then
               dos(:,jl) = CZERO
               if (rad_deriv) then
                  der_dos(:,jl) = CZERO
               endif
            endif
         enddo
      endif
!     ================================================================
   enddo LOOP_ic
!
   nullify( wfr_reg, der_wfr_reg, dos, der_dos, OmegaHat_mat, Omega_mat )
   nullify( green, der_green, p_tmp, q_tmp )
!
   end subroutine computeDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computePDOS(spin,site,atom)
!  ===================================================================
   use SystemSymmetryModule, only : getSymmetryFlags
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use RadialGridModule, only : getGrid
!
   implicit none
!
   integer (kind=Intkind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind) :: id, is, ic
   integer (kind=IntKind) :: ir, jl, kl, klc, l, m, kl1, kl1p, kl2, kl2c, kl2p, kl2pc, m2, m2p, ma, i
   integer (kind=IntKind) :: js, ia, nr, lm, jm, sz, loc1, loc2, m1
   integer (kind=IntKind) :: lmax_dos, jmax_dos
   integer (kind=IntKind) :: kmax_kkr_loc, kmax_phi_loc
   integer (kind=IntKind), pointer :: dos_flags(:)
!
   real (kind=RealKind) :: dos_cell, dos_mt
!
   complex (kind=CmplxKind), pointer :: pdos(:,:,:), p_tmp(:,:), sjm(:,:), pdos_kl(:,:)
!
   complex (kind=CmplxKind) :: cfac, cmat
   complex (kind=CmplxKind), pointer :: green(:,:)
   complex (kind=CmplxKind), pointer :: wfr_reg(:,:,:)
   complex (kind=CmplxKind), pointer :: Omega_mat(:,:)
   complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   if (present(spin)) then
      is = min(NumSpins,spin)
   else
      is = LocalSpin
   endif
!
   if (present(site)) then
      id = site
   else
      id = LocalIndex
   endif
!
   if (.not. allocated(wks_pdos)) then
      sz = 0
      do ia = 1, LocalNumSites
         nr = Scatter(ia)%numrs_cs
         lm = Scatter(ia)%lmax_green
         jm = (lm+1)*(lm+2)/2
         sz = sz + nr*jm*Scatter(ia)%kmax_kkr*NumSpins*Scatter(ia)%NumSpecies
      enddo
      allocate( wks_pdos(sz) )
      loc1 = 0
      do ia = 1, LocalNumSites
         nr = Scatter(ia)%numrs_cs
         lm = Scatter(ia)%lmax_green
         jm = (lm+1)*(lm+2)/2
         do js = 1, NumSpins
            do ic = 1, Scatter(ia)%NumSpecies
               loc2 = loc1 + nr*jm*Scatter(ia)%kmax_kkr
               p1c => wks_pdos(loc1+1:loc2)
               Scatter(ia)%Solutions(ic,js)%pdos =>                   &
                           aliasArray3_c(p1c,nr,jm,Scatter(ia)%kmax_kkr)
               loc1 = loc2
            enddo
         enddo
      enddo
   endif
!
   lmax_dos = Scatter(id)%lmax_green
   jmax_dos = (lmax_dos+1)*(lmax_dos+2)/2
   kmax_kkr_loc = Scatter(id)%kmax_kkr
   kmax_phi_loc = Scatter(id)%kmax_phi
   nr = Scatter(id)%numrs_cs
   Grid => getGrid(id)
!
   if (.not.allocated(space_integrated_dos_cell)) then
      allocate( space_integrated_dos_cell(MaxSpecies,NumSpins,LocalNumSites),       &
                space_integrated_dos_mt(MaxSpecies,NumSpins,LocalNumSites),         &
                space_integrated_pdos_cell(kmax_max_phi,MaxSpecies,NumSpins,LocalNumSites), &
                space_integrated_pdos_mt(kmax_max_phi,MaxSpecies,NumSpins,LocalNumSites) )
   endif
!
   LOOP_ic: do ic = 1, Scatter(id)%NumSpecies
      if (present(atom)) then
         if (atom > 0 .and. ic /= atom) then
            cycle LOOP_ic
         endif
      endif
      pdos => Scatter(id)%Solutions(ic,is)%pdos
      Omega_mat => Scatter(id)%Solutions(ic,is)%Omega_mat
      OmegaHat_mat => Scatter(id)%Solutions(ic,is)%OmegaHat_mat
      if (abs(aimag(energy)) < TEN2m8) then
         wfr_reg => Scatter(id)%Solutions(ic,is)%reg_sol
!        if ( SphericalSolver .or. isSphPotZeroOutsideRmt ) then
         if (isSphericalSolverOn) then
            do kl1 = 1, kmax_kkr_loc
               do i = 1, nj3(kl1,kl1)
                  kl = kj3(i,kl1,kl1); m = mofk(kl); l = lofk(kl)
                  if (m >= 0 .and. l<= lmax_dos) then
                     jl = (l+1)*(l+2)/2 - l + m
                     cfac = cgnt(i,kl1,kl1)*Omega_mat(kl1,kl1) ! Here uses Omega, not OmegaHat!!!
                     do ir = 1, nr
                        pdos(ir,jl,kl1) = cfac*wfr_reg(ir,kl1,kl1)**2
                     enddo
                  endif
               enddo
            enddo
            cfac = kappa/PI
!        else if (kmax_kkr_loc == kmax_phi_loc) then
         else if (.false.) then
            pdos = CZERO
            p_tmp => aliasArray2_c(wks_plhat,Scatter(id)%numrs,kmax_phi_loc)
            do kl2 = 1, kmax_kkr_loc
               m2 = mofk(kl2); kl2c = kl2-2*m2
!              =======================================================
!!             p_tmp = CZERO
!!             do kl1 = 1, kmax_kkr_loc
!!                do kl1p = 1, kmax_phi_loc
!!                   do ir = 1, Scatter(is,id)%numrs
!!                      p_tmp(ir,kl1p) = p_tmp(ir,kl1p) + OmegaHat_mat(kl1,kl2)*wfr_reg(ir,kl1p,kl1)
!!                   enddo
!!                enddo
!!             enddo
!              -------------------------------------------------------
               call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                          CONE,wfr_reg,Scatter(id)%numrs*kmax_phi_loc,     &
                          OmegaHat_mat(1:kmax_kkr_loc,kl2),1,CZERO,p_tmp,1)
!              -------------------------------------------------------
               do kl2p = 1, kmax_phi_loc
                  m2p = mofk(kl2p); kl2pc = kl2p-2*m2p
                  ma = abs(m2+m2p); ma = 1-2*mod(ma,2)
                  do kl1p = 1, kmax_phi_loc
                     do i = 1, nj3(kl2p,kl1p)
                        kl = kj3(i,kl2p,kl1p); m = mofk(kl); l = lofk(kl)
                        if (m >= 0 .and. l<= lmax_dos) then
                           cfac = ma*cgnt(i,kl2p,kl1p)
                           jl = (l+1)*(l+2)/2 - l + m
                           do ir = 1, nr
                              pdos(ir,jl,kl2) = pdos(ir,jl,kl2) + cfac*p_tmp(ir,kl1p)*wfr_reg(ir,kl2pc,kl2c)
                           enddo
                        endif
                     enddo
                  enddo
               enddo
            enddo
            cfac = kappa/PI
            nullify( p_tmp )
         else
            sjm => aliasArray2_c(wks_mtx1,kmax_kkr_loc,kmax_kkr_loc)
            do kl2 = 1, kmax_kkr_loc
               do kl1 = 1, kmax_kkr_loc
                  sjm(kl1,kl2)=kappa*OmegaHat_mat(kl1,kl2)-conjg(kappa*OmegaHat_mat(kl2,kl1))
               enddo
            enddo
            pdos = CZERO
            p_tmp => aliasArray2_c(wks_plhat,Scatter(id)%numrs,kmax_phi_loc)
            do kl2 = 1, kmax_kkr_loc
               m2 = mofk(kl2); kl2c = kl2-2*m2
!              -------------------------------------------------------
               call zgemv('n',Scatter(id)%numrs*kmax_phi_loc,kmax_kkr_loc, &
                          CONE,wfr_reg,Scatter(id)%numrs*kmax_phi_loc,     &
                          sjm(1:kmax_kkr_loc,kl2),1,CZERO,p_tmp,1)
!              -------------------------------------------------------
               do kl2p = 1, kmax_phi_loc
                  m2p = mofk(kl2p); kl2pc = kl2p-2*m2p
                  ma = abs(m2+m2p); ma = 1-2*mod(ma,2)
                  do kl1p = 1, kmax_phi_loc
                     do i = 1, nj3(kl2p,kl1p)
                        kl = kj3(i,kl2p,kl1p); m = mofk(kl); l = lofk(kl)
                        if (m >= 0 .and. l<= lmax_dos) then
                           cfac = ma*cgnt(i,kl2p,kl1p)
                           jl = (l+1)*(l+2)/2 - l + m
                           do ir = 1, nr
                              pdos(ir,jl,kl2) = pdos(ir,jl,kl2) + cfac*p_tmp(ir,kl1p)*wfr_reg(ir,kl2pc,kl2c)
                           enddo
                        endif
                     enddo
                  enddo
               enddo
            enddo
            cfac = sqrtm1/PI2
            nullify( p_tmp, sjm )
         endif
      else
         call ErrorHandler('computePDOS','Not implemented for non-real energy yet')
      endif
      pdos = cfac*pdos
!
!     ================================================================
!     Symmetrizing the DOS, if needed.  Added by Yang on 09-21-2018
!     ================================================================
      if (isDosSymmOn) then
         dos_flags => getSymmetryFlags(id)
         do kl = 1, kmax_kkr_loc
            do jl = 1, jmax_dos
               if (dos_flags(jl) == 0) then
                  pdos(:,jl,kl) = CZERO
               endif
            enddo
         enddo
      endif
!     ================================================================
!
      space_integrated_dos_cell(ic,is,id) = ZERO
      space_integrated_dos_mt(ic,is,id) = ZERO
      do kl = kmax_kkr_loc, 1, -1
         pdos_kl => pdos(:,:,kl)
!        -------------------------------------------------------------
         dos_cell = getVolumeIntegration(id, nr, Grid%r_mesh, jmax_dos, 2, pdos_kl, dos_mt)
!        -------------------------------------------------------------
         space_integrated_pdos_cell(kl,ic,is,id) = dos_cell
         space_integrated_pdos_mt(kl,ic,is,id) = dos_mt
         space_integrated_dos_cell(ic,is,id) = space_integrated_dos_cell(ic,is,id) + dos_cell
         space_integrated_dos_mt(ic,is,id) = space_integrated_dos_mt(ic,is,id) + dos_mt
      enddo
   enddo LOOP_ic
!
   nullify( wfr_reg, pdos, Omega_mat, OmegaHat_mat, pdos_kl, p_tmp )
!
   end subroutine computePDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntSphHankelSq0(lmax,rc,energy,fint)
!  ===================================================================
   use BesselModule, only : SphericalBessel, SphericalNeumann
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: rc
!
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind), intent(out) :: fint(0:lmax)
   complex (kind=CmplxKind) :: bessjl(0:lmax+1), bessnl(0:lmax+1), x, besshl, besshlp1, besshlm1
   complex (kind=CmplxKind) :: Ul, Ulp1
!
   x = rc*sqrt(energy)+SQRTm1*0.00001d0
!  -------------------------------------------------------------------
   call SphericalBessel(lmax+1,x,bessjl)
   call SphericalNeumann(lmax+1,x,bessnl)
!  -------------------------------------------------------------------
   besshl = bessjl(0)+SQRTm1*bessnl(0)
   besshlp1 = bessjl(1)+SQRTm1*bessnl(1)
   fint(0) = (besshl*besshlp1/x-besshl**2-besshlp1**2)*HALF
   do l = 1, lmax
      besshlm1 = besshl
      besshl = besshlp1
      besshlp1 = bessjl(l+1)+SQRTm1*bessnl(l+1)
!     fint(l) = ((2*l+1)*besshl*besshlp1/x-besshl**2-besshlp1**2)*HALF
      fint(l) = (besshlm1*besshlp1-besshl**2)*HALF
   enddo
!
!  ===================================================================
!  Another way of calculating fint is to use recursive relation.
!  Both ways have shown to give the same results.
!  ===================================================================
!! Ul = -SQRTm1*exp(SQRTm1*TWO*x)/(TWO*x**3)
!! write(6,'(a,2d16.8,2x,2d16.8)')'fint(0), Ul = ',fint(0), Ul
!! fint(0) = Ul
!! do l = 1, lmax
!!    besshlm1 = bessjl(l-1)+SQRTm1*bessnl(l-1)
!!    besshl = bessjl(l)+SQRTm1*bessnl(l)
!!    Ulp1 = (besshl**2+besshlm1**2+(2*l+1)*fint(l-1))/(2*l-1.0d0)
!!    write(6,'(a,2d16.8,2x,2d16.8)')'fint(l), Ulp1 = ',fint(l), Ulp1
!!    fint(l) = Ulp1
!! enddo
!  ===================================================================
!
   do l = 0, lmax
      fint(l) = fint(l)*rc**3
   enddo
!
   end subroutine calIntSphHankelSq0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntSphHankelSq1(id,lmax,rm,rc,energy,fint)
!  ===================================================================
   use BesselModule, only : SphericalBessel, SphericalNeumann
   use StepFunctionModule, only : getNumGaussRs, getGaussR, getGaussRWeight, &
                                  getSigmaLL
!
   implicit none
!
   integer (kind=intKind), intent(in) :: id, lmax
   integer (kind=intKind) :: ng, ig, kl, klp, kmax, l, lp
!     
   real (kind=RealKind), intent(in) :: rm, rc
   real (kind=RealKind), pointer :: rg(:), wg(:)
!  
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind), intent(out) :: fint(:,:)
   complex (kind=CmplxKind) :: bessjl(0:lmax+1), bessnl(0:lmax+1)
   complex (kind=CmplxKind), pointer :: sig_LL(:,:,:)
   complex (kind=CmplxKind) :: x, besshl, besshlp
   complex (kind=CmplxKind) :: Ul, Ulp1
!  
   kmax = (lmax+1)**2
   ng = getNumGaussRs(id)
   rg => getGaussR(id)
   wg => getGaussRWeight(id)
   sig_LL => getSigmaLL(id)
!  
   fint = CZERO
   do ig = 1, ng 
      if (rg(ig) < rm .or. rg(ig) > rc) then
         write(6,'(a,3f12.5)')'rg,rm,rc = ',rg(ig),rm,rc
         call ErrorHandler('calIntSphHankelSq','rg is out of range')
      endif
      x = rg(ig)*sqrt(energy)+SQRTm1*0.00001d0
!     ----------------------------------------------------------------
      call SphericalBessel(lmax+1,x,bessjl)
      call SphericalNeumann(lmax+1,x,bessnl)
!     ----------------------------------------------------------------
      do kl = 1, kmax
         l = lofk(kl)
         besshl = bessjl(l)+SQRTm1*bessnl(l)
         do klp = 1, kmax
            lp = lofk(klp)
            besshlp = bessjl(lp)+SQRTm1*bessnl(lp)
!           if (klp == kl) then
!              fint(klp,kl) = fint(klp,kl)+(sig_LL(klp,kl,ig)-CONE)*besshl*besshlp*wg(ig)*rg(ig)**2
!           else
!              fint(klp,kl) = fint(klp,kl)+sig_LL(klp,kl,ig)*besshl*besshlp*wg(ig)*rg(ig)**2
!           endif
            fint(klp,kl) = fint(klp,kl)+sig_LL(klp,kl,ig)*besshl*besshlp*wg(ig)*rg(ig)**2
         enddo
      enddo
   enddo
!  
   end subroutine calIntSphHankelSq1
!  ===================================================================
end module SSSolverModule
