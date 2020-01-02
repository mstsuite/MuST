module ScfDataModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicParamDefinitionsModule, only : MaxLenFileName
!
public :: initScfData,                 &
          endScfData,                  &
          isNonRelativisticCore,       &
          isNonRelativisticValence,    &
          isScalarRelativisticValence, &
          isRelativisticValence,       &
          getPotentialTypeParam,       &
          getSingleSiteSolverType,     &
          setSingleSiteSolverType,     &
          getSSSIntegrationMethod,     &
          getSSSolutionsMethod,        &
          getSingleSiteSolverMethod,   &
          getLmaxPotSolver,            &
          getLmaxSolution,             &
          getDOSrunid,                 &
          isSSIrregularSolOn,          &
          isFittedChargeDen,           &
          isSSPhaseShiftOn,            &
          isEfIterateOn,               &
          getCantedMomentTorqueFactor, &
          isReadKmesh,                 &
          isReadEmesh,                 &
          isLSMS,                      &
          isRealAxisSS,                &
          isSingleSiteEContour,        &
          isScreenKKR,                 &
          isKKR,                       &
          isKKRCPA,                    &
          isEmbeddedCluster,           &
          isScreenKKR_LSMS,            &
          isSingleSite,                &
          isFrozenCore,                &
          isLloyd,                     &
          getLloydMode,                &
          getKmeshFileName,            &
          getEmeshFileName,            &
          isExchangeParamNeeded,       &
          isSimpleMixing,              &
          isDGAMixing,                 &
          isBroydenMixing,             &
          isPotentialMixing,           &
          isChargeMixing,              &
          isDirectKSpaceSolver,        &
          isCheckSuperLUSolver,        &
          isInterstitialElectronPolarized, &
          isLdaCorrectionNeeded,       &
          isChargeSymm,                &
          getAdaptiveIntegrationMethod,&
          getUJfile,                   &
          setSCFMethod,                &
          getPoleSearchStep,           &
          retreiveEffectiveMediumParams,   &
          printScfData
!
public
   character (len=50) :: istop
   character (len=50) :: inputpath
   character (len=50) :: inputfile
   character (len=50) :: info_table
   character (len=50) :: excorr_name
!  character (len=50) :: info_path
!
   integer (kind=IntKind) :: nscf  = 0
   integer (kind=IntKind) :: nspin = 0
   integer (kind=IntKind) :: Tinv_alg = -1
   integer (kind=IntKind) :: Minv_alg = -1
   integer (kind=IntKind) :: ngaussr  = 0
   integer (kind=IntKind) :: ngaussq  = 0
   integer (kind=IntKind) :: nrelc  = -1
   integer (kind=IntKind) :: nrelv  = -1
   integer (kind=IntKind) :: harris = -1
   integer (kind=IntKind) :: potwrite = 0
   integer (kind=IntKind) :: movie    = 0
   integer (kind=IntKind) :: ntstep   = 0
   integer (kind=IntKind) :: ContourType = -1
   integer (kind=IntKind) :: NumEs = 0
   integer (kind=IntKind) :: NumExtraEs = 0
   integer (kind=IntKind) :: NumSS_IntEs = 0
   integer (kind=IntKind) :: OffsetE = -1
   integer (kind=IntKind) :: eGridType = 0
   integer (kind=IntKind) :: iterateEf = 0
   integer (kind=IntKind) :: RealAxisSS = 0
   integer (kind=IntKind) :: SingleSiteEContour = 0
   integer (kind=IntKind) :: KSpaceSolver = 0
!
   integer (kind=IntKind) :: i_vdif = 1
   integer (kind=IntKind) :: n_spin_cant = 0
   integer (kind=IntKind) :: n_spin_pola = 0
   integer (kind=IntKind) :: ss_solver_type = -1
   integer (kind=IntKind) :: lmax_pot_solver = -1
   integer (kind=IntKind) :: lmax_sol_cutoff = -1
   integer (kind=IntKind) :: ss_irrsol = 0
   integer (kind=IntKind) :: ss_phaseshift = 0
   integer (kind=IntKind) :: pot_type = -1
!
   integer (kind=IntKind) :: mixing_type = 0
   integer (kind=IntKind) :: mixing_quantity = 0
!
   real (kind=RealKind) :: Temperature
!
   real (kind=RealKind) :: EvBottom
   real (kind=RealKind) :: ErBottom
   real (kind=RealKind) :: ErTop
   real (kind=RealKind) :: EiBottom
   real (kind=RealKind) :: EiTop
!
   real (kind=RealKind) :: etol
   real (kind=RealKind) :: ptol
   real (kind=RealKind) :: eftol
   real (kind=RealKind) :: slutol
   real (kind=RealKind) :: sluchecktol
   real (kind=RealKind) :: rmstol
   real (kind=RealKind) :: tstep
   real (kind=RealKind) :: pole_step
!
   integer (kind=IntKind), parameter :: LeadPE = 0
!
   integer (kind=IntKind) :: NumKMeshs
   integer (kind=IntKind) :: kGenScheme
   integer (kind=IntKind), allocatable :: Kdiv(:,:)
   integer (kind=IntKind) :: ChargeSymmetry = 0
   integer (kind=IntKind) :: Symmetrize
   integer (kind=IntKind) :: Lloyd = 0
   integer (kind=IntKind) :: LloydMode = 1
   integer (kind=IntKind) :: TableID
!
   integer (kind=IntKind), parameter, private :: readEmesh   = 1
   integer (kind=IntKind), parameter, private :: readKmesh   = 1
   integer (kind=IntKind), parameter, private :: SingleSite  = -2
   integer (kind=IntKind), parameter, private :: ScreenKKR_LSMS = -1
   integer (kind=IntKind), parameter, private :: ScreenKKR = 0
   integer (kind=IntKind), parameter, private :: LSMS = 1
   integer (kind=IntKind), parameter, private :: KKR = 2
   integer (kind=IntKind), parameter, private :: KKRCPA = 3
   integer (kind=IntKind), parameter, private :: EmbeddedCluster = 4
!
   integer (kind=IntKind), private :: read_emesh = 0
   integer (kind=IntKind), private :: read_kmesh = 0
   integer (kind=IntKind), private :: scf_method = -10
   integer (kind=IntKind), private :: scf_method_save = -10
   integer (kind=IntKind), private :: LdaCorrectionType = 0
   integer (kind=IntKind), private :: DOSrunid   = 0
   integer (kind=IntKind), private :: RealEIntMethod = 0
!
   real (kind=RealKind), private :: ctq
!
   character (len=1),  private :: j_ij
   character (len=50), private :: UJfile
!
   integer (kind=IntKind), private :: ss_integral_method = -1
   integer (kind=IntKind), private :: ss_solution_method = -1
   integer (kind=IntKind), private :: sss_method = 2
!
   logical, private :: FrozenCore = .false.
   logical, private :: FrozenCoreFile_specified = .false.
   integer (kind=IntKind), private :: nFrozenCore = 0
   character (len=50), private :: FrozenCoreFile_name = ' '
!
!  ===================================================================
!  Effective medium parameters
!  ===================================================================
   integer (kind=IntKind), private :: EM_mix_type=2
   integer (kind=IntKind), private :: EM_max_iter = 30
   real (kind=RealKind), private ::   EM_mix_0 = 0.1d0
   real (kind=RealKind), private ::   EM_mix_1 = 0.01d0
   real (kind=RealKind), private ::   EM_tol = 0.0000001d0
   real (kind=RealKind), private ::   EM_switch = 0.003
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initScfData(tbl_id)
!  ===================================================================
   use MathParamModule, only : TEN2m6, TEN2m8, ZERO, TEN2m4, ONE
   use PublicParamDefinitionsModule, only : ASA, MuffinTin, MuffinTinASA
   use InputModule, only : getKeyValue
   implicit none
!
   integer (kind=IntKind), intent(in) :: tbl_id
!
   integer (kind=IntKind) :: rstatus, n
!
   character (len=50) :: s50
!
   real (kind=RealKind) :: rp(3)
!
   TableID = tbl_id
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Current File Path',inputpath)
   rstatus = getKeyValue(tbl_id,'Current File Name',inputfile)
!  rstatus = getKeyValue(tbl_id,'Path to Info Table File',info_path)
   rstatus = getKeyValue(tbl_id,'Info Table File Name',info_table)
   rstatus = getKeyValue(tbl_id,'No. Iterations (> 0)',nscf)
   rstatus = getKeyValue(tbl_id,'Method of SCF Calculation',scf_method_save)
   rstatus = getKeyValue(tbl_id,'K-space Solver Method',KSpaceSolver )
   rstatus = getKeyValue(tbl_id,'Stop-at Routine Name',istop)
   rstatus = getKeyValue(tbl_id,'No. Iter for Each Pot. Write',potwrite)
   rstatus = getKeyValue(tbl_id,'No. Iter for Each Movie',movie)
   rstatus = getKeyValue(tbl_id,'Calc. Harris Energy (H.E.)',harris)
   rstatus = getKeyValue(tbl_id,'No. Gauss Pts. along r',ngaussr)
   rstatus = getKeyValue(tbl_id,'No. Gauss Pts. along theta',ngaussq)
   rstatus = getKeyValue(tbl_id,'Valence Band Bottom Est.',EvBottom)
   rstatus = getKeyValue(tbl_id,'Temperature Parameter (K)',Temperature)
   rstatus = getKeyValue(tbl_id,'DOS Run ID',DOSrunid)
   rstatus = getKeyValue(tbl_id,'Energy (Ryd) Tol (> 0)',etol)
   rstatus = getKeyValue(tbl_id,'Potential Tol (> 0)',ptol)
   rstatus = getKeyValue(tbl_id,'Fermi Energy Tol (> 0)',eftol)
   rstatus = getKeyValue(tbl_id,'SuperLU Tol (> 0)',slutol)
   rstatus = getKeyValue(tbl_id,'K-space Check Tol (> 0)',sluchecktol)
   rstatus = getKeyValue(tbl_id,'Other RMS Tol (> 0)',rmstol)
   rstatus = getKeyValue(tbl_id,'Val. Electron Rel (>= 0)',nrelv)
   rstatus = getKeyValue(tbl_id,'Core Electron Rel (>= 0)',nrelc)
   rstatus = getKeyValue(tbl_id,'Charge Symmetry (>=0)',ChargeSymmetry)
   rstatus = getKeyValue(tbl_id,'SS Integration Method (>=0)',ss_integral_method)
   rstatus = getKeyValue(tbl_id,'Potential Type (>= 0)',pot_type)
!  rstatus = getKeyValue(tbl_id,'Single Site Solver (>= 0)',ss_solver_type)
   if (getKeyValue(tbl_id,'Single Site Solver (>= 0)',ss_solver_type,default_param=.false.) /= 0) then
      if (nrelv < 2) then
         if (pot_type == ASA .or. pot_type == MuffinTin .or. pot_type == MuffinTinASA) then
            ss_solver_type = 0
         else
            ss_solver_type = 2
         endif
      else
         if (pot_type == ASA .or. pot_type == MuffinTin .or. pot_type == MuffinTinASA) then
            ss_solver_type = 1
         else
            ss_solver_type = 3
         endif
      endif
   endif
!
   rstatus = getKeyValue(tbl_id,'SS Solutions Method (>=0)',ss_solution_method)
   if ( getKeyValue(tbl_id,'Single Site Solution Method (>=-1)',sss_method) /= 0) then
      if (ss_integral_method == 0) then
         sss_method = -1
      else if (ss_integral_method == 1 .and. ss_solution_method == 0) then
         sss_method = 0
      else if (ss_integral_method == 1 .and. ss_solution_method == 1) then
         sss_method = 1
      else
         sss_method = 2
      endif
   endif
   if ( getKeyValue(tbl_id,'SS Lmax Potential Solver',lmax_pot_solver) /= 0) then
      lmax_pot_solver = -1
   endif
   if ( getKeyValue(tbl_id,'Solutions Lmax Cutoff',lmax_sol_cutoff) /= 0) then
      lmax_sol_cutoff = -1
   endif
!
   if (getKeyValue(tbl_id,'Irregular Solutions (>=0)',ss_irrsol,default_param=.false.) /= 0) then
      if (getKeyValue(tbl_id,'Irregular Solutions (>=0)',ss_irrsol,default_param=.true.) == 0) then
         if (ss_irrsol == 0 .and. pot_type /= 3 .and. pot_type /= 5 .and. pot_type /= 6) then
            ss_irrsol = 1  ! In case of default case and muffin-tin potential, enable irr. sol. calc.
         endif
      else
         call ErrorHandler('initScfData','Input parameter is not found','Irregular Solutions (>=0)')
      endif
   endif
   if (ss_irrsol==1) then
!     necessary to account for the boundary conditions of irregular solutions
      ss_solution_method = 1
      sss_method = 2
   endif
   rstatus = getKeyValue(tbl_id,'Pole Search Step (>0.0)',pole_step)
   if ( rstatus /= 0 .or. pole_step < TEN2m6 .or. pole_step >= 0.05d0) then
      pole_step = 0.010d0
   endif
   rstatus = getKeyValue(tbl_id,'Compute Phase Shifts (>=0)',ss_phaseshift)
   rstatus = getKeyValue(tbl_id,'Exch-Corr. LDA Type (>= 0)',excorr_name)
   rstatus = getKeyValue(tbl_id,'Spin Index Param (>= 1)',nspin)
   if (nspin > 1) then
      rstatus = getKeyValue(tbl_id,'Interstitial Electron Spin',i_vdif)
   endif
   rstatus = getKeyValue(tbl_id,'Canted Moment Torque Coef.',ctq)
   rstatus = getKeyValue(tbl_id,'Calculate J_ij (y/n)',j_ij)
   rstatus = getKeyValue(tbl_id,'T-matrix inversion (>= 0)',Tinv_alg)
   rstatus = getKeyValue(tbl_id,'M-matrix inversion (>= 0)',Minv_alg)
   rstatus = getKeyValue(tbl_id,'No. Time Steps (>= 0)',ntstep)
   rstatus = getKeyValue(tbl_id,'Time Step',tstep)
   rstatus = getKeyValue(tbl_id,'Time Step',tstep)
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Read E-mesh from emeshs.inp',read_emesh)
   rstatus = getKeyValue(tbl_id,'Contour Type (>= 0)',ContourType)
   rstatus = getKeyValue(tbl_id,'Energy Grid Type (>= 0)',eGridType)
   rstatus = getKeyValue(tbl_id,'No. Energy Grids',NumEs)
   rstatus = getKeyValue(tbl_id,'No. Extra Energy Points',NumExtraEs)
   rstatus = getKeyValue(tbl_id,'SS Real Axis Int. Points',NumSS_IntEs)
   if (rstatus /= 0) then
      NumSS_IntEs = 500
   endif
   rstatus = getKeyValue(tbl_id,'SS Real Axis Int. Method',RealEIntMethod)
   if (rstatus /= 0) then
      RealEIntMethod = 0
   endif
   rstatus = getKeyValue(tbl_id,'Offset Energy Point',OffsetE)
   rstatus = getKeyValue(tbl_id,'Real Axis Bottom, erbot',ErBottom)
   rstatus = getKeyValue(tbl_id,'Real Axis Top, ertop',ErTop)
   rstatus = getKeyValue(tbl_id,'Imag Axis Bottom, eibot',EiBottom)
   rstatus = getKeyValue(tbl_id,'Imag Axis Top, eitop',EiTop)
   rstatus = getKeyValue(tbl_id,'Iterate Fermi energy',iterateEf)
   rstatus = getKeyValue(tbl_id,'Real Axis Single Site',RealAxisSS)
   rstatus = getKeyValue(tbl_id,'Single Site Contour',SingleSiteEContour)
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Read K-mesh from kmeshs.inp',read_kmesh)
   rstatus = getKeyValue(tbl_id,'Scheme to Generate K (>=0)',kGenScheme)
   rstatus = getKeyValue(tbl_id,'No. K Meshs in IBZ (> 0)',NumKMeshs)
   rstatus = getKeyValue(tbl_id,'Symmetrize BZ Integration',Symmetrize)
   rstatus = getKeyValue(tbl_id,'Mixing quantity type',mixing_quantity)
   rstatus = getKeyValue(tbl_id,'Mixing algorithm',mixing_type)
   rstatus = getKeyValue(tbl_id,'Lloyd correction',Lloyd)
   rstatus = getKeyValue(tbl_id,'Lloyd mode',LloydMode)
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'LDA Improvement Scheme',LdaCorrectionType)
!  -------------------------------------------------------------------
   if (NumKMeshs > 0) then
      allocate(Kdiv(3,NumKMeshs))
      rstatus = getKeyValue(tbl_id,'Kx, Ky, Kz Division (> 0)',3,Kdiv,NumKMeshs)
   endif
!
!  if (abs(ErTop) < TEN2m8) then
!     ErTop = -1.0d0
!  endif
!  -------------------------------------------------------------------
   if (nspin == 3) then
      n_spin_cant=2
      n_spin_pola=2
!     Symmetrize = 0
   else if (nspin == 2) then
      n_spin_pola=2
      n_spin_cant=1
   else if (nspin == 1) then
      n_spin_pola=1
      n_spin_cant=1
   else
      call ErrorHandler('initScfData','invalid nspin',nspin)
   endif
!
   if (i_vdif < 1 .or. i_vdif > 2) then
      call WarningHandler('initScfData','Unknown interstitial spin index', &
                          i_vdif)
      i_vdif = 1
   endif
!
   scf_method = scf_method_save
   if (scf_method==ScreenKKR_LSMS) then
      scf_method = LSMS ! first energy points are done with LSMS 
                        ! the switch to KKR is controled in the energy loop
   endif 
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Frozen-Core Calculation',n)
!  -------------------------------------------------------------------
   if (rstatus /= 0) then
      FrozenCore = .false.
      nFrozenCore = nscf
   else if (n < 1) then
      FrozenCore = .false.
      nFrozenCore = nscf
   else
      FrozenCore = .true.
      nFrozenCore = n
   endif
!
   if (getKeyValue(tbl_id,'Frozen-Core File Name',s50) == 0) then
      FrozenCoreFile_name = s50
      FrozenCoreFile_specified = .true.
   else
      FrozenCoreFile_name = ' '
      FrozenCoreFile_specified = .false.
   endif
!
   rstatus = getKeyValue(tbl_id,'Effective Medium Mixing Scheme',EM_mix_type)
   rstatus = getKeyValue(tbl_id,'Maximum Effective Medium Iterations',EM_max_iter)
   if ( getKeyValue(tbl_id,'Effective Medium Mixing Parameters',2,rp) == 0) then
      EM_mix_0 = rp(1); EM_mix_1 = rp(2)
   else
      call ErrorHandler('initScfData','Effective Medium Mixing Parameters are not found')
   endif
   rstatus = getKeyValue(tbl_id,'Effective Medium T-matrix Tol (> 0)',EM_tol)
   rstatus = getKeyValue(tbl_id,'Effective Medium Mixing eSwitch Value',EM_switch)
!
   end subroutine initScfData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endScfData()
!  ===================================================================
   implicit none
   if (NumKMeshs > 0) then
      deallocate(Kdiv)
   endif
   end subroutine endScfData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isRealAxisSS() result(rc)
!  ===================================================================
   implicit none
!
   logical rc
!
   if (RealAxisSS == 1) then
      rc = .true.
   else 
      rc = .false.
   endif
!
   end function isRealAxisSS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSingleSiteEContour() result(rc)
!  ===================================================================
   implicit none
!
   logical rc
!
   if (SingleSiteEContour == 1) then
      rc = .true.
   else 
      rc = .false.
   endif
!
   end function isSingleSiteEContour
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isFrozenCore(iter,fcf_name,fcf_exist) result(fc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: iter
!
   character (len=*), intent(inout), optional :: fcf_name
!
   logical, intent(out), optional :: fcf_exist
   logical :: fc
!
   if (present(fcf_name) .and. FrozenCoreFile_specified) then 
      fcf_name = FrozenCoreFile_name
   endif
!
   if (present(fcf_exist)) then
      if (present(fcf_name)) then
         inquire(file=fcf_name,exist=fcf_exist)
      else if (FrozenCoreFile_specified) then
         inquire(file=FrozenCoreFile_name,exist=fcf_exist)
      else
         call ErrorHandler('isFrozenCore','Frozen core file name is not specified')
      endif
   endif
!
   if (present(iter)) then
      if (iter > nFrozenCore) then
         fc = FrozenCore
      else
         fc = .false.
      endif
   else
      fc = FrozenCore
   endif
!
   end function isFrozenCore
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isNonRelativisticCore() result(rc)
!  ===================================================================
   implicit none
!
   logical rc
!
   if (nrelc == 0) then
      rc = .true.
   else if (nrelc == 1) then
      rc = .false.
   else
      call ErrorHandler('isNonRelativisticCore','invalid nrelc',nrelc)
   endif
!
   end function isNonRelativisticCore
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isNonRelativisticValence() result(rv)
!  ===================================================================
   implicit none
!
   logical rv
!
   if (nrelv == 0) then
      rv = .true.
   else if (nrelv == 1 .or. nrelv == 2) then
      rv = .false.
   else
      call ErrorHandler('isNonRelativisticValence','invalid nrelv',nrelv)
   endif
!
   end function isNonRelativisticValence
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isScalarRelativisticValence() result(rv)
!  ===================================================================
   implicit none
!
   logical rv
!
   if (nrelv == 0 .or. nrelv == 2) then
      rv = .false.
   else if (nrelv == 1) then
      rv = .true.
   else
      call ErrorHandler('isScalarRelativisticValence','invalid nrelv',nrelv)
   endif
!
   end function isScalarRelativisticValence
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isRelativisticValence() result(rv)
!  ===================================================================
   implicit none
!
   logical rv
!
   if (nrelv == 0 .or. nrelv == 1) then
      rv = .false.
   else if (nrelv == 2) then
      rv = .true.
   else
      call ErrorHandler('isRelativisticValence','invalid nrelv',nrelv)
   endif
!
   end function isRelativisticValence
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printScfData(fu_in)
!  ===================================================================
   use PublicParamDefinitionsModule, only : MuffinTin, ASA, MuffinTinASA, &
                                            MuffinTinTest
   implicit none
!
   integer (kind=IntKind), optional :: fu_in
   integer (kind=IntKind) :: fu
!
   if (present(fu_in)) then
      fu = fu_in
   else
      fu = 6
      write(fu,'(/,80(''-''))')
      write(fu,'(/,24x,a)')'********************************'
      write(fu,'( 24x,a )')'*   Output from printScfData   *'
      write(fu,'(24x,a,/)')'********************************'
      write(fu,'(/,80(''=''))')
   endif
   if ( .not.isSSIrregularSolOn() ) then
      write(fu,'(a)')'# Green function    : Integration along real energy axis for single scattering'
      write(fu,'(a)')'#                     Integration along complex energy contour for multiple scattering'
   endif
   if ( pot_type==MuffinTin .or. pot_type==ASA .or. pot_type==MuffinTinASA .or.  &
        pot_type==MuffinTinTest ) then
      write(fu,'(a)')'# Single Site Solver: Spherical'
   else
      write(fu,'(a)')'# Single Site Solver: Full-Potential'
      if ( sss_method == -1 ) then
         write(fu,'(a)')'#                     Regular solution is calculated without truncating potential'
         write(fu,'(a)')'#                     S and C are calculated using surface integration method'
      else if ( sss_method == 0 ) then
         write(fu,'(a)')'#                     Regular solution is calculated without truncating potential'
         write(fu,'(a)')'#                     S and C are calculated using volume integration and untruncated potential'
      else if ( sss_method == 1 ) then
         write(fu,'(a)')'#                     Regular solution is calculated without truncating potential'
         write(fu,'(a)')'#                     S and C are calculated using volume integration and truncated potential'
      else
         write(fu,'(a)')'#                     Regular solution is calculated using truncated potential'
         write(fu,'(a)')'#                     S and C are calculated using bounding sphere method'
      endif
   endif
   if ( scf_method == 0) then
      write(fu,'(a)')'# MST Method        : ScreenKKR_LSMS'
   else if ( scf_method == 1) then
      write(fu,'(a)')'# MST Method        : LSMS'
   else if ( scf_method == 2) then
      write(fu,'(a)')'# MST Method        : KKR'
   else if ( scf_method == 3) then
      write(fu,'(a)')'# MST Method        : KKRCPA'
   endif
   if ( nspin==1 ) then
      write(fu,'(a,i3)')'# Spin Parameter    : Non-magnetic -',nspin
   else if ( nspin==2) then
      write(fu,'(a,i3)')'# Spin Parameter    : Collinear -',nspin
   else if ( nspin>2) then
      write(fu,'(a,i3)')'# Spin Parameter    : Non-Collinear -',nspin
   endif
   if (fu==6) then
      write(fu,'(80(''=''))')
   endif
!
   end subroutine printScfData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isScreenKKR() result(md)
!  ===================================================================
   implicit none
   logical :: md
!
   if (scf_method == ScreenKKR) then
       md = .true.
   else
       md = .false.
   endif
   end function isScreenKKR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isKKR() result(md)
!  ===================================================================
   implicit none
   logical :: md
!
   if (scf_method == KKR) then
       md = .true.
   else
       md = .false.
   endif
   end function isKKR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLSMS() result(md)
!  ===================================================================
   implicit none
   logical :: md
!
   if (scf_method == LSMS) then
       md = .true.
   else
       md = .false.
   endif
   end function isLSMS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isKKRCPA() result(md)
!  ===================================================================
   implicit none
   logical :: md
!
   if (scf_method == KKRCPA) then
       md = .true.
   else
       md = .false.
   endif
   end function isKKRCPA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isEmbeddedCluster() result(md)
!  ===================================================================
   implicit none
   logical :: md
!
   if (scf_method == EmbeddedCluster) then
       md = .true.
   else
       md = .false.
   endif
   end function isEmbeddedCluster
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isScreenKKR_LSMS() result(md)
!  ===================================================================
   implicit none
   logical :: md
!
   if (scf_method_save == ScreenKKR_LSMS) then
       md = .true.
   else
       md = .false.
   endif
   end function isScreenKKR_LSMS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSingleSite() result(md)
!  ===================================================================
   implicit none
   logical :: md
!
   if (scf_method == SingleSite) then
       md = .true.
   else
       md = .false.
   endif
   end function isSingleSite
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSCFMethod(str_method)
!  ===================================================================
   implicit none
   character(len=*), intent(in) :: str_method
!
   if (str_method == "SingleSite") then
       scf_method = -2
   else if (str_method == "ScreenKKR") then
       scf_method = 0
   else if (str_method == "ScreenKKR_LSMS") then
       scf_method = -1
   else if (str_method == "LSMS") then
       scf_method = 1
   else if (str_method == "KKR") then
       scf_method = 2
   else if (str_method == "KKRCPA") then
       scf_method = 3
   else
       call ErrorHandler('serSCFMethod','unknown SCF method',scf_method)
   endif
   end subroutine setSCFMethod
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isReadEmesh() result(pt)
!  ===================================================================
   implicit none
   logical :: pt
   if (read_emesh == ReadEmesh) then
      pt = .true.
   else
      pt = .false.
   endif
   end function isReadEmesh
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isReadKmesh() result(pt)
!  ===================================================================
   implicit none
   logical :: pt
   if (read_kmesh == ReadKmesh) then
      pt = .true.
   else
      pt = .false.
   endif
   end function isReadKmesh
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKmeshFileName() result(fn)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: fn
   fn = trim(adjustl(inputpath))//'kmeshs.inp'
   end function getKmeshFileName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEmeshFileName() result(fn)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: fn
   fn = trim(adjustl(inputpath))//'emeshs.inp'
   end function getEmeshFileName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSingleSiteSolverType() result(sst)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: sst
!
   sst = ss_solver_type
!
   end function getSingleSiteSolverType
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSSIrregularSolOn() result (isIrrSolOn)
!  ===================================================================
   implicit none
!
   logical :: isIrrSolOn
!
   if (ss_irrsol == 0 ) then
      isIrrSolOn = .false.
   else
      isIrrSolOn = .true.
   endif
!
   end function isSSIrregularSolOn
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isFittedChargeDen() result (isChargeFitted)
!  ===================================================================
   implicit none
!
   logical :: isChargeFitted
!
   if (ss_irrsol == 0 .or. ss_irrsol == 2 ) then
      isChargeFitted = .false.
   else
      isChargeFitted = .true.
   endif
!
   end function isFittedChargeDen
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isEfIterateOn(pin) result (isEfIterOn)
!  ===================================================================
   implicit none
!
   logical, intent(out), optional :: pin
   logical :: isEfIterOn
!
   if (present(pin)) then
      if (iterateEf == 2) then ! pin the Ef to a fixed value
         pin = .true.
      else
         pin = .false.
      endif
   endif
!
   if (iterateEf == 1 ) then
      isEfIterOn = .false.
   else
      isEfIterOn = .true.
   endif
!
   end function isEfIterateOn
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSSPhaseShiftOn() result (isPShiftOn)
!  ===================================================================
   implicit none
!
   logical :: isPShiftOn
!
   if (ss_phaseshift == 0 ) then
      isPShiftOn = .false.
   else
      isPShiftOn = .true.
   endif
!
   end function isSSPhaseShiftOn
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSingleSiteSolverType(sst)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: sst
!
   ss_solver_type = sst
!
   end subroutine setSingleSiteSolverType
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSSSIntegrationMethod()               result(ssintegral)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: ssintegral
!
   ssintegral = ss_integral_method
!
   end function getSSSIntegrationMethod
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSSSolutionsMethod()               result(sssolution)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: sssolution
!
   sssolution = ss_solution_method
!
   end function getSSSolutionsMethod
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSingleSiteSolverMethod()            result(sssm)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: sssm
!
   sssm = sss_method
!
   end function getSingleSiteSolverMethod
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxPotSolver()                      result(lmax_pot)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: lmax_pot
!
   lmax_pot = lmax_pot_solver
!
   end function getLmaxPotSolver
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxSolution()                      result(lmax_sol)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: lmax_sol
!
   lmax_sol = lmax_sol_cutoff
!
   end function getLmaxSolution
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDOSrunid()                      result(Drid)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: Drid
!
   Drid = DOSrunid
!
   end function getDOSrunid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotentialTypeParam() result(pt)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: pt
!
   pt = pot_type
!
   end function getPotentialTypeParam
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCantedMomentTorqueFactor() result(c)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: c
!
   c = ctq
!
   end function getCantedMomentTorqueFactor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPoleSearchStep() result(step)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: step
!
   step = pole_step
!
   end function getPoleSearchStep

!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isExchangeParamNeeded() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if (j_ij == 'y' .or. j_ij == 'Y') then
      c = .true.
   else
      c = .false.
   endif
!
   end function isExchangeParamNeeded
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSimpleMixing() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( mixing_type == 0 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isSimpleMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isDGAMixing() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( mixing_type == 1 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isDGAMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isBroydenMixing() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( mixing_type == 2 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isBroydenMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isPotentialMixing() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( mixing_quantity == 1 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isPotentialMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isChargeMixing() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( mixing_quantity == 0 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isChargeMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isCheckSuperLUSolver()  result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( KSpaceSolver == 1 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isCheckSuperLUSolver
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLloyd()  result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( .not.isLSMS() .and. Lloyd /=0 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isLloyd
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLloydMode()  result(mlloyd)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind) :: mlloyd
!
   mlloyd = LloydMode
!
 end function getLloydMode
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isDirectKSpaceSolver()  result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( KSpaceSolver == 2 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isDirectKSpaceSolver
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isInterstitialElectronPolarized() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( i_vdif == 1 ) then
      c = .false.
   else
      c = .true.
   endif
!
   end function isInterstitialElectronPolarized
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLdaCorrectionNeeded() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( LdaCorrectionType == 0 ) then
      c = .false.
   else if (n_spin_pola ==2) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isLdaCorrectionNeeded
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isChargeSymm() result(c)
!  ===================================================================
   implicit none
!
   logical :: c
!
   if ( ChargeSymmetry == 1 ) then
      c = .true.
   else
      c = .false.
   endif
!
   end function isChargeSymm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAdaptiveIntegrationMethod() result(m)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: m
!
   m = RealEIntMethod
!
   end function getAdaptiveIntegrationMethod
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getUJfile() result(fn)
!  ===================================================================
   use InputModule, only : getKeyValue
   implicit none
!
   character (len=MaxLenFileName) :: fn
!
   integer (kind=IntKind) :: rstatus
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(TableID,'LDA+U Parameter File Name',UJfile)
!  -------------------------------------------------------------------
!
   UJfile = adjustl(UJfile)
   if ( nocaseCompare(UJfile,'None') ) then
      fn = 'None'
   else
      fn = trim(adjustl(inputpath))//UJfile
   endif
!
   end function getUJfile
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine retrieveEffectiveMediumParams(mix_type, max_iter,       &
                                            alpha_0, alpha_1, eSwitch, tol)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out) :: mix_type
   integer (kind=IntKind), intent(out) :: max_iter
!
   real (kind=RealKind), intent(out) ::  alpha_0
   real (kind=RealKind), intent(out) ::  alpha_1
   real (kind=RealKind), intent(out) ::  tol
   real (kind=RealKind), intent(out) ::  eSwitch
!
   mix_type = EM_mix_type
   max_iter = EM_max_iter
   alpha_0  = EM_mix_0
   alpha_1  = EM_mix_1
   tol      = EM_tol
   eSwitch  = EM_switch
!
   end subroutine retrieveEffectiveMediumParams
!  ===================================================================
end module ScfDataModule
