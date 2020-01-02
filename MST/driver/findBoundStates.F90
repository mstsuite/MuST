program findBoundStates
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use DataServiceCenterModule, only : isDataStorageExisting, &
                                       getDataStorage, RealMark
!
   use TimerModule, only : getTime
!
   use SystemModule, only : getNumAtoms, getAtomPosition
!
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop, EvBottom, isNonRelativisticCore
   use ScfDataModule, only : ngaussr, ngaussq
   use ScfDataModule, only : pole_step
!
   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies
!
   use MathParamModule, only : ZERO, TEN2m8, TEN2m6, THIRD, HALF, ONE, TWO, PI, FIVE, PI2, PI4, &  
                               CZERO, CONE, SQRTm1, Y0
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use MadelungModule, only : initMadelung, endMadelung
!
   use SpinRotationModule, only : initSpinRotation, endSpinRotation
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,   &
                                    calSymmetryFlags
!
   use CoreStatesModule, only : initCoreStates, calCoreStates, endCoreStates
   use CoreStatesModule, only : readCoreStates, printCoreStates, printCoreDensity
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : getVolume, getInscrSphVolume
!
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : getVolumeIntegration
!
   use OutputModule, only : getStandardOutputLevel
!
   use Atom2ProcModule, only : getLocalNumAtoms, getGlobalIndex
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential
!
   use SurfElementsModule, only : initSurfElements, endSurfElements,  &
                                  genGaussPoints, genGaussSphHarmonics
!
   use SSSolverModule, only : initSSSolver, endSSSolver
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : computeDOS, getDOS
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use RadialGridModule, only : getMaxNumRmesh, getGrid
!
   use WriteFunctionModule, only : writeFunction
!
   use IntegerFactorsModule, only : initIntegerFactors, endIntegerFactors, &
                                    mofj, lofj
   use MPPModule, only : syncAllPEs
!
   implicit   none
!
   logical :: non_zero
!
   character (len=6) :: state_string, jl_string
   character (len=8) :: app_string
   character (len=50) :: file_name
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: NumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: lmax_max, lmax_pot_max, lmax_kkr_max, jmax_rho, &
                             lmax_phi_max, lmax_rho_max, lmax_step_max, lmax_green_max
   integer (kind=IntKind) :: kmax_kkr, kmax_kkr_max, jmax_rho_max, max_bound_states
   integer (kind=IntKind) :: id, ig, is, i, ib, jl, l, m
   integer (kind=IntKind) :: RelativisticFlag
   integer (kind=IntKind) :: node_print_level
!
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_green(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: NumBoundPoles(:,:), NumResPoles(:,:)
   integer (kind=IntKind), allocatable :: NumBPDegens(:,:,:), NumRPDegens(:,:,:)
   integer (kind=IntKind) :: eGID, NumPEsInEGroup, MyPEInEGroup
   integer (kind=IntKind) :: kl, klp, kl1, kl2, kl3, ir, nr, occ, ie
!
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), allocatable :: AtomPosition(:,:), rho0(:)
   real (kind=RealKind) :: t0, t1, t2, ebot, etop, rfac
!
   real (kind=RealKind), allocatable :: BoundPoles(:,:,:), ResPoles(:,:,:), ResWidth(:,:,:)
   real (kind=RealKind), parameter :: ZeroIntHalfWidth = TEN2m6
!
   complex (kind=CmplxKind) :: e 
   complex (kind=CmplxKind), allocatable :: BmatResidual(:,:,:,:), RmatResidual(:,:,:,:)
   complex (kind=CmplxKind), allocatable :: Bdensity(:,:,:,:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
!
   type (GridStruct), pointer :: Grid
!
   interface
      subroutine calQuadraticPoles(LocalNumAtoms,lkkr,lphi,                        &
                                   n_spin_pola,eb,et,del,ldp,                      &
                                   NumBoundPoles,NumBPDegens,BoundPoles,Bmat,      &
                                   NumResPoles,NumRPDegens,ResPoles,ResWidth,Rmat, &
                                   CheckPoles,PanelOnZero)
      use KindParamModule, only : IntKind, RealKind, CmplxKind
      implicit none
      integer (kind=IntKind), intent(in) :: LocalNumAtoms
      integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
      integer (kind=IntKind), intent(in) :: lphi(1:LocalNumAtoms)
      integer (kind=IntKind), intent(in) :: n_spin_pola,ldp
      integer (kind=IntKind), intent(out) :: NumBoundPoles(n_spin_pola,LocalNumAtoms)
      integer (kind=IntKind), intent(out) :: NumBPDegens(ldp,n_spin_pola,LocalNumAtoms)
      integer (kind=IntKind), intent(out) :: NumResPoles(n_spin_pola,LocalNumAtoms)
      integer (kind=IntKind), intent(out) :: NumRPDegens(ldp,n_spin_pola,LocalNumAtoms)
!
      logical, optional, intent(in) :: CheckPoles
      logical, optional, intent(in) :: PanelOnZero
!
      real (kind=RealKind), intent(in) :: eb, et, del
      real (kind=RealKind), intent(out) :: BoundPoles(ldp,n_spin_pola,LocalNumAtoms)
      real (kind=RealKind), intent(out) :: ResPoles(ldp,n_spin_pola,LocalNumAtoms)
      real (kind=RealKind), intent(out) :: ResWidth(ldp,n_spin_pola,LocalNumAtoms)
!
      complex (kind=CmplxKind), intent(out) :: Bmat(ldp*ldp,ldp,n_spin_pola,LocalNumAtoms)
      complex (kind=CmplxKind), intent(out) :: Rmat(ldp*ldp,ldp,n_spin_pola,LocalNumAtoms)
!
      end subroutine calQuadraticPoles
   end interface
!
   interface
      subroutine calBoundStateDensity(LocalNumAtoms,lmax_kkr,lmax_rho, &
                                      n_spin_pola,max_nr,kmax_kkr_max, &
                                      jmax_rho_max,max_bound_states,   &
                                      NumBoundPoles, BoundPoles, Bmat, &
                                      Bdensity)
      use KindParamModule, only : IntKind, RealKind, CmplxKind
!
      implicit none
!
      integer (kind=IntKind), intent(in) :: LocalNumAtoms
      integer (kind=IntKind), intent(in) :: lmax_kkr(1:LocalNumAtoms)
      integer (kind=IntKind), intent(in) :: lmax_rho(1:LocalNumAtoms)
      integer (kind=IntKind), intent(in) :: n_spin_pola,kmax_kkr_max,jmax_rho_max,max_nr,max_bound_states
      integer (kind=IntKind), intent(in) :: NumBoundPoles(n_spin_pola,LocalNumAtoms)
!
      real (kind=RealKind), intent(in) :: BoundPoles(kmax_kkr_max,n_spin_pola,LocalNumAtoms)
!
      complex (kind=CmplxKind), intent(in) :: Bmat(kmax_kkr_max**2,kmax_kkr_max,n_spin_pola,LocalNumAtoms)
      complex (kind=CmplxKind), intent(out) :: Bdensity(max_nr,jmax_rho_max,max_bound_states,n_spin_pola,LocalNumAtoms)
!
      end subroutine calBoundStateDensity
   end interface
!
!  -------------------------------------------------------------------
   call startProcess()
!  -------------------------------------------------------------------
   t0 = getTime()
!  -------------------------------------------------------------------
   call initializeModules()
!  -------------------------------------------------------------------
   call initCoreStates(LocalNumAtoms,EvBottom,n_spin_pola,            &
                       isNonRelativisticCore(),istop,1)
!  -------------------------------------------------------------------
   call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, &
                     lmax_kkr, lmax_phi, lmax_pot, lmax_step, lmax_green,  &
                     n_spin_pola, n_spin_cant, RelativisticFlag, istop, atom_print_level)
!  -------------------------------------------------------------------
   call syncAllPEs()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Check the symmetry of single site density at an energy > 0
!  ===================================================================
   do id = 1, LocalNumAtoms
      Grid => getGrid(id)
      do is = 1,n_spin_pola
!        -------------------------------------------------------------
         call solveSingleScattering(is, id, (0.6d0,0.0d0), CZERO)
         call computeDOS()
         dos_r_jl => getDOS()
!        -------------------------------------------------------------
         jl = 0
         do l = 0, lmax_green(id)
            do m = 0, l
               jl = jl + 1
               non_zero = .false.
               LOOP_ir_2: do ir = 1, Grid%jend
                  if (abs(dos_r_jl(ir,jl)) > TEN2m6) then
                     non_zero = .true.
                     exit LOOP_ir_2
                  endif
               enddo LOOP_ir_2
               if (non_zero .and. node_print_level >= 0) then
                  write(6,'(a,i2,a,i2,a)')'Non-zero single site density component at (l,m) channel: (',l,',',m,')'
               endif
            enddo
         enddo
      enddo
   enddo
!
   write(6,'(/,a,f12.5,a)') 'Setup time: ',getTime() - t0,' sec.'
!
!  ===================================================================
!  read core states and density data of local atoms from file
!  ===================================================================
   call readCoreStates()
!  -------------------------------------------------------------------
   call calCoreStates()
!  -------------------------------------------------------------------
   do id = 1,LocalNumAtoms
      if (atom_print_level(id) >= 0) then
!        -------------------------------------------------------------
         call printCoreStates(id)
         call printCoreDensity(id)
!        -------------------------------------------------------------
      endif
   enddo
!
!  -------------------------------------------------------------------
   allocate( NumBoundPoles(n_spin_pola,LocalNumAtoms) )
   allocate( NumBPDegens(kmax_kkr_max,n_spin_pola,LocalNumAtoms) )
   allocate( NumResPoles(n_spin_pola,LocalNumAtoms) )
   allocate( NumRPDegens(kmax_kkr_max,n_spin_pola,LocalNumAtoms) )
   allocate( BoundPoles(kmax_kkr_max,n_spin_pola,LocalNumAtoms) )
   allocate( ResPoles(kmax_kkr_max,n_spin_pola,LocalNumAtoms) )
   allocate( ResWidth(kmax_kkr_max,n_spin_pola,LocalNumAtoms) )
   allocate( BmatResidual(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max,n_spin_pola,LocalNumAtoms) )
   allocate( RmatResidual(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max,n_spin_pola,LocalNumAtoms) )
!  -------------------------------------------------------------------
!
   t1 = getTime()
   ebot = -10.0d0
   etop = -TEN2m6
!  ------------------------------------------------------------------
   call calQuadraticPoles(LocalNumAtoms,lmax_kkr,lmax_phi,            &
                          n_spin_pola,ebot,etop,pole_step,            &
                          kmax_kkr_max,                               &
                          NumBoundPoles,NumBPDegens,BoundPoles,BmatResidual,      &
                          NumResPoles,NumRPDegens,ResPoles,ResWidth,RmatResidual, &
                          CheckPoles =.false.)
!  ------------------------------------------------------------------
   if (node_print_level >= 0) then
      write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
      do id = 1,LocalNumAtoms
         do is = 1,n_spin_pola
            write(6,'(a,i3,a,i3)')'Spin index: ',is,',  Atom index: ',id
            write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of bound states found within (',  &
                                             ebot,', ',ZERO,'): ',NumBoundPoles(is,id)
            do i = 1, NumBoundPoles(is,id)
               write(6,'(a,i2,5x,a,f20.12)')'Degeneracy = ',NumBPDegens(i,is,id), &
                                            ', Bound state energy = ',BoundPoles(i,is,id)
            enddo
            write(6,'(/,a,f13.8,a,f13.8,a,i5)')'The Number of resonance states found within (',  &
                                               ebot,', ',etop,'): ',NumResPoles(is,id)
            do i = 1, NumResPoles(is,id)
               write(6,'(a,i2,5x,a,f20.12,a,f20.12)')'Degeneracy = ',NumRPDegens(i,is,id), &
                                                     ', Resonance state energy = ',ResPoles(i,is,id), &
                                                     ', Width = ',ResWidth(i,is,id)
            enddo
         enddo
      enddo
      write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
   endif
!
   max_bound_states = 0 
   do id = 1,LocalNumAtoms
      do is = 1,n_spin_pola
         max_bound_states = max(max_bound_states,NumBoundPoles(is,id))
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call initIntegerFactors(lmax_kkr_max)
!  -------------------------------------------------------------------
!
   nr = getMaxNumRmesh()
   allocate( Bdensity(nr,jmax_rho_max,max_bound_states,n_spin_pola,LocalNumAtoms) )
   allocate( rho0(nr) )
   t2 = getTime()
!  ------------------------------------------------------------------
   call calBoundStateDensity(LocalNumAtoms,lmax_kkr,lmax_rho,        &
                             n_spin_pola,nr,kmax_kkr_max,            &
                             jmax_rho_max,max_bound_states,          &
                             NumBoundPoles, BoundPoles, BmatResidual,&
                             Bdensity)
!  ------------------------------------------------------------------
   if (node_print_level >= 0) then
      write(6,'(/,a,f12.5,a,/)') 'calBoundStateDensity time: ',getTime() - t1,' sec.'
      do id = 1, LocalNumAtoms
         Grid => getGrid(id)
         nr = Grid%jend
         r_mesh => Grid%r_mesh
         jmax_rho = (lmax_rho(id)+1)*(lmax_rho(id))/2
         do is = 1, n_spin_pola
            do ie = 1, NumBoundPoles(is,id)
               if (id < 10 .and. ie < 10) then
                  write(app_string,'(a,i1,a,i1,a,i1)')'a',id,'s',is,'e',ie
               else if (id >= 10 .and. ie < 10) then
                  write(app_string,'(a,i2,a,i1,a,i1)')'a',id,'s',is,'e',ie
               else if (id < 10 .and. ie >= 10) then
                  write(app_string,'(a,i1,a,i1,a,i2)')'a',id,'s',is,'e',ie
               else
                  write(app_string,'(a,i2,a,i1,a,i2)')'a',id,'s',is,'e',ie
               endif
               occ = (3-n_spin_pola)*NumBPDegens(ie,is,id)
               if (occ < 10) then
                  write(state_string,'(a,i1)')'Occ',occ
               else
                  write(state_string,'(a,i2)')'Occ',occ
               endif
!              =======================================================
!              Note: BmatResidual already contains a summation of the 
!                    contributions from the degenerate poles.
!              =======================================================
               rfac = Y0/real(NumBPDegens(ie,is,id),kind=RealKind)
               do ir = 1, nr
                  rho0(ir) = real(Bdensity(ir,1,ie,is,id),kind=RealKind)*rfac
               enddo
               file_name = 'BndState_'//trim(state_string)//app_string
!              -------------------------------------------------------
               call writeFunction(file_name,nr,r_mesh,rho0)
!              -------------------------------------------------------
               do jl = 2, jmax_rho
                  non_zero = .false.
                  LOOP_ir: do ir = 1, nr
                     if (abs(Bdensity(ir,jl,ie,is,id)) > TEN2m6) then
                        non_zero = .true.
                        exit LOOP_ir
                     endif
                  enddo LOOP_ir
                  if (non_zero) then
                     if (lofj(jl) < 10) then
                        write(jl_string,'(a,i1,a,i1)')'l',lofj(jl),'m',mofj(jl)
                     else if (mofj(jl) < 10) then
                        write(jl_string,'(a,i2,a,i1)')'l',lofj(jl),'m',mofj(jl)
                     else
                        write(jl_string,'(a,i2,a,i2)')'l',lofj(jl),'m',mofj(jl)
                     endif
                     file_name = 'BndState_'//trim(state_string)//trim(app_string)//jl_string
!                    -------------------------------------------------
                     call writeFunction(file_name,nr,r_mesh,Bdensity(:,jl,ie,is,id))
!                    -------------------------------------------------
                  endif
               enddo
            enddo
         enddo
      enddo
   endif
!
!  -------------------------------------------------------------------
   call endIntegerFactors()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  clean up the allocated spaces.
!  ===================================================================
   deallocate( NumBoundPoles, BoundPoles, NumResPoles, ResPoles, ResWidth )
   deallocate( NumBPDegens, NumRPDegens, BmatResidual, RmatResidual, Bdensity, rho0 )
!
   if (node_print_level >= 0) then
      write(6,*) " Run Ending! "
      call FlushFile(6)
   endif
!
!  ===================================================================
!  finish the process
!  -------------------------------------------------------------------
   call endSSSolver()
   call endCoreStates()
   call cleanModules()
   call finishProcess()
!  -------------------------------------------------------------------
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initializeModules()
!  ===================================================================
   use PublicTypeDefinitionsModule, only : GridStruct
   use ScfDataModule, only : ErBottom, ErTop, EiBottom, EiTop, Temperature
   use ScfDataModule, only : isNonRelativisticValence
   use ScfDataModule, only : isScalarRelativisticValence
   use ScfDataModule, only : getSingleSiteSolverMethod
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax, getStepFuncLmax
   use AtomModule, only : setTruncPotLmax, getTruncPotLmax, setPotLmax, getLocalEvecOld
   use PotentialTypeModule, only : isTestPotential, isFullPotential
   use TestPotentialModule, only : initTestPotential, endTestPotential,      &
                                   readTestPotential, getTestPotential,      &
                                   getTestPotEf, getTestV0
   use RadialGridModule, only : getGrid
   use PotentialModule, only : getPotEf, setV0, setPotential, isPotComponentZero, getSphPotr
   use PotentialModule, only : getPotComponentFlag, getTruncatedPotComponentFlag
   use PotentialModule, only : getPotential, getTruncatedPotential
   use ContourModule, only : getEPoint, setupContour
   use ProcMappingModule, only : getNumEsOnMyProc, getEnergyIndex
!
   implicit none
!
   integer (kind=IntKind) :: i, id, ig, is, l, m, jl, ie_glb, ie_loc
   integer (kind=IntKind) :: jmax_pot, jmax_rho, NumEsOnMyProc
   integer (kind=IntKind) :: ndivin, ndivout, nmult, lmax_max
!
   real (kind=RealKind) :: rmt, rend, rws, rinsc, v0, Efermi
   real (kind=RealKind), pointer :: potr_0(:)
   real (kind=RealKind), allocatable :: evec(:,:)
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!
   complex (kind=CmplxKind), pointer :: EPoint(:)
   complex (kind=CmplxKind), pointer :: pot_l(:), pot_jl(:,:)
   complex (kind=CmplxKind) :: energy
!
   type (GridStruct), pointer :: Grid
!
   logical :: non_zero
!
   integer (kind=IntKind), pointer :: pflag(:)
!
   NumAtoms = getNumAtoms()
   LocalNumAtoms=getLocalNumAtoms()
!
   allocate(AtomPosition(1:3,1:NumAtoms))
   allocate(GlobalIndex(LocalNumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
   enddo
!
   node_print_level = getStandardOutputLevel()
!
   allocate(atom_print_level(1:LocalNumAtoms))
   do i=1,LocalNumAtoms
      atom_print_level(i) = getStandardOutputLevel(i)
   enddo
!
!  ===================================================================
!  If it is not the full-potential calculation, set the Lmax for potential to 0             
!  ===================================================================
   if (.not.isFullPotential()) then
      do i=1,LocalNumAtoms
!        -------------------------------------------------------------
         call setPotLmax(i,0)
         call setTruncPotLmax(i,0)
!        -------------------------------------------------------------
      enddo
   endif
!
!  ===================================================================
!  setup the Lmax values.             
!  ===================================================================
   allocate(lmax_pot(LocalNumAtoms), lmax_rho(LocalNumAtoms))
   allocate(lmax_kkr(LocalNumAtoms), lmax_phi(LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms))
   allocate(lmax_green(LocalNumAtoms))
   lmax_kkr_max = 0
   lmax_phi_max = 0
   lmax_rho_max = 0
   lmax_pot_max = 0
   lmax_step_max = 0
   lmax_green_max = 0
   lmax_max = 0
   do i=1,LocalNumAtoms
      lmax_kkr(i) = getKKRLmax(i)
      lmax_phi(i) = getPhiLmax(i)
      lmax_rho(i) = getRhoLmax(i)
      lmax_pot(i) = getPotLmax(i)
      lmax_kkr_max = max(lmax_kkr_max,lmax_kkr(i))
      lmax_phi_max = max(lmax_phi_max,lmax_phi(i))
      lmax_rho_max = max(lmax_rho_max,lmax_rho(i))
      lmax_pot_max = max(lmax_pot_max,lmax_pot(i))
      lmax_step_max = max(lmax_step_max,lmax_step(i))
      if (isFullPotential()) then
         lmax_step(i)  = max(getStepFuncLmax(i),lmax_phi(i)*2)
         lmax_max = max( lmax_max, 2*lmax_step(i), 2*lmax_pot(i),     &
                         2*lmax_rho(i), lmax_phi(i) )
      else
         lmax_step(i)  = getStepFuncLmax(i)
         lmax_max = max( lmax_max, lmax_step(i), lmax_rho(i), lmax_kkr(i), lmax_phi(i) )
      endif
!     lmax_green(i) = min(2*lmax_phi(i),lmax_rho(i))
      lmax_green(i) = lmax_rho(i)
      lmax_green_max = max(lmax_green(i),lmax_green_max)
   enddo
!
   jmax_rho_max = (lmax_rho_max+1)*(lmax_rho_max+2)/2
   kmax_kkr_max = (lmax_kkr_max+1)**2
!
!  -------------------------------------------------------------------
!  call initSphericalHarmonics(4*lmax_max)
   call initSphericalHarmonics(2*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
!  call initGauntFactors(2*lmax_max,istop,0)
   call initGauntFactors(lmax_max,istop,0)
!  -------------------------------------------------------------------
!
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i)
      GlobalIndex(i)=ig
   enddo
!
!  -------------------------------------------------------------------
   call setupRadGridAndCell(LocalNumAtoms,lmax_max)
!  -------------------------------------------------------------------
!
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('testSineMatrixPole','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  initialize potential module
!  -------------------------------------------------------------------
   call initMadelung(LocalNumAtoms,NumAtoms,GlobalIndex,              &
                     lmax_rho_max,lmax_pot_max,bravais,AtomPosition,0)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call initSystemSymmetry( NumAtoms, LocalNumAtoms, lmax_pot, lmax_step, atom_print_level )
!  -------------------------------------------------------------------
   call calSymmetryFlags()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot,lmax_step,               &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
!
   if (isTestPotential()) then
!     ----------------------------------------------------------------
      call initTestPotential(LocalNumAtoms,n_spin_pola)
!     ----------------------------------------------------------------
      call readTestPotential(lmax_pot)
!     ----------------------------------------------------------------
      do id = 1, LocalNumAtoms
         jmax_pot = (lmax_pot(id)+1)*(lmax_pot(id)+2)/2
         do is = 1,n_spin_pola
            do jl=1,jmax_pot
!              -------------------------------------------------------
               pot_l => getTestPotential(id,is,jl)
!              -------------------------------------------------------
               call setPotential(id,1,is,jl,pot_l,1)
!              -------------------------------------------------------
            enddo
         enddo
      enddo
      do is = 1,n_spin_pola
         v0 = getTestV0(is)
         call setV0(is,v0)
      enddo
      Efermi=getTestPotEf()
!     ----------------------------------------------------------------
      call endTestPotential()
!     ----------------------------------------------------------------
   else
!     ================================================================
!     read potential data
!     ----------------------------------------------------------------
      call readPotential()
!     ----------------------------------------------------------------
      do id = 1, LocalNumAtoms
         Grid => getGrid(id)
         ig = GlobalIndex(id)
         do is = 1,n_spin_pola
            jl = 0
            do l = 0, lmax_pot(id)
               do m = 0, l
                  jl = jl + 1
                  if (.not.isPotComponentZero(id,jl) .and. node_print_level >= 0) then
                     write(6,'(a,i2,a,i2,a)')'Non-zero potential component at (l,m) channel: (',l,',',m,')'
                  endif
               enddo
            enddo
            potr_0 =>  getSphPotr(id,1,is)
         enddo
         write(6,'(/)')
         pflag => getPotComponentFlag(id)
         jl = 0
         do l = 0, lmax_pot(id)
            do m = 0, l
               jl = jl + 1
               if (pflag(jl) /= 0 .and. node_print_level >= 0) then
                  write(6,'(a,i2,a,i2,a)')'pflag<>0 component at (l,m) channel: (',l,',',m,')'
               endif
            enddo
         enddo
         write(6,'(/)')
         pflag => getTruncatedPotComponentFlag(id)
         jl = 0
         do l = 0, lmax_pot(id)
            do m = 0, l
               jl = jl + 1
               if (pflag(jl) /= 0 .and. node_print_level >= 0) then
                  write(6,'(a,i2,a,i2,a)')'Truncated pflag<>0 component at (l,m) channel: (',l,',',m,')'
               endif
            enddo
         enddo
         do is = 1,n_spin_pola
            write(6,'(/,a,i2)')'spin = ',is
            pot_jl => getPotential(id,1,is)
            jl = 0
            do l = 0, lmax_pot(id)
               do m = 0, l
                  jl = jl + 1
                  non_zero = .false.
                  LOOP_ir: do ir = 1, Grid%jend
                     if (abs(pot_jl(ir,jl)) > TEN2m6) then
                        non_zero = .true.
                        exit LOOP_ir
                     endif
                  enddo LOOP_ir
                  if (non_zero .and. node_print_level >= 0) then
                     write(6,'(a,i2,a,i2,a)')'Recheck non-zero potential component at (l,m) channel: (',l,',',m,')'
                  endif
               enddo
            enddo
            write(6,'(/)')
            pot_jl => getTruncatedPotential(id,1,is)
            jl = 0
            do l = 0, getTruncPotLmax(id)
               do m = 0, l
                  jl = jl + 1
                  non_zero = .false.
                  LOOP_ir_1: do ir = 1, Grid%jend-Grid%jmt+1
                     if (abs(pot_jl(ir,jl)) > TEN2m6) then
                        non_zero = .true.
                        exit LOOP_ir_1
                     endif
                  enddo LOOP_ir_1
                  if (non_zero .and. node_print_level >= 0) then
                     write(6,'(a,i2,a,i2,a)')'Recheck non-zero truncated potential component at (l,m) channel: (',l,',',m,')'
                  endif
               enddo
            enddo
         enddo
      enddo
   endif
!
   if (ErTop > ZERO) then
      Efermi = ErTop
   else
      Efermi=getPotEf()
   endif
!
   allocate(evec(3,LocalNumAtoms))
   do i = 1,LocalNumAtoms
      evec(1:3,i) = getLocalEvecOld(i)
   enddo
!  -------------------------------------------------------------------
   call initSpinRotation(LocalNumAtoms,evec)
!  -------------------------------------------------------------------
   call setupContour( ErBottom, Efermi, EiBottom, EiTop )
!  -------------------------------------------------------------------
!
   EPoint => getEPoint()
   NumEsOnMyProc = getNumEsOnMyProc()
   write(6,'(/,a)')'=================================================='
   write(6,'(  a)')'Energy #                    Energy Value'
   write(6,'(  a)')'--------------------------------------------------'
   do ie_loc = 1, NumEsOnMyProc
      ie_glb = getEnergyIndex(ie_loc)
      energy = EPoint(ie_glb)
      write(6,'(i5,12x,2d16.8)')ie_glb,energy
   enddo
   write(6,'(  a)')'=================================================='
!
   if (getSingleSiteSolverMethod() == -1) then
      call initSurfElements('none',-1)
      call genGaussPoints()
      call genGaussSphHarmonics(lmax_max)
   endif
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize single scattering solver
!  ===================================================================
   if (isNonRelativisticValence()) then
      RelativisticFlag = 0
   else if (isScalarRelativisticValence()) then
      RelativisticFlag = 1
   else
      RelativisticFlag = 2
   endif
   
   end subroutine initializeModules
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cleanModules()
!  ===================================================================
   use RadialGridModule, only : endRadialGrid
   implicit none
!
   deallocate( atom_print_level )
   deallocate( AtomPosition,lmax_kkr,lmax_phi,lmax_rho,lmax_pot )
   deallocate( lmax_step, lmax_green, GlobalIndex )
!
   call endSurfElements()
   call endSpinRotation()
   call endPotential()
   call endStepFunction()
   call endRadialGrid()
   call endSystemSymmetry()
   call endMadelung()
   call endGauntFactors()
   call endSphericalHarmonics()
!
   end subroutine cleanModules
!  ===================================================================
end program findBoundStates
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calQuadraticPoles(LocalNumAtoms,lkkr,lphi,                                    &
                                n_spin_pola,eb,et,Delta,ldp,                                &
                                NumBoundPoles, NumBPDegens, BoundPoles, BmatResidual,       &
                                NumResPoles, NumRPDegens, ResPoles, ResWidth, RmatResidual, &
                                CheckPoles,PanelOnZero)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                                Ten2m8, FOURTH, CZERO, TEN2m7
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : MyPE
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEInGroup
   use GroupCommModule, only : GlobalSumInGroup, bcastMessageInGroup
!
   use SSSolverModule, only : solveSingleScattering, getJostMatrix
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix, &
                                     solveQuadraticEquation, getEigenValue,   &
                                     solveLinearEquation, getEigenVector,     &
                                     getEigenMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: lphi(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: n_spin_pola,ldp
   integer (kind=IntKind), intent(out) :: NumBoundPoles(n_spin_pola,LocalNumAtoms)
   integer (kind=IntKind), intent(out) :: NumBPDegens(ldp,n_spin_pola,LocalNumAtoms)
   integer (kind=IntKind), intent(out) :: NumResPoles(n_spin_pola,LocalNumAtoms)
   integer (kind=IntKind), intent(out) :: NumRPDegens(ldp,n_spin_pola,LocalNumAtoms)
!
   integer (kind=IntKind) :: id, ie, is, iw, kmax_kkr, NumWindows, info
   integer (kind=IntKind) :: i, j, l, m, kl, lp, mp, klp, nv, nb, nr, je, ip
   integer (kind=IntKind) :: eGID, NumPEsInEGroup, MyPEInEGroup
   integer (kind=IntKind) :: MyNumWindows, nbr(2), nb0, nr0, ib, ir
   integer (kind=IntKind) :: bpdeg(ldp), rpdeg(ldp), degens(ldp)
!
   logical, optional, intent(in) :: CheckPoles
   logical, optional, intent(in) :: PanelOnZero
   logical :: isZeroInterval = .false.
   logical :: chkpole = .false.
   logical, parameter :: isGeneral = .false.
!
   real (kind=RealKind), intent(in) :: eb, et, Delta
   real (kind=RealKind) :: WindowWidth
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0, err
   real (kind=RealKind), intent(out) :: BoundPoles(ldp,n_spin_pola,LocalNumAtoms)
   real (kind=RealKind), intent(out) :: ResPoles(ldp,n_spin_pola,LocalNumAtoms)
   real (kind=RealKind), intent(out) :: ResWidth(ldp,n_spin_pola,LocalNumAtoms)
   real (kind=RealKind) :: bpe(ldp), ebr(ldp), bpe_prev, rpe_prev
!
   complex (kind=CmplxKind), intent(out) :: BmatResidual(ldp*ldp,ldp,n_spin_pola,LocalNumAtoms)
   complex (kind=CmplxKind), intent(out) :: RmatResidual(ldp*ldp,ldp,n_spin_pola,LocalNumAtoms)
   complex (kind=CmplxKind) :: e, a2l, a2lp, c, cde, ce0, det
   complex (kind=CmplxKind), pointer :: jost_mat(:,:)
   complex (kind=CmplxKind), allocatable :: s0(:,:), s1(:,:), s2(:,:)
   complex (kind=CmplxKind), allocatable :: sT(:,:), sTs(:,:), sm(:,:)
   complex (kind=CmplxKind), allocatable :: vt(:), diag(:)
   complex (kind=CmplxKind), pointer :: pv(:), evr(:), evl(:), em(:,:)
   complex (kind=CmplxKind) :: rpe(ldp), erc(ldp), bmat(ldp*ldp,ldp), rmat(ldp*ldp,ldp)
!
   if (eb > et) then
      call ErrorHandler('calQuadraticPoles','eb > et',eb,et)
   endif
!
   eGID = getGroupID('Energy Mesh')
   NumPEsInEGroup = getNumPEsInGroup(eGID)
   MyPEinEGroup = getMyPEinGroup(eGID)
   write(6,'(a,3i5)')'MyPE, MyPEinEGroup, NumPEsInEGroup = ',MyPE,MyPEinEGroup,NumPEsInEGroup
!
   WindowWidth = 4.0d0*Delta
!
   if (.not.present(CheckPoles)) then
      chkpole = .false.
   else
      chkpole = CheckPoles
   endif
!
   if (.not.present(PanelOnZero)) then
      isZeroInterval = .false.
   else
      isZeroInterval = PanelOnZero
   endif
!
   if(isZeroInterval) then
      NumWindows = 1
      MyNumWindows = 1
!     de = FOURTH*(et-eb)
      de = HALF*(et-eb)
   else
      NumWindows = int((et-eb)/WindowWidth)
      NumWindows = NumWindows - mod(NumWindows,4)
      if (NumWindows < NumPEsInEGroup) then
         NumWindows = NumPEsInEGroup
         MyNumWindows = 1
      else
         NumWindows = ceiling(NumWindows/real(NumPEsInEGroup))*NumPEsInEGroup
         MyNumWindows = NumWindows/NumPEsInEGroup
      endif
      WindowWidth = (et-eb)/real(NumWindows,kind=RealKind)
!     de = Delta
      de = WindowWidth/4.0d0
   endif
print *,'NumWindows = ',NumWindows
!
   de2 = de*TWO; dede2 = de*de*TWO
   cde = de
!
   kmax_kkr = (lkkr(1)+1)**2
   allocate( s0(1:kmax_kkr,1:kmax_kkr) )
   allocate( s1(1:kmax_kkr,1:kmax_kkr) )
   allocate( s2(1:kmax_kkr,1:kmax_kkr) )
   allocate( sm(1:kmax_kkr,1:kmax_kkr), vt(kmax_kkr), diag(kmax_kkr) )
   call initQuadraticMatrix(kmax_kkr,isGeneral)
!
   NumBoundPoles = 0; NumBPDegens = 0
   NumResPoles = 0; NumRPDegens = 0
   do id = 1,LocalNumAtoms
      if ((lkkr(id)+1)**2 /= kmax_kkr) then
         deallocate(s0, s1, s2, sT, sTs)
         kmax_kkr = (lkkr(id)+1)**2
         allocate( s0(1:kmax_kkr,1:kmax_kkr) )
         allocate( s1(1:kmax_kkr,1:kmax_kkr) )
         allocate( s2(1:kmax_kkr,1:kmax_kkr) )
         call endQuadraticMatrix()
         call initQuadraticMatrix(kmax_kkr)
      endif
!
      do is = 1,n_spin_pola
         nr = 0; nb = 0
         bpe = ZERO; rpe = CZERO
         bpe_prev = ZERO; rpe_prev = ZERO
         bpdeg = 0; rpdeg = 0
         bmat = CZERO; rmat = CZERO
         do iw = 1, MyNumWindows
            w0 = eb + (iw+MyPEInEGroup*MyNumWindows-1)*WindowWidth
            e0 = w0 + (HALF)*WindowWidth
            if (isZeroInterval) then
               e0 = ZERO
            else if ((abs(e0) < Ten2m6 .or. abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6)) then
               if (e0 < ZERO) then
                  e0 = e0 - HALF*de
               else
                  e0 = e0 + HALF*de
               endif
            endif
            ce0 = e0
!
            if(isZeroInterval) then
               s0(1:kmax_kkr,1:kmax_kkr) = CZERO
            else
               e = cmplx(e0,ZERO,kind=CmplxKind)
!              -------------------------------------------------------
               call solveSingleScattering(is, id, ce0, CZERO)
!              -------------------------------------------------------
               jost_mat => getJostMatrix()
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,jost_mat,1,s0,1)
!              -------------------------------------------------------
            endif
!
            e = cmplx(e0+de,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is, id, e, CZERO)
!           call solveSingleScattering(is, id, ce0, -cde)
!           ----------------------------------------------------------
            jost_mat => getJostMatrix()
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,jost_mat,1,s2,1)
!           ----------------------------------------------------------
!
            e = cmplx(e0-de,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is, id, e, CZERO)
!           call solveSingleScattering(is, id, ce0, cde)
!           ----------------------------------------------------------
            jost_mat => getJostMatrix()
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,jost_mat,1,s1,1)
!           ----------------------------------------------------------
!
            s1 = (s2 - jost_mat)/de2
            s2 = (s2 + jost_mat - TWO*s0)/dede2
!
            if(isZeroInterval) then
!              -------------------------------------------------------
               call solveLinearEquation(s1,s2,info)
!              -------------------------------------------------------
            else
!              -------------------------------------------------------
               call solveQuadraticEquation(s0,s1,s2,info)
!              -------------------------------------------------------
            endif
!
            if (info /= 0) then
               stop 'Error in s0, s1, s2'
            endif
!
!           ----------------------------------------------------------
            pv => getEigenValue(nv)
!           ----------------------------------------------------------
            do ie = 1, nv
               if (.false.) then ! change to .false. to turn off self-checking
!                 ====================================================
!                 Check eigenvalues and eigenvectors
!                 ====================================================
                  write(6,'(a,i5,2d15.8)')'Index, Eigenvalue = ',ie,pv(ie)
                  sm = s0 + s1*pv(ie) + s2*(pv(ie)*pv(ie))
                  evr => getEigenVector('R',ie)
                  vt = CZERO
                  do j = 1, kmax_kkr
                     do i = 1, kmax_kkr
                        vt(i) = vt(i) + sm(i,j)*evr(j)
                     enddo
                  enddo
                  do i = 1, kmax_kkr
                     err = abs(vt(i))
                     if (err > ten2m7) then
                        call ErrorHandler('calQuadraticPoles',        &
                                          'Right-side eigenvector error > 10^7',err)
                     endif
                  enddo
                  write(6,'(a)')'Right-side eigenvector passed!'
!
                  evl => getEigenVector('L',ie)
                  vt = CZERO
                  do j = 1, kmax_kkr
                     do i = 1, kmax_kkr
                        vt(j) = vt(j) + evl(i)*sm(i,j)
                     enddo
                  enddo
                  do i = 1, kmax_kkr
                     err = abs(vt(i))
                     if (err > ten2m7) then
                        call WarningHandler('calQuadraticPoles',        &
                                          'Left-side eigenvector error > 10^7',err)
                     endif
                  enddo
                  write(6,'(a)')'Left-side eigenvector passed!'
               endif
!              =======================================================
!              if (aimag(pv(ie)) > ZERO .and. real(pv(ie),kind=RealKind) + e0 > ZERO) then
               if (abs(aimag(pv(ie))) < Ten2m8 .and. real(pv(ie),kind=RealKind)+e0 < ZERO) then   ! Bound states
                  pe = real(pv(ie),kind=RealKind) + e0
                  if (pe >= w0 .and. pe <= w0+WindowWidth) then
!                    -------------------------------------------------
                     em => getEigenMatrix(ie) ! em is the residule matrix of
                                              ! integrating sm^{-1} around its eigenvalue
!                    -------------------------------------------------
                     if (size(em,1) /= kmax_kkr) then
                        call ErrorHandler('calQuadraticPoles','inconsistent matrix size',size(em,1),kmax_kkr)
                     endif
                     if (abs(pe-bpe_prev) > TEN2m6) then
                        nb = nb + 1
!                       BoundPoles(nb,is,id) = pe
                        bpe(nb) = pe
                        bpdeg(nb) = 1
!write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
!                       ----------------------------------------------
                        call zcopy(kmax_kkr*kmax_kkr,em,1,bmat(1,nb),1)
!                       ----------------------------------------------
                        bpe_prev = pe
                     else if (nb == 0) then
                        call ErrorHandler('calQuadraticPoles','bound state pe = ZERO',pe)
                     else ! In degeneracy case, em is added to bmat of the same energy
                        bpdeg(nb) = bpdeg(nb) + 1
!                       ----------------------------------------------
                        call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,bmat(1,nb),1)
!                       ----------------------------------------------
                     endif
!write(6,'(a,3i4,a,f12.8)')'MyPEinEGroup,nb,deg,pe = ',MyPEinEGroup,nb,bpdeg(nb),', ',pe
                  endif
               else if (aimag(sqrt(pv(ie)+e0)) < ZERO) then  ! Resonance states
                  pe = real(pv(ie),kind=RealKind) + e0
                  if (pe >= w0 .and. pe <= w0+WindowWidth .and. pe > ZERO) then
!                    -------------------------------------------------
                     em => getEigenMatrix(ie) ! em is the residule matrix of
                                              ! integrating sm^{-1} around its eigenvalue
!                    -------------------------------------------------
                     if (abs(pe-rpe_prev) > TEN2m6) then
                        nr = nr + 1
!                       ResPoles(nr,is,id) = pe
!                       ResWidth(nr,is,id) = aimag(sqrt(pv(ie)))**2
                        rpe(nr) = cmplx(pe,aimag(sqrt(pv(ie)))**2)
                        rpdeg(nb) = 1
!write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
!                       ----------------------------------------------
                        call zcopy(kmax_kkr*kmax_kkr,em,1,rmat(1,nr),1)
!                       ----------------------------------------------
                        rpe_prev = pe
                     else if (nr == 0) then
                        call ErrorHandler('calQuadraticPoles','resonance state pe = ZERO',pe)
                     else
                        rpdeg(nr) = rpdeg(nr) + 1
!                       ----------------------------------------------
                        call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,rmat(1,nr),1)
!                       ----------------------------------------------
                     endif
                  endif
               endif
            enddo
         enddo
!
         if (chkpole) then
            do je = 1, nb
!              e0 = BoundPoles(je,is,id)
               e0 = bpe(je)
               do ie = -10, 10
                  e = e0 + ie*0.001d0
                  call solveSingleScattering(is, id, e, CZERO)
                  jost_mat => getJostMatrix()
                  call calcDet(jost_mat,kmax_kkr,det,diag)
                  write(6,'(a,f12.5,2x,2d16.8)')'e,det = ',real(e),det
               enddo
               write(6,'(/)')
            enddo
            do je = 1, nr
!              e0 = ResPoles(je,is,id)
               e0 = rpe(je)
               do ie = -10, 10
                  e = e0 + ie*0.001d0
                  call solveSingleScattering(is, id, e, CZERO)
                  jost_mat => getJostMatrix()
                  call calcDet(jost_mat,kmax_kkr,det,diag)
                  write(6,'(a,f12.5,2x,2d16.8)')'e,det = ',real(e),det
               enddo
               write(6,'(/)')
            enddo
         endif
!
         ebr = ZERO; erc = CZERO
         do ip = 1, NumPEsInEGroup
            if (MyPEinEGroup == ip-1) then
               nbr(1) = nb
               nbr(2) = nr
            endif
            call bcastMessageInGroup(eGID,nbr,2,ip-1)
            if (nbr(1) > 0) then
               nb0 = NumBoundPoles(is,id)
               if (MyPEinEGroup == ip-1) then
                  ebr(1:nbr(1)) = bpe(1:nbr(1))
!                 ----------------------------------------------------
!                 call zcopy(ldp*ldp*ldp,bmat,1,brmat,1)
!                 do ib = 1, nbr(1)
!                    -------------------------------------------------
!                    call zcopy(ldp*ldp,bmat(1,ib),1,BmatResidual(1,nb0+ib,is,id),1)
!                    -------------------------------------------------
!                 enddo
!                 ----------------------------------------------------
                  call zcopy(ldp*ldp*nbr(1),bmat,1,BmatResidual(1,nb0+1,is,id),1)
!                 ----------------------------------------------------
                  degens(1:nbr(1)) = bpdeg(1:nbr(1))
               endif
!              -------------------------------------------------------
               call bcastMessageInGroup(eGID,ebr,nbr(1),ip-1)
!              call bcastMessageInGroup(eGID,brmat,ldp*ldp,ldp,ip-1)
               call bcastMessageInGroup(eGID,degens,nbr(1),ip-1)
!              -------------------------------------------------------
               do ib = 1, nbr(1)
                  BoundPoles(nb0+ib,is,id) = ebr(ib)
                  NumBPDegens(nb0+ib,is,id) = degens(ib)
!                 ----------------------------------------------------
!                 call zcopy(kmax_kkr*kmax_kkr,brmat(1,ib),1,          &
!                            BmatResidual(1,nb0+ib,is,id),1)
!                 ----------------------------------------------------
               enddo
!              -------------------------------------------------------
               call bcastMessageInGroup(eGID,BmatResidual(1:,nb0+1:,is,id),&
                                        ldp*ldp,nbr(1),ip-1)
!              -------------------------------------------------------
            endif
            if (nbr(2) > 0) then
               nr0 = NumResPoles(is,id)
               if (MyPEinEGroup == ip-1) then
                  erc(1:nbr(2)) = rpe(1:nbr(2))
!                 ----------------------------------------------------
!                 call zcopy(ldp*ldp*ldp,rmat,1,brmat,1)
                  call zcopy(ldp*ldp*nbr(2),rmat,1,RmatResidual(1,nr0+1,is,id),1)
!                 ----------------------------------------------------
                  degens(1:nbr(2)) = rpdeg(1:nbr(2))
               endif
!              -------------------------------------------------------
               call bcastMessageInGroup(eGID,erc,nbr(2),ip-1)
!              call bcastMessageInGroup(eGID,brmat,ldp*ldp,nbr(2),ip-1)
               call bcastMessageInGroup(eGID,degens,nbr(2),ip-1)
!              -------------------------------------------------------
               do ir = 1, nbr(2)
                  ResPoles(nr0+ir,is,id) = real(erc(ir),kind=RealKind)
                  ResWidth(nr0+ir,is,id) = aimag(erc(ir))
                  NumRPDegens(nr0+ir,is,id) = degens(ir)
!                 ----------------------------------------------------
!                 call zcopy(kmax_kkr*kmax_kkr,brmat(1,ir),1,          &
!                            RmatResidual(1,nr0+ir,is,id),1)
!                 ----------------------------------------------------
               enddo
!              -------------------------------------------------------
               call bcastMessageInGroup(eGID,RmatResidual(1:,nr0+1:,is,id),&
                                        ldp*ldp,nbr(2),ip-1)
!              -------------------------------------------------------
            endif
            NumBoundPoles(is,id) = NumBoundPoles(is,id) + nbr(1)
            NumResPoles(is,id) = NumResPoles(is,id) + nbr(2)
         enddo
      enddo
   enddo
!
   nullify(pv, evl, evr, em)
   call endQuadraticMatrix()
!
   deallocate( s0, s1, s2 )
   deallocate( sm, vt )
!
   end subroutine calQuadraticPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calBoundStateDensity(LocalNumAtoms,lmax_kkr,lmax_rho,   &
                                   n_spin_pola,max_nr,kmax_kkr_max,   &
                                   jmax_rho_max,max_bs,               &
                                   NumBoundPoles, BoundPoles, Bmat,   &
                                   Bdensity)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE, CZERO
!
   use MPPModule, only : syncAllPEs, MyPE
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEInGroup
   use GroupCommModule, only : GlobalSumInGroup
!
   use SSSolverModule, only : solveSingleScattering, getSineMatrix
   use SSSolverModule, only : getRegSolution, getSolutionRmeshSize
!  use SSSolverModule, only : getJostInvMatrix, getOmegaHatMatrix
!  use SSSolverModule, only : computeDOS, getDOS
!
   use RadialGridModule, only : getGrid
!
   use MatrixModule, only : computeAStarTInv
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk, bofk, m1m
!
   use GauntFactorsModule, only : getK3, getNumK3, getGauntFactor
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: lmax_kkr(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: lmax_rho(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: n_spin_pola,kmax_kkr_max,jmax_rho_max,max_nr,max_bs
   integer (kind=IntKind), intent(in) :: NumBoundPoles(n_spin_pola,LocalNumAtoms)
!
   integer (kind=IntKind) :: id, ie, is, nr, ib, info
   integer (kind=IntKind) :: kmax_kkr
   integer (kind=IntKind) :: jmax_rho, lmax_rho_max, kmax_rho_max, kmax_rho
   integer (kind=IntKind) :: kl, klp, klp_bar, kl1, kl2, kl3, kl3_bar, m3, mp, ir, i2, jl3
   integer (kind=IntKind) :: eGID, NumPEsInEGroup, MyPEInEGroup
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
!
   real (kind=RealKind), intent(in) :: BoundPoles(kmax_kkr_max,n_spin_pola,LocalNumAtoms)
   real (kind=RealKind), pointer :: cgnt(:,:,:)
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), allocatable :: gaunt(:,:,:)
!
   complex (kind=CmplxKind), intent(in) :: Bmat(kmax_kkr_max**2,kmax_kkr_max,n_spin_pola,LocalNumAtoms)
   complex (kind=CmplxKind), intent(out) :: Bdensity(max_nr,jmax_rho_max,max_bs,n_spin_pola,LocalNumAtoms)
   complex (kind=CmplxKind), allocatable, target :: wspace0(:), wspace1(:), wspace2(:), wspace3(:)
   complex (kind=CmplxKind), pointer :: sine_mat(:,:), smat_inv(:,:), BSinv(:,:)
   complex (kind=CmplxKind), pointer :: PhiLr(:,:,:), BPhiLr(:,:,:), PPr(:,:,:), PPrSum(:,:)
!  complex (kind=CmplxKind), pointer :: jost_inv(:,:)
!  complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind) :: e, cfac, cfac0, cfac1, cfac2, kappa
!
   type (GridStruct), pointer :: Grid
!
   logical :: non_zero
!
   lmax_rho_max = 0
   do id = 1,LocalNumAtoms
      lmax_rho_max = max(lmax_rho_max,lmax_rho(id))
   enddo
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   kmax_rho_max = (lmax_rho_max+1)**2
   allocate( gaunt(kmax_kkr_max,kmax_kkr_max,kmax_rho_max) )
   gaunt = ZERO
   do kl3 = 1, kmax_rho_max
      do kl1 = 1, kmax_kkr_max
         do i2 = 1, nj3(kl1,kl3)
            kl2 = kj3(i2,kl1,kl3)
            if (kl2 <= kmax_kkr_max) then
               gaunt(kl2,kl1,kl3) = cgnt(i2,kl1,kl3)
            endif
         enddo
      enddo
   enddo
!
!  ===================================================================
!  calculate the charge density associated with each bound state
!  ===================================================================
   allocate( wspace0(kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace1(kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace2(max_nr*kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace3(max_nr*kmax_kkr_max*kmax_kkr_max) )
   eGID = getGroupID('Energy Mesh')
   NumPEsInEGroup = getNumPEsInGroup(eGID)
   MyPEinEGroup = getMyPEinGroup(eGID)
   Bdensity = CZERO
   do id = 1,LocalNumAtoms
      kmax_kkr = (lmax_kkr(id)+1)**2
      jmax_rho = (lmax_rho(id)+1)*(lmax_rho(id)+2)/2
      kmax_rho = (lmax_rho(id)+1)**2
      smat_inv => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
      BSinv => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
!
      Grid => getGrid(id)
      r_mesh => Grid%r_mesh
      nr = getSolutionRmeshSize()
      BPhiLr => aliasArray3_c(wspace2,nr,kmax_kkr,kmax_kkr)
      PPrSum => aliasArray2_c(wspace2,nr,kmax_rho)
      PPr => aliasArray3_c(wspace3,nr,kmax_kkr,kmax_kkr)
!
      do is = 1,n_spin_pola
         do ib = MyPEinEGroup+1, NumBoundPoles(is,id), NumPEsInEGroup
            e = BoundPoles(ib,is,id)
!e = -e*0.2d0
            kappa = sqrt(e)
            cfac0 = HALF*kappa
!           ----------------------------------------------------------
            call solveSingleScattering(is, id, e, CZERO)
!           ----------------------------------------------------------
            sine_mat => getSineMatrix()
!jost_inv => getJostInvMatrix()
!           ==========================================================
!           calculate sine_mat^(-T*) and store the result in smat_inv
!           ----------------------------------------------------------
            call computeAStarTInv(sine_mat,kmax_kkr,kmax_kkr,smat_inv)
!           ----------------------------------------------------------
!
!           =============================================================
!           calculate BmatResidual*sine_mat^{-T*} and store the result in BSinv
!           ----------------------------------------------------------
            call zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,       &
                       Bmat(1,ib,is,id),kmax_kkr,smat_inv,kmax_kkr,   &
                       CZERO,BSinv,kmax_kkr)
!           ----------------------------------------------------------
!call zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,jost_inv,kmax_kkr,smat_inv,kmax_kkr,CZERO,BSinv,kmax_kkr)
            PhiLr => getRegSolution()
!           ----------------------------------------------------------
            call zgemm('n','n',nr*kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                       PhiLr,nr*kmax_kkr,BSinv,kmax_kkr,              &
                       CZERO,BPhiLr,nr*kmax_kkr)
!           ----------------------------------------------------------
            PPr = CZERO
            do klp = 1, kmax_kkr
               mp = mofk(klp)
               klp_bar = bofk(klp)
               cfac = m1m(mp)
                do kl2 = 1, kmax_kkr
                   do kl1 = 1, kmax_kkr
                     do ir = 1, nr
                        PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*BPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar)
                     enddo
                  enddo
               enddo
            enddo
            do jl3 = 1, jmax_rho
               m3 = mofj(jl3)
               kl3 = kofj(jl3)
               kl3_bar = bofk(kl3)
               do kl2 = 1, kmax_kkr
                  do kl1 = 1, kmax_kkr
                     cfac1 = cfac0*gaunt(kl1,kl2,kl3)
                     cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
                     do ir = 1, nr
                        Bdensity(ir,jl3,ib,is,id) = Bdensity(ir,jl3,ib,is,id) &
                                                  + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
!+ cfac1*PPr(ir,kl1,kl2) - conjg(cfac2*PPr(ir,kl1,kl2))
                     enddo
                  enddo
               enddo
!non_zero = .false.
!LOOP_ir: do ir = 1, nr
!if (abs(Bdensity(ir,jl3,ib,is,id)) > TEN2m6) then
!non_zero = .true.
!exit LOOP_ir
!endif
!enddo LOOP_ir
!if (non_zero) then
!write(6,'(a,3i5)')'None zero Bdensity component: ib,lp,mp = ',ib,lofj(jl3),mofj(jl3)
!endif
            enddo
!
!call computeDOS()
!dos_r_jl => getDOS()
!do jl3 = 1, jmax_rho
!non_zero = .false.
!LOOP_ir_3: do ir = 1, Grid%jend
!if (abs(dos_r_jl(ir,jl3)) > TEN2m6) then
!non_zero = .true.
!exit LOOP_ir_3
!endif
!enddo LOOP_ir_3
!if (non_zero) then
!write(6,'(a,3i5)')'None zero dos_r_jl component: ib,lp,mp = ',ib,lofj(jl3),mofj(jl3)
!write(6,'(a,2d15.8,2x,2d15.8)')'den,dos = ',Bdensity(1000,jl3,ib,is,id),dos_r_jl(1000,jl3)
!endif
!enddo
         enddo
!        ------------------------------------------------------------
         call GlobalSumInGroup(eGID,Bdensity(:,:,:,is,id),           &
                               max_nr,jmax_rho_max,max_bs)
!        ------------------------------------------------------------
      enddo
   enddo
!  ------------------------------------------------------------------
   call syncAllPEs()
!  ------------------------------------------------------------------
!
   deallocate( wspace0, wspace1, wspace2, wspace3, gaunt )
   nullify(sine_mat, BSinv, smat_inv, Grid, r_mesh, BPhiLr, PhiLr, PPr, PPrSum)
!
contains
!
   include '../lib/arrayTools.F90'
!
   end subroutine calBoundStateDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calcDet(mat,kmax,det,diag)
!  ===================================================================
   use KindParamModule, only : RealKind, CmplxKind, IntKind
!
   use MathParamModule, only : CONE
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kmax
   integer (kind=IntKind) :: kl
!
   complex (kind=CmplxKind), intent(out) :: det
   complex (kind=CmplxKind), intent(out) :: diag(kmax)
   complex (kind=CmplxKind), intent(in) :: mat(kmax,kmax)
   complex (kind=CmplxKind) :: matU(kmax,kmax)
!
!  ----------------------------------------------------------
   call zcopy(kmax*kmax,mat,1,matU,1)
!  ----------------------------------------------------------
   call GaussianElim(matU,kmax)
!  ----------------------------------------------------------
   det = CONE
   do kl = 1,kmax
      det = det*matU(kl,kl)
      diag(kl) = matU(kl,kl)
   enddo
!
   end subroutine calcDet
!  ===================================================================
