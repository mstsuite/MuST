program findResonance
!  *******************************************************************
!  *********** Headers of modules usage. No need to modify ***********
!  *******************************************************************
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8, TEN2m6, TEN2m3, FIVE, CZERO, CONE, TWO, PI, PI2, SQRTm1, HALF, &
                               PI4, THIRD
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use DataServiceCenterModule, only : isDataStorageExisting, &
                                       getDataStorage, RealMark
!
   use OutputModule, only : getStandardOutputLevel
!
!  ===================================================================
!  initPolyhedra and endPolyhedra are called in SystemVolumeModule, called
!  in startProcess()
!  ===================================================================
   use PolyhedraModule, only : getVolume, getInscrSphVolume
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
   use PolyhedraModule, only : getWignerSeitzRadius
   use PolyhedraModule, only : getNeighborDistance
   use PolyhedraModule, only : printPolyhedron
   use PolyhedraModule, only : printPolyhedronBoundary
!
   use SurfElementsModule, only : initSurfElements, endSurfElements,  &
                                  genGaussPoints, genGaussSphHarmonics
!
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : printStepFunction, testStepFunction
   use StepFunctionModule, only : getVolumeIntegration
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential, setPotentialOutsideMT, getPotential
   use PotentialModule, only : getPotEf, setV0, setPotential, isPotComponentZero, getSphPotr
!
   use ScfDataModule, only : ngaussr, ngaussq
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop
   use ScfDataModule, only : ErBottom, ErTop, EiBottom, EiTop, Temperature
   use ScfDataModule, only : isNonRelativisticValence
   use ScfDataModule, only : isScalarRelativisticValence
   use ScfDataModule, only : getSingleSiteSolverType, getSingleSiteSolverMethod
!
   use PotentialTypeModule, only : isASAPotential, isMuffinTinPotential,     &
                                   isMuffinTinASAPotential, isTestPotential, &
                                   isMuffinTinTestPotential, isFullPotential,&
                                   isTestPotential
!
   use TestPotentialModule, only : initTestPotential, endTestPotential,      &
                                   readTestPotential, getTestPotential,      &
                                   getTestPotEf, getTestV0
!
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,   &
                                    calSymmetryFlags
!
   use MPPModule, only : packMessage, unpackMessage, sendPackage, recvPackage, AnyPE, MyPE
   use MPPModule, only : setCommunicator, resetCommunicator, GlobalMax, bcastMessage
   use GroupCommModule, only : getMyPEinGroup,  getGroupID, getNumPEsInGroup, getGroupCommunicator
!
!  ===================================================================
!  initAtom2Proc and endAtom2Proc are called in SystemModule
!  ===================================================================
   use Atom2ProcModule, only : getGlobalIndex, getLocalNumAtoms, getLocalIndex
!
   use AtomModule, only : getStepFuncLmax, setTruncPotLmax, setPotLmax
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax
   use AtomModule, only : getGridData, getLocalEvecOld
   use AtomModule, only : getLocalNumSpecies, getLocalAtomicNumber
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
   use ContourModule, only : getNumEs, getEPoint, setupContour
!
   use ProcMappingModule, only : getNumEsOnMyProc, getEnergyIndex,    &
                                 getNumRedundantEsOnMyProc,           &
                                 getProcWithEnergyIndex
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use RadialGridModule, only : getGrid, endRadialGrid
!
   use TimerModule, only : getTime
!  *******************************************************************
!  ******  For different testing purposes, add necessary   *********** 
!  ******  module usage after the following line.          ***********
!  *******************************************************************
   use SSSolverModule, only : initSSSolver, endSSSolver, solveSingleScattering, &
                              computeDOS, computePDOS, getDOS, getPDOS, getJostMatrix
   use SSSolverModule, only : getSineMatrix, getCosineMatrix, getTMatrix, getSMatrix
!
   use ScfDataModule, only : pole_step
!
   implicit none
!
!  *******************************************************************
   integer (kind=IntKind) :: LocalNumAtoms, NumAtoms, NumEsOnMyProc, NumEs
   integer (kind=IntKind) :: RelativisticFlag
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_green(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
!
   integer (kind=IntKind), parameter :: funit = 11
   integer (kind=IntKind) :: lmax_step_max, lmax_kkr_max, lmax_phi_max, lmax_rho_max, lmax_pot_max, lmax_green_max
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: evec(:,:)
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind) :: Efermi
!
   complex (kind=CmplxKind) :: energy
!  *******************************************************************
!  ******  For different testing purposes, add necessary   *********** 
!  ******  variables after the following line.             ***********
!  *******************************************************************
!
   real (kind=RealKind) :: t0, t1, t2, t3
!
   complex (kind=CmplxKind), pointer :: p_emesh(:)
!
   integer (kind=IntKind) :: i, id, ig, is, ie_loc, ie_glb
   integer (kind=IntKind) :: eGID, aGID, NumPEsInGroup, comm
!
   integer (kind=IntKind), allocatable :: NumBoundPoles(:,:), NumResPoles(:,:)
   real (kind=RealKind), allocatable :: BoundPoles(:,:,:), ResPoles(:,:,:), ResWidth(:,:,:)
   real (kind=RealKind), parameter :: ZeroIntHalfWidth = TEN2m6
!
   interface
      subroutine calQuadraticPoles(LocalNumAtoms,lkkr,lphi,       &
                                   n_spin_pola,eb,et,del,ldp,     &
                                   NumBoundPoles,BoundPoles,      &
                                   NumResPoles,ResPoles,ResWidth, &
                                   CheckPoles,PanelOnZero)
      use KindParamModule, only : IntKind, RealKind, CmplxKind
      implicit none
      integer (kind=IntKind), intent(in) :: LocalNumAtoms
      integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
      integer (kind=IntKind), intent(in) :: lphi(1:LocalNumAtoms)
      integer (kind=IntKind), intent(in) :: n_spin_pola,ldp
      integer (kind=IntKind), intent(out) :: NumBoundPoles(LocalNumAtoms,n_spin_pola)
      integer (kind=IntKind), intent(out) :: NumResPoles(LocalNumAtoms,n_spin_pola)
!
      logical, optional, intent(in) :: CheckPoles
      logical, optional, intent(in) :: PanelOnZero
!
      real (kind=RealKind), intent(in) :: eb, et, del
      real (kind=RealKind), intent(out) :: BoundPoles(ldp,LocalNumAtoms,n_spin_pola)
      real (kind=RealKind), intent(out) :: ResPoles(ldp,LocalNumAtoms,n_spin_pola)
      real (kind=RealKind), intent(out) :: ResWidth(ldp,LocalNumAtoms,n_spin_pola)
!
      end subroutine calQuadraticPoles
   end interface
!
!  ===================================================================
!  Initialize modules and data
!  -------------------------------------------------------------------
   call startProcess()
!  -------------------------------------------------------------------
   t0 = getTime()
!  -------------------------------------------------------------------
   call initializeModules()
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, &
                     lmax_kkr, lmax_phi, lmax_pot, lmax_step, lmax_green,  &
                     n_spin_pola, n_spin_cant, RelativisticFlag, istop, atom_print_level)
!  -------------------------------------------------------------------
!
   write(6,'(/,a,f12.5,a)') 'Setup time: ',getTime() - t0,' sec.'
!
   aGID = getGroupID('Unit Cell')
   eGID = getGroupID('Energy Mesh')
!
   t1 = getTime()
!
!   call initFindPolesModule( LocalNumAtoms, MyPE, n_spin_pola, n_spin_cant, &
!                             lmax_kkr, lmax_phi, ZERO, Efermi, pole_step)
!  -------------------------------------------------------------------
!  call calScatteringPoles(LocalNumAtoms,lmax_kkr,lmax_phi,           &
!                          n_spin_pola,n_spin_cant,ZERO,Efermi,ZERO,nstep)
!   call findPoles()
!  -------------------------------------------------------------------
!   call endFindPolesModule()
!  write(6,'(/,a,f12.5,a)') 'Zooming technique time: ',getTime() - t1,' sec.'
!  -------------------------------------------------------------------
   allocate( NumBoundPoles(LocalNumAtoms,n_spin_pola) )
   allocate( NumResPoles(LocalNumAtoms,n_spin_pola) )
   allocate( BoundPoles((lmax_kkr_max+1)**2,LocalNumAtoms,n_spin_pola) )
   allocate( ResPoles((lmax_kkr_max+1)**2,LocalNumAtoms,n_spin_pola) )
   allocate( ResWidth((lmax_kkr_max+1)**2,LocalNumAtoms,n_spin_pola) )
   t1 = getTime()
!  ------------------------------------------------------------------
   call calQuadraticPoles(LocalNumAtoms,lmax_kkr,lmax_phi,            &
                          n_spin_pola,ErBottom,Efermi,pole_step,      &
                          (lmax_kkr_max+1)**2,                        &
                          NumBoundPoles,BoundPoles,                   &
                          NumResPoles,ResPoles,ResWidth,              &
                          CheckPoles =.false.)
!  ------------------------------------------------------------------
   write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
   do is = 1,n_spin_pola
      do id = 1,LocalNumAtoms
         write(6,'(a,i3,a,i3)')'Spin index: ',is,',  Atom index: ',id
         write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of bound states found within (',  &
                                          ErBottom,', ',ZERO,'): ',NumBoundPoles(id,is)
         do i = 1, NumBoundPoles(id,is)
            write(6,'(f20.12)')BoundPoles(i,id,is)
         enddo
         write(6,'(/,a,f13.8,a,f13.8,a,i5)')'The Number of resonance states found within (',  &
                                            ZERO,', ',Efermi,'): ',NumResPoles(id,is)
         do i = 1, NumResPoles(id,is)
            write(6,'(a,f20.12,a,f20.12)')'Resonance e = ',ResPoles(i,id,is),', Width = ',ResWidth(i,id,is)
         enddo
      enddo
   enddo
   write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
!
!  ===================================================================
!  clean up the allocated spaces.
!  ===================================================================
   deallocate( NumBoundPoles, BoundPoles, NumResPoles, ResPoles, ResWidth )
!
   if (node_print_level >= 0) then
      write(6,*) " Run Ending! "
      call FlushFile(6)
   endif
!
!  ===================================================================
!  clean up the allocated spaces.
!  ===================================================================
   nullify(p_emesh)
!
   call endSSSolver()
!
!  ===================================================================
!  finish the process
!  -------------------------------------------------------------------
   call finishProcess()
!  -------------------------------------------------------------------
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initializeModules()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, id, ig, is, l, m, jl, ie_glb
   integer (kind=IntKind) :: jmax_pot
   integer (kind=IntKind) :: ndivin, ndivout, nmult, kmax_kkr, lmax_max
!
   real (kind=RealKind) :: rmt, rend, rws, rinsc, v0
   real (kind=RealKind), pointer :: potr_0(:)
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!
   complex (kind=CmplxKind), pointer :: EPoint(:)
   complex (kind=CmplxKind), pointer :: pot_l(:)
!
   type (GridStruct), pointer :: Grid
!
   NumAtoms = getNumAtoms()
   LocalNumAtoms=getLocalNumAtoms()
!
   allocate(AtomPosition(1:3,1:NumAtoms), AtomicNumber(1:NumAtoms))
   allocate(GlobalIndex(LocalNumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
      AtomicNumber(i)=getAtomicNumber(i)
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
   kmax_kkr = (lmax_kkr_max + 1)**2
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
                  if (.not.isPotComponentZero(id,jl)) then
                     write(6,'(a,i2,a,i2,a)')'Non-zero potential component at (l,m) channel: (',l,',',m,')'
                  endif
               enddo
            enddo
            potr_0 =>  getSphPotr(id,1,is)
!           print *, 'vzero = ', (potr_0(1)+TWO*AtomicNumber(ig))/Grid%r_mesh(1)
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
   implicit none
!
   deallocate( atom_print_level )
   deallocate( AtomPosition,AtomicNumber,lmax_kkr,lmax_phi,lmax_rho,lmax_pot )
   deallocate( lmax_step, lmax_green, GlobalIndex, evec )
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
end program findResonance
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calQuadraticPoles(LocalNumAtoms,lkkr,lphi,              &
                                n_spin_pola,eb,et,Delta,ldp,          &
                                NumBoundPoles, BoundPoles,            &
                                NumResPoles, ResPoles, ResWidth,      &
                                CheckPoles,PanelOnZero)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                                Ten2m8, FOURTH, CZERO, TEN2m7
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use SSSolverModule, only : solveSingleScattering, getJostMatrix, getSineMatrix
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
   integer (kind=IntKind), intent(out) :: NumBoundPoles(LocalNumAtoms,n_spin_pola)
   integer (kind=IntKind), intent(out) :: NumResPoles(LocalNumAtoms,n_spin_pola)
!
   integer (kind=IntKind) :: id, ie, is, iw, kmax_kkr, NumWindows, info
   integer (kind=IntKind) :: i, j, l, m, kl, lp, mp, klp, nv, nb, nr, je
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
   real (kind=RealKind), intent(out) :: BoundPoles(ldp,LocalNumAtoms,n_spin_pola)
   real (kind=RealKind), intent(out) :: ResPoles(ldp,LocalNumAtoms,n_spin_pola)
   real (kind=RealKind), intent(out) :: ResWidth(ldp,LocalNumAtoms,n_spin_pola)
!
   complex (kind=CmplxKind) :: e, a2l, a2lp, c, cde, ce0, det
   complex (kind=CmplxKind), pointer :: jost_mat(:,:)
   complex (kind=CmplxKind), allocatable :: s0(:,:), s1(:,:), s2(:,:)
   complex (kind=CmplxKind), allocatable :: sT(:,:), sTs(:,:), sm(:,:)
   complex (kind=CmplxKind), allocatable :: vt(:), diag(:)
   complex (kind=CmplxKind), pointer :: pv(:), evr(:), evl(:), em(:,:)
!
   if (eb > et) then
      call ErrorHandler('calQuadraticPoles','eb > et',eb,et)
   endif
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
   else
      NumWindows = int((et-eb)/WindowWidth)
   endif
!
   if(isZeroInterval) then
!     de = FOURTH*(et-eb)
      de = HALF*(et-eb)
   else
      de = Delta
   endif
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
   do is = 1,n_spin_pola
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
         nr = 0; nb = 0
         do iw = 1, NumWindows
            w0 = eb + (iw-1)*WindowWidth
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
                     nb = nb + 1
                     BoundPoles(nb,id,is) = pe
!write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
                     em => getEigenMatrix(ie) ! em is the residule matrix of
                                              ! integrating sm^{-1} around its eigenvalue
                  endif
               else if (aimag(sqrt(pv(ie)+e0)) < ZERO) then  ! Resonance states
                  pe = real(pv(ie),kind=RealKind) + e0
                  if (pe >= w0 .and. pe <= w0+WindowWidth .and. pe > ZERO) then
                     nr = nr + 1
                     ResPoles(nr,id,is) = pe
                     ResWidth(nr,id,is) = aimag(sqrt(pv(ie)))**2
!write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
                     em => getEigenMatrix(ie) ! em is the residule matrix of
                                              ! integrating sm^{-1} around its eigenvalue
                  endif
               endif
            enddo
         enddo
!
         if (chkpole) then
            do je = 1, nb
               e0 = BoundPoles(je,id,is)
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
               e0 = ResPoles(je,id,is)
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
         NumBoundPoles(id,is) = nb
         NumResPoles(id,is) = nr
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
