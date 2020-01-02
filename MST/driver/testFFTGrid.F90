program testFFTGrid
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
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : getVolume, getInscrSphVolume
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
   use PolyhedraModule, only : getWignerSeitzRadius
   use PolyhedraModule, only : getNeighborDistance
   use PolyhedraModule, only : printPolyhedron
   use PolyhedraModule, only : printPolyhedronBoundary
   use PolyhedraModule, only : getPointLocationFlag
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
   use AtomModule, only : getGridData, getLocalEvecOld, getAtomMuffinTinRad
   use AtomModule, only : getLocalAtomPosition
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
   use RadialGridModule, only : getGrid
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
   use TimerModule, only : getTime
!  *******************************************************************
!  ******  For different testing purposes, add necessary   *********** 
!  ******  module usage after the following line.          ***********
!  *******************************************************************
   use SystemModule, only : getUniformGridParam, getBravaisLattice
   use PublicTypeDefinitionsModule, only : UniformGridStruct
   use Uniform3DGridModule, only : initUniform3DGrid, endUniform3DGrid, printUniform3DGrid
   use Uniform3DGridModule, only : createUniform3DGrid, insertAtomsInGrid
   use Uniform3DGridModule, only : getUniform3DGrid, createProcessorMesh
   use Uniform3DGridModule, only : distributeUniformGrid
!  use FFTTransformModule, only : initFFTTransform, endFFTTransform
!  use FFTTransformModule, only : getProcessorMesh
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
   integer (kind=IntKind) :: grid_start(3), grid_end(3)
!
   integer (kind=IntKind), parameter :: funit = 11
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: LocalAtomPosi(:,:)
   real (kind=RealKind), allocatable :: evec(:,:)
   real (kind=RealKind), pointer :: bravais(:,:)
!
   complex (kind=CmplxKind) :: energy
!  *******************************************************************
!  ******  For different testing purposes, add necessary   *********** 
!  ******  variables after the following line.             ***********
!  *******************************************************************
!
   real (kind=RealKind) :: t0, t1, t2, t3
!
   integer (kind=IntKind) :: i, id, ig, is, ie_loc, ie_glb
   integer (kind=IntKind) :: eGID, aGID, NumPEsInGroup, comm
!
   type (UniformGridStruct), pointer:: fftgp
   integer (kind=IntKind) :: ng_uniform(3)
   real (kind=RealKind), allocatable :: radius(:)
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
   call initUniform3DGrid(istop, node_print_level)
   ng_uniform(1:3) = getUniformGridParam()
   t1 = getTime()
!  -------------------------------------------------------------------
   call createUniform3DGrid('FFT', ng_uniform(1), ng_uniform(2), ng_uniform(3), bravais)
!  -------------------------------------------------------------------
!  call initFFTTransform(fftgp%nga, fftgp%ngb, fftgp%ngc, bravais)
!  -------------------------------------------------------------------
!  call createProcessorMesh('FFT', getProcessorMesh())
!  grid_start = 1; grid_end = ng_uniform
!  call distributeUniformGrid('FFT',grid_start,grid_end)
   call distributeUniformGrid('FFT')
!  -------------------------------------------------------------------
!
   allocate(radius(LocalNumAtoms))
   do i = 1, LocalNumAtoms
      radius(i) = getOutscrSphRadius(i)
      write(6,'(a,i5,2x,d15.8)')'i,radius = ',i,radius(i)
   enddo
!
!  -------------------------------------------------------------------
   call insertAtomsInGrid('FFT', LocalNumAtoms, LocalAtomPosi,        &
                          getPointLocationFlag, radius)
!  -------------------------------------------------------------------
   write(6,'(/,a,f12.5,a)') 'set FFT grid time: ',getTime() - t1,' sec.'
   if (node_print_level >= 0) then
      call printUniform3DGrid('FFT')
   endif
   fftgp => getUniform3DGrid('FFT')
!
!  ===================================================================
!  finish the process
!  -------------------------------------------------------------------
   deallocate(radius)
   call endUniform3DGrid()
!  call endFFTTransform()
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
   integer (kind=IntKind) :: lmax_step_max, lmax_kkr_max, lmax_phi_max, lmax_rho_max, lmax_pot_max, lmax_green_max
   integer (kind=IntKind) :: ndivin, ndivout, nmult, kmax_kkr, lmax_max
!
   real (kind=RealKind) :: rmt, rend, rws, rinsc, v0, Efermi
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
   allocate(LocalAtomPosi(3,LocalNumAtoms))
   do i=1,LocalNumAtoms
      atom_print_level(i) = getStandardOutputLevel(i)
      LocalAtomPosi(1:3,i)=getLocalAtomPosition(i)
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
!  ===================================================================
!  initialize radial grid
!  -------------------------------------------------------------------
   call initRadialGrid(LocalNumAtoms, istop, node_print_level)
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
!  -------------------------------------------------------------------
!   call initPolyhedra(NumAtoms,bravais,'main',0)
!  -------------------------------------------------------------------
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i)
      GlobalIndex(i)=ig
!     ----------------------------------------------------------------
      call getGridData(i,ndivin,ndivout,nmult)
!     ----------------------------------------------------------------
!     call genPolyhedron(i,ig,NumAtoms,AtomPosition)
!     ----------------------------------------------------------------
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printPolyhedron(i)
!        call printPolyhedronBoundary(i)
!        -------------------------------------------------------------
      endif
      rend =  getOutscrSphRadius(i)
      if (isMuffinTinPotential() .or. isMuffinTinTestPotential()) then
         rmt = getAtomMuffinTinRad(i)
         rinsc = getInscrSphRadius(i)
         if ( rmt < 0.010d0 ) then
            rmt = rinsc
         endif
         rws = getWignerSeitzRadius(i)
         if (getSingleSiteSolverType()==1) then
            rend=rws
         endif
!        -------------------------------------------------------------
         call genRadialGrid(i,xstart, rmt, rinsc, rend, ndivin)
!        -------------------------------------------------------------
      else if ( isASAPotential() ) then
         rend =  getWignerSeitzRadius(i)
         rmt = getAtomMuffinTinRad(i)
         rinsc = getWignerSeitzRadius(i)
         if ( rmt < 0.010d0 ) then
            rmt = rinsc
         endif
!        -------------------------------------------------------------
         call genRadialGrid(i,xstart, rmt, rinsc, rend, ndivin )
!        -------------------------------------------------------------
      else if (isMuffinTinASAPotential()) then
         rend =  getWignerSeitzRadius(i)
         rmt = getAtomMuffinTinRad(i)
         rinsc = getWignerSeitzRadius(i)
         if ( rmt < 0.010d0 ) then
            rmt = rinsc
         endif
!        -------------------------------------------------------------
         call genRadialGrid( i, xstart, rmt, rinsc, rend, ndivin )
!        -------------------------------------------------------------
      else
         if (getNeighborDistance(i,1)-getOutscrSphRadius(i) < TEN2m8) then
!           ----------------------------------------------------------
            call WarningHandler('testSSSolver',                       &
                     'Ill condition found: Neighbor distance <= Rcs', &
                     getNeighborDistance(i,1),getOutscrSphRadius(i))
!           ----------------------------------------------------------
         endif
         rmt = getAtomMuffinTinRad(i)
         rinsc = getInscrSphRadius(i)
         if ( rmt < 0.010d0 ) then
            rmt = getInscrSphRadius(i)
         endif
         rws = getWignerSeitzRadius(i)
!        -------------------------------------------------------------
         call genRadialGrid( i, rmt, rws, rinsc, rend, ndivin, ndivout, nmult)
!        call genRadialGrid(i,getInscrSphRadius(i),getOutscrSphRadius(i), &
!                           ndivin,ndivout,nmult)
!        -------------------------------------------------------------
      endif
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printRadialGrid(i)
!        -------------------------------------------------------------
      endif
   enddo
!
!  ===================================================================
!  initialize step function module
!  ===================================================================
   allocate( ngr(LocalNumAtoms), ngt(LocalNumAtoms) )
   do i=1,LocalNumAtoms
      ngr(i) = ngaussr
      ngt(i) = ngaussq
   enddo
!
!  -------------------------------------------------------------------
   call initStepFunction(LocalNumAtoms, lmax_max, lmax_step, ngr, ngt, &
                         istop,node_print_level)
!  -------------------------------------------------------------------
   deallocate( ngr, ngt )
!
   do i=1,LocalNumAtoms
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printStepFunction(i)
!        -------------------------------------------------------------
      endif
!     ----------------------------------------------------------------
      call testStepFunction(i)
!     ----------------------------------------------------------------
   enddo
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
end program testFFTGrid
