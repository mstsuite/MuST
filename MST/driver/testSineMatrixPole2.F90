program testSineMatrixPole
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8, TEN2m6, TEN2m3, FIVE
!
   use TimerModule, only : initTimer, getTime
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : initMPP, endMPP,MyPE
!
   use GroupCommModule, only : initGroupComm, getGroupID, endGroupComm
!
   use ParallelIOModule, only : initParallelIO, endParallelIO
!
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter,  &
                                       isDataStorageExisting, &
                                       getDataStorage, RealMark
!
   use InputModule, only : initInput, endInput, readInputData
   use InputModule, only : openInputFile, closeInputFile
   use InputModule, only : getKeyValue, printKeyNames, getTableIndex
!
   use OutputModule, only : initOutput, endOutput,  &
                            getDensityPrintFlag,    & 
                            getStandardOutputLevel
!
!  ===================================================================
!  initPolyhedra and endPolyhedra are called in SystemModule
!  ===================================================================
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
   use PolyhedraModule, only : getWignerSeitzRadius
   use PolyhedraModule, only : getNeighborDistance
   use PolyhedraModule, only : printPolyhedron
   use PolyhedraModule, only : printPolyhedronBoundary
!
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : printStepFunction, testStepFunction
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential
   use PotentialModule, only : getPotEf
!
   use ScfDataModule, only : ngaussr, ngaussq
   use ScfDataModule, only : initScfData, printScfData
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop
   use ScfDataModule, only : ErBottom, ErTop, pole_step
   use ScfDataModule, only : inputpath
   use ScfDataModule, only : isNonRelativisticValence
   use ScfDataModule, only : isScalarRelativisticValence
   use ScfDataModule, only : getPotentialTypeParam
   use ScfDataModule, only : getSingleSiteSolverType
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType,      &
                                   isASAPotential, isMuffinTinPotential,     &
                                   isMuffinTinASAPotential, isTestPotential, &
                                   isMuffinTinTestPotential, isFullPotential,&
                                   printPotentialType
!
   use SystemModule, only : initSystem, endSystem
   use Systemmodule, only : printSystem
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,   &
                                    calSymmetryFlags
!
   use SystemVolumeModule, only : initSystemVolume, endSystemVolume
!
   use ProcMappingModule, only : initProcMapping, endProcMapping,            &
                                 createParallelization
!
!  ===================================================================
!  initAtom2Proc and endAtom2Proc are called in SystemModule
!  ===================================================================
   use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc
   use Atom2ProcModule, only : getGlobalIndex
   use Atom2ProcModule, only : printAtom2ProcTable, getLocalNumAtoms
!
   use AtomModule, only : initAtom, endAtom
   use AtomModule, only : printAtom, getStepFuncLmax
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax
   use AtomModule, only : getGridData, getLocalEvecOld, getAtomMuffinTinRad
   use AtomModule, only : getLocalNumSpecies, getLocalAtomicNumber
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
!
   use MadelungModule, only : initMadelung, endMadelung
!
   use SpinRotationModule, only : initSpinRotation, endSpinRotation
!
   use SSSolverModule, only : initSSSolver, endSSSolver
!
!   use FindPolesModule, only : initFindPolesModule, endFindPolesModule, &
!                                findPoles
!
   implicit none
!
   character (len=8)  :: exec_date
   character (len=10) :: exec_time
!
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: i, j, k, id, ig, nstep, is, lmax_max
   integer (kind=IntKind) :: lmax_step_max, lmax_kkr_max, lmax_phi_max, lmax_rho_max, lmax_pot_max
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind) :: LocalNumAtoms, NumAtoms
   integer (kind=IntKind) :: RelativisticFlag
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: lmax_tmp(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: NumPoles(:,:)
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind), allocatable :: Poles(:,:,:)
   real (kind=RealKind), allocatable :: evec(:,:)
!
   real (kind=RealKind) :: Efermi
   real (kind=RealKind) :: t0, t1, t2, t3
   real (kind=RealKind) :: rmt, rend, rws, rinsc
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!  real (kind=RealKind), parameter :: xstart = -.913096740000D+01
   real (kind=RealKind), parameter :: ZeroIntHalfWidth = TEN2m6
!
!  -------------------------------------------------------------------
   call initTimer()
!  -------------------------------------------------------------------
   t0 = getTime()
!
!  -------------------------------------------------------------------
   call initMPP()
   call initGroupComm()
   call initDataServiceCenter()
   call initInput()
   call readInputs(def_id,info_id)
   call initScfData(def_id)
   call initPotentialType(getPotentialTypeParam())
   call initSystem(def_id)
!  -------------------------------------------------------------------
!  
   NumAtoms = getNumAtoms()
   if (NumAtoms < 1) then
      call ErrorHandler('main','invalid NumAtoms',NumAtoms)
   endif
!
!  ===================================================================
!  Initialize the processes mapping module that determines how the
!  parallization will be performed
!  -------------------------------------------------------------------
   call initProcMapping(NumAtoms, 1, 1, isFullPotential(), istop, 0, NumAtoms)
!  -------------------------------------------------------------------
   call createParallelization()
!  -------------------------------------------------------------------
   call initParallelIO(getGroupID('Unit Cell'),1) ! only the 1st cluster
                                                  ! in the group performs
                                                  ! writing potential data
!  -------------------------------------------------------------------
   call initAtom2Proc(NumAtoms, NumAtoms)
!  -------------------------------------------------------------------
   LocalNumAtoms=getLocalNumAtoms()
!  -------------------------------------------------------------------
   call initSystemVolume()
!  -------------------------------------------------------------------
!  ===================================================================
!
   allocate(AtomPosition(1:3,1:NumAtoms), AtomicNumber(1:NumAtoms))
   allocate(GlobalIndex(LocalNumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
      AtomicNumber(i)=getAtomicNumber(i)
   enddo
!
!  ===================================================================
!  set up print level
!  -------------------------------------------------------------------
   call initOutput(def_id)
   node_print_level = getStandardOutputLevel()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  call initAtom and setAtomData to setup Atom Module.................
!  -------------------------------------------------------------------
   call initAtom(info_id,istop,node_print_level)
!  -------------------------------------------------------------------
   if (node_print_level >= 0) then
!     ----------------------------------------------------------------
      call printPotentialType()
      call printScfData()
      call printSystem()
      call printAtom()
!     ----------------------------------------------------------------
   endif
!
   allocate(atom_print_level(1:LocalNumAtoms))
   do i=1,LocalNumAtoms
      atom_print_level(i) = getStandardOutputLevel(i)
   enddo
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  setup the Lmax values.             
!  ===================================================================
   allocate(lmax_pot(LocalNumAtoms), lmax_rho(LocalNumAtoms))
   allocate(lmax_kkr(LocalNumAtoms), lmax_phi(LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms), lmax_tmp(LocalNumAtoms))
   lmax_kkr_max = 0
   lmax_phi_max = 0
   lmax_rho_max = 0
   lmax_pot_max = 0
   lmax_step_max = 0
   lmax_max = 0
   do i=1,LocalNumAtoms
      lmax_kkr(i) = getKKRLmax(i)
      lmax_phi(i) = getPhiLmax(i)
      lmax_rho(i) = getRhoLmax(i)
      lmax_pot(i) = getPotLmax(i)
      lmax_step(i)  = getStepFuncLmax(i)
!      lmax_step(i) = max(2*lmax_phi(i), lmax_pot(i)+2*lmax_phi(i))
      lmax_kkr_max = max(lmax_kkr_max,lmax_kkr(i))
      lmax_phi_max = max(lmax_phi_max,lmax_phi(i))
      lmax_rho_max = max(lmax_rho_max,lmax_rho(i))
      lmax_pot_max = max(lmax_pot_max,lmax_pot(i))
      lmax_step_max = max(lmax_step_max,lmax_step(i))
      lmax_max = max( lmax_max, lmax_step(i)+lmax_pot(i), lmax_rho(i),&
                      2*lmax_phi(i) )
   enddo
!
   do i=1,LocalNumAtoms
      lmax_tmp(i) = max(lmax_rho(i), lmax_pot(i), lmax_step(i)+lmax_pot(i))
   enddo
!  -------------------------------------------------------------------
   call initSphericalHarmonics(4*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(2*lmax_max,istop,0)
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
            call WarningHandler('main',                               &
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
         call genRadialGrid( i, rmt, rinsc, rws, rend, ndivin, ndivout, nmult)
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
!  -------------------------------------------------------------------
   call initSystemSymmetry( NumAtoms, LocalNumAtoms, lmax_pot, lmax_step, atom_print_level )
!  -------------------------------------------------------------------
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  initialize potential module
!  -------------------------------------------------------------------
   call initMadelung(LocalNumAtoms,NumAtoms,GlobalIndex,              &
                     lmax_rho_max,lmax_max,bravais,AtomPosition,0)
!  -------------------------------------------------------------------
   call calSymmetryFlags()
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot,lmax_step,      &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  read potential data
!  -------------------------------------------------------------------
   call readPotential()
!  -------------------------------------------------------------------
   if (ErTop > ZERO) then
      Efermi = ErTop
   else
      Efermi=getPotEf()
   endif
!
!  -------------------------------------------------------------------
   call calSymmetryFlags()
!  -------------------------------------------------------------------
   allocate(evec(3,LocalNumAtoms))
   do i = 1,LocalNumAtoms
      evec(1:3,i) = getLocalEvecOld(i)
   enddo
!  -------------------------------------------------------------------
   call initSpinRotation(LocalNumAtoms,evec)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize Single Site Scatterer
!  ===================================================================
   if (isNonRelativisticValence()) then
      RelativisticFlag = 0
   else if (isScalarRelativisticValence()) then
      RelativisticFlag = 1
   else
      RelativisticFlag = 2
   endif
!  -------------------------------------------------------------------
   call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, lmax_kkr,    &
                            lmax_phi, lmax_pot, lmax_step, lmax_rho,  &
                            n_spin_pola, n_spin_cant,                 &
                            RelativisticFlag, istop, atom_print_level)
!  -------------------------------------------------------------------
   write(6,'(/,a,f12.5,a)') 'Setup time: ',getTime() - t0,' sec.'
!
   nstep = max(int((Efermi-ErBottom)/pole_step),1)
   write(6,'(/,a,f8.5,a,f8.5,a,i5)')'Ebot = ',ErBottom,',  Etop = ',Efermi,',  nstep = ',nstep
!  t1 = getTime()
!   call initFindPolesModule( LocalNumAtoms, MyPE, n_spin_pola, n_spin_cant, &
!                             lmax_kkr, lmax_phi, ErBottom, Efermi, pole_step)
!  -------------------------------------------------------------------
!  call calScatteringPoles(LocalNumAtoms,lmax_kkr,lmax_phi,           &
!                          n_spin_pola,n_spin_cant,ErBottom,Efermi,ZERO,nstep)
!   call findPoles()
!  -------------------------------------------------------------------
!   call endFindPolesModule()
!  write(6,'(/,a,f12.5,a)') 'Zooming technique time: ',getTime() - t1,' sec.'
!  -------------------------------------------------------------------
!  call plotSinDet(ErBottom, Efermi, 1000,                   &
!                  LocalNumAtoms, n_spin_pola, lmax_kkr)
!  -------------------------------------------------------------------
   allocate( NumPoles(LocalNumAtoms,n_spin_pola) )
   allocate( Poles((lmax_kkr_max+2)**2,LocalNumAtoms,n_spin_pola) )
!  if(ErBottom < -ZeroIntHalfWidth) then
!     t1 = getTime()
!     -------------------------------------------------------------------
!     call calQuadraticPoles(LocalNumAtoms,lmax_kkr,lmax_phi,            &
!                            n_spin_pola,ErBottom,-ZeroIntHalfWidth,     &
!                            (lmax_kkr_max+1)**2,NumPoles,Poles)
!     -------------------------------------------------------------------
!     write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
!     do is = 1,n_spin_pola
!        do id = 1,LocalNumAtoms
!           write(6,'(a,i3,a,i3)')'Spin index: ',is,',  Atom index: ',id
!           write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of poles found within (',ErBottom,', ', &
!                                             -ZeroIntHalfWidth,'): ',NumPoles(id,is)
!           do i = 1, NumPoles(id,is)
!              write(6,'(f20.12)')Poles(i,id,is)
!           enddo
!        enddo
!     enddo
!  endif
!
!  if (ErBottom <= ZeroIntHalfWidth) then
!     t1 = getTime()
!     -------------------------------------------------------------------
!     call calQuadraticPoles(LocalNumAtoms,lmax_kkr,lmax_phi,            &
!                            n_spin_pola,-ZeroIntHalfWidth,              &
!                            ZeroIntHalfWidth,(lmax_kkr_max+1)**2,       &
!                            NumPoles,Poles,.true.)
!     -------------------------------------------------------------------
!     write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
!     do is = 1,n_spin_pola
!        do id = 1,LocalNumAtoms
!           write(6,'(a,i3,a,i3)')'Spin index: ',is,',  Atom index: ',id
!           write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of poles found within (',  &
!                                             -ZeroIntHalfWidth,', ',               &
!                                             ZeroIntHalfWidth,'): ',NumPoles(id,is)
!           do i = 1, NumPoles(id,is)
!              write(6,'(f20.12)')Poles(i,id,is)
!           enddo
!        enddo
!     enddo
!
!     t1 = getTime()
!     -------------------------------------------------------------------
!     call calQuadraticPoles(LocalNumAtoms,lmax_kkr,lmax_phi,            &
!                            n_spin_pola,ZeroIntHalfWidth,Efermi,        &
!                            (lmax_kkr_max+1)**2,NumPoles,Poles)
!     -------------------------------------------------------------------
!     write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
!     do is = 1,n_spin_pola
!        do id = 1,LocalNumAtoms
!           write(6,'(a,i3,a,i3)')'Spin index: ',is,',  Atom index: ',id
!           write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of poles found within (',  &
!                                            ZeroIntHalfWidth,', ',                 &
!                                            Efermi,'): ',NumPoles(id,is)
!           do i = 1, NumPoles(id,is)
!              write(6,'(f20.12)')Poles(i,id,is)
!           enddo
!        enddo
!     enddo
!  else
!     t1 = getTime()
!     -------------------------------------------------------------------
!     call calQuadraticPoles(LocalNumAtoms,lmax_kkr,lmax_phi,            &
!                            n_spin_pola,ErBottom,Efermi,                &
!                            (lmax_kkr_max+1)**2,NumPoles,Poles)
!     -------------------------------------------------------------------
!     write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
!     do is = 1,n_spin_pola
!        do id = 1,LocalNumAtoms
!           write(6,'(a,i3,a,i3)')'Spin index: ',is,',  Atom index: ',id
!           write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of poles found within (',  &
!                                            ErBottom,', ',                 &
!                                            Efermi,'): ',NumPoles(id,is)
!           do i = 1, NumPoles(id,is)
!              write(6,'(f20.12)')Poles(i,id,is)
!           enddo
!        enddo
!     enddo
!  endif
   t1 = getTime()
!  -------------------------------------------------------------------
!  call calQuadraticPoles2(LocalNumAtoms,lmax_kkr,lmax_phi,            &
!                         n_spin_pola,ZeroIntHalfWidth,Efermi)
!  -------------------------------------------------------------------
   write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
!
!  -------------------------------------------------------------------
!  call calGoldenPoles(LocalNumAtoms,lmax_kkr,n_spin_pola,            &
!                      ZeroIntHalfWidth,Efermi)
!  -------------------------------------------------------------------
!
   t1 = getTime()
!  -------------------------------------------------------------------
   call calMinPoles(LocalNumAtoms,lmax_kkr,n_spin_pola,               &
                    ErBottom,-ZeroIntHalfWidth)
!  -------------------------------------------------------------------
   call calMinPoles(LocalNumAtoms,lmax_kkr,n_spin_pola,               &
                    ZeroIntHalfWidth,Efermi)
!  -------------------------------------------------------------------
   write(6,'(/,a,f12.5,a,/)') 'Minimization technique time: ',getTime() - t1,' sec.'
!
!  ===================================================================
!  clean up the allocated spaces.
!  ===================================================================
   deallocate( atom_print_level )
   deallocate( AtomPosition,AtomicNumber,lmax_kkr,lmax_phi,lmax_rho,lmax_pot )
   deallocate( lmax_step, GlobalIndex, Poles, evec )
!
   if (node_print_level >= 0) then
      write(6,*) " Run Ending! "
      call FlushFile(6)
   endif
!
   call endSSSolver()
   call endSpinRotation()
   call endPotential()
   call endStepFunction()
   call endRadialGrid()
   call endSystemSymmetry()
   call endMadelung()
   call endGauntFactors()
   call endSphericalHarmonics()
   call endPotentialType()
   call endSystemVolume()
   call endAtom()
   call endAtom2Proc()
   call endParallelIO()
   call endSystem()
   call endDataServiceCenter()
   call endGroupComm()
   call endMPP()
!
!  ===================================================================
   if (node_print_level >= 0) then
      call date_and_time(exec_date,exec_time)
      write(6,'(/,12a)')'Execution ends at ',                         &
           exec_time(1:2),':',exec_time(3:4),':',exec_time(5:6),', ', &
           exec_date(5:6),'-',exec_date(7:8),'-',exec_date(1:4)
      write(6,'(80(''-''))')
   endif
!
   call endOutput()
!
   if ( node_print_level >= 0 ) then
      stop 'End of the program. OK!'
   else
      stop
   endif
end program testSineMatrixPole
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calQuadraticPoles(LocalNumAtoms,lkkr,lphi,              &
                                n_spin_pola,eb,et,ldp,NumPoles,Poles, &
                                isZeroInterval)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                               Ten2m8, FOURTH, CZERO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getSineMatrix, solveSingleScattering
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix, &
                                     solveQuadraticEquation, getEigenValue,   &
                                     solveLinearEquation
!
   use WriteMatrixModule, only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: lphi(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: n_spin_pola,ldp
   integer (kind=IntKind), intent(out) :: NumPoles(LocalNumAtoms,n_spin_pola)
!
   integer (kind=IntKind) :: id, ie, is, iw, n, kmax_kkr, NumWindows, info
   integer (kind=IntKind) :: i, j, l, m, kl, lp, mp, klp, nv
!
   logical, optional :: isZeroInterval
   logical, parameter :: isGeneral = .false.
!
   real (kind=RealKind), parameter :: Delta = 0.002d0
   real (kind=RealKind), parameter :: WindowWidth = 4.05*Delta
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0
   real (kind=RealKind), intent(out) :: Poles(2*ldp,LocalNumAtoms,n_spin_pola)
!
   complex (kind=CmplxKind) :: e, kappa, a2l, a2lp
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), allocatable :: s0(:,:), s1(:,:), s2(:,:)
   complex (kind=CmplxKind), allocatable :: sT(:,:), sTs(:,:)
   complex (kind=CmplxKind), pointer :: pv(:)
!
   if (eb > et) then
      call ErrorHandler('calQuadraticPoles','eb > et',eb,et)
   endif
!
   if(.not.present(isZeroInterval)) isZeroInterval = .false.
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
!
   kmax_kkr = (lkkr(1)+1)**2
   allocate( s0(1:kmax_kkr,1:kmax_kkr) )
   allocate( s1(1:kmax_kkr,1:kmax_kkr) )
   allocate( s2(1:kmax_kkr,1:kmax_kkr) )
   allocate( sT(1:kmax_kkr,1:kmax_kkr) )
   allocate( sTs(1:kmax_kkr,1:kmax_kkr) )
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
            allocate( sT(1:kmax_kkr,1:kmax_kkr) )
            allocate( sTs(1:kmax_kkr,1:kmax_kkr) )
            call endQuadraticMatrix()
            call initQuadraticMatrix(kmax_kkr)
         endif
!
         n = 0
         do iw = 1, NumWindows
            w0 = eb + (iw-1)*WindowWidth
            e0 = w0 + (HALF)*WindowWidth
            if (.not.isZeroInterval.and.(abs(e0) < Ten2m6 .or.        &
                abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6)) then
               write(6,*)'e0 shifted'
               e0 = e0 - HALF*de
            else if(isZeroInterval) then
               e0 = ZERO
            endif
!
            if(isZeroInterval) then
               s0(1:kmax_kkr,1:kmax_kkr) = CZERO
            else
               e = cmplx(e0,ZERO,kind=CmplxKind)
               kappa = sqrt(e)
!              -------------------------------------------------------
               call solveSingleScattering(is,id,e,CZERO)
!              -------------------------------------------------------
               sin_mat => getSineMatrix(is,id)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!              -------------------------------------------------------
!
!              =======================================================
!              rescaling sine matrix
!              =======================================================
!              kl = 0
!              a2l = CONE/kappa
!              do l = 0, lkkr(id)
!                 do m = -l, l
!                    kl = kl + 1
!                    klp = 0
!                    a2lp = CONE
!                    do lp = 0, lkkr(id)
!                       do mp = -lp, lp
!                          klp = klp + 1
!                          s0(klp,kl) = s0(klp,kl)*a2l*a2lp
!                       enddo
!                       a2lp = a2lp*(2*lp+ONE)/kappa
!                    enddo
!                 enddo
!                 a2l = a2l*(2*l+ONE)/kappa
!              enddo
!
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!              -------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                          sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sTs,1,s0,1)
!              -------------------------------------------------------
            endif
!
            e = cmplx(e0+de,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
            sin_mat => getSineMatrix(is,id)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s2,1)
!           ----------------------------------------------------------
!
!           ==========================================================
!           rescaling sine matrix
!           ==========================================================
!           kl = 0
!           a2l = CONE/kappa
!           do l = 0, lkkr(id)
!              do m = -l, l
!                 kl = kl + 1
!                 klp = 0
!                 a2lp = CONE
!                 do lp = 0, lkkr(id)
!                    do mp = -lp, lp
!                       klp = klp + 1
!                       s2(klp,kl) = s2(klp,kl)*a2l*a2lp
!                    enddo
!                    a2lp = a2lp*(2*lp+ONE)/kappa
!                 enddo
!              enddo
!              a2l = a2l*(2*l+ONE)/kappa
!           enddo
!
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,s2,1,sT,1)
!           ----------------------------------------------------------
            call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,       &
                       sT,kmax_kkr,s2,kmax_kkr,CZERO,sTs,kmax_kkr)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sTs,1,s2,1)
!           ----------------------------------------------------------
!
            e = cmplx(e0-de,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
            sin_mat => getSineMatrix(is,id)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s1,1)
!           ----------------------------------------------------------
!
!           ==========================================================
!           rescaling sine matrix
!           ==========================================================
!           kl = 0
!           a2l = CONE/kappa
!           do l = 0, lkkr(id)
!              do m = -l, l
!                 kl = kl + 1
!                 klp = 0
!                 a2lp = CONE
!                 do lp = 0, lkkr(id)
!                    do mp = -lp, lp
!                       klp = klp + 1
!                       s1(klp,kl) = s1(klp,kl)*a2l*a2lp
!                    enddo
!                    a2lp = a2lp*(2*lp+ONE)/kappa
!                 enddo
!              enddo
!              a2l = a2l*(2*l+ONE)/kappa
!           enddo
!
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,s1,1,sT,1)
!           ----------------------------------------------------------
            call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,       &
                       sT,kmax_kkr,s1,kmax_kkr,CZERO,sTs,kmax_kkr)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sTs,1,sin_mat,1)
!           ----------------------------------------------------------
!
            s1 = (s2 - sin_mat)/de2
            s2 = (s2 + sin_mat - TWO*s0)/dede2
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
!           write(6,'(a,e20.12)')'Eigenvalues: e0=',e0
            do ie = 1, kmax_kkr*2
!              write(6,'(4x,2e20.12)')pv(ie)
               if (abs(aimag(pv(ie))) < Ten2m8) then
                  pe = real(pv(ie),kind=RealKind) + e0
                  if (pe >= w0 .and. pe <= w0+WindowWidth              &
                      .and. abs(pe) > Ten2m8) then
                     n = n + 1
                     Poles(n,id,is) = pe
!                    write(6,'(4x,a,e20.12)')'Pole Found:',pe
                  endif
               endif
            enddo
         enddo
         NumPoles(id,is) = n
      enddo
   enddo
!
   nullify(pv)
   call endQuadraticMatrix()
!
!  call writeMatrix('s0',s0,kmax_kkr,kmax_kkr)
   deallocate( sTs )
   deallocate( sT )
   deallocate( s1 )
   write(6,*)'OK2'
   call FlushFile(6)
   deallocate( s2 )
   write(6,*)'OK3'
   call FlushFile(6)
   deallocate( s0 ) !,s1,s2,sT,sTs )
   write(6,*)'OK1'
   call FlushFile(6)
!
   end subroutine calQuadraticPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calQuadraticPoles2(LocalNumAtoms,lkkr,lphi,              &
                                n_spin_pola,eb,et)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                               Ten2m8, CZERO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getSineMatrix,                    &
                                     solveSingleScattering
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix, &
                                     solveQuadraticEquation, getEigenValue,   &
                                     solveLinearEquation
!
   use WriteMatrixModule, only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: lphi(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: n_spin_pola
!
   integer (kind=IntKind) :: id, ie, is, iw, n, kmax_kkr, info, nv
   integer (kind=IntKind), parameter :: max_iter = 5
   integer (kind=IntKind), parameter :: NumWindows = 10
!
   logical, parameter :: isGeneral = .false.
!
   real (kind=RealKind), parameter :: Delta = TEN2m6, max_dist = 10.0d0
   real (kind=RealKind), parameter :: GoldenWidth = 0.01
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind) :: e0, e0new, de, de2, dede2, pe, temp
   real (kind=RealKind) :: WindowWidth, ebtemp, ettemp, sfac2inv
   real (kind=RealKind), allocatable :: ScalingFactors(:)
!
   complex (kind=CmplxKind) :: e, dets0, dets1, dets2
   complex (kind=CmplxKind) :: detderiv, detderiv2
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), allocatable :: s0(:,:), s1(:,:), s2(:,:)
   complex (kind=CmplxKind), allocatable :: sT(:,:), sTs(:,:)
   complex (kind=CmplxKind), pointer :: pv(:)
!
   interface
      subroutine calcDet(mat,kmax,det,diag)
      use KindParamModule, only : IntKind, CmplxKind
!
      integer (kind=IntKind), intent(in) :: kmax
!
      complex (kind=CmplxKind), intent(out) :: det
      complex (kind=CmplxKind), optional :: diag(kmax)
      complex (kind=CmplxKind), intent(in) :: mat(kmax,kmax)
      end subroutine calcDet
   end interface
!
   if (eb > et) then
      call ErrorHandler('calQuadraticPoles2','eb > et',eb,et)
   endif
!
   de = Delta
   de2 = de*TWO; dede2 = de*de*TWO
!
   WindowWidth = (et-eb)/NumWindows
!
   kmax_kkr = (lkkr(1)+1)**2
   allocate( s0(1:kmax_kkr,1:kmax_kkr) )
   allocate( s1(1:kmax_kkr,1:kmax_kkr) )
   allocate( s2(1:kmax_kkr,1:kmax_kkr) )
   allocate( sT(1:kmax_kkr,1:kmax_kkr) )
   allocate( sTs(1:kmax_kkr,1:kmax_kkr) )
   allocate(ScalingFactors(kmax_kkr))
!
   call initQuadraticMatrix(kmax_kkr,isGeneral)
!
   call calScalingFactors(ScalingFactors,kmax_kkr)
!
   write(6,'(a,6a20)')'SineMatDet::','Energy','det(s0)','det(s1)',     &
                      'det(s2)','d(det)/de','d2(det)/de2'
   do is = 1,n_spin_pola
      do id = 1,LocalNumAtoms
         if ((lkkr(id)+1)**2 /= kmax_kkr) then
            deallocate(s0, s1, s2, sT, sTs)
            deallocate(ScalingFactors)
            kmax_kkr = (lkkr(id)+1)**2
            allocate( s0(1:kmax_kkr,1:kmax_kkr) )
            allocate( s1(1:kmax_kkr,1:kmax_kkr) )
            allocate( s2(1:kmax_kkr,1:kmax_kkr) )
            allocate( sT(1:kmax_kkr,1:kmax_kkr) )
            allocate( sTs(1:kmax_kkr,1:kmax_kkr) )
            allocate(ScalingFactors(kmax_kkr))
            call endQuadraticMatrix()
            call initQuadraticMatrix(kmax_kkr)
         endif
!
         do iw = 1, NumWindows
            ebtemp = eb+(iw-1)*WindowWidth
            ettemp = eb+iw*WindowWidth
            e0 = ebtemp+HALF*WindowWidth
            e0new = ettemp
            write(6,'(a,i5,4x,a,f16.12)')'Window #',iw,'e0=',e0
            do n = 1, max_iter
               if (abs(e0) < Ten2m6 .or. abs(e0-de) < Ten2m6 .or.        &
                   abs(e0+de) < Ten2m6) then
                  write(6,*)'e0 shifted'
                  e0 = e0 + HALF*de
               endif
!         
               e = cmplx(e0,ZERO,kind=CmplxKind)
!              -------------------------------------------------------
               call solveSingleScattering(is,id,e,CZERO)
!              -------------------------------------------------------
               sin_mat => getSineMatrix(is,id)
!              -------------------------------------------------------
               call scaleMatrix(sin_mat,ScalingFactors,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!              -------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                          sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sTs,1,s0,1)
!              -------------------------------------------------------
               call calcDet(s0,kmax_kkr,dets0)
!              -------------------------------------------------------
!         
               e = cmplx(e0+de,ZERO,kind=CmplxKind)
!              -------------------------------------------------------
               call solveSingleScattering(is,id,e,CZERO)
!              -------------------------------------------------------
               sin_mat => getSineMatrix(is,id)
!              -------------------------------------------------------
               call scaleMatrix(sin_mat,ScalingFactors,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s2,1)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,s2,1,sT,1)
!              -------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                          sT,kmax_kkr,s2,kmax_kkr,CZERO,sTs,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sTs,1,s2,1)
!              -------------------------------------------------------
               call calcDet(s2,kmax_kkr,dets2)
!              -------------------------------------------------------
!         
               e = cmplx(e0-de,ZERO,kind=CmplxKind)
!              -------------------------------------------------------
               call solveSingleScattering(is,id,e,CZERO)
!              -------------------------------------------------------
               sin_mat => getSineMatrix(is,id)
!              -------------------------------------------------------
               call scaleMatrix(sin_mat,ScalingFactors,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s1,1)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,s1,1,sT,1)
!              -------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                          sT,kmax_kkr,s1,kmax_kkr,CZERO,sTs,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sTs,1,sin_mat,1)
!              -------------------------------------------------------
               call calcDet(sin_mat,kmax_kkr,dets1)
!              -------------------------------------------------------
!
               sfac2inv = ONE/(ScalingFactors(1)*ScalingFactors(1))
               detderiv = (dets2-dets1)/de2
               detderiv = sfac2inv*detderiv
               detderiv2 = (dets2+dets1-TWO*dets0)/(de*de)
               detderiv2 = sfac2inv*detderiv2
               s1 = (s2 - sin_mat)/de2
               s2 = (s2 + sin_mat - TWO*s0)/dede2
!              -------------------------------------------------------
               call calcDet(s1,kmax_kkr,dets1)
!              -------------------------------------------------------
               call calcDet(s2,kmax_kkr,dets2)
!              -------------------------------------------------------
               dets0 = sfac2inv*dets0
               dets1 = sfac2inv*dets1
               dets2 = sfac2inv*dets2
               write(6,'(a,6e20.12)')'SineMatDet::',e0,real(dets0),    &
                                      real(dets1),real(dets2),        &
                                      real(detderiv),real(detderiv2)
!         
!              -------------------------------------------------------
               call solveQuadraticEquation(s0,s1,s2,info)
!              -------------------------------------------------------
!         
               if (info /= 0) then
                  stop 'Error in s0, s1, s2'
               endif
!         
!              -------------------------------------------------------
               pv => getEigenValue(nv)
!              -------------------------------------------------------
               temp = max_dist
!              write(6,*)'Real Eigenvalues:'
               do ie = 1, kmax_kkr*2
                  if (abs(aimag(pv(ie))) < Ten2m8) then
                     pe = real(pv(ie),kind=RealKind) + e0
!                    write(6,'(4x,f16.12)')pe
                     if(abs(pe-e0) < temp .and. pe > ebtemp .and. pe < ettemp) then
                        temp = abs(pe-e0)
                        e0new = pe
                     endif
                  endif
               enddo
               write(6,'(a,f16.12)')'New e0=',e0new
               if(abs(e0-e0new) < TEN2m6) exit
               e0 = e0new
            enddo
            call calGoldenPoles(LocalNumAtoms,lkkr,n_spin_pola,       &
                                e0-GoldenWidth,e0+GoldenWidth)
         enddo
      enddo
   enddo
!
   nullify(pv)
   call endQuadraticMatrix()
!
   deallocate( s0,s1,s2,sT,sTs,ScalingFactors )
!
   end subroutine calQuadraticPoles2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine plotSinDet(eb, et, n, LocalNumAtoms, n_spin_pola,       &
                         lkkr)
!  ===================================================================
   use KindParamModule, only : RealKind, CmplxKind, IntKind
!
   use MathParamModule, only : ZERO, ONE, CONE, CZERO, HALF, TEN2m8
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use MPPModule, only : MyPE
!
   use SSSolverModule, only : getSineMatrix, solveSingleScattering
!
   implicit none
!
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind) :: e_step, e0
!
   complex (kind=CmplxKind) :: e, det, a2l, a2lp, kappa
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), allocatable :: SinU(:,:)
   complex (kind=CmplxKind), allocatable :: sT(:,:), sTs(:,:)
!
   integer (kind=IntKind), intent(in) :: n, LocalNumAtoms, n_spin_pola
   integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
   integer (kind=IntKind) :: i, l, m, kl, lp, mp, klp, id, is
   integer (kind=IntKind) :: kmax_kkr
!
   if(et < eb) then
      call ErrorHandler('plotSinDet','et < eb', et, eb)
   endif
!
   kmax_kkr = (lkkr(1)+1)**2
!
   allocate(SinU(1:kmax_kkr,1:kmax_kkr))
   allocate(sT(1:kmax_kkr,1:kmax_kkr))
   allocate(sTs(1:kmax_kkr,1:kmax_kkr))
!
   if(MyPE == 0) then
      open(unit=12,file='SineMatrixDet.dat',status='unknown',         &
           form='formatted')
   endif
!
   if(MyPE == 0) then
      write(12,'(3a20)')'Energy','Real(Det)'
   endif
!
   e_step = (et-eb)/n
!
   do id = 1, LocalNumAtoms
      if(MyPE == 0) then
         write(12,'(a,i5)')'Atom:',id
      endif
      do is = 1, n_spin_pola
         if(MyPE == 0) then
            write(12,'(a,i5)')'Spin:',is
         endif
!
         do i = 1, n
            e0 = eb+i*e_step
            if(e0 < TEN2m8) then
               e0 = e0+HALF*e_step
            endif
!
            e = cmplx(e0,ZERO,CmplxKind)
            kappa = sqrt(e)
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
            sin_mat => getSineMatrix(is,id)
            do kl = 1, kmax_kkr
               SinU(1:kmax_kkr,kl) = sin_mat(1:kmax_kkr,kl)
            enddo
!           ==========================================================
!           rescaling sine matrix
!           ==========================================================
!           kl = 0
!           a2l = CONE/kappa
!           do l = 0, lkkr(id)
!              do m = -l, l
!                 kl = kl + 1
!                 klp = 0
!                 a2lp = CONE
!                 do lp = 0, lkkr(id)
!                    do mp = -lp, lp
!                       klp = klp + 1
!                       SinU(klp,kl) = SinU(klp,kl)*a2l*a2lp
!                    enddo
!                    a2lp = a2lp*(2*lp+ONE)/kappa
!                 enddo
!              enddo
!              a2l = a2l*(2*l+ONE)/kappa
!           enddo
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,SinU,1,sT,1)
!           ----------------------------------------------------------
            call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,       &
                       sT,kmax_kkr,SinU,kmax_kkr,CZERO,sTs,kmax_kkr)
!           ----------------------------------------------------------
            do kl = 1, kmax_kkr
               SinU(1:kmax_kkr,kl) = sTs(1:kmax_kkr,kl)
            enddo
!           ----------------------------------------------------------
            call GaussianElim(SinU,kmax_kkr)
!           ----------------------------------------------------------
            det = CONE
            do kl = 1,kmax_kkr
               det = det*SinU(kl,kl)
            enddo
!
            write(12,'(2e20.12)')e0,real(det,RealKind)
!
         enddo
      enddo
   enddo
!
   close(12)
!
   deallocate(SinU, sT, sTs)
!
   end subroutine plotSinDet
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calScalingFactors(ScalingFactors,kmax)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, ONE, TEN2m8, TEN2m6, CZERO
!
   use SSSolverModule, only : getSineMatrix, solveSingleScattering
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kmax
   integer (kind=IntKind) :: i, j
!
   real (kind=RealKind), intent(out) :: ScalingFactors(kmax)
   real (kind=RealKind), parameter :: energy = 0.2
   real (kind=RealKind) :: row_min
!
   complex (kind=CmplxKind) :: e
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
!
   e = cmplx(energy, ZERO, CmplxKind)
!
!  -------------------------------------------------------------------
   call solveSingleScattering(1,1,e,CZERO)
!  -------------------------------------------------------------------
   sin_mat => getSineMatrix(1,1)
!  ===================================================================
!  Find the non-zero row minimum of the magnitudes of the elements of &
!  sin_mat and place inverse in ScalingFactors
!  ===================================================================
   write(6,*)'calScalingFactors ::'
   write(6,'(8x,2a5,a20)')'i','j','abs(sin_mat(i,j))'
   do i = 1, kmax
!     row_min = ONE
!     do j = 1, kmax
!        if(abs(sin_mat(i,j)) /= ZERO) then
!           write(6,'(8x,2i5,e20.12)')i,j,abs(sin_mat(i,j))
!        endif
!        if(abs(sin_mat(i,j)) > TEN2m8) then
!           row_min = min(row_min, abs(sin_mat(i,j)))
!        endif
!     enddo
!     ScalingFactors(i) = ONE/row_min
      ScalingFactors(i) = ONE
      write(6,'(4x,a,i3,a,e20.12)')'Scaling Factor #',i,'=',          &
                                   ScalingFactors(i)
   enddo
!
   end subroutine calScalingFactors
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine scaleMatrix(mat,ScalingFactors,kmax)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kmax
   integer (kind=IntKind) :: i, j
!
   real (kind=RealKind), intent(in) :: ScalingFactors(kmax)
!
   complex (kind=CmplxKind), intent(inout) :: mat(kmax, kmax)
!
   do i = 1, kmax
      do j = 1, kmax
         mat(i,j) = ScalingFactors(i)*mat(i,j)
      enddo
   enddo
!
   end subroutine scaleMatrix
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
   complex (kind=CmplxKind), optional :: diag(kmax)
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
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calGoldenPoles(LocalNumAtoms,lkkr,n_spin_pola,eb,et)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                               Ten2m8, CZERO, FIVE
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getSineMatrix,                    &
                                     solveSingleScattering
!
   use WriteMatrixModule, only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: n_spin_pola
!
   integer (kind=IntKind) :: id, ie, is, iw, n, kmax_kkr, info, kl
   integer (kind=IntKind) :: kbracket, l, m
   integer (kind=IntKind), parameter :: max_iter = 20
   integer (kind=IntKind), parameter :: NumWindows = 20
!
   logical, parameter :: isGeneral = .false.
!
   real (kind=RealKind), parameter :: Delta = TEN2m6, max_dist = 10.0d0
   real (kind=RealKind), parameter :: phi = (ONE+sqrt(FIVE))/TWO
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind) :: WindowWidth, sfac2inv, e0, e1, e2, e3
   real (kind=RealKind), allocatable :: ScalingFactors(:)
!
   complex (kind=CmplxKind) :: e, dets0, dets1, dets2, dets3
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), allocatable :: s0diag(:), s1diag(:)
   complex (kind=CmplxKind), allocatable :: s2diag(:), s3diag(:)
   complex (kind=CmplxKind), allocatable :: s0(:,:)
   complex (kind=CmplxKind), allocatable :: sT(:,:), sTs(:,:)
!
   interface
      subroutine calcDet(mat,kmax,det,diag)
      use KindParamModule, only : IntKind, CmplxKind
!
      integer (kind=IntKind), intent(in) :: kmax
!
      complex (kind=CmplxKind), intent(out) :: det
      complex (kind=CmplxKind), optional :: diag(kmax)
      complex (kind=CmplxKind), intent(in) :: mat(kmax,kmax)
      end subroutine calcDet
   end interface
!
   if (eb > et) then
      call ErrorHandler('calGoldenPoles','eb > et',eb,et)
   endif
!
   WindowWidth = (et-eb)/NumWindows
!
   kmax_kkr = (lkkr(1)+1)**2
   allocate( s0(1:kmax_kkr,1:kmax_kkr) )
   allocate( sT(1:kmax_kkr,1:kmax_kkr) )
   allocate( sTs(1:kmax_kkr,1:kmax_kkr) )
   allocate(ScalingFactors(kmax_kkr))
   allocate(s0diag(kmax_kkr),s1diag(kmax_kkr))
   allocate(s2diag(kmax_kkr),s3diag(kmax_kkr))
!
   call calScalingFactors(ScalingFactors,kmax_kkr)
   sfac2inv = ONE/(ScalingFactors(1)*ScalingFactors(1))
!
   write(6,*)'calGoldenPoles::'
   do is = 1,n_spin_pola
      do id = 1,LocalNumAtoms
         if ((lkkr(id)+1)**2 /= kmax_kkr) then
            deallocate(s0, sT, sTs)
            deallocate(s0diag, s1diag, s2diag, s3diag)
            kmax_kkr = (lkkr(id)+1)**2
            allocate( s0(1:kmax_kkr,1:kmax_kkr) )
            allocate( sT(1:kmax_kkr,1:kmax_kkr) )
            allocate( sTs(1:kmax_kkr,1:kmax_kkr) )
            allocate(s0diag(kmax_kkr),s1diag(kmax_kkr))
            allocate(s2diag(kmax_kkr),s3diag(kmax_kkr))
         endif
!
         do iw = 1, NumWindows
            e0 = eb+(iw-1)*WindowWidth
            e2 = eb+iw*WindowWidth
            e1 = (e2+phi*e0)/(ONE+phi)
            e3 = e0+e2-e1
            write(6,'(4x,a,i3,4x,a,f16.12)')'Window #',iw,'e0=',e0
            do n = 1, max_iter
               write(6,'(8x,a,i3,4(4x,a,f16.12))')'iter',n,'e0=',e0,  &
                                                  'e1=',e1,'e2=',e2,  &
                                                  'e3=',e3
!         
               e = cmplx(e0,ZERO,kind=CmplxKind)
!              -------------------------------------------------------
               call solveSingleScattering(is,id,e,CZERO)
!              -------------------------------------------------------
               sin_mat => getSineMatrix(is,id)
!              -------------------------------------------------------
               call scaleMatrix(sin_mat,ScalingFactors,kmax_kkr)
!              call writeMatrix('sin_mat0',sin_mat,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!              -------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                          sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!              call writeMatrix('sTs0',sTs,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call calcDet(sTs,kmax_kkr,dets0,s0diag)
!              -------------------------------------------------------
!
               e = cmplx(e1,ZERO,kind=CmplxKind)
!              -------------------------------------------------------
               call solveSingleScattering(is,id,e,CZERO)
!              -------------------------------------------------------
               sin_mat => getSineMatrix(is,id)
!              -------------------------------------------------------
               call scaleMatrix(sin_mat,ScalingFactors,kmax_kkr)
!              call writeMatrix('sin_mat1',sin_mat,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!              -------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                          sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!              call writeMatrix('sTs1',sTs,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call calcDet(sTs,kmax_kkr,dets1,s1diag)
!              -------------------------------------------------------
!
               e = cmplx(e2,ZERO,kind=CmplxKind)
!              -------------------------------------------------------
               call solveSingleScattering(is,id,e,CZERO)
!              -------------------------------------------------------
               sin_mat => getSineMatrix(is,id)
!              -------------------------------------------------------
               call scaleMatrix(sin_mat,ScalingFactors,kmax_kkr)
!              call writeMatrix('sin_mat2',sin_mat,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!              -------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                          sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!              call writeMatrix('sTs2',sTs,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call calcDet(sTs,kmax_kkr,dets2,s2diag)
!              -------------------------------------------------------
!              write(6,'(8x,3(4x,a,f16.12))')'dets0=',real(dets0),    &
!                                            'dets1=',real(dets1),    &
!                                            'dets2=',real(dets2)
               l = 0
               m = 0
               do kl = 1, kmax_kkr
                  write(6,'(a,i3,2e20.12)')'SineMatDiag:: l=',l,e0,    &
                                        real(s0diag(kl))      
                  write(6,'(a,i3,2e20.12)')'SineMatDiag:: l=',l,e1,    &
                                        real(s1diag(kl))      
                  write(6,'(a,i3,2e20.12)')'SineMatDiag:: l=',l,e2,    &
                                        real(s2diag(kl))
                  m = m+1
                  if(m == 2*l+1) then
                     l = l+1
                     m = 0
                  endif
               enddo
               do kl = 1, kmax_kkr
                  write(6,'(12x,a,i3,3(4x,a,e20.12))')'kl=',kl,       &
                                                's0diag=',            &
                                                real(s0diag(kl)),     &
                                                's1diag=',            &
                                                real(s1diag(kl)),     &
                                                's2diag=',            &
                                                real(s2diag(kl))
                  if(real(s1diag(kl)) < real(s0diag(kl)) .and.        &
                     real(s1diag(kl)) < real(s2diag(kl))) then
                     kbracket = kl
                     write(6,'(a,2e20.12)')'SineMatBrak::',e0,         &
                                           real(s0diag(kl))
                     write(6,'(a,2e20.12)')'SineMatBrak::',e1,         &
                                           real(s1diag(kl))
                     write(6,'(a,2e20.12)')'SineMatBrak::',e2,         &
                                           real(s2diag(kl))
                     exit
                  else
                     kbracket = 0
                  endif
               enddo
               if(kbracket == 0) exit
               write(6,'(8x,a,i3)')'kbracket=',kbracket
!
               e = cmplx(e3,ZERO,kind=CmplxKind)
!              -------------------------------------------------------
               call solveSingleScattering(is,id,e,CZERO)
!              -------------------------------------------------------
               sin_mat => getSineMatrix(is,id)
!              -------------------------------------------------------
               call scaleMatrix(sin_mat,ScalingFactors,kmax_kkr)
!              call writeMatrix('sin_mat3',sin_mat,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!              -------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,    &
                          sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!              call writeMatrix('sTs3',sTs,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call calcDet(sTs,kmax_kkr,dets3,s3diag)
!              -------------------------------------------------------
!              write(6,'(12x,a,f16.12)')'dets3=',real(dets3)
!
               write(6,'(8x,3(4x,a,e20.12))')'s3diag=',               &
                                             real(s3diag(kbracket))
               write(6,'(a,2e20.12)')'SineMatBrak::',e3,               &
                                     real(s3diag(kbracket))
               if(real(s3diag(kbracket)) < real(s1diag(kbracket))) then
                  if(e2-e1 > e1-e0) then
                     e0 = e1
                     e1 = (e2+phi*e0)/(ONE+phi)
                  else
                     e2 = e1
                     e1 = (e0+phi*e2)/(ONE+phi)
                  endif
               else
                  if(e2-e1 > e1-e0) then
                     e2 = e3
                     e1 = (e2+phi*e0)/(ONE+phi)
                  else
                     e0 = e3
                     e1 = (e0+phi*e2)/(ONE+phi)
                  endif
               endif
               e3 = e0+e2-e1
!
            enddo
         enddo
      enddo
   enddo
!
   deallocate( s0,sT,sTs,ScalingFactors )
   deallocate(s0diag, s1diag, s2diag, s3diag)
!
   end subroutine calGoldenPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calMinPoles(LocalNumAtoms,lkkr,n_spin_pola,ebot,etop)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, ONE, TEN2m6
!
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms, n_spin_pola
   integer (kind=IntKind), intent(in) :: lkkr(LocalNumAtoms)
   integer (kind=IntKind) :: id, is, ic, kl
   integer (kind=IntKind) :: kmax_kkr, numMax, numMin
   integer (kind=IntKind), allocatable :: MaxChannels(:)
   integer (kind=IntKind), allocatable :: MinChannels(:)
!
   real (kind=RealKind), intent(in) :: ebot, etop
   real (kind=RealKind), parameter :: eStep = 0.01
   real (kind=RealKind) :: eb, et
   real (kind=RealKind), allocatable :: MaxValues(:), MinValues(:)
!
   if (ebot > etop) then
      call ErrorHandler('calMinPoles','eb > et',ebot,etop)
   endif
!
   et = etop
   eb = ebot
   if(eb >= ZERO .and. et >= ZERO) then
      if(eb < TEN2m6) eb = TEN2m6
      if(et < TEN2m6) et = TEN2m6+eStep
   else if(eb < ZERO .and. et < ZERO) then
      if(et > -TEN2m6) et = -TEN2m6
      if(eb > -TEN2m6) eb = -TEN2m6-eStep
   else
      call ErrorHandler('calManPoles','eb and et opposite sign',  &
                        eb,et)
   endif
!
   kmax_kkr = (lkkr(1)+1)**2
!
   allocate(MaxChannels(kmax_kkr), MaxValues(kmax_kkr))
   allocate(MinChannels(kmax_kkr), MinValues(kmax_kkr))
!
   do kl = 1, kmax_kkr
      MaxChannels(kl) = 0
      MinChannels(kl) = 0
   enddo
!
   do id = 1, LocalNumAtoms
      do is = 1, n_spin_pola
         if((lkkr(id)+1)**2 /= kmax_kkr) then
            deallocate(MaxChannels, MaxValues, MinChannels, MinValues)
            kmax_kkr = (lkkr(id)+1)**2
            allocate(MaxChannels(kmax_kkr), MaxValues(kmax_kkr))
            allocate(MinChannels(kmax_kkr), MinValues(kmax_kkr))
         endif
!
         write(6,'(a,e20.12,a,e20.12,a)')'calMinPoles:: Energy range = (', &
                                       eb,',',et,')'
         write(6,'(2(4x,a,i3))')'Atom:',id,'Spin:',is
!        -------------------------------------------------------------
         call findMaxChannels(id, is, eb, et, kmax_kkr, MaxChannels,  &
                              MaxValues, numMax)
!        -------------------------------------------------------------
         call findGoldenMin(id, is, eb, et, kmax_kkr, MaxChannels,    &
                            MaxValues, numMax, MinChannels,           &
                            MinValues, numMin)
!        -------------------------------------------------------------
         call findNewtonMin(id, is, eb, et, kmax_kkr, MinChannels,    &
                            MinValues, numMin, numMax)
!        -------------------------------------------------------------
!
         write(6,'(8x,a,i4)')'Num Poles Found:',numMin
         do ic = 1, numMax
            kl = MinChannels(ic)
            if(kl == 0) cycle
            write(6,'(12x,a,i4,4x,a,f16.12)')'kl=',kl,'e=',           &
                                             MinValues(ic)
         enddo
      enddo
   enddo
!
   deallocate(MaxChannels, MaxValues, MinChannels, MinValues)
!
   end subroutine calMinPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine findMaxChannels(id, is, eb, et, kmax_kkr, MaxChannels,  &
                              MaxValues, numMax)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, ONE, TWO, FIVE, HALF, CZERO, CONE
   use MathParamModule, only : TEN2m6
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getSineMatrix, solveSingleScattering
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, kmax_kkr
   integer (kind=IntKind), intent(out) :: MaxChannels(kmax_kkr)
   integer (kind=IntKind), intent(out) :: numMax
   integer (kind=IntKind), parameter :: maxGold = 5
   integer (kind=IntKind) :: numEs
   integer (kind=IntKind) :: ie, ic, ig, kl
!
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind), intent(out) :: MaxValues(kmax_kkr)
   real (kind=RealKind), parameter :: eStep = 0.01
   real (kind=RealKind), parameter :: phi = (ONE+sqrt(FIVE))/TWO
   real (kind=RealKind) :: e1, e2, e3, e4, nFac
!
   complex (kind=CmplxKind) :: e, dets1, dets2, dets3, dets4
   complex (kind=CmplxKind) :: s1diag(kmax_kkr), s2diag(kmax_kkr)
   complex (kind=CmplxKind) :: s3diag(kmax_kkr), s4diag(kmax_kkr)
   complex (kind=CmplxKind) :: s0(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind) :: sT(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind) :: sTs(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
!
   interface
      subroutine calcDet(mat,kmax,det,diag)
      use KindParamModule, only : IntKind, CmplxKind
!
      integer (kind=IntKind), intent(in) :: kmax
!
      complex (kind=CmplxKind), intent(out) :: det
      complex (kind=CmplxKind), optional :: diag(kmax)
      complex (kind=CmplxKind), intent(in) :: mat(kmax,kmax)
      end subroutine calcDet
   end interface
   !
   if (eb > et) then
      call ErrorHandler('findMaxChannels','eb > et',eb,et)
   endif
!
   if(eb >= ZERO .and. et >= ZERO) then
      nFac = ONE
   else if(eb < ZERO .and. et < ZERO) then
      nFac = -ONE
   else
      call ErrorHandler('findMaxChannels','eb and et opposite sign',  &
                        eb,et)
   endif
!
   numEs = floor((et-eb)/eStep, IntKind)
!
!  ==================================================================
!  Find any channels of the SineT*Sine Matrix that have a maximum
!  ==================================================================
   ic = 0
   if(nFac == ONE) then
      e1 = eb
      e2 = eb+eStep
      e3 = eb+TWO*eStep
   else
      e1 = et
      e2 = et-eStep
      e3 = et-TWO*eStep
   endif
!
   e = cmplx(e1,ZERO,kind=CmplxKind)
!  ----------------------------------------------------------------
   call solveSingleScattering(is,id,e,CZERO)
!  ----------------------------------------------------------------
   sin_mat => getSineMatrix(is,id)                                  
!  ----------------------------------------------------------------
   call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)                    
!  ----------------------------------------------------------------
   call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)                         
!  ----------------------------------------------------------------
   call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
              sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)          
!  ----------------------------------------------------------------
   call calcDet(nFac*sTs,kmax_kkr,dets1,s1diag)                         
!  ----------------------------------------------------------------
!                                                                  
   e = cmplx(e2,ZERO,kind=CmplxKind)                               
!  ----------------------------------------------------------------
   call solveSingleScattering(is,id,e,CZERO)
!  ----------------------------------------------------------------
   sin_mat => getSineMatrix(is,id)                                  
!  ----------------------------------------------------------------
   call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)                    
!  ----------------------------------------------------------------
   call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)                         
!  ----------------------------------------------------------------
   call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
              sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)          
!  ----------------------------------------------------------------
   call calcDet(nFac*sTs,kmax_kkr,dets2,s2diag)                         
!  ----------------------------------------------------------------
!                                                                  
   e = cmplx(e3,ZERO,kind=CmplxKind)                               
!  ----------------------------------------------------------------
   call solveSingleScattering(is,id,e,CZERO)                             
!  ----------------------------------------------------------------
   sin_mat => getSineMatrix(is,id)                                  
!  ----------------------------------------------------------------
   call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)                    
!  ----------------------------------------------------------------
   call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)                         
!  ----------------------------------------------------------------
   call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
              sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)          
!  ----------------------------------------------------------------
   call calcDet(nFac*sTs,kmax_kkr,dets3,s3diag)                         
!  ----------------------------------------------------------------
!
   do kl = 1, kmax_kkr
      write(6,'(12x,a,i4,3(4x,a,e20.12))')'kl=',kl,'s1diag=',       &
                                       real(s1diag(kl)),'s2diag=', &
                                       real(s2diag(kl)),'s3diag=', &
                                       real(s3diag(kl))
      if(real(s2diag(kl)) > real(s1diag(kl)) .and.                 &
         real(s2diag(kl)) > real(s3diag(kl))) then
         ic = ic+1
!
         if(ic > kmax_kkr) then
            call ErrorHandler('findMaxChannels', 'ic > kmax_kkr',  &
                              ic, kmax_kkr)
         endif
!
         MaxChannels(ic) = kl
         MaxValues(ic) = e2
         write(6,'(4x,a,i4,4x,a,f16.12)')'Max Bracketed: kl=',kl,   &
                                        'e2=',e2
      endif
   enddo
   do ie = 3, numEs
      e1 = e2
      e2 = e3
      e3 = e2+nFac*eStep
      write(6,'(4x,3(4x,a,f16.12))')'e1=',e1,'e2=',e2,'e3=',e3
!
!     ----------------------------------------------------------------
      call zcopy(kmax_kkr,s2diag,1,s1diag,1)                    
!     ----------------------------------------------------------------
      call zcopy(kmax_kkr,s3diag,1,s2diag,1)                         
!     ----------------------------------------------------------------
!                                                                     
      e = cmplx(e3,ZERO,kind=CmplxKind)                               
!     ----------------------------------------------------------------
      call solveSingleScattering(is,id,e,CZERO)                             
!     ----------------------------------------------------------------
      sin_mat => getSineMatrix(is,id)                                  
!     ----------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)                    
!     ----------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)                         
!     ----------------------------------------------------------------
      call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
                 sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)          
!     ----------------------------------------------------------------
      call calcDet(nFac*sTs,kmax_kkr,dets3,s3diag)                         
!     ----------------------------------------------------------------
!
      do kl = 1, kmax_kkr
         write(6,'(12x,a,i4,3(4x,a,e20.12))')'kl=',kl,'s1diag=',       &
                                          real(s1diag(kl)),'s2diag=', &
                                          real(s2diag(kl)),'s3diag=', &
                                          real(s3diag(kl))
         if(real(s2diag(kl)) > real(s1diag(kl)) .and.                 &
            real(s2diag(kl)) > real(s3diag(kl))) then
            ic = ic+1
!
            if(ic > kmax_kkr) then
               call ErrorHandler('findMaxChannels', 'ic > kmax_kkr',  &
                                 ic, kmax_kkr)
            endif
!
            MaxChannels(ic) = kl
            MaxValues(ic) = e2
            write(6,'(4x,a,i4,4x,a,f16.12)')'Max Bracketed: kl=',kl,   &
                                           'e2=',e2
         endif
      enddo
   enddo
   numMax = ic
!
!  ==================================================================
!  For each bracketed maximum, do a golden maximization
!  ==================================================================
   do ic = 1, numMax
      kl = MaxChannels(ic)
      e2 = MaxValues(ic)
      e1 = (e2*(ONE+phi)-nFac*TWO*eStep)/(ONE+phi)
      e3 = e2+(e2+(nFac*TWO*eStep-e2*(ONE+phi))/(ONE+phi))*phi
      e4 = e1+e3-e2
      do ig = 1, maxGold
         e = cmplx(e1,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets1,s1diag)
!        ------------------------------------------------------------
!      
         e = cmplx(e2,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets2,s2diag)
!        ------------------------------------------------------------
!      
         e = cmplx(e3,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets3,s3diag)
!        ------------------------------------------------------------
!
         if(real(s2diag(kl)) <= real(s1diag(kl)) .or.                 &
            real(s2diag(kl)) <= real(s3diag(kl))) exit
!
         e = cmplx(e4,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets4,s4diag)
!        ------------------------------------------------------------
!
         if(real(s4diag(kl)) > real(s2diag(kl))) then
            if(nFac*(e3-e2) > nFac*(e2-e1)) then
               e1 = e2
               e2 = (e3+phi*e1)/(ONE+phi)
            else
               e3 = e2
               e2 = (e1+phi*e3)/(ONE+phi)
            endif
         else
            if(nFac*(e3-e2) > nFac*(e2-e1)) then
               e3 = e4
               e2 = (e3+phi*e1)/(ONE+phi)
            else
               e1 = e4
               e2 = (e1+phi*e3)/(ONE+phi)
            endif
         endif
         write(6,'(8x,a,i4,4x,a,f16.12)')'kl=',kl,'new e2:',e2
         e4 = e1+e3-e2
      enddo
      MaxValues(ic) = e2
   enddo
!
   end subroutine findMaxChannels
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine findGoldenMin(id, is, eb, et, kmax_kkr, MaxChannels,    &
                            MaxValues, numMax, MinChannels,           &
                            MinValues, numMin)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, ONE, TWO, FIVE, CZERO, CONE
   use MathParamModule, only : TEN2m6
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getSineMatrix, solveSingleScattering
!
   implicit none
!
   logical :: foundMin, e3Switch
!
   integer (kind=IntKind), intent(in) :: id, is, kmax_kkr
   integer (kind=IntKind), intent(in) :: MaxChannels(kmax_kkr)
   integer (kind=IntKind), intent(in) :: numMax
   integer (kind=IntKind), intent(out) :: MinChannels(kmax_kkr)
   integer (kind=IntKind), intent(out) :: numMin
   integer (kind=IntKind), parameter :: maxGold = 10
   integer (kind=IntKind) :: numEs
   integer (kind=IntKind) :: ie, ic, ig, kl
!
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind), intent(in) :: MaxValues(kmax_kkr)
   real (kind=RealKind), intent(out) :: MinValues(kmax_kkr)
   real (kind=RealKind), parameter :: eStep = 0.02
   real (kind=RealKind), parameter :: phi = (ONE+sqrt(FIVE))/TWO
   real (kind=RealKind), parameter :: phiInv = ONE/phi
   real (kind=RealKind) :: e1, e2, e3, e4, ebt, nFac
!
   complex (kind=CmplxKind) :: e, dets1, dets2, dets3, dets4
   complex (kind=CmplxKind) :: s1diag(kmax_kkr), s2diag(kmax_kkr)
   complex (kind=CmplxKind) :: s3diag(kmax_kkr), s4diag(kmax_kkr)
   complex (kind=CmplxKind) :: s0(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind) :: sT(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind) :: sTs(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
!
   interface
      subroutine calcDet(mat,kmax,det,diag)
      use KindParamModule, only : IntKind, CmplxKind
!
      integer (kind=IntKind), intent(in) :: kmax
!
      complex (kind=CmplxKind), intent(out) :: det
      complex (kind=CmplxKind), optional :: diag(kmax)
      complex (kind=CmplxKind), intent(in) :: mat(kmax,kmax)
      end subroutine calcDet
   end interface
!
   if (eb > et) then
      call ErrorHandler('findGoldenMin','eb > et',eb,et)
   endif
!
   if(eb >= ZERO .and. et >= ZERO) then
      nFac = ONE
   else if(eb < ZERO .and. et < ZERO) then
      nFac = -ONE
   else
      call ErrorHandler('findGoldenMin','eb and et opposite sign',  &
                        eb,et)
   endif
!  =====================================================================
!  For each channel with a maximum, do a golden minimization
!  =====================================================================
   numMin = 0
   do ic = 1, numMax
      kl = MaxChannels(ic)
      ebt = MaxValues(ic)
      numEs = floor((et-ebt)/eStep, IntKind)
!     ==================================================================
!     Bracket the minimum of the SineT*Sine Matrix
!     ==================================================================
      foundMin = .false.
      e1 = ebt
      e3 = ebt+nFac*TWO*eStep
      e2 = (e3+phi*e1)/(ONE+phi)
      e3Switch = .false.
!     
      e = cmplx(e1,ZERO,kind=CmplxKind)
!     -------------------------------------------------------------
      call solveSingleScattering(is,id,e,CZERO)
!     -------------------------------------------------------------
      sin_mat => getSineMatrix(is,id)
!     -------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!     -------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!     -------------------------------------------------------------
      call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,          &
                 sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!     -------------------------------------------------------------
      call calcDet(nFac*sTs,kmax_kkr,dets1,s1diag)
!     -------------------------------------------------------------
!     
      e = cmplx(e2,ZERO,kind=CmplxKind)
!     -------------------------------------------------------------
      call solveSingleScattering(is,id,e,CZERO)
!     -------------------------------------------------------------
      sin_mat => getSineMatrix(is,id)
!     -------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!     -------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!     -------------------------------------------------------------
      call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,          &
                 sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!     -------------------------------------------------------------
      call calcDet(nFac*sTs,kmax_kkr,dets2,s2diag)
!     -------------------------------------------------------------
!     
      e = cmplx(e3,ZERO,kind=CmplxKind)
!     -------------------------------------------------------------
      call solveSingleScattering(is,id,e,CZERO)
!     -------------------------------------------------------------
      sin_mat => getSineMatrix(is,id)
!     -------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!     -------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!     -------------------------------------------------------------
      call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,          &
                 sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!     -------------------------------------------------------------
      call calcDet(nFac*sTs,kmax_kkr,dets3,s3diag)
!     -------------------------------------------------------------
!
      write(6,'(12x,a,i4,3(4x,a,e20.12))')'kl=',kl,'s1diag=',       &
                                       real(s1diag(kl)),'s2diag=', &
                                       real(s2diag(kl)),'s3diag=', &
                                       real(s3diag(kl))
      if(real(s2diag(kl)) < real(s1diag(kl)) .and.                 &
         real(s2diag(kl)) < real(s3diag(kl))) then
         foundMin = .true.
         numMin = numMin+1
         MinChannels(ic) = kl
         MinValues(ic) = e2
         write(6,'(4x,a,i4,4x,a,f16.12)')'Min Bracketed: kl=',kl,  &
                                        'e2=',e2
         e4 = e1+e3-e2
!
         e = cmplx(e4,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets4,s4diag)
!        ------------------------------------------------------------
!
         if(real(s4diag(kl)) < real(s2diag(kl))) then
            if(nFac*(e3-e2) > nFac*(e2-e1)) then
               e1 = e2
               e2 = (e3+phi*e1)/(ONE+phi)
            else
               e3 = e2
               e2 = (e1+phi*e3)/(ONE+phi)
            endif
         else
            if(nFac*(e3-e2) > nFac*(e2-e1)) then
               e3 = e4
               e2 = (e3+phi*e1)/(ONE+phi)
            else
               e1 = e4
               e2 = (e1+phi*e3)/(ONE+phi)
            endif
         endif
         write(6,'(8x,a,i4,4x,a,f16.12)')'kl=',kl,'new e2:',e2
      endif
      do ie = 3, numEs
         e1 = e2
         e2 = e3
         if(e3Switch) then
            e3 = e2+(e2+(nFac*TWO*eStep-e2*(ONE+phi))/(ONE+phi))*phi
         else
            e3 = e2+(e2+(nFac*TWO*eStep-e2*(ONE+phiInv))/             &
                                           (ONE+phiInv))*phiInv
         endif
         e3Switch = .not.e3Switch
!     
!        -------------------------------------------------------------
         call zcopy(kmax_kkr,s2diag,1,s1diag,1)
!        -------------------------------------------------------------
         call zcopy(kmax_kkr,s3diag,1,s2diag,1)
!        -------------------------------------------------------------
!     
         e = cmplx(e3,ZERO,kind=CmplxKind)
!        -------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        -------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        -------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        -------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        -------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,          &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        -------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets3,s3diag)
!        -------------------------------------------------------------
!
         write(6,'(12x,a,i4,3(4x,a,e20.12))')'kl=',kl,'s1diag=',       &
                                          real(s1diag(kl)),'s2diag=', &
                                          real(s2diag(kl)),'s3diag=', &
                                          real(s3diag(kl))
         if(real(s2diag(kl)) < real(s1diag(kl)) .and.                 &
            real(s2diag(kl)) < real(s3diag(kl))) then
            foundMin = .true.
            numMin = numMin+1
            MinChannels(ic) = kl
            MinValues(ic) = e2
            write(6,'(4x,a,i4,4x,a,f16.12)')'Min Bracketed: kl=',kl,  &
                                           'e2=',e2
            e4 = e1+e3-e2
!
            e = cmplx(e4,ZERO,kind=CmplxKind)
!           ------------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ------------------------------------------------------------
            sin_mat => getSineMatrix(is,id)
!           ------------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!           ------------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!           ------------------------------------------------------------
            call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                       sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!           ------------------------------------------------------------
            call calcDet(nFac*sTs,kmax_kkr,dets4,s4diag)
!           ------------------------------------------------------------
!
            if(real(s4diag(kl)) < real(s2diag(kl))) then
               if(nFac*(e3-e2) > nFac*(e2-e1)) then
                  e1 = e2
                  e2 = (e3+phi*e1)/(ONE+phi)
               else
                  e3 = e2
                  e2 = (e1+phi*e3)/(ONE+phi)
               endif
            else
               if(nFac*(e3-e2) > nFac*(e2-e1)) then
                  e3 = e4
                  e2 = (e3+phi*e1)/(ONE+phi)
               else
                  e1 = e4
                  e2 = (e1+phi*e3)/(ONE+phi)
               endif
            endif
            write(6,'(8x,a,i4,4x,a,f16.12)')'kl=',kl,'new e2:',e2
            e4 = e1+e3-e2
            exit
         endif
      enddo
      if(.not.foundMin) cycle
      do ig = 1, maxGold
         e = cmplx(e1,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets1,s1diag)
!        ------------------------------------------------------------
!      
         e = cmplx(e2,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets2,s2diag)
!        ------------------------------------------------------------
!      
         e = cmplx(e3,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets3,s3diag)
!        ------------------------------------------------------------
!
         if(real(s2diag(kl)) >= real(s1diag(kl)) .or.                 &
            real(s2diag(kl)) >= real(s3diag(kl))) exit
!
         e = cmplx(e4,ZERO,kind=CmplxKind)
!        ------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,         &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ------------------------------------------------------------
         call calcDet(nFac*sTs,kmax_kkr,dets4,s4diag)
!        ------------------------------------------------------------
!
         if(real(s4diag(kl)) < real(s2diag(kl))) then
            if(nFac*(e3-e2) > nFac*(e2-e1)) then
               e1 = e2
               e2 = (e3+phi*e1)/(ONE+phi)
            else
               e3 = e2
               e2 = (e1+phi*e3)/(ONE+phi)
            endif
         else
            if(nFac*(e3-e2) > nFac*(e2-e1)) then
               e3 = e4
               e2 = (e3+phi*e1)/(ONE+phi)
            else
               e1 = e4
               e2 = (e1+phi*e3)/(ONE+phi)
            endif
         endif
         write(6,'(8x,a,i4,4x,a,f16.12)')'kl=',kl,'new e2:',e2
         e4 = e1+e3-e2
      enddo
      MinValues(ic) = e2
   enddo
!
   end subroutine findGoldenMin
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine findNewtonMin(id, is, eb, et, kmax_kkr, MinChannels,    &
                            MinValues, numMin, numMax)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TWO, HALF, CZERO, CONE
   use MathParamModule, only : TEN2m6, TEN2m8
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getSineMatrix, solveSingleScattering
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, kmax_kkr
   integer (kind=IntKind), intent(in) :: numMax
   integer (kind=IntKind), intent(inout) :: MinChannels(kmax_kkr)
   integer (kind=IntKind), intent(inout) :: numMin
   integer (kind=IntKind), parameter :: maxNewt = 10
   integer (kind=IntKind) :: n, ic, kl
!
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind), intent(inout) :: MinValues(kmax_kkr)
   real (kind=RealKind), parameter :: de = TEN2m6
   real (kind=RealKind) :: e0, eStep, s1, s2, s3
!
   complex (kind=CmplxKind) :: e, dets1, dets2, dets3
   complex (kind=CmplxKind) :: s1diag(kmax_kkr), s2diag(kmax_kkr)
   complex (kind=CmplxKind) :: s3diag(kmax_kkr)
   complex (kind=CmplxKind) :: s0(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind) :: sT(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind) :: sTs(kmax_kkr,kmax_kkr)
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
!
   interface
      subroutine calcDet(mat,kmax,det,diag)
      use KindParamModule, only : IntKind, CmplxKind
!
      integer (kind=IntKind), intent(in) :: kmax
!
      complex (kind=CmplxKind), intent(out) :: det
      complex (kind=CmplxKind), optional :: diag(kmax)
      complex (kind=CmplxKind), intent(in) :: mat(kmax,kmax)
      end subroutine calcDet
   end interface
!
   if (eb > et) then
      call ErrorHandler('findNewtonMin','eb > et',eb,et)
   endif
!
   do ic = 1, numMax
      kl = MinChannels(ic)
      if(kl == 0) cycle
      e0 = MinValues(ic)
!
      do n = 1, maxNewt
         e = cmplx(e0-de,ZERO,kind=CmplxKind)
!        ---------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ---------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ---------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ---------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ---------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,            &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ---------------------------------------------------------------
         call calcDet(sTs,kmax_kkr,dets1,s1diag)
!        ---------------------------------------------------------------
!      
         e = cmplx(e0,ZERO,kind=CmplxKind)
!        ---------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ---------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ---------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ---------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ---------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,            &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ---------------------------------------------------------------
         call calcDet(sTs,kmax_kkr,dets2,s2diag)
!        ---------------------------------------------------------------
!      
         e = cmplx(e0+de,ZERO,kind=CmplxKind)
!        ---------------------------------------------------------------
         call solveSingleScattering(is,id,e,CZERO)
!        ---------------------------------------------------------------
         sin_mat => getSineMatrix(is,id)
!        ---------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!        ---------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,s0,1,sT,1)
!        ---------------------------------------------------------------
         call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,            &
                    sT,kmax_kkr,s0,kmax_kkr,CZERO,sTs,kmax_kkr)
!        ---------------------------------------------------------------
         call calcDet(sTs,kmax_kkr,dets3,s3diag)
!        ---------------------------------------------------------------
!
         s1 = real(s1diag(kl), RealKind)
         s2 = real(s2diag(kl), RealKind)
         s3 = real(s3diag(kl), RealKind)
         write(6,'(4x,a,3(4x,a,e20.12))')'Newton:','s1=',s1,'s2=',    &
                                         s2,'s3=',s3
         eStep = HALF*de*(s3-s1)/(s3-TWO*s2+s1)
         write(6,'(4x,a,e20.12)')'Newton: eStep=',eStep
         if(abs(eStep) < TEN2m8) exit
         e0 = e0-eStep
         write(6,'(8x,a,i4,4x,a,f16.12)')'kl=',kl,'new e0:',e0
      enddo
!
      if(e0 < eb .or. e0 > et) then
         MinChannels(ic) = 0
         numMin = numMin-1
      endif
!
      MinValues(ic) = e0
   enddo
!
   end subroutine findNewtonMin
