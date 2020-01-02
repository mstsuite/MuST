program testParallelIO
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8
!
   use TimerModule, only : initTimer, getTime
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : initMPP, endMPP, GlobalSum, MyPE, NumPEs
!
   use GroupCommModule, only : initGroupComm, endGroupComm, getGroupID
!
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter,  &
                                       isDataStorageExisting, &
                                       createDataStorage,     &
                                       setDataStorage2Value,  &
                                       setDataStorageLDA,     &
                                       getDataStorage, RealMark, RealType
!
   use InputModule, only : initInput, endInput, readInputData
   use InputModule, only : openInputFile, closeInputFile
   use InputModule, only : getKeyValue, printKeyNames, getTableIndex
!
   use OutputModule, only : initOutput, endOutput,  &
                            getStandardOutputLevel
!
!  ===================================================================
!  initPolyhedra and endPolyhedra are called in SystemModule
!  ===================================================================
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
   use PotentialModule, only : readPotential, writePotential
   use PotentialModule, only : getPotEf
!
   use ScfDataModule, only : ngaussr, ngaussq
   use ScfDataModule, only : initScfData, printScfData
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop, EvBottom, isNonRelativisticCore
   use ScfDataModule, only : inputpath
   use ScfDataModule, only : getPotentialTypeParam
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
   use SystemVolumeModule, only : initSystemVolume, endSystemVolume,         &
                                  printSystemVolume, updateSystemVolume,     &
                                  setSystemVolumeMT, getAtomicVPVolume
!
   use ProcMappingModule, only : initProcMapping, endProcMapping,            &
                                 createParallelization
!
   use ParallelIOModule, only : initParallelIO, endParallelIO,               &
                                isInputProc, isOutputProc
!
!  ===================================================================
!  initAtom2Proc and endAtom2Proc are called in SystemModule
!  ===================================================================
   use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc
   use Atom2ProcModule, only : getGlobalIndex
   use Atom2ProcModule, only : printAtom2ProcTable, getLocalNumAtoms
!
   use AtomModule, only : initAtom, endAtom
   use AtomModule, only : printAtom
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax
   use AtomModule, only : getGridData, getLocalNumSpecies
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid, getNumRmesh
!
   use MadelungModule, only : initMadelung, endMadelung
!
   use CoreStatesModule, only : initCoreStates, calCoreStates, endCoreStates
   use CoreStatesModule, only : readCoreStates
   use CoreStatesModule, only : printCoreStates
!
   implicit none
!
   character (len=8)  :: exec_date
   character (len=10) :: exec_time
!
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: i, j, k, id, ig, nstep, is, nr, ns, ip
   integer (kind=IntKind) :: lmax_step_max, lmax_kkr_max, lmax_phi_max, lmax_rho_max, lmax_pot_max
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind) :: LocalNumAtoms, NumAtoms, NumIOProcs(2)
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind) :: MaxSpecies
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: DataSize(:)
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind), pointer :: rho0_old(:,:)
   real (kind=RealKind), pointer :: mom0_old(:,:)
   real (kind=RealKind), pointer :: rho0_new(:,:)
   real (kind=RealKind), pointer :: mom0_new(:,:)
!
   real (kind=RealKind) :: Efermi
   real (kind=RealKind) :: t0, t1, t2, t3
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!  real (kind=RealKind), parameter :: xstart = -.913096740000D+01
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
   call initProcMapping(NumAtoms, 1, 1, isFullPotential(), istop, -1, NumAtoms)
!  -------------------------------------------------------------------
   call createParallelization()
!  -------------------------------------------------------------------
   call initAtom2Proc(NumAtoms, NumAtoms)
!  -------------------------------------------------------------------
   call initSystemVolume()
!  -------------------------------------------------------------------
!
   LocalNumAtoms=getLocalNumAtoms()
!
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
      call printSystemVolume()
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
   allocate(lmax_step(LocalNumAtoms))
   lmax_kkr_max = 0
   lmax_phi_max = 0
   lmax_rho_max = 0
   lmax_pot_max = 0
   lmax_step_max = 0
   do i=1,LocalNumAtoms
      lmax_kkr(i) = getKKRLmax(i)
      lmax_phi(i) = getPhiLmax(i)
      lmax_rho(i) = getRhoLmax(i)
      lmax_pot(i) = getPotLmax(i)
      lmax_step(i) = max(2*lmax_phi(i), lmax_pot(i)+2*lmax_phi(i))
      lmax_kkr_max = max(lmax_kkr_max,lmax_kkr(i))
      lmax_phi_max = max(lmax_phi_max,lmax_phi(i))
      lmax_rho_max = max(lmax_rho_max,lmax_rho(i))
      lmax_pot_max = max(lmax_pot_max,lmax_pot(i))
      lmax_step_max = max(lmax_step_max,lmax_step(i))
   enddo
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(max(lmax_kkr_max+lmax_phi_max,2*lmax_step_max))
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(lmax_step_max,istop,0)
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
      call ErrorHandler('testParallelIO','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i)
      GlobalIndex(i)=ig
!     ----------------------------------------------------------------
      call getGridData(i,ndivin,ndivout,nmult)
!     ----------------------------------------------------------------
      call genPolyhedron(i,ig,NumAtoms,AtomPosition)
!     ----------------------------------------------------------------
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printPolyhedron(i)
!        call printPolyhedronBoundary(i)
!        -------------------------------------------------------------
      endif
      if (isMuffinTinPotential() .or. isMuffinTinTestPotential()) then
!        -------------------------------------------------------------
         call genRadialGrid(i,xstart,getInscrSphRadius(i),            &
                            getInscrSphRadius(i),                     &
                            getOutscrSphRadius(i),ndivin)
!        -------------------------------------------------------------
      else if (isASAPotential() .or. isMuffinTinASAPotential()) then
!        -------------------------------------------------------------
         call genRadialGrid(i,xstart,getWignerSeitzRadius(i),         &
                            getWignerSeitzRadius(i),                  &
                            getOutscrSphRadius(i),ndivin)
!                           getWignerSeitzRadius(i),ndivin)
!        -------------------------------------------------------------
      else
         if (getNeighborDistance(i,1)-getOutscrSphRadius(i) < TEN2m8) then
!           ----------------------------------------------------------
            call WarningHandler('main',                               &
                     'Ill condition found: Neighbor distance <= Rcs', &
                     getNeighborDistance(i,1),getOutscrSphRadius(i))
!           ----------------------------------------------------------
         endif
!        -------------------------------------------------------------
         call genRadialGrid(i,xstart,getInscrSphRadius(i),            &
                            getInscrSphRadius(i),                     &
                            getOutscrSphRadius(i),ndivin)
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
   call initStepFunction(LocalNumAtoms, lmax_step_max, lmax_step, ngr, ngt, &
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
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  initialize potential module
!  -------------------------------------------------------------------
   call initMadelung(LocalNumAtoms,NumAtoms,GlobalIndex,              &
                     lmax_rho_max,lmax_pot_max,bravais,AtomPosition,node_print_level)
!  -------------------------------------------------------------------
   call initSystemSymmetry( NumAtoms, LocalNumAtoms, lmax_pot, lmax_step, atom_print_level )
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot,lmax_step,      &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize core states module
!  -------------------------------------------------------------------
   call initCoreStates(LocalNumAtoms,EvBottom,n_spin_pola,         &
                       isNonRelativisticCore(),istop,node_print_level)
!  -------------------------------------------------------------------
!
   allocate( DataSize(LocalNumAtoms) )
   MaxSpecies = 0
   do id = 1,LocalNumAtoms
      MaxSpecies = max(MaxSpecies, getLocalNumSpecies(id))
   enddo    
!           
   do id = 1,LocalNumAtoms
      DataSize(id)  = (getNumRmesh(id)+1)*MaxSpecies
   enddo
!
   ip = 1
   do while (ip <= NumPEs) 
!     ----------------------------------------------------------------
      call initParallelIO(getGroupID('Unit Cell'),1,ip,ip)
!     ----------------------------------------------------------------
      if (isInputProc()) then
         NumIOProcs(1) = 1
      else
         NumIOProcs(1) = 0
      endif
      if (isOutputProc()) then
         NumIOProcs(2) = 1
      else
         NumIOProcs(2) = 0
      endif
!     ----------------------------------------------------------------
      call GlobalSum(NumIOProcs,2)
!     ----------------------------------------------------------------
!
!     ================================================================
!     read potential data
!     ----------------------------------------------------------------
      t1 = getTime()
!     ----------------------------------------------------------------
      call readPotential()
!     ----------------------------------------------------------------
      if (MyPE == 0) then
         write(6,'(/,a)')       '============================================='
         write(6,'(a,i5)')      'Number of processors reading potential data: ',NumIOProcs(1)
         write(6,'(a,f12.5,a)') 'Parallel input  time: ',getTime() - t1,' sec.'
         write(6,'(a,/)')       '============================================='
      endif
!
!     ================================================================
!     read core states and density data of local atoms from file
!     ================================================================
      call readCoreStates()
!     ----------------------------------------------------------------
      call calCoreStates()
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         if (atom_print_level(id) >= 0) then
!           ----------------------------------------------------------
            call printCoreStates(id)
!           ----------------------------------------------------------
         endif
      enddo
!
      if (.not.isDataStorageExisting('NewSphericalElectronDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'NewSphericalElectronDensity',  &
                                DataSize(1:LocalNumAtoms),RealType)  ! en extra space is used
!        -------------------------------------------------------------
         call setDataStorage2Value('NewSphericalElectronDensity',ZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)+1
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'NewSphericalElectronDensity',nr)
!           ----------------------------------------------------------
         enddo
      endif
!
      if (n_spin_pola == 2) then
         if (.not.isDataStorageExisting('NewSphericalMomentDensity')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'NewSphericalMomentDensity', &
                                   DataSize(1:LocalNumAtoms),RealType) ! en extra space is used
!           ----------------------------------------------------------
            call setDataStorage2Value('NewSphericalMomentDensity',ZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)+1
!              -------------------------------------------------------
               call setDataStorageLDA(id,'NewSphericalMomentDensity',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
!
      do id = 1,LocalNumAtoms
         nr = getNumRmesh(id)+1
         ns = getLocalNumSpecies(id)
         rho0_old => getDataStorage( id, 'OldSphericalElectronDensity',  &
                                     nr, ns, RealMark )
         rho0_new => getDataStorage( id, 'NewSphericalElectronDensity',  &
                                     nr, ns, RealMark )
         call dcopy(nr*ns,rho0_old,1,rho0_new,1)
         if ( n_spin_pola==2 ) then
            mom0_old => getDataStorage(id,'OldSphericalMomentDensity',   &
                                       nr, ns, RealMark )
            mom0_new => getDataStorage(id,'NewSphericalMomentDensity',   &
                                       nr, ns, RealMark )
            call dcopy(nr*ns,mom0_old,1,mom0_new,1)
         endif
      enddo
!
      t1 = getTime()
!     ----------------------------------------------------------------
      call writePotential()
!     ----------------------------------------------------------------
      if (MyPE == 0) then
         write(6,'(/,a)')       '============================================='
         write(6,'(a,i5)')      'Number of processors writing potential data: ',NumIOProcs(2)
         write(6,'(a,f12.5,a)') 'Parallel output time: ',getTime() - t1,' sec.'
         write(6,'(a,/)')       '============================================='
         call FlushFile(6)
      endif
!
!     ----------------------------------------------------------------
      call endParallelIO()
!     ----------------------------------------------------------------
      ip = ip*2
   enddo
!
   call endCoreStates()
   call endPotential()
   call endStepFunction()
   call endRadialGrid()
   call endSystemSymmetry()
   call endMadelung()
   call endGauntFactors()
   call endSphericalHarmonics()
   call endPotentialType()
   call endAtom()
   call endSystemVolume()
   call endAtom2Proc()
   call endProcMapping()
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
end program testParallelIO
