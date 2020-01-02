program main
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use TimerModule, only : initTimer, getTime
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : initMPP, endMPP, initParent,startChildTasks
   use MPPModule, only : getMyPE, getNumPEs, GlobalSum
!
   use InputModule, only : initInput, endInput, readInputData
   use InputModule, only : openInputFile, closeInputFile
   use InputModule, only : getKeyValue, printKeyNames
!
   use OutputModule, only : initOutput, endOutput, getStandardOutputLevel
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
   use PolyhedraModule, only : printPolyhedron
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential, writePotential
   use PotentialModule, only : getPotEf
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType
!
   use ScfDataModule, only : initScfData, nscf, printScfData
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop
   use ScfDataModule, only : EvBottom
   use ScfDataModule, only : isReadEmesh, isReadKmesh
   use ScfDataModule, only : inputpath, Temperature
   use ScfDataModule, only : isNonRelativisticCore
   use ScfDataModule, only : isNonRelativisticValence
   use ScfDataModule, only : isScalarRelativisticValence
   use ScfDataModule, only : getKmeshFileName
   use ScfDataModule, only : NumKMeshs, kGenScheme, Kdiv, Symmetrize
   use ScfDataModule, only : getEmeshFileName
   use ScfDataModule, only : ContourType, eGridType
   use ScfDataModule, only : NumEs, ErBottom, ErTop, EiBottom, EiTop
   use ScfDataModule, only : isKKR, isKKRCPA, isLSMS
   use ScfDataModule, only : getPotentialTypeParam
!
   use SystemModule, only : initSystem, endSystem
   use Systemmodule, only : printSystem, getBravaisLattice
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber
!
   use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc, getGlobalIndex
   use Atom2ProcModule, only : printAtom2ProcTable, getLocalNumAtoms
!
   use AtomModule, only : initAtom, endAtom
   use AtomModule, only : getMaxLmax, printAtom
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax
   use AtomModule, only : getGridData
   use AtomModule, only : getLocalAtomName, getLocalAtomicNumber
   use AtomModule, only : getLocalAtomPosition, getLocalEvecOld
   use AtomModule, only : getInPotFileName, getInPotFileForm
   use AtomModule, only : getOutPotFileName, getOutPotFileForm
   use AtomModule, only : getInValDenFileName, getInValDenFileForm
   use AtomModule, only : getOutValDenFileName, getOutValDenFileForm
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use NeighborModule, only : initNeighbor, endNeighbor, getNumNeighbor
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
!
   use CoreStatesModule, only : initCoreStates, calCoreStates, endCoreStates
   use CoreStatesModule, only : readCoreStates, writeCoreStates
   use CoreStatesModule, only : printCoreStates
!
   use ContourModule, only : initContour, setupContour, printContour
   use ContourModule, only : getNumEs, getEPoint, getEWeight
!
   use BZoneModule, only : initBZone, printBZone, endBZone
!
!   use ValenceStatesModule, only : initValenceStates, printValenceStates
!   use ValenceStatesModule, only : calValenceStates
!   use ValenceStatesModule, only : writeValenceDensity, readValenceDensity
!
   use MultipleScatteringModule, only : initMultipleScattering,       &
                                        endMultipleScattering,        &
                                        getTauijLocalSpinBlock
!
   use WriteMatrixModule, only : writeMatrix
!
   implicit none
!
   character (len=80) :: info_table(1), info_path(1)
!
   integer (kind=IntKind) :: MyPE, NumPEs
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: i, j, id, ig, ia, isi, isj, nn
   integer (kind=IntKind) :: lmax_max, NumAtoms
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: rel
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: LocalAtomicNumber(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: fit_pot(:)
   integer (kind=IntKind), allocatable :: fit_rho(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: LocalAtomPosi(:,:)
   real (kind=RealKind), allocatable :: LocalEvec(:,:)
   real (kind=RealKind) :: rmt, rend
   real (kind=RealKind) :: Efermi
   real (kind=RealKind) :: TotalEnergy
   real (kind=RealKind) :: q_ws, q_mt, m_ws, m_mt, evec(3)
!
   complex (kind=CmplxKind), pointer :: Ewght(:)
   complex (kind=CmplxKind), pointer :: Epoint(:)
   complex (kind=CmplxKind), pointer :: Tauij(:,:)
   complex (kind=CmplxKind) :: energy
!
   character (len=10) :: hosts(2)
   character (len=7) :: exec(2)
   integer (kind=IntKind) :: num_tasks(2)
!
!   hosts(1) = 'dy.theory'
!   hosts(2) = 'lu'
   exec(1) = 'exe_dy'
   exec(2) = 'exe_nb'
   num_tasks = 1
!
!  -------------------------------------------------------------------
   call initTimer()
!  -------------------------------------------------------------------
   call initMPP()
!   call initParent()
!   call startChildTasks(1,hosts,exec,num_tasks)
!  -------------------------------------------------------------------
   MyPE = getMyPE()
   NumPEs = getNumPEs()

   write(6,*)" MPPmodule :: MyPE, NumPEs",MyPE,NumPEs
!  -------------------------------------------------------------------
!
   if (MyPE == 0) then
   write(6,'(/)')
   write(6,'(12x,a)')'********************************************************'
   write(6,'(12x,a)')'*                                                      *'
   write(6,'(12x,a)')'*    Full-Potential Multiple Scattering Theory Based   *'
   write(6,'(12x,a)')'*                                                      *'
   write(6,'(12x,a)')'*  Ab Initio Electronic Structure Calculation Package  *'
   write(6,'(12x,a)')'*                                                      *'
   write(6,'(12x,a)')'*                  Version 2.0 Beta 1                  *'
   write(6,'(12x,a)')'*                                                      *'
   write(6,'(12x,a)')'********************************************************'
   write(6,'(/,80(''=''))')
   endif
!
!  ===================================================================
!  call readInput to obtain input data................................
!      Proc 0: read input data file and broadcast to other processors
!      Other : wait and receive the data from Proc 0.
!  -------------------------------------------------------------------
   call initInput()
!  -------------------------------------------------------------------
   call readInputData(5,def_id)
!  -------------------------------------------------------------------
   call getKeyValue(def_id,'Info Table File Name',info_table,1)
   call getKeyValue(def_id,'Path to Info Table File',info_path,1)
!  -------------------------------------------------------------------
   call openInputFile(10,trim(adjustl(info_path(1)))//adjustl(info_table(1)))
   call readInputData(10,info_id)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call initOutput(def_id)
   print_level = getStandardOutputLevel()
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call initSystem(def_id)
   call initScfData(def_id)
!  -------------------------------------------------------------------
   NumAtoms = getNumAtoms()
!  -------------------------------------------------------------------
   call initAtom2Proc(NumAtoms)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  call initAtom and setAtomData to setup Atom Module.................
!  -------------------------------------------------------------------
   call initAtom(info_id,istop,print_level)
!  -------------------------------------------------------------------
   if (print_level >= 0) then
!     ----------------------------------------------------------------
      call printSystem()
      call printScfData()
      call printAtom()
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  setup factors
!  -------------------------------------------------------------------
   lmax_max=getMaxLmax()
!  -------------------------------------------------------------------
   call initGauntFactors(lmax_max, istop, print_level)
!  -------------------------------------------------------------------
!
   LocalNumAtoms=getLocalNumAtoms(MyPE)
!
   allocate(AtomPosition(3,NumAtoms), AtomicNumber(NumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
      AtomicNumber(i)=getAtomicNumber(i)
   enddo
!
!  ===================================================================
!  initialize LIZ
!  -------------------------------------------------------------------
   call initNeighbor(LocalNumAtoms)
!  -------------------------------------------------------------------
   call setupLizNeighbor()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize radial grid
!  -------------------------------------------------------------------
   call initRadialGrid(LocalNumAtoms, istop, print_level)
!  -------------------------------------------------------------------
!
   bravais(1:3,1:3)=getBravaisLattice()
!  -------------------------------------------------------------------
   call initPolyhedra(LocalNumAtoms,bravais,istop,print_level)
!  -------------------------------------------------------------------
!
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i,MyPE)
!     ----------------------------------------------------------------
      call getGridData(i,ndivin,ndivout,nmult)
!     ----------------------------------------------------------------
      call genPolyhedron(i,ig,NumAtoms,AtomPosition)
!     ----------------------------------------------------------------
      if (print_level >= 0) then
!        -------------------------------------------------------------
         call printPolyhedron(i)
!        -------------------------------------------------------------
      endif
!     ----------------------------------------------------------------
      call genRadialGrid(i,getInscrSphRadius(i),getOutscrSphRadius(i), &
                         ndivin,ndivout,nmult)
!     ----------------------------------------------------------------
      call printRadialGrid(i)
!     ----------------------------------------------------------------
   enddo
!
!  ===================================================================
!  initialize potential module
!  ===================================================================
   allocate(fit_pot(LocalNumAtoms), lmax_pot(LocalNumAtoms))
   allocate(fit_rho(LocalNumAtoms), lmax_rho(LocalNumAtoms))
   allocate(lmax_kkr(LocalNumAtoms), lmax_phi(LocalNumAtoms))
   do i=1,LocalNumAtoms
      fit_pot(i)=0                 ! = 0, no fit;  = 1, fit.
      fit_rho(i)=0                 ! = 0, no fit;  = 1, fit.
      lmax_kkr(i)=getKKRLmax(i)
      lmax_phi(i)=getPhiLmax(i)
      lmax_rho(i)=getRhoLmax(i)
      lmax_pot(i)=getPotLmax(i)
   enddo
!
!  -------------------------------------------------------------------
   call initPotentialType(getPotentialTypeParam())
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,fit_pot,lmax_pot,n_spin_pola,     &
                      istop,print_level)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  read potential data
!  ===================================================================
   do i=1,LocalNumAtoms
!     ----------------------------------------------------------------
      call readPotential(i,getLocalAtomName(i),                       &
                         getInPotFileName(i),getInPotFileForm(i))
!     ----------------------------------------------------------------
   enddo
!
!  ===================================================================
!  initialize core states module
!  -------------------------------------------------------------------
   call initCoreStates(LocalNumAtoms,EvBottom,n_spin_pola,            &
                       isNonRelativisticCore(),istop,print_level)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  read core states and density data of local atoms from file
!  ===================================================================
   do i=1,LocalNumAtoms
!     ----------------------------------------------------------------
      call readCoreStates(i,getLocalAtomName(i),                      &
                          getInPotFileName(i),getInPotFileForm(i))
!     ----------------------------------------------------------------
   enddo
!
   if (isReadEmesh()) then
!     ----------------------------------------------------------------
      call initContour(getEmeshFileName(),istop,print_level)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call initContour(ContourType,eGridType,NumEs,Temperature,istop,print_level)
!     ----------------------------------------------------------------
   endif
!
   if (isNonRelativisticValence()) then
      rel = 0
   else if (isScalarRelativisticValence()) then
      rel = 1
   else 
      rel = 2
   endif
!
   allocate(LocalAtomPosi(3,LocalNumAtoms), LocalAtomicNumber(LocalNumAtoms))
   allocate(GlobalIndex(LocalNumAtoms), LocalEvec(3,LocalNumAtoms))
   do id=1,LocalNumAtoms
      LocalAtomPosi(1:3,id)=getLocalAtomPosition(id)
      LocalEvec(1:3,id)=getLocalEvecOld(id)
      GlobalIndex(id)=getGlobalIndex(id,MyPE)
      LocalAtomicNumber(id)=getLocalAtomicNumber(id)
   enddo
!
!  -------------------------------------------------------------------
   call initMultipleScattering( LocalNumAtoms,LocalAtomicNumber,GlobalIndex, &
                                lmax_kkr,lmax_phi,lmax_pot,LocalAtomPosi,    &
                                LocalEvec,n_spin_pola,n_spin_cant,rel, &
                                istop,print_level)
!  -------------------------------------------------------------------
!
   energy = (0.753825D+00,0.100000D-02)
   if ( MyPE==0 ) then
      write(6,*)"     energy :: ",energy
   endif
   do i = 1,LocalNumAtoms
      do j= 1, LocalNumAtoms
         print *,"before call getTauijLocalSpinBlock" 
         Tauij => getTauijLocalSpinBlock( i, energy, 1, 1, i, j )
         if (MyPE == 0) then
            write(6,*)"   LocalAtom ::", i
            call writeMatrix( "Tauij",Tauij,n_spin_cant*(lmax_kkr(i)+1)**2,    &
                              n_spin_cant*(lmax_kkr(i)+1)**2 )
         endif
            nullify(Tauij)
      enddo
   enddo
!
   call endCoreStates()
   call endPotential()
   call endPolyhedra()
   call endRadialGrid()
   call endGauntFactors()
   call endNeighbor()
   call endPotentialType()
   call endAtom()
   call endAtom2Proc()
   call endSystem()
   call endOutput()
   call endMPP()
!  -------------------------------------------------------------------
   stop 'End of the program. OK!'
end program main
