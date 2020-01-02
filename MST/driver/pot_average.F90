program pot_average
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8
!
   use TimerModule, only : initTimer, getTime
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : initMPP, endMPP, syncAllPEs, GlobalMax
   use MPPModule, only : getMyPE, getNumPEs, GlobalSum
!
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter
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
!  use PolyhedraModule, only : initPolyhedra, endPolyhedra
!  use PolyhedraModule, only : genPolyhedron
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
   use PotentialModule, only : readPotential, getPotential
!
   use ScfDataModule, only : ngaussr, ngaussq
   use ScfDataModule, only : initScfData
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop
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
   use Systemmodule, only : printSystem, getBravaisLattice
   use SystemModule, only : getNumAtoms, getAtomPosition
   use SystemModule, only : getNumAtomTypes
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,   &
                                    calSymmetryFlags ! ,getSymmetryFlags
!
!  ===================================================================
!  initAtom2Proc and endAtom2Proc are called in SystemModule
!  ===================================================================
!  use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc
   use Atom2ProcModule, only : getGlobalIndex
   use Atom2ProcModule, only : printAtom2ProcTable, getLocalNumAtoms
!
   use AtomModule, only : initAtom, endAtom
   use AtomModule, only : getMaxLmax, printAtom
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax
   use AtomModule, only : getGridData
   use AtomModule, only : getLocalAtomPosition
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid, getGridRadius
!
!  ===================================================================
   use SurfElementsModule, only : initSurfElements, endSurfElements
   use SurfElementsModule, only : genGaussPoints, genGaussSphHarmonics
   use SurfElementsModule, only : getSurfAverage
!  ===================================================================
!
   implicit none
!
   logical :: FileNamed = .false.
!
   character (len=12) :: str_stdin= 'i_lsms_stdin'
   character (len=80) :: info_table, info_path
   character (len=160) :: itname
   character (len=2)  :: anm
   character (len=8)  :: exec_date
   character (len=10) :: exec_time
   character (len=200) :: FileName
!
   integer (kind=IntKind) :: MyPE, NumPEs
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: i, id, ig, ia, ir, nr, is
   integer (kind=IntKind) :: lmax_max, lmax_kkr_max, NumAtoms, jmax_pot
   integer (kind=IntKind) :: lmax_rho_max, lmax_pot_max
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: node_print_level, lprint
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind) :: MaxVal_Integer(2)
   integer (kind=IntKind) :: fstatus
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind) :: rmt, rend
   real (kind=RealKind) :: t0, t1, t2, t3
   real (kind=RealKind) :: pot_aver(2)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
!
!  -------------------------------------------------------------------
   call initTimer()
!  -------------------------------------------------------------------
   t0 = getTime()
!
!  -------------------------------------------------------------------
   call initMPP()
!  -------------------------------------------------------------------
   MyPE = getMyPE()
   NumPEs = getNumPEs()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  call initDataServiceCenter to startup the Data Storage Service
!  -------------------------------------------------------------------
   call initDataServiceCenter()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  call readInput to obtain input data................................
!      Proc 0: read input data file and broadcast to other processors
!      Other : wait and receive the data from Proc 0.
!  -------------------------------------------------------------------
   call initInput()
!  -------------------------------------------------------------------
   inquire(unit=5,name=FileName,named=FileNamed)
   if (FileNamed) then
!     ----------------------------------------------------------------
      call readInputData(5,def_id)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call openInputFile(7,str_stdin)
      call readInputData(7,def_id)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   fstatus = getKeyValue(def_id,'Info Table File Name',info_table)
   fstatus = getKeyValue(def_id,'Current File Path',info_path)
!  -------------------------------------------------------------------
   itname = trim(info_path)//info_table
   info_id = getTableIndex(itname)
   if (info_id < 1) then
!     ----------------------------------------------------------------
      call openInputFile(10,itname)
      call readInputData(10,info_id)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  initialize ScfData module
!  -------------------------------------------------------------------
   call initScfData(def_id)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize Potential Type
!  -------------------------------------------------------------------
   call initPotentialType(getPotentialTypeParam())
!  -------------------------------------------------------------------
!  ===================================================================
!
!  ===================================================================
   call initSystem(def_id)
!  -------------------------------------------------------------------
   NumAtoms = getNumAtoms()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  set up print level
!  -------------------------------------------------------------------
   call initOutput(def_id)
   node_print_level = getStandardOutputLevel()
!  -------------------------------------------------------------------
!
   if (node_print_level >= 0) then
   call date_and_time(exec_date,exec_time)
   write(6,'(13a,i5,a)')'Execution starts at ',                               &
                   exec_time(1:2),':',exec_time(3:4),':',exec_time(5:6),', ', &
                   exec_date(5:6),'-',exec_date(7:8),'-',exec_date(1:4),      &
                   ' on ',NumPEs,' nodes'
   endif
!
!  ===================================================================
!  call initAtom and setAtomData to setup Atom Module.................
!  -------------------------------------------------------------------
   call initAtom(info_id,istop,node_print_level)
!  -------------------------------------------------------------------
   if (node_print_level >= 0) then
!     ----------------------------------------------------------------
      call printPotentialType()
      call printSystem()
      call printAtom()
!     ----------------------------------------------------------------
   endif
!
   bravais(1:3,1:3)=getBravaisLattice()
!
!  ===================================================================
!  setup factors
!  -------------------------------------------------------------------
   lmax_max=getMaxLmax()
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_max)
!  -------------------------------------------------------------------
   call initGauntFactors(lmax_max, istop, node_print_level)
!  -------------------------------------------------------------------
!
   LocalNumAtoms=getLocalNumAtoms(MyPE)
!
   allocate(atom_print_level(1:LocalNumAtoms))
   do i=1,LocalNumAtoms
      atom_print_level(i) = getStandardOutputLevel(i)
   enddo
!
   allocate(AtomPosition(3,NumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
   enddo
!
!  ===================================================================
!  initialize atomic polyhedra
!  -------------------------------------------------------------------
!  call initPolyhedra(LocalNumAtoms,bravais,istop,node_print_level)
!  -------------------------------------------------------------------
!
!  *******************************************************************
!
!  ===================================================================
!  initialize radial grid
!  -------------------------------------------------------------------
   call initRadialGrid(LocalNumAtoms, istop, node_print_level)
!  -------------------------------------------------------------------
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i,MyPE)
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
      if (isMuffinTinPotential() .or. isMuffinTinTestPotential()) then
!        -------------------------------------------------------------
!         call genRadialGrid(i,xstart,getInscrSphRadius(i),           &
!                            getWignerSeitzRadius(i),ndivin)
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
!
!  *******************************************************************
!
!  ===================================================================
!  setup the Lmax values.             
!  ===================================================================
   allocate(lmax_pot(LocalNumAtoms))
   allocate(lmax_rho(LocalNumAtoms))
   allocate(lmax_kkr(LocalNumAtoms), lmax_phi(LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms))
   lmax_kkr_max = 0
   lmax_rho_max = 0
   lmax_pot_max = 0
   do i=1,LocalNumAtoms
      lmax_kkr(i)=getKKRLmax(i)
      lmax_phi(i)=getPhiLmax(i)
      lmax_rho(i)=getRhoLmax(i)
      lmax_pot(i)=getPotLmax(i)
      lmax_kkr_max = max(lmax_kkr_max,lmax_kkr(i))
      lmax_rho_max = max(lmax_rho_max,lmax_rho(i))
      lmax_pot_max = max(lmax_pot_max,lmax_pot(i))
      lmax_step(i) = max(lmax_rho(i), lmax_phi(i),                    &
                         lmax_pot(i)+lmax_phi(i)+lmax_kkr(i))
   enddo
!  ===================================================================
   MaxVal_Integer(1) = lmax_rho_max
   MaxVal_Integer(2) = lmax_pot_max
!  -------------------------------------------------------------------
   call GlobalMax(MaxVal_Integer,2)
!  -------------------------------------------------------------------
   lmax_rho_max = MaxVal_Integer(1)
   lmax_pot_max = MaxVal_Integer(2)
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  initialize step function module
!  ===================================================================
   allocate( ngr(LocalNumAtoms), ngt(LocalNumAtoms) )
   do i=1,LocalNumAtoms
      ngr(i) = ngaussr
      ngt(i) = ngaussq
   enddo
!  -------------------------------------------------------------------
   call initStepFunction(LocalNumAtoms,lmax_max,lmax_step,ngr,ngt,istop,node_print_level)
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
   call initSystemSymmetry( LocalNumAtoms, LocalNumAtoms, lmax_pot, lmax_step, &
                            atom_print_level )
!  -------------------------------------------------------------------
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  initialize potential module
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot,lmax_step,               &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  read potential data
!  -------------------------------------------------------------------
   call readPotential()
!  -------------------------------------------------------------------
   call syncAllPEs()
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  calculate the average of potential
!  -------------------------------------------------------------------
   call initSurfElements(istop,-1)
   call genGaussPoints()
   call genGaussSphHarmonics(lmax_pot_max)
!  -------------------------------------------------------------------
   do is = 1, n_spin_pola
      pot_aver(is) = ZERO
      do i = 1, LocalNumAtoms
         nr = getNumRmesh(i)
         r_mesh => getRmesh(i)
         pot_l => getPotential(i,1,is)
         pot_aver(is) = pot_aver(is)                                  &
                       +getSurfAverage(i,nr,lmax_pot(i),r_mesh,pot_l)
      enddo
   enddo
!  -------------------------------------------------------------------
   call GlobalSum(pot_aver,2)
!  -------------------------------------------------------------------
   pot_aver(1:2) = pot_aver(1:2)/real(NumAtoms,kind=RealKind)
   write(6,'(/,a,2f12.5,/)')'Surface average potential:',pot_aver(1:n_spin_pola)
!
!  ===================================================================
!  clean up the allocated spaces.
!  ===================================================================
   nullify(r_mesh, pot_l)
   deallocate( atom_print_level )
   deallocate(AtomPosition,lmax_kkr,lmax_phi,lmax_pot,lmax_rho)
   deallocate( lmax_step )
!
   call endSurfElements()
   call endPotential()
   call endStepFunction()
   call endRadialGrid()
   call endSystemSymmetry()
   call endGauntFactors()
   call endSphericalHarmonics()
   call endPotentialType()
   call endAtom()
   call endSystem()
   call endDataServiceCenter()
   call endMPP()
!
!  ===================================================================
   if (node_print_level >= 0) then
   call date_and_time(exec_date,exec_time)
   write(6,'(/,12a)')'Execution ends at ',                            &
                   exec_time(1:2),':',exec_time(3:4),':',exec_time(5:6),', ', &
                   exec_date(5:6),'-',exec_date(7:8),'-',exec_date(1:4)
   write(6,'(80(''-''))')
   endif
!
   call endOutput()
!
   if (node_print_level >= 0) then
      stop 'End of the program. OK!'
   else
      stop
   endif
end program pot_average
