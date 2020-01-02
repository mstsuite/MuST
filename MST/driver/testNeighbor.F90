program testNeighbor
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8
!
   use ChemElementModule, only : getZval
!
   use TimerModule, only : initTimer, getTime
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : initMPP, endMPP, syncAllPEs, GlobalMax, GlobalMin
   use MPPModule, only : getMyPE
!
   use GroupCommModule, only : initGroupComm
!
   use ProcMappingModule, only : initProcMapping, endProcMapping,            &
                                 createParallelization
!
   use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc
   use Atom2ProcModule, only : getGlobalIndex, getLocalNumAtoms
!
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, isDataStorageExisting
!
   use InputModule, only : initInput, endInput, readInputData
   use InputModule, only : openInputFile, closeInputFile
   use InputModule, only : getKeyValue, printKeyNames, getTableIndex
!
   use OutputModule, only : initOutput, endOutput, getStandardOutputLevel
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
!
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : printStepFunction, testStepFunction
!
   use ScfDataModule, only : ngaussr, ngaussq
   use ScfDataModule, only : initScfData, nscf, printScfData, potwrite, movie
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop
   use ScfDataModule, only : EvBottom, ErTop
   use ScfDataModule, only : isReadKmesh
   use ScfDataModule, only : inputpath, Temperature
   use ScfDataModule, only : isNonRelativisticCore
   use ScfDataModule, only : isNonRelativisticValence
   use ScfDataModule, only : isScalarRelativisticValence
   use ScfDataModule, only : getKmeshFileName
   use ScfDataModule, only : NumKMeshs, kGenScheme, Kdiv, Symmetrize
   use ScfDataModule, only : isKKR, isKKRCPA, isLSMS
   use ScfDataModule, only : getPotentialTypeParam
   use ScfDataModule, only : isChargeMixing
   use ScfDataModule, only : eftol, etol, ptol, rmstol
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType,      &
                                   isASAPotential, isMuffinTinPotential,     &
                                   isMuffinTinASAPotential, isTestPotential, &
                                   isMuffinTinTestPotential, isFullPotential,&
                                   printPotentialType
!
   use SystemModule, only : initSystem, endSystem
   use Systemmodule, only : printSystem, getBravaisLattice
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber,   &
                            getNumVacancies
   use SystemModule, only : getUniformGridParam
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
   use AtomModule, only : getLocalAtomName, getLocalAtomicNumber
   use AtomModule, only : getLocalAtomPosition, getLocalEvecOld
   use AtomModule, only : getMixingParam4Evec
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use NeighborModule, only : initNeighbor, endNeighbor, getNumNeighbors, &
                              printNeighbor
!
   implicit none
!
   character (len=80) :: info_table, info_path
   character (len=160) :: itname
   character (len=2) :: anm
!
   integer (kind=IntKind) :: MaxLIZ, MinLIZ
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: i, id, ig, ia, is, jl, ir, nr, iscf
   integer (kind=IntKind) :: n_movie, n_potwrite, funit
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind) :: NumRotations
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind) :: ng_uniform(3)
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: LocalAtomicNumber(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: fit_pot(:)
   integer (kind=IntKind), allocatable :: fit_rho(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
!
   real (kind=RealKind), allocatable :: LocalNumValenceElectrons(:)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: LocalAtomPosi(:,:)
   real (kind=RealKind), allocatable :: LocalEvec(:,:)
!
   real (kind=RealKind), allocatable :: rho_rms(:,:),pot_rms(:,:)
   real (kind=RealKind) :: ef_diff, tote_diff, max_rms(4)
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind) :: rmt, rend
   real (kind=RealKind) :: Efermi
   real (kind=RealKind) :: TotalEnergy
   real (kind=RealKind) :: q_ws, q_mt, m_ws, m_mt, evec(3)
   real (kind=RealKind) :: v0, Beta(3), Alpha(3), val
   real (kind=RealKind) :: t0
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!
   complex (kind=CmplxKind), pointer :: pot_l(:)
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
!
!  ===================================================================
!  Initialize the processes mapping module that determines how the
!  parallization will be performed
!  -------------------------------------------------------------------
   call initProcMapping(NumAtoms, 1, 1, isFullPotential(), istop, 0, NumAtoms)
!  -------------------------------------------------------------------
   call createParallelization()
!  -------------------------------------------------------------------
   call initAtom2Proc(NumAtoms, NumAtoms)
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
   write(6,'(/)')
   write(6,'(12x,a)')'********************************************************'
   write(6,'(12x,a)')'*                                                      *'
   write(6,'(12x,a)')'*    Full-Potential Multiple Scattering Theory Based   *'
   write(6,'(12x,a)')'*                                                      *'
   write(6,'(12x,a)')'*  Ab Initio Electronic Structure Calculation Package  *'
   write(6,'(12x,a)')'*                                                      *'
   write(6,'(12x,a)')'*                    Version 2.0d                      *'
   write(6,'(12x,a)')'*                                                      *'
   write(6,'(12x,a)')'********************************************************'
   write(6,'(/,80(''=''))')
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
      call printScfData()
      call printSystem()
      call printAtom()
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   bravais(1:3,1:3)=getBravaisLattice()
!  -------------------------------------------------------------------
!
   LocalNumAtoms=getLocalNumAtoms()
!
   allocate(atom_print_level(1:LocalNumAtoms))
   do i=1,LocalNumAtoms
      atom_print_level(i) = getStandardOutputLevel(i)
   enddo
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
   call initNeighbor(LocalNumAtoms,atom_print_level)
!  -------------------------------------------------------------------
   call setupLizNeighbor(atom_print_level)
!  -------------------------------------------------------------------
!
   MinLIZ = 1000000
   MaxLIZ = 0
   do i = 1, LocalNumAtoms
      MinLIZ = min(MinLIZ,getNumNeighbors(i))
      MaxLIZ = max(MaxLIZ,getNumNeighbors(i))
   enddo
!
   call GlobalMax(MaxLIZ)
   call GlobalMin(MinLIZ)
!
   if (getMyPE() == 0) then
      print *,'Max LIZ = ',MaxLIZ
      print *,'Min LIZ = ',MinLIZ
   endif
!
   if (MaxLIZ /= MinLIZ) then
      do i=1, LocalNumAtoms
         if (MinLIZ == getNumNeighbors(i)) then
            print *,'Global Index = ',getGlobalIndex(i)
            call printNeighbor(i)
         else if (MaxLIZ == getNumNeighbors(i)) then
            print *,'Global Index = ',getGlobalIndex(i)
            call printNeighbor(i)
         endif
      enddo
   endif
   call syncAllPEs()
!
   call endNeighbor()
   call endPotentialType()
   call endAtom()
!  call endAtom2Proc()
   call endOutput()
   call endSystem()
   call endDataServiceCenter()
   call endMPP()
   stop
end program testNeighbor
