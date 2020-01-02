program mst2
!  use ISO_FORTRAN_ENV, only : compiler_version, compiler_options
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8, TEN2m6, TEN2m5, TEN2m3, PI4, THIRD
!
   use ChemElementModule, only : getZval
!
   use TimerModule, only : initTimer, getTime
!
   use ErrorHandlerModule, only : setErrorOutput, ErrorHandler, WarningHandler
!
   use MPPModule, only : initMPP, endMPP, syncAllPEs, GlobalMax
   use MPPModule, only : getMyPE, getNumPEs, getMaxWaits, setMaxWaits
!
   use GroupCommModule, only : initGroupComm, endGroupComm, getGroupID, getMyPEinGroup
!
   use ParallelIOModule, only : initParallelIO, endParallelIO
!
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, &
                                       getDataStorage, RealMark
!
   use PublicTypeDefinitionsModule, only : MixListRealStruct,  &
                                           MixListCmplxStruct, &
                                           UniformGridStruct
!
   use PublicParamDefinitionsModule, only : PrintDOSswitchOff, ButterFly
   use PublicParamDefinitionsModule, only : StandardInputFile
!
   use InputModule, only : initInput, endInput, readInputData
   use InputModule, only : openInputFile, closeInputFile
   use InputModule, only : getKeyValue, printKeyNames, getTableIndex
!
   use OutputModule, only : initOutput, endOutput,  &
                            isStandardOutToScreen,  &
                            getDensityPrintFlag,    &
                            getPotentialPrintFlag,  &
                            getStandardOutputLevel, isOutputAtomBased
!
!  ===================================================================
!  initPolyhedra and endPolyhedra are called in SystemModule
!  ===================================================================
!  use PolyhedraModule, only : initPolyhedra, endPolyhedra
!  use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
   use PolyhedraModule, only : getWignerSeitzRadius
   use PolyhedraModule, only : getNeighborDistance
   use PolyhedraModule, only : getPointLocationFlag
   use PolyhedraModule, only : printPolyhedron
   use PolyhedraModule, only : printPolyhedronBoundary
!
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : printStepFunction, testStepFunction
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential, writePotential
   use PotentialModule, only : getPotEf, setV0, printPot_L
   use PotentialModule, only : setPotential, setPotEf
   use PotentialModule, only : isSphericalInputFile, setPotentialOutsideMT
!
   use TestPotentialModule, only : initTestPotential,  &
                                   endTestPotential,   &
                                   readTestPotential,  &
                                   getTestPotential,   &
                                   getTestPotEf,       &
                                   getTestV0,          &
                                   getTestValenceNum
!
   use ScfDataModule, only : ngaussr, ngaussq
   use ScfDataModule, only : initScfData, nscf, printScfData, potwrite, movie
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : ntstep, tstep
   use ScfDataModule, only : istop, iterateEf
   use ScfDataModule, only : SingleSiteEContour, RealAxisSS
   use ScfDataModule, only : inputpath, Temperature
   use ScfDataModule, only : isNonRelativisticCore
   use ScfDataModule, only : isRelativisticValence
   use ScfDataModule, only : EvBottom, ErTop, ErBottom, EiTop, EiBottom
   use ScfDataModule, only : NumEs, ContourType, eGridType
   use ScfDataModule, only : isReadEmesh, getEmeshFileName
   use ScfDataModule, only : isReadKmesh, getKmeshFileName
   use ScfDataModule, only : NumKMeshs, kGenScheme, Kdiv, Symmetrize
   use ScfDataModule, only : isScreenKKR, isKKRCPA, isLSMS, isScreenKKR_LSMS, &
                             isKKR, isEmbeddedCluster, isSingleSite, setScfMethod
   use ScfDataModule, only : isExchangeParamNeeded
   use ScfDataModule, only : getPotentialTypeParam
   use ScfDataModule, only : isChargeMixing, isPotentialMixing
   use ScfDataModule, only : excorr_name, eftol, etol, ptol, rmstol
   use ScfDataModule, only : isLdaCorrectionNeeded, getUJfile
   use ScfDataModule, only : getSingleSiteSolverType, getDOSrunID
   use ScfDataModule, only : NumSS_IntEs, isSSIrregularSolOn
   use ScfDataModule, only : isFrozenCore
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType,      &
                                   isASAPotential, isMuffinTinPotential,     &
                                   isMuffinTinASAPotential, isTestPotential, &
                                   isMuffinTinTestPotential, isFullPotential,&
                                   printPotentialType, isSphericalPotential, &
                                   isMuffinTinFullPotential, setIsSphericalPotential
!
   use SystemModule, only : initSystem, endSystem, printSysMovie
   use Systemmodule, only : printSystem, getBravaisLattice, getLatticeConstant
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber,   &
                            getNumVacancies
   use SystemModule, only : getUniformGridParam
   use SystemModule, only : getSystemID, getNumAtomTypes
   use SystemModule, only : updateSystem, setRMSInfo
   use SystemModule, only : writeMomentDirectionData, writeAtomPositionData
   use SystemModule, only : writeMomentMovie, writeAtomPositionMovie
   use SystemModule, only : updateSystemMovie, resetSystemMovie
   use SystemModule, only : getLmaxKKR, getLmaxRho, getLmaxPot, getLmaxPhi
!
   use MediumHostModule, only : initMediumHost, endMediumHost, printMediumHost
!
   use LatticeModule, only : initLattice, endLattice, getLatticeType
!
   use ContourModule, only : initContour, endContour, getNumEs
!
   use ProcMappingModule, only : initProcMapping, endProcMapping,            &
                                 createParallelization, getNumEsOnMyProc,    &
                                 getNumKsOnMyProc, getNumRedundantEsOnMyProc,&
                                 getNumRedundantKsOnMyProc
!
   use SystemVolumeModule, only : initSystemVolume, endSystemVolume,         &
                                  printSystemVolume, updateSystemVolume,     &
                                  setSystemVolumeMT, getAtomicVPVolume
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,&
                                    calSymmetryFlags ! ,getSymmetryFlags
!
!  ===================================================================
!  initAtom2Proc and endAtom2Proc are called in main now!!!
!  ===================================================================
   use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc
   use Atom2ProcModule, only : getGlobalIndex, getAtom2ProcInGroup
   use Atom2ProcModule, only : printAtom2ProcTable, getLocalNumAtoms
!
   use AtomModule, only : initAtom, endAtom
   use AtomModule, only : getMaxLmax, printAtom, getStepFuncLmax
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax
   use AtomModule, only : getTruncPotLmax
   use AtomModule, only : getGridData, getAtomMuffinTinRad
   use AtomModule, only : getLocalAtomName, getLocalAtomicNumber
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
   use AtomModule, only : getLocalAtomNickName, printAtomMomentInfo
!   use AtomModule, only : setAtomVolMT, setAtomVolVP
!   use AtomModule, only : setAtomVolWS, setAtomVolINSC
   use AtomModule, only : getLocalAtomPosition, getLocalEvecOld
   use AtomModule, only : getLocalConstrainField
   use AtomModule, only : getMixingParam4Evec ! , resetCompFlags, getCompFlags
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use NeighborModule, only : initNeighbor, endNeighbor
   use SendRecvTmatModule, only : initSendRecvTmat, endSendRecvTmat, getNumProcWiseSends
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid
!
   use CoreStatesModule, only : initCoreStates, calCoreStates, endCoreStates
   use CoreStatesModule, only : readCoreStates, readCoreDensity, writeCoreDensity
   use CoreStatesModule, only : printCoreStates, printCoreDensity
!
   use BZoneModule, only : initBZone, printBZone, endBZone, getNumRotations, &
                           getNumKs
!
   use IBZRotationModule, only : initIBZRotation, endIBZRotation, computeRotationMatrix, &
                                 printIBZRotationMatrix, checkCrystalSymmetry
!
   use ValenceDensityModule, only : initValenceDensity, getFermiEnergy, endValenceDensity
   use ValenceDensityModule, only : printRho_L
   use ValenceDensityModule, only : setValenceVPCharge, setValenceVPMomentSize
   use ValenceDensityModule, only : getExchangeEnergy
!
   use GFMethodModule, only : initGFMethod, printValenceStates, &
                              endGFMethod
   use GFMethodModule, only : calValenceStates, calValenceDOS, writeSS_DOS, writeDOS, calPartialDOS
   use GFMethodModule, only : printExchangeParam
!
   use Uniform3DGridModule, only : initUniform3DGrid, endUniform3DGrid
   use Uniform3DGridModule, only : isUniform3DGridInitialized, createUniform3DGrid, printUniform3DGrid
   use Uniform3DGridModule, only : createProcessorMesh, insertAtomsInGrid, getUniform3DGrid
   use Uniform3DGridModule, only : distributeUniformGrid
!
   use IsoparametricIntegrationModule, only : initIsoparametricIntegration, &
                                              endIsoparametricIntegration
!
   use MadelungModule, only : initMadelung, endMadelung, printMadelungMatrix
!
   use ChargeDistributionModule, only : initChargeDistribution,     &
                                        endChargeDistribution,      &
                                        updateChargeDistribution,   &
                                        printChargeDistribution,    &
                                        getGlobalVPCellMomentTable
!
   use PotentialGenerationModule, only : initPotentialGeneration,  &
                                         endPotentialGeneration,   &
                                         computeNewPotential,      &
                                         printPotentialGeneration, &
                                         printMadelungShiftTable
   use PotentialGenerationModule, only : printNewPot_L => printPot_L
!
   use TotalEnergyModule, only : initTotalEnergy, endTotalEnergy,     &
                                 computeEnergyFunctional, printTotalEnergy, &
                                 getEnergyPerAtom, getPressurePerAtom
!
   use MixingModule, only : initMixing, endMixing, mixValues
!
   use ExchCorrFunctionalModule, only : initExchCorrFunctional,   &
                                        endExchCorrFunctional,    &
                                        isLDAFunctional,          &
                                        isGGAFunctional,          &
                                        isMGGAFunctional,         &
                                        isHybridFunctional
!
   use ConvergenceCheckModule, only : initConvergenceCheck, &
                                      endConvergenceCheck,  &
                                      checkConvergence
!
   use ConstrainLocalMomentModule, only : initConstrainLM,            &
                                          endConstrainLM,             &
                                          printConstrainMoment,       &
                                          calConstrainLM,             &
                                          getFieldRms,                &
                                          updateConstrainLM
!
   use SpinRotationModule, only: printSpinRotation
!
   use BookKeepingModule, only : initBookKeeping, endBookKeeping
!
!v2.0  use ExchangeInteractionModule, only : initExchangeInteraction, &
!v2.0                   endExchangeInteraction, calExchgInteractions, &
!v2.0                   printExchange
!
   use ChargeDensityModule, only : initChargeDensity,         &
                                   endChargeDensity,          &
                                   constructChargeDensity,    &
                                   isChargeComponentZero,     &
                                   printChargeDensity,        &
                                   printChargeDensity_L,      &
                                   printMomentDensity_L,      &
                                   updateTotalDensity
!
   use LdaCorrectionModule, only : initLdaCorrection, endLdaCorrection, &
                                   insertCorrectionOrbitals, getTransformMatrix, &
                                   checkLdaCorrection
!
   use SurfElementsModule, only : printSurfaceElements
!
   use ForceModule, only : initForce, endForce, calForce, printForce
!
   use KreinModule, only : initKrein, endKrein
!
   implicit none
!
   logical :: ScfConverged = .false.
   logical :: FileNamed = .false.
   logical :: isRmtExternal = .false.
   logical :: initSystemMovie=.true.
   logical :: IsoParamVINT = .false.
   logical :: FrozenCoreFileExist = .false.
   logical :: StandardInputExist = .false.
   logical :: isDOSCalculationOnly = .false.
!
   character (len=80) :: info_table, info_path
   character (len=160) :: itname, cmd
   character (len=12) :: anm
   character (len=8)  :: exec_date
   character (len=10) :: exec_time
   character (len=200) :: FileName
   character (len=50) :: StorageKey
   character (len=50) :: FrozenCoreFileName = 'FrozenCoreDensity.dat'
!
   integer (kind=IntKind) :: MyPE, NumPEs
   integer (kind=IntKind) :: funit, funit_sysmov, en_movie
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: i, id, ig, is, jl, nk, ne, n, na, ia
   integer (kind=IntKind) :: iscf, itstep, niter, sdstep_new
   integer (kind=IntKind) :: n_chgtab, n_madtab, n_potwrite, n_visual, n_sysmov
   integer (kind=IntKind) :: lmax_max, lmax_kkr_max, GlobalNumAtoms, jmax_pot
   integer (kind=IntKind) :: lmax_rho_max, lmax_pot_max, lmax_gaunt
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind) :: NumRotations
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: node_print_level
!  integer (kind=IntKind) :: ng_uniform(3)
   integer (kind=IntKind) :: InitMode
   integer (kind=IntKind), pointer :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: lmax_tmp(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: lmax_green(:)
   integer (kind=IntKind), allocatable :: lmax_mad(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind) :: NumMix(1)
   integer (kind=IntKind) :: MaxVal_Integer(2)
!  Needed for L-SIC:
   integer (kind=IntKind) :: lsic_mode
   integer (kind=IntKind) :: vna, vnb, vnc
!
   real (kind=RealKind), pointer :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: LocalNumValenceElectrons(:)
   real (kind=RealKind), allocatable :: LocalAtomPosi(:,:)
   real (kind=RealKind), allocatable :: LocalEvec(:,:)
   real (kind=RealKind), allocatable :: radius(:)
!
   real (kind=RealKind) :: rho_rms_av(2), pot_rms_av(2)
   real (kind=RealKind), allocatable :: rho_rms(:,:), pot_rms(:,:)
   real (kind=RealKind), allocatable :: evec_rms(:), bcon_rms(:)
   real (kind=RealKind), allocatable :: bcon_sd(:,:)
   real (kind=RealKind) :: ef_diff, tote_diff
   real (kind=RealKind) :: max_rms(8), keep_rms(4), rms, rms_sys(4)
!
   real (kind=RealKind) :: vcell(3,3), vorigin(3)
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind) :: alat
   real (kind=RealKind) :: rmt, rinsc, rend, rws
   real (kind=RealKind) :: Efermi, volume, cfac
   real (kind=RealKind) :: v0, val, evb
   real (kind=RealKind) :: t0, t1, t2, t3
!
   real (kind=RealKind), pointer :: mom_table(:)
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!  real (kind=RealKind), parameter :: xstart = -.913096740000D+01
!
   complex (kind=CmplxKind), pointer :: pot_l(:)
!
!  ===================================================================
!  Mixing quantity
!  ===================================================================
   type (MixListRealStruct), target :: RealArrayList
   type (MixListCmplxStruct), target :: CmplxArrayList
!
!  type (UniformGridStruct), pointer:: fftgp
!
   interface
!
!      subroutine setupLizNeighbor(print_level,isScreened)
      subroutine setupLizNeighbor(print_level)
         use KindParamModule, only : IntKind
         integer (kind=IntKind), intent(in) :: print_level(*)
!         logical, optional,intent(in)       :: isScreened
      end subroutine setupLizNeighbor
!
      subroutine setupMixRealArrayList(NLA, nsp, RAList, r_rms, p_rms)
         use KindParamModule, only : IntKind, RealKind
         use PublicTypeDefinitionsModule, only : MixListRealStruct
         implicit none
         integer (kind=IntKind), intent(in) :: NLA,nsp
         real (kind=RealKind), intent(in) :: r_rms(:,:)
         real (kind=RealKind), intent(in) :: p_rms(:,:)
         type (MixListRealStruct), target :: RAList
      end subroutine setupMixRealArrayList
!
      subroutine setupMixCmplxArrayList(NLA, nsp, CAList, r_rms, p_rms)
         use KindParamModule, only : IntKind, RealKind
         use PublicTypeDefinitionsModule, only : MixListCmplxStruct
         implicit none
         integer (kind=IntKind), intent(in) :: NLA,nsp
         real (kind=RealKind), intent(in) :: r_rms(:,:)
         real (kind=RealKind), intent(in) :: p_rms(:,:)
         type (MixListCmplxStruct), target :: CAList
      end subroutine setupMixCmplxArrayList
!
      subroutine updateMixRealValues( NLA, nsp, RAList )
         use KindParamModule, only : IntKind, RealKind
         use PublicTypeDefinitionsModule, only : MixListRealStruct
         implicit none
         integer (kind=IntKind), intent(in) :: NLA,nsp
         type (MixListRealStruct), target :: RAList
      end subroutine updateMixRealValues
!
      subroutine updateMixCmplxValues( NLA, nsp, CAList)
         use KindParamModule, only : IntKind, RealKind
         use PublicTypeDefinitionsModule, only : MixListCmplxStruct
         implicit none
         integer (kind=IntKind), intent(in) :: NLA,nsp
         type (MixListCmplxStruct), target :: CAList
      end subroutine updateMixCmplxValues
!
      function getValBandEnergy(na,ia,is) result(be)
         use KindParamModule, only : IntKind, RealKind
         integer (kind=IntKind), intent(in) :: na, ia, is
         real (kind=RealKind) :: be
      end function getValBandEnergy
!
      function getSphValDensity(na,ia,is) result(rho)
         use KindParamModule, only : IntKind, RealKind
         integer (kind=IntKind), intent(in) :: na, ia, is
         real (kind=RealKind), pointer :: rho(:)
      end function getSphValDensity
!
   end interface
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
   call date_and_time(exec_date,exec_time)
!
   if (MyPE == 0) then
      write(6,'(13a,i5,a)')'Execution starts at ',                            &
                   exec_time(1:2),':',exec_time(3:4),':',exec_time(5:6),', ', &
                   exec_date(5:6),'-',exec_date(7:8),'-',exec_date(1:4),      &
                   ' on ',NumPEs,' nodes'
!!    write(6,'(80(''-''),/)')
!!    write(6,'(/)')
!!    write(6,'(12x,a)')'********************************************************'
!!    write(6,'(12x,a)')'*                                                      *'
!!    write(6,'(12x,a)')'*    Full-Potential Multiple Scattering Theory Based   *'
!!    write(6,'(12x,a)')'*                                                      *'
!!    write(6,'(12x,a)')'*  Ab Initio Electronic Structure Calculation Package  *'
!!    write(6,'(12x,a)')'*                                                      *'
!!    write(6,'(12x,a)')'*                    Version 2.4.x.x                   *'
!!    write(6,'(12x,a)')'*                                                      *'
!!    write(6,'(12x,a)')'********************************************************'
!!    write(6,'(/,80(''=''))')
!     ----------------------------------------------------------------
      call print_version(6)
!     ----------------------------------------------------------------
   endif
!
!
!  -------------------------------------------------------------------
   call initGroupComm()
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
!  FileName = 'None'
!  inquire(unit=5,name=FileName,named=FileNamed)
   inquire(unit=5,name=FileName,exist=StandardInputExist)
!  write(6,*) "main:: Input file open: ",trim(FileName)
!  if (FileNamed) then
   if (StandardInputExist) then
!     ----------------------------------------------------------------
      call readInputData(5,def_id)
!     ----------------------------------------------------------------
   else ! The standard input can be taken from a file named as StandardInputFile
!     ----------------------------------------------------------------
      call openInputFile(7,StandardInputFile)
      call readInputData(7,def_id)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  Check if we are performing a self interaction corrected calculation
!  -------------------------------------------------------------------
   if (getKeyValue(def_id,'Local SIC',lsic_mode) /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('main','Input parameter is not found!')
!     ----------------------------------------------------------------
   endif
   if(lsic_mode.ne.0) then
     write(6,*) 'Local SIC mode ',lsic_mode
     call ErrorHandler('main','L-SIC not implemented yet! Be patient ...')
   end if
!
!  ===================================================================
!  get the name and path of the system information data table file
!  -------------------------------------------------------------------
   if (getKeyValue(def_id,'Current File Path',info_path) /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('main','Input parameter for Current File Path is not found!')
!     ----------------------------------------------------------------
   endif
!
   if (getKeyValue(def_id,'Info Table File Name',info_table) == 0) then
      itname = trim(info_path)//info_table
      info_id = getTableIndex(itname)
      if (info_id < 1) then
!        -------------------------------------------------------------
         call openInputFile(10,itname)
         call readInputData(10,info_id)
!        -------------------------------------------------------------
      endif
   else
      info_id = def_id
   endif
!
!  ===================================================================
!  initialize SCF-calculation related data
!  -------------------------------------------------------------------
   call initScfData(def_id)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize Potential Type
!  -------------------------------------------------------------------
   call initPotentialType(getPotentialTypeParam())
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize system (unit cell and related data)
!  ===================================================================
   call initSystem(def_id)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Check data consistency
!  ===================================================================
   if (getLmaxKKR() /= getLmaxPhi()) then
      call WarningHandler('main',                                  &
           'Recommended setting lmax_kkr = lmax_phi in info_table file')
   else if (isFullPotential()) then
      if (getLmaxRho() /= 2*getLmaxKKR()) then
         call WarningHandler('main',                                  &
              'Recommended setting lmax_rho = 2*lmax_kkr in info_table file')
      else if (getLmaxRho() /= getLmaxPot()) then
         call WarningHandler('main',                                  &
              'Recommended setting lmax_rho = lmax_pot in info_table file')
      endif
   endif
!
!  ===================================================================
!  Determine if using Z*Tau*Z - Z*J formula for the Green function
!  ===================================================================
   if (.not.isSSIrregularSolOn() .and. NumSS_IntEs < 2) then
      call ErrorHandler('main','No. of mesh points on real energy axis < 2',NumSS_IntEs)
   endif
!
!  ===================================================================
!  Initialize the contour in energy complex plane to find the total
!  number of energy mesh needed
!  ===================================================================
   if (isReadEmesh()) then
!     ----------------------------------------------------------------
      call initContour(getEmeshFileName(), istop, -1)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call initContour( ContourType, eGridType, NumEs, Temperature,   &
                        istop, -1, .true.)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   ne = getNumEs()
!  ne = max(NumSS_IntEs,getNumEs())
!  -------------------------------------------------------------------
!
!  ===================================================================
!  After the number of energy points is obtained, we need to call endContour
!  since initContour will be called again in GFMethod within SCF loop
!  -------------------------------------------------------------------
   call endContour()
!  -------------------------------------------------------------------
!
!  *******************************************************************
!
!  ===================================================================
   bravais(1:3,1:3)=getBravaisLattice()
!
   GlobalNumAtoms = getNumAtoms()
!
   AtomPosition => getAtomPosition()
   AtomicNumber => getAtomicNumber()
!
!  ===================================================================
!  Initialize the the lattice system
!  -------------------------------------------------------------------
   call initLattice(bravais)
!  -------------------------------------------------------------------
!
!  *******************************************************************
!
!  ===================================================================
!  Initialize the Brillouin zone mesh for k-space integration
!  ===================================================================
   if (isKKR() .or. isScreenKKR_LSMS() .or. isKKRCPA()) then
      if (isReadKmesh()) then
!        -------------------------------------------------------------
         call initBZone(getKmeshFileName(),istop,-1)
!        -------------------------------------------------------------
      else if (NumKMeshs > 0) then
!        -------------------------------------------------------------
         call initBZone(NumKMeshs,kGenScheme,Kdiv,Symmetrize,bravais, &
                        GlobalNumAtoms,AtomPosition,AtomicNumber,istop,-MyPE)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call WarningHandler('main','No K mesh is initialized')
!        -------------------------------------------------------------
      endif
      nk = getNumKs()
   else if (isFullPotential() .and. getLmaxKKR() > 2) then
      nk = (getLmaxKKR()+1)**2
   else
      nk = 1
   endif
!  ===================================================================
!
!  ===================================================================
!  Initialize the processes mapping module that determines how the
!  parallization will be performed
!  ===================================================================
   na = getNumAtoms()
!  -------------------------------------------------------------------
   call initProcMapping(na, ne, nk, isFullPotential(), istop, 0)
!  -------------------------------------------------------------------
#ifdef USE_SCALAPACK
!  -------------------------------------------------------------------
   call createParallelization()
!  -------------------------------------------------------------------
#else
!  -------------------------------------------------------------------
   call createParallelization(isLSMS())
!  -------------------------------------------------------------------
#endif
   call initParallelIO(getGroupID('Unit Cell'),1) ! only the 1st cluster
                                                  ! in the group performs
                                                  ! writing potential data
!  -------------------------------------------------------------------
   call initAtom2Proc(GlobalNumAtoms)
!  -------------------------------------------------------------------
!
   LocalNumAtoms=getLocalNumAtoms()
!
!  *******************************************************************
!
!  ===================================================================
!  initialize system volume and the volume of each atom in the system
!  -------------------------------------------------------------------
   call initSystemVolume()
!  -------------------------------------------------------------------
!
!  *******************************************************************
!
!  ===================================================================
!  set up print level
!  -------------------------------------------------------------------
!  call initOutput(def_id,-MyPE)
   call initOutput(def_id)
!  -------------------------------------------------------------------
   node_print_level = getStandardOutputLevel()
   call setErrorOutput(1)
!  -------------------------------------------------------------------
!
   allocate(atom_print_level(1:LocalNumAtoms))
   if (isOutputAtomBased()) then
      do i=1,LocalNumAtoms
         atom_print_level(i) = getStandardOutputLevel(i)
         node_print_level = max(atom_print_level(i), node_print_level)
      enddo
   else
      do i=1,LocalNumAtoms
         atom_print_level(i) = node_print_level
      enddo
   endif
!
   if (node_print_level >= 0 .and. .not.isStandardOutToScreen()) then
!     ----------------------------------------------------------------
      call getcwd(cmd)
!     ----------------------------------------------------------------
      write(6,'(a,a)')'Current working directory: ',trim(cmd)
!     ----------------------------------------------------------------
      call get_command(cmd)
!     ----------------------------------------------------------------
      write(6,'(a,a)')'Execution command: ',trim(cmd)
!     ----------------------------------------------------------------
!     write(6,'(a,a)') 'This code was compiled by ', compiler_version(), &
!                      ' using the options ', compiler_options()
!     ----------------------------------------------------------------
      n = command_argument_count()
!     ----------------------------------------------------------------
      do i = 1, n
!        -------------------------------------------------------------
         call get_command_argument(i,cmd)
!        -------------------------------------------------------------
         if (len_trim(cmd) > 0) then
            write(6,'(a,a)')'Command argument: ',trim(cmd)
         endif
      enddo
!     ----------------------------------------------------------------
      call get_environment_variable("HOSTNAME", cmd)
!     ----------------------------------------------------------------
      if (len(trim(cmd)) == 0) then
!        -------------------------------------------------------------
         call get_environment_variable("HOST", cmd)
!        -------------------------------------------------------------
      endif
      write(6,'(13a,i5,2a)')'Execution starts at ',                           &
                   exec_time(1:2),':',exec_time(3:4),':',exec_time(5:6),', ', &
                   exec_date(5:6),'-',exec_date(7:8),'-',exec_date(1:4),      &
                   ' on ',NumPEs,' nodes of ',trim(cmd)
!!   write(6,'(80(''-''),/)')
!!   write(6,'(/)')
!!   write(6,'(12x,a)')'********************************************************'
!!   write(6,'(12x,a)')'*                                                      *'
!!   write(6,'(12x,a)')'*    Full-Potential Multiple Scattering Theory Based   *'
!!   write(6,'(12x,a)')'*                                                      *'
!!   write(6,'(12x,a)')'*  Ab Initio Electronic Structure Calculation Package  *'
!!   write(6,'(12x,a)')'*                                                      *'
!!   write(6,'(12x,a)')'*                    Version 2.4                       *'
!!   write(6,'(12x,a)')'*                                                      *'
!!   write(6,'(12x,a)')'********************************************************'
!!   write(6,'(/)')
!    -----------------------------------------------------------------
     call print_version(6)
!    -----------------------------------------------------------------
     write(6,'(12x,a,i5)')'Number of atoms on each processor:              ', &
                          LocalNumAtoms
     write(6,'(12x,a,i5)')'Number of k-points on each processor:           ', &
                          getNumKsOnMyProc()
     write(6,'(12x,a,i5)')'Number of redundant k-points on each processor: ', &
                          getNumRedundantKsOnMyProc()
     write(6,'(12x,a,i5)')'Number of energies on each processor:           ', &
                          getNumEsOnMyProc()
     write(6,'(12x,a,i5)')'Number of redundant energies on each processor: ', &
                          getNumRedundantEsOnMyProc()
     write(6,'(/,80(''=''))')
!    -----------------------------------------------------------------
!    call force_openmp     ! use this only if necessary
     call print_threads(6)
!    -----------------------------------------------------------------
     if ( isKKR() ) then
        write(6,'(/,14x,a,/)')'::::  KKR Electronic Structure Calculation ::::'
     else if ( isLSMS() ) then
        write(6,'(/,14x,a,/)')'::::  LSMS Electronic Structure Calculation ::::'
     else if ( isScreenKKR() ) then
        write(6,'(/,14x,a,/)')'::::  Screend-KKR Electronic Structure Calculation ::::'
     else if ( isKKRCPA() ) then
        write(6,'(/,14x,a,/)')'::::  KKR-CPA Electronic Structure Calculation ::::'
     else if ( isScreenKKR_LSMS() ) then
        write(6,'(/,14x,a,/)')'::::  Screened-KKR-LSMS Electronic Structure Calculation ::::'
     else if ( isEmbeddedCluster() ) then
        write(6,'(/,14x,a,/)')'::::  Embedded-LSMS Electronic Structure Calculation ::::'
     else if ( isSingleSite() ) then
        write(6,'(/,14x,a,/)')'::::  Single Site Electronic Structure Calculation ::::'
     endif
   endif
!
!  ===================================================================
!  call initAtom and setAtomData to setup Atom Module.................
!  -------------------------------------------------------------------
   call initAtom(info_id,istop,node_print_level)
!  -------------------------------------------------------------------
   do i = 1,LocalNumAtoms
      rmt   = getAtomMuffinTinRad(i)
      rinsc = getInscrSphRadius(i)
      if ( rinsc-rmt>Ten2m3 .and. rmt>TEN2m6 ) then
!        -------------------------------------------------------------
         call setSystemVolumeMT(i,rmt)
!        -------------------------------------------------------------
         isRmtExternal = .true.
      endif
   enddo
!
   if ( isRmtExternal ) then
!     ----------------------------------------------------------------
      call updateSystemVolume()
!     ----------------------------------------------------------------
   endif
!
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
   if (getDOSrunID() == PrintDOSswitchOff) then
      isDOSCalculationOnly = .false.
   else
      isDOSCalculationOnly = .true.
   endif
!
!  ===================================================================
!  initialize medium system if needed
!  ===================================================================
   if (isKKRCPA() .or. isEmbeddedCluster()) then
!     ----------------------------------------------------------------
      call initMediumHost(def_id)
!     ----------------------------------------------------------------
      if (NumPEs > 8) then
!        -------------------------------------------------------------
         call printMediumHost(0)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call printMediumHost(-1)
!        -------------------------------------------------------------
      endif
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  initialize LIZ
!  -------------------------------------------------------------------
   call initNeighbor(LocalNumAtoms,atom_print_level)
!  -------------------------------------------------------------------
!   if ( isScreenKKR_LSMS() ) then
!      call setupLizNeighbor(atom_print_level,isScreenKKR_LSMS())
!   else
      call setupLizNeighbor(atom_print_level)
!   endif
!  -------------------------------------------------------------------
   if ( (isKKR() .or. isKKRCPA()) .and. node_print_level >= 0) then
!     ----------------------------------------------------------------
      call printBZone()
!     ----------------------------------------------------------------
   endif
   call initSendRecvTmat(LocalNumAtoms,node_print_level)
   n = getNumProcWiseSends()
   if (n > getMaxWaits()) then
!     ----------------------------------------------------------------
      call setMaxWaits(n)
!     ----------------------------------------------------------------
   endif
!  ===================================================================
!
   allocate(LocalAtomPosi(3,LocalNumAtoms))
   allocate(GlobalIndex(LocalNumAtoms), LocalEvec(3,LocalNumAtoms))
!
   do id=1,LocalNumAtoms
      LocalAtomPosi(1:3,id)=getLocalAtomPosition(id)
      LocalEvec(1:3,id)=getLocalEvecOld(id)
      GlobalIndex(id)=getGlobalIndex(id)
   enddo
!
!
!  *******************************************************************
!
!  ===================================================================
!  initialize atomic polyhedra
!  -------------------------------------------------------------------
!   call initPolyhedra(LocalNumAtoms,bravais,istop,node_print_level)
!  -------------------------------------------------------------------
!
!  *******************************************************************
!
!  ===================================================================
!  initialize radial grid
   allocate(radius(LocalNumAtoms))
   radius = TEN2m6
!  -------------------------------------------------------------------
   call initRadialGrid(LocalNumAtoms, istop, node_print_level)
!  -------------------------------------------------------------------
   do i=1,LocalNumAtoms
      ig=GlobalIndex(i)
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
!
      rend =  getOutscrSphRadius(i)
!
      if (isMuffinTinPotential() .or. isMuffinTinTestPotential()) then
         rmt = getAtomMuffinTinRad(i)
         rinsc = getInscrSphRadius(i)
         if ( rmt < 0.010d0 ) then
            rmt = rinsc
         endif
         rws = getWignerSeitzRadius(i)
!         volume = PI4*THIRD*(rws**3)
!         call setAtomVolWS(i,volume)
!         volume = PI4*THIRD*(rmt**3)
!         call setAtomVolMT(i,volume)
!         volume = PI4*THIRD*(rinsc**3)
!         call setAtomVolINSC(i,volume)
!         volume = getAtomicVPVolume(i)
!         call setAtomVolVP(i,volume)
!         call genRadialGrid( i, rinsc, rinsc, rws, ndivin, ndivout, nmult)
!
!
!        ========================================================
!???     Temporary fix due to inconsistent definitions of pzz/pzj
!???     between the relativistic and non-relativistic solvers
!        ========================================================
         if (getSingleSiteSolverType()==1) then
            rend=rws
         endif
!        ========================================================
!        -------------------------------------------------------------
         call genRadialGrid( i, xstart, rmt, rinsc, rws, rend, ndivin)
!        -------------------------------------------------------------
      else if (isASAPotential() ) then
!        rend =  getWignerSeitzRadius(i)
!        rmt = getAtomMuffinTinRad(i)
         rinsc = getWignerSeitzRadius(i)
         rmt = rinsc
!        if ( rmt < 0.010d0 ) then
!           rmt = rinsc
!        endif
         ndivin = ndivin+11
!        -------------------------------------------------------------
         call genRadialGrid( i, xstart, rmt, rinsc, rinsc, rend, ndivin )
!        -------------------------------------------------------------
      else if ( isMuffinTinASAPotential() ) then
         rend =  getWignerSeitzRadius(i)
         rmt = getAtomMuffinTinRad(i)
         rinsc = getWignerSeitzRadius(i)
         if ( rmt < 0.010d0 ) then
            rmt = rinsc
         endif
!         volume = PI4*THIRD*(rend**3)
!         call setAtomVolWS(i,volume)
!         volume = PI4*THIRD*(rmt**3)
!         call setAtomVolMT(i,volume)
!         volume = PI4*THIRD*(rinsc**3)
!         call setAtomVolINSC(i,volume)
!         volume = getAtomicVPVolume(i)
!         call setAtomVolVP(i,volume)
!         call genRadialGrid( i, rmt, rinsc, rinsc, ndivin, ndivout, nmult)
!        -------------------------------------------------------------
         call genRadialGrid( i, xstart, rmt, rinsc, rend, rend, ndivin )
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
!         volume = PI4*THIRD*(rws**3)
!         call setAtomVolWS(i,volume)
!         volume = PI4*THIRD*(rmt**3)
!         call setAtomVolMT(i,volume)
!         volume = PI4*THIRD*(rinsc**3)
!         call setAtomVolINSC(i,volume)
!         volume = getAtomicVPVolume(i)
!         call setAtomVolVP(i,volume)
!        -------------------------------------------------------------
!!!       call genRadialGrid( i, xstart, rmt, rinsc, rws, rend, ndivin )
         call genRadialGrid( i, rmt, rinsc, rws, rend, ndivin, ndivout, nmult)
!        -------------------------------------------------------------
      endif
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printRadialGrid(i)
!        -------------------------------------------------------------
      endif
      radius(i) = getOutscrSphRadius(i)
   enddo
!
!  ===================================================================
!  initialize and setup the Visual grid
!  -------------------------------------------------------------------
   if ( getDensityPrintFlag() >= 1 .or. getPotentialPrintFlag() >= 1) then
      if (.not.isUniform3DGridInitialized()) then
         call initUniform3DGrid(istop, node_print_level)
      endif
!     ----------------------------------------------------------------
      call fetchVisualDomainParameters(vna, vnb, vnc, vcell, vorigin)
!     ----------------------------------------------------------------
      call createUniform3DGrid('Visual',vna, vnb, vnc, vcell, cell_origin=vorigin)
!     ----------------------------------------------------------------
      call distributeUniformGrid('Visual')
      call insertAtomsInGrid('Visual', LocalNumAtoms, LocalAtomPosi,  &
                             getPointLocationFlag, radius)
!     ----------------------------------------------------------------
      if (node_print_level > 0) then
         call printUniform3DGrid('Visual')
      endif
   endif
   deallocate(radius)
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  setup the Lmax values.
!  ===================================================================
!  -------------------------------------------------------------------
   lmax_max=getMaxLmax()
!  -------------------------------------------------------------------
   allocate(lmax_pot(LocalNumAtoms))
   allocate(lmax_rho(LocalNumAtoms))
   allocate(lmax_kkr(LocalNumAtoms))
   allocate(lmax_phi(LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms))
   allocate(lmax_green(LocalNumAtoms))
   allocate(lmax_mad(LocalNumAtoms))
   allocate(lmax_tmp(LocalNumAtoms))
   lmax_kkr_max = 0
   lmax_rho_max = 0
   lmax_pot_max = 0
   do i=1,LocalNumAtoms
      lmax_kkr(i)   = getKKRLmax(i)
      lmax_phi(i)   = getPhiLmax(i)
      lmax_rho(i)   = getRhoLmax(i)
      lmax_pot(i)   = getPotLmax(i)
!
!     lmax_step(i)  = max(2*lmax_phi(i),getStepFuncLmax(i)) ! This change is
      lmax_step(i)  = getStepFuncLmax(i)                    ! made on 4/12/2015
!
      lmax_green(i) = min(2*lmax_phi(i),lmax_rho(i))
      lmax_mad(i) = max(lmax_rho(i),lmax_pot(i))
!new  lmax_max = max( lmax_max, lmax_step(i)+lmax_pot(i), lmax_rho(i), lmax_pot(i),&
      lmax_max = max( lmax_max, getTruncPotLmax(i), lmax_rho(i), lmax_pot(i), &
                      lmax_phi(i), lmax_kkr(i), lmax_green(i) )
      lmax_kkr_max  = max(lmax_kkr_max,lmax_kkr(i))
      lmax_rho_max  = max(lmax_rho_max,lmax_rho(i))
      lmax_pot_max  = max(lmax_pot_max,lmax_pot(i))
   enddo
!  ===================================================================
   MaxVal_Integer(1) = lmax_rho_max
   MaxVal_Integer(2) = lmax_pot_max
!  -------------------------------------------------------------------
   call GlobalMax(MaxVal_Integer,2)
!  -------------------------------------------------------------------
   lmax_rho_max = MaxVal_Integer(1)
   lmax_pot_max = MaxVal_Integer(2)
   lmax_gaunt = max(lmax_max,2*lmax_pot_max,2*lmax_rho_max)
!  ===================================================================
!  setup factors
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_gaunt)
!  -------------------------------------------------------------------
   call initGauntFactors(lmax_gaunt, istop, node_print_level)
!  -------------------------------------------------------------------
   if (IsoParamVINT) then
!     ----------------------------------------------------------------
      call initIsoparametricIntegration(lmax_max)
!     ----------------------------------------------------------------
   endif
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
!
!  -------------------------------------------------------------------
!  call initStepFunction(LocalNumAtoms,lmax_step,istop,node_print_level)
   call initStepFunction(LocalNumAtoms,lmax_max,lmax_step,ngr,ngt,istop,node_print_level)
!  -------------------------------------------------------------------
   do i=1,LocalNumAtoms
      lmax_tmp(i) = 2*max(lmax_rho(i), lmax_pot(i))
   enddo
   deallocate( ngr, ngt )
!
   do i=1,LocalNumAtoms
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printStepFunction(i)
!        -------------------------------------------------------------
         call testStepFunction(i)
!        -------------------------------------------------------------
      endif
   enddo
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  setup Madelung matrix
!  ===================================================================
   t2 = getTime()
!  -------------------------------------------------------------------
   call initMadelung(LocalNumAtoms,GlobalNumAtoms,GlobalIndex,        &
                     lmax_rho_max,lmax_pot_max,bravais,AtomPosition,  &
                     node_print_level)
!  -------------------------------------------------------------------
   if (MyPE == 0) then
      write(6,'(/,a,f10.5,/)')'Time:: initMadelung: ',getTime()-t2
   endif
!
   if (node_print_level >= 0) then
!     ----------------------------------------------------------------
      call printMadelungMatrix()
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call initSystemSymmetry( GlobalNumAtoms, LocalNumAtoms,            &
                            lmax_mad, lmax_step, atom_print_level )
!  -------------------------------------------------------------------
   call calSymmetryFlags()
!  -------------------------------------------------------------------
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  initialize potential module
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot, lmax_step,      &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
!
   allocate(LocalNumValenceElectrons(LocalNumAtoms))
!
   if (isTestPotential()) then
!     ----------------------------------------------------------------
      call initTestPotential(LocalNumAtoms,n_spin_pola)
!     ----------------------------------------------------------------
      call readTestPotential(lmax_pot)
!     ----------------------------------------------------------------
      do i = 1,LocalNumAtoms
         jmax_pot = (lmax_pot(i)+1)*(lmax_pot(i)+2)/2
         val = ZERO
         do is = 1,n_spin_pola
            do jl=1,jmax_pot
!              -------------------------------------------------------
               pot_l => getTestPotential(i,is,jl)
!              -------------------------------------------------------
               call setPotential(i,1,is,jl,pot_l,1)
!              -------------------------------------------------------
            enddo
!           ----------------------------------------------------------
            val = val + getTestValenceNum(i,is)
!           ----------------------------------------------------------
         enddo
         LocalNumValenceElectrons(i) = val
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
!     ================================================================
!     ----------------------------------------------------------------
      call readPotential()
!     ----------------------------------------------------------------
      Efermi=getPotEf()
!
!     ================================================================
!     initialize core states module
!     ----------------------------------------------------------------
      call initCoreStates(LocalNumAtoms,ErBottom,n_spin_pola,         &
                          isNonRelativisticCore(),istop,node_print_level)
!     ----------------------------------------------------------------
!
!     ================================================================
!     read core states and density data of local atoms from file
!     ================================================================
      call readCoreStates()
!     ----------------------------------------------------------------
      do i = 1,LocalNumAtoms
         LocalNumValenceElectrons(i) = ZERO
         do ia = 1, getLocalNumSpecies(i)
            LocalNumValenceElectrons(i) = LocalNumValenceElectrons(i)  &
                 + getLocalSpeciesContent(i,ia)*getZval( getLocalAtomicNumber(i,ia) )
         enddo
      enddo
!
      if (isFrozenCore(fcf_name=FrozenCoreFileName,fcf_exist=FrozenCoreFileExist)) then
         if (FrozenCoreFileExist) then
!           ----------------------------------------------------------
            call readCoreDensity(FrozenCoreFileName)
!           ----------------------------------------------------------
         endif
      endif
   endif
!  -------------------------------------------------------------------
   call syncAllPEs()
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   if (isKKR() .or. isScreenKKR_LSMS() .or. isKKRCPA()) then
!     ================================================================
!     initialize IBZ rotation matrix module
!     ----------------------------------------------------------------
      call initIBZRotation(isRelativisticValence(),getLatticeType(),  &
                           lmax_kkr_max,Symmetrize)
      call computeRotationMatrix(bravais,GlobalNumAtoms,AtomPosition,anum=AtomicNumber)
!     ----------------------------------------------------------------
      if( checkCrystalSymmetry(bravais,GlobalNumAtoms,AtomPosition,   &
                               anum=AtomicNumber) ) then
         if (node_print_level >= 0) then
            write(6,'(/,a,/)')'The crystal system has the point group symmetry!'
         endif
      else
         call ErrorHandler('main',                                    &
                           'The crystal system does not have the point group symmetry')
      endif
      if (node_print_level >= 0) then
!        -------------------------------------------------------------
         call printIBZRotationMatrix(Rot3D_Only=.true.)
!        -------------------------------------------------------------
      endif
!
!     if ( Symmetrize<0 .or. Symmetrize==1 ) then
!        NumRotations = getNumRotations()
!        -------------------------------------------------------------
!        call initSpecKInteg(lmax_kkr_max, NumRotations)
!        -------------------------------------------------------------
!     endif
!
!     ================================================================
!     Llody Module
!     ----------------------------------------------------------------
      call initKrein(LocalNumAtoms, n_spin_cant, n_spin_pola, lmax_kkr, &
                     getNumKs(), istop)
!     ----------------------------------------------------------------
   endif
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   if (abs(ErTop) > TEN2m8 .and. getDOSrunID() == 0) then
      Efermi = ErTop
      if (node_print_level >= 0) then
         write(6,'(/,a,f12.8)')' Initial Fermi energy is set by the input:',Efermi
      endif
   else
      if (node_print_level >= 0) then
         write(6,'(/,a,f12.8)')' Initial Fermi energy read from the potential:',Efermi
      endif
   endif
!  ===================================================================
!  initialize exchange-correlation scheme, potential generation, total
!  energy calculation, and charge distribution book-keeping
!  -------------------------------------------------------------------
   call initExchCorrFunctional(n_spin_pola,excorr_name,node_print_level)
!  -------------------------------------------------------------------
   if (node_print_level >= 0) then
      if (isLDAFunctional()) then
         write(6,'(/,2x,a)')'========================================='
         write(6,'(2x,a)')  'Exchange-Correlation Functional Type: LDA'
         write(6,'(2x,a,/)')'========================================='
      else if (isGGAFunctional()) then
         write(6,'(/,2x,a)')'========================================='
         write(6,'(2x,a)')  'Exchange-Correlation Functional Type: GGA'
         write(6,'(2x,a,/)')'========================================='
      else if (isHybridFunctional()) then
         write(6,'(/,2x,a)')'============================================'
         write(6,'(2x,a)')  'Exchange-Correlation Functional Type: Hybrid'
         write(6,'(2x,a,/)')'============================================'
      else if (isMGGAFunctional()) then
         write(6,'(/,2x,a)')'=========================================='
         write(6,'(2x,a)')  'Exchange-Correlation Functional Type: MGGA'
         write(6,'(2x,a,/)')'=========================================='
      else
         call ErrorHandler('main','Unknown exchange-correlation functional type')
      endif
   endif
!  ===================================================================
!
!  *******************************************************************
!
!  -------------------------------------------------------------------
   call initValenceDensity(LocalNumAtoms,LocalNumValenceElectrons,    &
                           lmax_rho,n_spin_pola,n_spin_cant,Efermi,   &
                           istop,atom_print_level,                    &
                           isGGA = isGGAFunctional().and.(.not.isDOSCalculationOnly))
!  -------------------------------------------------------------------
   call initGFMethod(LocalNumAtoms,GlobalIndex,LocalAtomPosi,         &
                     lmax_kkr,lmax_phi,lmax_pot,lmax_step,lmax_green, &
                     n_spin_pola,n_spin_cant,istop,atom_print_level,  &
                     isGGA = isGGAFunctional().and.(.not.isDOSCalculationOnly))
!  -------------------------------------------------------------------
   call initChargeDensity(LocalNumAtoms, GlobalIndex, lmax_rho,       &
                          n_spin_pola, n_spin_cant, atom_print_level, &
                          isGGA = isGGAFunctional().and.(.not.isDOSCalculationOnly))
!  -------------------------------------------------------------------
   call initChargeDistribution(LocalNumAtoms,GlobalNumAtoms,n_spin_pola)
!  -------------------------------------------------------------------
   call initPotentialGeneration(LocalNumAtoms,GlobalNumAtoms,lmax_pot,&
                                lmax_rho,n_spin_pola,istop,atom_print_level,&
                                isGGA = isGGAFunctional().and.(.not.isDOSCalculationOnly))
!  -------------------------------------------------------------------
   call initTotalEnergy(LocalNumAtoms,GlobalNumAtoms,getNumVacancies(),&
                        n_spin_pola,istop,atom_print_level,            &
                        isGGA = isGGAFunctional().and.(.not.isDOSCalculationOnly))
!  -------------------------------------------------------------------
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  Setup the mixing
!  ===================================================================
   NumMix(1) = 0
   do i = 1,LocalNumAtoms
      NumMix(1) = NumMix(1) + getLocalNumSpecies(i)
   enddo
   if ( .not.isFullPotential() ) then
!     ----------------------------------------------------------------
      call initMixing( 1, NumMix, RealArrayList )
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call initMixing( 1, NumMix, CmplxArrayList )
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call setupMixingScheme(LocalNumAtoms,n_spin_pola)
!  -------------------------------------------------------------------
   call initConstrainLM( LocalNumAtoms, n_spin_cant, n_spin_pola,     &
                         tstep, atom_print_level )
!  -------------------------------------------------------------------
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  If LDA correction is needed, call initLdaCorrection
!  ===================================================================
   if (isLdaCorrectionNeeded()) then
!     ----------------------------------------------------------------
      call initLdaCorrection(LocalNumAtoms,GlobalIndex,               &
                             getLocalAtomicNumber,getLocalNumSpecies)
!     ----------------------------------------------------------------
      call insertCorrectionOrbitals(getUJfile())
!     ----------------------------------------------------------------
   endif
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   call initConvergenceCheck(LocalNumAtoms,n_spin_pola,atom_print_level)
!  -------------------------------------------------------------------
   n = 0
   do id = 1, LocalNumAtoms
      n = n + getLocalNumSpecies(id)
   enddo
   allocate( rho_rms(n_spin_pola,n) )
   allocate( pot_rms(n_spin_pola,n) )
   allocate( evec_rms(LocalNumAtoms) )
   allocate( bcon_rms(LocalNumAtoms) )
   allocate( bcon_sd(3,LocalNumAtoms) )
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  Setup the force calculation
!  ===================================================================
   call initForce(LocalNumAtoms,istop,node_print_level)
!  ===================================================================
!
!  ===================================================================
!  Start SCF iterations.
!  ===================================================================
   if (node_print_level >= 0) then
      write(6,'(/,80(''=''))')
      write(6,'(''Time:: Start-up   '',5x,'' :'',f12.5,''Sec'')') &
                    getTime()-t0
      write(6,'(80(''-''))')
   endif
!
   if (GlobalNumAtoms > 100) then
      funit = 21     ! if tooo many atoms, charge table is written to a file
   else
      funit = 6
   endif
   n_chgtab = 0
   n_madtab = 0
   n_potwrite = 0
   n_visual = 0
   n_sysmov = 0
!
   if (MyPE == 0) then
!     ----------------------------------------------------------------
      call initBookKeeping(info_path,getSystemID(),MyPE)
!     ----------------------------------------------------------------
   endif
!
   if ( printSysMovie==1 ) then
!     ----------------------------------------------------------------
      call updateSystem('LIZ Size',getAtom2ProcInGroup)
!     ----------------------------------------------------------------
   endif
!
   niter = 0
   SD_LOOP: do itstep = 1,ntstep
!     ================================================================
!     reset the Initialization for calculating charge density on uniform grids
!     This needs to be done every time the atom positions are modified
!     ================================================================
      InitMode = 0
!
      sdstep_new = 1
      rho_rms(:,:) = ZERO
      pot_rms(:,:) = ZERO
      evec_rms(:)  = ZERO
      bcon_rms(:)  = ZERO
      do i = 1,LocalNumAtoms
         bcon_sd(1:3,i) = getLocalConstrainField(i)
      enddo
      rms = ZERO
      if (node_print_level >= 0) then
         write(6,*) "SD_Loop::", itstep
         call printAtomMomentInfo()
      endif
!
      SCF_LOOP: do iscf = 1, nscf
!
         if (isDOSCalculationOnly) then
!           ----------------------------------------------------------
            call calValenceDOS()
!           ----------------------------------------------------------
            call writeDOS(inputpath,getSystemID(),getDOSrunID())
            call writeSS_DOS(inputpath,getSystemID(),getDOSrunID())
!           ----------------------------------------------------------
            exit SCF_LOOP
         endif
!
         if (iscf < 10) then
            write(anm,'(a,i1,a)')'scf',iscf,'pre_'
         else if (iscf < 100) then
            write(anm,'(a,i2,a)')'scf',iscf,'pre_'
         else if (iscf < 1000) then
            write(anm,'(a,i3,a)')'scf',iscf,'pre_'
         else
            write(anm,'(a,i4,a)')'scf',iscf,'pre_'
         endif
!
         if ( isFullPotential() ) then
            if ( getPotentialPrintFlag() >= 0 .and. node_print_level >= 1 ) then
               do id = 1,LocalNumAtoms
                  call printPot_L(id,aux_name=anm)
               enddo
            endif
         endif
!
         niter = niter+1
         if (node_print_level >= 0) then
            write(6,'(//,80(''*''))')
            write(6,'(28x,a,i3)')'SCF Iteration Number: ',iscf
            write(6,'(80(''*''),/)')
            call FlushFile(6)
         endif
!        -------------------------------------------------------------
         call resetSystemMovie()
!        -------------------------------------------------------------
!        =============================================================
!        sync all processors.
!        -------------------------------------------------------------
         call syncAllPEs()
!        -------------------------------------------------------------
!        =============================================================
!        calculate the core states.
!        =============================================================
         if (.not.isTestPotential() .and. .not.isFrozenCore(iter=iscf)) then
!           ----------------------------------------------------------
            call calCoreStates(evb)
!           ----------------------------------------------------------
            if (evb < ErBottom - TEN2m6) then
               ErBottom = evb
            endif
            do id = 1,LocalNumAtoms
               if (atom_print_level(id) >= 0) then
!                 ----------------------------------------------------
                  call printCoreStates(id)
!                 call printCoreDensity(id)
!                 ----------------------------------------------------
               endif
            enddo
         else if (isFrozenCore(iter=iscf)) then
            if (.not.isFrozenCore(iter=iscf-1) .and. .not.FrozenCoreFileExist) then
               call writeCoreDensity(FrozenCoreFileName)
               FrozenCoreFileExist = .true.
            endif
            if (node_print_level >= 0) then
               write(6,'(//,80(''=''))')
               write(6,'(18x,a)')'Frozen core: core states are not recalculated'
               write(6,'(80(''=''),/)')
               call FlushFile(6)
            endif
         endif
!        =============================================================
!        generate new evec direction and set constrain field to the
!        potential.
!        =============================================================
         if ( n_spin_cant == 2 ) then
!           ----------------------------------------------------------
            call calConstrainLM(itstep)
!           ----------------------------------------------------------
         endif
!
         t1 = getTime()
!
!        =============================================================
!        call calValenceStates to calculation the valence density of
!        states, electron density, and Fermi energy
!        =============================================================
         t3 = getTime()
!        -------------------------------------------------------------
         call calValenceStates()
!        -------------------------------------------------------------
         if (node_print_level >= 0) then
            write(6,'(''Time:: calValenceStates :'',f12.5,''Sec'')')getTime()-t3
!           ----------------------------------------------------------
            call printValenceStates()
!           ----------------------------------------------------------
         endif
!
!        =============================================================
!        update the constrain field, moment direction, reset
!        the potential to the local direction.
!        =============================================================
         if ( n_spin_cant == 2 ) then
!           ----------------------------------------------------------
            call updateConstrainLM()
!           ----------------------------------------------------------
         endif
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
!        Construct the total electron density (and moment density)
!        -------------------------------------------------------------
         call constructChargeDensity()
!        -------------------------------------------------------------
!
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
!        update the charge distribution table.
!        -------------------------------------------------------------
         call updateChargeDistribution(getExchangeEnergy)
!        -------------------------------------------------------------
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
!        Calculate the new DFT-LDA potential from the new density
!        =============================================================
#ifdef EPrint_MT
         if ( isFullPotential() ) then
!           ----------------------------------------------------------
            call computeNewPotential(isMT=.true.)
!           ----------------------------------------------------------
            call computeEnergyFunctional(isMT=.true.)
!           ----------------------------------------------------------
            if ( node_print_level >= 0 ) then
!              -------------------------------------------------------
               call printTotalEnergy(isMT=.true.)
!              -------------------------------------------------------
            endif
         endif
#endif
         t2 = getTime()
!        -------------------------------------------------------------
         call computeNewPotential()
!        -------------------------------------------------------------
         if (MyPE == 0) then
            write(6,'(/,a,f10.5,/)')'Time:: computeNewPotential: ',getTime()-t2
         endif
!
         if (iscf < 10) then
            write(anm,'(a,i1,a)')'scf',iscf,'post_'
         else if (iscf < 100) then
            write(anm,'(a,i2,a)')'scf',iscf,'post_'
         else if (iscf < 1000) then
            write(anm,'(a,i3,a)')'scf',iscf,'post_'
         else
            write(anm,'(a,i4,a)')'scf',iscf,'post_'
         endif
!
         if (node_print_level >= 1) then
            do id = 1,LocalNumAtoms
               call printRho_L(id,aux_name=anm)
            enddo
         endif
!
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
!        calculate the DFT total energy.
!        -------------------------------------------------------------
         call computeEnergyFunctional()
!        -------------------------------------------------------------
         if ( node_print_level >= 0 ) then
!           ----------------------------------------------------------
            call printTotalEnergy()
!           ----------------------------------------------------------
         endif
!
!        =============================================================
!        Calculate the force...
!        =============================================================
         if ( isFullPotential() ) then
!           ----------------------------------------------------------
            call calForce()
!           ----------------------------------------------------------
            if (node_print_level >= 0) then
!              -------------------------------------------------------
               call printForce()
!              -------------------------------------------------------
            endif
         endif
!        =============================================================
!        check for convergence
!        -------------------------------------------------------------
         call checkConvergence( rho_rms, pot_rms, ef_diff, tote_diff, &
                                evec_rms, bcon_rms, getFermiEnergy(), &
                                getEnergyPerAtom() )
!        -------------------------------------------------------------
         if ( node_print_level >= 0 .and. GlobalNumAtoms <= 64 ) then
!           ----------------------------------------------------------
            call printChargeDistribution(iscf,funit)
!           ----------------------------------------------------------
         endif
!
!        =============================================================
         max_rms(:) = ZERO
         rms_sys=ZERO
         n = 0
         do id = 1, LocalNumAtoms
            rho_rms_av = 0
            pot_rms_av = 0
            do ia = 1, getLocalNumSpecies(id)
               n = n + 1
               cfac = getLocalSpeciesContent(id,ia)
               rho_rms_av(1:n_spin_pola) = rho_rms_av(1:n_spin_pola) + &
                                           cfac*rho_rms(1:n_spin_pola,n)
               pot_rms_av(1:n_spin_pola) = pot_rms_av(1:n_spin_pola) + &
                                           cfac*pot_rms(1:n_spin_pola,n)
            enddo
            rms_sys(1) = maxval(rho_rms_av(1:n_spin_pola))
            rms_sys(2) = maxval(pot_rms_av(1:n_spin_pola))
            if ( n_spin_cant==2 ) then
               rms_sys(3:4) = getFieldRms(id)
            endif
            max_rms(1) = max(max_rms(1),rho_rms_av(1))
            max_rms(3) = max(max_rms(3),pot_rms_av(1))
            if (n_spin_pola==2) then
               max_rms(2) = max(max_rms(2),rho_rms_av(2))
               max_rms(4) = max(max_rms(4),pot_rms_av(2))
            endif
            max_rms(5) = max(max_rms(5),ef_diff)
            max_rms(6) = max(max_rms(6),tote_diff)
            call setRMSInfo(getGlobalIndex(id),rms_sys)
         enddo
         if ( n_spin_cant==2 ) then
            max_rms(7:8) = getFieldRms()
         endif
!        -------------------------------------------------------------
         call GlobalMax(max_rms,8)
!        -------------------------------------------------------------
         rms = maxval(max_rms)
         if ( (max(max_rms(1),max_rms(n_spin_pola)) <= ptol .and.     &
               max(max_rms(3),max_rms(2+n_spin_pola)) <= ptol .and.   &
               max_rms(5) <= eftol .and. max_rms(6) <= etol) .or.     &
              (max_rms(7) <= rmstol .and. itstep>1 &
               .and. n_spin_cant==2)) then
            if (node_print_level >= 0) then
               write(6,'(//,80(''=''))')
               write(6,'(a,78x,a)')'#','#'
               write(6,'(a,24x,a,24x,a)')'#','SCF Convergence is reached !!!','#'
               write(6,'(a,78x,a)')'#','#'
               write(6,'(80(''=''))')
               write(6,'(//,a)')                                         &
'     Iteration   RMS of Rho      RMS of Pot        Diff Ef        Diff Etot      Diff Evec'
               write(6,'(80(''-''))')
               write(6,'(5x,i5,6(4x,f12.8))')iscf,max(max_rms(1),max_rms(2)), &
                      max(max_rms(3),max_rms(4)),max_rms(5:7)
               write(6,'(80(''=''),/)')
            endif
            ScfConverged = .true.
         endif
!
         if ( printSysMovie==1 .and. ( (n_sysmov+1)==movie .or.       &
                                  iscf==nscf .or. ScfConverged) ) then
!           ==========================================================
!           ----------------------------------------------------------
            call updateSystem('RMS Info',getAtom2ProcInGroup)
!           ----------------------------------------------------------
            call updateSystem('Energy Pressure',getAtom2ProcInGroup)
!           ----------------------------------------------------------
            if ( isFullPotential() ) then
!              -------------------------------------------------------
               call updateSystem('Force Data',getAtom2ProcInGroup)
!              -------------------------------------------------------
            endif
!           ----------------------------------------------------------
            if (n_spin_cant == 2) then
!              -------------------------------------------------------
               call updateSystem('Old Moment Direction',getAtom2ProcInGroup)
               call updateSystem('Moment Direction',getAtom2ProcInGroup)
               call updateSystem('Constrain Field',getAtom2ProcInGroup)
!              -------------------------------------------------------
            endif
         endif
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
!        output potential and electron density
!        =============================================================
         n_potwrite = n_potwrite + 1
         if ( potwrite > 0 .and. (n_potwrite==potwrite .or. iscf==nscf &
                                  .or. ScfConverged) ) then
!           ----------------------------------------------------------
            call writePotential()
!           ----------------------------------------------------------
            n_potwrite = 0
         endif
!        output the charge density and magnetic moment (added by xianglin)
         if ( isFullPotential() ) then
            if ( getDensityPrintFlag()==0 .and. (iscf==nscf .or. &
                                                ScfConverged) ) then
               do id = 1,LocalNumAtoms
!                 ----------------------------------------------------------
                  call printChargeDensity_L(id,"TotalNew")
                  call printChargeDensity_L(id,"Valence")
!                 ----------------------------------------------------------
                  if ( n_spin_pola==2 ) then
                     call printMomentDensity_L(id,"Valence")
                     call printMomentDensity_L(id,"TotalNew")
                  endif
!                 ----------------------------------------------------------
               enddo
            endif
         endif
!
         if ( isFullPotential() ) then
!           if ( node_print_level >=0 .and. getDensityPrintFlag()>=0 ) then
            if ( getPotentialPrintFlag() >= 0 .and. node_print_level >= 1 ) then
               do id = 1,LocalNumAtoms
!                 call printNewPot_L(id,-1,"Coulomb",aux_name=anm)
                  call printNewPot_L(id,-1,"Total",aux_name=anm)
!                 call printNewPot_L(id,-1,"Exchg",aux_name=anm)
                  if ( atom_print_level(id) >= 0 ) then
!                    call printNewPot_L(id,-1,"Total",1,aux_name=anm)
!                    call printNewPot_L(id,-1,"Coulomb",1,aux_name=anm)
!                    call printNewPot_L(id,-1,"Exchg",1,aux_name=anm)
!                    ----------------------------------
!                    Components of Coulomb potential
!                    ----------------------------------
!                    call printNewPot_L(id,-1,"Tilda",aux_name=anm)
!                    call printNewPot_L(id,-1,"Madelung",aux_name=anm)
!                    if (.not.isMuffinTinFullPotential()) then
!                       call printNewPot_L(id,-1,"Pseudo",aux_name=anm)
!                    endif
!                    call printNewPot_L(id,-1,"Tilda",1,aux_name=anm)
!                    call printNewPot_L(id,-1,"Madelung",1,aux_name=anm)
!                    if (.not.isMuffinTinFullPotential()) then
!                       call printNewPot_L(id,-1,"Pseudo",1,aux_name=anm)
!                    endif
                  endif
               enddo
            endif
         endif
!
         n_visual = n_visual + 1
         if ( n_visual==movie .or. iscf==nscf .or. ScfConverged ) then
!
            t2 = getTime()
!
            if ( node_print_level >=1 ) then
               do id = 1,LocalNumAtoms
                  if ( atom_print_level(id) >= 0 ) then
                     call printChargeDensity_L(id,"Valence")
                     call printChargeDensity_L(id,"TotalNew")
                     if ( n_spin_pola==2 ) then
                        call printMomentDensity_L(id,"Valence")
                        call printMomentDensity_L(id,"TotalNew")
                     endif
                     if ( isFullPotential().and.atom_print_level(id)>=1) then
                        call printChargeDensity_L(id,"Valence",1)
                        if (.not.isMuffinTinFullPotential()) then
                           call printChargeDensity_L(id,"Pseudo")
                           call printChargeDensity_L(id,"Pseudo",1)
                        endif
                     endif
                  endif
               enddo
            endif
!
            if ( getDensityPrintFlag()>=1 .and. MyPE == 0 ) then
               call printDensityOnGrid(LocalNumAtoms,node_print_level)
            endif
!
            if ( getPotentialPrintFlag()>=1 .and. MyPE == 0 ) then
               call printPotentialOnGrid(LocalNumAtoms,node_print_level)
            endif
!
            n_visual = 0
            if ( node_print_level >= 0 ) then
               write(6,'(/,80(''=''))')
               write(6,'(''Time:: IO for VisualGrid :'',f12.5,''Sec'')') &
                     getTime() - t2
               write(6,'(80(''=''))')
            endif
         endif
!
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
!        setup the quantities for mixing.
!        =============================================================
         if ( .not. isFullPotential() ) then
!           ----------------------------------------------------------
            call setupMixRealArrayList( LocalNumAtoms, n_spin_pola,    &
                                     RealArrayList, rho_rms, pot_rms )
!           ----------------------------------------------------------
            call mixValues(RealArrayList)
!           ----------------------------------------------------------
            call updateMixRealValues( LocalNumAtoms, n_spin_pola,      &
                                      RealArrayList )
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call setupMixCmplxArrayList( LocalNumAtoms, n_spin_pola,   &
                                    CmplxArrayList, rho_rms, pot_rms )
!           ----------------------------------------------------------
            if ( .not.( (isSphericalInputFile() .or. isChargeMixing()) .and. iscf==1) ) then
               call mixValues(CmplxArrayList )
            endif
!           ----------------------------------------------------------
            call updateMixCmplxValues(LocalNumAtoms,n_spin_pola,CmplxArrayList)
!           ----------------------------------------------------------
         endif
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
         if ( isChargeMixing() ) then
!           ==========================================================
!           update the potential
!           ==========================================================
!           ----------------------------------------------------------
            call updatePotential(LocalNumAtoms,n_spin_pola)
!           ----------------------------------------------------------
!           ==========================================================
!           update the charge distribution table.
!           ==========================================================
!           ----------------------------------------------------------
            call updateChargeDistribution(getExchangeEnergy)
!           ----------------------------------------------------------
            call computeNewPotential()
!           ----------------------------------------------------------
         endif
!
         if ( node_print_level >= 0 ) then
!           ----------------------------------------------------------
            call printPotentialGeneration()
!           ----------------------------------------------------------
         endif
!
         n_madtab = n_madtab + 1
         if ( (printSysMovie==1 .and. n_madtab == movie) .or.         &
              iscf == nscf .or. ScfConverged) then
            call syncAllPEs()
            if ( node_print_level >= 0 ) then
!              -------------------------------------------------------
               call printMadelungShiftTable(iscf,funit)
!              -------------------------------------------------------
               if (n_spin_pola == 2) then
!                 ----------------------------------------------------
                  call printMomentVsCoreSplit(iscf,funit)
!                 ----------------------------------------------------
               endif
            endif
            n_madtab = 0
         endif
!
         n_chgtab = n_chgtab + 1
         if (n_chgtab==movie .or. iscf==nscf .or. ScfConverged) then
!            if (node_print_level >= 0) then
!              -------------------------------------------------------
!               call printChargeDistribution(iscf,funit)
!              -------------------------------------------------------
!            endif
!            if (n_spin_cant == 2) then
!              -------------------------------------------------------
!               mom_table => getGlobalVPCellMomentTable()
!              -------------------------------------------------------
!               call updateSystemMovie(mom_table)
!              -------------------------------------------------------
!              call writeMomentMovie(iscf,itstep)
!              -------------------------------------------------------
!            endif
            n_chgtab = 0
         endif
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
!        update potentials and desities
!        =============================================================
!        -------------------------------------------------------------
         call updatePotential(LocalNumAtoms,n_spin_pola)
!        -------------------------------------------------------------
         call updateTotalDensity(setValenceVPCharge, setValenceVPMomentSize)
!        -------------------------------------------------------------
!
         Efermi = getFermiEnergy()
         call setPotEf(Efermi)
!        =============================================================
!
!        *************************************************************
!
!        =============================================================
!        write data to the book-keeping file
!        -------------------------------------------------------------
         if ( MyPE == 0 ) then
            keep_rms(1) = max(max_rms(1),max_rms(2))
            keep_rms(2) = max(max_rms(3),max_rms(4))
            keep_rms(3:4) = max_rms(5:6)
!           ----------------------------------------------------------
            call keep(sdstep_new,niter,keep_rms,getFermiEnergy(),     &
                      getEnergyPerAtom(), getPressurePerAtom())
!           ----------------------------------------------------------
            sdstep_new = 0
         endif
#ifdef DEBUG
         write(6,*) "Main:: End SCF"
         call printAtomMomentInfo()
#endif
         t1 = getTime() - t1
         if ( node_print_level >= 0 ) then
            write(6,'(/,80(''=''))')
            write(6,'(''Time:: Iteration  '',i5,'' :'',f12.5,''Sec'')') &
                    iscf, t1
            write(6,'(''Time:: Cummulative'',5x,'' :'',f12.5,''Sec'')') &
                    getTime()-t0
            write(6,'(80(''=''))')
         endif
!
         if ( printSysMovie==1 ) then
            if (  MyPE==0 .and. initSystemMovie ) then
               call driverSystemMovie( initSystemMovie, getEnergyPerAtom(), &
                                       funit_sysmov, en_movie )
               initSystemMovie=.false.
            endif
            n_sysmov = n_sysmov+1
            if ( n_sysmov==movie .or. (iscf==nscf .or. ScfConverged) ) then
               if ( MyPE==0 ) then
                  call printSystemMovie(iscf,itstep,funit_sysmov,en_movie)
               endif
               n_sysmov = 0
            endif
         endif
!
         if ( ScfConverged .or. nscf == 1 ) then
            if ( .not.isRelativisticValence() ) then
!              changed by xianglin because calPartialDOS has not been implemented for REL
!              -------------------------------------------------------
               call calPartialDOS(getFermiEnergy())
!              -------------------------------------------------------
            endif
            exit SCF_LOOP
         endif
!
      enddo SCF_LOOP
!
      ScfConverged = .false.
      do i = 1,LocalNumAtoms
         bcon_sd(1:3,i) = getLocalConstrainField(i)-bcon_sd(1:3,i)
         bcon_rms(i) = sqrt( bcon_sd(1,i)**2 + &
                             bcon_sd(2,i)**2 + &
                             bcon_sd(3,i)**2 )
      enddo
      max_rms(8) = maxval(bcon_rms)
      call GlobalMax(max_rms(8))
      if ( max_rms(8) < rmstol .and. max_rms(8)/=ZERO ) then
         if (node_print_level >= 0) then
            write(6,'(/,a)')' SD-Convergence is reached!!!'
            write(6,'(/,80(''=''))')
         endif
         exit SD_LOOP
      endif
   enddo SD_LOOP
!
   if ( MyPE==0 .and. printSysMovie==1) then
      call driverSystemMovie( initSystemMovie, getEnergyPerAtom(), &
                              funit_sysmov, en_movie )
      initSystemMovie=.true.
   endif
!
   deallocate(LocalEvec)
   deallocate(LocalNumValenceElectrons)
!
   if (n_spin_cant == 2) then
!     ----------------------------------------------------------------
      call writeMomentDirectionData()
!     ----------------------------------------------------------------
   endif
!
   if ( isExchangeParamNeeded() .and. n_spin_cant == 2 ) then
      call printExchangeParam(trim(inputpath)//'OneSiteExchangeParam.dat', &
                              trim(inputpath)//'PairExchangeParam.dat')
   endif
!
   if ( isExchangeParamNeeded() .and. (isKKR() .or. isScreenKKR_LSMS() &
                                       .or. isKKRCPA()) ) then
      if ( isScreenKKR_LSMS() ) then
         call setScfMethod("ScreenKKR")
      endif
      if ( Symmetrize/=0 ) then
         Symmetrize = 0
      endif
!     ----------------------------------------------------------------
      call endBZone()
!     ----------------------------------------------------------------
      if (isReadKmesh()) then
!        -------------------------------------------------------------
         call initBZone(getKmeshFileName(),istop,node_print_level)
!        -------------------------------------------------------------
      else if (NumKMeshs > 0) then
!        -------------------------------------------------------------
         call initBZone(NumKMeshs,kGenScheme,Kdiv,Symmetrize,bravais, &
                       GlobalNumAtoms,AtomPosition,AtomicNumber,istop,&
                       node_print_level)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call WarningHandler('main','No K mesh is initialized')
!        -------------------------------------------------------------
      endif
!
! The following lines of code were from version 2.0 and are commented out
! In order to calcuculate the exchange interaction, ExchangeInteractionModule
! needs to be modified
!     ----------------------------------------------------------------
!v2.0 call initExchangeInteraction( LocalNumAtoms, lmax_kkr, lmax_phi,    &
!v2.0                               lmax_pot, LocalAtomPosi, n_spin_pola, &
!v2.0                               n_spin_cant, getFermiEnergy(),        &
!v2.0                               istop, atom_print_level )
!     ----------------------------------------------------------------
!v2.0 call calExchgInteractions()
!     ----------------------------------------------------------------
!v2.0 call printExchange()
!     ----------------------------------------------------------------
!v2.0 call endExchangeInteraction()
!     ----------------------------------------------------------------
stop 'Under construction...'
   endif
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
!  clean up the allocated spaces.
!  ===================================================================
   nullify(AtomPosition,AtomicNumber)
   deallocate( GlobalIndex, LocalAtomPosi, atom_print_level )
   deallocate( lmax_kkr, lmax_phi, lmax_pot, lmax_rho, lmax_green, lmax_mad)
   deallocate( lmax_step, rho_rms, pot_rms, lmax_tmp )
!
!  -------------------------------------------------------------------
   call endInput()
!  -------------------------------------------------------------------
   if ( MyPE== 0) then
!     ----------------------------------------------------------------
      call endBookKeeping()
!     ----------------------------------------------------------------
   endif
   if (node_print_level >= 0) then
      write(6,'(/,80(''=''))')
      write(6,'(''Time:: Job total  '',5x,'' :'',f12.5,''Sec'')') &
                    getTime()-t0
      write(6,'(80(''-''))')
   endif
!
   if (isKKR() .or. isScreenKKR_LSMS() .or. isKKRCPA()) then
!     ----------------------------------------------------------------
      call endBZone()
      call endIBZRotation()
!     ----------------------------------------------------------------
!     if ( Symmetrize==0 .or. Symmetrize==1 ) then
!        call endSpecKInteg()
!     endif
   endif
!
   if (.not.isTestPotential()) then
      call endCoreStates()
   endif
!
   if (node_print_level >= 0) then
      write(6,*) " Run Ending! "
      call FlushFile(6)
   endif
!
   if (isLdaCorrectionNeeded()) then
!     ----------------------------------------------------------------
      call endLdaCorrection()
!     ----------------------------------------------------------------
   endif
   call endConstrainLM()
   call endForce()
   call endConvergenceCheck()
   call endTotalEnergy()
   call endChargeDensity()
   call endChargeDistribution()
   call endMadelung()
   call endPotentialGeneration()
   call endExchCorrFunctional()
   call endValenceDensity()
   call endGFMethod()
   call endPotential()
   call endStepFunction()
   call endRadialGrid()
   if ( IsoParamVINT ) then
      call endIsoparametricIntegration()
   endif
   if (isUniform3DGridInitialized()) then
      call endUniform3DGrid()
   endif
   call endKrein()
   call endSystemSymmetry()
   call endLattice()
   call endGauntFactors()
   call endSphericalHarmonics()
   call endSendRecvTmat()
   call endNeighbor()
   call endPotentialType()
   call endAtom()
   call endSystemVolume()
   call endAtom2Proc()
   call endParallelIO()
   call endProcMapping()
   if (isKKRCPA() .or. isEmbeddedCluster()) then
!     ----------------------------------------------------------------
      call endMediumHost()
!     ----------------------------------------------------------------
   endif
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
      write(6,'(a)')'End of a successful run......!'
   endif
!
   call endOutput()
!
end program mst2
