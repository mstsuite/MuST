program kubo
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8, TEN2m6, TEN2m5, TEN2m3, ONE, PI4, THIRD
!
   use ChemElementModule, only : getZval
!
   use TimerModule, only : initTimer, getTime, fetchStoredTime
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
   use ProcMappingModule, only : initProcMapping, endProcMapping,            &
                                 createParallelization, getNumEsOnMyProc,    &
                                 getNumKsOnMyProc, getNumRedundantEsOnMyProc,&
                                 getNumRedundantKsOnMyProc
!
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, &
                                       getDataStorage, RealMark
!
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

   use KuboDataModule, only : initKuboData, includeVertexCorrections, endKuboData
   use KuboDataModule, only : printTildeMatrices
 
   use ScfDataModule, only : n_spin_pola, n_spin_cant 
   use ScfDataModule, only : ngaussr, ngaussq, isRelativisticValence
   use ScfDataModule, only : initScfData, endScfData, isLSMS, getPotentialTypeParam
   use ScfDataModule, only : NumEs, ContourType, eGridType
   use ScfDataModule, only : Temperature
   use ScfDataModule, only : istop, NumKMeshs, kGenScheme, Kdiv, Symmetrize
   use ScfDataModule, only : isKKRCPA, isKKRCPASRO, getSingleSiteSolverType
   use ScfDataModule, only : pole_step

   use LatticeModule, only : initLattice, endLattice, getLatticeType

   use ContourModule, only : initContour, endContour, getNumEs

   use PotentialTypeModule, only : initPotentialType, endPotentialType,      &    
                                   isASAPotential, isMuffinTinPotential,     &    
                                   isMuffinTinASAPotential, isTestPotential, &
                                   isMuffinTinTestPotential, isFullPotential,&
                                   printPotentialType, isSphericalPotential, &
                                   isMuffinTinFullPotential, setIsSphericalPotential

   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
   use PolyhedraModule, only : getWignerSeitzRadius, getNeighborDistance

   use SystemModule, only : initSystem, getBravaisLattice, getLatticeConstant
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber
   use SystemModule, only : endSystem, getLmaxKKR, getLmaxRho, getLmaxPot, getLmaxPhi
!
   use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc
   use Atom2ProcModule, only : getGlobalIndex, getAtom2ProcInGroup
   use Atom2ProcModule, only : printAtom2ProcTable, getLocalNumAtoms

   use BZoneModule, only : initBZone, printBZone, endBZone, getNumRotations, &
                           getNumKs
  
   use SystemVolumeModule, only : initSystemVolume, endSystemVolume,         &
                                  printSystemVolume, updateSystemVolume,     &
                                  setSystemVolumeMT, getAtomicVPVolume

   use AtomModule, only : initAtom, endAtom
   use AtomModule, only : getRadialGridData, getMuffinTinRadius, setMuffinTinRadius
   use AtomModule, only : getAtomCoreRad, setAtomCoreRad
   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies
   use AtomModule, only : getLocalAtomPosition, getLocalEvec
   use AtomModule, only : getMaxLmax, printAtom, getStepFuncLmax, getTruncPotLmax
   use AtomModule, only : getPotLmax, getKKRLmax, getPhiLmax, getRhoLmax

   use MediumHostModule, only : initMediumHost, endMediumHost, printMediumHost 
   use SendRecvTmatModule, only : initSendRecvTmat, endSendRecvTmat, getNumProcWiseSends
   use NeighborModule, only : initNeighbor, endNeighbor

   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
   use RadialGridModule, only : getRadialGridRadius
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : printStepFunction, testStepFunction
   use MadelungModule, only : initMadelung, endMadelung, printMadelungMatrix
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,&
                                    calSymmetryFlags
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential, writePotential
   use PotentialModule, only : getPotEf, setV0, printPot_L
   use PotentialModule, only : setPotential, setPotEf
   use PotentialModule, only : isSphericalInputFile

   use IBZRotationModule, only : initIBZRotation, endIBZRotation, computeRotationMatrix, &
                                 printIBZRotationMatrix, checkCrystalSymmetry

   use SSSolverModule, only : initSSSolver, endSSSolver
   use MSSolverModule, only : initMSSolver, endMSSolver
   use ConductivityModule, only : initConductivity, calConductivity, endConductivity
   use ConductivityModule, only : sigma, sigmatilde, sigmatilde2, sigmatilde3, sigmatilde4

   use WriteMatrixModule, only : writeMatrix
   use MatrixInverseModule, only : MtxInv_LU

   use SineMatrixZerosModule, only : initSineMatrixZeros, endSineMatrixZeros
   use SineMatrixZerosModule, only : findSineMatrixZeros, printSineMatrixZerosInfo, &
                                     getNumSineMatrixZeros

!  ===================================================================
   implicit none
!
   logical :: StandardInputExist = .false.
   logical :: isRmtUpdated = .false.
   logical :: vc = .true.

   character (len=80) :: info_table, info_path
   character (len=160) :: itname
   character (len=8)  :: exec_date
   character (len=10) :: exec_time
   character (len=200) :: FileName
!
   integer (kind=IntKind) :: MyPE, NumPEs
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: lmax_max, lmax_kkr_max
   integer (kind=IntKind) :: lmax_rho_max, lmax_pot_max, lmax_gaunt
   integer (kind=IntKind) :: GlobalNumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: i, j, k, ig, id, n, na, ne, nk, is, ia, ib
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind), pointer :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:), atom_print_level(:)  
   integer (kind=IntKind), allocatable :: lmax_tmp(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: lmax_green(:)
   integer (kind=IntKind), allocatable :: lmax_mad(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind), allocatable :: LocalNumSpecies(:)
   integer (kind=IntKind) :: MaxVal_Integer(2)
!
   real (kind=RealKind) :: Efermi, t0, t2, t3, t_inp, t_outp
   real (kind=RealKind) :: rmt, rinsc, rend, rws, hin, rmt_grid, rc, ei
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), pointer :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: final_sigma(:,:,:)
   real (kind=RealKind), allocatable :: LocalAtomPosi(:,:)
   real (kind=RealKind), allocatable :: LocalEvec(:,:)
   real (kind=RealKind), allocatable :: radius(:)
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!  real (kind=RealKind), parameter :: xstart = -.913096740000D+01

   interface
!
!      subroutine setupLizNeighbor(print_level,isScreened)
      subroutine setupLizNeighbor(print_level)
         use KindParamModule, only : IntKind
         integer (kind=IntKind), intent(in) :: print_level(*)
!         logical, optional,intent(in)       :: isScreened
      end subroutine setupLizNeighbor
!
   end interface
!
!  -------------------------------------------------------------------
   call initTimer()
!  -------------------------------------------------------------------
   t0 = getTime()
   t_inp = ZERO; t_outp = ZERO
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
      write(6,'(80(''-''),/)')
      write(6,'(/)')
      write(6,'(12x,a)')'********************************************************'
      write(6,'(12x,a)')'*                                                      *'
      write(6,'(12x,a)')'*                        KUBO                          *'
      write(6,'(12x,a)')'*                                                      *'
      write(6,'(12x,a)')'*    Muffin-Tin Multiple Scattering Theory Based       *'
      write(6,'(12x,a)')'*                                                      *'
      write(6,'(12x,a)')'*    Electrical Conductivity Calculation Package       *'
      write(6,'(12x,a)')'*                                                      *'
      write(6,'(12x,a)')'*                    Version 1.0.0                     *'
      write(6,'(12x,a)')'*                                                      *'
      write(6,'(12x,a)')'********************************************************'
      write(6,'(/,80(''=''))')
   endif
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
   inquire(unit=5,name=FileName,exist=StandardInputExist)
   t3 = getTime()
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
   t_inp = t_inp + (getTime() - t3)
!  ===================================================================
!  Enter KUBO code here

   call initScfData(def_id)
   call initKuboData(def_id)
   call initPotentialType(getPotentialTypeParam())
   call initSystem(def_id)

   t_inp = t_inp + fetchStoredTime()

   allocate(final_sigma(3,3,n_spin_pola))
   final_sigma = 0

   bravais(1:3,1:3)=getBravaisLattice()
!
   GlobalNumAtoms = getNumAtoms()
!
   AtomPosition => getAtomPosition()
   AtomicNumber => getAtomicNumber()
!
   call initContour( ContourType, eGridType, NumEs, Temperature,   &
                        istop, -1, .true.)
   ne = getNumEs()
   call endContour()

!  ===================================================================
!  Initialize the the lattice system
!  -------------------------------------------------------------------
   call initLattice(bravais)
!  -------------------------------------------------------------------
   call initBZone(NumKMeshs,kGenScheme,Kdiv,Symmetrize,bravais, &
                        GlobalNumAtoms,AtomPosition,AtomicNumber,istop,-MyPE)
!  -------------------------------------------------------------------
   nk = getNumKs()
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
!  ------------------------------------------------------------------
   call initAtom2Proc(GlobalNumAtoms)
!  ------------------------------------------------------------------
!
   LocalNumAtoms=getLocalNumAtoms()
!  ------------------------------------------------------------------   
   call initSystemVolume()
!  ------------------------------------------------------------------
   call initOutput(def_id)
   node_print_level = getStandardOutputLevel()
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

   call setErrorOutput(e_print_level=1,w_print_level=node_print_level,&
                       m_print_level=node_print_level)

   call initAtom(info_id,istop,node_print_level)

   do i = 1,LocalNumAtoms
      rmt   = getMuffinTinRadius(i)
      rinsc = getInscrSphRadius(i)
      if (rmt < ONE) then
         if (abs(rmt) < TEN2m6) then
!           ==========================================================
!           The muffin-tin radius is set to be the inscribed radius.
!           ==========================================================
            rmt = rinsc
         else
!           ==========================================================
!           In this case, rmt is treated as a scaling factor for rinsc
!           The muffin-tin radius is set to be the inscribed radius
!           multiplied by the scaling factor
!           ==========================================================
            rmt = rmt*rinsc
         endif
!        -------------------------------------------------------------
         call setMuffinTinRadius(i,rmt)
!        -------------------------------------------------------------
         if (.not.isFullPotential()) then
!           ==========================================================
!           For Muffin-tin, ASA, or Muffin-tin-ASA calculations, since
!           potential outside rmt is 0, core radius is set to rmt
!           ----------------------------------------------------------
            call setAtomCoreRad(i,rmt)
!           ----------------------------------------------------------
         endif
      endif
!
      if ( abs(rinsc-rmt) > TEN2m6) then
!        =============================================================
!        The muffin-tin radius is set to be other than the inscribed radius.
!        -------------------------------------------------------------
         call setSystemVolumeMT(i,rmt)
!        -------------------------------------------------------------
         isRmtUpdated = .true.
      endif
   enddo
   if ( isRmtUpdated ) then
!     ----------------------------------------------------------------
      call updateSystemVolume()
!     ----------------------------------------------------------------
   endif

   if (isKKRCPA() .or. isKKRCPASRO()) then
!     ----------------------------------------------------------------
      call initMediumHost(def_id)
!     ----------------------------------------------------------------
      call printMediumHost(0)
!     ----------------------------------------------------------------
   endif

!  ===================================================================
!  initialize LIZ
!  -------------------------------------------------------------------
   call initNeighbor(LocalNumAtoms,atom_print_level)
!  -------------------------------------------------------------------
    call setupLizNeighbor(atom_print_level)
    call initSendRecvTmat(LocalNumAtoms,node_print_level)
    n = getNumProcWiseSends()
    if (n > getMaxWaits()) then
!     ----------------------------------------------------------------
      call setMaxWaits(n)
!     ----------------------------------------------------------------
    endif
!   ===================================================================
!
    allocate(LocalAtomPosi(3,LocalNumAtoms))
    allocate(GlobalIndex(LocalNumAtoms), LocalEvec(3,LocalNumAtoms))
!
    do id=1,LocalNumAtoms
      LocalAtomPosi(1:3,id)=getLocalAtomPosition(id)
      LocalEvec(1:3,id)=getLocalEvec(id,'old')
      GlobalIndex(id)=getGlobalIndex(id)
    enddo

   allocate(radius(LocalNumAtoms))
   radius = TEN2m6
!  -------------------------------------------------------------------
   call initRadialGrid(LocalNumAtoms, istop, node_print_level)
!  -------------------------------------------------------------------
   isRmtUpdated = .false.

   do i=1,LocalNumAtoms
      ig=GlobalIndex(i)
!     ----------------------------------------------------------------
      call getRadialGridData(i,ndivin,ndivout,nmult,hin)
!     ----------------------------------------------------------------
      rend = max(getOutscrSphRadius(i),getAtomCoreRad(i))
      rmt = getMuffinTinRadius(i)
      if ( rmt < 0.010d0 ) then
         call ErrorHandler('main','rmt < 0.01',rmt)
      endif
!
      if (isMuffinTinPotential()) then
         rinsc = getInscrSphRadius(i)
         if (rmt - rinsc > TEN2m6) then
            call WarningHandler('main','rinsc < rmt, rmt is reset to rinsc',rinsc,rmt)
            rmt = rinsc
         endif
         rws = getWignerSeitzRadius(i)
!
!        ========================================================
!???     Temporary fix due to inconsistent definitions of pzz/pzj
!???     between the relativistic and non-relativistic solvers
!        ========================================================
         if (getSingleSiteSolverType()==1) then
            rend=rws
         endif
!        ========================================================
      else if (isASAPotential() ) then
         rend =  getWignerSeitzRadius(i)
         rinsc = getWignerSeitzRadius(i)
         rmt = rinsc
         ndivin = ndivin+29
      else if ( isMuffinTinASAPotential() ) then
         rend =  getWignerSeitzRadius(i)
         rinsc = getWignerSeitzRadius(i)
      else
         if (getNeighborDistance(i,1)-getOutscrSphRadius(i) < TEN2m8) then
!           ----------------------------------------------------------
            call WarningHandler('main',                               &
                     'Ill condition found: Neighbor distance <= Rcs', &
                     getNeighborDistance(i,1),getOutscrSphRadius(i))
!           ----------------------------------------------------------
         endif
         rinsc = getInscrSphRadius(i)
         rws = getWignerSeitzRadius(i)
      endif
      rc = getAtomCoreRad(i)
      if (hin > TEN2m6 .or. rc < ONE) then
!        -------------------------------------------------------------
         call genRadialGrid(id=i, rmt=rmt, rinsc=rinsc, rend=rend,    &
                            ndivin=ndivin, xstep=hin)
!        -------------------------------------------------------------
      else if (isASAPotential()) then
!        -------------------------------------------------------------
         call genRadialGrid(id=i, rmt=rmt, rinsc=rinsc, rend=rend,    &
                            ndivin=ndivin)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call genRadialGrid(id=i, rmt=rmt, rinsc=rinsc, rend=rend,    &
                            ndivin=ndivin, rfix=rc)
!        -------------------------------------------------------------
      endif
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printRadialGrid(i)
!        -------------------------------------------------------------
      endif
      rmt_grid = getRadialGridRadius(i,MT=.true.)
      if (abs(rmt_grid-rmt) > TEN2m6) then
!        -------------------------------------------------------------
         call setMuffinTinRadius(i,rmt_grid)
!        -------------------------------------------------------------
         if (.not.isFullPotential()) then
!           ==========================================================
!           For Muffin-tin, ASA, or Muffin-tin-ASA calculations, since
!           potential outside rmt is 0, core radius is set to rmt
!           ----------------------------------------------------------
            call setAtomCoreRad(i,rmt_grid)
!           ----------------------------------------------------------
         endif
!        -------------------------------------------------------------
         call setSystemVolumeMT(i,rmt_grid)
!        -------------------------------------------------------------
         isRmtUpdated = .true.
      endif
      radius(i) = getOutscrSphRadius(i)
   enddo

   lmax_max=getMaxLmax()
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
      lmax_step(i)  = getStepFuncLmax(i) 
      lmax_green(i) = min(2*lmax_phi(i),lmax_rho(i))
      lmax_mad(i) = max(lmax_rho(i),lmax_pot(i))
      lmax_max = max( lmax_max, getTruncPotLmax(i), lmax_rho(i), lmax_pot(i), &
                      lmax_phi(i), lmax_kkr(i), lmax_green(i) )
      lmax_kkr_max  = max(lmax_kkr_max,lmax_kkr(i))
      lmax_rho_max  = max(lmax_rho_max,lmax_rho(i))
      lmax_pot_max  = max(lmax_pot_max,lmax_pot(i))
   enddo
   MaxVal_Integer(1) = lmax_rho_max
   MaxVal_Integer(2) = lmax_pot_max
!  -------------------------------------------------------------------
   call GlobalMax(MaxVal_Integer,2)
!  -------------------------------------------------------------------
   lmax_rho_max = MaxVal_Integer(1)
   lmax_pot_max = MaxVal_Integer(2)
   lmax_gaunt = max(lmax_max,2*lmax_pot_max,2*lmax_rho_max)
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_gaunt)
!  -------------------------------------------------------------------
   call initGauntFactors(lmax_gaunt, istop, node_print_level)
!  -------------------------------------------------------------------
   allocate( ngr(LocalNumAtoms), ngt(LocalNumAtoms) )
   do i=1,LocalNumAtoms
      ngr(i) = ngaussr
      ngt(i) = ngaussq
   enddo
!  -------------------------------------------------------------------
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

   t2 = getTime()
!  -------------------------------------------------------------------
   call initMadelung(LocalNumAtoms,GlobalNumAtoms,GlobalIndex,        &
                     lmax_rho_max,lmax_pot_max,bravais,AtomPosition,  &
                     node_print_level)
!  -------------------------------------------------------------------
   if (node_print_level >= 0) then
     write(6,'(/,a,f10.5,/)')'Time:: initMadelung: ',getTime()-t2
   endif
   call initSystemSymmetry( GlobalNumAtoms, LocalNumAtoms,            &
                            lmax_mad, lmax_step, atom_print_level )
!  -------------------------------------------------------------------
   call calSymmetryFlags()
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot, lmax_step,      &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------

   if (isKKRCPA() .or. isKKRCPASRO()) then
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
!       -------------------------------------------------------------
        call printIBZRotationMatrix(Rot3D_Only=.true.)
!       -------------------------------------------------------------
      endif
   endif
   call readPotential()
   Efermi=getPotEf()
   if (node_print_level >= 0) then
     write(6,'(/,a,f12.8)')' Fermi energy read from the potential:',Efermi
   endif
   call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, &
                     lmax_kkr, lmax_phi, lmax_pot, lmax_step, lmax_green, &
                     n_spin_pola, n_spin_cant, 0,   &
                     istop, atom_print_level, derivative=.true.)

   call initMSSolver(LocalNumAtoms, GlobalIndex,                        &
                     lmax_kkr, lmax_phi, lmax_green, AtomPosition,            &
                        n_spin_pola, n_spin_cant, 0,      &
                        istop, atom_print_level, derivative=.true.)

   vc = includeVertexCorrections()
   call initConductivity(LocalNumAtoms, lmax_kkr, lmax_phi, lmax_green, &
           n_spin_pola, n_spin_cant, 0, istop, atom_print_level, vc)

!=====================================================================
!  Determine if there are sine matrix zeros within (0.0, Ef)
!*********************************************************************
   if (getKeyValue(def_id,'Imaginary energy shift',ei) /= 0) then
      ei = 0.001d0
   endif
   allocate(LocalNumSpecies(LocalNumAtoms))
   do id = 1, LocalNumAtoms
      LocalNumSpecies(id) = getLocalNumSpecies(id)
   enddo
!  -------------------------------------------------------------------
   call initSineMatrixZeros(LocalNumAtoms,n_spin_pola,LocalNumSpecies,  &
                            lmax_kkr,lmax_rho,node_print_level)
!  -------------------------------------------------------------------
   do id =  1, LocalNumAtoms
      do is = 1, n_spin_pola
         do ia = 1, LocalnumSpecies(id)
!           ----------------------------------------------------------
            call findSineMatrixZeros(id,ia,is,ZERO,Efermi,EiBound=ei,Delta=pole_step)
!           ----------------------------------------------------------
            if (node_print_level >= 0 .and. getNumSineMatrixZeros(id,ia,is) > 0) then
!              -------------------------------------------------------
               call printSineMatrixZerosInfo(id,ia,is)
!              -------------------------------------------------------
            endif
         enddo
      enddo
   enddo
!  -------------------------------------------------------------------
   call endSineMatrixZeros()
!  -------------------------------------------------------------------
   deallocate(LocalNumSpecies)
!*********************************************************************
!  End of finding sine matrix poles
!=====================================================================

   call calConductivity(Efermi, LocalNumAtoms, n_spin_pola)
  
   do i = 1, 3
     do j = 1, 3
       do k = 1, n_spin_pola
         final_sigma(i, j, k) = (2.30384174*real(sigma(i, j, k)))/100.0
       enddo
     enddo
   enddo  
   
   do k = 1, n_spin_pola
     call MtxInv_LU(final_sigma(:,:,k), 3)
   enddo

   if (node_print_level >= 0) then
     if (n_spin_pola == 1) then
       write(6,'(5x,a)')'                                                        '
       write(6,'(5x,a)')'********************************************************'
       write(6,'(5x,a)')'--------------------------------------------------------'
       write(6,'(5x,a)')'             RESISTIVITY TENSOR (muOhm-cm)          '
       write(6,'(5x,a)')'--------------------------------------------------------'
       write(6,'(3x,f10.5,14x,f10.5,14x,f10.5)') final_sigma(1,1,1),final_sigma(1,2,1),final_sigma(1,3,1)
       write(6,'(3x,f10.5,14x,f10.5,14x,f10.5)') final_sigma(2,1,1),final_sigma(2,2,1),final_sigma(2,3,1)
       write(6,'(3x,f10.5,14x,f10.5,14x,f10.5)') final_sigma(3,1,1),final_sigma(3,2,1),final_sigma(3,3,1)
       write(6,'(5x,a)')'********************************************************'
     else if (n_spin_pola == 2) then
       write(6,'(5x,a)')'                                                        '
       write(6,'(5x,a)')'********************************************************'
       write(6,'(5x,a)')'--------------------------------------------------------'
       write(6,'(/,a,/)')'       RESISTIVITY TENSOR  (SPIN UP) (muOhm-cm)            '
       write(6,'(5x,a)')'--------------------------------------------------------'
       write(6,'(3x,f10.5,15x,f10.5,14x,f10.5)') final_sigma(1,1,1),final_sigma(1,2,1),final_sigma(1,3,1)
       write(6,'(3x,f10.5,15x,f10.5,14x,f10.5)') final_sigma(2,1,1),final_sigma(2,2,1),final_sigma(2,3,1)
       write(6,'(3x,f10.5,15x,f10.5,14x,f10.5)') final_sigma(3,1,1),final_sigma(3,2,1),final_sigma(3,3,1)
       write(6,'(5x,a)')'--------------------------------------------------------'
       write(6,'(/,a,/)')'      RESISTIVITY TENSOR  (SPIN DOWN) (muOhm-cm)          '
       write(6,'(5x,a)')'--------------------------------------------------------'
       write(6,'(3x,f10.5,14x,f10.5,14x,f10.5)') final_sigma(1,1,2),final_sigma(1,2,2),final_sigma(1,3,2)
       write(6,'(3x,f10.5,14x,f10.5,14x,f10.5)') final_sigma(2,1,2),final_sigma(2,2,2),final_sigma(2,3,2)
       write(6,'(3x,f10.5,14x,f10.5,14x,f10.5)') final_sigma(3,1,2),final_sigma(3,2,2),final_sigma(3,3,2)
       write(6,'(5x,a)')'********************************************************'
     endif
   endif

   if (printTildeMatrices() .and. node_print_level >= 0) then
     write(6,'(a)')'                                       '
     write(6,'(a)')'******************************************************************************'
     write(6,'(a)')' Raw conductivity matrices (complex), in Atomic Units, printed for debug purposes'
     write(6,'(a)')' i,j indices represent direction and k represents the spin index                '
     write(6,'(a)')' To see the meaning of sigma-tilde(2,3,4), refer to Phys. Rev. B 31, 3260 (1985)'
     write(6,'(a)')'******************************************************************************'
     write(6,'(a)')'                                       '
     call writeMatrix('sigma', sigma, 3, 3, n_spin_pola)
     call writeMatrix('sigma-tilde1', sigmatilde, 3, 3, n_spin_pola)
     call writeMatrix('sigma-tilde2', sigmatilde2, 3, 3, n_spin_pola)
     call writeMatrix('sigma-tilde3', sigmatilde3, 3, 3, n_spin_pola)
     call writeMatrix('sigma-tilde4', sigmatilde4, 3, 3, n_spin_pola)
   endif

   if (node_print_level >= 0) then
      write(6,'(/,80(''=''))')
      t3 = getTime()-t0
      t2 = t_inp + t_outp
      write(6,'(''Time:: Job total including IO :'',f12.5,''Sec'')') t3
      write(6,'(''       Job total excluding IO :'',f12.5,''Sec'')') t3 - t2
      write(6,'(80(''-''))')
      write(6,*) " Run Ending! "
      call FlushFile(6)
   endif

   nullify(AtomPosition,AtomicNumber)
   deallocate( GlobalIndex, LocalAtomPosi, atom_print_level )
   deallocate( lmax_kkr, lmax_phi, lmax_pot, lmax_rho, lmax_green, lmax_mad)
   deallocate( lmax_step, lmax_tmp , final_sigma)
!
   call endConductivity()
   call endInput()
   call endMadelung()
   call endPotential()
   call endStepFunction()
   call endRadialGrid()
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
   if (isKKRCPA() .or. isKKRCPASRO()) then
      call endMediumHost()
   endif
   call endSystem()
   call endDataServiceCenter()
   call endGroupComm()
   call endMPP()
   call endScfData
   call endKuboData()
   call endBZone()
   if (isKKRCPA() .or. isKKRCPASRO()) then
      call endIBZRotation()
   endif

!  ==================================================================
   call date_and_time(exec_date,exec_time)
   if (node_print_level >= 0) then
     write(6,'(/,12a)')'Execution ends at ',                         &
           exec_time(1:2),':',exec_time(3:4),':',exec_time(5:6),', ', &
           exec_date(5:6),'-',exec_date(7:8),'-',exec_date(1:4)
     write(6,'(80(''-''))')
     write(6,'(a)')'End of a successful run......!'
   endif

   call endOutput()
!
end program kubo
