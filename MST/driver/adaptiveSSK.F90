program adaptiveSSK
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8, TEN2m6, TEN2m3, FIVE, CZERO, CONE, TWO, PI, PI2, SQRTm1, HALF, &
                               PI4, THIRD, FOUR, ONE
!
   use MatrixDeterminantModule, only : MtxDet
!
   use TimerModule, only : getTime
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
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
   use ScfDataModule, only : istop, ErTop, ErBottom
   use ScfDataModule, only : NumSS_IntEs, getAdaptiveIntegrationMethod
   use ScfDataModule, only : getSingleSiteSolverType
!  use ScfDataModule, only : ss_integral_method, ss_solution_method, ss_irrsol
!
   use PotentialTypeModule, only : isASAPotential, isMuffinTinPotential,     &
                                   isMuffinTinASAPotential, isFullPotential
!
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber,   &
                            getBravaisLattice
!
   use InputModule, only : getKeyValue
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
   use ProcMappingModule, only : getNumEsOnMyProc, getEnergyIndex,    &
                                 getNumRedundantEsOnMyProc,           &
                                 getProcWithEnergyIndex
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use RadialGridModule, only : getGrid
!
   use SSSolverModule, only : initSSSolver, endSSSolver
!
   use WriteMatrixModule,  only : writeMatrix
   use MatrixInverseModule, only : MtxInv_LU
!
   use IntegerFactorsModule, only : lofk, mofk, m1m
!
   use IntegrationModule, only : calIntegration
!
   use DerivativeModule, only : derv5
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,   &
                                    calSymmetryFlags
!
   use MadelungModule, only : initMadelung, endMadelung
!
   use AdaptIntegrationModule, only : initAdaptIntegration, endAdaptIntegration
   use AdaptIntegrationModule, only : setupAdaptMesh, getAdaptIntegration
   use AdaptIntegrationModule, only : getAuxDataAdaptIntegration
   use AdaptIntegrationModule, only : getUniformIntegration, getUniMeshIntegration
!
   implicit none
!
   character (len=2)  :: sform
   character (len=8)  :: exec_date
   character (len=10) :: exec_time
   character (len=10) :: string_rb
   character (len=35), allocatable :: filename(:,:)
   character (len=50), allocatable :: filename1(:,:), filename2(:,:)
!
   logical :: NeedPDOS = .false.
   logical, allocatable :: checked(:,:,:)
!
   integer (kind=IntKind) :: def_id, info_id, eGID, aGID, NumPEsInGroup, comm
   integer (kind=IntKind) :: i, j, k, id, ig, is, lmax_max, l, m, jl, kl, n, n1, n2
   integer (kind=IntKind) :: lmax_step_max, lmax_kkr_max, lmax_phi_max, lmax_rho_max, lmax_pot_max, lmax_green_max
   integer (kind=IntKind) :: ndivin, ndivout, nmult, kmax_kkr, klp, mp, nm
   integer (kind=IntKind) :: jmax_dos, kmax_dos, jmax_pot, lmax_gaunt
   integer (kind=IntKind) :: LocalNumAtoms, NumAtoms, NumEsOnMyProc, NumEs
   integer (kind=IntKind) :: RelativisticFlag
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind) :: info(4)
   integer (kind=IntKind), pointer :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_green(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: lmax_tmp(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
!
   integer (kind=IntKind), parameter :: funit = 11
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), pointer :: AtomPosition(:,:)
!
   real (kind=RealKind) :: re, vol, de
   real (kind=RealKind) :: t0, t1, t2, t3
   real (kind=RealKind) :: rmt, rend, rws, rinsc, Rb, ps, IDOS_cell, ssdos_int
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!
   complex (kind=CmplxKind), pointer :: p_aux(:)
!
   type (GridStruct), pointer :: Grid
!
   interface
      function computeSingleSiteDOS(info,e,aux,rfac,redundant) result(dos)
         use KindParamModule, only : intKind, realKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: e
         real (kind=RealKind), intent(in), optional :: rfac
         complex (kind=CmplxKind), intent(out), target :: aux(:)
         real (kind=RealKind) :: dos
         logical, intent(in), optional :: redundant
      end function computeSingleSiteDOS
   end interface
!
   interface
      function computeSingleSitePS(info,e) result(tps)
         use KindParamModule, only : intKind, realKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: e
         real (kind=RealKind) :: tps
      end function computeSingleSitePS
   end interface
!
!  ===================================================================
!  Initialize modules and data
!  -------------------------------------------------------------------
   call startProcess(1)
!  -------------------------------------------------------------------
   t0 = getTime()
   NumAtoms = getNumAtoms()
   LocalNumAtoms=getLocalNumAtoms()
   bravais = getBravaisLattice()
   AtomPosition => getAtomPosition()
   AtomicNumber => getAtomicNumber()
!
   node_print_level = getStandardOutputLevel()
!
   allocate(GlobalIndex(LocalNumAtoms), atom_print_level(1:LocalNumAtoms))
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
   allocate(lmax_step(LocalNumAtoms), lmax_tmp(LocalNumAtoms))
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
      lmax_green(i) = lmax_rho(i)
      lmax_green_max = max(lmax_green(i),lmax_green_max)
   enddo
   lmax_gaunt = max(lmax_max,2*lmax_pot_max,2*lmax_rho_max)
!
   kmax_kkr = (lmax_kkr_max + 1)**2
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_gaunt)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(lmax_gaunt,istop,0)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize radial grid
!  -------------------------------------------------------------------
   call initRadialGrid(LocalNumAtoms, istop, node_print_level)
!  -------------------------------------------------------------------
!
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i)
      GlobalIndex(i)=ig
!     ----------------------------------------------------------------
      call getGridData(i,ndivin,ndivout,nmult)
!     ----------------------------------------------------------------
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printPolyhedron(i)
!        -------------------------------------------------------------
      endif
      rend =  getOutscrSphRadius(i)
      if (isMuffinTinPotential()) then
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
         call genRadialGrid( i, rmt, rinsc, rws, rend, ndivin, ndivout, nmult)
!        -------------------------------------------------------------
      endif
      if (atom_print_level(i) >= 0) then
!        -------------------------------------------------------------
         call printRadialGrid(i)
!        -------------------------------------------------------------
      endif
   enddo
   if (MyPE == 0) then
      if (getKeyValue(1,'Large sphere radius (a.u.)',Rb) > 0) then
         Rb = 500.0d0
         call WarningHandler('testSSSolver','No input for Rb. It is set to default',Rb)
      endif
      write(6,'(a,f12.5)')'You entered Rb = ',Rb
      write(string_rb,'(f10.1)')10000000.0+Rb
      string_rb(1:3)='_Rb'
   endif
   call bcastMessage(Rb,0)
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
   do i=1,LocalNumAtoms
!     lmax_tmp(i) = max(lmax_rho(i), lmax_pot(i), lmax_step(i)+lmax_pot(i))
      lmax_tmp(i) = max(lmax_rho(i), lmax_pot(i))
   enddo
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
!  ===================================================================
!  read potential data
!  -------------------------------------------------------------------
   call readPotential()
!  -------------------------------------------------------------------
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
      enddo
   enddo
!
   NumEs = NumSS_IntEs
   NumEsOnMyProc = getNumEsOnMyProc()
!
!  ===================================================================
!  The following parameters are enforced...
!  ===================================================================
   RelativisticFlag = 0
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
!  comm = getGroupCommunicator(eGID)
!  NumPEsInGroup = getNumPEsInGroup(eGID)
!  call setCommunicator(comm,getMyPEinGroup(eGID),NumPEsInGroup)
!
   t1 = getTime()
!
   if (NumSS_IntEs < 1) then
      call ErrorHandler('main','NumSS_IntEs < 1',NumSS_IntEs)
   else if (ErTop <= ErBottom) then
      call ErrorHandler('main','ErTop <= ErBottom',ErTop,ErBottom)
   else if (ErTop <= ZERO) then
      call ErrorHandler('main','ErTop <= 0.0',ErTop)
   else if (ErBottom <= ZERO) then
      call WarningHandler('main','ErBottom <= 0.0, and is reset to 0.000001',ErBottom)
      ErBottom = 0.000001d0
   endif
!
   n = 0
   do id = 1,LocalNumAtoms
      Grid => getGrid(id)
      jmax_dos = (lmax_green(id)+1)*(lmax_green(id)+2)/2
      n = max(n,Grid%jend_plus_n*jmax_dos)
   enddo
!  -------------------------------------------------------------------
   call initAdaptIntegration(n+4,eGID)
!  -------------------------------------------------------------------
!
   do is = 1, n_spin_pola
      do id = 1, LocalNumAtoms
         if ( node_print_level >= 0) then
            write(6,'(/,2(a,i2))')'is = ',is,', id = ',id
            write(6,'(a,2(f11.8,a))')'In the real energy interval:  [',ErBottom,', ',ErTop,']'
            write(6,'(a)')   '=========================================================================================='
            write(6,'(4x,a)')'Energy    Single_Site_DOS_ws   Single_Site_DOS_mt    DOS_Outside     Total_Phase_Shift'
            write(6,'(a)')   '------------------------------------------------------------------------------------------'
         endif
         info(1) = is; info(2) = id; info(3) = node_print_level; info(4) = eGID
         if (getAdaptiveIntegrationMethod() == 0) then
            ssdos_int = getUniMeshIntegration(NumSS_IntEs,ErBottom,ErTop,info,computeSingleSiteDOS,nm)
         else if (getAdaptiveIntegrationMethod() == 1) then
            call setupAdaptMesh(ErBottom,ErTop,NumSS_IntEs,info,computeSingleSitePS)
            ssdos_int = getAdaptIntegration(info,computeSingleSiteDOS,nm)
         else
            ssdos_int = getUniformIntegration(NumSS_IntEs,ErBottom,ErTop,info,computeSingleSiteDOS,nm)
         endif
         if ( node_print_level >= 0) then
            write(6,'(a)')   '=========================================================================================='
            write(6,'(a,i4)')'Number of mesh points for the integration: ',nm
         endif
         p_aux => getAuxDataAdaptIntegration(); n = size(p_aux); n = n - 4
         ps = computeSingleSitePS(info,ErTop)
         if ( node_print_level >= 0) then
            write(6,'(a,f12.8,a,d15.8)')'At energy = ',ErTop,     &
                                        ': Single site IDOS given by Green function  =',ssdos_int
            write(6,'(26x,a,d15.8)')      'Integrated DOS outside the atomic cell    =',real(p_aux(n+4),kind=RealKind)
!           ==========================================================
            write(6,'(26x,a,d15.8)')      'Single site sum. of partial phase shifts  =',ps
            IDOS_cell = (2/n_spin_pola)*(ps+getVolume(id)*sqrt(ErTop**3)/(6.0d0*PI))/PI - real(p_aux(n+4),kind=RealKind)
            write(6,'(26x,a,d15.8)')      'Single site {phase shift sum-OutsideIDOS} =',IDOS_cell
!           ==========================================================
         endif
      enddo
   enddo
!
!  ===================================================================
!  clean up the allocated spaces.
!  ===================================================================
   nullify(AtomPosition, AtomicNumber, p_aux, Grid)
   deallocate( atom_print_level )
   deallocate( lmax_kkr,lmax_phi,lmax_rho,lmax_pot )
   deallocate( lmax_step, lmax_green, GlobalIndex )
!
   call endAdaptIntegration()
   call endSSSolver()
   call endPotential()
   call endStepFunction()
   call endRadialGrid()
   call endSystemSymmetry()
   call endMadelung()
   call endGauntFactors()
   call endSphericalHarmonics()
!
!  ===================================================================
!  finish the process
!  -------------------------------------------------------------------
   call finishProcess()
!  -------------------------------------------------------------------
end program adaptiveSSK
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function computeSingleSiteDOS(info,e,aux,rfac,redundant) result(dos)
!  ===================================================================
!  
!  This function returns single site DOS for a given energy on the real
!  energy axis.
!  
!  The function also returns aux, a data set that can be integrated
!  on an energy grid 
!  ===================================================================
   use KindParamModule, only : intKind, realKind, CmplxKind
   use MathParamModule, only : ZERO, TWO, CZERO
   use PublicTypeDefinitionsModule, only : GridStruct
   use PotentialModule, only : getVdif
   use StepFunctionModule, only : getVolumeIntegration
   use SSSolverModule, only : solveSingleScattering, computeDOS, getDOS, getOutsideDOS
   use ScfDataModule, only : isInterstitialElectronPolarized
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use RadialGridModule, only : getGrid
   use GroupCommModule, only : getMyPEinGroup, getNumPEsInGroup, GlobalSumInGroup
   implicit none
!  
   integer (kind=IntKind), intent(in) :: info(*)
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind), intent(in), optional :: rfac
   complex (kind=CmplxKind), intent(out), target :: aux(:)
!
   logical, intent(in), optional :: redundant
   logical :: red
!  
   integer (kind=IntKind) :: is, id
   integer (kind=IntKind) :: jmax_dos, iend, n, jl, ir
   integer (kind=IntKind) :: node_print_level, eGID, MyPEinEGroup, NumPEsinEGroup
!  
   real (kind=RealKind) :: dos, dos_out
   real (kind=RealKind) :: sfac, vdif, dos_mt, tps
   real (kind=RealKind), pointer :: p_vdif(:)
   real (kind=RealKind) :: msgbuf(5,100) ! Hard coded here!!!
!  
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
!  
   type (GridStruct), pointer :: Grid
!
   interface
      function computeSingleSitePS(info,e) result(tps)
         use KindParamModule, only : intKind, realKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: info(*)
         real (kind=RealKind), intent(in) :: e
         real (kind=RealKind) :: tps
      end function computeSingleSitePS
   end interface
!
   if (present(redundant)) then
      red = redundant
   else
      red = .false.
   endif
!  
   is = info(1)
   id = info(2)
   node_print_level = info(3)
   eGID = info(4)
   MyPEinEGroup = getMyPEinGroup(eGID)
   NumPEsInEGroup = getNumPEsInGroup(eGID)
!
   sfac= TWO/real(n_spin_pola,kind=RealKind)
   Grid => getGrid(id)
!  
   if (isInterstitialElectronPolarized()) then
      p_vdif => getVdif()
      vdif = p_vdif(1)
   else
      vdif = ZERO
   endif
!  
   energy = e-(is-1)*vdif
!  -------------------------------------------------------------------
   call solveSingleScattering(is,id,energy,CZERO)
   call computeDOS()
   dos_r_jl => getDOS()
!  -------------------------------------------------------------------
   iend = size(dos_r_jl,1); jmax_dos = size(dos_r_jl,2)
!  -------------------------------------------------------------------
   dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,     &
                                    jmax_dos, 2, dos_r_jl, dos_mt )
!  -------------------------------------------------------------------
   dos_mt = sfac*dos_mt
!  -------------------------------------------------------------------
   dos_out = sfac*getOutsideDOS()
   tps = computeSingleSitePS(info,e)
!  -------------------------------------------------------------------
!
   msgbuf = ZERO
   msgbuf(1,MyPEinEGroup+1) = real(energy)
   msgbuf(2,MyPEinEGroup+1) = dos
   msgbuf(3,MyPEinEGroup+1) = dos_mt
   msgbuf(4,MyPEinEGroup+1) = dos_out
   msgbuf(5,MyPEinEGroup+1) = tps
!  -------------------------------------------------------------------
   call GlobalSumInGroup(eGID,msgbuf,5,NumPEsInEGroup)
!  -------------------------------------------------------------------
   if ( node_print_level >= 0 .and. .not.red) then
      do n = 1, NumPEsInEGroup
         write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')msgbuf(1:5,n)
      enddo
   endif
!
   do jl = 1, jmax_dos
      n = (jl-1)*iend
      do ir = 1, iend
         aux(n+ir) = sfac*dos_r_jl(ir,jl)
      enddo
   enddo
!  n = iend*jmax_dos
   n = size(aux) - 4
   aux(n+1) = dos
   aux(n+2) = dos_mt
   aux(n+3) = dos*energy
   aux(n+4) = dos_out
!  
   end function computeSingleSiteDOS
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function computeSingleSitePS(info,e) result(tps)
!  ===================================================================
   use KindParamModule, only : intKind, realKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO
   use ScfDataModule, only : isInterstitialElectronPolarized
   use PotentialModule, only : getVdif
   use SSSolverModule, only : solveSingleScattering, computePhaseShift, getPhaseShift
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(*)
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind) :: vdif
   real (kind=RealKind), pointer :: p_vdif(:)
!
   integer (kind=IntKind) :: kl, kmax_phi, is, id
!
   complex (kind=CmplxKind) :: energy
!
   real (kind=RealKind) :: tps
   real (kind=RealKind), pointer :: ps(:)
!
   is = info(1); id = info(2)
!
   if (isInterstitialElectronPolarized()) then
      p_vdif => getVdif()
      vdif = p_vdif(1)
   else
      vdif = ZERO
   endif
!  
   energy = e-(is-1)*vdif
!  -------------------------------------------------------------------
   call solveSingleScattering(is,id,energy,CZERO)
   call computePhaseShift()
   ps => getPhaseShift()
!  -------------------------------------------------------------------
   kmax_phi = size(ps)
   tps = ZERO
   do kl = 1, kmax_phi
      tps = tps + ps(kl)
   enddo
!
   end function computeSingleSitePS
!  ===================================================================
