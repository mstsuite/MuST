program testSSSolver
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8, TEN2m6, TEN2m3, FIVE, CZERO, CONE, TWO, PI, PI2, SQRTm1, HALF, &
                               PI4, THIRD
!
   use MatrixDeterminantModule, only : MtxDet
!
   use TimerModule, only : getTime
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
   use InputModule, only : getKeyValue
!
   use MPPModule, only : packMessage, unpackMessage, sendPackage, recvPackage, AnyPE, MyPE
   use MPPModule, only : setCommunicator, resetCommunicator, GlobalMax, bcastMessage, syncAllPEs
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
   use RadialGridModule, only : getGrid
!
   use SSSolverModule, only : initSSSolver, endSSSolver, solveSingleScattering, &
                              computeDOS, computePDOS, getDOS, getPDOS, getJostMatrix
   use SSSolverModule, only : getSineMatrix, getCosineMatrix, getTMatrix, getSMatrix
!
   use SCPolesModule, only : initSCPoles, endSCPoles, computeSCPoles
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
   implicit none
!
   character (len=2)  :: sform
   character (len=8)  :: exec_date
   character (len=10) :: exec_time
   character (len=10) :: string_rb
   character (len=35), allocatable :: filename(:), filename1(:), filename2(:)
   character (len=50), allocatable :: filename3(:), filename4(:)
!
   logical :: Add_HighL_Contr = .false.
   logical :: NeedPDOS = .false.
   logical :: Search4Poles = .false.
   logical :: ModifyPot = .false.
   logical, allocatable :: checked(:,:,:)
!
   integer (kind=IntKind) :: def_id, info_id, eGID, aGID, NumPEsInGroup, comm
   integer (kind=IntKind) :: i, j, k, id, ig, nstep, is, lmax_max, nw, ldp, ie_loc, ie_glb, jl, l, m, kl, nfit, je, n
   integer (kind=IntKind) :: lmax_step_max, lmax_kkr_max, lmax_phi_max, lmax_rho_max, lmax_pot_max, lmax_green_max
   integer (kind=IntKind) :: ndivin, ndivout, nmult, kmax_kkr, klp, mp, INFO
   integer (kind=IntKind) :: jmax_dos, kmax_dos, jmax_pot
   integer (kind=IntKind) :: LocalNumAtoms, NumAtoms, NumEsOnMyProc, NumEs
   integer (kind=IntKind) :: RelativisticFlag
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind) :: num_poles
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_green(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: lmax_tmp(:)
   integer (kind=IntKind), allocatable :: ngr(:), ngt(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: IPVT(:)
!
   integer (kind=IntKind), parameter :: funit = 11
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: evec(:,:)
   real (kind=RealKind), allocatable :: dos_jl(:)
   real (kind=RealKind), allocatable :: dos_e(:,:,:), dos_mt_e(:,:,:), dos_int(:,:,:), dos_int_fe(:,:)
   real (kind=RealKind), allocatable :: pdos_e(:,:,:,:), pdos_mt_e(:,:,:,:)
   real (kind=RealKind), allocatable :: e_mesh(:), val_tmp(:), dos_corr(:,:,:), dos_corr_rb(:,:,:)
   real (kind=RealKind), allocatable :: kappa_mesh(:)
   real (kind=RealKind), allocatable, target :: deriv_dos_int(:,:), integr_dos_e(:,:), deriv_feidos(:,:)
   real (kind=RealKind), allocatable, target :: integr_dos_corr(:,:), integr_dos_corr_rb(:,:)
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind), pointer :: r_mesh(:), p_val(:)
   real (kind=RealKind), pointer :: potr_0(:)
!
   real (kind=RealKind) :: Efermi, v0, re, vol
   real (kind=RealKind) :: t0, t1, t2, t3
   real (kind=RealKind) :: rmt, rend, rws, rinsc, Rb, si, sr
   real (kind=RealKind) :: dos, sfac, dos_mt, dos_tmp, intdos, intdos0, intdos1, corr, corr_rb, phase
   real (kind=RealKind), allocatable :: diags(:), diags_e(:,:,:,:), phase_e(:,:,:)
!
   real (kind=RealKind), parameter :: xstart = -.1113096740000D+02
!
   complex (kind=CmplxKind) :: energy, a0, delta_p, deto, dets, smdexp, eiphi
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind), pointer :: EPoint(:)
   complex (kind=CmplxKind), pointer :: pdos(:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl_kl(:,:,:)
   complex (kind=CmplxKind), pointer :: pot_l(:)
   complex (kind=CmplxKind), pointer :: jost_mat(:,:), sin_mat(:,:), ci_mat(:,:), xi_mat(:,:), cos_mat(:,:)
   complex (kind=CmplxKind), pointer :: t_mat(:,:), cp_val(:)
   complex (kind=CmplxKind), pointer :: S_matrix(:,:), Sc_matrix(:,:), U_matrix(:,:), xic_mat(:,:)
   complex (kind=CmplxKind), allocatable, target :: space_t1(:), space_t2(:), space_t3(:)
   complex (kind=CmplxKind), allocatable, target :: space_t4(:), space_t5(:), space_t6(:)
   complex (kind=CmplxKind), allocatable :: fint(:), gint(:,:), detS_e(:,:,:), cval_tmp(:)
!
   type (GridStruct), pointer :: Grid
!
!  ===================================================================
!  Initialize modules and data
!  -------------------------------------------------------------------
   call startProcess()
!  -------------------------------------------------------------------
   t0 = getTime()
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
!     lmax_green(i) = min(2*lmax_phi(i),lmax_rho(i))
      lmax_green(i) = lmax_rho(i)
      lmax_green_max = max(lmax_green(i),lmax_green_max)
   enddo
!
   kmax_kkr = (lmax_kkr_max + 1)**2
!
   allocate( IPVT(kmax_kkr) )
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
!
!     ================================================================
!     If necessary, use the following code to set up a model full-potential
!     ================================================================
      if (ModifyPot) then
         do id = 1, LocalNumAtoms
            Grid => getGrid(id)
            do is = 1,n_spin_pola
!              -------------------------------------------------------
!ywg           call setPotentialOutsideMT(id,1,is,nfit=3)
!              -------------------------------------------------------
               jl = 0
               do l = 0, lmax_pot(id)
                  do m = 0, l
                     jl = jl + 1
                     if (.not.isPotComponentZero(id,jl)) then
                        if (l > 0) then
                           pot_l => getPotential(id,1,is,jl)
                           rmt = getInscrSphRadius(id)
                           a0 = pot_l(Grid%jmt)/rmt**l
                           do i = 1, Grid%jmt
                              pot_l(i) = a0*Grid%r_mesh(i)**l
                           enddo
                        endif
                     endif
                  enddo
               enddo
            enddo
         enddo
      endif
!     ================================================================
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
!  -------------------------------------------------------------------
   call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, &
                     lmax_kkr, lmax_phi, lmax_pot, lmax_step, lmax_green,  &
                     n_spin_pola, n_spin_cant, RelativisticFlag, istop, atom_print_level)
!  -------------------------------------------------------------------
   if (n_spin_pola == 1) then
      sfac = 2.0d0
   else
      sfac = 1.0d0
   endif
!
   write(6,'(/,a,f12.5,a)') 'Setup time: ',getTime() - t0,' sec.'
!
   allocate( dos_jl((lmax_green_max+1)*(lmax_green_max+2)/2) )
!
   aGID = getGroupID('Unit Cell')
   eGID = getGroupID('Energy Mesh')
   if (getMyPEinGroup(eGID) == 0) then
      NumEs = getNumEs()
      allocate(filename(LocalNumAtoms), filename1(LocalNumAtoms), filename2(LocalNumAtoms), filename3(LocalNumAtoms), &
               filename4(LocalNumAtoms))
      do id = 1, LocalNumAtoms
         ig = GlobalIndex(id)
         write(filename(id),'(5i6)')100000+lmax_step(id),100000+lmax_kkr(id),100000+lmax_phi(id),100000+ig,100000+NumEs
         filename(id)(1:4)='_stp'; filename(id)(7:10)='_kkr'; filename(id)(13:16)='_phi'; filename(id)(19:21)='_at'
         filename(id)(25:26)='_e'
         filename1(id) = 'fed'//trim(filename(id))
         filename2(id) = 'smd'//trim(filename(id))
         filename3(id) = 'del'//trim(filename(id))//trim(string_rb)
         filename4(id) = 'DOS'//trim(filename(id))//trim(string_rb)
         filename(id) = 'dos'//trim(filename(id))
         open(unit=funit+id,file=filename(id),form='formatted',status='unknown')
         open(unit=100+funit+id,file=filename1(id),form='formatted',status='unknown')
         open(unit=200+funit+id,file=filename2(id),form='formatted',status='unknown')
         open(unit=300+funit+id,file=filename3(id),form='formatted',status='unknown')
         open(unit=400+funit+id,file=filename4(id),form='formatted',status='unknown')
         write(funit+id,'(5x,2a)')  &
            'Energy(Ryd)              DOS               DOS_mt          Integrated_DOS      Integrate_n(E)      ', &
            'Derivative_N(E)'
         write(100+funit+id,'(6x,2a)')  &
            'Energy(Ryd)         Krein_N(E)            FE_N(E)          d[KN(E)]/dE         d[FN(E)]/dE         ', &
            'Anal-FE_n(E)'
         write(200+funit+id,'(5x,a,$)')'Energy(Ryd)        total-phase'
         do l = 0, lmax_kkr(id)
            if (l < 10) then
               do m = -l, l
                  if (m < -9) then
                     write(200+funit+id,'(9x,a,i1,a,i3,a,$)')'(',l,',',m,')  '
                  else if (m < 0 .or. m > 9) then
                     write(200+funit+id,'(9x,a,i1,a,i2,a,$)')'(',l,',',m,')   '
                  else
                     write(200+funit+id,'(9x,a,i1,a,i1,a,$)')'(',l,',',m,')    '
                  endif
               enddo
            else
               do m = -l, l
                  if (m < -9) then
                     write(200+funit+id,'(9x,a,i2,a,i3,a,$)')'(',l,',',m,')  '
                  else if (m < 0 .or. m > 9) then
                     write(200+funit+id,'(9x,a,i2,a,i2,a,$)')'(',l,',',m,')   '
                  else
                     write(200+funit+id,'(9x,a,i2,a,i1,a,$)')'(',l,',',m,')    '
                  endif
               enddo
            endif
         enddo
         write(200+funit+id,'(a)')' '
         write(300+funit+id,'(5x,3a)')  &
            'Energy(Ryd)        Krein_DOS         FE_DOS         Delta_DOS       Delta_DOS_Rb      Green_DOS        ', &
            'Krein_Cell_DOS      Krein_IDOS       ', &
            'FE_IDOS         Delta_IDOS        Delta_IDOS_Rb       Green_IDOS        Krein_Cell_IDOS'
         write(400+funit+id,'(5x,2a)')  &
            'Energy(Ryd)       Krein_DOS+FE_DOS-Delta_DOS    Krein_DOS+FE_DOS-Delta_DOS_Rb     Green_DOS    ', &
            'Krein_IDOS+FE_IDOS-Delta_IDOS     Krein_IDOS+FE_IDOS-Delta_IDOS_Rb    Green_IDOS'
      enddo
      allocate( dos_mt_e(NumAtoms,NumEs,n_spin_pola), dos_e(NumAtoms,NumEs,n_spin_pola), &
                dos_int(NumAtoms,NumEs,n_spin_pola) )
      allocate( dos_corr(NumAtoms,NumEs,n_spin_pola), dos_corr_rb(NumAtoms,NumEs,n_spin_pola) )
      allocate( checked(NumAtoms,NumEs,n_spin_pola) ); checked = .false.
      allocate( phase_e(NumAtoms,NumEs,n_spin_pola), diags_e(kmax_kkr,NumAtoms,NumEs,n_spin_pola) )
      allocate( detS_e(NumAtoms,NumEs,n_spin_pola) )
      if (NeedPDOS) then
         allocate( pdos_mt_e(kmax_kkr,NumAtoms,NumEs,n_spin_pola), pdos_e(kmax_kkr,NumAtoms,NumEs,n_spin_pola) )
      endif
   endif
!
   if (Search4Poles) then
!     ----------------------------------------------------------------
      call initSCPoles(LocalNumAtoms, lmax_kkr, lmax_phi, n_spin_cant)
!     ----------------------------------------------------------------
      call computeSCPoles(ErBottom, Efermi)
!     ----------------------------------------------------------------
   endif
!
   t1 = getTime()
!
   allocate( space_t1(kmax_kkr*kmax_kkr), space_t2(kmax_kkr*kmax_kkr) )
   allocate( space_t3(kmax_kkr*kmax_kkr) )
   allocate( space_t4(kmax_kkr*kmax_kkr) )
   allocate( space_t5(kmax_kkr*kmax_kkr) )
   allocate( space_t6(kmax_kkr*kmax_kkr) )
   allocate( fint(0:lmax_kkr_max), gint(kmax_kkr,kmax_kkr) )
   allocate( diags(1:kmax_kkr) )
!
   do is = 1, n_spin_pola
      do ie_loc = 1, NumEsOnMyProc
         ie_glb = getEnergyIndex(ie_loc)
         energy = EPoint(ie_glb)
         kappa = sqrt(energy)
         if (getMyPEinGroup(eGID) > 0) then
            call packMessage(ie_glb)
         endif
         do id = 1, LocalNumAtoms
            ig = GlobalIndex(id)
            t0 = getTime()
!           ----------------------------------------------------------
            call solveSingleScattering(is, id, energy, CZERO)
!           ----------------------------------------------------------
            if (lmax_kkr(id) == lmax_phi(id)) then
               kmax_kkr = (lmax_kkr(id)+1)**2
!              -------------------------------------------------------
!              jost_mat => getJostMatrix()
               sin_mat => getSineMatrix()
               cos_mat => getCosineMatrix()
!              -------------------------------------------------------
               ci_mat => aliasArray2_c(space_t1,kmax_kkr,kmax_kkr)
               ci_mat = cos_mat
!              -------------------------------------------------------
               call MtxInv_LU(ci_mat,kmax_kkr)
!              -------------------------------------------------------
               xi_mat => aliasArray2_c(space_t2,kmax_kkr,kmax_kkr)
               xi_mat = CZERO
               do kl = 1, kmax_kkr
                  xi_mat(kl,kl) = CONE
               enddo
!              -------------------------------------------------------
               call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, SQRTm1,       &
                          sin_mat, kmax_kkr, ci_mat, kmax_kkr, CONE, xi_mat, kmax_kkr)
!              -------------------------------------------------------
               call MtxDet(kmax_kkr,xi_mat,deto)
!              call MtxDet(kmax_kkr,jost_mat,delta_p)
!              -------------------------------------------------------
               dets = conjg(deto)
!              intdos = aimag(log(deto)-log(dets))/(n_spin_pola*PI)
!              =======================================================
!              Determine the integrated DOS based on Krein's formula
!              Note: Fortran Im[log(r*e^{i*x})] gives -PI <= x <= PI
!              =======================================================
               intdos0 = aimag(log(deto/dets))/(TWO*PI)
!              if (intdos0 < ZERO) then
!                 intdos = sfac*(intdos0+HALF)
!              else
                  intdos = sfac*intdos0
!              endif
!
!              =======================================================
!              Calculate the S-matrix
!              =======================================================
               xic_mat => aliasArray2_c(space_t3,kmax_kkr,kmax_kkr)
               xic_mat = CZERO
               do kl = 1, kmax_kkr
                  xic_mat(kl,kl) = CONE
               enddo
!              -------------------------------------------------------
               call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, -SQRTm1,      &
                          sin_mat, kmax_kkr, ci_mat, kmax_kkr, CONE, xic_mat, kmax_kkr)
!              -------------------------------------------------------
               call MtxInv_LU(xic_mat,kmax_kkr)
!              -------------------------------------------------------
               S_matrix => aliasArray2_c(space_t4,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, CONE,         &
                          xi_mat, kmax_kkr, xic_mat, kmax_kkr, CZERO, S_matrix, kmax_kkr)
!              -------------------------------------------------------
   !!          S_matrix => getSMatrix()
!
!              =======================================================
!              An alternative way to calculate the integrated DOS
!              delta_p = det[S-matrix]
!              =======================================================
               call MtxDet(kmax_kkr,S_matrix,delta_p)
               intdos = sfac*aimag(log(delta_p))/(TWO*PI)
!              write(6,'(a,f12.5,2x,2d16.8)')'energy, detS = ',real(energy,RealKind),delta_p
!              =======================================================
!
!              =======================================================
!              =======================================================
!              =======================================================
!              A correction term needs to be added to the Krein formula's DOS
!              due to the DOS outside the atomic cell.
!              -------------------------------------------------------
!              t_mat => getTMatrix()
!
!              Check different ways of calculting the t-matrix if needed
!              -------------------------------------------------------
!              call checkTmatrix(kmax_kkr,energy,sin_mat,cos_mat,t_mat,S_matrix)
!              -------------------------------------------------------
               t_mat => aliasArray2_c(space_t5,kmax_kkr,kmax_kkr)
               t_mat = -S_matrix
               do kl = 1, kmax_kkr
                  t_mat(kl,kl) = CONE + t_mat(kl,kl)
               enddo
               t_mat = t_mat/(TWO*kappa*SQRTm1)
!
!              Firstly, add the contribution from a bounding sphere to infinity
!              -------------------------------------------------------
!              if (isMuffinTinPotential()) then
!                 rend = getInscrSphRadius(id)
!              else if (isASAPotential() .or. isMuffinTinASAPotential()) then
!                 rend = getWignerSeitzRadius(id)
!              else
!                 rend =  getOutscrSphRadius(id)
!              endif
               if (isASAPotential() .or. isMuffinTinASAPotential()) then
                  rend = getWignerSeitzRadius(id)
               else  ! In both muffin-tin and full-potential cases, the integration is carried out 
                     ! from muffin-tin radius to infinity
                  rend = getInscrSphRadius(id)
               endif
!              -------------------------------------------------------
               call calIntSphHankelSq0(lmax_kkr(id),rend,energy,fint)
!              -------------------------------------------------------
               corr = ZERO
               do kl = kmax_kkr, 1, -1
                  l = lofk(kl)
!!                corr = corr + TWO*real(kappa,kind=RealKind)*aimag(t_mat(kl,kl)*fint(l))
                  corr = corr + real((S_matrix(kl,kl)-CONE)*fint(l),kind=RealKind)
               enddo
!
!              Now subtract the contribution from the region inside the cell
!              but outside the muffin-tin sphere
!              -------------------------------------------------------
               if (isFullPotential()) then
                  rmt =  getInscrSphRadius(id)
                  rend =  getOutscrSphRadius(id)
!                 ----------------------------------------------------
                  call calIntSphHankelSq1(id,lmax_kkr(id),rmt,rend,energy,gint)
!                 ----------------------------------------------------
                  do kl = 1, kmax_kkr
                     do klp = 1, kmax_kkr
!!                      corr = corr - TWO*real(kappa,kind=RealKind)*aimag(t_mat(klp,kl)*gint(klp,kl))
                        if (klp == kl) then
                           corr = corr - real((S_matrix(kl,kl)-CONE)*gint(kl,kl),kind=RealKind)
                        else
                           corr = corr - real(S_matrix(klp,kl)*gint(klp,kl),kind=RealKind)
                        endif
                     enddo
                  enddo
!                 write(6,'(a,2f12.5)')'t*int[hl*hl] for full pot = ',real(energy,kind=RealKind),corr
               else
!                 write(6,'(a,2f12.5)')'t*int[hl*hl] for sphe pot = ',real(energy,kind=RealKind),corr
               endif
!
               corr = sfac*real(kappa,kind=RealKind)*corr/(TWO*PI)
!              write(6,'(a,2f12.5)')'energy, corr = ',real(energy,kind=RealKind),corr
!
!              now check for the correction term up to a large radius
!              -------------------------------------------------------
               call calIntSphHankelSq0(lmax_kkr(id),Rb,energy,fint)
!              -------------------------------------------------------
               corr_rb = ZERO
               do kl = 1, kmax_kkr
                  l = lofk(kl)
!!                corr_rb = corr_rb + aimag(t_mat(kl,kl)*fint(l))
                  corr_rb = corr_rb + real((S_matrix(kl,kl)-CONE)*fint(l),kind=RealKind)
               enddo
!!             corr_rb = corr - sfac*real(energy,kind=RealKind)*corr_rb/PI
               corr_rb = corr - sfac*real(kappa,kind=RealKind)*corr_rb/(TWO*PI)
!              =======================================================
!              =======================================================
!              =======================================================
!
!
!              =======================================================
!              Test the unitarity of the S-matrix
!                  Note: space_t5 was used by t_mat, now it is used by
!                        Sc_matrix, since t_mat is no longer needed
!                        until the next cycle step of the loop
!              =======================================================
               Sc_matrix => aliasArray2_c(space_t5,kmax_kkr,kmax_kkr)
               do kl = 1, kmax_kkr
                  do klp = 1, kmax_kkr
                     Sc_matrix(klp,kl) = conjg(S_matrix(klp,kl))
                  enddo
               enddo
               U_matrix => aliasArray2_c(space_t6,kmax_kkr,kmax_kkr)
!              -------------------------------------------------------
               call zgemm( 'n', 't', kmax_kkr, kmax_kkr, kmax_kkr, CONE,         &
                          S_matrix, kmax_kkr, Sc_matrix, kmax_kkr, CZERO, U_matrix, kmax_kkr)
!              -------------------------------------------------------
!              call writeMatrix('S*S^{*T}',U_matrix,kmax_kkr,kmax_kkr,1.0d-8)
!              -------------------------------------------------------
!
!              =======================================================
!              Triangularize the S_matrix
!              -------------------------------------------------------
               call ZGETRF_nopivot(kmax_kkr, kmax_kkr, S_matrix, kmax_kkr, IPVT, INFO)
!              -------------------------------------------------------
!
!              =======================================================
!              Determine the generialized partial phase shift based on
!              triangularized S_matrix (Due to Balazs Gyorffy)
!              =======================================================
               phase = ZERO
               do kl = 1, kmax_kkr
                  smdexp = log(S_matrix(kl,kl))
                  si = aimag(S_matrix(kl,kl)); sr = real(S_matrix(kl,kl), kind=RealKind)
                  diags(kl) = aimag(smdexp)/TWO
                  if (lofk(kl) == 2 .and. si < ZERO) then
                     diags(kl) = diags(kl) + PI
                  else if (lofk(kl) < 2 .and. (si > ZERO .and. sr < ZERO)) then
                     diags(kl) = diags(kl) - PI
                  endif
                  phase = phase + diags(kl)
               enddo
            else
               intdos = ZERO; phase = ZERO; diags = ZERO
            endif
            write(6,'(/,a,f12.5,a,i3)') 'solveSS time: ',getTime() - t0,' sec., ie = ',ie_loc
            t0 = getTime()
!           ----------------------------------------------------------
            call computeDOS(Add_HighL_Contr)
!           ----------------------------------------------------------
            write(6,'(/,a,f12.5,a)') 'compDOS time: ',getTime() - t0,' sec.'
            dos_r_jl => getDOS()
            Grid => getGrid(id)
            jmax_dos = (lmax_green(id)+1)*(lmax_green(id)+2)/2
            kmax_dos = (lmax_green(id)+1)**2
!           t0 = getTime()
!           ----------------------------------------------------------
            dos = sfac*getVolumeIntegration( id, Grid%jend, Grid%r_mesh,   &
                                             jmax_dos, 2, dos_r_jl, dos_mt)
!           ----------------------------------------------------------
            dos_mt = sfac*dos_mt
            if (isMuffinTinPotential()) then
               dos = dos_mt
            endif
!           write(6,'(/,a,f12.5,a)') 'getVInt time: ',getTime() - t0,' sec.'
            if (getMyPEinGroup(eGID) > 0) then
               call packMessage(ig)
               call packMessage(dos_mt)
               call packMessage(dos)
               call packMessage(intdos)
               call packMessage(corr)
               call packMessage(corr_rb)
               call packMessage(phase)
               call packMessage(delta_p)
               call packMessage(diags,kmax_kkr)
            else
               checked(ig,ie_glb,is)=.true.
               dos_e(ig,ie_glb,is) = dos
               dos_mt_e(ig,ie_glb,is) = dos_mt
               dos_int(ig,ie_glb,is) = intdos
               dos_corr(ig,ie_glb,is) = corr
               dos_corr_rb(ig,ie_glb,is) = corr_rb
               phase_e(ig,ie_glb,is) = phase
               detS_e(ig,ie_glb,is) = delta_p
               diags_e(1:kmax_kkr,ig,ie_glb,is) =  diags(1:kmax_kkr)
            endif
!
            if (NeedPDOS) then
!              -------------------------------------------------------
               call computePDOS()
!              -------------------------------------------------------
               dos_r_jl_kl => getPDOS()
               kl = 0
               do l = 0, lmax_kkr(id)
                  do m = -l, l
                     kl = kl + 1
                     pdos => dos_r_jl_kl(:,:,kl)
!                    -------------------------------------------------
                     dos = sfac*getVolumeIntegration( id, Grid%jend, Grid%r_mesh,   &
                                                      jmax_dos, 2, pdos, dos_mt )
!                    -------------------------------------------------
                     dos_mt = sfac*dos_mt
                     if (getMyPEinGroup(eGID) > 0) then
                        call packMessage(dos_mt)
                        call packMessage(dos)
                     else
                        pdos_e(kl,ig,ie_glb,is)=dos
                        pdos_mt_e(kl,ig,ie_glb,is)=dos_mt
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo
!
   comm = getGroupCommunicator(eGID)
   NumPEsInGroup = getNumPEsInGroup(eGID)
   call setCommunicator(comm,getMyPEinGroup(eGID),NumPEsInGroup)
   if (getMyPEinGroup(eGID) == 0) then
      do i = 2, NumPEsInGroup
         call recvPackage(101,AnyPE)
         do is = 1, n_spin_pola
            do ie_loc = 1, NumEsOnMyProc
               call unpackMessage(ie_glb)
               do id = 1, LocalNumAtoms
                  kmax_kkr = (lmax_kkr(id)+1)**2
                  call unpackMessage(ig)
                  if (checked(ig,ie_glb,is)) then
                     call unpackMessage(dos_mt)
                     call unpackMessage(dos)
                     call unpackMessage(intdos)
                     call unpackMessage(corr)
                     call unpackMessage(corr_rb)
                     call unpackMessage(phase)
                     call unpackMessage(delta_p)
                     call unpackMessage(diags,kmax_kkr)
                     if (NeedPDOS) then
                        do kl = 0, (lmax_kkr(id)+1)**2
                           call unpackMessage(dos_mt)
                           call unpackMessage(dos)
                        enddo
                     endif
                  else
                     call unpackMessage(dos_mt_e(ig,ie_glb,is))
                     call unpackMessage(dos_e(ig,ie_glb,is))
                     call unpackMessage(dos_int(ig,ie_glb,is))
                     call unpackMessage(dos_corr(ig,ie_glb,is))
                     call unpackMessage(dos_corr_rb(ig,ie_glb,is))
                     call unpackMessage(phase_e(ig,ie_glb,is))
                     call unpackMessage(detS_e(ig,ie_glb,is))
                     call unpackMessage(diags_e(1:kmax_kkr,ig,ie_glb,is),kmax_kkr)
                     if (NeedPDOS) then
                        do kl = 0, (lmax_kkr(id)+1)**2
                           call unpackMessage(pdos_mt_e(kl,ig,ie_glb,is))
                           call unpackMessage(pdos_e(kl,ig,ie_glb,is))
                        enddo
                     endif
                     checked(ig,ie_glb,is)=.true.
                  endif
               enddo
            enddo
         enddo
      enddo
   else
      call sendPackage(101,0)
   endif
   call resetCommunicator()
!
   if (getMyPEinGroup(eGID) == 0) then
!
      allocate( e_mesh(NumEs), deriv_dos_int(NumEs,LocalNumAtoms), integr_dos_e(NumEs,LocalNumAtoms) )
      allocate( kappa_mesh(NumEs) )
      allocate( val_tmp(NumEs), deriv_feidos(NumEs,LocalNumAtoms), dos_int_fe(NumEs,LocalNumAtoms) )
      allocate( integr_dos_corr(NumEs,LocalNumAtoms), integr_dos_corr_rb(NumEs,LocalNumAtoms), cval_tmp(NumEs) )
!
      do ie_glb = 1, NumEs
         e_mesh(ie_glb) = real(EPoint(ie_glb),kind=RealKind)
         kappa_mesh(ie_glb) = sqrt(e_mesh(ie_glb))
      enddo
!
      do id = 1, LocalNumAtoms
         if (isMuffinTinPotential() .or. isMuffinTinTestPotential()) then
            rinsc = getInscrSphRadius(id)
            vol = PI4*rinsc**3*THIRD
         else if (isASAPotential() .or. isMuffinTinASAPotential()) then
            rinsc = getWignerSeitzRadius(id)
            vol = PI4*rinsc**3*THIRD
         else
            vol = getVolume(id)
         endif
         do ie_glb = 1, NumEs
            dos_int_fe(ie_glb,id) = sfac*vol*sqrt(abs(e_mesh(ie_glb))**3)/(6.0d0*PI**2)
         enddo
      enddo
!
      if (getMyPEinGroup(aGID) == 0) then
         write(6,'(/,a)')'========================================================================='
         write(6,'(  a)')'Spin  Atom        Energy Value              DOS            Integrated DOS'
         write(6,'(  a)')'-------------------------------------------------------------------------'
      endif
      do is = 1, n_spin_pola
!        =============================================================
!        The scheme here for setting up integrated DOS needs to be further tested.
!        Note: both dos_e and dos_int contain a factor of 2/n_spin_pola
!              dos_int is calculated based on Krein's formula and is
!              relative to the free electron IDOS
!        =============================================================
         do id = 1, LocalNumAtoms
            ig = GlobalIndex(id)
!
!           ==========================================================
!           Add free-electron integrated DOS to the Krein theorem formula
!           ==========================================================
!           do ie_glb = 1, NumEs
!              dos_int(ig,ie_glb,is) = dos_int(ig,ie_glb,is) + dos_int_fe(ie_glb,id)
!           enddo
!
            do ie_glb = 1, NumEs
               val_tmp(ie_glb) = dos_corr(ig,ie_glb,is)
            enddo
            p_val => integr_dos_corr(1:NumEs,id)
!           ----------------------------------------------------------
            call calIntegration(0,NumEs,e_mesh,val_tmp,p_val)
!           ----------------------------------------------------------
            do ie_glb = 1, NumEs
               val_tmp(ie_glb) = dos_corr_rb(ig,ie_glb,is)
            enddo
            p_val => integr_dos_corr_rb(1:NumEs,id)
!           ----------------------------------------------------------
            call calIntegration(0,NumEs,e_mesh,val_tmp,p_val)
!           ----------------------------------------------------------
!
!           ==========================================================
!           Determine the integrated DOS
!           Another approach is to set dos_int = sfac*phase_e/PI
!           ==========================================================
!           intdos0 = ZERO
            intdos0 = dos_int(ig,1,is)
            do ie_glb = 1, NumEs
               dos_int(ig,ie_glb,is) = dos_int(ig,ie_glb,is)-intdos0
            enddo
            do ie_glb = 2, NumEs
               if (dos_int(ig,ie_glb-1,is)+dos_int_fe(ie_glb-1,id) - &
                   (dos_int(ig,ie_glb,is)+dos_int_fe(ie_glb,id)) > 0.00001d0 ) then
                  intdos1 = (dos_e(ig,ie_glb-1,is)-dos_int_fe(ie_glb-1,id)+dos_corr(ig,ie_glb-1,is)+ &
                             dos_e(ig,ie_glb,is)-dos_int_fe(ie_glb,id)+                              &
                             dos_corr(ig,ie_glb,is))*abs(EPoint(ie_glb)-EPoint(ie_glb-1))/TWO
                  n = 0
                  do while (dos_int(ig,ie_glb,is)-dos_int(ig,ie_glb-1,is)+n+0.5 < intdos1)
                     n = n + 2
                  enddo
!                 if (n > 0) then
!                    write(6,'(a,3f12.7,i5)') 'e, dos, n = ',e_mesh(ie_glb),dos_int(ig,ie_glb-1,is), &
!                                             dos_int(ig,ie_glb,is),n
!                 endif
                  do je = ie_glb, NumEs
                     dos_int(ig,je,is) = dos_int(ig,je,is)+n
                  enddo
               endif
            enddo
         enddo
!        =============================================================
!
!        =============================================================
!        Perform integration on dos_e and derivative on dos_int so to
!        compare the results with dos_int and dos_e, respectively.
!        Note: integr_dos_e (integration of dos_e) is IDOS calculated
!              from Green function
!              deriv_dos_int (derivative of dos_int) is DOS calculated
!              from Krein's formula and is relative to the free electron DOS
!        =============================================================
         do id = 1, LocalNumAtoms
            ig = GlobalIndex(id)
            do ie_glb = 1, NumEs
               val_tmp(ie_glb) = dos_e(ig,ie_glb,is)
            enddo
            p_val => integr_dos_e(1:NumEs,id)
!           ----------------------------------------------------------
            call calIntegration(0,NumEs,e_mesh,val_tmp,p_val)
!           ----------------------------------------------------------
!
!           ==========================================================
!           Calculate the derivative of the det[S-matrix] and calculate
!           the derivative of the integrated DOS in Krein's formula
!           ==========================================================
            do ie_glb = 1, NumEs
               val_tmp(ie_glb) = dos_int(ig,ie_glb,is)
!Nov5          val_tmp(ie_glb) = dos_int(ig,ie_glb,is)/kappa_mesh(ie_glb)
            enddo
            p_val => deriv_dos_int(1:NumEs,id)
!           ----------------------------------------------------------
!           call newder(val_tmp,p_val,e_mesh,NumEs)
!Nov2       call derv5(val_tmp,p_val,e_mesh,NumEs)
            call derv5(val_tmp,p_val,kappa_mesh,NumEs)
!Nov5       call newder(val_tmp,p_val,kappa_mesh,NumEs)
            do ie_glb = 1, NumEs
!Nov5          deriv_dos_int(ie_glb,id) = HALF*(e_mesh(ie_glb)*deriv_dos_int(ie_glb,id) + &
!Nov5                                           dos_int(ig,ie_glb,is))/e_mesh(ie_glb)
               deriv_dos_int(ie_glb,id) = HALF*deriv_dos_int(ie_glb,id)/kappa_mesh(ie_glb)
            enddo
!           ----------------------------------------------------------
!!          call derv5(detS_e(1:NumEs,ie_glb,is),cval_tmp,e_mesh,NumEs)
!           ----------------------------------------------------------
!!          do ie_glb = 1, NumEs
!!             deriv_dos_int(ie_glb,id) = sfac*aimag(cval_tmp(ie_glb)/detS_e(ig,ie_glb,is))/PI2 &
!!             deriv_dos_int(ie_glb,id) = deriv_dos_int(ie_glb,id) - dos_corr(ig,ie_glb,is)
!!          enddo
!
!           ==========================================================
!           This is to check the derivative of the free electron integrated DOS
!           ----------------------------------------------------------
!           call newder(dos_int_fe(1:NumEs,id),deriv_feidos(1:NumEs,id),e_mesh,NumEs)
            call derv5(dos_int_fe(1:NumEs,id),deriv_feidos(1:NumEs,id),e_mesh,NumEs)
!           ----------------------------------------------------------
         enddo
!
!        =============================================================
!        Note: in the following print statement, dos_int and deriv_dos_int
!              from the Krein' formula does NOT include the free electron
!              contribution
!        =============================================================
         do ie_glb = 1, NumEs
            energy = EPoint(ie_glb)
            do id = 1, LocalNumAtoms
               ig = GlobalIndex(id)
               write(funit+id,'(2x,d16.8,5(4x,d16.8))')e_mesh(ie_glb),dos_e(ig,ie_glb,is),                         &
                                                       dos_mt_e(ig,ie_glb,is),dos_int(ig,ie_glb,is),   &
                                                       integr_dos_e(ie_glb,id),deriv_dos_int(ie_glb,id)
               write(100+funit+id,'(2x,d16.8,5(4x,d16.8))')e_mesh(ie_glb),                                         &
                              dos_int(ig,ie_glb,is)+dos_int_fe(ie_glb,id), dos_int_fe(ie_glb,id),  &
                              deriv_dos_int(ie_glb,id)+deriv_feidos(ie_glb,id),              &
                              deriv_feidos(ie_glb,id),sfac*getVolume(id)*sqrt(abs(e_mesh(ie_glb)))/(4.0d0*PI**2)
!
               kmax_kkr = (lmax_kkr(id)+1)**2
               write(200+funit+id,'(2x,d16.8,2x,d16.8,$)')real(energy,kind=RealKind),phase_e(ig,ie_glb,is)
               do kl = 1, kmax_kkr-1
                  write(200+funit+id,'(2x,d16.8,$)') diags_e(kl,ig,ie_glb,is)
               enddo
               write(200+funit+id,'(2x,d16.8)') diags_e(kmax_kkr,ig,ie_glb,is)
!
               write(300+funit+id,'(2x,d16.8,5(1x,d16.8),2x,6(2x,d16.8),5x,d16.8)')e_mesh(ie_glb),            &
                     deriv_dos_int(ie_glb,id), sfac*getVolume(id)*sqrt(abs(e_mesh(ie_glb)))/(4.0d0*PI**2),    &
                     dos_corr(ig,ie_glb,is), dos_corr_rb(ig,ie_glb,is), dos_e(ig,ie_glb,is),                  &
                     deriv_dos_int(ie_glb,id)+sfac*getVolume(id)*sqrt(abs(e_mesh(ie_glb)))/(4.0d0*PI**2)-     &
                     dos_corr(ig,ie_glb,is), dos_int(ig,ie_glb,is), dos_int_fe(ie_glb,id), integr_dos_corr(ie_glb,id),&
                     integr_dos_corr_rb(ie_glb,id), integr_dos_e(ie_glb,id),                                  &
                     dos_int(ig,ie_glb,is)+dos_int_fe(ie_glb,id)-integr_dos_corr(ie_glb,id)
!
               write(400+funit+id,'(2x,d16.8,1(8x,d16.8),8x,3(8x,d16.8),12x,2(8x,d16.8))')e_mesh(ie_glb),     &
                     deriv_dos_int(ie_glb,id)+sfac*getVolume(id)*sqrt(abs(e_mesh(ie_glb)))/(4.0d0*PI**2)-     &
                     dos_corr(ig,ie_glb,is), deriv_dos_int(ie_glb,id)+sfac*getVolume(id)*                     &
                     sqrt(abs(e_mesh(ie_glb)))/(4.0d0*PI**2)-dos_corr_rb(ig,ie_glb,is), dos_e(ig,ie_glb,is),  &
                     dos_int(ig,ie_glb,is)+dos_int_fe(ie_glb,id)-integr_dos_corr(ie_glb,id),                  &
                     dos_int(ig,ie_glb,is)+dos_int_fe(ie_glb,id)-integr_dos_corr_rb(ie_glb,id),integr_dos_e(ie_glb,id)
!
               if (getMyPEinGroup(aGID) == 0) then
                  write(6,'(1x,i2,4x,i2,6x,d16.8,6x,d16.8,4x,d16.8)')is,ig,real(energy,kind=RealKind),        &
                        dos_e(ig,ie_glb,is),dos_int(ig,ie_glb,is)
                  if (NeedPDOS) then
                     dos_tmp = ZERO
                     kl = 0
                     do l = 0, lmax_kkr(id)
                        do m = -l, l
                           kl = kl + 1
                           write(6,'(34x,i3,1x,i3,f14.8)')l,m,pdos_e(kl,ig,ie_glb,is)
                           dos_tmp = dos_tmp + pdos_e(kl,ig,ie_glb,is)
                        enddo
                     enddo
                     write(6,'(34x,a,f14.8)')'Total =',dos_tmp
                  endif
               endif
            enddo
         enddo
      enddo
      write(6,'(/,a)')'========================================================================='
!
      do ig = 1, NumAtoms
         close(unit=funit+ig)
         close(unit=100+funit+ig)
         close(unit=200+funit+ig)
         close(unit=300+funit+ig)
      enddo
!
      deallocate( dos_mt_e, dos_e, dos_int, checked, dos_corr, dos_corr_rb, fint, gint )
      deallocate( e_mesh, deriv_dos_int, integr_dos_e, val_tmp, deriv_feidos, dos_int_fe )
      deallocate( phase_e, detS_e, diags_e, integr_dos_corr, integr_dos_corr_rb )
      if (NeedPDOS) then
         deallocate( pdos_mt_e, pdos_e )
      endif
!
      write(6,'(/,a,f12.5,a)') 'Compute time: ',getTime() - t1,' sec.'
   endif
!
   call syncAllPEs()
!
!  ===================================================================
!  clean up the allocated spaces.
!  ===================================================================
   nullify(EPoint, dos_r_jl)
   if (getMyPEinGroup(eGID) == 0) then
      deallocate( filename, filename1, filename2, filename3, filename4 )
   endif
   deallocate( atom_print_level )
   deallocate( AtomPosition,AtomicNumber,lmax_kkr,lmax_phi,lmax_rho,lmax_pot )
   deallocate( lmax_step, lmax_green, GlobalIndex, evec, dos_jl )
   deallocate( IPVT)
   deallocate( space_t1, space_t2, space_t3, space_t4, space_t5, space_t6 )
   deallocate( diags )
!
   if (Search4Poles) then
      call endSCPoles()
   endif
   call endSSSolver()
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
!  ===================================================================
!  finish the process
!  -------------------------------------------------------------------
   call finishProcess()
!  -------------------------------------------------------------------
!
contains
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function computeSingleSiteDOS(is,id,e,m,aux) result(dos)
!  ===================================================================
!  
!  This function returns single site DOS for a given energy on the real
!  energy axis.
!  
!  The function also returns aux, a data set that can be integrated
!  on an energy grid 
!  ===================================================================
   use PotentialModule, only : getVdif
   use StepFunctionModule, only : getVolumeIntegration
   use SSSolverModule, only : solveSingleScattering, computeDOS,  getDOS
   use ScfDataModule, only : isInterstitialElectronPolarized
   use RadialGridModule, only : getGrid
   implicit none
!  
   integer (kind=IntKind), intent(in) :: is, id
!  
   real (kind=RealKind), intent(in) :: e
!  
   integer (kind=IntKind), intent(out) :: m   ! this is the actual size of aux
   integer (kind=IntKind) :: jmax_dos, iend, n, js
!  
   real (kind=RealKind) :: dos(n_spin_cant)
   real (kind=RealKind) :: sfac, vdif, dos_mt
   real (kind=RealKind), pointer :: p_vdif(:)
!  
   complex (kind=CmplxKind), intent(out) :: aux(*)
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
!  
   type (GridStruct), pointer :: Grid
!  
   jmax_dos = (lmax_green(id)+1)*(lmax_green(id)+2)/2
   sfac= TWO/real(n_spin_pola,kind=RealKind)
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
!  -------------------------------------------------------------------
   call computeDOS()
!  -------------------------------------------------------------------
!  
   Grid => getGrid(id)
   m = 0
   do js = 1, n_spin_cant
!     ----------------------------------------------------------------
      dos_r_jl => getDOS(js,id)
!     ----------------------------------------------------------------
      iend = size(dos_r_jl,1); jmax_dos = size(dos_r_jl,2); n = iend*jmax_dos
!     ----------------------------------------------------------------
      dos(js) = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,     &
                                           jmax_dos, 2, dos_r_jl, dos_mt )
      dos_mt = sfac*dos_mt
!     ----------------------------------------------------------------
      call zcopy(n,dos_r_jl,1,aux(m+1),1)
!     ----------------------------------------------------------------
      m = m + n
      aux(m+1) = dos(js)
      aux(m+2) = dos_mt
      aux(m+3) = dos(js)*energy
      m = m + 3
   enddo
!  
   end function computeSingleSiteDOS
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntSphHankelSq0(lmax,rc,energy,fint)
!  ===================================================================
      use KindParamModule, only : IntKind, RealKind, CmplxKind
      use MathParamModule, only : Half, SQRTm1
      use BesselModule, only : SphericalBessel, SphericalNeumann
!
      implicit none
!
      integer (kind=IntKind), intent(in) :: lmax
      integer (kind=IntKind) :: l
!
      real (kind=RealKind), intent(in) :: rc
!
      complex (kind=CmplxKind), intent(in) :: energy
      complex (kind=CmplxKind), intent(out) :: fint(0:lmax)
      complex (kind=CmplxKind) :: bjl(0:lmax+1), bnl(0:lmax+1), x, bhl, bhlp1, bhlm1
      complex (kind=CmplxKind) :: Ul, Ulp1
!
      x = rc*sqrt(energy)+SQRTm1*0.00001d0
!     ----------------------------------------------------------------
      call SphericalBessel(lmax+1,x,bjl)
      call SphericalNeumann(lmax+1,x,bnl)
!     ----------------------------------------------------------------
      bhl = bjl(0)+SQRTm1*bnl(0)
      bhlp1 = bjl(1)+SQRTm1*bnl(1)
      fint(0) = (bhl*bhlp1/x-bhl**2-bhlp1**2)*HALF
      do l = 1, lmax
         bhlm1 = bhl
         bhl = bhlp1
         bhlp1 = bjl(l+1)+SQRTm1*bnl(l+1)
!        fint(l) = ((2*l+1)*bhl*bhlp1/x-bhl**2-bhlp1**2)*HALF
         fint(l) = (bhlm1*bhlp1-bhl**2)*HALF
      enddo
!
!     ================================================================
!     Another way of calculating fint is to use recursive relation.
!     Both ways have shown to give the same results.
!     ================================================================
!!    Ul = -SQRTm1*exp(SQRTm1*TWO*x)/(TWO*x**3)
!!    write(6,'(a,2d16.8,2x,2d16.8)')'fint(0), Ul = ',fint(0), Ul
!!    fint(0) = Ul
!!    do l = 1, lmax
!!       bhlm1 = bjl(l-1)+SQRTm1*bnl(l-1)
!!       bhl = bjl(l)+SQRTm1*bnl(l)
!!       Ulp1 = (bhl**2+bhlm1**2+(2*l+1)*fint(l-1))/(2*l-1.0d0)
!!       write(6,'(a,2d16.8,2x,2d16.8)')'fint(l), Ulp1 = ',fint(l), Ulp1
!!       fint(l) = Ulp1
!!    enddo
!     ================================================================
!
      do l = 0, lmax
         fint(l) = fint(l)*rc**3
      enddo
!
   end subroutine calIntSphHankelSq0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntSphHankelSq1(id,lmax,rm,rc,energy,fint)
!  ===================================================================
      use KindParamModule, only : IntKind, RealKind, CmplxKind
      use MathParamModule, only : Half, SQRTm1
      use BesselModule, only : SphericalBessel, SphericalNeumann
      use StepFunctionModule, only : getNumGaussRs, getGaussR, getGaussRWeight, &
                                     getSigmaLL
!
      implicit none
!
      integer (kind=intKind), intent(in) :: id, lmax
      integer (kind=intKind) :: ng, ig, kl, klp, kmax, l, lp
!
      real (kind=RealKind), intent(in) :: rm, rc
      real (kind=RealKind), pointer :: rg(:), wg(:)
!
      complex (kind=CmplxKind), intent(in) :: energy
      complex (kind=CmplxKind), intent(out) :: fint(:,:)
      complex (kind=CmplxKind) :: bjl(0:lmax+1), bnl(0:lmax+1)
      complex (kind=CmplxKind), pointer :: sig_LL(:,:,:)
      complex (kind=CmplxKind) :: x, bhl, bhlp
      complex (kind=CmplxKind) :: Ul, Ulp1
!
      kmax = (lmax+1)**2
      ng = getNumGaussRs(id)
      rg => getGaussR(id)
      wg => getGaussRWeight(id)
      sig_LL => getSigmaLL(id)
!
      fint = CZERO
      do ig = 1, ng
         if (rg(ig) < rm .or. rg(ig) > rc) then
            write(6,'(a,3f12.5)')'rg,rm,rc = ',rg(ig),rm,rc
            call ErrorHandler('calIntSphHankelSq','rg is out of range')
         endif
         x = rg(ig)*sqrt(energy)+SQRTm1*0.00001d0
!        -------------------------------------------------------------
         call SphericalBessel(lmax+1,x,bjl)
         call SphericalNeumann(lmax+1,x,bnl)
!        -------------------------------------------------------------
         do kl = 1, kmax
            l = lofk(kl)
            bhl = bjl(l)+SQRTm1*bnl(l)
            do klp = 1, kmax
               lp = lofk(klp)
               bhlp = bjl(lp)+SQRTm1*bnl(lp)
!              if (klp == kl) then
!                 fint(klp,kl) = fint(klp,kl)+(sig_LL(klp,kl,ig)-CONE)*bhl*bhlp*wg(ig)*rg(ig)**2
!              else
!                 fint(klp,kl) = fint(klp,kl)+sig_LL(klp,kl,ig)*bhl*bhlp*wg(ig)*rg(ig)**2
!              endif
               fint(klp,kl) = fint(klp,kl)+sig_LL(klp,kl,ig)*bhl*bhlp*wg(ig)*rg(ig)**2
            enddo
         enddo
      enddo
!
   end subroutine calIntSphHankelSq1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkTmatrix(kmax_kkr,energy,sm,cm,t0_mat,S_mat)
!  ===================================================================
!
!  Check different ways of calculting the t-matrix
!
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: kmax_kkr
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind), intent(in) :: sm(:,:), cm(:,:), t0_mat(:,:), S_mat(:,:)
!
   integer (kind=IntKind) :: kl
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind), allocatable :: t1_mat(:,:), t2_mat(:,:), jm(:,:), cmi(:,:), scmi(:,:)
!
   kappa = sqrt(energy)
   allocate(t1_mat(kmax_kkr,kmax_kkr), t2_mat(kmax_kkr,kmax_kkr), jm(kmax_kkr,kmax_kkr), cmi(kmax_kkr,kmax_kkr), &
            scmi(kmax_kkr,kmax_kkr))
!   
!  ===================================================================
   jm = SQRTm1*sm - cm
   call MtxInv_LU(jm,kmax_kkr)
   call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, CONE/kappa,    &
               sm, kmax_kkr, jm, kmax_kkr, CZERO, t1_mat, kmax_kkr)
!  ===================================================================
!
!  ===================================================================
   cmi = cm
   call MtxInv_LU(cmi,kmax_kkr)
   call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, CONE,          &
               sm, kmax_kkr, cmi, kmax_kkr, CZERO, jm, kmax_kkr)
   call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, CONE,          &
               sm, kmax_kkr, cmi, kmax_kkr, CZERO, scmi, kmax_kkr)
   jm = SQRTm1*scmi
   do kl = 1, kmax_kkr
      jm(kl,kl) = jm(kl,kl) - CONE
   enddo
   call MtxInv_LU(jm,kmax_kkr)
   call zgemm( 'n', 'n', kmax_kkr, kmax_kkr, kmax_kkr, CONE/kappa,    &
               scmi, kmax_kkr, jm, kmax_kkr, CZERO, t2_mat, kmax_kkr)
!  ===================================================================
! 
   do kl = 1, kmax_kkr
      write(6,'(a,i3,2x,2d16.8,2x,2d16.8,2x,2d16.8)')'kl,t0-mat,t1-mat,t2-mat:',kl,t0_mat(kl,kl),t1_mat(kl,kl), &
            t2_mat(kl,kl)
      write(6,'(a,i3,2x,2d16.8,2x,2d16.8)')'kl, S-mat, 1-2*k*i*t1-mat:',kl,S_mat(kl,kl), &
            CONE-TWO*kappa*SQRTm1*t1_mat(kl,kl)
!     write(6,'(a,2d16.8,2x,a,2d16.8)')'S-mat:',S_mat(kl,kl),', 1-2*k*i*t2-mat:',CONE-TWO*kappa*SQRTm1*t2_mat(kl,kl)
!     write(6,'(a,2d16.8,2x,a,2d16.8)')'(1-S-mat)/(2*i*kappa):',(CONE-S_mat(kl,kl))/(TWO*kappa*SQRTm1),', t-mat:', &
!                                      t1_mat(kl,kl)
!     write(6,'(a,2d16.8,2x,a,2d16.8)')'(1-S-mat)/(2*i*kappa):',(CONE-S_mat(kl,kl))/(TWO*kappa*SQRTm1),', t-mat:', &
!                                      t2_mat(kl,kl)
   enddo

   deallocate(t1_mat, t2_mat, jm, cmi, scmi)
!
   end subroutine checkTmatrix
!  ===================================================================
end program testSSSolver
