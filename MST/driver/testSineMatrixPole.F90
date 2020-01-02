program testSineMatrixPole
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m8
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
   integer (kind=IntKind) :: i, j, k, id, ig, nstep, is, lmax_max, nw, ldp
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
   t1 = getTime()
!   call initFindPolesModule( LocalNumAtoms, MyPE, n_spin_pola, n_spin_cant, &
!                             lmax_kkr, lmax_phi, ErBottom, Efermi, pole_step)
!  -------------------------------------------------------------------
   call calScatteringPoles(LocalNumAtoms,lmax_kkr,lmax_phi,           &
                           n_spin_pola,n_spin_cant,ErBottom,Efermi,ZERO,nstep)
!   call findPoles()
!  -------------------------------------------------------------------
!   call endFindPolesModule()
   write(6,'(/,a,f12.5,a)') 'Zooming technique time: ',getTime() - t1,' sec.'
!
   nw = int((Efermi-ErBottom)/(4*.025d0))
   ldp = nw*(lmax_kkr_max+1)**2
   allocate( NumPoles(LocalNumAtoms,n_spin_pola) )
   allocate( Poles(ldp,LocalNumAtoms,n_spin_pola) )
   t1 = getTime()
!  -------------------------------------------------------------------
   call calQuadraticPoles(LocalNumAtoms,lmax_kkr,lmax_phi,            &
                          n_spin_pola,ErBottom,Efermi,ldp,NumPoles,Poles)
!  -------------------------------------------------------------------
   write(6,'(/,a,f12.5,a,/)') 'Quadratic technique time: ',getTime() - t1,' sec.'
   do is = 1,n_spin_pola
      do id = 1,LocalNumAtoms
         write(6,'(a,i3,a,i3)')'Spin index: ',is,',  Atom index: ',id
         write(6,'(a,f10.5,a,f10.5,a,i5)')'The Number of poles found within (',ErBottom,', ', &
                                           Efermi,'): ',NumPoles(id,is)
         do i = 1, NumPoles(id,is)
            write(6,'(f20.12)')Poles(i,id,is)
         enddo
      enddo
   enddo
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
                                n_spin_pola,eb,et,ldp,NumPoles,Poles)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE, Ten2m8, CZERO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getSineMatrix,                    &
                                     solveSingleScattering
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix, &
                                     solveQuadraticEquation, getEigenValue
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: lphi(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: n_spin_pola,ldp
   integer (kind=IntKind), intent(out) :: NumPoles(LocalNumAtoms,n_spin_pola)
!
   integer (kind=IntKind) :: id, ie, is, iw, n, kmax_kkr, NumWindows, info, nv
!
   real (kind=RealKind), parameter :: Delta = 0.025d0
   real (kind=RealKind), parameter :: WindowWidth = 4*Delta
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind), intent(out) :: Poles(ldp,LocalNumAtoms,n_spin_pola)
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0
!
   complex (kind=CmplxKind) :: e
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), allocatable :: s0(:,:), s1(:,:), s2(:,:)
   complex (kind=CmplxKind), pointer :: pv(:)
!
   if (eb > et) then
      call ErrorHandler('calQuadraticPoles','eb > et',eb,et)
   endif
!
   de = Delta; de2 = de*TWO; dede2 = de*de*TWO
   NumWindows = int((et-eb)/WindowWidth)
!
   kmax_kkr = (lkkr(1)+1)**2
   allocate( s0(1:kmax_kkr,1:kmax_kkr) )
   allocate( s1(1:kmax_kkr,1:kmax_kkr) )
   allocate( s2(1:kmax_kkr,1:kmax_kkr) )
   call initQuadraticMatrix(kmax_kkr)
!
   do is = 1,n_spin_pola
      do id = 1,LocalNumAtoms
         if ((lkkr(id)+1)**2 /= kmax_kkr) then
            deallocate(s0, s1, s2)
            kmax_kkr = (lkkr(id)+1)**2
            allocate( s0(1:kmax_kkr,1:kmax_kkr) )
            allocate( s1(1:kmax_kkr,1:kmax_kkr) )
            allocate( s2(1:kmax_kkr,1:kmax_kkr) )
            call endQuadraticMatrix()
            call initQuadraticMatrix(kmax_kkr)
         endif
!
         n = 0
         do iw = 1, NumWindows
            w0 = eb + (iw-1)*WindowWidth
            e0 = w0 + HALF*WindowWidth
            if (abs(e0) < Ten2m6 .or. abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6) then
               e0 = e0 + HALF*de
            endif
!
            e = cmplx(e0,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
            sin_mat => getSineMatrix(is,id)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0,1)
!           ----------------------------------------------------------
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
            e = cmplx(e0-de,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
            sin_mat => getSineMatrix(is,id)
!           ----------------------------------------------------------
!
            s1 = (s2 - sin_mat)/de2
            s2 = (s2 + sin_mat - TWO*s0)/dede2
!
!           ----------------------------------------------------------
            call solveQuadraticEquation(s0,s1,s2,info)
!           ----------------------------------------------------------
!
            if (info /= 0) then
               stop 'Error in s0, s1, s2'
            endif
!
!           ----------------------------------------------------------
            pv => getEigenValue(nv)
!           ----------------------------------------------------------
            do ie = 1, kmax_kkr*2
               if (abs(aimag(pv(ie))) < Ten2m8) then
                  pe = real(pv(ie),kind=RealKind) + e0
                  if (pe >= w0 .and. pe <= w0+WindowWidth) then
                     n = n + 1
                     Poles(n,id,is) = pe
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
   deallocate( s0, s1, s2 )
!
   end subroutine calQuadraticPoles
