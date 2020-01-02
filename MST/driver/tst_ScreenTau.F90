program test_ScreenTau
   use KindParamModule, only : CmplxKind, RealKind, IntKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use MPPModule, only : initMPP, endMPP
   use MPPModule, only : getMyPE, getNumPEs
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
   use PotentialModule, only : readPotential
!
   use ScfDataModule, only : initScfData, nscf, printScfData
   use ScfDataModule, only : pot_type, n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop
   use ScfDataModule, only : ContourType, eGridType
   use SCFDataModule, only : Kdiv,NumKMeshs,kGenScheme,Symmetrize
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
   use AtomModule, only : getScreenPotential, getScreenRcut
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
!
   use SingleScattererModule, only : initSingleScatterer, getTMatrixLocalSpinBlock
   use SingleScattererModule, only : getSinMatrix, getCosMatrix
!
   use SpecKIntegrationModule, only : getBZirr
!
   use TauScreenKKRModule, only : initTauScreenKKRModule, endTauScreenKKRModule
   use TauScreenKKRModule, only : calScreenTau0J, printScreenTau0J, getDeltaM
   use TauScreenKKRModule, only : SpinSpace
!
   use BZoneModule, only : initBZone, getNumKs, getKPoint, getWeight
!
   use MathParamModule, only : czero, cone, sqrtm1
   use WriteMatrixModule, only : writeMatrix
   use MatrixInverseModule, only : MtxInv_GE
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType
!
   implicit none
!
   logical :: isRelativistic
!
   character (len=80) :: info_table(1), info_path(1)
!
   integer (kind=IntKind) :: MyPE, NumPEs
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: i, ig, k, l_max
   integer (kind=IntKind) :: is, ns, rel, ptype, sdim
   integer (kind=IntKind) :: lmax_max,kmax, NumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: ndivin, ndivout, nmult
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind), allocatable :: AtomicNumber(:), AtomType(:)
   integer (kind=IntKind), allocatable :: lmax(:),lkkr(:),lphi(:)
   integer (kind=IntKind), allocatable :: lpot(:),lrho(:)
   integer (kind=IntKind), allocatable :: fit_pot(:),fit_rho(:)
!
   real (kind=RealKind) :: bravais(3,3), kvec(3)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: rmt(:)
   real (kind=RealKind), allocatable :: v_swiss(:)
   real (kind=RealKind), allocatable :: rcut_swiss(:)
   real (kind=RealKind) :: NumKpoints, weight
!
   complex (kind=CmplxKind), pointer :: tau_0j_swiss(:,:,:)
   complex (kind=CmplxKind) :: vshift
   complex (kind=CmplxKind) :: epsilon
   complex (kind=CmplxKind), pointer :: DeltaM(:,:),Tmatrix(:,:)
   complex (kind=CmplxKind), pointer :: TauK(:,:)
!
!  -------------------------------------------------------------------
   call initMPP()
!  -------------------------------------------------------------------
   MyPE = getMyPE()
   NumPEs = getNumPEs()
!
   if (MyPE == 0) then
      write(6,*)'  ******************************************************'
      write(6,*)'  *'
      write(6,*)'  *    Test Screen Structure Constant                  *'
      write(6,*)'  *'
      write(6,*)'  ******************************************************'
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
   if (print_level >= 0) then
!     ----------------------------------------------------------------
      call printSystem()
      call printAtom()
!     ----------------------------------------------------------------
   endif
!  ===================================================================
!  setup factors
!  -------------------------------------------------------------------
   lmax_max=getMaxLmax()
!  -------------------------------------------------------------------
!
   bravais(1:3,1:3)=getBravaisLattice()
   LocalNumAtoms=getLocalNumAtoms(MyPE)
!  ===================================================================
!  initialize radial grid
!  -------------------------------------------------------------------
   call initRadialGrid(LocalNumAtoms, istop, print_level)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call initPolyhedra(LocalNumAtoms,bravais,istop,print_level)
!  -------------------------------------------------------------------
   NumAtoms = getNumAtoms()
!
   allocate(AtomPosition(3,NumAtoms), AtomicNumber(NumAtoms), lmax(NumAtoms))
   allocate(rmt(NumAtoms),rcut_swiss(NumAtoms),v_swiss(NumAtoms))
   allocate(AtomType(NumAtoms))
!
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
      AtomicNumber(i)=getAtomicNumber(i)
   enddo
!
   do i=1,LocalNumAtoms
      ig=getGlobalIndex(i,MyPE)
!     ----------------------------------------------------------------
      call genPolyhedron(i,ig,NumAtoms,AtomPosition)
!     ----------------------------------------------------------------
      call getGridData(i,ndivin,ndivout,nmult)
!     ----------------------------------------------------------------
      if (print_level >= 2) then
!        -------------------------------------------------------------
         call printPolyhedron(i)
!        -------------------------------------------------------------
      endif
      write(6,*)'GridData ::'
      write(6,*)'id, rmt, rend, ndivin, ndivout, nmult :: '
      write(6,*) i, getInscrSphRadius(i), getOutscrSphRadius(i),        &
                ndivin,ndivout,nmult
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
   allocate(fit_pot(LocalNumAtoms), lpot(LocalNumAtoms))
   allocate(fit_rho(LocalNumAtoms), lrho(LocalNumAtoms))
   allocate(lkkr(LocalNumAtoms), lphi(LocalNumAtoms))
   do i=1,LocalNumAtoms
      fit_pot(i)=0                 ! = 0, no fit;  = 1, fit.
      fit_rho(i)=0                 ! = 0, no fit;  = 1, fit.
      lkkr(i)=getKKRLmax(i)
      lphi(i)=getPhiLmax(i)
      lrho(i)=getRhoLmax(i)
      lpot(i)=getPotLmax(i)
   enddo
!
!  -------------------------------------------------------------------
   call initPotentialType(getPotentialTypeParam())
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,fit_pot,lpot,n_spin_pola,         &
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
   do i=1,NumAtoms
      rmt(i) = getInscrSphRadius(i)
   enddo
   do i = 1,NumAtoms
      rcut_swiss(i) = getScreenRcut(i)
      lmax(i) = getScreenLmax(i)
      v_swiss(i) = getScreenPotential(i)
   enddo
!  -------------------------------------------------------------------
   call initSingleScattering( NumAtoms, AtomicNumber, lmax, lphi, lpot,   &
                              n_spin_cant, 0, istop, print_level )
!  -------------------------------------------------------------------
!   rcut_swiss(:) = 4.8d0
   write(*,*)'lmax   ::',lmax
   write(*,*)'v_swiss::',v_swiss
   write(*,*)'rcut   ::',rcut_swiss
   write(*,*)'rmt    ::',rmt
   write(*,*)'bravais::',bravais
   write(*,*)'AtomPosition::',AtomPosition
!  -------------------------------------------------------------------
!   write(6,*)"tst_Screen :: initScreenStrConstModule"
   call initTauScreenKKRModule(bravais,NumAtoms,n_spin_cant,AtomPosition(:,:), &
         lmax(:), rmt(:), v_swiss(:), rcut_swiss(:), print_level,istop)
!  -------------------------------------------------------------------
!
   epsilon = (0.75382452000000D+00,0.10000000000000D-02)
   isRelativistic = .false.
!  -------------------------------------------------------------------
!  Calculates Tau0j blocks for atom 1
!
   call calScreenTau0J(1,epsilon,isRelativistic)
!  ===================================================================
!  retrieve Tmatrix of the system and 
!  calculates DeltaM = [T_sys^(-1)-Tswiss^(-1)]
!  ===================================================================
!  define the size of the matrices
!
   vshift = czero
   kmax = 0
   l_max = 0
   do i =1,NumAtoms
      kmax = max( kmax,(lmax(i)+1)**2 )
      l_max = max( l_max,lmax(i) )
   enddo
   sdim = n_spin_cant*kmax
!
   allocate( DeltaM(sdim,sdim), Tmatrix(sdim,sdim)) ! , TauK(sdim,sdim) )
!  ===================================================================
!  get the single scatter T-matrix
!  ===================================================================
   Tmatrix = czero
   do i = 1,n_spin_cant
      Tmatrix((i-1)*kmax+1:i*kmax,(i-1)*kmax+1:i*kmax) =               &
                             getTMatrixLocalSpinBlock(1,i,i)
   enddo
   call writeMatrix('Tmatrix', Tmatrix,sdim,sdim )
!  -------------------------------------------------------------------
!  needed for charge density calculations
!  -------------------------------------------------------------------
   DeltaM = getDeltaM( 1 )
   call printScreenTau0J(n_spin_cant,2)
!  -------------------------------------------------------------------
!  initialize BZ
!  -------------------------------------------------------------------
   AtomType = 1
   call initBZone( NumKMeshs, kGenScheme, Kdiv, Symmetrize, bravais,  &
                NumAtoms, AtomPosition, AtomType, istop, print_level )
!  -------------------------------------------------------------------
!  get the Tau00 - integration over special k-points
!  -------------------------------------------------------------------
!   call initScreenTauBZ(NumAtoms,n_spin_cant,1,1)
   TauK => getBZirr(NumKMeshs,NumAtoms,n_spin_cant,l_max,kmax,.false.)
   do i = 1,NumAtoms
      write(*,*)"  TauK :: block",i
      call writeMatrix('TauK', TauK((i-1)*sdim+1:i*sdim,(i-1)*sdim+1:i*sdim), &
                       sdim,sdim)
   enddo
   call writeMatrix('TauK', TauK, sdim*NumAtoms, sdim*NumAtoms)
!  -------------------------------------------------------------------
!   DeltaM = getDelta(n_spin_cant,1)
!  ===================================================================
!  final output
!  ===================================================================
   call printScreenTau0J(n_spin_cant,2)
!   call printScreenTau00(n_spin_cant)
!  -------------------------------------------------------------------
   deallocate( AtomPosition, AtomicNumber, AtomType )
   deallocate( lmax, rmt, rcut_swiss, v_swiss )
   deallocate( fit_pot, lpot, fit_rho, lrho, lkkr, lphi )
   deallocate( DeltaM, TauK, Tmatrix )
   nullify( DeltaM, TauK, Tmatrix )
!  -------------------------------------------------------------------
   call endAtom()
   call endSystem()
   call endTauScreenKKRModule()
   call endPotentialType()
!   call endScreenTauBZ()
   call endOutput()
   call endMPP()
!  -------------------------------------------------------------------
   stop 'Ok'
end program test_ScreenTau
