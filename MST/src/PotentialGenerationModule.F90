Module PotentialGenerationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
   use MathParamModule, only : CZERO, ZERO, ONE, TWO, THREE, FOUR,    &
                               FIVE, SIX, THIRD, HALF, SQRTM1, PI,    &
                               PI2, PI4, SQRT_PI, TEN2m6, TEN2m7,     &
                               TEN2m8, TEN2m9, TEN2m10, TEN2m12, TEN2m13, Y0
   use IntegerFactorsModule, only : lofk, mofk, lofj, mofj, kofj, m1m
   use PublicTypeDefinitionsModule, only : GridStruct, UniformGridStruct
   use TimerModule, only : getTime
   use ChebyshevModule, only : ChebyshevStruct
!
   implicit none
!
public :: initPotentialGeneration, &
          endPotentialGeneration,  &
          computeNewPotential,     &
          getPotential,            &
          getPotentialAtPoint,     &
          getPotLmax,              &
          isPotComponentZero,      &
          getPotComponentFlag,     &
          getNumPotComponents,     &
          getSphPotr,              &
          getVdif,                 &
          setVdif,                 &
          getVshift,               &
          getVcoulomb_R0,          &
          getPseudoDipoleField,    &
          printMadelungShiftTable, &
          printPot_L,              &
          printPotentialGeneration
!
   interface getPotential
      module procedure getPot_L, getPot_Lj
   end interface
!
private
   type NewPotentialStruct
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind) :: lmax
      integer (kind=IntKind) :: jmax
      integer (kind=IntKind) :: jmt
      integer (kind=IntKind) :: jend
      integer (kind=IntKind) :: n_Rpts
      integer (kind=IntKind) :: ifit_XC
      integer (kind=IntKind) :: NumFlagJl
      real (kind=RealKind) :: Madelung_Shift
      real (kind=RealKind) :: VcoulombR0(3)
!
      type (GridStruct), pointer :: Grid
!
      real (kind=RealKind), pointer :: potr_sph(:,:,:)
!
      complex (kind=CmplxKind), pointer :: potL(:,:,:,:) ! = potL_Coulomb + potL_Exch
      complex (kind=CmplxKind), pointer :: potL_Tilda(:,:,:) ! Intra site potential
      complex (kind=CmplxKind), pointer :: potL_Madelung(:,:)
      complex (kind=CmplxKind), pointer :: potL_Pseudo(:,:)
      complex (kind=CmplxKind), pointer :: potL_Coulomb(:,:,:) ! = potL_Tilda + potL_Madelung + potL_Pseudo
!
      complex (kind=CmplxKind), pointer :: potL_Exch(:,:,:,:)
      complex (kind=CmplxKind), pointer :: enL_Exch(:,:,:,:)
!
      complex (kind=CmplxKind), pointer :: potL_XCHat(:,:,:,:) ! Seems not used
      complex (kind=CmplxKind), pointer :: enL_XCHat(:,:,:,:) ! Seems not used
!
      integer (kind=IntKind), pointer :: PotCompFlag(:)
      integer (kind=IntKind), pointer :: potL_Tilda_flag(:)
      integer (kind=IntKind), pointer :: potL_Madelung_flag(:)
      integer (kind=IntKind), pointer :: potL_Pseudo_flag(:)
      integer (kind=IntKind), pointer :: potL_Coulomb_flag(:)
      integer (kind=IntKind), pointer :: potL_Exch_flag(:)
      integer (kind=IntKind), pointer :: potL_XCHat_flag(:) ! Seems not used
!
   end type NewPotentialStruct
!
   type (NewPotentialStruct), allocatable, target :: Potential(:)
!
   logical :: Initialized = .false.
   logical :: EmptyTable = .true.
   logical :: isFullPot = .false.
   logical :: isMTFP = .false.
   logical :: isSphPot = .false.
   logical :: isChargeSymmOn = .false.
   logical :: isFirstExchg = .true.
   logical :: vshift_switch_on = .false.
   logical :: gga_functional = .false.
!
   character (len=50) :: StopRoutine
!
   integer (kind=IntKind) :: jend_max
   integer (kind=IntKind) :: lmax_max
   integer (kind=IntKind) :: jmax_max
   integer (kind=IntKind) :: lmax_rho_max
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: PotentialType
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind) :: n_inter = 5
   integer (kind=IntKind) :: InitMode = 0
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: Print_Level(:)
   integer (kind=IntKind), allocatable :: indrl_fit(:,:)
!
   real (kind=RealKind), allocatable :: alpha_mad(:)
   real (kind=RealKind), allocatable :: MadelungShiftTable(:)
   real (kind=RealKind) :: rhoint_sav, MadelungSum
   real (kind=RealKind), target :: vdif(1)
   real (kind=RealKind) :: v0(2)
   real (kind=RealKind) :: v_shift(2)
   real (kind=RealKind) :: MuffinTinZeroCoulumb
   real (kind=RealKind) :: MuffinTinZeroExchCorr(2)
!
   real (kind=RealKind), allocatable :: sqrt_r(:), V0_inter(:)
   real (kind=RealKind), allocatable, target :: rho_tmp_r(:), mom_tmp_r(:)
   real (kind=RealKind), allocatable, target :: rho_tmp_i(:), mom_tmp_i(:)
   real (kind=RealKind), allocatable, target :: V1_r(:), V2_r(:)
   real (kind=RealKind), allocatable, target :: E1_r(:), E2_r(:)
   real (kind=RealKind), allocatable, target :: V1_i(:), V2_i(:)
!
   complex (kind=CmplxKind), allocatable :: pXC_r(:,:),pXC_rphi(:,:)
   complex (kind=CmplxKind), allocatable :: eXC_r(:,:),eXC_rphi(:,:)
   complex (kind=CmplxKind), allocatable, target :: rho_tmp_c(:), mom_tmp_c(:)
   complex (kind=CmplxKind), allocatable, target :: V1_c(:), V2_c(:)
   complex (kind=CmplxKind), allocatable, target :: E1_c(:), E2_c(:)
!
   real (kind=RealKind), allocatable, target :: DF_Pseudo(:,:)
!
   real (kind=RealKind), parameter :: rmt_fraction = 0.70d0
   real (kind=RealKind), parameter :: pot_tol = TEN2m8
!
   real (kind=RealKind), allocatable :: vmt1(:)
!
   type AngularIntegrationData
      integer (kind=IntKind) :: kmax
      integer (kind=IntKind) :: n_theta
      integer (kind=IntKind) :: n_phi
      integer (kind=IntKind) :: ngl
      real (kind=RealKind), pointer:: radial_data(:,:,:,:,:)
      real (kind=RealKind), pointer :: theta(:), phi(:)
      real (kind=RealKind), pointer :: wght_the(:), wght_phi(:)
      complex (kind=CmplxKind), pointer :: ylm_ngl(:,:)
   end type AngularIntegrationData
!
   type (AngularIntegrationData) :: AngularData
!
   integer (kind=IntKind) :: numk_local
   integer (kind=IntKind), parameter :: nr_int_max = 3
   integer (kind=IntKind), parameter :: n_interp_max = 100
   integer (kind=IntKind), allocatable :: iparam(:,:)
!
   complex(kind=CmplxKind), pointer :: fft_c(:)
   complex(kind=CmplxKind), allocatable :: Ylm(:)
   complex(kind=CmplxKind), allocatable, target :: v_interp(:,:)
   real (kind=RealKind), allocatable, target ::  LocalAtomPosi(:,:)
   type (ChebyshevStruct), allocatable, target :: chebv_struct(:)
!
   interface
      subroutine constructDataOnGrid(grid_name, value_name, value_type, getData, den, lmax, spin)
         use KindParamModule, only : IntKind, RealKind
         use PublicTypeDefinitionsModule, only : UniformGridStruct
         implicit none
         character (len=*), intent(in) :: grid_name
         character (len=*), intent(in) :: value_name
         character (len=*), intent(in) :: value_type
         real (kind=RealKind), intent(out) :: den(:)
         integer (kind=IntKind), intent(in), optional :: lmax, spin
!
         interface
            function getData( dname, id, ia, r, jmax_in, n, grad ) result(v)
               use KindParamModule, only : IntKind, RealKind
               implicit none
               character (len=*), intent(in) :: dname
               integer (kind=IntKind), intent(in) :: id, ia
               real (kind=RealKind), intent(in) :: r(3)
               real (kind=RealKind), intent(out), optional :: grad(3)
               integer (kind=IntKind), intent(in), optional :: jmax_in, n
               real (kind=RealKind) :: v
            end function getData
         end interface
!
      end subroutine constructDataOnGrid
   end interface
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
contains
!  ===================================================================
   include '../lib/arrayTools.F90'
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initPotentialGeneration(nlocal,num_atoms,lmax_p,lmax_r, &
                                      npola,istop,iprint,isGGA)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalSumInGroup, GlobalMinInGroup, GlobalMaxInGroup
!
   use ParallelFFTModule, only : initParallelFFT, allocateFunctionSpace
   use ParallelFFTModule, only : getGlobalGridIndexRange, getProcessorMesh
!
   use RadialGridModule, only : getGrid, getNumRmesh, getGridRadius
!
   use PolyhedraModule, only : getVolume
   use PolyhedraModule, only : getPointLocationFlag, getOutscrSphRadius
!
   use MadelungModule, only : getMadelungMatrix
!
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       setDataStorage2Value,  &
                                       setDataStorageLDA,     &
                                       getDataStorageLDA,     &
                                       isDataStorageExisting, &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark
!
   use PotentialTypeModule, only : isFullPotential, isSphericalPotential, &
                                   isMuffinTinFullPotential
#ifdef POT_DEBUG
!
   use SurfElementsModule, only : initSurfElements
   use SurfElementsModule, only : genGaussPoints, genGaussSphHarmonics
#endif
!
   use SystemModule, only : getLmaxRho, getLmaxMax
   use SystemModule, only : getUniformGridParam, getBravaisLattice
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use AtomModule, only : getLocalNumSpecies
   use AtomModule, only : getLocalAtomPosition
!
   use ScfDataModule, only : isChargeSymm
!
   use Uniform3DGridModule, only : initUniform3DGrid
   use Uniform3DGridModule, only : isUniform3DGridInitialized, createUniform3DGrid, printUniform3DGrid
   use Uniform3DGridModule, only : createProcessorMesh, getUniform3DGrid
   use Uniform3DGridModule, only : distributeUniformGrid, insertAtomsInGrid
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   logical, intent(in), optional :: isGGA
!
   integer (kind=IntKind), intent(in) :: nlocal,num_atoms
   integer (kind=IntKind), intent(in) :: lmax_p(nlocal)
   integer (kind=IntKind), intent(in) :: lmax_r(nlocal)
   integer (kind=IntKind), intent(in) :: iprint(nlocal)
   integer (kind=IntKind), intent(in) :: npola
!
   integer (kind=IntKind), allocatable :: DataSize(:), DataSizeL(:)
   integer (kind=IntKind), allocatable :: SpeciesDataSize(:), SpeciesDataSizeL(:)
!
   type (GridStruct), pointer :: Grid
!
   integer (kind=IntKind) :: jl, kl, l, m, ig, id, nr, ir
   integer (kind=IntKind) :: kmax_max, lmax_pot, jmax_pot
   integer (kind=IntKind) :: jmt, jend, num_species
   integer (kind=IntKind) :: grid_start(3), grid_end(3), gir(3,3), ng_uniform(3)
!
   real (kind=RealKind) :: alpha, r, alpha_min, alpha_max
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), pointer :: madmat(:), r_mesh(:)
   real (kind=RealKind), allocatable :: radius(:)
!
   type (NewPotentialStruct), pointer :: p_Pot
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   LocalNumAtoms = nlocal
   GlobalNumAtoms = num_atoms
   n_spin_pola = npola
   isFullPot = isFullPotential()
   isMTFP = isMuffinTinFullPotential()
   isSphPot = isSphericalPotential()
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   allocate(Print_Level(nlocal))
   Print_Level(1:nlocal) = iprint(1:nlocal)
   allocate( DataSize(LocalNumAtoms), SpeciesDataSize(LocalNumAtoms))
   allocate( DataSizeL(LocalNumAtoms), SpeciesDataSizeL(LocalNumAtoms) )
!
   StopRoutine = istop
!
   if (present(isGGA)) then
      gga_functional = isGGA
   else
      gga_functional = .false.
   endif
!
   jend_max = 0 
   lmax_max = 0
   node_print_level = -1
   do id = 1, LocalNumAtoms
      lmax_max=max(lmax_max, lmax_p(id), lmax_r(id))
      DataSize(id) = getNumRmesh(id)
      SpeciesDataSize(id) = DataSize(id)*getLocalNumSpecies(id)
      DataSizeL(id) = getNumRmesh(id)*(((lmax_p(id)+1)*(lmax_p(id)+2))/2)
      SpeciesDataSizeL(id) = DataSizeL(id)*getLocalNumSpecies(id)
      jend_max=max(jend_max,DataSize(id))
      node_print_level = max(node_print_level,Print_Level(id))
   enddo
   lmax_rho_max = getLmaxRho()
   lmax_max = max(lmax_max,lmax_rho_max)
   lmax_max = max(lmax_max,getLmaxMax())
   jmax_max = ((lmax_max+1)*(lmax_max+2))/2
   kmax_max = (lmax_max+1)*(lmax_max+1)
!
!  ------------------------------------------------------------------
   call initIntegerFactors(lmax_max)
!  ------------------------------------------------------------------
   allocate(vmt1(LocalNumAtoms))
   if ( isFullPot ) then
      allocate(V0_inter(n_spin_pola) )
   endif
!
!  ==================================================================
!  generate Storage space for Potentials
!  ==================================================================
   if (.not.isDataStorageExisting('NewSphericalPotential')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'NewSphericalPotential',    &
                             SpeciesDataSize*n_spin_pola,RealType)
!     ----------------------------------------------------------------
      call setDataStorage2Value('NewSphericalPotential',ZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'NewSphericalPotential',getNumRmesh(id))
!        -------------------------------------------------------------
      enddo
   endif
   if (.not.isDataStorageExisting('L-Potential')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'L-Potential',             &
                             SpeciesDataSizeL*n_spin_pola,ComplexType)
!     ----------------------------------------------------------------
      call setDataStorage2Value('L-Potential',CZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'L-Potential',getNumRmesh(id))
!        -------------------------------------------------------------
      enddo
   endif
   if ( isFullPot ) then
      if (.not.isDataStorageExisting('TildaPotential')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'TildaPotential',        &
                                SpeciesDataSizeL,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('TildaPotential',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'TildaPotential',getNumRmesh(id))
!           ----------------------------------------------------------
         enddo
      endif
      if (.not.isDataStorageExisting('MadelungPotential')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'MadelungPotential',     &
                                DataSizeL,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('MadelungPotential',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'MadelungPotential',getNumRmesh(id))
!           ----------------------------------------------------------
         enddo
      endif
      if (.not.isDataStorageExisting('CoulombPotential')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'CoulombPotential',      &
                                SpeciesDataSizeL,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('CoulombPotential',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'CoulombPotential',getNumRmesh(id))
!           ----------------------------------------------------------
         enddo
      endif
      if (.not.isDataStorageExisting('XchgCorrPotential')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'XchgCorrPotential',     &
                                SpeciesDataSizeL*n_spin_pola,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('XchgCorrPotential',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'XchgCorrPotential',getNumRmesh(id))
!           ----------------------------------------------------------
         enddo
      endif
      if (.not.isDataStorageExisting('XchgCorrEnergy')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'XchgCorrEnergy',        &
                                SpeciesDataSizeL*n_spin_pola,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('XchgCorrEnergy',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'XchgCorrEnergy',getNumRmesh(id))
!           ----------------------------------------------------------
         enddo
      endif
      if (.not.isMTFP) then
         if (.not.isDataStorageExisting('PseudoPotential')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'PseudoPotential',       &
                                   DataSizeL,ComplexType)
!           ----------------------------------------------------------
            call setDataStorage2Value('PseudoPotential',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
!              -------------------------------------------------------
               call setDataStorageLDA(id,'PseudoPotential',getNumRmesh(id))
!              -------------------------------------------------------
            enddo
         endif
         if (.not.isDataStorageExisting('XchgCorrHatPotential')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'XchgCorrHatPotential',  &
                                   SpeciesDataSizeL*n_spin_pola,ComplexType)
!           ----------------------------------------------------------
            call setDataStorage2Value('XchgCorrHatPotential',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
!              -------------------------------------------------------
               call setDataStorageLDA(id,'XchgCorrHatPotential',getNumRmesh(id))
!              -------------------------------------------------------
            enddo
         endif
         if (.not.isDataStorageExisting('XchgCorrHatEnergy')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'XchgCorrHatEnergy',     &
                                   SpeciesDataSizeL*n_spin_pola,ComplexType)
!           ----------------------------------------------------------
            call setDataStorage2Value('XchgCorrHatEnergy',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
!              -------------------------------------------------------
               call setDataStorageLDA(id,'XchgCorrHatEnergy',getNumRmesh(id))
!              -------------------------------------------------------
            enddo
         endif
      endif
!
!ywg 11/114/18
!ywg  call setAngularData(60,40,lmax_max,npola,jend_max)
  call setAngularData(80,50,lmax_max,npola,jend_max)
!......
      allocate( indrl_fit(jmax_max,npola) )
      allocate( pXC_r(jmax_max,npola),pXC_rphi(jmax_max,npola) )
      allocate( eXC_r(jmax_max,npola),eXC_rphi(jmax_max,npola) )
!
   endif
!
   jend_max = 0
   allocate( Potential(LocalNumAtoms) )
!
   isChargeSymmOn = isChargeSymm()
!
   do id = 1,LocalNumAtoms
      num_species = getLocalNumSpecies(id)
      Potential(id)%NumSpecies = num_species
      p_Pot => Potential(id)
      lmax_pot = lmax_p(id)
      Grid => getGrid(id)
      p_Pot%Grid   => Grid
      jmax_pot = ((lmax_pot+1)*(lmax_pot+2))/2
      jend = Grid%jend
      jend_max = max(jend_max,jend)
      jmt = Grid%jmt
      r_mesh => Grid%r_mesh(1:jend)
      r = rmt_fraction*r_mesh(jmt)
!     ----------------------------------------------------------------
      call hunt(jmt,r_mesh(1:jmt),r,ir)
!     ----------------------------------------------------------------
      p_Pot%ifit_XC = ir
!     ================================================================
      p_Pot%potr_sph => getDataStorage( id, 'NewSphericalPotential',  &
                                        jend, n_spin_pola, num_species, RealMark )
!     ================================================================
      if ( isFullPot ) then
         p_Pot%potL => getDataStorage( id, 'L-Potential',             &
                                       jend, jmax_pot, n_spin_pola, num_species, ComplexMark )
         p_Pot%potL_Tilda => getDataStorage( id, 'TildaPotential',    &
                                             jend, jmax_pot, num_species, ComplexMark )
         p_Pot%potL_Madelung => getDataStorage( id, 'MadelungPotential',&
                                                jend, jmax_pot, ComplexMark )
         p_Pot%potL_Coulomb => getDataStorage( id, 'CoulombPotential',  &
                                               jend, jmax_pot, num_species, ComplexMark )
         p_Pot%potL_Exch => getDataStorage( id, 'XchgCorrPotential',  &
                                            jend, jmax_pot, n_spin_pola, num_species, ComplexMark )
         p_Pot%enL_Exch => getDataStorage( id, 'XchgCorrEnergy',      &
                                           jend, jmax_pot, n_spin_pola, num_species, ComplexMark )
         if (.not.isMTFP) then
            p_Pot%potL_Pseudo => getDataStorage( id, 'PseudoPotential',    &
                                                 jend, jmax_pot, ComplexMark )
            p_Pot%potL_XCHat => getDataStorage( id, 'XchgCorrHatPotential',&
                                                jend, jmax_pot, n_spin_pola, num_species, ComplexMark )
            p_Pot%enL_XCHat  => getDataStorage( id, 'XchgCorrHatEnergy',   &
                                                jend, jmax_pot, n_spin_pola, num_species, ComplexMark )
         endif
      else
         p_Pot%potL => getDataStorage( id, 'L-Potential',             &
                                       jend, jmax_pot, n_spin_pola, num_species, ComplexMark )
         nullify( p_Pot%potL_Tilda, p_Pot%potL_Madelung, p_Pot%potL_Exch, &
                  p_Pot%potL_Pseudo, p_Pot%potL_Coulomb, p_Pot%potL_XCHat )
         nullify( p_Pot%enL_Exch, p_Pot%enL_XCHat )
      endif
!
      allocate( p_Pot%PotCompFlag(jmax_pot) )
      allocate( p_Pot%potL_Tilda_flag(jmax_pot) )
      allocate( p_Pot%potL_Madelung_flag(jmax_pot) )
      allocate( p_Pot%potL_Coulomb_flag(jmax_pot) )
      allocate( p_Pot%potL_Exch_flag(jmax_pot) )
      allocate( p_Pot%potL_Pseudo_flag(jmax_pot) )
      allocate( p_Pot%potL_XCHat_flag(jmax_pot) )
!
      p_Pot%PotCompFlag = 0
      p_Pot%potL_Tilda_flag = 0
      p_Pot%potL_Madelung_flag = 0
      p_Pot%potL_Coulomb_flag = 0 
      p_Pot%potL_Exch_flag = 0
!
      p_Pot%PotCompFlag(1) = 1
      p_Pot%potL_Tilda_flag(1) = 1
      p_Pot%potL_Madelung_flag(1) = 1
      p_Pot%potL_Coulomb_flag(1) = 1 
      p_Pot%potL_Exch_flag(1) = 1
!
      p_Pot%potL_Pseudo_flag = 0
      p_Pot%potL_XCHat_flag = 0
      p_Pot%potL_Pseudo_flag(1) = 1
      p_Pot%potL_XCHat_flag(1) = 1
!
      p_Pot%lmax   = lmax_pot
      p_Pot%jmax   = jmax_pot
      p_Pot%jmt    = jmt
      p_Pot%jend   = jend
      p_Pot%n_Rpts = jend
   enddo
!
   allocate( alpha_mad(LocalNumAtoms), GlobalIndex(LocalNumAtoms) )
   allocate( LocalAtomPosi(3,LocalNumAtoms) )
   allocate( MadelungShiftTable(GlobalNumAtoms) )
   allocate( radius(LocalNumAtoms) )
!
   alpha_min =  1.0d+20
   alpha_max = -1.0d+20
   MadelungSum = ZERO
   radius = TEN2m6
   do id = 1, LocalNumAtoms
      madmat => getMadelungMatrix(id)
      alpha=ZERO
      do ig=1,GlobalNumAtoms
         alpha=alpha+madmat(ig)
      enddo
      alpha_min = min(alpha_min,alpha)
      alpha_max = max(alpha_max,alpha)
      MadelungSum = MadelungSum + alpha
      alpha_mad(id)=TWO*alpha*getVolume(id)
      GlobalIndex(id) = getGlobalIndex(id)
      LocalAtomPosi(1:3,id)=getLocalAtomPosition(id)
!     radius(id) = getOutscrSphRadius(id)
      radius(id) = getGridRadius(id)
   enddo
!
!  ===================================================================
!  Check if the summation of each Madelung column is a constant.
!      If yes, the total Madelung shift in the unit cell will be zero,
!      one can take advantage of this fact to correct the numerical errors.
!  -------------------------------------------------------------------
   call GlobalMinInGroup(GroupID,alpha_min)
   call GlobalMaxInGroup(GroupID,alpha_max)
!  -------------------------------------------------------------------
   if (abs(alpha_min-alpha_max) > TEN2m6) then
      vshift_switch_on = .false.
   else
      vshift_switch_on = .true.
      MadelungSum = MadelungSum/real(GlobalNumAtoms,kind=RealKind)
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,MadelungSum)
!     ----------------------------------------------------------------
   endif
!
   deallocate(DataSize, SpeciesDataSize)
   deallocate( DataSizeL, SpeciesDataSizeL )
   if ( isFullPot ) then
#ifdef POT_DEBUG
      call initSurfElements(istop,-1)
      call genGaussPoints()
      call genGaussSphHarmonics(lmax_max)
#endif
   endif
!
   allocate( sqrt_r(0:jend_max) )
   allocate( rho_tmp_r(0:jend_max), mom_tmp_r(0:jend_max) )
   allocate( rho_tmp_i(0:jend_max), mom_tmp_i(0:jend_max) )
   allocate( V1_r(0:jend_max), V2_r(0:jend_max) )
   allocate( E1_r(0:jend_max), E2_r(0:jend_max) )
   allocate( V1_i(0:jend_max), V2_i(0:jend_max) )
   allocate( rho_tmp_c(0:jend_max), mom_tmp_c(0:jend_max) )
   allocate( V1_c(0:jend_max), V2_c(0:jend_max) )
   allocate( E1_c(0:jend_max), E2_c(0:jend_max) )
   if ( isFullPot .and. .not.isMTFP) then
      allocate( DF_Pseudo(3,LocalNumAtoms) )
   endif
!
   if ( isFullPot ) then
      if (.not.isUniform3DGridInitialized()) then
!        -------------------------------------------------------------
         call initUniform3DGrid(istop, node_print_level)
!        -------------------------------------------------------------
      endif
      ng_uniform(1:3) = getUniformGridParam()
      bravais = getBravaisLattice()
!     ----------------------------------------------------------------
      call createUniform3DGrid('FFT', ng_uniform(1), ng_uniform(2), ng_uniform(3), bravais)
!     ----------------------------------------------------------------
      call initParallelFFT()
!     ----------------------------------------------------------------
      call createProcessorMesh('FFT', getProcessorMesh())
!     ----------------------------------------------------------------
      gir = getGlobalGridIndexRange('R')
!     ----------------------------------------------------------------
      grid_start(1:3) = gir(1:3,1)
      grid_end(1:3) = gir(1:3,2)
!     ----------------------------------------------------------------
      call distributeUniformGrid('FFT',grid_start,grid_end)
!     ----------------------------------------------------------------
      call insertAtomsInGrid('FFT', LocalNumAtoms, LocalAtomPosi,     &
                             getPointLocationFlag, radius)
!     ----------------------------------------------------------------
      if (node_print_level >= 0) then
         call printUniform3DGrid('FFT')
      endif
!     ----------------------------------------------------------------
      call allocateFunctionSpace(fft_c,numk_local)
!     ----------------------------------------------------------------
!
      allocate( Ylm(kmax_max) )
      allocate( iparam(4+nr_int_max,LocalNumAtoms) )
      allocate( v_interp(n_interp_max*nr_int_max*jmax_max,LocalNumAtoms) )
      allocate( chebv_struct(jmax_max) )
   endif
!
   deallocate( radius )
!
   Initialized = .true.
   EmptyTable = .true.
   isFirstExchg = .true.
!
   end subroutine initPotentialGeneration
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endPotentialGeneration()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
#ifdef POT_DEBUG
   use SurfElementsModule, only : endSurfElements
#endif
!
   use ParallelFFTModule, only : endParallelFFT
!
   implicit none
!
   integer (kind=IntKind) :: n
!
   do n = 1,LocalNumAtoms
      deallocate( Potential(n)%PotCompFlag )
      deallocate( Potential(n)%potL_Tilda_flag )
      deallocate( Potential(n)%potL_Madelung_flag )
      deallocate( Potential(n)%potL_Pseudo_flag )
      deallocate( Potential(n)%potL_Coulomb_flag )
      deallocate( Potential(n)%potL_Exch_flag )
      deallocate( Potential(n)%potL_XCHat_flag )
      nullify( Potential(n)%potr_sph)
      nullify( Potential(n)%Grid )
      nullify( Potential(n)%potL )
      if ( isFullPot ) then ! .or. .not.isSphPot ) then
         nullify( Potential(n)%potL_Tilda )
         nullify( Potential(n)%potL_Madelung )
         nullify( Potential(n)%potL_Coulomb )
         nullify( Potential(n)%potL_Exch )
         nullify( Potential(n)%enL_Exch )
         if (.not.isMTFP) then
            nullify( Potential(n)%potL_Pseudo )
            nullify( Potential(n)%potL_XCHat )
            nullify( Potential(n)%enL_XCHat )
         endif
      endif
   enddo
!
   deallocate( MadelungShiftTable )
   deallocate( Potential, alpha_mad, GlobalIndex, Print_Level )
!
   deallocate( sqrt_r )
   deallocate( rho_tmp_i, mom_tmp_i )
   deallocate( rho_tmp_r, mom_tmp_r )
   deallocate( V1_r, V2_r, V1_i, V2_i )
   deallocate( rho_tmp_c, mom_tmp_c, V1_c, V2_c )
   deallocate( vmt1 )
   if ( isFullPot ) then
      deallocate(V0_inter)
      if (.not.isMTFP) then
         deallocate( DF_Pseudo )
      endif
      deallocate( Ylm )
      deallocate( v_interp )
      deallocate( chebv_struct )
   endif
!
#ifdef POT_DEBUG
!  -------------------------------------------------------------------
   call endSurfElements()
!  -------------------------------------------------------------------
#endif
!  -------------------------------------------------------------------
   call endIntegerFactors()
!  -------------------------------------------------------------------
   if ( isFullPot ) then
!     ----------------------------------------------------------------
      call clearAngularData()
      call endParallelFFT()
!     ----------------------------------------------------------------
      deallocate( fft_c ); nullify( fft_c )
      deallocate( indrl_fit )
      deallocate( pXC_r, pXC_rphi, eXC_r, eXC_rphi )
   endif
   Initialized = .false.
   EmptyTable = .true.
   isChargeSymmOn = .false.
!
   end subroutine endPotentialGeneration
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphPotr(id,ia,is)                         result(potr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind) :: npts
!
   real (kind=RealKind), pointer :: potr(:)
!
   if ( .not.Initialized ) then
      call ErrorHandler("getSphrPotr",                                &
             "The PotentialGenerationModule has to be initialized first")
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphPotr','invalid id',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSphPotr','invalid spin index',is)
   else if (ia < 1 .or. ia > Potential(id)%NumSpecies) then
      call ErrorHandler('getSphPotr','invalid species index',ia)
   endif
!
   npts = Potential(id)%Grid%jend
   potr => Potential(id)%potr_sph(1:npts,is,ia)
!
   end function getSphPotr
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotLmax(id)                                result(lmax)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: lmax
!
   if ( .not.Initialized ) then
      call ErrorHandler("getPotLmax",                                &
             "The PotentialGenerationModule has to be initialized first")
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPotLmax','invalid id',id)
   endif
!
   lmax = Potential(id)%lmax
!
   end function getPotLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVdif() result(vd)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: vd(:)
!
   vd => vdif
!
   end function getVdif
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVshift(is) result(vs)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind) :: vs
!
   vs = v_shift(is)
!
   end function getVshift
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVcoulomb_R0(id) result(vc_0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: vc_0(3)
!
   vc_0 = Potential(id)%VcoulombR0
!
   end function getVcoulomb_R0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setVdif(vd)
!  ===================================================================
   implicit none
   real (kind=RealKind), intent(in) :: vd
!
   vdif(1) = vd
!
   end subroutine setVdif
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isPotComponentZero(id,jl) result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: jl
!
   logical :: flag
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('isPotComponentZero','invalid id',id)
   else if (jl < 1 .or. jl > Potential(id)%jmax) then
      flag = .true.
   else if (Potential(id)%PotCompFlag(jl) == 0) then
      flag = .true.
   else
      flag = .false.
   endif
!
   end function isPotComponentZero
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotComponentFlag(id)                    result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer, pointer :: flag(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPotComponentZero','invalid id',id)
   endif
!
   flag => Potential(id)%PotCompFlag(:)
!
   end function getPotComponentFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setPotComponentFlag()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind):: id, is, nr, jl, ia
!
   integer, pointer :: flag(:)
!
   do id = 1,LocalNumAtoms
      flag => Potential(id)%PotCompFlag(1:Potential(id)%jmax)
      flag = 0
      do ia = 1, Potential(id)%NumSpecies
         do is = 1,n_spin_pola
            do jl = 1,Potential(id)%jmax
               LOOP_NR: do nr = 1,Potential(id)%n_Rpts
!                 if ( abs(Potential(id)%potL(nr,jl,is)) > ten2m10) then
                  if ( abs(Potential(id)%potL(nr,jl,is,ia)) > pot_tol) then
                     flag(jl) = 1
                     cycle LOOP_NR
                  endif
               enddo LOOP_NR
            enddo
         enddo
      enddo
   enddo
!
   end subroutine setPotComponentFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumPotComponents(id)                   result(ncomp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: ncomp
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getNumPotComponents','invalid id',id)
   endif
!
   ncomp = Potential(id)%NumFlagJl
!
   end function getNumPotComponents
!  ===================================================================

!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotentialAtPosi(potentialType,site,atom,is,posi) result(pot)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: PotentialType
!
   integer (kind=IntKind), intent(in) :: site, atom, is
!
   real (kind=RealKind), intent(in) :: posi(3)
   real (kind=RealKind) :: pot
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
!
   interface
      function getValueAtPosi(r_mesh,posi,val_l) result(val)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         real (kind=RealKind), intent(in) :: r_mesh(:), posi(3)
         real (kind=RealKind) :: val
         complex (kind=CmplxKind), intent(in) :: val_l(:,:)
      end function getValueAtPosi
   end interface
!
   if (site < 1 .or. site > LocalNumAtoms) then
      call ErrorHandler('getPotential','invalid site',site)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getPotential','invalid is',is)
   else if (atom < 1 .or. atom > Potential(site)%NumSpecies) then
      call ErrorHandler('getPotential','invalid species index',atom)
   else if (isMTFP .and. PotentialType=="Pseudo" ) then
      call ErrorHandler('getPotential','invalid potential type','Pseudo')
   endif
!
   if ( PotentialType=="Total" .or. PotentialType=="Potential" ) then
      if ( sqrt(posi(1)**2+posi(2)**2+posi(3)**2) < TEN2m8 ) then
         pot = -1.0d12
         return
      endif
      pot_l => Potential(site)%potL(:,:,is,atom)
   else if ( PotentialType=="Pseudo" ) then
      pot_l => Potential(site)%potL_Pseudo(:,:)
   else if ( PotentialType=="Tilda" ) then
      pot_l => Potential(site)%potL_Tilda(:,:,atom)
   else if ( PotentialType=="Coulomb" ) then
      if ( sqrt(posi(1)**2+posi(2)**2+posi(3)**2) < TEN2m8 ) then
         pot = -1.0d12
         return
      endif
      pot_l => Potential(site)%potL_Coulomb(:,:,atom)
   else if ( PotentialType=="Madelung" ) then
      pot_l => Potential(site)%potL_Madelung(:,:)
   else if ( PotentialType=="Exchg" ) then
      pot_l => Potential(site)%potL_Exch(:,:,is,atom)
   else if ( PotentialType=="En_Exchg" ) then
      pot_l => Potential(site)%enL_Exch(:,:,is,atom)
   else
      call ErrorHandler("getPot_L","Undefined potential type")
   endif
!
!  -------------------------------------------------------------------
   pot = getValueAtPosi(Potential(site)%Grid%r_mesh,posi,pot_l)
!  -------------------------------------------------------------------
!
   end function getPotentialAtPosi
!  ===================================================================

!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPot_Lj(potentialType,site,atom,is,jl) result(pot_l)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: PotentialType
!
   integer (kind=IntKind), intent(in) :: site, atom
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind) :: iend
!
   complex (kind=CmplxKind), pointer :: pot_l(:)
!
   if (site < 1 .or. site > LocalNumAtoms) then
      call ErrorHandler('getPotential','invalid site',site)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getPotential','invalid is',is)
   else if (atom < 1 .or. atom > Potential(site)%NumSpecies) then
      call ErrorHandler('getPotential','invalid species index',atom)
   else if (jl < 1 .or. jl > Potential(site)%jmax) then
      call ErrorHandler('getPotential','invalid jl',jl)
   else if (isMTFP .and. nocaseCompare(PotentialType,"Pseudo") ) then
      call ErrorHandler('getPot_Lj','invalid potential type','Pseudo')
   endif
!
   iend = Potential(site)%Grid%jend
   if ( nocaseCompare(PotentialType,"Total") ) then
      pot_l => Potential(site)%potL(1:iend,jl,is,atom)
   else if ( nocaseCompare(PotentialType,"Pseudo") ) then
      pot_l => Potential(site)%potL_Pseudo(1:iend,jl)
   else if ( nocaseCompare(PotentialType,"Tilda") ) then
      pot_l => Potential(site)%potL_Tilda(1:iend,jl,atom)
   else if ( nocaseCompare(PotentialType,"Madelung") ) then
      pot_l => Potential(site)%potL_Madelung(1:iend,jl)
   else if ( nocaseCompare(PotentialType,"Coulomb") ) then
      pot_l => Potential(site)%potL_Coulomb(1:iend,jl,atom)
   else if ( nocaseCompare(PotentialType,"Exchg") ) then
      pot_l => Potential(site)%potL_Exch(1:iend,jl,is,atom)
   else if ( nocaseCompare(PotentialType,"En_Exchg") ) then
      pot_l => Potential(site)%enL_Exch(1:iend,jl,is,atom)
   else
      call ErrorHandler("getPot_Lj","Undefined potential type")
   endif
!
   end function getPot_Lj
!  ===================================================================

!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPot_L(potentialType,site,atom,is) result(pot_l)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: PotentialType
!
   integer (kind=IntKind), intent(in) :: site, atom
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind) :: iend, jmax
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
!
   if (site < 1 .or. site > LocalNumAtoms) then
      call ErrorHandler('getPotential','invalid site',site)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getPotential','invalid is',is)
   else if (atom < 1 .or. atom > Potential(site)%NumSpecies) then
      call ErrorHandler('getPotential','invalid species index',atom)
   else if (isMTFP .and. nocaseCompare(PotentialType,"Pseudo") ) then
      call ErrorHandler('getPot_L','invalid potential type','Pseudo')
   endif
!
   iend = Potential(site)%Grid%jend
   jmax = Potential(site)%jmax
   if ( nocaseCompare(PotentialType,"Total") ) then
      pot_l => Potential(site)%potL(1:iend,1:jmax,is,atom)
   else if ( nocaseCompare(PotentialType,"Pseudo") ) then
      pot_l => Potential(site)%potL_Pseudo(1:iend,1:jmax)
   else if ( nocaseCompare(PotentialType,"Tilda") ) then
      pot_l => Potential(site)%potL_Tilda(1:iend,1:jmax,atom)
   else if ( nocaseCompare(PotentialType,"Coulomb") ) then
      pot_l => Potential(site)%potL_Coulomb(1:iend,1:jmax,atom)
   else if ( nocaseCompare(PotentialType,"Madelung") ) then
      pot_l => Potential(site)%potL_Madelung(1:iend,1:jmax)
   else if ( nocaseCompare(PotentialType,"Exchg") ) then
      pot_l => Potential(site)%potL_Exch(1:iend,1:jmax,is,atom)
   else if ( nocaseCompare(PotentialType,"En_Exchg") ) then
      pot_l => Potential(site)%enL_Exch(1:iend,1:jmax,is,atom)
   else
      call ErrorHandler("getPot_L","Undefined potential type")
   endif
!
   end function getPot_L
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeNewPotential(isMT)
!  ===================================================================
   use MPPModule, only : MyPE
   use GroupCommModule, only : GlobalSumInGroup
   use Atom2ProcModule, only : getMaxLocalNumAtoms
   use Atom2ProcModule, only : getLocalIndex, getAtom2ProcInGroup
   use AtomModule, only : getLocalSpeciesContent
!
#ifdef POT_DEBUG
   use SurfElementsModule, only : getSurfAverage
#endif
!
   use PotentialTypeModule, only : isMuffintinASAPotential,            &
                                   isMuffintinPotential, isASAPotential
!
   use ChargeDistributionModule, only : getInterstitialElectronDensity
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use SystemVolumeModule, only : getSystemVolume, getTotalInterstitialVolume
!
   use SystemSymmetryModule, only : getSymmetryFlags
!
   use PolyhedraModule, only : getVolume, getInscrSphVolume
!
   implicit   none
!
   logical, optional :: isMT
   logical :: isMTon = .false.
!
   integer (kind=IntKind) :: id, ia, ig, ip, is, jl, jmt, ir, nRpts, jend, l
   integer (kind=IntKind) :: jmax, lmax, kmax
   integer (kind=IntKind) :: MaxLocalAtoms
   integer (kind=IntKind), pointer :: flag_jl(:), p_flags(:)
!
   real (kind=RealKind) :: pot_r, pot_i, vol_ints, vtmp, vs, vs_mt, rfac
#ifdef POT_DEBUG
   real (kind=RealKind) :: pot_aver(2), vol_int(6)
   real (kind=RealKind) :: potaver_tilda, potaver_hat,                 &
                           potaver_madelung, potaver_Coulomb
   real (kind=RealKind) :: membuf(6)
#endif
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind) :: t0, t1
!
   complex (kind=CmplxKind) :: cfact
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
!
#ifdef TIMING
   t0 = getTime()
#endif
   cfact = -SQRTm1
   if ( present(isMT) ) then
      isMTon = isMT
   else
      isMTon = .false.
   endif
   if ( .not.isFullPot .or. isMTon ) then
!
      rhoint_sav = getInterstitialElectronDensity()
!
!     ================================================================
!     calculate hartree and exchange-correlation potential
!     ----------------------------------------------------------------
      call calSphPotential()
!     ----------------------------------------------------------------
!
!     ================================================================
!     calculate zero potential
!     ================================================================
      if ( isMuffintinPotential() ) then
!        -------------------------------------------------------------
         call calZeroPotential_MT()
!        -------------------------------------------------------------
      else if ( isASAPotential() ) then
!        -------------------------------------------------------------
         call calZeroPotential_ASA()
!        -------------------------------------------------------------
      else if ( isMuffintinASAPotential() ) then
!        -------------------------------------------------------------
         call calZeroPotential_MTASA()
!        -------------------------------------------------------------
      endif
!     ================================================================
!     calculate new potential
!     ----------------------------------------------------------------
      call constructNewPotential()
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         nRpts = Potential(id)%n_Rpts
         r_mesh => Potential(id)%Grid%r_mesh(1:nRpts)
         do ia = 1, Potential(id)%NumSpecies
            do is = 1,n_spin_pola
               pot_l => Potential(id)%potL(1:nRpts,1:Potential(id)%jmax,is,ia)
               do ir = 1,nRpts
                  pot_l(ir,1) = Potential(id)%potr_sph(ir,is,ia)/(Y0*r_mesh(ir))
               enddo
            enddo
         enddo
      enddo
!
      MadelungShiftTable = ZERO
      do id = 1,LocalNumAtoms
         ig = GlobalIndex(id)
         MadelungShiftTable(ig) = Potential(id)%Madelung_Shift
      enddo
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,MadelungShiftTable,GlobalNumAtoms)
!     ----------------------------------------------------------------
      EmptyTable = .false.
#ifdef TIMING
      t0 = getTime()-t0
      write(6,'(a,f10.5)')' PotentialGeneration :: Time ',t0
#endif
   else 
      V0_inter = ZERO
      do id = 1,LocalNumAtoms
         Potential(id)%VcoulombR0 = ZERO
      enddo
      if (isMTFP) then
!        -------------------------------------------------------------
         call calIntraPot()
!        -------------------------------------------------------------
         call calInterPlusMadPot()
!        -------------------------------------------------------------
!        call calZeroPotential_MTFP()  ! Needs to look into...
         call calZeroPotential_MT()    ! We use calZeroPotential_MT instead
!        -------------------------------------------------------------
      else
         t1 = getTime()
!        -------------------------------------------------------------
         call calIntraPot()
!        -------------------------------------------------------------
         if (MyPE == 0) then
            write(6,'(/,a,f10.5,/)')'Time:: calIntraPot: ',getTime()-t1
         endif
         t1 = getTime()
!        -------------------------------------------------------------
         call calInterPlusMadPot()
!        -------------------------------------------------------------
         if (MyPE == 0) then
            write(6,'(/,a,f10.5,/)')'Time:: calInterPlusMadPot: ',getTime()-t1
         endif
         t1 = getTime()
!        -------------------------------------------------------------
         call calFFTPseudoPot()
!        -------------------------------------------------------------
         if (MyPE == 0) then
            write(6,'(/,a,f10.5,/)')'Time:: calFFTPseudoPot: ',getTime()-t1
         endif
      endif
!
      t1 = getTime()
      do id = 1,LocalNumAtoms
         do ia = 1, Potential(id)%NumSpecies
!           ----------------------------------------------------------
            call calExchangeJl(id,ia,-1)
!           ----------------------------------------------------------
         enddo
      enddo
      if (MyPE == 0) then
         write(6,'(/,a,f10.5,/)')'Time:: calExchangeJl: ',getTime()-t1
      endif
!
      t1 = getTime()
      do id = 1,LocalNumAtoms
         nRpts = Potential(id)%n_Rpts
         do ia = 1, Potential(id)%NumSpecies
            Potential(id)%potL = CZERO
            Potential(id)%potr_sph = ZERO
            do is = 1,n_spin_pola
               do jl = 1,Potential(id)%jmax
                  if ( is==1 ) then
                     if (isMTFP) then
                        Potential(id)%potL_Coulomb(1:nRpts,jl,ia) =          &
                                   Potential(id)%potL_Tilda(1:nRpts,jl,ia) + &
                                   Potential(id)%potL_Madelung(1:nRpts,jl)
                     else
                        Potential(id)%potL_Coulomb(1:nRpts,jl,ia) =          &
                                   Potential(id)%potL_Tilda(1:nRpts,jl,ia) + &
                                   Potential(id)%potL_Madelung(1:nRpts,jl) + &
                                   Potential(id)%potL_Pseudo(1:nRpts,jl)
                     endif
                  endif
                  Potential(id)%potL(1:nRpts,jl,is,ia) =                     &
                                Potential(id)%potL_Coulomb(1:nRpts,jl,ia) +  &
                                Potential(id)%potL_Exch(1:nRpts,jl,is,ia)
               enddo
!
               r_mesh => Potential(id)%Grid%r_mesh(1:nRpts)
               do ir = 1,nRpts
                  pot_r = real(Potential(id)%potL(ir,1,is,ia),kind=RealKind)
                  Potential(id)%potr_sph(ir,is,ia) = Y0*r_mesh(ir)*pot_r
               enddo
            enddo
         enddo
      enddo
      if (MyPE == 0) then
         write(6,'(/,a,f10.5,/)')'Time:: computeNewPotentail:: Loop_id 1: ',getTime()-t1
      endif
!
#ifdef POT_DEBUG
      vol_int = ZERO
!
      pot_aver(1:2) = ZERO
      potaver_tilda    = ZERO
      potaver_madelung = ZERO
      potaver_hat      = ZERO
      potaver_Coulomb  = ZERO
      do id = 1,LocalNumAtoms
         nRpts = Potential(id)%n_Rpts
         r_mesh => Potential(id)%Grid%r_mesh(1:nRpts)
         lmax = Potential(id)%lmax
         jmax = Potential(id)%jmax
         jmt  = Potential(id)%jmt
         do ia = 1, Potential(id)%NumSpecies
            do is = 1,n_spin_pola
               pot_l => Potential(id)%potL(:,:,is,ia)
               pot_aver(is) = pot_aver(is) + &
                              getSurfAverage(id,nRpts,lmax,r_mesh,pot_l)
            enddo
            pot_l => Potential(id)%potL_Tilda(:,:,ia)
            pot_r = getSurfAverage(id,nRpts,lmax,r_mesh,pot_l)
            potaver_tilda = potaver_tilda + pot_r*getLocalSpeciesContent(id,ia)
            pot_l => Potential(id)%potL_Coulomb(:,:,ia)
            pot_r = getSurfAverage(id,nRpts,lmax,r_mesh,pot_l)
            potaver_Coulomb = potaver_Coulomb + pot_r*getLocalSpeciesContent(id,ia)
         enddo
         pot_l => Potential(id)%potL_Madelung(:,:)
         pot_r = getSurfAverage(id,nRpts,lmax,r_mesh,pot_l)
         potaver_madelung = potaver_madelung + pot_r
         if (.not.isMTFP) then
            pot_l => Potential(id)%potL_Pseudo(:,:)
            pot_r = getSurfAverage(id,nRpts,lmax,r_mesh,pot_l)
            potaver_hat = potaver_hat + pot_r*getLocalSpeciesContent(id,ia)
         endif
      enddo
!
      write(6,'(/,a,/)')   'Surface average potential before global sum:'
      write(6,'(a,2f13.6)') "             Total    :", pot_aver(1:n_spin_pola)
      write(6,'(a,f13.6)')  '             Tilda    :', potaver_tilda
      write(6,'(a,f13.6)')  '             Madelung :', potaver_madelung
      write(6,'(a,f13.6)')  '             Hat      :', potaver_hat
      write(6,'(a,f13.6)')  '             Coulomb  :', potaver_Coulomb
!
      membuf(1) = pot_aver(1)
      membuf(2) = pot_aver(2)
      membuf(3) = potaver_tilda
      membuf(4) = potaver_madelung
      membuf(5) = potaver_hat
      membuf(6) = potaver_Coulomb
!     ---------------------------------------------------------------
      call GlobalSumInGroup(GroupID,membuf,6)
!     ---------------------------------------------------------------
      pot_aver(1:2) = membuf(1:2)/real(GlobalNumAtoms,kind=RealKind)
      potaver_tilda = membuf(3)/real(GlobalNumAtoms,kind=RealKind)
      potaver_madelung = membuf(4)/real(GlobalNumAtoms,kind=RealKind)
      potaver_hat = membuf(5)/real(GlobalNumAtoms,kind=RealKind)
      potaver_Coulomb = membuf(6)/real(GlobalNumAtoms,kind=RealKind)
!
      write(6,'(/,a,/)')   'Surface average potential after global sum:'
      write(6,'(a,2f13.6)') "             Total    :", pot_aver(1:n_spin_pola)
      write(6,'(a,f13.6)')  '             Tilda    :', potaver_tilda
      write(6,'(a,f13.6)')  '             Madelung :', potaver_madelung
      write(6,'(a,f13.6)')  '             Hat      :', potaver_hat
      write(6,'(a,f13.6)')  '             Coulomb  :', potaver_Coulomb
#endif
!
      t1 = getTime()
      do id = 1,LocalNumAtoms
         nRpts = Potential(id)%n_Rpts
         jmax = Potential(id)%jmax
         lmax = Potential(id)%lmax
         kmax = (lmax+1)*(lmax+1)
         jend = Potential(id)%jend
         vol_ints = getVolume(id)-getInscrSphVolume(id)
         flag_jl => Potential(id)%PotCompFlag(1:jmax)
         if ( isChargeSymmOn ) then
            p_flags => getSymmetryFlags(id)
            flag_jl(1:jmax) = p_flags(1:jmax)
            do ia = 1, Potential(id)%NumSpecies
               do is = 1,n_spin_pola
                  do jl = 1,jmax
                     if ( flag_jl(jl) ==0 ) then
                        Potential(id)%potL(1:nRpts,jl,is,ia) = CZERO
                        Potential(id)%potL_Exch(1:nRpts,jl,is,ia) = CZERO
                        Potential(id)%enL_Exch(1:nRpts,jl,is,ia) = CZERO
#ifdef UNCOMMENT
                     else if ( flag_jl(jl) == 1 ) then
                        do ir = 1,nRpts
                           pot_r = real(Potential(id)%potL(ir,jl,is,ia),kind=RealKind)
                           Potential(id)%potL(ir,jl,is,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                           pot_r = real(Potential(id)%potL_Exch(ir,jl,is,ia),kind=RealKind)
                           Potential(id)%potL_Exch(ir,jl,is,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                           pot_r = real(Potential(id)%enL_Exch(ir,jl,is,ia),kind=RealKind)
                           Potential(id)%enL_Exch(ir,jl,is,ia)  = cmplx(pot_r,ZERO,kind=CmplxKind)
                        enddo
                     else if ( flag_jl(jl) == 2 ) then
                        do ir = 1,nRpts
                           pot_r = real( cfact*Potential(id)%potL(ir,jl,is,ia),kind=RealKind )
                           Potential(id)%potL(ir,jl,is,ia) = cmplx(ZERO,pot_r,kind=CmplxKind)
                           pot_r = real( cfact*Potential(id)%potL_Exch(ir,jl,is,ia),kind=RealKind )
                           Potential(id)%potL_Exch(ir,jl,is,ia) = cmplx(ZERO,pot_r,kind=CmplxKind)
                           pot_r = real( cfact*Potential(id)%enL_Exch(ir,jl,is,ia),kind=RealKind )
                           Potential(id)%enL_Exch(ir,jl,is,ia)  = cmplx(ZERO,pot_r,kind=CmplxKind)
                        enddo
#endif
                     endif
                  enddo
               enddo
               do jl = 1,jmax
                  if ( flag_jl(jl) ==0 ) then
                     Potential(id)%potL_Coulomb(1:nRpts,jl,ia)  = CZERO
                     Potential(id)%potL_Tilda(1:nRpts,jl,ia)    = CZERO
#ifdef UNCOMMENT
                  else if ( flag_jl(jl) == 1 ) then
                     do ir = 1,nRpts
                        pot_r = real(Potential(id)%potL_Coulomb(ir,jl,ia),kind=RealKind)
                        Potential(id)%potL_Coulomb(ir,jl,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                        pot_r = real(Potential(id)%potL_Tilda(ir,jl,ia),kind=RealKind)
                        Potential(id)%potL_Tilda(ir,jl,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                     enddo
                  else if ( flag_jl(jl) == 2 ) then
                     do ir = 1,nRpts
                        pot_r = real( cfact*Potential(id)%potL_Coulomb(ir,jl,ia),kind=RealKind )
                        Potential(id)%potL_Coulomb(ir,jl,ia) = cmplx(ZERO,pot_r,kind=CmplxKind)
                        pot_r = real( cfact*Potential(id)%potL_Tilda(ir,jl,ia),kind=RealKind )
                        Potential(id)%potL_Tilda(ir,jl,ia) = cmplx(ZERO,pot_r,kind=CmplxKind)
                     enddo
#endif
                  endif
               enddo
            enddo
            do jl = 1,jmax
               if ( flag_jl(jl) ==0 ) then
                  Potential(id)%potL_Madelung(1:nRpts,jl) = CZERO
#ifdef UNCOMMENT
               else if ( flag_jl(jl) == 1 ) then
                  do ir = 1,nRpts
                     pot_r = real(Potential(id)%potL_Madelung(ir,jl),kind=RealKind)
                     Potential(id)%potL_Madelung(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                  enddo
               else if ( flag_jl(jl) == 2 ) then
                  do ir = 1,nRpts
                     pot_r = real( cfact*Potential(id)%potL_Madelung(ir,jl),kind=RealKind )
                     Potential(id)%potL_Madelung(ir,jl)  = cmplx(ZERO,pot_r,kind=CmplxKind)
                  enddo
#endif
               endif
            enddo
            if (.not.isMTFP) then
               do jl = 1,jmax
                  if ( flag_jl(jl) ==0 ) then
                     Potential(id)%potL_Pseudo(1:nRpts,jl) = CZERO
#ifdef UNCOMMENT
                  else if ( flag_jl(jl) == 1 ) then
                     do ir = 1,nRpts
                        pot_r = real(Potential(id)%potL_Pseudo(ir,jl),kind=RealKind)
                        Potential(id)%potL_Pseudo(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                     enddo
                  else if ( flag_jl(jl) == 2 ) then
                     do ir = 1,nRpts
                        pot_r = real( cfact*Potential(id)%potL_Pseudo(ir,jl),kind=RealKind )
                        Potential(id)%potL_Pseudo(ir,jl)  = cmplx(ZERO,pot_r,kind=CmplxKind)
                     enddo
#endif
                  endif
               enddo
            endif
         else
            flag_jl(1) = 1
            do ia = 1, Potential(id)%NumSpecies
               do is = 1,n_spin_pola
                  do jl = 2,jmax
                     flag_jl(jl) = 0
                     l = lofj(jl)
                     LOOP_nr0: do ir = 1,nRpts
                        if ( abs(Potential(id)%potL(ir,jl,is,ia)) > ONE/(10.0d0**(l+4)) ) then
!                       if ( abs(Potential(id)%potL(ir,jl,is,ia)) > pot_tol ) then
                           flag_jl(jl) = 1
                           exit LOOP_nr0
                        endif
                     enddo LOOP_nr0
                  enddo
               enddo
            enddo
            do ia = 1, Potential(id)%NumSpecies
               do is = 1,n_spin_pola
                  do jl = 2,jmax
                     if ( flag_jl(jl) == 0 ) then
                        Potential(id)%potL(1:nRpts,jl,is,ia) = CZERO
                        Potential(id)%potL_Exch(1:nRpts,jl,is,ia)  = CZERO
                        Potential(id)%enL_Exch(1:nRpts,jl,is,ia)   = CZERO
                     else if (mofj(jl) == 0) then
                        do ir = 1,nRpts
                           pot_r = real(Potential(id)%potL(ir,jl,is,ia),kind=RealKind)
                           Potential(id)%potL(ir,jl,is,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                           pot_r = real(Potential(id)%potL_Exch(ir,jl,is,ia),kind=RealKind)
                           Potential(id)%potL_Exch(ir,jl,is,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                           pot_r = real(Potential(id)%enL_Exch(ir,jl,is,ia),kind=RealKind)
                           Potential(id)%enL_Exch(ir,jl,is,ia)  = cmplx(pot_r,ZERO,kind=CmplxKind)
                        enddo
                     else
                        do ir = 1,nRpts
                           pot_r = real(Potential(id)%potL(ir,jl,is,ia),kind=RealKind)
                           pot_i = real(cfact*Potential(id)%potL(ir,jl,is,ia),kind=RealKind)
                           if ( abs(pot_r) /= ZERO  .and. abs(pot_i)/abs(pot_r)<TEN2m7 ) then
                              flag_jl(jl) = 1
                              Potential(id)%potL(ir,jl,is,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                              pot_r = real(Potential(id)%potL_Exch(ir,jl,is,ia),kind=RealKind)
                              Potential(id)%potL_Exch(ir,jl,is,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                              pot_r = real(Potential(id)%enL_Exch(ir,jl,is,ia),kind=RealKind)
                              Potential(id)%enL_Exch(ir,jl,is,ia)  = cmplx(pot_r,ZERO,kind=CmplxKind)
                           else if ( abs(pot_i)/=ZERO .and. abs(pot_r)/abs(pot_i)< TEN2m7 ) then
                              Potential(id)%potL(ir,jl,is,ia) = cmplx(ZERO,pot_i,kind=CmplxKind)
                              flag_jl(jl) = 2
                              pot_r = real(cfact*Potential(id)%potL_Exch(ir,jl,is,ia),kind=RealKind )
                              Potential(id)%potL_Exch(ir,jl,is,ia) = cmplx(ZERO,pot_r,kind=CmplxKind)
                              pot_r = real(cfact*Potential(id)%enL_Exch(ir,jl,is,ia),kind=RealKind )
                              Potential(id)%enL_Exch(ir,jl,is,ia)  = cmplx(ZERO,pot_r,kind=CmplxKind)
                           else
                              flag_jl(jl) = 3
                           endif
                        enddo
                     endif
#ifdef UNCOMMENT
                     l = lofj(jl)
                     LOOP_NR5: do ir = 1,nRpts
                        if ( abs(Potential(id)%potL_Exch(ir,jl,is,ia)) > ONE/(10.0d0**(l+4)) ) then
                           Potential(id)%potL_Exch_flag(jl) = 1
                           exit LOOP_NR5
                        endif
                     enddo LOOP_NR5
                     LOOP_NR6: do ir = 1,nRpts
                        if ( abs(Potential(id)%potL_XCHat(ir,jl,is,ia)) > ONE/(10.0d0**(l+4)) ) then
                           Potential(id)%potL_XCHat_flag(jl) = 1
                           exit LOOP_NR6
                        endif
                     enddo LOOP_NR6
#endif
                  enddo
               enddo
               do jl = 2,jmax
                  if ( flag_jl(jl) == 0 ) then
                     Potential(id)%potL_Coulomb(1:nRpts,jl,ia) = CZERO
                     Potential(id)%potL_Tilda(1:nRpts,jl,ia) = CZERO
                  else if ( flag_jl(jl) == 1 ) then
                     do ir = 1,nRpts
                        pot_r = real(Potential(id)%potL_Coulomb(ir,jl,ia),kind=RealKind)
                        Potential(id)%potL_Coulomb(ir,jl,ia) = cmplx(pot_r,ZERO,kind=RealKind)
                        pot_r = real(Potential(id)%potL_Tilda(ir,jl,ia),kind=RealKind)
                        Potential(id)%potL_Tilda(ir,jl,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                     enddo
#ifdef UNCOMMENT
                  else if (mofj(jl) == 0) then
                     do ir = 1,nRpts
                        pot_r = real(Potential(id)%potL_Coulomb(ir,jl,ia),kind=RealKind)
                        Potential(id)%potL_Coulomb(ir,jl,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                        pot_r = real(Potential(id)%potL_Tilda(ir,jl,ia),kind=RealKind)
                        Potential(id)%potL_Tilda(ir,jl,ia) = cmplx(pot_r,ZERO,kind=CmplxKind)
                     enddo
                  else if (flag_jl(jl) == 2) then
                     do ir = 1,nRpts
                        pot_r = real(cfact*Potential(id)%potL_Coulomb(ir,jl,ia),kind=RealKind )
                        Potential(id)%potL_Coulomb(ir,jl,ia) = cmplx(ZERO,pot_r,kind=CmplxKind)
                        pot_r = real(cfact*Potential(id)%potL_Tilda(ir,jl,ia),kind=RealKind )
                        Potential(id)%potL_Tilda(ir,jl,ia) = cmplx(ZERO,pot_r,kind=CmplxKind)
                     enddo
#endif
                  endif
#ifdef UNCOMMENT
                  l = lofj(jl)
                  LOOP_NR4: do ir = 1,nRpts
                     if ( abs(Potential(id)%potL_Coulomb(ir,jl,ia)) > ONE/(10.0d0**(l+4)) ) then
                        Potential(id)%potL_Coulomb_flag(jl) = 1 
                        exit LOOP_NR4
                     endif
                  enddo LOOP_NR4
                  LOOP_NR1: do ir = 1,nRpts
                     if ( abs(Potential(id)%potL_Tilda(ir,jl,ia)) > ONE/(10.0d0**(l+4)) ) then
                        Potential(id)%potL_Tilda_flag(jl) = 1
                        exit LOOP_NR1
                     endif
                  enddo LOOP_NR1
#endif
               enddo
            enddo
            do jl = 2,jmax
               if ( flag_jl(jl) == 0 ) then
                  Potential(id)%potL_Madelung(1:nRpts,jl) = CZERO
               else if ( flag_jl(jl) == 1 ) then
                  do ir = 1,nRpts
                     pot_r = real(Potential(id)%potL_Madelung(ir,jl),kind=RealKind)
                     Potential(id)%potL_Madelung(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                  enddo
#ifdef UNCOMMENT
               else if (mofj(jl) == 0) then
                  do ir = 1,nRpts
                     pot_r = real(Potential(id)%potL_Madelung(ir,jl),kind=RealKind)
                     Potential(id)%potL_Madelung(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                  enddo
#endif
               else if ( flag_jl(jl) == 2 ) then
                  do ir = 1,nRpts
                     pot_r = real( cfact*Potential(id)%potL_Madelung(ir,jl),kind=RealKind )
                     Potential(id)%potL_Madelung(ir,jl)  = cmplx(ZERO,pot_r,kind=CmplxKind)
                  enddo
               endif
#ifdef UNCOMMENT
               l = lofj(jl)
               LOOP_NR2: do ir = 1,nRpts
                  if ( abs(Potential(id)%potL_Madelung(ir,jl)) > ONE/(10.0d0**(L+4)) ) then
                     Potential(id)%potL_Madelung_flag(jl) = 1
                     exit LOOP_NR2
                  endif
               enddo LOOP_NR2
#endif
            enddo
            if (.not.isMTFP) then
               do jl = 2,jmax
                  if ( flag_jl(jl) == 0 ) then
                     Potential(id)%potL_Pseudo(1:nRpts,jl) = CZERO
                  else if ( flag_jl(jl) == 1 ) then
                     do ir = 1,nRpts
                        pot_r = real(Potential(id)%potL_Pseudo(ir,jl),kind=RealKind)
                        Potential(id)%potL_Pseudo(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                     enddo
#ifdef UNCOMMENT
                  else if (mofj(jl) == 0) then
                     do ir = 1,nRpts
                        pot_r = real(Potential(id)%potL_Pseudo(ir,jl),kind=RealKind)
                        Potential(id)%potL_Pseudo(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                     enddo
#endif
                  else if ( flag_jl(jl) == 2 ) then
                     do ir = 1,nRpts
                        pot_r = real( cfact*Potential(id)%potL_Pseudo(ir,jl),kind=RealKind )
                        Potential(id)%potL_Pseudo(ir,jl) = cmplx(ZERO,pot_r,kind=CmplxKind)
                     enddo
                  endif
#ifdef UNCOMMENT
                  l = lofj(jl)
                  LOOP_NR3: do ir = 1,nRpts
                     if ( abs(Potential(id)%potL_Pseudo(ir,jl)) > ONE/(10.0d0**(l+4)) ) then
                        Potential(id)%potL_Pseudo_flag(jl) = 1
                        exit LOOP_NR3
                     endif
                  enddo LOOP_NR3
#endif
               enddo
            endif
         endif
!
         do ia = 1, Potential(id)%NumSpecies
            do is = 1,n_spin_pola
               do jl = 1,jmax
                  Potential(id)%potL(jend+1:nRpts,jl,is,ia) = CZERO
               enddo
               Potential(id)%potr_sph(jend+1:nRpts,is,ia) = ZERO
!
               if (.not.isMTFP) then
                  r_mesh => Potential(id)%Grid%r_mesh(1:nRpts)
                  pot_l => Potential(id)%potL(1:nRpts,1:jmax,is,ia)
!                 ====================================================
!                 Integration of potential over the intestitial region
!                 ====================================================
                  vs = getVolumeIntegration( id, nRpts, r_mesh,         &
                                            kmax, jmax, 0, pot_l, vs_mt)
                  V0_inter(is) = V0_inter(is) + getLocalSpeciesContent(id,ia)*(vs - vs_mt)
               endif
            enddo
         enddo
!
#ifdef POT_DEBUG
         do ia = 1, Potential(id)%NumSpecies
            rfac = getLocalSpeciesContent(id,ia)
            do is = 1,n_spin_pola
               pot_l => Potential(id)%potL(1:nRpts,1:jmax,is,ia)
               vol_int(1) = vol_int(1) + rfac*getVolumeIntegration( id, nRpts, &
                                            r_mesh, kmax, jmax, 0, pot_l, vs_mt )
               pot_l => Potential(id)%potL_Exch(1:nRpts,1:jmax,is,ia)
               vol_int(5) = vol_int(5) + rfac*getVolumeIntegration( id, nRpts, &
                                            r_mesh, kmax, jmax, 0, pot_l, vs_mt )
            enddo
            pot_l => Potential(id)%potL_Tilda(1:nRpts,1:jmax,ia)
            vol_int(2) = vol_int(2) + rfac*getVolumeIntegration( id, nRpts, &
                                            r_mesh, kmax, jmax, 0, pot_l, vs_mt )
            pot_l => Potential(id)%potL_Coulomb(1:nRpts,1:jmax,ia)
            vol_int(6) = vol_int(6) + rfac*getVolumeIntegration( id, nRpts, &
                                            r_mesh, kmax, jmax, 0, pot_l, vs_mt )
         enddo
         pot_l => Potential(id)%potL_Madelung(1:nRpts,1:jmax)
         vol_int(3) = vol_int(3) + rfac*getVolumeIntegration( id, nRpts, &
                                            r_mesh, kmax, jmax, 0, pot_l, vs_mt )
         if (isMTFP) then
            vol_int(4) = ZERO
         else
            pot_l => Potential(id)%potL_Pseudo(1:nRpts,1:jmax)
            vol_int(4) = vol_int(4) + rfac*getVolumeIntegration( id, nRpts, &
                                              r_mesh, kmax, jmax, 0, pot_l, vs_mt )
         endif
#endif
!
         Potential(id)%Madelung_Shift = ZERO
      enddo
      if (MyPE == 0) then
         write(6,'(/,a,f10.5,/)')'Time:: computeNewPotentail:: Loop_id 2: ',getTime()-t1
      endif
!
      if (isMTFP) then
         V0_inter = v0 ! v0(is) is calculated in calZeroPotential_MTFP
      else
         call GlobalSumInGroup(GroupID,V0_inter,n_spin_pola)
         V0_inter = V0_inter/getTotalInterstitialVolume()
      endif
      vtmp = HALF*(V0_inter(1)+V0_inter(n_spin_pola))
      if ( maxval(print_level) >= 0 ) then
         write(6,'(/,a,2x,2d16.8)')"Potential shift V0_inter =", V0_inter
      endif
!
!     ================================================================
!     substract muffin-tin zero from the potential.....
!     ================================================================
      t1 = getTime()
      do id = 1,LocalNumAtoms
         jend = Potential(id)%jend
         r_mesh => Potential(id)%Grid%r_mesh
         do ia = 1, Potential(id)%NumSpecies
            do is = 1,n_spin_pola
               do ir = 1,jend
                  Potential(id)%potL(ir,1,is,ia)= Potential(id)%potL(ir,1,is,ia) - &
                                                  cmplx(V0_inter(is)/Y0,Zero,kind=CmplxKind)
                  pot_r = real(Potential(id)%potL(ir,1,is,ia),kind=RealKind)
                  Potential(id)%potr_sph(ir,is,ia) = Y0*r_mesh(ir)*pot_r
               enddo
            enddo
         enddo
      enddo
      if (MyPE == 0) then
         write(6,'(/,a,f10.5,/)')'Time:: computeNewPotentail:: Loop_id 3: ',getTime()-t1
      endif
      v_shift(1:n_spin_pola) = V0_inter(1:n_spin_pola)
!
#ifdef POT_DEBUG
      call GlobalSumInGroup(GroupID,vol_int,6)
      vol_int(1:6) = vol_int(1:6)/getSystemVolume()
!
      write(6,'(/,a,/)')   'Volume average potential before shift:'
      write(6,'(a,2f13.6)') "             Total    :", vol_int(1)
      write(6,'(a,f13.6)')  '             Tilda    :', vol_int(2)
      write(6,'(a,f13.6)')  '             Madelung :', vol_int(3)
      write(6,'(a,f13.6)')  '             Pseudo   :', vol_int(4)
      write(6,'(a,f13.6)')  '             Exchg    :', vol_int(5)
      write(6,'(a,f13.6,/)')'             Coulomb  :', vol_int(6)
#endif
!
      vdif(1) = V0_inter(n_spin_pola) - V0_inter(1)
      MadelungShiftTable = ZERO
      EmptyTable = .false.
!
!     ================================================================
!     The following 2 lines are added to take care full-potential case
!     ----------------------------------------------------------------
      rhoint_sav = ZERO
      MuffinTinZeroCoulumb = ZERO
      MuffinTinZeroExchCorr = ZERO
!     ================================================================
#ifdef TIMING
      t0 = getTime()-t0
      write(6,'(a,f10.5)')' PotentialGeneration ::  Time:',t0
#endif
   endif
!
   nullify(p_flags, flag_jl)
!
   end subroutine computeNewPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine constructNewPotential()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: j, is, ir, jmt, ia
!
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), pointer :: NewSphPotr(:)
!
   vdif(1) = v0(n_spin_pola)-v0(1)
!
   do j = 1,LocalNumAtoms
      jmt = Potential(j)%jmt
      r_mesh => Potential(j)%Grid%r_mesh(1:jmt)
      do ia = 1, Potential(j)%NumSpecies
         do is = 1,n_spin_pola
            if (abs(v0(is)) > TEN2m8) then
               NewSphPotr => Potential(j)%potr_sph(1:jmt,is,ia)
               do ir = 1, jmt
                  NewSphPotr(ir) = NewSphPotr(ir) - v0(is)*r_mesh(ir)
               enddo
            endif
         enddo
      enddo
   enddo
!
   do is = 1,n_spin_pola
      v0(is) = ZERO
   enddo
!
   end subroutine constructNewPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSphPotential()
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use InterpolationModule, only : FitInterp
   use IntegrationModule, only : calIntegration
!
   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies
   use AtomModule, only : getLocalSpeciesContent
!
   use SystemModule, only : getAtomicNumber, getAlloyElementContent, &
                            getNumAlloyElements
!
   use SystemVolumeModule, only : getAtomicVPVolume
   use SystemVolumeModule, only : getSystemVolume
!
   use PolyhedraModule, only : getVolume
!
   use MadelungModule, only : getMadelungMatrix
!
   use ChargeDistributionModule, only : getGlobalOnSiteElectronTableOld,   &
                                        getGlobalMTSphereElectronTableOld, &
                                        getInterstitialElectronDensityOld, &
                                        getGlobalTableLine
!
   use ChargeDensityModule, only : getSphChargeDensity, getSphMomentDensity
!
   use ExchCorrFunctionalModule, only : calSphExchangeCorrelation
   use ExchCorrFunctionalModule, only : getExchCorrPot
!
   use PotentialTypeModule, only : isMuffintinASAPotential
!
   implicit   none
!
   real (kind=RealKind) :: rhoint
   real (kind=RealKind), pointer :: Q_Table(:)
   real (kind=RealKind), pointer :: Qmt_Table(:)
!
   integer (kind=IntKind) :: na, j, ia, lig
   integer (kind=IntKind) :: ir
   integer (kind=IntKind) :: is
   integer (kind=IntKind) :: jmt, jmt_max
   integer (kind=IntKind) :: n_Rpts
   integer (kind=IntKind), pointer :: global_table_line(:)
!
   real (kind=RealKind) :: dps
   real (kind=RealKind) :: ztotss
   real (kind=RealKind) :: surfamt
   real (kind=RealKind) :: rmt
   real (kind=RealKind) :: qsub, dq, dq_mt
   real (kind=RealKind) :: V2rmt
   real (kind=RealKind) :: sums(2)
   real (kind=RealKind) :: vsum, vmad_corr
   real (kind=RealKind), parameter :: PI8 = PI4*TWO
!
   real (kind=RealKind), pointer :: madmat(:)
   real (kind=RealKind), pointer :: Vexc(:)
   real (kind=RealKind), pointer :: rho0_p(:), mom0_p(:)
   real (kind=RealKind), pointer :: der_rho0_p(:), der_mom0_p(:)
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), allocatable, target :: der_rho_ws(:), der_mom_ws(:)
!
!   real (kind=RealKind), allocatable :: vmt1(:)
!
   rhoint = getInterstitialElectronDensityOld()
   global_table_line => getGlobalTableLine()
   Q_Table => getGlobalOnSiteElectronTableOld()
   Qmt_Table => getGlobalMTSphereElectronTableOld()
!
!  ===================================================================
!  vmt1_i = pi4*rho_0*Rmt_i^2 + 2*sum_j[madmat(j,i)*qsub_j], where:
!
!       qsub_i = dQ_i + rho_0*Omega_i
!         dQ_i = Z_i - (Qmt_i + rho_0*Omega0_i)
!
!       For Muffin-tin-ASA, an extra term is added to vmt1_i:
!
!    -TWO*dq_i/rmt-sum_i[dq_i*(vmt1_i+TWO*(Z_i-Qmt_i-dq_i)/rmt)]/dq_int
!
!         dq_i = Q_i - Qmt_i - rho_0*Omega0_i = dQ_i - Z_i + Q_i
!         dq_int = sum_i[dq_i] = sum_i[Q_i - Z_i]
!  ===================================================================
!
   vsum = ZERO
   jmt_max = 0
   do na=1, LocalNumAtoms
      rmt = Potential(na)%Grid%rmt
      surfamt=PI4*rmt*rmt
      madmat => getMadelungMatrix(na)
      vmt1(na) = ZERO
      do j=1,GlobalNumAtoms
         qsub = ZERO
         do ia = 1, getNumAlloyElements(j)
            lig = global_table_line(j) + ia
            qsub = qsub + getAlloyElementContent(j,ia)*(getAtomicNumber(j,ia)-Q_Table(lig))
         enddo
         qsub = qsub + rhoint*getAtomicVPVolume(j)
         vmt1(na) = vmt1(na) + madmat(j)*qsub
      enddo
      vsum = vsum + TWO*vmt1(na)
      vmt1(na) = TWO*vmt1(na) + rhoint*surfamt
      jmt_max = max(jmt_max, Potential(na)%jmt)
   enddo
   allocate(der_rho_ws(jmt_max), der_mom_ws(jmt_max))
   der_rho_ws = ZERO; der_mom_ws = ZERO
!
!  ===================================================================
!  If vshift_switch_on is false, the unit cell summation of the Madelung
!  shift is not zero. We will calculate the Madelung potential in
!  conventional way without the correction term for minimizing the truncation
!  and/or round-off error.
!  ===================================================================
   if (vshift_switch_on) then
      vsum = vsum/real(GlobalNumAtoms,kind=RealKind)
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,vsum)
!     ----------------------------------------------------------------
      vmad_corr = vsum - TWO*MadelungSum*rhoint*getSystemVolume()/    &
                         real(GlobalNumAtoms,kind=RealKind)
      do na = 1, LocalNumAtoms
         vmt1(na) = vmt1(na) - vmad_corr
      enddo
   endif
!
   if (isMuffintinASAPotential()) then
      sums(1:2) = ZERO
      do na=1, LocalNumAtoms
         rmt = Potential(na)%Grid%rmt
         j = GlobalIndex(na)
         dq = getDQinter(j,rhoint)
         vmt1(na)=vmt1(na)-TWO*dq/rmt
         dq_mt = ZERO
         do ia = 1, getLocalNumSpecies(na)
            lig = global_table_line(j) + ia
            dq_mt = dq_mt +                                           &
                  getLocalSpeciesContent(na,ia)*(getLocalAtomicNumber(na,ia)-Qmt_Table(lig))
         enddo
         sums(1) = sums(1) + dq*(vmt1(na)+TWO*(dq_mt-dq)/rmt)
         sums(2) = sums(2) + dq
      enddo
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,sums,2)
!     ----------------------------------------------------------------
!
      do na=1, LocalNumAtoms
         vmt1(na)=vmt1(na)-sums(1)/sums(2)
      enddo
   endif
!
   nullify( madmat )
!
!  ===================================================================
!  calculate hartree and exchange-correlation potential ...........
!  then subtract muffin tin zero...................................
!  ===================================================================
!
   do na = 1, LocalNumAtoms
      jmt = Potential(na)%jmt
      n_Rpts = Potential(na)%n_Rpts
      r_mesh => Potential(na)%Grid%r_mesh(1:n_Rpts)
      sqrt_r(0) = ZERO
      do ir = 1, jmt
         sqrt_r(ir) = sqrt(r_mesh(ir))
      enddo
!
      do ia = 1, Potential(na)%NumSpecies
         ztotss = real( getLocalAtomicNumber(na,ia),kind=RealKind)
!        ============================================================
!        retrieve the charge density and the moment density
!        ============================================================
         if (gga_functional) then
            rho0_p => getSphChargeDensity("TotalNew",na,ia,der_rho0_p)
         else
            rho0_p => getSphChargeDensity("TotalNew",na,ia)
            der_rho0_p => der_rho_ws
         endif
         
#ifdef No_BLAS
         rho_tmp_r(1:jmt) = rho0_p(1:jmt)
#else
!        ------------------------------------------------------------
         call dcopy(jmt,rho0_p(1:jmt),1,rho_tmp_r(1:jmt),1)
!        ------------------------------------------------------------
#endif
         if (n_spin_pola == 2) then
            if (gga_functional) then
               mom0_p => getSphMomentDensity("TotalNew",na,ia,der_mom0_p)
            else
               mom0_p => getSphMomentDensity("TotalNew",na,ia)
               der_mom0_p => der_mom_ws
            endif
#ifdef No_BLAS
            mom_tmp_r(1:jmt) = mom0_p(1:jmt)
#else
!           ---------------------------------------------------------
            call dcopy(jmt,mom0_p(1:jmt),1,mom_tmp_r(1:jmt),1)
!           ---------------------------------------------------------
#endif
         else
            mom_tmp_r = ZERO
            mom0_p => mom_tmp_r
            der_mom0_p => der_mom_ws
         endif
!
!        ============================================================
!        Coulumb potential: vh0+vh1+vh2-vmt1,   r <= Rmt_i, where:
!
!             vh0(r) = -2Z/r
!             vh1(r) = (2/r)*Integ[4*pi*r1^2*rho(r1)]dr1 ;   r1<r
!             vh2(r) = 2*Integ[4*pi*r1*rho(r1)]dr1       ; r<r1<Rmt_i
!
!             V_1(r) = Integ[r1^2*rho(r1)*dr1];             r1 <= r
!             V_2(r) = Integ[r1*rho(r1)*dr1];               r1 <= r
!
!             vmt1 = pi4*rho_0*Rmt_i^2 + 2*sum_j[madmat(j,i)*qsub_j]
!                  qsub_i = dQ_i + rho_0*Omega_i
!                  dQ_i   = Z_i - (Qmt_i + rho_0*Omega0_i)
!        ------------------------------------------------------------
         call FitInterp(4,sqrt_r(1:4),rho_tmp_r(1:4),ZERO,rho_tmp_r(0),dps)
!        ------------------------------------------------------------
         call calIntegration(jmt+1,sqrt_r(0:jmt),rho_tmp_r(0:jmt),V1_r(0:jmt),5)
!        ------------------------------------------------------------
         call calIntegration(jmt+1,sqrt_r(0:jmt),rho_tmp_r(0:jmt),V2_r(0:jmt),3)
!        ------------------------------------------------------------
         if (abs(rmt-r_mesh(jmt)) > TEN2m8) then
!           ---------------------------------------------------------
            call FitInterp(jmt+1,sqrt_r(0:jmt),V2_r(0:jmt),sqrt(rmt),V2rmt,dps)
!           ---------------------------------------------------------
         else
            V2rmt = V2_r(jmt)
         endif
         if (gga_functional) then
            if (n_spin_pola == 1) then
!              ------------------------------------------------------
               call calSphExchangeCorrelation(jmt,rho0_p,der_rho_den=der_rho0_p)
!              ------------------------------------------------------
            else
!              ------------------------------------------------------
               call calSphExchangeCorrelation(jmt,rho0_p,der_rho0_p,mom0_p,der_mom0_p)
!              ------------------------------------------------------
            endif
         else
            if (n_spin_pola == 1) then
!              ------------------------------------------------------
               call calSphExchangeCorrelation(jmt,rho0_p)
!              ------------------------------------------------------
            else
!              ------------------------------------------------------
               call calSphExchangeCorrelation(jmt,rho0_p,mag_den=mom0_p)
!              ------------------------------------------------------
            endif
         endif
         do is =1, n_spin_pola
            Vexc => getExchCorrPot(jmt,is)
            do ir = 1, jmt
               Potential(na)%potr_sph(ir,is,ia) =                                  &
                     TWO*(-ztotss+PI8*(V1_r(ir)+(V2rmt-V2_r(ir))*r_mesh(ir))) + &
                     (Vexc(ir) - vmt1(na))*r_mesh(ir)
            enddo
            do ir = jmt+1, n_Rpts
               Potential(na)%potr_sph(ir,is,ia) = ZERO
            enddo
         enddo
      enddo
      Potential(na)%Madelung_Shift = -vmt1(na)
   enddo
!
!  deallocate( vmt1 )
   deallocate(der_rho_ws, der_mom_ws)
   nullify( r_mesh, Vexc)
!
   end subroutine calSphPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calZeroPotential_MT()
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use PolyhedraModule, only : getVolume
!
   use SystemModule, only : getAtomicNumber, getNumAlloyElements,     &
                            getAlloyElementContent
!
   use SystemVolumeModule, only : getAtomicVPVolume, getTotalInterstitialVolume
!
   use ScfDataModule, only : isInterstitialElectronPolarized
!
   use MadelungModule, only : getMadelungMatrix
!
   use ChargeDistributionModule, only : getGlobalOnSiteElectronTableOld,   &
                                        getInterstitialElectronDensityOld, &
                                        getInterstitialMomentDensity,      &
                                        getGlobalTableLine
!
   use ExchCorrFunctionalModule, only : calSphExchangeCorrelation
   use ExchCorrFunctionalModule, only : getExchCorrPot
!
   implicit   none
!
   real (kind=RealKind) :: rhoint
   real (kind=RealKind) :: mdenint
   real (kind=RealKind), pointer :: Q_Table(:)
!
   integer (kind=IntKind) :: j, na, is, lig, ia
   integer (kind=IntKind), pointer :: global_table_line(:)
!
   real (kind=RealKind), pointer :: madmat(:)
!
   real (kind=RealKind) :: rmt
   real (kind=RealKind) :: surfamt
   real (kind=RealKind) :: omegmt
   real (kind=RealKind) :: qsub, term3, vxout(2)
   real (kind=RealKind) :: vmt
!
   real (kind=RealKind), parameter :: FIFTH = ONE/FIVE
!
!  ===================================================================
!  calculate muffin-tin zero potential for Muffin-tin case
!  ===================================================================
   rhoint = getInterstitialElectronDensityOld()
   if ( isInterstitialElectronPolarized() ) then
      mdenint = getInterstitialMomentDensity()
   else
      mdenint = ZERO
   endif
   global_table_line => getGlobalTableLine()
   Q_Table => getGlobalOnSiteElectronTableOld()
!
!  ===================================================================
!  vmt = [1/5*rhoint*sum_i(surf_i*tau_i) + sum_i(qsub_i*surf_i)
!         + 2*sum_ij(madmat(i,j)*tau_i*qsub_j)]/Omega_0
!  where:
!       surf_i = pi4*Rmt_i^2
!       qsub_i = dQ_i + rho_0*Omega_i
!         dQ_i = Z_i - (Qmt_i + rho_0*Omega0_i)
!       Omega_0 = sum_i(Omega0_i)
!  ===================================================================
!
   vmt=zero
   do na = 1, LocalNumAtoms
      rmt = Potential(na)%Grid%rmt
      surfamt=PI4*rmt*rmt
      omegmt=surfamt*rmt*THIRD
!
      madmat => getMadelungMatrix(na)
      term3 = ZERO
      do j = 1, GlobalNumAtoms
         qsub = ZERO
         do ia = 1, getNumAlloyElements(j)
            lig = global_table_line(j) + ia
            qsub = qsub + getAlloyElementContent(j,ia)*(getAtomicNumber(j,ia)-Q_Table(lig))
         enddo
         qsub = qsub + rhoint*getAtomicVPVolume(j)
         term3=term3+madmat(j)*qsub
      enddo
!
      j = GlobalIndex(na)
      qsub = ZERO
      do ia = 1, getNumAlloyElements(j)
         lig = global_table_line(j) + ia
         qsub = qsub + getAlloyElementContent(j,ia)*(getAtomicNumber(j,ia)-Q_Table(lig))
      enddo
      qsub = qsub + rhoint*getAtomicVPVolume(j)
      vmt = vmt + (surfamt*(FIFTH*omegmt*rhoint+qsub) + term3*TWO*omegmt)
   enddo
   vmt = vmt/getTotalInterstitialVolume()
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,vmt)
!  -------------------------------------------------------------------
!
   nullify( madmat )
!
   MuffinTinZeroCoulumb = vmt
   if (gga_functional) then
      if (n_spin_pola == 1) then
!        -------------------------------------------------------------
         call calSphExchangeCorrelation(rhoint,der_rho_den=ZERO)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call calSphExchangeCorrelation(rhoint,der_rho_den=ZERO,      &
                                        mag_den=mdenint,der_mag_den=ZERO)
!        -------------------------------------------------------------
      endif
   else
      if (n_spin_pola == 1) then
!        -------------------------------------------------------------
         call calSphExchangeCorrelation(rhoint)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call calSphExchangeCorrelation(rhoint,mag_den=mdenint)
!        -------------------------------------------------------------
      endif
   endif
   do is = 1,n_spin_pola
      vxout(is) = getExchCorrPot(is)
      v0(is) = vmt + vxout(is)
      MuffinTinZeroExchCorr(is) = vxout(is)
#ifdef DEBUG
      if (is==1) then
          write(6,'(1x,a,3(1x,f11.8))')"Muffin-tin Zero potential::", &
               v0(is),vmt,vxout(is)
!
      else
         write(6,'(28x,3(1x,f11.8))') v0(is),vmt,vxout(is)
      endif
#endif
   enddo
   nullify(global_table_line)
!
   end subroutine calZeroPotential_MT
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calZeroPotential_MTFP()
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use PolyhedraModule, only : getVolume
!
   use SystemVolumeModule, only : getSystemVolume, getTotalInterstitialVolume
!
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
!
   use MadelungModule, only : getMadelungMatrix
!
   use ChargeDensityModule, only : getNeutralChargeDensity
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use ExchCorrFunctionalModule, only : calSphExchangeCorrelation
   use ExchCorrFunctionalModule, only : getExchCorrPot
!
   implicit   none
!
   integer (kind=IntKind) :: id, nRpts, jmax

   integer (kind=IntKind) :: j, na, is, ia
!
   real (kind=RealKind), pointer :: madmat(:)
!
   real (kind=RealKind) :: vol_int, volmt_int, vtmp, vmt
   real (kind=RealKind) :: rhoint, mdenint
   real (kind=RealKind) :: vxout(2)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
!
   vol_int = ZERO
   do id = 1,LocalNumAtoms
      nRpts = Potential(id)%n_Rpts
      r_mesh => Potential(id)%Grid%r_mesh(1:nRpts)
      jmax = Potential(id)%jmax
      do ia = 1, getLocalNumSpecies(id)
         pot_l => Potential(id)%potL_Tilda(1:nRpts,1:jmax,ia)
         vtmp = getVolumeIntegration( id, nRpts, r_mesh, jmax, 0, pot_l, volmt_int)
         vol_int = vol_int + getLocalSpeciesContent(id,ia)*(vtmp - volmt_int)
      enddo
      pot_l => Potential(id)%potL_Madelung(1:nRpts,1:jmax)
      vtmp = getVolumeIntegration( id, nRpts, r_mesh, jmax, 0, pot_l, volmt_int)
      vol_int = vol_int + vtmp - volmt_int
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,vol_int)
!  -------------------------------------------------------------------
   vmt = vol_int/getTotalInterstitialVolume()
!
!  ===================================================================
!  calculate muffin-tin zero potential for Muffin-tin FP case:
!      vmt = -int_{mt volume} V_coul d^3r/Omega_interstial
!      v0 = vmt + vxc[rhoint,mdenint]
!  ===================================================================
   rhoint = getNeutralChargeDensity()
   mdenint = ZERO ! Since moment density is not calculated in ChargeDensityModule, 
                  ! for now, we assume the interstitial region is not spin-polarized.
                  ! This can be modified in the future.
!
   MuffinTinZeroCoulumb = vmt
   if (gga_functional) then
      if (n_spin_pola == 1) then
!        -------------------------------------------------------------
         call calSphExchangeCorrelation(rhoint,der_rho_den=ZERO)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call calSphExchangeCorrelation(rhoint,der_rho_den=ZERO,      &
                                        mag_den=mdenint,der_mag_den=ZERO)
!        -------------------------------------------------------------
      endif
   else
      if (n_spin_pola == 1) then
!        -------------------------------------------------------------
         call calSphExchangeCorrelation(rhoint)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call calSphExchangeCorrelation(rhoint,mag_den=mdenint)
!        -------------------------------------------------------------
      endif
   endif
   do is = 1,n_spin_pola
      vxout(is) = getExchCorrPot(is)
      v0(is) = vmt + vxout(is)
      MuffinTinZeroExchCorr(is) = vxout(is)
   enddo
!
   end subroutine calZeroPotential_MTFP
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calZeroPotential_ASA()
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use InterpolationModule, only : FitInterp
!
   use ScfDataModule, only : isInterstitialElectronPolarized
!
   use AtomModule, only : getLocalAtomicNumber, getLocalSpeciesContent
!
   use ChargeDistributionModule, only : getGlobalMTSphereElectronTableOld
   use ChargeDistributionModule, only : getGlobalTableLine
!
   use ChargeDensityModule, only : getSphChargeDensity, getSphMomentDensity
!
   use ExchCorrFunctionalModule, only : calSphExchangeCorrelation
   use ExchCorrFunctionalModule, only : getExchCorrPot
!
   implicit   none
!
   real (kind=RealKind), pointer :: Qmt_Table(:)
!
   integer (kind=IntKind) :: na, is, ir, n_Rpts, jmt, j, ia, lig
   integer (kind=IntKind), pointer :: global_table_line(:)
!
   real (kind=RealKind) :: rmt, dq_mt, rhoint_tmp, mdenint_tmp
   real (kind=RealKind) :: dps, vxout(2), sums(3)
   real (kind=RealKind) :: vmt, rho_rmt, mden_rmt
   real (kind=RealKind), pointer :: rho0(:), mom0(:)
   real (kind=RealKind), pointer :: r_mesh(:)
!
!  ===================================================================
!  calculate muffin-tin zero potential for ASA case
!  ===================================================================
!
!  ===================================================================
!  vmt = 2*sum_i(dQmt_i)/N
!  where:
!       dQmt_i = Qmt_i - Z_i
!  ===================================================================
   global_table_line => getGlobalTableLine()
   Qmt_Table => getGlobalMTSphereElectronTableOld()
   vmt=ZERO
   do na = 1, LocalNumAtoms
      j     = GlobalIndex(na)
      rmt   = Potential(na)%Grid%rmt
      dq_mt = ZERO
      do ia = 1, Potential(na)%NumSpecies
         lig = global_table_line(j) + ia
         dq_mt = dq_mt +                                              &
           getLocalSpeciesContent(na,ia)*(Qmt_Table(lig)-getLocalAtomicNumber(na,ia))
      enddo
      vmt = vmt + dq_mt/rmt
   enddo
   vmt = TWO*vmt/real(GlobalNumAtoms,kind=RealKind)
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,vmt)
!  -------------------------------------------------------------------
!
   sums(1:3) = ZERO
   do na = 1, LocalNumAtoms
      n_Rpts = Potential(na)%n_Rpts
      jmt = Potential(na)%Grid%jmt
      rmt = Potential(na)%Grid%rmt
      r_mesh => Potential(na)%Grid%r_mesh(1:n_Rpts)
!
      rhoint_tmp = ZERO
      mdenint_tmp = ZERO
      do ia = 1, Potential(na)%NumSpecies
         rho0 => getSphChargeDensity("TotalNew",na,ia)
         if (abs(rmt-r_mesh(jmt)) > TEN2m8) then
            do ir = 1, n_Rpts
               sqrt_r(ir) = sqrt(r_mesh(ir))
            enddo
!           ----------------------------------------------------------
            call FitInterp( n_Rpts, sqrt_r(1:n_Rpts), rho0(1:n_Rpts),     &
                            sqrt(rmt), rho_rmt, dps )
!           ----------------------------------------------------------
         else
            rho_rmt = rho0(jmt)
         endif
         rhoint_tmp = rhoint_tmp + getLocalSpeciesContent(na,ia)*rho_rmt
!
         if (n_spin_pola == 2 .and. isInterstitialElectronPolarized()) then
            mom0 => getSphMomentDensity("TotalNew",na,ia)
            if (abs(rmt-r_mesh(jmt)) > TEN2m8) then
!              -------------------------------------------------------
               call FitInterp( n_Rpts, sqrt_r(1:n_Rpts), mom0(1:n_Rpts),  &
                               sqrt(rmt), mden_rmt, dps )
!              -------------------------------------------------------
            else
               mden_rmt = mom0(jmt)
            endif
         else
            mden_rmt = ZERO
         endif
         mdenint_tmp = mdenint_tmp + getLocalSpeciesContent(na,ia)*mden_rmt
      enddo
!
      if (gga_functional) then
         if (n_spin_pola == 1) then
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(rhoint_tmp,der_rho_den=ZERO)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(rhoint_tmp,der_rho_den=ZERO, &
                                           mag_den=mdenint_tmp,der_mag_den=ZERO)
!           ----------------------------------------------------------
         endif
      else
         if (n_spin_pola == 1) then
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(rhoint_tmp)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(rhoint_tmp,mag_den=mdenint_tmp)
!           ----------------------------------------------------------
         endif
      endif
      do is = 1,n_spin_pola
         sums(is) = sums(is) + getExchCorrPot(is)
      enddo
      sums(3) = sums(3) + rhoint_tmp*getLocalSpeciesContent(na,ia)
   enddo
!
   sums(1:3) = sums(1:3)/real(GlobalNumAtoms,kind=RealKind)
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,sums,3)
!  -------------------------------------------------------------------
   do is=1,n_spin_pola
      vxout(is) = sums(is)
   enddo
   rhoint_tmp = sums(3)
!
   MuffinTinZeroCoulumb = vmt
   do is = 1,n_spin_pola
      v0(is) = vmt + vxout(is)
      MuffinTinZeroExchCorr(is) = vxout(is)
   enddo
!
!  rhoint_sav = rhoint_tmp
!
   nullify(global_table_line)
   nullify( rho0, mom0)
!
   end subroutine calZeroPotential_ASA
!  ===================================================================
!
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calZeroPotential_MTASA()
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use InterpolationModule, only : FitInterp
!
   use PolyhedraModule, only : getVolume, getInscrSphRadius
!
   use SystemModule, only : getAtomicNumber, getNumAlloyElements, &
                            getAlloyElementContent
!
   use SystemVolumeModule, only : getAtomicVPVolume, getTotalInterstitialVolume
!
   use MadelungModule, only : getMadelungMatrix
!
   use ChargeDistributionModule, only : getGlobalOnSiteElectronTableOld, &
                                        getInterstitialElectronDensityOld
   use ChargeDistributionModule, only : getGlobalTableLine
!
   use ChargeDensityModule, only : getSphChargeDensity, getSphMomentDensity
!
   use ExchCorrFunctionalModule, only : calSphExchangeCorrelation
   use ExchCorrFunctionalModule, only : getExchCorrPot
!
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent, &
                          getLocalAtomicNumber
!
   implicit   none
!
   real (kind=RealKind) :: rhoint
   real (kind=RealKind), pointer :: Q_Table(:)
!
   integer (kind=IntKind) :: j, na, is, ir, jmt, n_Rpts, ia, lig
   integer (kind=IntKind), pointer :: global_table_line(:)
!
   real (kind=RealKind), pointer :: madmat(:)
   real (kind=RealKind), pointer :: rho0(:)
   real (kind=RealKind), pointer :: mom0(:)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind) :: rmt, sqrmt, width
   real (kind=RealKind) :: surfamt
   real (kind=RealKind) :: omegmt
   real (kind=RealKind) :: qsub, term3, dps, vxout(2), dq
   real (kind=RealKind) :: rhoint_tmp, mdenint_tmp
   real (kind=RealKind) :: vmt, rho_rmt, mden_rmt
   real (kind=RealKind) :: sums(3)
!
   real (kind=RealKind), parameter :: FIFTH = ONE/FIVE
!
!  ===================================================================
!  calculate muffin-tin zero potential for Muffin-tin case
!
!  vmt = [1/5*rhoint*sum_i(surf_i*tau_i) + sum_i(qsub_i*surf_i)
!         + 2*sum_ij(madmat(i,j)*tau_i*qsub_j)]/Omega_0
!  where:
!       surf_i = pi4*Rmt_i^2
!       qsub_i = dQ_i + rho_0*Omega_i
!         dQ_i = Z_i - (Qmt_i + rho_0*Omega0_i)
!       Omega_0 = sum_i(Omega0_i)
!  ===================================================================
!
   rhoint = getInterstitialElectronDensityOld()
   global_table_line => getGlobalTableLine()
   Q_Table => getGlobalOnSiteElectronTableOld()
   vmt = zero
   do na = 1, LocalNumAtoms
      rmt = Potential(na)%Grid%rmt
      surfamt = PI4*rmt*rmt
      omegmt  = surfamt*rmt*THIRD
!
      madmat => getMadelungMatrix(na)
      term3 = ZERO
      do j = 1, GlobalNumAtoms
         qsub = ZERO
         do ia = 1, getNumAlloyElements(j)
            lig = global_table_line(j) + ia
            qsub = qsub + getAlloyElementContent(j,ia)*(getAtomicNumber(j,ia)-Q_Table(lig))
         enddo
         qsub = qsub + rhoint*getAtomicVPVolume(j)
         term3 = term3 + madmat(j)*qsub
      enddo
!
      j = GlobalIndex(na)
      qsub = ZERO
      do ia = 1, getLocalNumSpecies(na)
         lig = global_table_line(j) + ia
         qsub = qsub + getLocalSpeciesContent(na,ia)*(getLocalAtomicNumber(na,ia)-Q_Table(lig))
      enddo
      qsub = qsub + rhoint*getAtomicVPVolume(j)
      vmt = vmt + (surfamt*(FIFTH*omegmt*rhoint+qsub) + term3*TWO*omegmt)
   enddo
   vmt = vmt/getTotalInterstitialVolume()
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,vmt)
!  -------------------------------------------------------------------
!
   nullify( madmat )
!
   sums(1:3) = ZERO
   do na = 1, LocalNumAtoms
      j   = GlobalIndex(na)
      jmt = Potential(na)%Grid%jmt
      rmt = Potential(na)%Grid%rmt
      n_Rpts =  Potential(na)%n_Rpts
      r_mesh => Potential(na)%Grid%r_mesh(1:n_Rpts)
!
      do ir = 1, n_Rpts
         sqrt_r(ir) = sqrt(r_mesh(ir))
      enddo
!
      width = sqrt_r(jmt)-sqrt_r(jmt-1)
      sqrmt = sqrt(getInscrSphRadius(na))
!
      rhoint_tmp = ZERO
      mdenint_tmp = ZERO
      do ia = 1, Potential(na)%NumSpecies
         do is = 1, n_spin_pola
             do ir = 1, jmt
               Potential(na)%potr_sph(ir,is,ia) = Potential(na)%potr_sph(ir,is,ia) &
                                            +vmt/(ONE+exp((sqrmt-sqrt_r(ir))/width))
            enddo
         enddo
!
         rho0 => getSphChargeDensity("TotalNew",na,ia)
         if (abs(rmt-r_mesh(jmt)) > TEN2m8) then
!           ----------------------------------------------------------
            call FitInterp(n_Rpts,sqrt_r(1:n_Rpts),rho0(1:n_Rpts),sqrt(rmt),    &
                           rho_rmt,dps)
!           ----------------------------------------------------------
         else
            rho_rmt = rho0(jmt)
         endif
         rhoint_tmp = rhoint_tmp + getLocalSpeciesContent(na,ia)*rho_rmt
!
         if (n_spin_pola == 2) then
            mom0 => getSphMomentDensity("TotalNew",na,ia)
            if (abs(rmt-r_mesh(jmt)) > TEN2m8) then
!              -------------------------------------------------------
               call FitInterp(n_Rpts,sqrt_r(1:n_Rpts),mom0(1:n_Rpts),sqrt(rmt), &
                              mden_rmt,dps)
!              -------------------------------------------------------
            else
               mden_rmt = mom0(jmt)
            endif
         else
            mden_rmt = ZERO
         endif
         mdenint_tmp = mdenint_tmp + getLocalSpeciesContent(na,ia)*mden_rmt
      enddo
!
      dq = getDQinter(j,rhoint_tmp)
      sums(3) = sums(3) + dq
      if (gga_functional) then
         if (n_spin_pola == 1) then
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(rhoint_tmp,der_rho_den=ZERO)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(rhoint_tmp,der_rho_den=ZERO, &
                                           mag_den=mdenint_tmp,der_mag_den=ZERO)
!           ----------------------------------------------------------
         endif
      else
         if (n_spin_pola == 1) then
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(rhoint_tmp)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(rhoint_tmp,mag_den=mdenint_tmp)
!           ----------------------------------------------------------
         endif
      endif
      do is = 1,n_spin_pola
         sums(is) = sums(is) + dq*getExchCorrPot(is)
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,sums,3)
!  -------------------------------------------------------------------
   do is = 1,n_spin_pola
      vxout(is) = sums(is)/sums(3)
   enddo
!
   MuffinTinZeroCoulumb = vmt
   do is = 1,n_spin_pola
      v0(is) = vmt + vxout(is)
      MuffinTinZeroExchCorr(is) = vxout(is)
   enddo
!
   rhoint_sav = rhoint_tmp
!
   nullify(global_table_line)
   nullify( rho0, mom0)
!
   end subroutine calZeroPotential_MTASA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printMadelungShiftTable(iter,fu)
!  ===================================================================
   use ChargeDistributionModule, only : getGlobalMTSphereElectronTableOld, &
                                        getGlobalOnSiteElectronTableOld, &
                                        getGlobalTableLine
!
   use SystemModule, only : getAtomicNumber, getAtomName, getAtomEnergy
   use SystemModule, only : getNumAlloyElements
!
   implicit none
!
   logical :: FileExist
!
   integer (kind=IntKind), intent(in) :: iter
   integer (kind=IntKind), intent(in) :: fu
   integer (kind=IntKind) :: ig, ia, na, lig
   integer (kind=IntKind), pointer :: global_table_line(:)
!
   real (kind=RealKind), pointer :: Qmt_Table(:)
   real (kind=RealKind), pointer :: Q_Table(:)
   real (kind=RealKind), pointer :: atom_en(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('printMadelungShiftTable',                    &
                        'Need to initialize PotentialGenerationModule first')
   else if (EmptyTable) then
      call ErrorHandler('printMadelungShiftTable','Table is empty')
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,20x,a)')'***************************************'
   write(6,'( 20x,a )')'* Output from printMadelungShiftTable *'
   write(6,'(20x,a,/)')'***************************************'
   write(6,'(80(''=''))')
   write(6,'(2x,''Interstitial charge density'',t40,''='',f18.11)')rhoint_sav
!
   if (fu /= 6) then
      inquire(file='MadelungShiftTable',exist=FileExist)
      if (FileExist) then
         open(unit=fu,file='MadelungShiftTable',form='formatted',     &
              status='old',position='append')
      else
         open(unit=fu,file='MadelungShiftTable',form='formatted',     &
              status='unknown')
         write(fu,'(/,80(''-''))')
         write(fu,'(/,20x,a)')'***************************************'
         write(fu,'( 20x,a )')'* Output from printMadelungShiftTable *'
         write(fu,'(20x,a,/)')'***************************************'
      endif
   endif
!
   write(fu,'(/,a,i5)')'# ITERATION :',iter
   if (fu == 6) then
      write(fu,'(80(''=''))')
   else
      write(fu,'(60(''=''))')
   endif
!
   na = 1
   LOOP_ig: do ig = 1, GlobalNumAtoms
      na = getNumAlloyElements(ig)
      if (na > 1) then
         exit LOOP_ig
      endif
   enddo LOOP_ig
!
   if (na == 1) then
      write(fu,'(a)')' Atom   Index      Q          Qmt        dQ        Vmad        Local Energy'
   else
      write(fu,'(a)')' Atom  Site  Species    Q          Qmt        dQ        Vmad        Local Energy'
   endif
   global_table_line => getGlobalTableLine()
   Q_Table => getGlobalOnSiteElectronTableOld()
   Qmt_Table => getGlobalMTSphereElectronTableOld()
   atom_en => getAtomEnergy()
   do ig = 1, GlobalNumAtoms
      do ia = 1, getNumAlloyElements(ig)
!        write(fu,'(2x,a3,2x,i6,4(2x,f20.16))')getAtomName(ig,ia), ig, ia,     &
         if (na == 1) then
            write(fu,'(2x,a3,2x,i6,4(2x,f9.5),6x,f12.5)')getAtomName(ig), ig,  &
                                                 Q_Table(ig),Qmt_Table(ig),    &
                                                 Q_Table(ig)-getAtomicNumber(ig), &
                                                 MadelungShiftTable(ig),atom_en(1,ig)
         else
            lig = global_table_line(ig) + ia
            write(fu,'(2x,a3,2x,i6,2x,i2,4(2x,f9.5),6x,f12.5)')getAtomName(ig,ia), ig, ia,  &
                                                 Q_Table(lig),Qmt_Table(lig),         &
                                                 Q_Table(lig)-getAtomicNumber(ig,ia), &
                                                 MadelungShiftTable(ig),atom_en(1,ig)
         endif
      enddo
   enddo
   if (fu == 6) then
      write(fu,'(80(''=''),/)')
   else
      close(unit=fu)
   endif
   nullify(global_table_line)
!
   end subroutine printMadelungShiftTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPotentialGeneration()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: na
!
   real (kind=RealKind) :: rmt, omegmt
!
   if (.not.Initialized) then
      call ErrorHandler('printPotentialGeneration',                   &
                        'Need to initialize PotentialGenerationModule first')
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,20x,a)')'****************************************'
   write(6,'( 20x,a )')'* Output from printPotentialGeneration *'
   write(6,'(20x,a,/)')'****************************************'
   write(6,'(80(''=''))')
   write(6,'(2x,''Interstitial charge density'',t40,''='',f18.11)')rhoint_sav
!
   do na = 1, LocalNumAtoms
      if (Print_Level(na) >= 0) then
         write(6,'(2x,''For local atom'',t40,'':'',i5)')na
         rmt = Potential(na)%Grid%rmt
         omegmt = PI4*THIRD*rmt*rmt*rmt
         write(6,'(2x,''MT sphere volume'',t40,''='',f18.11)')omegmt
         write(6,'(2x,a,t40,''='',f18.11)')                            &
                 'Muffin-tin zero potential from Coul',MuffinTinZeroCoulumb
         if (n_spin_pola == 1) then
            write(6,'(2x,a,t40,''='',f18.11)')                         &
                 'Muffin-tin zero potential from ExC',MuffinTinZeroExchCorr(1)
            write(6,'(2x,a,t40,''='',f18.11)')                         &
                 'Muffin-tin zero potential',                          &
                 MuffinTinZeroCoulumb + MuffinTinZeroExchCorr(1)
         else
            write(6,'(2x,a,t40,''='',2f18.11)')                        &
                 'Muffin-tin zero potential from ExC',                 &
                 MuffinTinZeroExchCorr(1), MuffinTinZeroExchCorr(2)
            write(6,'(2x,a,t40,''='',f18.11)')                         &
                 'Muffin-tin zero for spin up',                        &
                 MuffinTinZeroCoulumb + MuffinTinZeroExchCorr(1)
            write(6,'(2x,a,t40,''='',f18.11)')                         &
                 'Muffin-tin zero for spin down',                      &
                 MuffinTinZeroCoulumb + MuffinTinZeroExchCorr(2)
         endif
!        write(6,'(10x,''Madelung constant shift'',t40,''='',f18.11)') &
!              alpha_mad(na)*rhoint_sav
!        write(6,'(10x,''Madelung site shift'',t40,''='',f18.11)')     &
!              -alpha_mad(na)*rhoint_sav+Potential(na)%Madelung_Shift
         write(6,'(10x,''Madelung site shift'',t40,''='',f18.11)')     &
               Potential(na)%Madelung_Shift
         write(6,'(80(''=''))')
      endif
   enddo
!
   end subroutine printPotentialGeneration
!  ===================================================================
!
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDQinter(ig,rhoint) result(dq)
!  ===================================================================
   use SystemVolumeModule, only : getAtomicMTVolume, getAtomicVPVolume
!
   use SystemModule, only : getNumAlloyElements, getAlloyElementContent
!
   use PotentialTypeModule, only : isASAPotential
!
   use ChargeDistributionModule, only : getGlobalMTSphereElectronTableOld, &
                                        getGlobalVPCellElectronTableOld, &
                                        getGlobalTableLine
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ig
   integer (kind=IntKind) :: ia, lig
   integer (kind=IntKind), pointer :: global_table_line(:)
!
   real (kind=RealKind), intent(in) :: rhoint
   real (kind=RealKind), pointer :: Qmt_Table(:)
   real (kind=RealKind), pointer :: Qvp_Table(:)
   real (kind=RealKind) :: dq
!
!  *******************************************************************
!  dq is called qint1 in old code, but there is a slight difference:
!  in old code, zsemss+zcorss is used instead of getCoreVPCharg
!
!     dq_i = Q_i - Qmt_i - rho_0*Omega0_i = dQ_i - Z_i + Q_i
!  *******************************************************************
   global_table_line => getGlobalTableLine()
   Qmt_Table => getGlobalMTSphereElectronTableOld()
   Qvp_Table => getGlobalVPCellElectronTableOld()
   if (isASAPotential()) then
      dq = ZERO
      do ia = 1, getNumAlloyElements(ig)
         lig = global_table_line(ig) + ia
         dq = dq + getAlloyElementContent(ig,ia)*(Qvp_Table(lig)-Qmt_Table(lig))
      enddo
   else
      dq = ZERO
      do ia = 1, getNumAlloyElements(ig)
         lig = global_table_line(ig) + ia
         dq = dq + getAlloyElementContent(ig,ia)*(Qvp_Table(lig)-Qmt_Table(lig))
      enddo
      dq = dq -rhoint*(getAtomicVPVolume(ig)-getAtomicMTVolume(ig))
   endif
   nullify(global_table_line)
!
   end function getDQinter
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntraPot()
!  ===================================================================
!
   use MathParamModule, only : CZERO, ZERO, PI4, TWO, Y0inv
   use ErrorHandlerModule, only : ErrorHandler


   use PotentialModule, only: getPotential, getPotLmax
   use ChargeDensityModule, only: getChargeDensity,          &
                                  getMomentDensity,          &
                                  getRhoLmax,                &
                                  getPseudoNumRPts
   use AtomModule, only : getMaxLmax, getLocalAtomicNumber
   use InterpolationModule, only : FitInterp
   use IntegrationModule, only : calIntegration
!
   implicit   none
!
   integer (kind=IntKind) :: ir, id, jl, l, m, lmax, jmax, jmin, ia
   integer (kind=IntKind) :: lmax_rho, jmax_rho
   integer (kind=IntKind) :: nRpts, nRpts_ps, flag_jlrho
!
   real (kind=RealKind) :: Rc, Zi, dps
   real (kind=RealKind) :: V2R0_r, V2R0_i
   real (kind=RealKind), pointer :: r_mesh(:)
#ifdef TIMING
   real (kind=RealKind) :: t0
#endif
!
   complex (kind=CmplxKind), pointer :: rhol(:), potl(:)
!
   type (GridStruct), pointer :: Grid
!
#ifdef TIMING
   t0 = getTime()
#endif
   do id = 1,LocalNumAtoms
!
      lmax = Potential(id)%lmax
      jmax = Potential(id)%jmax
      lmax_rho = getRhoLmax(id)
      jmax_rho = ((lmax_rho+1)*(lmax_rho+2))/2
      jmin = min(jmax,jmax_rho)
!
      Grid => Potential(id)%Grid
      r_mesh => Grid%r_mesh
      nRpts = Potential(id)%jend
      nRpts_ps = getPseudoNumRPts(id)
      Rc = r_mesh(nRpts)
      sqrt_r(0) = ZERO
      do ir = 1, nRpts
         sqrt_r(ir) = sqrt(r_mesh(ir))
      enddo
!
      Potential(id)%potL_Tilda = CZERO  ! set ZERO here
      do ia = 1, Potential(id)%NumSpecies
         Zi = getLocalAtomicNumber(id,ia)
         LoopJl:  do jl = 1,jmin
            l = lofj(jl)
            m = mofj(jl)
            potl => Potential(id)%potL_Tilda(1:nRpts,jl,ia)
!           potl = CZERO
            rhol => getChargeDensity("Tilda",id, ia, jl)
            if ( l == 0 ) then
               do ir = 1,nRpts
                  rho_tmp_r(ir) = real(rhol(ir),kind=RealKind)
               enddo
               rho_tmp_i(0:nRpts) = ZERO
            else if ( m==0 ) then
               do ir = 1,nRpts
                  rho_tmp_r(ir) = real(rhol(ir),kind=RealKind)
               enddo
               rho_tmp_i = ZERO
            else
               do ir = 1,nRpts
                  rho_tmp_r(ir) = real(rhol(ir),kind=RealKind)
                  rho_tmp_i(ir) = real(-sqrtm1*rhol(ir),kind=RealKind)
               enddo
            endif
!
            flag_jlrho=0
            LoopRho: do ir = nRpts,1,-1
               if ( abs(rho_tmp_r(ir)) > TEN2m10 .or. &
                    abs(rho_tmp_i(ir)) > TEN2m10 ) then
                  flag_jlrho = 1
                  exit LoopRho
               endif
            enddo LoopRho
!
            if ( flag_jlrho == 0 ) then
               cycle LoopJl
            endif
!           ----------------------------------------------------------
            call FitInterp( 4, sqrt_r(1:4), rho_tmp_r(1:4), ZERO,         &
                           rho_tmp_r(0), dps )
            call FitInterp( 4, sqrt_r(1:4), rho_tmp_i(1:4), ZERO,         &
                           rho_tmp_i(0), dps )
!           ----------------------------------------------------------
            call calIntegration( nRpts+1, sqrt_r(0:nRpts), rho_tmp_r(0:nRpts),    &
                                 V1_r(0:nRpts), 5+2*l)
            call calIntegration( nRpts+1, sqrt_r(0:nRpts), rho_tmp_i(0:nRpts),    &
                                 V1_i(0:nRpts), 5+2*l)
!           ----------------------------------------------------------
            if ( l <= 1 ) then
               call calIntegration( nRpts_ps+1, sqrt_r(0:nRpts_ps),       &
                                    rho_tmp_r(0:nRpts_ps), V2_r(0:nRpts_ps), 3-2*l)
               call calIntegration( nRpts_ps+1, sqrt_r(0:nRpts_ps),       &
                                    rho_tmp_i(0:nRpts_ps), V2_i(0:nRpts_ps), 3-2*l)
            else
               do ir = 1,nRpts_ps
                  rho_tmp_r(ir) = rho_tmp_r(ir)/(r_mesh(ir)**(l-1))
                  rho_tmp_i(ir) = rho_tmp_i(ir)/(r_mesh(ir)**(l-1))
               enddo
               call calIntegration( nRpts_ps+1, sqrt_r(0:nRpts_ps),       &
                                    rho_tmp_r(0:nRpts_ps), V2_r(0:nRpts_ps), 1 )
               call calIntegration( nRpts+1, sqrt_r(0:nRpts_ps),          &
                                    rho_tmp_i(0:nRpts_ps), V2_i(0:nRpts_ps), 1 )
            endif
!
            V2R0_r = V2_r(nRpts_ps)
            V2R0_i = V2_i(nRpts_ps)
!
            if ( l == 0 ) then
               do ir = 1, nRpts_ps
                  potl(ir) = cmplx( (FOUR*PI4*V1_r(ir)-TWO*Zi*Y0inv)/r_mesh(ir) +  &
                                    FOUR*PI4*(V2R0_r-V2_r(ir)), ZERO,kind=CmplxKind )
               enddo
               do ir = nRpts_ps+1, nRpts
                  potl(ir) = cmplx( (FOUR*PI4*V1_r(ir)-TWO*Zi*Y0inv)/r_mesh(ir),ZERO,kind=CmplxKind )
               enddo
!              =======================================================
!              Note: VcoulombR0 will be added to the total energy
!                    calculation in the rho*v_coul term.
!              =======================================================
               Potential(id)%VcoulombR0(1) = FOUR*PI4*V2R0_r*Y0
            else if ( m == 0 ) then
               do ir = 1, nRpts_ps
                  potl(ir) = cmplx( (FOUR*PI4/(2*l+1))*                     &
                                    ( V1_r(ir)/(r_mesh(ir)**(l+1)) +        &
                                      (V2R0_r-V2_r(ir))*(r_mesh(ir)**l) ),  &
                                    ZERO,kind=CmplxKind)
               enddo
               do ir = nRpts_ps+1, nRpts
                  potl(ir) = cmplx( (FOUR*PI4/(2*l+1))*                     &
                                    ( V1_r(ir)/(r_mesh(ir)**(l+1))),        &
                                    ZERO,kind=CmplxKind)
               enddo
            else
               do ir = 1,nRpts_ps
                  potl(ir) = cmplx( (FOUR*PI4/(2*l+1))*                     &
                                    ( V1_r(ir)/(r_mesh(ir)**(l+1)) +        &
                                      (V2R0_r-V2_r(ir))*(r_mesh(ir)**l)),   &
                                    (FOUR*PI4/(2*l+1))*                     &
                                    ( V1_i(ir)/(r_mesh(ir)**(l+1)) +        &
                                      (V2R0_i-V2_i(ir))*(r_mesh(ir)**l)),   &
                                    kind=CmplxKind)
               enddo
               do ir = nRpts_ps+1,nRpts
                  potl(ir) = cmplx( (FOUR*PI4/(2*l+1))*                  &
                                    ( V1_r(ir)/(r_mesh(ir)**(l+1)) )  ,  &
                                    (FOUR*PI4/(2*l+1))*                  &
                                    ( V1_i(ir)/(r_mesh(ir)**(l+1)) )  ,  &
                                    kind=CmplxKind)
               enddo
            endif
!
         enddo LoopJl
      enddo
   enddo
!
#ifdef TIMING
      t0 = getTime()-t0
      write(6,'(a,f10.5)')' calIntraPot :: Time ',t0
#endif
!
   end subroutine calIntraPot
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calInterPlusMadPot()
!  ===================================================================
   use AtomModule, only : getLocalAtomicNumber
   use SystemModule, only : getAtomicNumber, getAdditionalElectrons,  &
                            getNumAlloyElements, getAlloyElementContent
   use SystemVolumeModule, only : getSystemVolume
   use MadelungModule, only : getDLMatrix, getDLFactor, getMadelungMatrix
   use ChargeDensityModule, only : getMultipoleMoment,                &
                                   isChargeComponentZero,             &
                                   getNeutralChargeDensity,           &
                                   getPseudoNumRPts
   use ChargeDistributionModule, only : getInterstitialElectronDensity
   use PolyhedraModule, only : getVolume, getInscrSphVolume, getInscrSphRadius
   use StepFunctionModule, only : getVolumeIntegration
   use InterpolationModule, only : FitInterp
!
   implicit none
!
   integer (kind=IntKind) :: jl_pot, jmax_pot, kl_pot, m_pot, l_pot
   integer (kind=IntKind) :: kl_rho, kmax_rho, m_rho, l_rho, jl_rho
   integer (kind=IntKind) :: mdif, lsum,  jmax_rho
   integer (kind=IntKind) :: local_id, ig, ia, kl, ir, n_RPts, jmt, lig
   integer (kind=IntKind), pointer :: table_line(:)
!
   real (kind=RealKind)   :: Zi, Zj, a0, rho_neutral, rho_nj, alpha_mad
   real (kind=RealKind), target :: dlf_0(1:1)
   real (kind=RealKind), pointer :: dlf(:), r_mesh(:), mad(:)
#ifdef TIMING
   real (kind=RealKind) :: t0
#endif
!
   complex (kind=CmplxKind) :: sumjl, sumat, ql_av
   complex (kind=CmplxKind), allocatable, target :: dlm_0(:,:)
   complex (kind=CmplxKind), pointer :: dlm(:,:)
   complex (kind=CmplxKind), pointer :: v_tilt(:,:), ql(:)
!
#ifdef TIMING
   t0 = getTime()
#endif
   allocate(dlm_0(1:GlobalNumAtoms,1))
!
   rho_neutral = (getNeutralChargeDensity()-                           &
                  getAdditionalElectrons()/getSystemVolume())/Y0
   rho_nj = rho_neutral*getSystemVolume()/(GlobalNumAtoms*Y0)
!
   kmax_rho = (lmax_rho_max+1)**2
   jmax_rho = ((lmax_rho_max+1)*(lmax_rho_max+2))/2
   do local_id = 1,LocalNumAtoms
!
      jmax_pot = Potential(local_id)%jmax
      n_Rpts = Potential(local_id)%jend
      jmt = Potential(local_id)%jmt
      v_tilt => Potential(local_id)%potL_Madelung
      v_tilt = CZERO
      Zi = getLocalAtomicNumber(local_id)
      r_mesh => Potential(local_id)%Grid%r_mesh(1:n_Rpts)
!
      if ( jmax_pot==1 .and. lmax_rho_max==0 ) then
         mad => getMadelungMatrix(local_id)
         dlm_0(1:GlobalNumAtoms,1) = cmplx(mad,ZERO,kind=RealKind)
         a0 = ONE
         dlm => dlm_0
         dlf_0(1) = ONE
      else
         dlm => getDLMatrix(local_id,a0)
      endif
!
      LoopJl: do jl_pot = 1,jmax_pot
         m_pot  = mofj(jl_pot)
         l_pot  = lofj(jl_pot)
         kl_pot = kofj(jl_pot)
!
         if ( jmax_pot==1 .and. lmax_rho_max==0 ) then
            dlf => dlf_0
         else
            dlf => getDLFactor(jl_pot)
         endif
!
         sumjl  = czero
         do kl_rho = kmax_rho, 1, -1
            m_rho = mofk(kl_rho)
            l_rho = lofk(kl_rho)
            jl_rho = ((l_rho+1)*(l_rho+2))/2-l_rho+abs(m_rho)
!
            mdif = m_rho-m_pot
            lsum = l_pot+l_rho
            kl = (lsum+1)*(lsum+1)-lsum+mdif
!
            ql => getMultipoleMoment(jl_rho,table_line)
!
            if (maxval(print_level) >= 1) then
               do ig = 1,GlobalNumAtoms
                  do ia = 1, getNumAlloyElements(ig)
                     lig = table_line(ig)+ia
                     if (l_rho == 0) then
                        write(6,'(a,i4,a,i4,a,2f18.14)')                          &
                           'ig = ',ig,', ia = ',ia,', Q0*Y0-Zi, Q0*Y0 = ',        &
                           real(ql(lig),kind=RealKind)*Y0-getAtomicNumber(ig,ia), &
                           real(ql(lig),kind=RealKind)*Y0
                        write(6,'(14x,a,f18.14)')   'Zi-Q0*Y0+rho0*Omega_mt = ',  &
                           getAtomicNumber(ig,ia)-real(ql(lig),kind=RealKind)*Y0+ &
                           getInterstitialElectronDensity()*getInscrSphVolume(ig)
                     endif
                  enddo
                  write(6,'(14x,a,2f18.14)')        'rho0, rho_neutral      = ',  &
                           getInterstitialElectronDensity(),rho_neutral*Y0      
                  if (m_rho < 0) then
                     write(6,'(a,6i4,a,2f18.14)')'i,j,L,Lp = ',ig,local_id, &
                          l_pot,m_pot,l_rho,m_rho,', alpha_M = ',           &
                          dlf(kl_rho)*m1m(m_rho)*dlm(ig,kl)/a0**lsum
                  else
                     write(6,'(a,6i4,a,2f18.14)')'i,j,L,Lp = ',ig,local_id, &
                          l_pot,m_pot,l_rho,m_rho,', alpha_M = ',dlf(kl_rho)*dlm(ig,kl)/a0**lsum
                  endif
               enddo
            endif
!
            sumat  = czero
            if ( m_rho<0 ) then
               do ig = 1,GlobalNumAtoms
                  ql_av = CZERO
                  do ia = 1, getNumAlloyElements(ig)
                     lig = table_line(ig)+ia
                     ql_av = ql_av + getAlloyElementContent(ig,ia)*ql(lig)
                  enddo
                  sumat = sumat + dlm(ig,kl)*m1m(m_rho)*conjg(ql_av)/a0**lsum
               enddo
            else
               if ( l_rho==0 ) then
                  do ig = 1,GlobalNumAtoms
                     ql_av = CZERO
                     do ia = 1, getNumAlloyElements(ig)
                        lig = table_line(ig)+ia
                        Zj = getAtomicNumber(ig,ia)
                        ql_av = ql_av + getAlloyElementContent(ig,ia)*(ql(lig)-Zj/Y0)
                     enddo
                     sumat = sumat + dlm(ig,kl)*ql_av/a0**lsum
                  enddo
               else
                  do ig = 1,GlobalNumAtoms
                     ql_av = CZERO
                     do ia = 1, getNumAlloyElements(ig)
                        lig = table_line(ig)+ia
                        ql_av = ql_av + getAlloyElementContent(ig,ia)*ql(lig)
                     enddo
                     sumat = sumat + dlm(ig,kl)*ql_av/a0**lsum
                  enddo
               endif
            endif
            sumjl = sumjl + dlf(kl_rho)*sumat
         enddo
!        sumjl = ZERO
!
         if ( l_pot==0 ) then
            do ir = 1,n_RPts
               v_tilt(ir,jl_pot) = v_tilt(ir,jl_pot)  +               &
                     cmplx( real( (TWO*sumjl -                        &
                                PI4*THIRD*rho_neutral*r_mesh(ir)**2), &
                                kind=RealKind), ZERO,kind=CmplxKind)
            enddo
            Potential(local_id)%VcoulombR0(2) = TWO*sumjl*Y0
         else if ( m_pot == 0 ) then
            do ir = 1,n_RPts
               v_tilt(ir,jl_pot) = v_tilt(ir,jl_pot) +                &
                 cmplx( real(TWO*m1m(l_pot)*sumjl*(r_mesh(ir)**l_pot),&
                             kind=RealKind), ZERO, kind=CmplxKind )
            enddo
         else
            do ir = 1,n_RPts
               v_tilt(ir,jl_pot) = v_tilt(ir,jl_pot) +                &
                    TWO*m1m(l_pot)*sumjl*(r_mesh(ir)**l_pot)
            enddo
         endif
!ywg     if ( isChargeComponentZero(local_id,jl_pot) ) then
!ywg        v_tilt(1:n_Rpts,jl_pot) = CZERO
!ywg     endif
!
      enddo LoopJl
!
   enddo
!
   nullify(table_line)
   deallocate( dlm_0 )
!
#ifdef TIMING
      t0 = getTime()-t0
      write(6,'(a,f10.5)')' calInterPlusMadPot :: Time ',t0
#endif
!
   end subroutine calInterPlusMadPot
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calFFTPseudoPot()
!  ===================================================================
   use MPPModule, only : NumPEs, MyPE
   use MPPModule, only : GlobalSum, setCommunicator, resetCommunicator
   use Uniform3DGridModule, only : getNumGridPoints
   use ParallelFFTModule, only : performTransformR2C, allocateFunctionSpace
   use ParallelFFTModule, only : getParaFFTCommunicator
   use ChargeDensityModule, only : isChargeComponentZero, getRhoLmax, &
                                   getChargeDensityAtPoint
   use InterpolationModule, only : FitInterp
   use MathParamModule, only : SQRT2
!
!use Uniform3DGridModule, only : getUniform3DGrid, getGridIndex, getGridPosition
!use ParallelFFTModule, only : getGridPointCoord
!use PublicTypeDefinitionsModule, only : UniformGridStruct
!
   implicit none
!
!character (len=4) :: ac
!
   integer (kind=IntKind) :: id, lmax, jmax, ir, nr, jl, jmin
   integer (kind=IntKind) :: lmax_rho, jmax_rho
   integer (kind=IntKind) :: ng, i, j, k, nlocal, comm, p, n
!
   real (kind=RealKind) :: pot_r, pseudo_fft, vc_0, dps
   real (kind=RealKind) :: dfp(3), dfp_total(3), dfp_corr(3), buf(2)
   real (kind=RealKind), pointer :: r_mesh(:)
!  real (kind=RealKind), pointer :: p_den(:,:,:)
   real (kind=RealKind), pointer :: p_den(:)
!
   complex (kind=CmplxKind), pointer :: v_jl(:,:)
!
   real (kind=RealKind) :: t0, t1, t2
!
!type (UniformGridStruct), pointer :: gp
!
#ifdef TIMING
   t0 = getTime()
#endif
!  ng  = fft_grid%NumLocalGridPoints
   ng  = getNumGridPoints('FFT',local_grid=.true.)
   call allocateFunctionSpace( p_den, nlocal )
!
!  gather the charge density on FFT grid
!
!  -------------------------------------------------------------------
   call constructDataOnGrid( 'FFT', 'Charge', 'Pseudo', getChargeDensityAtPoint, p_den )
!  -------------------------------------------------------------------
#ifdef TIMING
   t1 = getTime()
   write(6,*) "calFFTPseudoPot:: Time in UniformGrid: ",t1-t0
#endif
!
!  p_den => getDataStorage(StorageKey, nga, ngb, ngc, RealMark)
!
!  This summation needs to be carried out over processors with grid
!  distribution. - 05/15/18
   pseudo_fft = ZERO
   do i = 1, ng
      pseudo_fft = pseudo_fft + p_den(i)
   enddo
!  write(6,'(a,i5,2x,3d16.8)')'MyPE,p_den = ',MyPE,p_den(1),p_den(ng),pseudo_fft
!
   buf(1) = ng
   buf(2) = pseudo_fft
   comm = getParaFFTCommunicator(MyProc=p,NumProcs=n)
   if (comm > 0) then ! For serial FFT, comm = -1.
      call setCommunicator(comm,p,n)
      call GlobalSum(buf,2)
      call resetCommunicator()
   endif
!  The following buf(1) should be the total number of FFT grids.
   pseudo_fft = buf(2)/buf(1)
!
   if ( maxval(print_level) >= 0 ) then
      write(6,'(a,d16.8)') "Sum of pseudo charge density on uniform grid =",pseudo_fft
   endif
!
   pseudo_fft = ZERO ! In priciple, pseudo_fft should be zero since the pseudo charge
!                    ! is neutral. Here I set it to zero -ywg
   do i = 1, ng
      p_den(i) = p_den(i) - pseudo_fft
   enddo
!gp => getUniform3DGrid('FFT')
!write(ac,'(i4)')1000+NumPEs*100+MyPE
!ac(1:1) = 'p'
!open(unit=111,file='den-'//ac//'.dat',status='unknown',form='formatted')
!do i = 1, ng
!j = getGridIndex(gp,i)
!write(111,'(i8,2x,3f15.8,2x,3f15.8,2x,d16.8)')j,getGridPosition(gp,j),getGridPointCoord('R',i),p_den(i)
!enddo
!do i = 14, 15
!j = getGridIndex(gp,i)
!write(6,'(i8,2x,3f15.8,2x,3f15.8,2x,d16.8)')j,getGridPosition(gp,j),getGridPointCoord('R',i),p_den(i)
!enddo
! open(unit=111,file='den-2001-new.dat',status='unknown',form='formatted')
! write(111,'(5d16.8)')(p_den(i),i=1,ng)
! close(unit=111)
!nullify(gp)
!
!  -------------------------------------------------------------------
   call performTransformR2C(p_den,fft_c)
!  -------------------------------------------------------------------
!open(unit=112,file='kden-new.dat',status='unknown',form='formatted')
!do k = 1, numk_local
!write(112,'(i5,2x,3f15.8,2x,2d15.8)')k,getGridPointCoord('K',k),fft_c(k)
!enddo
!close(112)
!
#ifdef TIMING
   t2 = getTime()
   t1 = t2 - t1
   write(6,*) "calFFTPseudoPot:: Time in performTransformR2C: ",t1
#endif
!  -------------------------------------------------------------------
   call calPseudoDipoleField(fft_c, DF_Pseudo)
!  -------------------------------------------------------------------
   dfp_total = ZERO
   do id = 1,LocalNumAtoms
      dfp_total(1:3) = dfp_total(1:3) + DF_Pseudo(1:3,id)
!     write(6,'(a,3d16.8)'),'DF_Pseudo = ',DF_Pseudo(1:3,id)
   enddo
!
!  Correction is made to ensure total force acting on the unit cell is zero.
!  But I don't think this is necessary. -Yang
!  dfp_corr(1:3) = -dfp_total(1:3)/real(LocalNumAtoms,kind=RealKind)
   dfp_corr = ZERO
   do id = 1,LocalNumAtoms
      DF_Pseudo(1:3,id) =  DF_Pseudo(1:3,id) + dfp_corr(1:3)
   enddo
#ifdef TIMING
   t1 = getTime()
   t2 = t1 - t2
   write(6,*) "calFFTPseudoPot:: Time in Dipole: ",t2
#endif
!
!  -------------------------------------------------------------------
   call calRadialInterpolation(fft_c,-2,60,iparam,v_interp)
!  -------------------------------------------------------------------
   do id = 1,LocalNumAtoms
      lmax   = Potential(id)%lmax
      lmax_rho = getRhoLmax(id)
      jmax = ((lmax+1)*(lmax+2))/2
      jmax_rho = ((lmax_rho+1)*(lmax_rho+2))/2
      jmin = min(jmax,jmax_rho)
      nr     = Potential(id)%jend
      r_mesh => Potential(id)%Grid%r_mesh(1:nr)
      v_jl => Potential(id)%potL_Pseudo(1:nr,1:jmax)
      v_jl = CZERO
!     ----------------------------------------------------------------
      call calRadialProjection(nr,r_mesh,iparam(:,id),v_interp(:,id),v_jl)
!     ----------------------------------------------------------------
      do jl = 1,jmax
         if ( isChargeComponentZero(id,jl) .and. jl<=jmax_rho ) then
            v_jl(1:nr,jl) = CZERO
         else if ( lofj(jl)==0 ) then
!ywg     if ( lofj(jl)==0 ) then
            do ir = 1,nr
               pot_r = real(v_jl(ir,jl),kind=RealKind)
               v_jl(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
               V1_r(ir) = pot_r
            enddo
            call FitInterp(4,r_mesh(1:4),V1_r(1:4),ZERO,vc_0,dps)
            Potential(id)%VcoulombR0(3) = vc_0*Y0
!ywg        Potential(id)%VcoulombR0 = Potential(id)%VcoulombR0 +    &
!ywg                                   real(V_L0,kind=RealKind)*Y0
            if ( print_level(1) >=0 ) then
               write(6,'(a,d15.8,a,d15.8)')'Extrapolation: r = ',ZERO,     'pot_pseudo_0 = ',vc_0
               write(6,'(a,d15.8,a,d15.8)')'               r = ',r_mesh(1),'pot_pseudo_0 = ',V1_r(1)
               write(6,'(a,d15.8,a,d15.8)')'               r = ',r_mesh(2),'pot_pseudo_0 = ',V1_r(2)
               write(6,'(a,d15.8,a,d15.8)')'               r = ',r_mesh(3),'pot_pseudo_0 = ',V1_r(3)
               write(6,'(a,d15.8,a,d15.8)')'               r = ',r_mesh(4),'pot_pseudo_0 = ',V1_r(4)
               write(6,'(a,d15.8,a,d15.8)')'Pseudo Pot at Site = ',vc_0*Y0, &
                                           ', with extrapolation error = ',dps
!ywg           write(6,'(a,d15.8,a,d15.8)')'Old Pseudo V0 = ',vc_0*Y0,  &
!ywg                   ', New Pseudo V0 = ',real(V_L0,kind=RealKind)*Y0
            endif
         else if ( mofj(jl)==0 ) then
            do ir = 1,nr
               pot_r = real(v_jl(ir,jl),kind=RealKind)
               v_jl(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
            enddo
         endif
      enddo
   enddo
!
   deallocate(p_den); nullify(p_den)
#ifdef TIMING
   t2 = getTime()
   t1 = t2 - t1
   write(6,*) "calFFTPseudoPot:: Time in BackProjection: ",t1
   write(6,*) "calFFTPseudoPot:: Time : ",t2-t0
#endif
!
   end subroutine calFFTPseudoPot
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calPseudoDipoleField(p_fft_c, dfp)
!  ===================================================================
   use MathParamModule, only : SQRTm1
!
   use MPPModule, only : GlobalSum, setCommunicator, resetCommunicator
!
   use ParallelFFTModule, only : getNumLocalGridPoints,               &
                                 getLocalIndexOfGridOrigin,           &
                                 getGridPointCoord, getParaFFTCommunicator
!
   use SystemModule, only : getAtomPosition
!
   implicit none
!
   real (kind=RealKind), intent(out), target :: dfp(3,LocalNumAtoms)
   real (kind=RealKind), target :: dfp_buf(3,GlobalNumAtoms)
   real (kind=RealKind), pointer :: dfp_p(:,:)
   real (kind=RealKind) :: posi(3)
   real (kind=RealKind) :: kv2, kvec(3)
!
   complex (kind=CmplxKind), intent(in) :: p_fft_c(:)
   complex (kind=CmplxKind) :: expikR, cfact
!
   integer (kind=IntKind) :: idk, idk0, ng, comm, p, n, id, ig, na
!
   idk0 = getLocalIndexOfGridOrigin('K')
   if (idk0 > 1 .or. idk0 < 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('calPseudoDipoleField',                       &
                        'Grid origin does not correspond to the first local index',idk0)
!     ----------------------------------------------------------------
   endif
!
!  cfact = -SQRTm1
   cfact =  SQRTm1    ! 09/26/2018
!
   comm = getParaFFTCommunicator(MyProc=p,NumProcs=n)
   if (comm > 0) then ! For serial FFT, comm = -1.
      na = GlobalNumAtoms
      dfp_p => dfp_buf
   else
      na = LocalNumAtoms
      dfp_p => dfp
   endif
!
   dfp_p = ZERO
   do ig = 1, na
      if (comm > 0) then ! For serial FFT, comm = -1.
         posi = getAtomPosition(ig)
      else
         posi = LocalAtomPosi(:,ig)
      endif
      do idk =idk0+1, numk_local
         kvec = getGridPointCoord('K',idk)
         kv2 = kvec(1)**2+kvec(2)**2+kvec(3)**2
!        expikR  = exp( -SQRTm1*(kvec(1)*posi(1) + kvec(2)*posi(2) + kvec(3)*posi(3)) )
         expikR  = exp(  SQRTm1*(kvec(1)*posi(1) + kvec(2)*posi(2) + kvec(3)*posi(3)) ) ! 09/26/2018
         dfp_p(:,ig) = dfp_p(:,ig) + real(cfact*expikR*p_fft_c(idk),kind=RealKind)*kvec/kv2
      enddo
   enddo
   if (comm > 0) then ! For serial FFT, comm = -1.
!     ----------------------------------------------------------------
      call setCommunicator(comm,p,n)
      call GlobalSum(dfp_buf,3,GlobalNumAtoms)
      call resetCommunicator()
!     ----------------------------------------------------------------
      do id = 1, LocalNumAtoms
         ig = GlobalIndex(id)
         dfp(:,id) = dfp_p(:,ig)
      enddo
   endif
!
!  ===================================================================
!  dfp is the gradient of V:
!  dfp(1:3) = 4pi*e^2*sum_K [theta_K/K^2 * i*K * e^{iK*posi}]
!           =-4pi*e^2 * sum_K [ Im[theta_K*e^{iK*posi}] * K /K^2 ]
!           = 4pi*e^2 * sum_K [ Re[i*theta_K*e^{iK*posi}] * K /K^2 ]
!  ===================================================================
!  dfp = TWO*PI4*dfp/real(ng,kind=RealKind)
   dfp = TWO*PI4*dfp  ! The parallel FFT transformation already includes
                      ! a factor of 1/ng
!  ===================================================================
!
   end subroutine calPseudoDipoleField
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPseudoDipoleField(id) result(pdf)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), pointer :: pdf(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPseudoDipoleField','invalid id',id)
   endif
!
   pdf => DF_Pseudo(1:3,id)
!
   end function getPseudoDipoleField
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchangeJl(id,ia,lmax_in)
!  ===================================================================
!  use SystemSymmetryModule, only : getSymmetryFlags ! The symmetry break is likely 
!                                                      caused by the calculation
!                                                      of the ex-corr potential
   use MPPModule, only : MyPE
   use ChargeDensityModule, only : getChargeDensity, getMomentDensity, &
                                   getChargeDensityAtPoint,           &
                                   getMomentDensityAtPoint,           &
                                   getRhoLmax
   use ExchCorrFunctionalModule, only : calSphExchangeCorrelation
   use ExchCorrFunctionalModule, only : calExchangeCorrelation
   use ExchCorrFunctionalModule, only : getExchCorrEnDen
   use ExchCorrFunctionalModule, only : getExchCorrPot
!
   implicit none
!
   logical :: isNegCharge
!
   integer (kind=IntKind), intent(in) :: id, ia, lmax_in
!
   integer (kind=IntKind) :: ir, it, ip, is, jend, ir_lastNeg, ir_fit
   integer (kind=IntKind) :: l, m, kl, jl, lmax, jmax, kmax, n0
   integer (kind=IntKind) :: ngl_the, ngl_phi, ngl, ing, ip0, it0
!
!  integer (kind=IntKind), pointer :: flags_jl(:)
!
   real (kind=RealKind) :: posi(3), fact, sint, cost, sinp, cosp, r
   real (kind=RealKind) :: grad_rho(3), grad_mom(3)
   real (kind=RealKind) :: rho, mom, pXC_rtp, eXC_rtp, pXC_r0(2), eXC_r0(2)
   real (kind=RealKind), pointer :: r_mesh(:), V_tmp(:)
   real (kind=RealKind), pointer :: theta(:), phi(:)
   real (kind=RealKind), pointer :: wght_the(:), wght_phi(:)
   real (kind=RealKind), pointer :: chgmom_data(:,:,:,:,:)
   real (kind=RealKind), allocatable :: der_rho_tmp(:), der_mom_tmp(:)
!#ifdef TIMING
   real (kind=RealKind) :: t0, t1, t2, ts, tc1, tc2, tc3
!#endif
!
   complex (kind=CmplxKind) :: a, b, y1, y2, x1, x2
   complex (kind=CmplxKind), pointer :: rhol(:,:), der_rhol(:,:)
   complex (kind=CmplxKind), pointer :: moml(:,:), der_moml(:,:)
   complex (kind=CmplxKind), pointer :: pylm(:), V_exc(:)
   complex (kind=CmplxKind), pointer :: potL_Exch(:,:,:), enL_Exch(:,:,:)
   complex (kind=CmplxKind), pointer :: ylm_ngl(:,:)
!
!  -------------------------------------------------------------------
!  generate points on the sphere to be integrated
!
!  integrate over the surface
!                       ->     ->      ^   ^
!  Vxc_L(r)=int{Vxc[rho(r),mom(r)]*Y_L(r)}dr
!                       ->     ->      ^   ^
!  Exc_L(r)=int{Exc[rho(r),mom(r)]*Y_L(r)}dr
!  -------------------------------------------------------------------
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
#ifdef TIMING
   t0 = getTime()
#endif
!
   if ( lmax_in<0 ) then
      lmax = Potential(id)%lmax
      jmax = Potential(id)%jmax
      kmax = (lmax+1)*(lmax+1)
   else
      lmax = min(lmax_in,Potential(id)%lmax)
      jmax = ((lmax+1)*(lmax+2))/2
      kmax = (lmax+1)*(lmax+1)
   endif
   jend  = Potential(id)%jend
   r_mesh => Potential(id)%Grid%r_mesh
!
   potL_Exch => Potential(id)%potL_Exch(:,:,:,ia)
   potL_Exch = CZERO
   enL_Exch  => Potential(id)%enL_Exch(:,:,:,ia)
   enL_Exch  = CZERO
!
!  Early exit
!
   if ( lmax_in==0 ) then
      allocate(der_rho_tmp(jend), der_mom_tmp(jend))
      if (gga_functional) then
         rhol => getChargeDensity("TotalNew",id,ia,der_rhol)
         do ir = 1,jend
            der_rho_tmp(ir) = real(der_rhol(ir,1),kind=RealKind)*Y0
         enddo
      else
         rhol => getChargeDensity("TotalNew",id,ia)
         der_rho_tmp = ZERO
      endif
      do ir = 1,jend
         rho_tmp_r(ir) = real(rhol(ir,1),kind=RealKind)*Y0
      enddo
      if ( n_spin_pola == 2 ) then
         if (gga_functional) then
            moml => getMomentDensity("TotalNew",id,ia,der_moml)
            do ir = 1,jend
               der_mom_tmp(ir) = real(der_moml(ir,1),kind=RealKind)*Y0
            enddo
         else
            moml => getMomentDensity("TotalNew",id,ia)
            der_mom_tmp = ZERO
         endif
         do ir = 1,jend
            mom_tmp_r(ir) = real(moml(ir,1),kind=RealKind)*Y0
         enddo
         if (gga_functional) then
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(jend,rho_tmp_r(1:),        &
                                           der_rho_den=der_rho_tmp,   &
                                           mag_den=mom_tmp_r(1:),     &
                                           der_mag_den=der_mom_tmp)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(jend,rho_tmp_r(1:),        &
                                           mag_den=mom_tmp_r(1:))
!           ----------------------------------------------------------
         endif
      else
         if (gga_functional) then
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(jend,rho_tmp_r(1:),        &
                                           der_rho_den=der_rho_tmp)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call calSphExchangeCorrelation(jend,rho_tmp_r(1:))
!           ----------------------------------------------------------
         endif
      endif
      do is = 1,n_spin_pola
         V_exc => Potential(id)%potL_Exch(1:jend,1,is,ia)
         V_tmp => getExchCorrPot(jend,is)
         do ir = 1,jend
            V_exc(ir) = cmplx(V_tmp(ir)/Y0,ZERO,kind=CmplxKind)
         enddo
         V_exc => Potential(id)%enL_Exch(1:jend,1,is,ia)
         V_tmp => getExchCorrEnDen(jend)
         do ir = 1,jend
            V_exc(ir) = cmplx(V_tmp(ir)/Y0,ZERO,kind=CmplxKind)
         enddo
      enddo
      deallocate(der_rho_tmp, der_mom_tmp)
!
      return
   endif
!
   indrl_fit = 0
!
   ngl_the = AngularData%n_theta
   ngl_phi = AngularData%n_phi
   ngl = AngularData%ngl
!
   phi      => AngularData%phi(1:ngl_phi)
   theta    => AngularData%theta(1:ngl_the)
   wght_phi => AngularData%wght_phi(1:ngl_phi)
   wght_the => AngularData%wght_the(1:ngl_the)
   ylm_ngl  => AngularData%ylm_ngl(1:kmax,1:ngl)
!
   chgmom_data => AngularData%radial_data(1:n_spin_pola,1:2,1:ngl_the,1:ngl_phi,1:jend)
!
#ifdef TIMING
    t1 = getTime()
    write(6,*) "calExchangeJl:: Init Time : ",t1-t0
#endif
!
   tc1 = ZERO; tc2 = ZERO; tc3 = ZERO
   ts = getTime()
   ir_lastNeg = 0
   do ir = 1,jend
      isNegCharge = .false.
!     pXC_r(1:jmax,1:n_spin_pola) = CZERO
!     eXC_r(1:jmax,1:n_spin_pola) = CZERO
      pXC_r = CZERO; eXC_r = CZERO
      Loop_phi01: do ip0 = 1,ngl_phi
         if ( mod(ip0,2)==0 ) then
            ip = ngl_phi-ip0/2+1
         else
            ip = (ip0+1)/2
         endif
         cosp  = cos(phi(ip))
         sinp  = sin(phi(ip))
!        pXC_rphi(1:jmax,1:n_spin_pola) = CZERO
!        eXC_rphi(1:jmax,1:n_spin_pola) = CZERO
         pXC_rphi = CZERO; eXC_rphi = CZERO
         Loop_the01: do it0 = 1,ngl_the
            if ( mod(it0,2)==0 ) then
               it = ngl_the-it0/2+1
            else
               it = (it0+1)/2
            endif
!
! The following code appears alright: theta corresonds to
! Phi, and phi corresponds to Theta. Just a bad way in naming the variables. -Yang Wang
!
            sint = sin(theta(it))
            cost = cos(theta(it))
            posi(1) = sint*sinp
            posi(2) = cost*sinp
!ywg 11/14/18
  posi(1) = cost*sinp
  posi(2) = sint*sinp
!
            posi(3) = cosp
            posi(1:3) = r_mesh(ir)*posi(1:3)
            t2 = getTime()
            if (gga_functional) then
               rho = getChargeDensityAtPoint( 'TotalNew', id, ia, posi, grad=grad_rho )
            else
               rho = getChargeDensityAtPoint( 'TotalNew', id, ia, posi )
               grad_rho = ZERO
            endif
            tc1 = tc1 + (getTime()-t2)  ! cummulating time on getChargeDensityAtPoint
            if ( rho <= ZERO ) then
               cycle Loop_the01
            endif
            t2 = getTime()
            if ( n_spin_pola==2 ) then
               if (gga_functional) then
                  mom = getMomentDensityAtPoint( 'TotalNew', id, ia, posi, grad=grad_mom )
!                 ----------------------------------------------------
                  call calExchangeCorrelation(rho,grad_rho,mom,grad_mom)
!                 ----------------------------------------------------
               else
                  mom = getMomentDensityAtPoint( 'TotalNew', id, ia, posi )
!                 ----------------------------------------------------
                  call calExchangeCorrelation(rho,mag_den=mom)
!                 ----------------------------------------------------
               endif
            else
               if (gga_functional) then
!                 ----------------------------------------------------
                  call calExchangeCorrelation(rho,grad_rho_den=grad_rho)
!                 ----------------------------------------------------
               else
!                 ----------------------------------------------------
                  call calExchangeCorrelation(rho)
!                 ----------------------------------------------------
               endif
            endif
            tc2 = tc2 + (getTime()-t2)  ! cummulating time on calExchangeCorrelation
            t2 = getTime()
            do is = 1,n_spin_pola
!              -------------------------------------------------------
               pXC_rtp = getExchCorrPot(is)
               eXC_rtp = getExchCorrEnDen()
!              -------------------------------------------------------
               chgmom_data(is,1,it0,ip0,ir) = pXC_rtp
               chgmom_data(is,2,it0,ip0,ir) = eXC_rtp
               pXC_rtp = pXC_rtp*wght_the(it)
               pXC_rphi(1,is) = pXC_rphi(1,is) + Y0*cmplx(pXC_rtp,ZERO,Kind=CmplxKind)
               eXC_rtp = eXC_rtp*wght_the(it)
               eXC_rphi(1,is) = eXC_rphi(1,is) + Y0*cmplx(eXC_rtp,ZERO,Kind=CmplxKind)
            enddo
            tc3 = tc3 + (getTime()-t2)  ! cummulating time on getExchCorrPot and getExchCorrEnDen
         enddo Loop_the01
         fact = sinp*wght_phi(ip)
         do is = 1,n_spin_pola
            pXC_r(1,is) = pXC_r(1,is) + fact*pXC_rphi(1,is)
            eXC_r(1,is) = eXC_r(1,is) + fact*eXC_rphi(1,is)
         enddo
      enddo Loop_phi01
!
      do is = 1,n_spin_pola
         potL_Exch(ir,1,is) = pXC_r(1,is)
         enL_Exch(ir,1,is)  = eXC_r(1,is)
!
#ifdef DEBUG
         if ( ir==999 ) then
            write(6,'(a,i4,a,4d16.8)') "id:",id,"Pot-En_Exchange MT-1:", pXC_r(1,is),eXC_r(1,is)
         endif
#endif
!
      enddo
!
   enddo
   if (MyPE == 0) then
      write(6,'(/,a,f10.5)')'Time:: calExchangeJl::tc1: ',tc1
      write(6,'(  a,f10.5)')'Time:: calExchangeJl::tc2: ',tc2
      write(6,'(  a,f10.5)')'Time:: calExchangeJl::tc3: ',tc3
      write(6,'(  a,f10.5/)')'Time:: calExchangeJl::Loop_ir 1: ',getTime()-ts
   endif
!
#ifdef TIMING
   t2 = getTime()
   write(6,*) "calExchangeJl:: l=0 Time : ",t2-t1
#endif
!
   if (jmax>1) then
!     flags_jl => getSymmetryFlags(id)
      r = 0.0005d0
!     ----------------------------------------------------------------
      call hunt(jend,r_mesh(1:jend),r,ir_fit)
!     ----------------------------------------------------------------
      ir_lastNeg = 0
      indrl_fit = ir_fit
!
      ts = getTime()
!
      do ir = jend,ir_fit,-1
         ing=0
         if ( ir < Potential(id)%jend ) then
            do is = 1,n_spin_pola
               pXC_r0(is) = Y0*real(potL_Exch(ir,1,is),kind=RealKind)
               eXC_r0(is) = Y0*real(enL_Exch(ir,1,is),kind=RealKind)
            enddo
         else
            pXC_r0 = ZERO
            eXC_r0 = ZERO
         endif
         pXC_r(1:jmax,1:n_spin_pola) = CZERO
         eXC_r(1:jmax,1:n_spin_pola) = CZERO
         Loop_phi1: do ip0 = 1,ngl_phi
            if ( mod(ip0,2)==0 ) then
               ip = ngl_phi-ip0/2+1
            else
               ip = (ip0+1)/2
            endif
            cosp  = cos(phi(ip))
            sinp  = sin(phi(ip))
            pXC_rphi(1:jmax,1:n_spin_pola) = CZERO
            eXC_rphi(1:jmax,1:n_spin_pola) = CZERO
            Loop_the1: do it0 = 1,ngl_the
               if ( mod(it0,2)==0 ) then
                  it = ngl_the-it0/2+1
               else
                  it = (it0+1)/2
               endif
               ing = (ip-1)*ngl_the+it
               sint = sin(theta(it))
               cost = cos(theta(it))
               posi(1) = sint*sinp
               posi(2) = cost*sinp
               posi(3) = cosp
               posi(1:3) = r_mesh(ir)*posi(1:3)
               pylm => ylm_ngl(1:kmax,ing)
!
               do is = 1,n_spin_pola
                  pXC_rtp = chgmom_data(is,1,it0,ip0,ir)
                  eXC_rtp = chgmom_data(is,2,it0,ip0,ir)
!
                  pXC_rtp = (pXC_rtp -pXC_r0(is))*wght_the(it)
                  eXC_rtp = (eXC_rtp -eXC_r0(is))*wght_the(it)
!                  pXC_rtp = pXC_rtp*wght_the(it)
!                  eXC_rtp = eXC_rtp*wght_the(it)
!
                  do jl = 2,jmax
!     if (flags_jl(jl) > 0) then   ! Added by Yang @ Dec 27, 2014
                     l = lofj(jl)
                     m = mofj(jl)
                     kl = (l+1)*(l+1)-l+mofj(jl)
                     fact = ONE
                     pXC_rphi(jl,is) = pXC_rphi(jl,is) + fact*          &
!                     cmplx(pXC_rtp,ZERO,Kind=CmplxKind)*m1m(m)*conjg(pylm(kl))/&
                     cmplx(pXC_rtp,ZERO,Kind=CmplxKind)*conjg(pylm(kl))/&
                     r_mesh(ir)**l
                     eXC_rphi(jl,is) = eXC_rphi(jl,is) + fact*          &
!                     cmplx(eXC_rtp,ZERO,Kind=CmplxKind)*m1m(m)*conjg(pylm(kl))/&
                     cmplx(eXC_rtp,ZERO,Kind=CmplxKind)*conjg(pylm(kl))/&
                     r_mesh(ir)**l
!     endif
                  enddo
               enddo
!
            enddo Loop_the1
!
            fact = sinp*wght_phi(ip)
            do is = 1,n_spin_pola
               do jl = 2,jmax
                  pXC_r(jl,is) = pXC_r(jl,is) + fact*pXC_rphi(jl,is)
                  eXC_r(jl,is) = eXC_r(jl,is) + fact*eXC_rphi(jl,is)
              enddo
            enddo
!
         enddo Loop_phi1
!
         do is = 1,n_spin_pola
            do jl = 2,jmax
               l = lofj(jl)
               pXC_r(jl,is) = pXC_r(jl,is)*r_mesh(ir)**l
               eXC_r(jl,is) = eXC_r(jl,is)*r_mesh(ir)**l
            enddo
         enddo
!
         do is = 1,n_spin_pola
            do jl = 2,jmax
               potL_Exch(ir,jl,is) = pXC_r(jl,is)
               enL_Exch(ir,jl,is)  = eXC_r(jl,is)
            enddo
         enddo
!
         do is = 1,n_spin_pola
            do jl = 2,jmax
               l = lofj(jl)
               if ( r_mesh(ir)<0.1d0 .and. abs(pXC_r(jl,is))<5.0d0*ten2m12 &
                    .and. ( indrl_fit(jl,is)<=ir_fit+1 ) ) then
                  indrl_fit(jl,is) = ir+1
               endif
            enddo
         enddo
!
      enddo
      if (MyPE == 0) then
         write(6,'(/,a,f10.5/)')'Time:: calExchangeJl::Loop_ir 2: ',getTime()-ts
      endif
!
      ts = getTime()
      n0 = 4
      do is = 1,n_spin_pola
         do jl = 2,jmax
            if ( indrl_fit(jl,is) /= 0 ) then
               l = lofj(jl)
               ir = indrl_fit(jl,is)
               x1 = r_mesh(ir+1)
               x2 = r_mesh(ir+n0/2)
!
               y1 = potL_Exch(ir+1,jl,is)/(r_mesh(ir+1)**l)
               y2 = potL_Exch(ir+n0/2,jl,is)/(r_mesh(ir+n0/2)**l)
               a  = (y2-y1)/(x2-x1)
               b  = (y1*x2-y2*x1)/(x2-x1)
               do ir = indrl_fit(jl,is)+n0/2-1,1,-1
                  potL_Exch(ir,jl,is) = (a*r_mesh(ir)+b)*(r_mesh(ir)**l)
               enddo
!
               ir = indrl_fit(jl,is)
               y1 = enL_Exch(ir+1,jl,is)/(r_mesh(ir+1)**l)
               y2 = enL_Exch(ir+n0/2,jl,is)/(r_mesh(ir+n0/2)**l)
               a  = (y2-y1)/(x2-x1)
               b  = (y1*x2-y2*x1)/(x2-x1)
               do ir = indrl_fit(jl,is)+n0/2-1,1,-1
                  enL_Exch(ir,jl,is) = (a*r_mesh(ir)+b)*(r_mesh(ir)**l)
               enddo
!
            endif
         enddo
      enddo
      if (MyPE == 0) then
         write(6,'(/,a,f10.5/)')'Time:: calExchangeJl::Loop_is : ',getTime()-ts
      endif
!
   endif
!
#ifdef TIMING
   t0 = getTime()-t0
   write(6,*) "calExchangeJl:: Time : ",t0
#endif
   end subroutine calExchangeJl
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setAngularData(nt,np,ll,ns,nr)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nt, np, ll, ns, nr
!
   integer (kind=IntKind) :: ngl, ing, ip, it, kl
!
   real (kind=RealKind) :: sint, cost, sinp, cosp, posi(3)
   real (kind=RealKind), pointer :: theta(:), phi(:)
   real (kind=RealKind), pointer :: wght_the(:), wght_phi(:)
   complex (kind=CmplxKind), pointer :: ylm_ngl(:,:)
!
   kl = (ll+1)*(ll+1)
   AngularData%kmax    = kl
   AngularData%n_theta = nt
   AngularData%n_phi   = np
   AngularData%ngl     = nt*np
   ngl = AngularData%ngl 
!
   allocate( AngularData%theta(nt), AngularData%phi(np),              &
             AngularData%wght_phi(np), AngularData%wght_the(nt) )
   allocate( AngularData%ylm_ngl(1:AngularData%kmax,1:ngl))
!
   phi   => AngularData%phi
   theta => AngularData%theta
   wght_the => AngularData%wght_the
   wght_phi => AngularData%wght_phi
   ylm_ngl => AngularData%ylm_ngl
!
   do it = 1,nt
      wght_the(it) = PI2/nt
      theta(it)    = (it-1)*wght_the(it)+HALF*wght_the(it)
   enddo
!   do ip = 1,np
!      wght_phi(ip) = PI/np
!      phi(ip)    = (ip-1)*wght_phi(ip)+HALF*wght_phi(ip)
!   enddo
   call  gauleg(ZERO, PI, phi(1:np), wght_phi(1:np), np)
!   call  gauleg(ZERO, PI/2, phi(1:np/2), wght_phi(1:np/2), np/2)
!   call  gauleg(Pi/2, PI, phi(np/2+1:np), wght_phi(np/2+1:np), np/2)
!   call  gauleg(ZERO, PI, theta(1:nt/2), wght_the(1:nt/2), nt/2)
!   call  gauleg(PI, PI2, theta(nt/2+1:nt), wght_the(nt/2+1:nt), nt/2)
!
   if ( kl>1 ) then!
      ing=0
      do ip = 1,np
         cosp  = cos(phi(ip))
         sinp  = sin(phi(ip))
         do it = 1,nt
            ing = ing+1
            sint = sin(theta(it))
            cost = cos(theta(it))
            posi(1) = sint*sinp
            posi(2) = cost*sinp
!
!ywg 11/14/18
  posi(1) = cost*sinp
  posi(2) = sint*sinp
!
            posi(3) = cosp
!           ----------------------------------------------------
            call calYlm(posi, ll, ylm_ngl(1:kl,ing) )
!           ----------------------------------------------------
         enddo
      enddo
   endif
!
   allocate(AngularData%radial_data(ns,2,nt,np,nr))
!
   end subroutine setAngularData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine clearAngularData()
!  ===================================================================
   implicit none
! 
   deallocate( AngularData%theta, AngularData%phi,                    &
               AngularData%wght_phi, AngularData%wght_the )
   deallocate( AngularData%ylm_ngl, AngularData%radial_data )
!
   end subroutine clearAngularData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPot_L( id, jmax_in, pot_type, flag_rl, aux_name )
!  ===================================================================
   use MathParamModule, only: Ten2m12
   use MPPModule, only : MyPE
!
   implicit none
!
   character (len=*), intent(in), optional :: pot_type, aux_name
!
   integer (kind=IntKind), intent(in), optional :: flag_rl
   integer (kind=IntKind), intent(in) :: id, jmax_in
!
   character(len=50) :: file_potl
   character(len=9)  :: char_lm
   integer (kind=IntKind) :: l, m, jl, is, jmax, ir, NumRs, ia
   integer (kind=IntKind) :: funit, ns, nc, len_potname
   integer (kind=IntKind) :: offset    = 100000
   integer (kind=IntKind) :: offset_at = 1000
   integer (kind=IntKind) :: offset_lm = 100
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: potl4(:,:,:,:), potl3(:,:,:), potl2(:,:)
   complex (kind=CmplxKind), pointer :: potl(:,:)
!
   if (isMTFP .and. nocaseCompare(pot_type,"Pseudo") ) then
      call ErrorHandler('printPot_L','invalid potential type','Pseudo')
   else if (isMTFP .and. nocaseCompare(pot_type,"XCHat") ) then
      call ErrorHandler('printPot_L','invalid potential type','XCHat')
   else if (isMTFP .and. nocaseCompare(pot_type,"En_XCHat") ) then
      call ErrorHandler('printPot_L','invalid potential type','En_XCHat')
   endif
!
   file_potl = " "
   if ( .not.present(pot_type) ) then
      if ( present(flag_rl) ) then
         if (flag_rl==0) then
            nc = 6
            file_potl(1:nc)  = "PotNL_"
         else
            nc = 7
            file_potl(1:nc)  = "PotNRL_"
         endif
      else
         nc = 6
         file_potl(1:nc)  = "PotNL_"
      endif
   else
      len_potname = len(pot_type)
      if ( present(flag_rl) ) then
         if (flag_rl==0) then
            nc = 5+len_potname
            file_potl(4+len_potname:nc)  = "L_"
         else
            nc = 5+len_potname+1
            file_potl(4+len_potname:nc)  = "RL_"
         endif
      else
         nc = 5+len_potname
         file_potl(4+len_potname:nc)  = "L_"
      endif
      file_potl(1:3)  = "Pot"
      file_potl(4:3+len_potname)  = pot_type(1:len_potname)
   endif
   if ( present(aux_name) ) then
      len_potname = len(trim(aux_name))
      file_potl(nc+1:nc+len_potname) = aux_name(1:len_potname)
      nc = nc + len_potname
   endif
!
   write(file_potl(nc+1:nc+6),'(i6)') offset+MyPE
   file_potl(nc+1:nc+1) = 'n'
   file_potl(nc+7:nc+7) = '_'
   write(file_potl(nc+8:nc+11),'(i4)') offset_at+id
   file_potl(nc+8:nc+8) = 'a'
   funit = 55+MyPE+id
   open(unit=funit,file=trim(file_potl),status='unknown')
!
   NumRs = Potential(id)%Grid%jend
   r_mesh => Potential(id)%Grid%r_mesh(1:NumRs)
   if ( jmax_in>0 ) then
      jmax =  min(jmax_in,Potential(id)%jmax)
   else
      jmax = Potential(id)%jmax
   endif
!
   if ( present(pot_type) ) then
      if ( nocaseCompare(pot_type,"Total") ) then
         potl4 => Potential(id)%potL
         flag_jl => Potential(id)%PotCompFlag
         ns = 4
      else if ( nocaseCompare(pot_type,"Tilda") ) then
         potl3 => Potential(id)%potL_Tilda
         flag_jl => Potential(id)%potL_Tilda_flag
         ns = 3
      else if ( nocaseCompare(pot_type,"Coulomb") ) then
         potl3 => Potential(id)%potL_Coulomb
         flag_jl => Potential(id)%potL_Coulomb_flag
         ns = 3
      else if ( nocaseCompare(pot_type,"Pseudo") ) then
         potl2 => Potential(id)%potL_Pseudo
         flag_jl => Potential(id)%potL_Pseudo_flag
         ns = 2
      else if ( nocaseCompare(pot_type,"Madelung") ) then
         potl2 => Potential(id)%potL_Madelung
         flag_jl => Potential(id)%potL_Madelung_flag
         ns = 2
      else if ( nocaseCompare(pot_type,"Exchg") ) then
         potl4 => Potential(id)%potL_Exch
         flag_jl => Potential(id)%potL_Exch_flag
         ns = 4
      else if ( nocaseCompare(pot_type,"XCHat") ) then
         potl4 => Potential(id)%potL_XCHat
         flag_jl => Potential(id)%potL_XCHat_flag
         ns = 4
      else if ( nocaseCompare(pot_type,"En_Exchg") ) then
         potl4 => Potential(id)%enL_Exch
         flag_jl => Potential(id)%potL_Exch_flag
         ns = 4
      else if ( nocaseCompare(pot_type,"En_XCHat") ) then
         potl4 => Potential(id)%enL_XCHat
         flag_jl => Potential(id)%potL_XCHat_flag
         ns = 4
      else
         potl4 => Potential(id)%potL
         flag_jl => Potential(id)%PotCompFlag
         ns = 4
      endif
   else
      potl4 => Potential(id)%potL
      flag_jl => Potential(id)%PotCompFlag
      ns = 4
   endif
!
   do is = 1,n_spin_pola
      write(funit,'(a,$)') "# Ind       r_mesh"
      do jl = 1,jmax
         if ( flag_jl(jl) /= 0 ) then
            l = lofj(jl)
            m = mofj(jl)
            char_lm(1:3) = "ReL"
            write(char_lm(4:6),'(i3)') offset_lm + l
            char_lm(4:4) = '_'
            write(char_lm(7:9),'(i3)') offset_lm + m
            char_lm(7:7) = '_'
            write(funit,'(9x,a9,$)') char_lm
            char_lm(1:3) = "ImL"
            write(funit,'(8x,a9,$)') char_lm
         endif
      enddo
      write(funit,'(a)') "  "
      do ia = 1, Potential(id)%NumSpecies
         if (ns == 4) then
            potl => potl4(:,:,is,ia)
         else if (ns == 3) then
            potl => potl3(:,:,ia)
         else if (ns == 2) then
            potl => potl2(:,:)
         endif
         do ir = 1, NumRs
            write(funit,'(i4,1x,i4,1x,d16.8,$)') id,ia,r_mesh(ir)
            JlLoop: do jl = 1,jmax
               if ( flag_jl(jl) /= 0 ) then
                  l = lofj(jl)
                  if ( present(flag_rl) ) then
                     if (flag_rl==0) then
                         write(funit,'(1x,(2(1x,d16.8)),$)') potl(ir,jl) ! /(r_mesh(ir)**l)
                     else
                         write(funit,'(1x,(2(1x,d16.8)),$)') potl(ir,jl)/(r_mesh(ir)**l)
                     endif
                  else
                     write(funit,'(1x,(2(1x,d16.8)),$)') potl(ir,jl) ! /(r_mesh(ir)**l)
                  endif
               endif
            enddo JlLoop
            write(funit,'(a)') " "
         enddo
      enddo
   enddo
!
   nullify(flag_jl)
   close(funit)
!
   end subroutine printPot_L
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDensity_r(denType, id, ia, nrs, r_mesh, den_r)
!  ===================================================================
   use MathParamModule, only: Ten2m8
   use MPPModule, only : MyPE
!
   implicit none
!
   character (len=*), intent(in) :: denType
!
   integer (kind=IntKind), intent(in) :: id, ia, nrs
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), pointer :: den_r(:)
!
   character(len=20) :: file_den
   integer (kind=IntKind) :: lenDenName
   integer (kind=IntKind) :: ir, funit
   integer (kind=IntKind) :: offset = 100000
!
   if (MyPE /= 0) then
     return
   endif
   lenDenName = len(denType)+1
   file_den(1:20) = "                    "
   file_den(1:lenDenName-1) = trim(denType)
   file_den(lenDenName:lenDenName) = "_"
!
   if ( id == 1 ) then
      write(file_den(lenDenName+1:lenDenName+6),'(i6)') offset+MyPE+id*10+ia
      file_den(lenDenName+1:lenDenName+1) = 'n'
      funit = 55+MyPE+id*10+ia
      open(unit=funit,file=trim(file_den),status='unknown')
      write(funit,'(a)') "#Ind   atom    r_mesh   (lm):"
      write(funit,'(a)') "  "
      do ir = 1, nrs
         write(funit,'(2i5,2(1x,d16.8))') id, ia, r_mesh(ir), den_r(ir)
      enddo
      write(funit,'(a)') " "
      close(funit)
   endif
!
   end subroutine printDensity_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotentialAtPoint( potType, id, ia, posi, jmax_in, spin, grad ) result(pot)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   use InterpolationModule, only : PolyInterp
!
   use PotentialTypeModule, only : isMuffinTinPotential
!
   use AtomModule, only : getLocalAtomicNumber
!
   use PolyhedraModule, only: isExternalPoint
!
   implicit none
!
   character (len=*), intent(in) :: potType
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in), optional :: jmax_in
   integer (kind=IntKind), intent(in), optional :: spin
!
   real (kind=RealKind), intent(in) :: posi(3)
   real (kind=RealKind), intent(out), optional :: grad(3) ! So far it is not used
!
   integer (kind=IntKind) :: iend, irp, ir, l, kl, jl, ns, jmt, jmax
!
   real (kind=RealKind) :: r, pot, err
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: pot_in
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPotentialAtPoint','invalid site index',id)
   else if (isMTFP) then
      if ( nocaseCompare(potType,"Pseudo") .or. &
           nocaseCompare(potType,"XCHat") .or.  &
           nocaseCompare(potType,"En_XCHat") ) then
         call ErrorHandler('getPotentialAtPoint','invalid potential type for MTFP case',potType)
      endif
   else if (ia < 1 .or. ia > Potential(id)%NumSpecies) then
      call ErrorHandler('getPotentialAtPoint','invalid atom index on the site',ia,id)
   endif
!
   if (present(jmax_in)) then
      if ( jmax_in<1 .or. jmax_in>Potential(id)%jmax ) then
         call ErrorHandler('getPotentialAtPoint','invalid jmax_in',jmax_in)
      endif
      jmax = jmax_in
   else
      jmax = Potential(id)%jmax
   endif
!
   pot = ZERO
!
   if ( isExternalPoint(id,posi(1),posi(2),posi(3)) ) then
      return
   endif
!
   iend = Potential(id)%Grid%jend
   r = sqrt( posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
   if ( r > Potential(id)%Grid%r_mesh(iend) ) then
      return
   endif
!
   r_mesh => Potential(id)%Grid%r_mesh(1:iend)
!  -------------------------------------------------------------------
   call hunt(iend,r_mesh(1:iend),r,ir)
!  -------------------------------------------------------------------
   if (ir > iend-(n_inter-1)/2) then
      irp=iend-n_inter+1
   else if (2*ir+1 > n_inter) then
      irp=ir-(n_inter-1)/2
   else
      irp=1
   endif
!
   if ( isMuffinTinPotential() .and. Potential(id)%jmax==1 ) then
      jmt = Potential(id)%Grid%jmt
      if ( ir > jmt ) then
         return
      else if ( ir>jmt-(n_inter-1)/2 ) then
         irp = jmt-n_inter+1
      endif
   endif
!
!  -------------------------------------------------------------------
   call calYlm(posi,Potential(id)%lmax,Ylm)
!  -------------------------------------------------------------------
!
   if ( .not.present(spin) ) then
      do ns = 1,n_spin_pola
         if ( nocaseCompare(potType,"Tilda") ) then
            pot_l => Potential(id)%potL_Tilda(1:iend,1:jmax,ia)
         else if ( nocaseCompare(potType,"Madelung") ) then
            pot_l => Potential(id)%potL_Madelung(1:iend,1:jmax)
         else if ( nocaseCompare(potType,"Total") ) then
            pot_l => Potential(id)%potL(1:iend,1:jmax,ns,ia)
         else if ( nocaseCompare(potType,"Pseudo") ) then
            pot_l => Potential(id)%potL_Pseudo(1:iend,1:jmax)
         else if ( nocaseCompare(potType,"Exchg") ) then
            pot_l => Potential(id)%potL_Exch(1:iend,1:jmax,ns,ia)
         else if ( nocaseCompare(potType,"Coulomb") ) then
            pot_l => Potential(id)%potL_Coulomb(1:iend,1:jmax,ia)
         else if ( nocaseCompare(potType,"XCHat") ) then
            pot_l => Potential(id)%potL_XCHat(1:iend,1:jmax,ns,ia)
         else if ( nocaseCompare(potType,"En_Exchg") ) then
            pot_l => Potential(id)%enL_Exch(1:iend,1:jmax,ns,ia)
         else if ( nocaseCompare(potType,"En_XCHat") ) then
            pot_l => Potential(id)%enL_XCHat(1:iend,1:jmax,ns,ia)
         else
            call ErrorHandler('getPotentialAtPoint','invalid potential type',potType)
!           pot_l => Potential(id)%potL(1:iend,1:jmax,ns)
         endif
!
         do jl = jmax,1,-1
            l = lofj(jl)
!           ----------------------------------------------------------
            call PolyInterp(n_inter, r_mesh(irp:irp+n_inter-1),       &
                         pot_l(irp:irp+n_inter-1,jl), r, pot_in, err)
!           ----------------------------------------------------------
            kl = (l+1)*(l+1)-l+mofj(jl)
            if (mofj(jl) == 0) then
               pot = pot + real(pot_in*Ylm(kl),RealKind)
            else
               pot = pot + TWO*real(pot_in*Ylm(kl),RealKind)
            endif
         enddo
      enddo
      pot = pot/(ONE+real(n_spin_pola-1,kind=RealKind))
   else
      if ( spin<1 .or. spin>2 ) then
         call ErrorHandler('getPotentialAtPoint','invalid spin index',spin)
      endif
      if ( nocaseCompare(potType,"Tilda") ) then
         pot_l => Potential(id)%potL_Tilda(1:iend,1:jmax,ia)
      else if ( nocaseCompare(potType,"Madelung") ) then
         pot_l => Potential(id)%potL_Madelung(1:iend,1:jmax)
      else if ( nocaseCompare(potType,"Total") ) then
         pot_l => Potential(id)%potL(1:iend,1:jmax,spin,ia)
      else if ( nocaseCompare(potType,"Pseudo") ) then
         pot_l => Potential(id)%potL_Pseudo(1:iend,1:jmax)
      else if ( nocaseCompare(potType,"Exchg") ) then
         pot_l => Potential(id)%potL_Exch(1:iend,1:jmax,spin,ia)
      else if ( nocaseCompare(potType,"Coulomb") ) then
         pot_l => Potential(id)%potL_Coulomb(1:iend,1:jmax,ia)
      else if ( nocaseCompare(potType,"XCHat") ) then
         pot_l => Potential(id)%potL_XCHat(1:iend,1:jmax,spin,ia)
      else if ( nocaseCompare(potType,"En_Exchg") ) then
         pot_l => Potential(id)%enL_Exch(1:iend,1:jmax,spin,ia)
      else if ( nocaseCompare(potType,"En_XCHat") ) then
         pot_l => Potential(id)%enL_XCHat(1:iend,1:jmax,spin,ia)
      else
         call ErrorHandler('getPotentialAtPoint','invalid potential type',potType)
!        pot_l => Potential(id)%potL(1:iend,1:jmax,spin)
      endif
!
      do jl = jmax,1,-1
         l = lofj(jl)
!        -------------------------------------------------------------
         call PolyInterp(n_inter, r_mesh(irp:irp+n_inter-1),       &
                         pot_l(irp:irp+n_inter-1,jl), r, pot_in, err)
!        -------------------------------------------------------------
         kl = (l+1)*(l+1)-l+mofj(jl)
         if (mofj(jl) == 0) then
            pot = pot + real(pot_in*Ylm(kl),RealKind)
         else
            pot = pot + TWO*real(pot_in*Ylm(kl),RealKind)
         endif
      enddo
   endif
!
   if (present(grad)) then
      grad = ZERO
   endif
!
   end function getPotentialAtPoint
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calRadialInterpolation(p_fft_c,kpow,ninterp,            &
                                     iparam_out,w_interp)
!  ===================================================================
   use BesselModule, only : SphericalBessel
!
   use SphericalHarmonicsModule, only : calYlmConjg
!
   use ChebyshevModule, only : calChebGrid
!
   use ParallelFFTModule, only : getNumLocalGridPoints,      &
                                 getNumGlobalGridPoints,     &
                                 getLocalIndexOfGridOrigin,  &
                                 getGridPointCoord, getParaFFTCommunicator
!
   use MPPModule, only : GlobalSum, setCommunicator, resetCommunicator
   use GroupCommModule, only : getGroupID, GlobalSumInGroup
!
   use SystemModule, only : getAtomPosition
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kpow, ninterp
   integer (kind=IntKind) :: n_rmesh, idk, idk0, nnr, comm, p, n, nmax, aGID
   integer (kind=IntKind) :: lmax, jmax, kl, jl, l, m, i, nr
   integer (kind=IntKind) :: ir_lsq, na, ia, id, ig
!
   real (kind=RealKind) :: dummy, k2, rmt, kfact_r, kfact_i, r0
   real (kind=RealKind) :: posi(3), kvec(3)
   real (kind=RealKind), pointer :: r_interp(:)
   real (kind=RealKind), target :: r_interp_global(n_interp_max*nr_int_max,GlobalNumAtoms)
   real (kind=RealKind) :: dr_int(nr_int_max)
   real (kind=RealKind) :: kri(n_interp_max*nr_int_max)
   real (kind=RealKind), allocatable, target :: Bj_l(:)
   real (kind=RealKind), pointer :: rmesh(:)
   real (kind=RealKind), pointer ::  p_Bj_l(:,:)
!
   complex (kind=CmplxKind), intent(in) :: p_fft_c(:)
   complex(kind=CmplxKind) :: expikr, kfact, i2l
   complex(kind=CmplxKind), intent(out), target :: w_interp(:,:)
   complex(kind=CmplxKind), target, allocatable :: w_interp_global(:,:)
   complex(kind=CmplxKind), pointer :: pv_interp(:,:), p1(:)
!
!  ===================================================================
!  Interpolation parameters
!  ===================================================================
   integer (kind=IntKind), intent(out) :: iparam_out(:,:)
   integer (kind=IntKind) :: n_interp, n_interp_in, nr_int
   integer (kind=IntKind) :: d_ir(nr_int_max)
   integer (kind=IntKind) :: iparam_global(4+nr_int_max,GlobalNumAtoms)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
!  ===================================================================
!  Back FFT to compute the radial jl components of potential
!  ===================================================================
   if (kpow < 0) then
      idk0 = getLocalIndexOfGridOrigin('K')
      if ( idk0 > 1 .or. idk0 < 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calRadialInterpolation',                  &
                           'Grid origin does not correspond to the first local index',idk0)
!        -------------------------------------------------------------
      endif
   else
      idk0 = 0
   endif
!
   n_interp  = ninterp
   r_interp_global = ZERO
   do id = 1, LocalNumAtoms
      lmax   = Potential(id)%lmax
      n_rmesh = Potential(id)%jend
      rmesh => Potential(id)%Grid%r_mesh
      rmt = rmesh(Potential(id)%jmt)
!
      ig = GlobalIndex(id)
      r_interp => r_interp_global(:,ig)
!
!     ================================================================
!     With increase of l_max, Gibs fenomena occurs at larger r.
!     For different values of l_max define ranges to apply the Chebyshev fit,
!     and use a linear square fit( order 3 polynomial )
!     to extrapolate the Chebyshev fit to small r
!     ================================================================
      if (lmax > 8) then
         nr_int = 3
         dr_int(1) = rmesh(1)
         dr_int(2) = HALF*rmt
         dr_int(3) = rmt
      else if (lmax > 4) then
         nr_int = 2
         dr_int(1) = rmesh(1)
         dr_int(2) = HALF*rmt
      else
         nr_int = 1
         dr_int(1) = rmesh(1)
      endif
!
      do i = 1, nr_int
         r0 = dr_int(i)
!        -------------------------------------------------------------
         call hunt(n_rmesh,rmesh,r0,ir_lsq)
!        -------------------------------------------------------------
         d_ir(i) = ir_lsq
         dr_int(i) = rmesh(ir_lsq)
!        -------------------------------------------------------------
         call calChebGrid( rmesh(ir_lsq), rmesh(n_rmesh), n_interp,   &
                           r_interp((i-1)*n_interp+1:i*n_interp) )
!        -------------------------------------------------------------
         if ( n_rmesh-ir_lsq+1 <= 10 ) then
            n_interp_in = n_rmesh-ir_lsq+1
         else
            n_interp_in = 10
         endif
      enddo
      iparam_out(1,id) = lmax
      iparam_out(2,id) = n_interp
      iparam_out(3,id) = n_interp_in
      iparam_out(4,id) = nr_int
      iparam_out(5:4+nr_int_max,id) = d_ir(1:nr_int_max)
   enddo
!
   comm = getParaFFTCommunicator(MyProc=p,NumProcs=n)
!
   if (comm > 0) then
      na = GlobalNumAtoms
      iparam_global = 0
      do id = 1, LocalNumAtoms
         ig = GlobalIndex(id)
         iparam_global(:,ig) = iparam_out(:,id)
      enddo
!     ================================================================
!     Receive the data stored in iparam_out for the atoms distributed
!     on other processors.
!     ================================================================
      aGID = getGroupID('Unit Cell')
!     ----------------------------------------------------------------
      call GlobalSumInGroup(aGID,iparam_global,                       &
                            4+nr_int_max,GlobalNumAtoms)
      call GlobalSumInGroup(aGID,r_interp_global,                     &
                            n_interp_max*nr_int_max,GlobalNumAtoms)
!     ----------------------------------------------------------------
      lmax = 0
      do ig = 1, GlobalNumAtoms
         lmax = max(lmax,iparam_global(1,ig))
      enddo
      jmax = (lmax+1)*(lmax+2)/2
      nmax = n_interp_max*nr_int_max*jmax
      allocate(w_interp_global(nmax,GlobalNumAtoms))
      w_interp_global = CZERO
   else ! For serial FFT, comm = -1.
      na = LocalNumAtoms
      nmax = 0
      lmax = 0
      do id = 1, LocalNumAtoms
         lmax = max(lmax,iparam_out(1,id))
      enddo
      w_interp = CZERO
   endif
   allocate(Bj_l(n_interp_max*nr_int_max*(lmax+1)))
!
   if ( kpow == -2) then
!     dummy = TWO*PI4*PI4/real(ng,kind=RealKind)
      dummy = TWO*PI4*PI4 ! parallel FFT implementatin already include 1/ng
   else if (kpow == 0) then
!     dummy = PI4/real(ng,kind=RealKind)
      dummy = PI4 ! parallel FFT implementatin already include 1/ng
   else
!     ----------------------------------------------------------------
      call ErrorHandler('calRadialInterpolation','Unconventional situation', &
                        'kpow <> -2 && kpow <> 0. Needs checking!')
!     ----------------------------------------------------------------
   endif
!
   do ia = 1, na
      if (comm > 0) then
!        =============================================================
!        In this case, pv_interp is to be calculated for all the atoms
!        on the local processor
!        =============================================================
         posi = getAtomPosition(ia)
         lmax = iparam_global(1,ia)
         n_interp = iparam_global(2,ia)
         n_interp_in = iparam_global(3,ia)
         nr_int = iparam_global(4,ia)
         d_ir(1:nr_int_max) = iparam_global(5:4+nr_int_max,ia)
!
         r_interp => r_interp_global(:,ia)
!
         nnr = n_interp*nr_int
         jmax = (lmax+1)*(lmax+2)/2
         p1 => w_interp_global(:,ia)
         pv_interp  => aliasArray2_c(p1,nnr,jmax)
      else ! For serial FFT, comm = -1.
!        =============================================================
!        In this case, pv_interp is to be calculated for the local atoms
!        on the local processor
!        =============================================================
         posi = LocalAtomPosi(:,ia)
         lmax = iparam_out(1,ia)
         n_interp = iparam_out(2,ia)
         n_interp_in = iparam_out(3,ia)
         nr_int = iparam_out(4,ia)
         d_ir(1:nr_int_max) = iparam_out(5:4+nr_int_max,ia)
!
         ig = GlobalIndex(ia)
         r_interp => r_interp_global(:,ig)
!
         nnr = n_interp*nr_int
         jmax = (lmax+1)*(lmax+2)/2
         p1 => w_interp(:,ia)
         pv_interp  => aliasArray2_c(p1,nnr,jmax)
      endif
      p_Bj_l(1:,0:) => aliasArray2_r(Bj_l, nnr, lmax+1)
!
      do idk = idk0+1, numk_local
         kvec = getGridPointCoord('K',idk)
         k2 = kvec(1)**2 + kvec(2)**2 + kvec(3)**2
         expikR  = exp(  ( kvec(1)*posi(1) + kvec(2)*posi(2) + kvec(3)*posi(3) )*sqrtm1 )
         do nr = 1,nnr
            kri(nr) = sqrt(k2)*r_interp(nr)
         enddo
!
         if (kpow == -2) then
            kfact = (expikR*p_fft_c(idk))/k2
         else if (kpow == 0) then
            kfact = expikR*p_fft_c(idk)
         else 
            kfact = (expikR*p_fft_c(idk))*sqrt(k2)**kpow
         endif
         kfact_r = real(kfact,kind=RealKind)
         kfact_i = -real(sqrtm1*kfact,kind=RealKind)
!
!        -------------------------------------------------------------
         call calYlmConjg( kvec, lmax, Ylm )
         call SphericalBessel( lmax, nnr, kri, p_Bj_l(1:nnr,0:lmax) )
!        -------------------------------------------------------------
!
         do l = 0,lmax
            if (mod(l,2) == 0) then
               do m = 0,l
                  kl = (l+1)*(l+1) -l + m
                  jl = ((l+1)*(l+2))/2 -l + m
                  pv_interp(1:nnr,jl) = pv_interp(1:nnr,jl) + kfact_r*Ylm(kl)*p_Bj_l(1:nnr,l)
               enddo
            else
               do m = 0,l
                  kl = (l+1)*(l+1) -l + m
                  jl = ((l+1)*(l+2))/2 -l + m
                  pv_interp(1:nnr,jl) = pv_interp(1:nnr,jl) + sqrtm1*kfact_i*Ylm(kl)*p_Bj_l(1:nnr,l)
               enddo
            endif
         enddo
      enddo
!
      i2l = dummy
      do l = 0,lmax
         do m = 0,l
            jl = ((l+1)*(l+2))/2 - l + m
            do i =1,nnr
               pv_interp(i,jl) = i2l*pv_interp(i,jl)/(r_interp(i)**l)
            enddo
         enddo
!        i2l = i2l*(-sqrtm1)  ! Mar 27, 2017 - Yang
         i2l = i2l*sqrtm1
      enddo
   enddo
   nullify(rmesh, r_interp, p_Bj_l, pv_interp)
!
!  ===================================================================
!  In parallel FFT case, a summation over all k-points, which are
!  distributed across pocessors, needs to be performed.
!  ===================================================================
   if (comm > 0) then
!     ================================================================
!     In this case, w_interp_global is summed across all processors, 
!     which, in effect, sums over the reciprocal uniform grid points.
!     ----------------------------------------------------------------
      call setCommunicator(comm,p,n)
      call GlobalSum(w_interp_global,nmax,GlobalNumAtoms)
      call resetCommunicator()
!     ----------------------------------------------------------------
      do id = 1, LocalNumAtoms
         ig = GlobalIndex(id)
!        -------------------------------------------------------------
         call zcopy(nmax,w_interp_global(1,ig),1,w_interp(1,id),1)
!        -------------------------------------------------------------
      enddo
      deallocate(w_interp_global)
   endif
   deallocate(Bj_l)
!
   end subroutine calRadialInterpolation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calRadialProjection(n_rmesh,rmesh,iparam_in,w_interp,pv_jl)
!  ===================================================================
   use IntegerFactorsModule, only : lofj
!
   use ChebyshevModule, only : initChebyshevSeries, ChebyshevEval,    &
                               endChebyshevSeries
!
   use InterpolationModule, only : LeastSqFitInterp, FitInterp
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n_rmesh
   integer (kind=IntKind), intent(in) :: iparam_in(:)
!
   integer (kind=IntKind) :: nnr
   integer (kind=IntKind) :: lmax, jmax, jl, l, i
   integer (kind=IntKind) :: ir_lsq, ir_lsqm1
!
   real (kind=RealKind), intent(in) :: rmesh(n_rmesh)
   real (kind=RealKind) :: r0
!
   complex(kind=CmplxKind), intent(in) :: w_interp(:)
   complex(kind=CmplxKind), intent(out), target :: pv_jl(:,:)
   complex(kind=CmplxKind), pointer :: pv_interp(:,:)
   complex(kind=CmplxKind), pointer :: pv(:), pv0(:)
!
!  ===================================================================
!  Interpolation parameters from subroutine calRadialInterpolation
!  ===================================================================
   integer (kind=IntKind) :: n_interp, n_interp_in, nr_int
   integer (kind=IntKind) :: d_ir(nr_int_max)
!
   lmax = iparam_in(1)
   n_interp = iparam_in(2)
   n_interp_in = iparam_in(3)
   nr_int = iparam_in(4)
   d_ir(1:nr_int_max) = iparam_in(5:4+nr_int_max)
!
!  ===================================================================
!  With increase of l_max, Gibs fenomena occurs at larger r.
!  For different values of l_max define ranges to apply the Chebyshev fit,
!  and use a linear square fit( order 3 polynomial )
!  to extrapolate the Chebyshev fit to small r
!  ===================================================================
   jmax = (lmax+1)*(lmax+2)/2
   nnr = n_interp*nr_int
   pv_interp  => aliasArray2_c(w_interp,nnr,jmax)
!
   pv_jl = CZERO
   do jl = 1, jmax
      l = lofj(jl)
      if (l <= 4) then
         pv => pv_interp(1:n_interp,jl)
         r0 = rmesh(d_ir(1))
         ir_lsq = d_ir(1)
      else if (l > 4 .and. l <= 8) then
         pv => pv_interp(n_interp+1:2*n_interp,jl)
         r0 = rmesh(d_ir(2))
         ir_lsq = d_ir(2)
      else
         pv => pv_interp(2*n_interp+1:3*n_interp,jl)
         r0 = rmesh(d_ir(3))
         ir_lsq = d_ir(3)
      endif
!     ----------------------------------------------------------------
      call initChebyshevSeries( n_interp, r0, rmesh(n_rmesh), &
                               pv, chebv_struct(jl) )
!     ----------------------------------------------------------------
      pv => pv_jl(ir_lsq:n_rmesh,jl)
!     ----------------------------------------------------------------
      call ChebyshevEval( chebv_struct(jl), n_interp, n_rmesh-ir_lsq+1, &
                          rmesh(ir_lsq:n_rmesh), pv )
      call endChebyshevSeries(chebv_struct(jl))
!     ----------------------------------------------------------------
   enddo
!
   if (nr_int > 1) then
      do jl = 1,jmax
         l = lofj(jl)
         if ( l<=4 ) then
            ir_lsq = d_ir(1)
            ir_lsqm1 = ir_lsq - 1
         else if (l > 4 .and. l <= 8) then
            ir_lsq = d_ir(2)
            ir_lsqm1 = ir_lsq - 1
         else
            ir_lsq = d_ir(3)
            ir_lsqm1 = ir_lsq - 1
         endif
         pv => pv_jl(1:ir_lsqm1,jl)
         pv0 => pv_jl(ir_lsq:ir_lsq+n_interp_in-1,jl)
!        -------------------------------------------------------------
         call LeastSqFitInterp( n_interp_in, rmesh(ir_lsq:n_rmesh),   &
                                pv0, ir_lsqm1, rmesh(1:ir_lsqm1), pv )
!        -------------------------------------------------------------
         do i = 1, n_rmesh
            pv_jl(i,jl) = pv_jl(i,jl)*(rmesh(i)**l)
         enddo
      enddo
   else
      do jl = 1,jmax
         l = lofj(jl)
         do i = 1, n_rmesh
            pv_jl(i,jl) = pv_jl(i,jl)*(rmesh(i)**l)
         enddo
      enddo
   endif
!
   nullify( pv_interp, pv, pv0 )
!
   end subroutine calRadialProjection
!  ===================================================================
end Module PotentialGenerationModule
