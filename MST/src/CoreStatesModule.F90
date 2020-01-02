module CoreStatesModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
   use MathParamModule, only : ZERO, HALF, ONE, TWO, THREE, FOUR, THIRD
   use MathParamModule, only : TEN, TEN2m5, TEN2m6, TEN2m8, TEN2m10, SQRTm1, &
                               TEN2p10,PI4,TEN2m12, THIRD, PI4
   use PhysParamModule, only : LightSpeed
   use IntegrationModule, only : calIntegration
   use DerivativeModule, only : derv5
   use PublicTypeDefinitionsModule, only : GridStruct
!
public :: initCoreStates,     &
          endCoreStates,      &
          calCoreStates,      &
          readCoreStates,     &
          readCoreDensity,    &
          printCoreStates,    &
          printCoreDensity,   &
          getDeepCoreDensity, &
          getSemiCoreDensity, &
          getDeepCoreDensityDerivative,&
          getSemiCoreDensityDerivative,&
          getDeepCoreEnergy,  &
          getDeepCoreKineticEnergy,  &
          getSemiCoreEnergy,  &
          getSemiCoreKineticEnergy,  &
          getCoreVPCharge,    &
          getCoreMTCharge,    &
          getCoreVPMoment,    &
          getCoreMTMoment,    &
          getInterstitialSemiCoreDensity,  &
          getInterstitialDeepCoreDensity,  &
          writeCoreStates,    &
          writeCoreDensity,   &
          getCoreSplitTable,     &
          getCoreNumStatesTable, &
          getCoreDescriptionTable
!
   interface getDeepCoreEnergy
      module procedure getDeepCoreE0, getDeepCoreE1
   end interface
!
   interface getDeepCoreDensity
      module procedure getDCD0, getDCD1, getDCD2
   end interface
!
   interface getSemiCoreEnergy
      module procedure getSemiCoreE0, getSemiCoreE1
   end interface
!
   interface getSemiCoreDensity
      module procedure getSCD0, getSCD1, getSCD2
   end interface
!
   interface getDeepCoreDensityDerivative
      module procedure getDCDDer0, getDCDDer1, getDCDDer2
   end interface
!
   interface getSemiCoreDensityDerivative
      module procedure getSCDDer0, getSCDDer1, getSCDDer2
   end interface
!
   interface getCoreVPCharge
      module procedure getCoreVPC0, getCoreVPC1
   end interface
!
   interface getCoreMTCharge
      module procedure getCoreMTC0, getCoreMTC1
   end interface
!
   interface getCoreVPMoment
      module procedure getCoreVPM0, getCoreVPM1
   end interface
!
   interface getCoreMTMoment
      module procedure getCoreMTM0, getCoreMTM1
   end interface
!
private
   type CoreStruct
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind) :: MaxNumc
      integer (kind=IntKind) :: rsize
      integer (kind=IntKind) :: jcore
      integer (kind=IntKind), pointer :: numc_below(:)
      integer (kind=IntKind), pointer :: numc(:)
      integer (kind=IntKind), pointer :: nc(:,:)
      integer (kind=IntKind), pointer :: lc(:,:)
      integer (kind=IntKind), pointer :: kc(:,:)
      real (kind=RealKind), pointer :: ec(:,:,:)
      real (kind=RealKind), pointer :: ecorv(:,:)
      real (kind=RealKind), pointer :: esemv(:,:)
      real (kind=RealKind), pointer :: corden(:,:,:) ! deep core density, without r^2
      real (kind=RealKind), pointer :: semden(:,:,:) ! semi core density, without r^2
      real (kind=RealKind), pointer :: dcorden(:,:,:) ! derivative of deep core density, without r^2
      real (kind=RealKind), pointer :: dsemden(:,:,:) ! derivative of semi core density, without r^2
      real (kind=RealKind), pointer :: OldCore(:,:,:)
      real (kind=RealKind), pointer :: r_mesh(:)
      real (kind=RealKind), pointer :: core_ke(:,:)
      real (kind=RealKind), pointer :: semi_ke(:,:)
      real (kind=RealKind) :: vol_core
      real (kind=RealKind) :: vol_coreint
      real (kind=RealKind) :: rcore_mt
      real (kind=RealKind), pointer :: ztotss(:)
      real (kind=RealKind), pointer :: zsemss(:)
      real (kind=RealKind), pointer :: zcorss(:)
      real (kind=RealKind), pointer :: qcpsc_mt(:)
      real (kind=RealKind), pointer :: qcpsc_ws(:)
      real (kind=RealKind), pointer :: mcpsc_mt(:)
      real (kind=RealKind), pointer :: mcpsc_ws(:)
      real (kind=RealKind), pointer :: qsemmt(:)
      real (kind=RealKind), pointer :: qcormt(:)
      real (kind=RealKind), pointer :: qsemws(:)
      real (kind=RealKind), pointer :: qcorws(:)
      real (kind=RealKind), pointer :: qcorout(:)
      type (GridStruct), pointer :: Grid
   end type CoreStruct
   type (CoreStruct), allocatable, target :: Core(:)
!
   logical :: Initialized = .false.
   logical :: isNonRelativistic
!
   character (len = 50) :: stop_routine
!
   integer (kind=IntKind), parameter :: integration_order = 3
   integer (kind=IntKind), parameter :: ipdeq = 5
   integer (kind=IntKind), parameter :: ipdeq2 = 2*ipdeq
   integer (kind=IntKind), parameter :: nitmax = 50
   integer (kind=IntKind), parameter :: LargeRs = 1101
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind), parameter :: n_inter = 5 ! order of polynomial
                                                    ! interpolation
!
   real (kind=RealKind), parameter :: tol = TEN2m12
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: ndeepz
   integer (kind=IntKind) :: GroupID
   integer (kind=IntKind) :: NumPEsInGroup
   integer (kind=IntKind) :: MyPEinGroup
   integer (kind=IntKind) :: LargeRsMult
   integer (kind=IntKind) :: lmax_core = 0
!
   real (kind=RealKind) :: MyLightSpeed
   real (kind=RealKind) :: cinv
   real (kind=RealKind) :: c2inv
   real (kind=RealKind) :: evbot
   real (kind=RealKind) :: etopcor
   real (kind=RealKind) :: rhoint_deep
   real (kind=RealKind) :: rhoint_semi
   real (kind=RealKind) :: TotalInterstitialCoreVol
   real (kind=RealKind), allocatable :: tmp(:)
   real (kind=RealKind), allocatable :: qmp(:)
   real (kind=RealKind), allocatable :: ftmp(:)
   real (kind=RealKind), allocatable :: dftmp(:)
!
   logical :: GlobalTable_Allocated = .false.
   integer (kind=IntKind) :: maxnc_save = 0
   character (len=2), allocatable, target :: DescriptionTable(:,:)
   integer (kind=IntKind), allocatable :: TmpDesTable(:,:)
   integer (kind=IntKind), allocatable, target :: NumStatesTable(:)
   real (kind=RealKind), allocatable, target :: SplitTable(:,:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initCoreStates(na,evb,ns,isNonRel,istop,iprint)
!  ===================================================================
   use GroupCommModule, only : getMyPEinGroup, getGroupID, getNumPEsInGroup
!
   use SystemModule, only : getNumAtoms, getNumAlloyElements
   use RadialGridModule, only : getGrid
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       setDataStorage2Value,  &
                                       setDataStorageLDA,     &
                                       isDataStorageExisting, &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark
   use AtomModule, only : getAtomCoreRad, setAtomCoreRad
   use PolyhedraModule, only : getVolume
!
   use GroupCommModule, only : GlobalSumInGroup
!
   use Atom2ProcModule, only : getLocalIndex, getAtom2ProcInGroup
   use Atom2ProcModule, only : getMaxLocalNumAtoms, getGlobalIndex
!
   implicit none
!
   logical, intent(in) :: isNonRel
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind), intent(in) :: ns
   integer (kind=IntKind) :: id, jcore, jend, ig, nlrs, ntmp, last
   integer (kind=IntKind) :: g2last
   integer (kind=IntKind), allocatable :: DataSize(:)
!
   real (kind=RealKind), intent(in) :: evb
!
   real (kind=RealKind), pointer :: r_mesh(:)
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
   isNonRelativistic = isNonRel
!
   LocalNumAtoms = na
   n_spin_pola = ns
   if (isNonRelativistic) then
      MyLightSpeed=LightSpeed*TEN2p10
      cinv=ONE/MyLightSpeed
!     cinv = ZERO  ! 02/11/18 -Yang
   else
      MyLightSpeed=LightSpeed
      cinv=ONE/MyLightSpeed
   endif
   c2inv=cinv*cinv
   allocate( Core(LocalNumAtoms), DataSize(LocalNumAtoms) )
!
   stop_routine = istop
   print_level = iprint
!
   TotalInterstitialCoreVol = ZERO
!
   do id = 1,LocalNumAtoms
      ig = getGlobalIndex(id)
      Core(id)%Grid => getGrid(id)
      if ( Core(id)%Grid%nmult>1 ) then
         ntmp = (Core(id)%Grid%jend_plus_n-Core(id)%Grid%jinsc)/Core(id)%Grid%nmult
         if ( LargeRs-Core(id)%Grid%jinsc>ntmp) then
            nlrs = Core(id)%Grid%jinsc+(LargeRs-Core(id)%Grid%jinsc)*Core(id)%Grid%nmult
            LargeRsMult = nlrs
         else
            LargeRsMult = Core(id)%Grid%jend_plus_n
         endif
      else
         LargeRsMult = LargeRs
      endif
      last=max(Core(id)%Grid%jend_plus_n, LargeRsMult)
      r_mesh => Core(id)%Grid%r_mesh
      Core(id)%NumSpecies = getNumAlloyElements(ig)
      Core(id)%rsize = last
      jcore = Core(id)%Grid%jmt
      if ( getAtomCoreRad(id) < 0.10d0 ) then
         Core(id)%rcore_mt = Core(id)%Grid%r_mesh(jcore)
         Core(id)%jcore    = jcore
         call setAtomCoreRad(id,Core(id)%rcore_mt)
      else
         Core(id)%rcore_mt = getAtomCoreRad(id)
!        Find the first point on the radial grid that is smaller
!        than the input rcore_mt
         jend = Core(id)%Grid%jend
         call hunt(jend,r_mesh,Core(id)%rcore_mt,jcore)
         if (r_mesh(jcore)>Core(id)%rcore_mt) then
            Core(id)%rcore_mt = r_mesh(jcore-1)
            Core(id)%jcore = jcore-1
         else
            Core(id)%rcore_mt = r_mesh(jcore)
            Core(id)%jcore = jcore
         endif
!
      endif
      Core(id)%vol_core = PI4*THIRD*(Core(id)%rcore_mt**3)
      Core(id)%vol_coreint = getVolume(id)-Core(id)%vol_core
      DataSize(id)  = last*n_spin_pola*Core(id)%NumSpecies
      TotalInterstitialCoreVol = TotalInterstitialCoreVol + &
                                 Core(id)%vol_coreint
!
      allocate( Core(id)%numc_below(Core(id)%NumSpecies),            &
                Core(id)%numc(Core(id)%NumSpecies),                  &
                Core(id)%ztotss(Core(id)%NumSpecies),                &
                Core(id)%zsemss(Core(id)%NumSpecies),                &
                Core(id)%zcorss(Core(id)%NumSpecies),                &
                Core(id)%qcpsc_mt(Core(id)%NumSpecies),              &
                Core(id)%qcpsc_ws(Core(id)%NumSpecies),              &
                Core(id)%mcpsc_mt(Core(id)%NumSpecies),              &
                Core(id)%mcpsc_ws(Core(id)%NumSpecies),              &
                Core(id)%qsemmt(Core(id)%NumSpecies),                &
                Core(id)%qcormt(Core(id)%NumSpecies),                &
                Core(id)%qsemws(Core(id)%NumSpecies),                &
                Core(id)%qcorws(Core(id)%NumSpecies),                &
                Core(id)%qcorout(Core(id)%NumSpecies) )
   enddo
!
!  ==================================================================
!  generate Storage space for charge densities
!  ==================================================================
   if (.not.isDataStorageExisting('CoreDensity')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'CoreDensity',  &
                             DataSize, RealType)
      call createDataStorage(LocalNumAtoms,'DeriCoreDensity',  &
                             DataSize, RealType)
!     ----------------------------------------------------------------
      call setDataStorage2Value('CoreDensity',ZERO)
      call setDataStorage2Value('DeriCoreDensity',ZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         last = Core(id)%rsize
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'CoreDensity',last)
         call setDataStorageLDA(id,'DeriCoreDensity',last)
!        -------------------------------------------------------------
      enddo
   endif
   if (.not.isDataStorageExisting('OldCoreDensity')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'OldCoreDensity',  &
                             DataSize, RealType)
!     ----------------------------------------------------------------
      call setDataStorage2Value('OldCoreDensity',ZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         last = Core(id)%rsize
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'OldCoreDensity',last)
!        -------------------------------------------------------------
      enddo
   endif
   if (.not.isDataStorageExisting('SemiCoreDensity')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'SemiCoreDensity',         &
                             DataSize, RealType)
      call createDataStorage(LocalNumAtoms,'DeriSemiCoreDensity',     &
                             DataSize, RealType)
!     ----------------------------------------------------------------
      call setDataStorage2Value('SemiCoreDensity',ZERO)
      call setDataStorage2Value('DeriSemiCoreDensity',ZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         last = Core(id)%rsize
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'SemiCoreDensity',last)
         call setDataStorageLDA(id,'DeriSemiCoreDensity',last)
!        -------------------------------------------------------------
      enddo
   endif
!
   evbot = evb
   g2last = 0
   do id =1,LocalNumAtoms
      last = Core(id)%rsize
      Core(id)%corden => getDataStorage( id, 'CoreDensity', &
                                         last, n_spin_pola, Core(id)%NumSpecies, RealMark )
      Core(id)%semden => getDataStorage( id, 'SemiCoreDensity', &
                                         last, n_spin_pola, Core(id)%NumSpecies, RealMark )
      Core(id)%dcorden => getDataStorage( id, 'DeriCoreDensity', &
                                          last, n_spin_pola, Core(id)%NumSpecies, RealMark )
      Core(id)%dsemden => getDataStorage( id, 'DeriSemiCoreDensity', &
                                          last, n_spin_pola, Core(id)%NumSpecies, RealMark )
      Core(id)%OldCore => getDataStorage( id, 'OldCoreDensity', &
                                         last, n_spin_pola, Core(id)%NumSpecies, RealMark )
!
      allocate( Core(id)%r_mesh(1:last) )
      allocate( Core(id)%ecorv(1:n_spin_pola,1:Core(id)%NumSpecies) )
      allocate( Core(id)%esemv(1:n_spin_pola,1:Core(id)%NumSpecies) )
      allocate( Core(id)%core_ke(1:n_spin_pola,1:Core(id)%NumSpecies) )
      allocate( Core(id)%semi_ke(1:n_spin_pola,1:Core(id)%NumSpecies) )
      g2last=max( last, g2last )
   enddo
!
   GlobalNumAtoms = getNumAtoms()
   allocate( tmp(0:g2last), qmp(0:g2last), ftmp(0:g2last), dftmp(1:g2last) )
!
   deallocate( DataSize )
!
   rhoint_semi = ZERO
   rhoint_deep = ZERO
   Initialized = .true.
   GlobalTable_Allocated = .false.
   maxnc_save = 0
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   call GlobalSumInGroup(GroupID,TotalInterstitialCoreVol)
!
   end subroutine initCoreStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readCoreStates()
!  ===================================================================
   use ChemElementModule, only : MaxLenOfAtomName
   use ChemElementModule, only : getZtot, getZcor, getZsem
   use ChemElementModule, only : getNumCoreStates
   use ChemElementModule, only : getCoreStateN
   use ChemElementModule, only : getCoreStateL
   use ChemElementModule, only : getCoreStateKappa
!
   use AtomModule, only : getInPotFileName, getInPotFileForm,         &
                          getLocalAtomNickName
!
   implicit none
!
   character (len=MaxLenOfAtomName) :: an
   character (len=10) :: acc
!
   integer (kind=IntKind) :: id, ia
   integer (kind=IntKind) :: MaxNumc, ic
   integer (kind=IntKind), parameter :: funit=92
!
   do id = 1, LocalNumAtoms
      MaxNumc = 0
      do ia = 1, Core(id)%NumSpecies
         an = getLocalAtomNickName(id,ia)
         MaxNumc = max(MaxNumc,getNumCoreStates(an))
      enddo
      Core(id)%MaxNumc = MaxNumc
      allocate( Core(id)%nc(1:MaxNumc,Core(id)%NumSpecies) )
      allocate( Core(id)%lc(1:MaxNumc,Core(id)%NumSpecies) )
      allocate( Core(id)%kc(1:MaxNumc,Core(id)%NumSpecies) )
      allocate( Core(id)%ec(1:MaxNumc,1:n_spin_pola,Core(id)%NumSpecies) )
      Core(id)%nc = 0
      Core(id)%lc = 0
      Core(id)%kc = 0
      Core(id)%ec = ZERO
!
      if (getInPotFileForm(id) == 'FORMATTED') then
         acc = 'SEQUENTIAL'
      else
         acc = 'DIRECT'
      endif
      do ia = 1, Core(id)%NumSpecies
         if (getInPotFileForm(id) == 'FORMATTED') then
!           ----------------------------------------------------------
            open(unit=funit,file=getInPotFileName(id,ia),form='FORMATTED',  &
                 access=acc)
!           ----------------------------------------------------------
            call readFormattedData(funit,id,ia)
!           ----------------------------------------------------------
            close(unit=funit)
!           ----------------------------------------------------------
         else
!           ----------------------------------------------------------
            call readUnformattedData(funit,id,ia)
!           ----------------------------------------------------------
         endif
      enddo
!
      do ia = 1, Core(id)%NumSpecies
         an = getLocalAtomNickName(id,ia)
         Core(id)%ztotss(ia)=getZtot(an)
         Core(id)%zsemss(ia)=getZsem(an)
         Core(id)%zcorss(ia)=getZcor(an)
         Core(id)%numc(ia)=getNumCoreStates(an)
      enddo
!
!     do ia = 1, Core(id)%NumSpecies
!        an = getLocalAtomNickName(id,ia)
!        do ic=1, Core(id)%numc(ia)
!           Core(id)%nc(ic,ia)=getCoreStateN(an,ic)
!           Core(id)%lc(ic,ia)=getCoreStateL(an,ic)
!           Core(id)%kc(ic,ia)=getCoreStateKappa(an,ic)
!        enddo
!     enddo
   enddo
!
   end subroutine readCoreStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readFormattedData(funit,id,ia)
!  ===================================================================
   use InterpolationModule, only : PolyInterp
!
   use ChemElementModule, only : MaxLenOfAtomName
   use ChemElementModule, only : getZtot, getZcor
   use ChemElementModule, only : getNumCoreStates, setNumCoreStates
!
   use PublicParamDefinitionsModule, only : MaxLenFileName
!
   use MPPModule, only : MyPE
!
   use AtomModule, only : getLocalAtomNickName
!
   implicit none
!
   logical :: nmd
!
   character (len=1) :: dummy
   character (len=MaxLenFileName) :: inp
   character (len=MaxLenOfAtomName) :: an
   character (len=17), parameter :: sname='readFormattedData'
!
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: is
   integer (kind=IntKind) :: j_inter, irp
   integer (kind=IntKind) :: ns,j,jmt,nrv,nrrho,nrcor
   integer (kind=IntKind) :: numc
   integer (kind=IntKind) :: nc, lc, kc
!
   real (kind=RealKind) :: ztss,alat,zcss,edum
   real (kind=RealKind) :: xstart,xmt,err,h, xvalws
   real (kind=RealKind), allocatable :: corden(:), x_mesh(:)
!
   type (CoreStruct), pointer :: ThisC
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
   ThisC=>Core(id)
   an = getLocalAtomNickName(id,ia)
!
   read(funit,'(a)') dummy
   read(funit,*) ns
   if (ns /= n_spin_pola .and. MyPE == 0) then
      inquire(unit=funit,named=nmd,name=inp)
      write(6,'(''File Name: '',a)')trim(inp)
      call WarningHandler(sname,'Spin index in file <> n_spin_pola', &
                          ns,n_spin_pola)
   endif
   do is=1,n_spin_pola
      read(funit,'(a)')dummy
      read(funit,*)ztss,alat,zcss
      if(int(ztss) /= getZtot(an)) then
         call ErrorHandler(sname,'Inconsistent atomic Z',ztss)
      else if (int(zcss) /= getZcor(an)) then
         call WarningHandler(sname,'The deep core charge is different from the database',zcss)
      endif
      read(funit,'(17x,2d20.13,i5)') xstart,xmt,jmt
!
!     ================================================================
!     read in formatted one-electron LDA potential..............
!     ================================================================
      nrv=ceiling(jmt/4.0)
      read(funit,'(a)') (dummy,j=1,nrv)
      read(funit,'(a)') dummy
!
!     ================================================================
!     read in formatted total charge density....................
!     ================================================================
      read(funit,'(i5,e20.13)') nrrho,xvalws
      nrrho=ceiling(nrrho/4.0)       ! data is stored 4 numbers per line
      read(funit,'(a)')(dummy,j=1,nrrho)
!
!     ================================================================
!     read in formatted core state information..................
!     ================================================================
      read(funit,'(2i5)') numc,nrcor
      do j=1,numc
         if (j <= getNumCoreStates(an)) then
            read(funit,'(3i5,f12.5,1x,a)')ThisC%nc(j,ia),ThisC%lc(j,ia), &
                                          ThisC%kc(j,ia),ThisC%ec(j,is,ia),dummy
!           if (ThisC%nc(j,ia) /= nc) then
!              call ErrorHandler(sname,'nc <> ThisC%nc',nc, ThisC%nc(j,ia))
!           else if (ThisC%lc(j,ia) /= lc) then
!              call ErrorHandler(sname,'lc <> ThisC%lc',lc, ThisC%lc(j,ia))
!           else if (ThisC%kc(j,ia) /= kc) then
!              call ErrorHandler(sname,'kc <> ThisC%kc',kc, ThisC%kc(j,ia))
!           endif
         else
            read(funit,'(3i5,f12.5,1x,a)')nc,lc,kc,edum,dummy
         endif
      enddo
!
      if (numc > getNumCoreStates(an)) then
         if (print_level >= 0) then
            call WarningHandler(sname,'numc > getNumCoreStates()',numc,&
                                getNumCoreStates(an))
         endif
         call setNumCoreStates(an,numc)
      else if (numc < getNumCoreStates(an) ) then
         if (print_level >= 0) then
            call WarningHandler(sname,'numc < getNumCoreStates', numc, &
                                getNumCoreStates(an))
         endif
         call setNumCoreStates(an,numc)
      endif
!
      if (nrcor > 0) then
         allocate(corden(nrcor), x_mesh(nrcor))
         read(funit,'(4e20.13)') (corden(j),j=1,nrcor)
!
!        =============================================================
!        generate grid for Core density read-in ......................
!        =============================================================
         h=(xmt-xstart)/dble(jmt-1)
         do j=1,nrcor
            x_mesh(j)=xstart+(j-1)*h
         enddo
!
!        =============================================================
!        interpolate chg density onto GRID mesh
!        =============================================================
         if (ThisC%Grid%jend_plus_n < n_inter) then
!           ----------------------------------------------------------
            call ErrorHandler(sname,'jend_plus_n < n_inter',          &
                              ThisC%Grid%jend_plus_n,n_inter)
!           ----------------------------------------------------------
         endif
         j_inter=1
         do j=1,ThisC%Grid%jend_plus_n
!           ----------------------------------------------------------
            call hunt(nrcor,x_mesh,ThisC%Grid%x_mesh(j),j_inter)
!           ----------------------------------------------------------
            if (j_inter > nrcor-(n_inter-1)/2) then
               irp=nrcor-n_inter+1
            else if (2*j_inter+1 > n_inter) then
               irp=j_inter-(n_inter-1)/2
            else
               irp=1
            endif
!           ----------------------------------------------------------
            call PolyInterp(n_inter,x_mesh(irp:irp+n_inter-1),        &
                            corden(irp:irp+n_inter-1), &
                            ThisC%Grid%x_mesh(j),ThisC%OldCore(j,is,ia),err)
!           ----------------------------------------------------------
         enddo
         deallocate(x_mesh, corden)
      endif
!
      if (ns < n_spin_pola) then
         do j=1,numc
            if (j <= getNumCoreStates(an)) then
               ThisC%ec(j,n_spin_pola,ia) = ThisC%ec(j,1,ia)
            endif
         enddo
         do j=1,ThisC%Grid%jend_plus_n
            ThisC%OldCore(j,1,ia) = HALF*ThisC%OldCore(j,1,ia)
         enddo
         do j=1,ThisC%Grid%jend_plus_n
            ThisC%OldCore(j,n_spin_pola,ia) = ThisC%OldCore(j,1,ia)
         enddo
         exit
      else if (ns > n_spin_pola) then
         do j=1,ThisC%Grid%jend_plus_n
            ThisC%OldCore(j,1,ia) = TWO*ThisC%OldCore(j,1,ia)
         enddo
      endif
   enddo
   end subroutine readFormattedData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readUnformattedData(funit,id,ia)
!  ===================================================================
   use DataServiceCenterModule, only : getDataStorage,        &
                                       isDataStorageExisting, &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark, &
                                       IntegerMark, IntegerType
!
   use AtomModule, only : getLocalAtomNickName
!
   use ChemElementModule, only : MaxLenOfAtomName
   use ChemElementModule, only : getNumCoreStates
   use ChemElementModule, only : getCoreStateN, getCoreStateL,&
                                 getCoreStateKappa, getCoreStateSymbol
!
   implicit none
!
   character (len=MaxLenOfAtomName) :: an
   character (len=19), parameter :: sname='readUnFormattedData'
!
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: id, ia
!
   integer (kind=IntKind) :: ic, is, data_sz, numc, n
   integer (kind=IntKind), pointer :: pc0(:,:,:)
!
   real (kind=RealKind), pointer :: ec0(:,:,:)
!
   if ( isDataStorageExisting('OldEstimatedCoreEnergy') ) then
      an = getLocalAtomNickName(id,ia)
      numc = getNumCoreStates(an)
      if (numc > 0) then
         data_sz = Core(id)%MaxNumc*n_spin_pola*Core(id)%NumSpecies
         ec0 => getDataStorage( id, 'OldEstimatedCoreEnergy',         &
                                Core(id)%MaxNumc, n_spin_pola,        &
                                Core(id)%NumSpecies, RealMark )
         pc0 => getDataStorage( id, 'OldEstimatedCoreStates',         &
                                Core(id)%MaxNumc, n_spin_pola,        &
                                Core(id)%NumSpecies, IntegerMark )
         do is = 1,n_spin_pola
            do ic = 1, numc
               n = pc0(ic,is,ia)
               Core(id)%nc(ic,ia)=getCoreStateN(an,n)
               Core(id)%lc(ic,ia)=getCoreStateL(an,n)
               Core(id)%kc(ic,ia)=getCoreStateKappa(an,n)
               Core(id)%ec(ic,is,ia) = ec0(ic,is,ia)
            enddo
         enddo
      endif
      nullify(ec0)
  else
      call ErrorHandler("readUnformatedDataCore",                     &
        'The data storage for Core Energy have not been created yet!')
   endif
!
   end subroutine readUnformattedData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeCoreStates(id,filename,fileform)
!  ===================================================================
   implicit none
!
   character (len=10) :: acc
   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: fileform
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), parameter :: funit=93
!
   call ErrorHandler('writeCoreStates','Not implemented')
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('writeCoreStates','invalid local atom index',id)
   endif
!
   if (fileform == 'FORMATTED') then
      acc = 'SEQUENTIAL'
   else
      acc = 'DIRECT'
   endif
   if (fileform == 'FORMATTED') then
!     ----------------------------------------------------------------
      open(unit=funit,file=filename,form=fileform,access=acc)
!     ----------------------------------------------------------------
      call writeFormattedData(funit,id)
   else if (fileform == 'UNFORMATTED') then
      call writeUnformattedData(funit,id)
   endif
!
   close(unit=funit)
   end subroutine writeCoreStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFormattedData(funit,id)
!  ===================================================================
   implicit none
   character (len=18), parameter :: sname='writeFormattedData'
!
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: id
!
   call ErrorHandler(sname,'Not implemented')
!
   end subroutine writeFormattedData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeUnformattedData(funit,id)
!  ===================================================================
   implicit none
   character (len=20), parameter :: sname='writeUnformattedData'
!
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: id
!
   call ErrorHandler(sname,'Not implemented')
!
   end subroutine writeUnformattedData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCoreStates()
!  ===================================================================
   integer (kind=IntKind) :: id
!
   do id = 1, LocalNumAtoms
      deallocate( Core(id)%nc, Core(id)%lc, Core(id)%kc, Core(id)%ec )
      deallocate( Core(id)%numc, Core(id)%numc_below )
!      deallocate( Core(id)%corden, Core(id)%semden, Core(id)%OldCore)
      nullify( Core(id)%corden, Core(id)%semden, Core(id)%OldCore)
      nullify( Core(id)%dcorden, Core(id)%dsemden )
      deallocate( Core(id)%ecorv, Core(id)%esemv, Core(id)%r_mesh )
      deallocate( Core(id)%core_ke, Core(id)%semi_ke )
      nullify(Core(id)%Grid)
   enddo
   deallocate(qmp, tmp, Core, ftmp, dftmp)
!
   if (GlobalTable_Allocated) then
      deallocate(NumStatesTable)
      deallocate(TmpDesTable,DescriptionTable)
      deallocate(SplitTable)
   endif
   GlobalTable_Allocated = .false.
   maxnc_save = 0
!
   end subroutine endCoreStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readCoreDensity(fname)
!  ===================================================================
   use InterpolationModule, only : PolyInterp
!
   use MPPMOdule, only : setCommunicator, resetCommunicator,          &
                         sendMessage, recvMessage, syncAllPEs
!
   use ParallelIOModule, only : isInputProc, getMyInputProc,          &
                                getNumInputClients, getInputClient,   &
                                getIOCommunicator, getMyPEinIOGroup,  &
                                getNumPEsInIOGroup
!
   use Atom2ProcModule, only : getGlobalIndex
!
   implicit none
!
   character (len=*), intent(in) :: fname
!
   integer (kind=IntKind) :: integer4_size,real8_size
   integer (kind=IntKind) :: cunit, i, id, ia, ic, ig, is, fp_pos
   integer (kind=IntKind) :: isize, fsize, imsgbuf_size, fmsgbuf_size
   integer (kind=IntKind) :: num_clients, proc_client, present_atom
!
   integer (kind=IntKind) :: MyPEinIOGroup, NumPEsInIOGroup, io_comm
   integer (kind=IntKind) :: msgid1, msgid2
!
!  Assuming that the size of the core density file does not exceed 2**31 bytes.
   integer (kind=IntKind) :: file_loc(GlobalNumAtoms+1) 
   integer (kind=IntKind) :: block_size(2,GlobalNumAtoms)
!
   integer (kind=IntKind), allocatable :: imsgbuf(:)
   integer (kind=IntKind) :: rsize, jcore, jinsc
   integer (kind=IntKind) :: j_inter, irp, jrp, j
!
   real (kind=RealKind), allocatable :: fmsgbuf(:) , x_mesh(:)
   real (kind=RealKind) :: xstart, hin, hout, vol_core, vol_coreint, rcore_mt
   real (kind=RealKind) :: err
!
   logical :: Interp
!
!  -------------------------------------------------------------------
   MyPEinIOGroup = getMyPEinIOGroup()
   NumPEsInIOGroup = getNumPEsInIOGroup()
   io_comm = getIOCommunicator()
   call setCommunicator(io_comm,MyPEinIOGroup,NumPEsInGroup,sync=.true.)
!  -------------------------------------------------------------------
   call c_dtsize(integer4_size,real8_size)
!  -------------------------------------------------------------------
!
   if ( isInputProc() ) then
!     ----------------------------------------------------------------
      call c_gopen(cunit,trim(adjustl(fname)),len_trim(adjustl(fname)), &
                   'READ',4,'XDR',3)
!     ----------------------------------------------------------------
      fp_pos=1
!     ----------------------------------------------------------------
      call c_fseek(cunit,fp_pos,0)
      call c_read_integer(cunit,block_size,2*GlobalNumAtoms)
!     ----------------------------------------------------------------
      num_clients = getNumInputClients()
      do i = 1, num_clients
!        -------------------------------------------------------------
         proc_client = getInputClient(i)
         call sendMessage(block_size,2,GlobalNumAtoms,212232,proc_client)
!        -------------------------------------------------------------
      enddo
   else
!     ----------------------------------------------------------------
      call recvMessage(block_size,2,GlobalNumAtoms,212232,getMyInputProc())
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  Determine the address of the core density data for each atomic site
!  in the core density file.
!  ===================================================================
   imsgbuf_size = 0
   fmsgbuf_size = 0
   do ig = 1, GlobalNumAtoms
      file_loc(ig) = block_size(1,ig)*integer4_size + block_size(2,ig)*real8_size
      imsgbuf_size = max(block_size(1,ig),imsgbuf_size)
      fmsgbuf_size = max(block_size(2,ig),fmsgbuf_size)
   enddo
   do ig = 2, GlobalNumAtoms
      file_loc(ig) = file_loc(ig) + file_loc(ig-1)
   enddo
   do ig = GlobalNumAtoms+1, 2, -1
      file_loc(ig) = file_loc(ig-1) + 1
   enddo
   file_loc(1) = 1
   file_loc = file_loc + 2*GlobalNumAtoms*integer4_size
!
   allocate(imsgbuf(imsgbuf_size), fmsgbuf(fmsgbuf_size))
   do id = 1, LocalNumAtoms
      ig = getGlobalIndex(id)
      if ( isInputProc() ) then
!        =============================================================
!        read in the core density data...................................
!        =============================================================
         fp_pos=file_loc(ig)
!        -------------------------------------------------------------
         call c_fseek(cunit,fp_pos,0)
         call c_read_integer(cunit,imsgbuf,block_size(1,ig))
!        -------------------------------------------------------------
         fp_pos = fp_pos + block_size(1,ig)*integer4_size
!        -------------------------------------------------------------
         call c_fseek(cunit,fp_pos,0)
         call c_read_double(cunit,fmsgbuf,block_size(2,ig))
!        -------------------------------------------------------------
!
         num_clients = getNumInputClients()
         do i = 1, num_clients
            proc_client = getInputClient(i)
            present_atom = getGlobalIndex(id,proc_client)
            fp_pos=file_loc(present_atom)
!           ----------------------------------------------------------
            call c_fseek(cunit,fp_pos,0)
            call c_read_integer(cunit,imsgbuf,block_size(1,present_atom))
!           ----------------------------------------------------------
            fp_pos = fp_pos + block_size(1,present_atom)*integer4_size
!           ----------------------------------------------------------
            call c_fseek(cunit,fp_pos,0)
            call c_read_double(cunit,fmsgbuf,block_size(2,present_atom))
!           ----------------------------------------------------------
            call sendMessage(imsgbuf,block_size(1,present_atom),212233,proc_client)
            call sendMessage(fmsgbuf,block_size(2,present_atom),212234,proc_client)
!           ----------------------------------------------------------
         enddo
      else
!        -------------------------------------------------------------
         call recvMessage(imsgbuf,block_size(1,ig),212233,getMyInputProc())
         call recvMessage(fmsgbuf,block_size(2,ig),212234,getMyInputProc())
!        -------------------------------------------------------------
      endif
!
!     ================================================================
!     Restore the data read in from the core density file.
!     ================================================================
      if (imsgbuf(1) /= ig) then
         call ErrorHandler('readCoreDensity','Inconsistent global index',imsgbuf(1),ig)
      endif
      if (Core(id)%NumSpecies /= imsgbuf(2)) then
         call ErrorHandler('readCoreDensity','Inconsistent number of species',imsgbuf(2),Core(id)%NumSpecies)
      endif
      Core(id)%MaxNumc = imsgbuf(3)
      rsize = imsgbuf(4)
      jcore = imsgbuf(5)
      jinsc = imsgbuf(6)
      isize = 6
!
      xstart = fmsgbuf(1)
      hin = fmsgbuf(2)
      hout = fmsgbuf(3)
      vol_core = fmsgbuf(4)
      vol_coreint = fmsgbuf(5)
      rcore_mt = fmsgbuf(6)
      fsize = 6
!
!
!     ================================================================
!     In case the radial grid associated with the core density data differs
!     from the current radial grid, we need to interpolate the data onto
!     the current radial grid. 
!     ================================================================
      if (abs(xstart-Core(id)%Grid%xstart) > TEN2m6 .or.              &
          abs(hin-Core(id)%Grid%hin) > TEN2m6 .or.                    &
          abs(hout-Core(id)%Grid%hout) > TEN2m6 .or.                  &
          abs(rcore_mt-Core(id)%rcore_mt) > TEN2m6 .or.               &
          jinsc /= Core(id)%Grid%jinsc .or.                           &
          rsize /= Core(id)%rsize) then
!        =============================================================
!        generate radial grid for the Core density read-in ...........
!        =============================================================
         allocate(x_mesh(rsize))
!        -------------------------------------------------------------
         call dcopy(rsize,fmsgbuf(fsize+1),1,x_mesh,1)
!        -------------------------------------------------------------
         do j=1, rsize
            x_mesh(j)=log(x_mesh(j))
         enddo
!
!        =============================================================
!        interpolate core density onto the existing radial mesh ......
!        =============================================================
         if (rsize < n_inter) then
!           ----------------------------------------------------------
            call ErrorHandler('readCoreDensity','rsize < n_inter',rsize,n_inter)
!           ----------------------------------------------------------
         endif
         Interp = .true.
      else
!        -------------------------------------------------------------
         call dcopy(rsize,fmsgbuf(fsize+1),1,Core(id)%r_mesh,1)
!        -------------------------------------------------------------
         fsize = fsize + rsize
         Interp = .false.
      endif
!
      do ia = 1, Core(id)%NumSpecies
         Core(id)%numc_below(ia) = imsgbuf(isize+1)
         Core(id)%numc(ia) = imsgbuf(isize+2)
         isize = isize + 2
         do ic = 1, Core(id)%numc(ia)
            Core(id)%nc(ic,ia) = imsgbuf(isize+1)
            Core(id)%lc(ic,ia) = imsgbuf(isize+2)
            Core(id)%kc(ic,ia) = imsgbuf(isize+3)
            isize = isize + 3
         enddo
!
         do is = 1, n_spin_pola
            Core(id)%ecorv(is,ia) = fmsgbuf(fsize+1)
            Core(id)%esemv(is,ia) = fmsgbuf(fsize+2)
            Core(id)%core_ke(is,ia) = fmsgbuf(fsize+3)
            Core(id)%semi_ke(is,ia) = fmsgbuf(fsize+4)
            fsize = fsize + 4
         enddo
         Core(id)%ztotss(ia) = fmsgbuf(fsize+1)
         Core(id)%zsemss(ia) = fmsgbuf(fsize+2)
         Core(id)%zcorss(ia) = fmsgbuf(fsize+3)
         Core(id)%qcpsc_mt(ia) = fmsgbuf(fsize+4)
         Core(id)%qcpsc_ws(ia) = fmsgbuf(fsize+5)
         Core(id)%mcpsc_mt(ia) = fmsgbuf(fsize+6)
         Core(id)%mcpsc_ws(ia) = fmsgbuf(fsize+7)
         Core(id)%qsemmt(ia) = fmsgbuf(fsize+8)
         Core(id)%qcormt(ia) = fmsgbuf(fsize+9)
         Core(id)%qsemws(ia) = fmsgbuf(fsize+10)
         Core(id)%qcorws(ia) = fmsgbuf(fsize+11)
         Core(id)%qcorout(ia) = fmsgbuf(fsize+12)
         fsize = fsize + 12
         do is = 1, n_spin_pola
            do ic = 1, Core(id)%numc(ia)
               Core(id)%ec(ic,is,ia) = fmsgbuf(fsize+1)
               fsize = fsize + 1
            enddo
            if (Interp) then
               j_inter=1
               do j=1, Core(id)%rsize
!                 -------------------------------------------------------
                  call hunt(rsize,x_mesh,Core(id)%Grid%x_mesh(j),j_inter)
!                 -------------------------------------------------------
                  if (j_inter > rsize-(n_inter-1)/2) then
                     irp=rsize-n_inter+1
                  else if (2*j_inter+1 > n_inter) then
                     irp=j_inter-(n_inter-1)/2
                  else
                     irp=1
                  endif
                  jrp = fsize + irp
!                 ----------------------------------------------------
                  call PolyInterp(n_inter,x_mesh(irp:irp+n_inter-1),  &
                                  fmsgbuf(jrp:jrp+n_inter-1),         &
                                  Core(id)%Grid%x_mesh(j),Core(id)%corden(j,is,ia),err)
!                 ----------------------------------------------------
                  jrp = jrp + rsize
!                 ----------------------------------------------------
                  call PolyInterp(n_inter,x_mesh(irp:irp+n_inter-1),  &
                                  fmsgbuf(jrp:jrp+n_inter-1),         &
                                  Core(id)%Grid%x_mesh(j),Core(id)%semden(j,is,ia),err)
!                 ----------------------------------------------------
               enddo
               fsize = fsize + 2*rsize
            else
!              -------------------------------------------------------
               call dcopy(rsize,fmsgbuf(fsize+1),1,Core(id)%corden(1,is,ia),1)
!              -------------------------------------------------------
               fsize = fsize + rsize
!              -------------------------------------------------------
               call dcopy(rsize,fmsgbuf(fsize+1),1,Core(id)%semden(1,is,ia),1)
!              -------------------------------------------------------
               fsize = fsize + rsize
            endif
         enddo
      enddo
!
      if (Interp) then
         deallocate(x_mesh)
      endif
   enddo
!
   if( isInputProc() ) then
      call c_close(cunit)
   endif
!
   call resetCommunicator()
!
   deallocate(imsgbuf, fmsgbuf)
!
   end subroutine readCoreDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeCoreDensity(fname)
!  ===================================================================
   use MPPModule, only : MyPE, GlobalSum, setCommunicator, resetCommunicator
   use MPPModule, only : nbsendMessage, recvMessage, waitMessage
!
   use ParallelIOModule, only : isOutputProc, getMyOutputProc,        &
                                getNumOutputClients, getOutputClient, &
                                getIOCommunicator, getMyPEinIOGroup,  &
                                getNumPEsInIOGroup, isOutputHeadProc
!
   use Atom2ProcModule, only : getGlobalIndex
!
   implicit none
!
   character (len=*), intent(in) :: fname
!
   integer (kind=IntKind) :: integer4_size,real8_size
   integer (kind=IntKind) :: cunit, i, id, ia, ic, ig, is, fp_pos
   integer (kind=IntKind) :: isize, fsize, imsgbuf_size, fmsgbuf_size
   integer (kind=IntKind) :: num_clients, proc_client, present_atom
!
   integer (kind=IntKind) :: MyPEinIOGroup, NumPEsInIOGroup, io_comm
   integer (kind=IntKind) :: msgid1, msgid2
!
!  Assuming that the size of the core density file does not exceed 2**31 bytes.
   integer (kind=IntKind) :: file_loc(GlobalNumAtoms+1) 
   integer (kind=IntKind) :: block_size(2,GlobalNumAtoms)
!
   integer (kind=IntKind), allocatable :: imsgbuf(:)
!
   real (kind=RealKind), allocatable :: fmsgbuf(:) 
!
   if (getMyOutputProc() < 0) then
      return
   endif
!
   MyPEinIOGroup = getMyPEinIOGroup()
   NumPEsInIOGroup = getNumPEsInIOGroup()
   io_comm = getIOCommunicator()
   call setCommunicator(io_comm,MyPEinIOGroup,NumPEsInGroup,sync=.true.)
!
!  -------------------------------------------------------------------
   call c_dtsize(integer4_size,real8_size)
!  -------------------------------------------------------------------
!
   if ( isOutputProc() ) then
!     ----------------------------------------------------------------
      call c_gopen(cunit,trim(adjustl(fname)),len_trim(adjustl(fname)), &
                   'WRITE',5,'XDR',3)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  Determine the total bytes of the core density data contributed by 
!  each atomic sites...
!  ===================================================================
   block_size = 0
   imsgbuf_size = 0
   fmsgbuf_size = 0
   do id = 1, LocalNumAtoms
!     ================================================================
!     The following data are to be stored
!
!     integer :: GlobalIndex
!     integer :: Core(id)%NumSpecies
!     integer :: Core(id)%MaxNumc
!     integer :: Core(id)%rsize
!     integer :: Core(id)%jcore
!     integer :: Core(id)%Grid%jinsc
!     real (kind=RealKind) :: Core(id)%Grid%xstart
!     real (kind=RealKind) :: Core(id)%Grid%hin
!     real (kind=RealKind) :: Core(id)%Grid%hout
!     real (kind=RealKind) :: Core(id)%vol_core
!     real (kind=RealKind) :: Core(id)%vol_coreint
!     real (kind=RealKind) :: Core(id)%rcore_mt
!     ================================================================
      ig = getGlobalIndex(id)
!     file_loc(ig) = file_loc(ig) + 6*integer4_size + 6*real8_size
      isize = 6
      fsize = 6
!
!     ================================================================
!     The following data are to be stored
!     real (kind=RealKind) :: Core(id)%r_mesh(:)
!     ================================================================
      fsize = fsize + Core(id)%rsize
!
      do ia = 1, Core(id)%NumSpecies
!        =============================================================
!        The following data are to be stored
!
!        integer :: Core(id)%numc_below(ia)
!        integer :: Core(id)%numc(ia)
!        real (kind=RealKind) :: Core(id)%ecorv(:,:)
!        real (kind=RealKind) :: Core(id)%esemv(:,:)
!        real (kind=RealKind) :: Core(id)%core_ke(:,:)
!        real (kind=RealKind) :: Core(id)%semi_ke(:,:)
!        real (kind=RealKind) :: Core(id)%ztotss(:)
!        real (kind=RealKind) :: Core(id)%zsemss(:)
!        real (kind=RealKind) :: Core(id)%zcorss(:)
!        real (kind=RealKind) :: Core(id)%qcpsc_mt(:)
!        real (kind=RealKind) :: Core(id)%qcpsc_ws(:)
!        real (kind=RealKind) :: Core(id)%mcpsc_mt(:)
!        real (kind=RealKind) :: Core(id)%mcpsc_ws(:)
!        real (kind=RealKind) :: Core(id)%qsemmt(:)
!        real (kind=RealKind) :: Core(id)%qcormt(:)
!        real (kind=RealKind) :: Core(id)%qsemws(:)
!        real (kind=RealKind) :: Core(id)%qcorws(:)
!        real (kind=RealKind) :: Core(id)%qcorout(:)
!        =============================================================
!        file_loc(ig) = file_loc(ig) + 2*integer4_size + (n_spin_pola*4+12)*real8_size
         isize = isize + 2
         fsize = fsize + (n_spin_pola*4+12)
!
         do ic = 1, Core(id)%numc(ia)
!           ==========================================================
!           The following data are to be stored
!
!           integer :: Core(id)%nc(:,:)
!           integer :: Core(id)%lc(:,:)
!           integer :: Core(id)%kc(:,:)
!           real (kind=RealKind) :: Core(id)%ec(:,:,:)
!           ==========================================================
!           file_loc(ig) = file_loc(ig) + 3*integer4_size + n_spin_pola*real8_size
            isize = isize + 3
            fsize = fsize + n_spin_pola
         enddo
!        =============================================================
!        The following data are to be stored
!
!        real (kind=RealKind) :: Core(id)%corden(:,:,:)
!        real (kind=RealKind) :: Core(id)%semden(:,:,:)
!        =============================================================
!        file_loc(ig) = file_loc(ig) + 2*Core(id)%rsize*n_spin_pola*real8_size
         fsize = fsize + 2*Core(id)%rsize*n_spin_pola
      enddo
      imsgbuf_size = max(isize,imsgbuf_size)
      fmsgbuf_size = max(fsize,fmsgbuf_size)
      block_size(1,ig) = isize
      block_size(2,ig) = fsize
   enddo
!  -------------------------------------------------------------------
   call GlobalSum(block_size,2,GlobalNumAtoms)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  The first piece of information is written which records the size of 
!  the core density data for each atomic site. This information will be
!  used to determine the location of the core density data of an atomic 
!  site when the file is read.
!  ===================================================================
   if ( isOutputHeadProc() ) then
      fp_pos = 1
!     ----------------------------------------------------------------
      call c_fseek(cunit,fp_pos,0)
      call c_write_integer(cunit,block_size,2*GlobalNumAtoms)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  Determine the address of the core density data for each atomic site
!  in the core density file.
!  ===================================================================
   do ig = 1, GlobalNumAtoms
      file_loc(ig) = block_size(1,ig)*integer4_size + block_size(2,ig)*real8_size
   enddo
   do ig = 2, GlobalNumAtoms
      file_loc(ig) = file_loc(ig) + file_loc(ig-1)
   enddo
   do ig = GlobalNumAtoms+1, 2, -1
      file_loc(ig) = file_loc(ig-1) + 1
   enddo
   file_loc(1) = 1
   file_loc = file_loc + 2*GlobalNumAtoms*integer4_size
!  ===================================================================
!
!  ===================================================================
!  Store integer datat in imsgbuf, real data in fmsgbuf, and write them
!  to the core density file.
!  ===================================================================
   allocate(imsgbuf(imsgbuf_size), fmsgbuf(fmsgbuf_size))
   do id = 1, LocalNumAtoms
      ig = getGlobalIndex(id)
      imsgbuf(1) = ig
      imsgbuf(2) = Core(id)%NumSpecies
      imsgbuf(3) = Core(id)%MaxNumc
      imsgbuf(4) = Core(id)%rsize
      imsgbuf(5) = Core(id)%jcore
      imsgbuf(6) = Core(id)%Grid%jinsc
      isize = 6
!
      fmsgbuf(1) = Core(id)%Grid%xstart
      fmsgbuf(2) = Core(id)%Grid%hin
      fmsgbuf(3) = Core(id)%Grid%hout
      fmsgbuf(4) = Core(id)%vol_core
      fmsgbuf(5) = Core(id)%vol_coreint
      fmsgbuf(6) = Core(id)%rcore_mt
      fsize = 6
!
!     ----------------------------------------------------------------
      call dcopy(Core(id)%rsize,Core(id)%r_mesh,1,fmsgbuf(fsize+1),1)
!     ----------------------------------------------------------------
      fsize = fsize + Core(id)%rsize
!
      do ia = 1, Core(id)%NumSpecies
         imsgbuf(isize+1) = Core(id)%numc_below(ia)
         imsgbuf(isize+2) = Core(id)%numc(ia)
         isize = isize + 2
         do ic = 1, Core(id)%numc(ia)
            imsgbuf(isize+1) = Core(id)%nc(ic,ia)
            imsgbuf(isize+2) = Core(id)%lc(ic,ia)
            imsgbuf(isize+3) = Core(id)%kc(ic,ia)
            isize = isize + 3
         enddo
!
         do is = 1, n_spin_pola
            fmsgbuf(fsize+1) = Core(id)%ecorv(is,ia)
            fmsgbuf(fsize+2) = Core(id)%esemv(is,ia)
            fmsgbuf(fsize+3) = Core(id)%core_ke(is,ia)
            fmsgbuf(fsize+4) = Core(id)%semi_ke(is,ia)
            fsize = fsize + 4
         enddo
         fmsgbuf(fsize+1) = Core(id)%ztotss(ia)
         fmsgbuf(fsize+2) = Core(id)%zsemss(ia)
         fmsgbuf(fsize+3) = Core(id)%zcorss(ia)
         fmsgbuf(fsize+4) = Core(id)%qcpsc_mt(ia)
         fmsgbuf(fsize+5) = Core(id)%qcpsc_ws(ia)
         fmsgbuf(fsize+6) = Core(id)%mcpsc_mt(ia)
         fmsgbuf(fsize+7) = Core(id)%mcpsc_ws(ia)
         fmsgbuf(fsize+8) = Core(id)%qsemmt(ia)
         fmsgbuf(fsize+9) = Core(id)%qcormt(ia)
         fmsgbuf(fsize+10) = Core(id)%qsemws(ia)
         fmsgbuf(fsize+11) = Core(id)%qcorws(ia)
         fmsgbuf(fsize+12) = Core(id)%qcorout(ia)
         fsize = fsize + 12
         do is = 1, n_spin_pola
            do ic = 1, Core(id)%numc(ia)
               fmsgbuf(fsize+1) = Core(id)%ec(ic,is,ia)
               fsize = fsize + 1
            enddo
!           ----------------------------------------------------------
            call dcopy(Core(id)%rsize,Core(id)%corden(1,is,ia),1,fmsgbuf(fsize+1),1)
!           ----------------------------------------------------------
            fsize = fsize + Core(id)%rsize
!           ----------------------------------------------------------
            call dcopy(Core(id)%rsize,Core(id)%semden(1,is,ia),1,fmsgbuf(fsize+1),1)
!           ----------------------------------------------------------
            fsize = fsize + Core(id)%rsize
         enddo
      enddo
      if (isize > imsgbuf_size) then
!        -------------------------------------------------------------
         call ErrorHandler('writeCoreDensity','integer data size > buffer size')
!        -------------------------------------------------------------
      else if (fsize > fmsgbuf_size) then
!        -------------------------------------------------------------
         call ErrorHandler('writeCoreDensity','Reral data size > buffer size')
!        -------------------------------------------------------------
      else if (isize /= block_size(1,ig)) then
!        -------------------------------------------------------------
         call ErrorHandler('writeCoreDensity','Inconsistent integer data size')
!        -------------------------------------------------------------
      else if (fsize /= block_size(2,ig)) then 
!        -------------------------------------------------------------
         call ErrorHandler('writeCoreDensity','Inconsistent real data size')
!        -------------------------------------------------------------
      else if (fsize*real8_size+isize*integer4_size /= file_loc(ig+1)-file_loc(ig)) then
!        -------------------------------------------------------------
         call ErrorHandler('writeCoreDensity','Inconsistent data location and size')
!        -------------------------------------------------------------
      endif
!
      if ( isOutputProc() ) then
!        =============================================================
!        write out imsgbuf of the present local atom..................
!        =============================================================
         fp_pos=file_loc(ig)
!        -------------------------------------------------------------
         call c_fseek(cunit,fp_pos,0)
         call c_write_integer(cunit,imsgbuf,isize)
!        -------------------------------------------------------------
         fp_pos = fp_pos + isize*integer4_size
!        -------------------------------------------------------------
         call c_fseek(cunit,fp_pos,0)
         call c_write_double(cunit,fmsgbuf,fsize)
!        -------------------------------------------------------------
!
         num_clients = getNumOutputClients()
         do i = 1, num_clients
            proc_client = getOutputClient(i)
            present_atom = getGlobalIndex(id,proc_client)
!           ==========================================================
!           Receive clients data
!           ----------------------------------------------------------
            call recvMessage(imsgbuf,block_size(1,present_atom),112233,proc_client)
            call recvMessage(fmsgbuf,block_size(2,present_atom),112234,proc_client)
!           ----------------------------------------------------------
!
!           ==========================================================
!           Write clients data
!           ==========================================================
            fp_pos=file_loc(present_atom)
!           ----------------------------------------------------------
            call c_fseek(cunit,fp_pos,0)
            call c_write_integer(cunit,imsgbuf,isize)
!           ----------------------------------------------------------
            fp_pos = fp_pos + isize*integer4_size
!           ----------------------------------------------------------
            call c_fseek(cunit,fp_pos,0)
            call c_write_double(cunit,fmsgbuf,fsize)
!           ----------------------------------------------------------
         enddo
      else
!        -------------------------------------------------------------
         msgid1 = nbsendMessage(imsgbuf,isize,112233,getMyOutputProc())
         msgid2 = nbsendMessage(fmsgbuf,fsize,112234,getMyOutputProc())
         call waitMessage(msgid1)
         call waitMessage(msgid2)
!        -------------------------------------------------------------
      endif
   enddo
!
   if( isOutputProc() ) then
      call c_close(cunit)
   endif
!
   call resetCommunicator()
!
   deallocate(imsgbuf, fmsgbuf)
!
   end subroutine writeCoreDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calCoreStates(evb)
!  ===================================================================
   use GroupCommModule, only : GlobalMaxInGroup, GlobalSumInGroup
!
   use Atom2ProcModule, only : getLocalIndex, getAtom2ProcInGroup,     &
                               getMaxLocalNumAtoms
!
   use AtomModule, only : getLocalSpeciesContent
!
   use SystemVolumeModule, only : getTotalInterstitialMTVolume,        &
                                  getSystemVolume, getAtomicMTVolume
!
   use PolyhedraModule, only : getVolume
!
   use PotentialTypeModule, only : isMuffinTinPotential, isASAPotential, &
                                   isMuffinTinASAPotential, isFullPotential
!
   use PotentialModule, only : getSphPotr
!
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       isDataStorageExisting, &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark, &
                                       IntegerMark, IntegerType
!
   use ChemElementModule, only : getCoreStateIndex
!
   use InterpolationModule, only : FitInterp
!
   implicit   none
!
   character (len=13), parameter ::  sname='calCoreStates'
!
   integer (kind=IntKind), parameter :: formula = 0
!
   integer (kind=IntKind) :: id, ia
   integer (kind=IntKind) :: is, ig, i, j, ir
   integer (kind=IntKind) :: nr, nws, nend, last, last2, nmult, jmt, jend_plus_n
   integer (kind=IntKind) :: DataSize(LocalNumAtoms)
   integer (kind=IntKind), pointer :: pc0(:,:,:)
!
   real (kind=RealKind), optional, intent(out) :: evb
   real (kind=RealKind), allocatable :: wrk1(:),wrk2(:)
   real (kind=RealKind), parameter :: tolch=TEN2m5
   real (kind=RealKind) :: h, hout, efact, dps
   real (kind=RealKind) :: qint_semi, qint_deep, msgbuf(2)
   real (kind=RealKind), allocatable :: vr_t(:)
   real (kind=RealKind), allocatable :: memtemp(:,:)
   real (kind=RealKind), pointer :: ec0(:,:,:)
   real (kind=RealKind), pointer :: vp(:)
   real (kind=RealKind), parameter :: PI8 = PI4*TWO
   real (kind=RealKind) :: msemmt, mcormt, msemws, mcorws
!
   real (kind=RealKind), allocatable :: sqrt_r(:)
   real (kind=RealKind), allocatable :: r_mesh_t(:)
   real (kind=RealKind), allocatable :: sqrt_r_t(:)
!  ===================================================================
!  setup single grids for r_mesh...................................
!  I chhose a large last value to make this code agree with the older version
!  ===================================================================
   last = 0
   do id =1, LocalNumAtoms
      last = max(last,Core(id)%rsize)
      DataSize(id) = Core(id)%MaxNumc*n_spin_pola*Core(id)%NumSpecies
      if (Core(id)%MaxNumc < 1) then
         DataSize(id) = n_spin_pola
      endif
   enddo
   allocate(vr_t(1:last))
   allocate(wrk1(0:last), wrk2(0:last), sqrt_r(0:last), sqrt_r_t(0:last), r_mesh_t(1:last))
!
   if (.not.isDataStorageExisting('NewEstimatedCoreEnergy')) then
!     ---------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'NewEstimatedCoreEnergy',  &
                             DataSize,RealType)
      call createDataStorage(LocalNumAtoms,'NewEstimatedCoreStates',  &
                             DataSize,IntegerType)
!     ---------------------------------------------------------------
   endif
!
   msgbuf(1:2) = ZERO
   etopcor=-10.0d+20
   do id =1, LocalNumAtoms
      jend_plus_n=Core(id)%Grid%jend_plus_n
      jmt=Core(id)%jcore
      nws=Core(id)%Grid%jws
      last=Core(id)%rsize      ! Definition of last: the final mesh point.
!
!     ----------------------------------------------------------------
      call dcopy(jend_plus_n,Core(id)%Grid%r_mesh,1,Core(id)%r_mesh,1)
!     ----------------------------------------------------------------
      hout = Core(id)%Grid%hout
      do ir = jend_plus_n+1,last
         Core(id)%r_mesh(ir) = exp(hout*(ir-Core(id)%Grid%jmt)+Core(id)%Grid%xmt)
      enddo
!
      sqrt_r(0) = ZERO
      do ir=1,last
         sqrt_r(ir)=sqrt(Core(id)%r_mesh(ir))
      enddo
!
      h = Core(id)%Grid%hin
      nmult=Core(id)%Grid%nmult
      if (nmult > 1) then
         r_mesh_t(1:Core(id)%Grid%jmt) = Core(id)%r_mesh(1:Core(id)%Grid%jmt)
         do ir = Core(id)%Grid%jmt+1, last
            r_mesh_t(ir) = exp(h*(ir-Core(id)%Grid%jmt)+Core(id)%Grid%xmt)
         enddo
      else
         r_mesh_t(1:last) = Core(id)%r_mesh(1:last)
      endif
!
      sqrt_r_t(0) = ZERO
      do ir = 1,last
         sqrt_r_t(ir) = sqrt(r_mesh_t(ir))
      enddo
!
!     ================================================================
!     calculate the core and semi-core densities......................
!     Definition of last2: the last mesh point at the cell boundary, or
!                          on the Wigner-Seitz sphere surface.
!     Definition of nend:  the last mesh point of uniform mesh of step hin
!                          at (or somewhat further than) the cell boundary, 
!                          or on (or somewhat further than) the Wigner-Seitz 
!                          sphere surface.
!     Notice the difference in the definition of last2 between this version 
!     and earlier versions. This difference causes the difference in the semi-core
!     kinetic energy and ultimately the difference in the total energy.
!     ================================================================
      if (isMuffinTinPotential() .or. isFullPotential() ) then
         last2 = Core(id)%Grid%jend   !  Definition of last2:
      else
         last2 = nws
      endif
      if (nmult > 1) then
         nend = last2 - mod(last2-Core(id)%Grid%jmt,nmult) + nmult  ! In case hout /= hin.
      else
         nend = last2
      endif
!
      do ia = 1, Core(id)%NumSpecies
         ndeepz=int(Core(id)%zcorss(ia)/n_spin_pola+HALF,IntKind)
!
         Core(id)%corden(:,:,ia)=ZERO
         Core(id)%semden(:,:,ia)=ZERO
         Core(id)%dcorden(:,:,ia)=ZERO
         Core(id)%dsemden(:,:,ia)=ZERO
         Core(id)%ecorv(1:n_spin_pola,ia)=ZERO
         Core(id)%esemv(1:n_spin_pola,ia)=ZERO
         Core(id)%core_ke(1:n_spin_pola,ia)=ZERO
         Core(id)%semi_ke(1:n_spin_pola,ia)=ZERO
         Core(id)%qsemmt(ia)=ZERO
         Core(id)%qcormt(ia)=ZERO
         Core(id)%qsemws(ia)=ZERO
         Core(id)%qcorws(ia)=ZERO
         Core(id)%mcpsc_mt(ia)=ZERO
         Core(id)%mcpsc_ws(ia)=ZERO
         Core(id)%qcorout(ia)=ZERO
!
         if(Core(id)%numc(ia).gt.0) then
            msemmt=ZERO; mcormt=ZERO; msemws=ZERO; mcorws=ZERO
            do is=1,n_spin_pola
!              =======================================================
!              setup single grids for vr_t............................
!              =======================================================
               vp=>getSphPotr(id,ia,is)
               if (nmult > 1) then
                  vr_t = ZERO
                  do ir = 1,Core(id)%Grid%jmt
                     vr_t(ir) = vp(ir)
                  enddo
!
                  ir = Core(id)%Grid%jmt+1
                  nr = Core(id)%Grid%jmt+nmult
                  do while (nr <= Core(id)%Grid%jend)
                     vr_t(ir) = vp(nr)
                     nr = nr+nmult
                     ir = ir+1
                  enddo
               else
                  vr_t = ZERO
                  call dcopy(Core(id)%Grid%jend,vp,1,vr_t,1)
               endif
!              -------------------------------------------------------
               call corslv(id,ia,is,nws,last,nend,h,r_mesh_t,sqrt_r_t(0:last),vr_t)
!              -------------------------------------------------------
!
!              =======================================================
!              The following two lines were added by Yang on 11/23/2018
!              =======================================================
               Core(id)%semden(last2+1:,is,ia)=ZERO
               Core(id)%corden(last2+1:,is,ia)=ZERO
!              =======================================================
!
!              =======================================================
!              adjust bottom of cotour: above highest semi-core state
!              =======================================================
               etopcor=max(Core(id)%ec(Core(id)%numc_below(ia),is,ia),etopcor)
!
!              =======================================================
!              check the charge of semi-cores.........................
!              -------------------------------------------------------
               wrk2(0)=ZERO
               do j = 1, last2
                  wrk2(j)=Core(id)%semden(j,is,ia)*Core(id)%r_mesh(j)
               enddo
!              -------------------------------------------------------
               call FitInterp( 4, sqrt_r(1:4), wrk2(1:4), ZERO, wrk2(0), dps )
!              -------------------------------------------------------
!              call calIntegration(0,last2,Core(id)%r_mesh,wrk2,wrk1)
               call calIntegration(last2+1,sqrt_r(0:last2),wrk2(0:last2),wrk1(0:last2),3)
!              -------------------------------------------------------
!              wrk1(jmt) = (wrk1(jmt) + HALF*wrk2(1)*Core(id)%r_mesh(1))*PI4
!              wrk1(last2) = (wrk1(last2) + HALF*wrk2(1)*Core(id)%r_mesh(1))*PI4
               Core(id)%qsemmt(ia)=Core(id)%qsemmt(ia)+wrk1(jmt)*PI8
               Core(id)%qsemws(ia)=Core(id)%qsemws(ia)+wrk1(last2)*PI8
               msemmt=msemmt+(n_spin_pola+1-2*is)*wrk1(jmt)*PI8
               msemws=msemws+(n_spin_pola+1-2*is)*wrk1(last2)*PI8
!              =======================================================
!
!              =======================================================
!              check the charge of deep-cores.........................
!              -------------------------------------------------------
               wrk2(0)=ZERO
               do j = 1, last2
                  wrk2(j)=Core(id)%corden(j,is,ia)*Core(id)%r_mesh(j)
               enddo
!              -------------------------------------------------------
               call FitInterp( 4, sqrt_r(1:4), wrk2(1:4), ZERO, wrk2(0), dps )
!              -------------------------------------------------------
!              call calIntegration(0,last2,Core(id)%r_mesh,wrk2,wrk1)
               call calIntegration(last2+1,sqrt_r(0:last2),wrk2(0:last2),wrk1(0:last2),3)
!              -------------------------------------------------------
!              wrk1(jmt) = (wrk1(jmt) + HALF*wrk2(1)*Core(id)%r_mesh(1))*PI4
!              wrk1(last2) = (wrk1(last2) + HALF*wrk2(1)*Core(id)%r_mesh(1))*PI4
               Core(id)%qcormt(ia)=Core(id)%qcormt(ia)+wrk1(jmt)*PI8
               Core(id)%qcorws(ia)=Core(id)%qcorws(ia)+wrk1(last2)*PI8
               mcormt=mcormt+(n_spin_pola+1-2*is)*wrk1(jmt)*PI8
               mcorws=mcorws+(n_spin_pola+1-2*is)*wrk1(last2)*PI8
!
!              =======================================================
!              Compute the kinetic energy
!              =======================================================
               wrk2(0)=ZERO
               do j = 1, last2
                  wrk2(j)=Core(id)%semden(j,is,ia)*vp(j)
               enddo
!              -------------------------------------------------------
               call FitInterp( 4, sqrt_r(1:4), wrk2(1:4), ZERO, wrk2(0), dps )
!              -------------------------------------------------------
               call calIntegration(last2+1,sqrt_r(0:last2),wrk2(0:last2),wrk1(0:last2),3)
!              -------------------------------------------------------
!write(6,'(a,2d15.8)')'Semi core energy and int[rho*v] =',Core(id)%esemv(is,ia),wrk1(last2)*PI8
               Core(id)%semi_ke(is,ia)=Core(id)%esemv(is,ia)-wrk1(last2)*PI8
!
               if (formula == 0) then
                  wrk2(0)=ZERO
                  do j = 1, last2
                     wrk2(j)=Core(id)%corden(j,is,ia)*vp(j)
                  enddo
!                 ----------------------------------------------------
                  call FitInterp( 4, sqrt_r(1:4), wrk2(1:4), ZERO, wrk2(0), dps )
!                 ----------------------------------------------------
                  call calIntegration(last2+1,sqrt_r(0:last2),wrk2(0:last2),wrk1(0:last2),3)
!                 ----------------------------------------------------
!write(6,'(a,2d15.8)')'Deep core energy and int[rho*v] =',Core(id)%ecorv(is,ia),wrk1(last2)*PI8
                  Core(id)%core_ke(is,ia)=Core(id)%ecorv(is,ia)-wrk1(last2)*PI8
               else
!                 ---------------------------------------------------
                  call newder(vp(1:jmt),wrk1(1:jmt),sqrt_r(1:jmt),jmt)
!                 ---------------------------------------------------
                  do j = 1, jmt
                     wrk2(j)=Core(id)%r_mesh(j)*core(id)%corden(j,is,ia)*(HALF*sqrt_r(j)*wrk1(j)-vp(j))
                  enddo
!                 ---------------------------------------------------
                  call FitInterp(4,sqrt_r(1:4),wrk2(1:4),ZERO,wrk2(0),dps)
                  call calIntegration(jmt+1,sqrt_r(0:jmt),wrk2(0:jmt),wrk1(0:jmt),1)
!                 ---------------------------------------------------
                  Core(id)%core_ke(is,ia)=wrk1(jmt)*PI4
               endif
            enddo
            Core(id)%mcpsc_mt(ia)= msemmt+mcormt
            Core(id)%mcpsc_ws(ia)= msemws+mcorws
         endif
!
         Core(id)%qcorout(ia)=(Core(id)%zsemss(ia)-Core(id)%qsemmt(ia))+ &
                              (Core(id)%zcorss(ia)-Core(id)%qcormt(ia))
         Core(id)%qcpsc_mt(ia)=Core(id)%qcormt(ia)+Core(id)%qsemmt(ia)
         Core(id)%qcpsc_ws(ia)=Core(id)%qcorws(ia)+Core(id)%qsemws(ia)
!
!        ============================================================
!        check that semi-core and deep-core charge is not being lost.
!        ============================================================
         if(abs(Core(id)%qsemws(ia)-Core(id)%zsemss(ia)) > tolch     &
            .and. print_level >= 0) then
            write(6,'(/,10x,''Lost semi-core charge'')')
            write(6,'(10x,''z='',d17.8,'' q='',f17.8)') &
                             Core(id)%zsemss(ia), Core(id)%qsemws(ia)
         endif
         if(abs(Core(id)%qcorws(ia)-Core(id)%zcorss(ia)) > tolch     &
            .and. print_level >= 0) then
            write(6,'(/,10x,''Lost      core charge'')')
            write(6,'(10x,''z='',d17.8,'' q='',f17.8)') &
                             Core(id)%zcorss(ia), Core(id)%qcorws(ia)
         endif
         msgbuf(1) = msgbuf(1) + (Core(id)%zsemss(ia)-Core(id)%qsemmt(ia))*&
                     getLocalSpeciesContent(id,ia)/real(GlobalNumAtoms,kind=Realkind)
         msgbuf(2) = msgbuf(2) + (Core(id)%zcorss(ia)-Core(id)%qcormt(ia))*&
                     getLocalSpeciesContent(id,ia)/real(GlobalNumAtoms,kind=Realkind)
      enddo
   enddo
!
   do id = 1,LocalNumAtoms
      if (Core(id)%MaxNumc >= 1) then
         ec0  => getDataStorage( id, 'NewEstimatedCoreEnergy',        &
                                 Core(id)%MaxNumc,n_spin_pola,        &
                                 Core(id)%NumSpecies,RealMark )
         pc0  => getDataStorage( id, 'NewEstimatedCoreStates',        &
                                 Core(id)%MaxNumc,n_spin_pola,        &
                                 Core(id)%NumSpecies,IntegerMark )
      endif
      do ia = 1, Core(id)%NumSpecies
         do is = 1, n_spin_pola
            do i = 1,Core(id)%numc(ia)
               ec0(i,is,ia) = Core(id)%ec(i,is,ia)
               pc0(i,is,ia) = getCoreStateIndex(Core(id)%nc(i,ia),    &
                                                Core(id)%lc(i,ia),    &
                                                Core(id)%kc(i,ia))
            enddo
         enddo
      enddo
   enddo
!
!  ===================================================================
   deallocate(vr_t)
   deallocate(wrk1, wrk2, sqrt_r, sqrt_r_t, r_mesh_t)
   nullify(ec0, vp)
!  ===================================================================
!
   msgbuf(1:2) = msgbuf(1:2)
!  ===================================================================
!  Find the uppermost core state across all the nodes/atoms........
!  -------------------------------------------------------------------
   call GlobalMaxInGroup(GroupID,etopcor)
   call GlobalSumInGroup(GroupID,msgbuf,2)
!  -------------------------------------------------------------------
   qint_semi = msgbuf(1)
   qint_deep = msgbuf(2)
!
!  ==================================================================
!  check that uppermost core state is well below evbot............
!  ==================================================================
   if(etopcor+0.1d0 .gt. evbot) then
!     call ErrorHandler('calCoreStates','etopcor+0.1 > evbot',etopcor,evbot)
      call WarningHandler('calCoreStates','reset evbot to a new value',etopcor-0.2d0)
      evbot = etopcor - 0.2d0
   endif
!
   if (present(evb)) then
      evb = evbot
   endif
!
   if(abs(qint_semi) < TEN2m8) then
      qint_semi = ZERO
   else if(qint_semi < ZERO .and. abs(qint_semi) < TEN2m6) then
      if (print_level >= 0) then
         call WarningHandler('calCoreStates','qint_semi is set to 0',qint_semi)
      endif
      qint_semi = ZERO
   else if (qint_semi < ZERO) then
      call WarningHandler('calCoreStates','qint_semi<0, likely caused by the valence band bottom being set too high',qint_semi)
   endif
!
   if(abs(qint_deep) < TEN2m8) then
       qint_deep = ZERO
   else if(qint_deep < ZERO .and. abs(qint_deep) < TEN2m6) then
       if (print_level >= 0) then
          call WarningHandler('calCoreStates','qint_deep is set to 0',qint_deep)
       endif
       qint_deep = ZERO
   else if(qint_deep < ZERO) then
       call WarningHandler('calCoreStates','qint_deep<0, likely caused by the valence band bottom being set too high',qint_deep)
   endif
!
   if(isASAPotential()) then
       rhoint_semi =                                                  &
          real(GlobalNumAtoms,kind=Realkind)*qint_semi/(getSystemVolume()*n_spin_pola)
       rhoint_deep =                                                  &
          real(GlobalNumAtoms,kind=Realkind)*qint_deep/(getSystemVolume()*n_spin_pola)
!
       if (rhoint_semi > TEN2m8) then
          do id = 1, LocalNumAtoms
             jmt = Core(id)%jcore
             do ia = 1, Core(id)%NumSpecies
                do is = 1, n_spin_pola
                   do j=1,jmt
                      Core(id)%semden(j,is,ia) = Core(id)%semden(j,is,ia) + rhoint_semi
                   enddo
                enddo
             enddo
          enddo
       endif
!
       if (rhoint_deep > TEN2m8) then
          do id = 1, LocalNumAtoms
             jmt = Core(id)%jcore
             do ia = 1, Core(id)%NumSpecies
                do is = 1, n_spin_pola
                   do j=1,jmt
                      Core(id)%corden(j,is,ia) = Core(id)%corden(j,is,ia) + rhoint_deep
                   enddo
                enddo
             enddo
          enddo
       endif
!
       rhoint_semi = ZERO
       rhoint_deep = ZERO
   else
       if (TotalInterstitialCoreVol < TEN2m6) then
           call ErrorHandler('calCoreStates','no interstitial volume', &
                             TotalInterstitialCoreVol)
       endif
       rhoint_semi =                                                   &
          real(GlobalNumAtoms,kind=Realkind)*qint_semi/TotalInterstitialCoreVol
       rhoint_deep =                                                   &
          real(GlobalNumAtoms,kind=Realkind)*qint_deep/TotalInterstitialCoreVol
   endif
!
!   if ( isASAPotential() .or. isMuffinTinPotential() .or.              &
!        isFullPotential() ) then
!      do id = 1, LocalNumAtoms
!         jmt = Core(id)%jcore
!         do ia = 1, Core(id)%NumSpecies
!            do is = 1, n_spin_pola
!               if (abs(Core(id)%semden(jmt+1,is,ia)) < TEN2m8) then
!                  Core(id)%semden(jmt+1,is,ia) = Core(id)%semden(jmt,is,ia)
!               endif
!               if (abs(Core(id)%semden(jmt+2,is,ia)) < TEN2m8) then
!                  Core(id)%semden(jmt+2,is,ia) = Core(id)%semden(jmt+1,is,ia)
!               endif
!               if (abs(Core(id)%corden(jmt+1,is,ia)) < TEN2m8) then
!                  Core(id)%corden(jmt+1,is,ia) = Core(id)%corden(jmt,is,ia)
!               endif
!               if (abs(Core(id)%corden(jmt+2,is,ia)) < TEN2m8) then
!                  Core(id)%corden(jmt+2,is,ia) = Core(id)%corden(jmt+1,is,ia)
!               endif
!            enddo
!         enddo
!      enddo
!   else if ( isMuffinTinASAPotential() ) then
!      do id = 1, LocalNumAtoms
!         jmt = Core(id)%jcore
!         do ia = 1, Core(id)%NumSpecies
!            do is = 1, n_spin_pola
!               if (abs(Core(id)%semden(jmt+1,is,ia)) < TEN2m8) then
!                  Core(id)%semden(jmt+1,is,ia) = Core(id)%semden(jmt,is,ia)
!               endif
!               if (abs(Core(id)%semden(jmt+2,is,ia)) < TEN2m8) then
!                  Core(id)%semden(jmt+2,is,ia) = Core(id)%semden(jmt+1,is,ia)
!               endif
!               if (abs(Core(id)%corden(jmt+1,is,ia)) < TEN2m8) then
!                  Core(id)%corden(jmt+1,is,ia) = Core(id)%corden(jmt,is,ia)
!               endif
!               if (abs(Core(id)%corden(jmt+2,is,ia)) < TEN2m8) then
!                  Core(id)%corden(jmt+2,is,ia) = Core(id)%corden(jmt+1,is,ia)
!               endif
!            enddo
!         enddo
!      enddo
!   endif
!
   call updateGlobalCoreStatesTable()
!
   if(stop_routine.eq.sname) then
      do id =1, LocalNumAtoms
         call printCoreStates(id)
      enddo
      call StopHandler(sname)
   endif
!
   end subroutine calCoreStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printCoreStates(id)
!  ===================================================================
   implicit none
!
   character (len=15), parameter :: sname='printCoreStates'
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: ia,is,i
!
   if (id > LocalNumAtoms .or. id < 1) then
      call ErrorHandler(sname,'Invalid local atom index',id)
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)') '*********************************'
   write(6,'(  24x,a )')'*  Output from printCoreStates  *'
   write(6,'( 24x,a,/)')'*********************************'
!
   write(6,'(''Local Atom Index: '',i6)')id
   write(6,'(''jend, jcore     : '',2i6)')Core(id)%Grid%jend,Core(id)%jcore
   write(6,'(''Upper limit of core level'',t31,''='',f20.11)')etopcor
   do ia = 1, Core(id)%NumSpecies
      write(6,'(80(''=''))')
      if (Core(id)%NumSpecies > 1) then
         write(6,'(''Species Index   : '',i6)')ia
      endif
      write(6,'(''ztotss          : '',i6)')int(Core(id)%ztotss(ia))
      do is=1,n_spin_pola
         write(6,'(''Spin Index      : '',i6)')is
         write(6,                                                        &
            '(''Eigenvalues: '',t20,''n'',t25,''l'',t30,''k'',t41,''energy'')')
         do i=1,Core(id)%numc(ia)
            if (Core(id)%ec(i,is,ia) <= evbot) then
               write(6,'(t18,i3,t23,i3,t28,i3,t32,f20.11)')                 &
                  Core(id)%nc(i,ia),Core(id)%lc(i,ia),Core(id)%kc(i,ia),Core(id)%ec(i,is,ia)
            else
               write(6,'(t18,i3,t23,i3,t28,i3,t32,f20.11,a)')               &
                  Core(id)%nc(i,ia),Core(id)%lc(i,ia),Core(id)%kc(i,ia),    &
                  Core(id)%ec(i,is,ia),'--- treated as valence state ---'
            endif
         enddo
         write(6,'(/)')
         write(6,'(''Total Core Energy'',t31,''='',f20.11)')Core(id)%ecorv(is,ia)
         write(6,'(''Total Semicore Energy'',t31,''='',f20.11)')Core(id)%esemv(is,ia)
         write(6,'(/)')
         write(6,'(''Deep Core Kinetic Energy'',t31,''='',f20.11)')Core(id)%core_ke(is,ia)
         write(6,'(''Semi Core Kinetic Energy'',t31,''='',f20.11)')Core(id)%semi_ke(is,ia)
!
!        do j=1,last2,50
!           write(6,'(''j,corden,semden: '',i5,3f20.11)')              &
!              j,Core(id)%corden(j,is,ia),Core(id)%semden(j,is,ia)
!        enddo
         write(6,'(80(''=''))')
      enddo
!
      write(6,'(''Muffin-tin core charge'',t31,''='',f20.11)') Core(id)%qcormt(ia)
      write(6,'(''Muffin-tin semicore charge'',t31,''='',f20.11)')     &
                Core(id)%qsemmt(ia)
      write(6,'(''Muffin-tin core+semi moment'',t31,''='',f20.11)')    &
                Core(id)%mcpsc_mt(ia)
      write(6,'(''Wigner-Seitz core+semi moment'',t31,''='',f20.11)')  &
                Core(id)%mcpsc_ws(ia)
      write(6,'(''Interstitial charge core'',t31,''='',f20.11)')       &
                Core(id)%qcorout(ia)
      write(6,'(/,''Lost semi-core charge'',t31,''='',f20.11)')        &
                Core(id)%qsemws(ia)-Core(id)%zsemss(ia)
      write(6,'(  ''Lost deep-core charge'',t31,''='',f20.11)')        &
                Core(id)%qcorws(ia)-Core(id)%zcorss(ia)
   enddo
   write(6,'(80(''=''))')
!
   end subroutine printCoreStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine corslv(id,ia,is,nws,last,last2,h,r,sqrt_r,rv)
!  ===================================================================
   use InterpolationModule, only: getInterpolation
!
   use BesselModule, only : IntegrateSphHankelSq
!
   use WriteFunctionModule, only : writeFunction
!
   use ChemElementModule, only : MaxLenOfAtomName
   use ChemElementModule, only : setConfiguration
   use ChemElementModule, only : getZtot, getZcor, getZsem
   use ChemElementModule, only : getNumCoreStates
!
   use AtomModule, only : getLocalAtomNickName
!
   implicit   none
!
   character (len=MaxLenOfAtomName) :: an
   character (len=6), parameter :: sname='corslv'
!
   character (len=4) :: state_string
   character (len=10) :: spin_string
!
   integer, intent(in) :: nws, last, last2
!
   integer, intent(in) :: id,ia,is
   integer :: i
   integer :: j
   integer :: ndeep
   integer :: kappa
   integer :: jmt, num_cores, num_electrons, occ, npt, lpt
!
   real (kind=RealKind), intent(in) :: h
   real (kind=RealKind), intent(in) :: r(last)
   real (kind=RealKind), intent(in) :: sqrt_r(0:last)
   real (kind=RealKind), intent(in) :: rv(last)
!
   real (kind=RealKind) :: fac1
   real (kind=RealKind) :: gnrm, norm_frac, h2nrm
   real (kind=RealKind) :: ftmp_int, dftmp_int, err
!
   integer, parameter :: solver = 1
   integer, parameter :: nr = 4
!
!  ===================================================================
!     this routine drives the radial equations solution for the
!     core states. it must be called for each component(sublattice)
!     separately.
!
!     calls:  deepst(outws,inws)
!             semden(outws,hank)
!
!     r= radial grid, rv=potential times r
!     ec= core levels guesses: to be updated
!     g,f upper and lower components times r
!     corden= core charge density
!     nc= principal q. nos., lc= orbital q. nos., kc= kappa q.nos.
!     numc= no. of core states, z= atomic no.
!     ndeepz= number of deep core electrons based on zcor
!     nws=top r index
!  ===================================================================
!
   Core(id)%esemv(is,ia)=ZERO
   Core(id)%ecorv(is,ia)=ZERO
!
   if (Core(id)%numc(ia) <= 0) then
      return
   endif
!
   if(solver /= 1) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'solver <> 1 is not implemented',solver)
!     ----------------------------------------------------------------
!ywg  allocate( potc(nr,last) )
!     ----------------------------------------------------------------
!ywg  call fitpot(r,rv,potc,nr,jmt,last)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  call deepst for each core state....................................
!  ===================================================================
   jmt=Core(id)%jcore
   num_cores = 0; num_electrons = 0
   ndeep = 0
   npt = -1; lpt = -1
   LOOP_i: do i = 1, Core(id)%numc(ia)
      fac1=(3-n_spin_pola)*abs(Core(id)%kc(i,ia))
      ndeep=ndeep+fac1
      if(solver.eq.1) then
!        -------------------------------------------------------------
         call deepst(Core(id)%nc(i,ia),Core(id)%lc(i,ia),Core(id)%kc(i,ia), &
                     Core(id)%ec(i,is,ia),rv,r,sqrt_r(0:last),ftmp(1:last), &
                     dftmp(1:last),h2nrm,h,Core(id)%ztotss(ia),nws,last)
!        -------------------------------------------------------------
!
         if (ndeep > ndeepz) then
!           ----------------------------------------------------------
            call semcst(Core(id)%nc(i,ia),Core(id)%lc(i,ia),Core(id)%kc(i,ia),  &
                        Core(id)%ec(i,is,ia),rv,r,sqrt_r(0:last),ftmp(1:last2), &
                        dftmp(1:last2),h2nrm,h,Core(id)%ztotss(ia),jmt,nws,last2)
!           ----------------------------------------------------------
         endif
      else
         if(isNonRelativistic) then
            kappa=Core(id)%lc(i,ia)
         else
            kappa=Core(id)%kc(i,ia)
         endif
!        -------------------------------------------------------------
!ywg     call srcore(nrelc.eq.0,ecore(i),ftmp(1),nc(i),kappa,potc,nr, &
!ywg                 r,nitmax,last2,1.d0)
!        -------------------------------------------------------------
!ywg     deallocate( potc )
      endif
!
!     ================================================================
!     Check if the core state energy is above the bottom of the valence contour
!     ================================================================
      if (Core(id)%ec(i,is,ia) > evbot) then
!        if (is == 2) then
!           call ErrorHandler(sname,                                  &
!                             'The bottom of the contour is between the spin up/down states of core electron', &
!                             Core(id)%ec(i,1,ia),Core(id)%ec(i,2,ia))
!        endif
         cycle LOOP_i
      else
         num_cores = i
         num_electrons = ndeep
      endif
!
!     ================================================================
!     normalize the wavefunctions
!     ================================================================
      ftmp(0) = ZERO
      do j=1,last2
         ftmp(j)=ftmp(j)/r(j)  ! get rid of factor r
      enddo
!     ----------------------------------------------------------------
!     call calIntegration(0,last2,r(1:last2),ftmp(1:last2),qmp(1:last2))
!     gnrm=ONE/(qmp(last2)+HALF*ftmp(1)*r(1))
!     ----------------------------------------------------------------
      call calIntegration(last2+1,sqrt_r(0:last2),ftmp(0:last2),qmp(0:last2),3)
!     ----------------------------------------------------------------
#ifdef CORE_NORM2INFINITY
      call IntegrateSphHankelSq(Core(id)%lc(i,ia),r(last2),Core(id)%ec(i,is,ia),norm_frac)
!     ----------------------------------------------------------------
      if (print_level >= 0) then
         write(6,'(a,3i4,2(a,d15.8),a,2d15.8)')'nc, lc, kc = ',       &
               Core(id)%nc(i,ia),Core(id)%lc(i,ia),Core(id)%kc(i,ia), &
               ', Int[rho] beyond Rc = ', PI4*norm_frac*h2nrm,        &
               ', Int[Rho] within Rc = ',TWO*PI4*qmp(last2),          &
               ', norm_frac, h2nrm = ',norm_frac,h2nrm
      endif
      gnrm=ONE/(TWO*PI4*qmp(last2)+PI4*norm_frac*h2nrm) ! Normalized to R = infinity
#else
      gnrm=ONE/(TWO*PI4*qmp(last2)) ! Normalized to R = r(last2)
#endif
      do j=1,last2
         ftmp(j)=ftmp(j)*gnrm/r(j)    ! get rid off another factor r so that
                                      ! at this stage, ftmp is just density
      enddo
      do j=1,last2
         dftmp(j)=dftmp(j)*gnrm/r(j)**2 ! get rid off factor r**2 so that
                                        ! at this stage, dftmp is just density
                                        ! derivative
      enddo
      if (print_level > 0 .and. (npt /= Core(id)%nc(i,ia) .or. lpt /= Core(id)%lc(i,ia))) then
         if (Core(id)%lc(i,ia) > 3) then
!           ----------------------------------------------------------
            call ErrorHandler(sname,'Core state l > 3',Core(id)%lc(i,ia))
!           ----------------------------------------------------------
         endif
         if(isNonRelativistic) then
            occ = (3-n_spin_pola)*(2*Core(id)%lc(i,ia)+1)
         else
            occ = (3-n_spin_pola)*abs(Core(id)%kc(i,ia))
         endif
         lpt = Core(id)%lc(i,ia)
         npt = Core(id)%nc(i,ia)
         if (occ < 10) then
            if (Core(id)%lc(i,ia) == 0) then
               write(state_string,'(i1,a,i1)')Core(id)%nc(i,ia),'s',occ
            else if (Core(id)%lc(i,ia) == 1) then
               write(state_string,'(i1,a,i1)')Core(id)%nc(i,ia),'p',occ
            else if (Core(id)%lc(i,ia) == 2) then
               write(state_string,'(i1,a,i1)')Core(id)%nc(i,ia),'d',occ
            else if (Core(id)%lc(i,ia) == 3) then
               write(state_string,'(i1,a,i1)')Core(id)%nc(i,ia),'f',occ
            endif
         else
            if (Core(id)%lc(i,ia) == 0) then
               write(state_string,'(i1,a,i2)')Core(id)%nc(i,ia),'s',occ
            else if (Core(id)%lc(i,ia) == 1) then
               write(state_string,'(i1,a,i2)')Core(id)%nc(i,ia),'p',occ
            else if (Core(id)%lc(i,ia) == 2) then
               write(state_string,'(i1,a,i2)')Core(id)%nc(i,ia),'d',occ
            else if (Core(id)%lc(i,ia) == 3) then
               write(state_string,'(i1,a,i2)')Core(id)%nc(i,ia),'f',occ
            endif
         endif
         if (n_spin_pola == 1) then
            spin_string = '_NonPola'
         else if (is == 1) then
            spin_string = '_SpinUp'
         else
            spin_string = '_SpinDown'
         endif
!        -------------------------------------------------------------
         call writeFunction('Core_'//trim(state_string)//spin_string,last2,r,ftmp(1:),2)
!        -------------------------------------------------------------
      endif
      if(ndeep > ndeepz)then
         do j=1,Core(id)%Grid%jmt
            Core(id)%semden(j,is,ia)=Core(id)%semden(j,is,ia) + fac1*ftmp(j)
         enddo
         do j=1,Core(id)%Grid%jmt
            Core(id)%dsemden(j,is,ia)=Core(id)%dsemden(j,is,ia) + fac1*dftmp(j)
         enddo
         if ( Core(id)%Grid%nmult>1 ) then
            do j=Core(id)%Grid%jmt+1,Core(id)%Grid%jmt+(last2-Core(id)%Grid%jmt)*Core(id)%Grid%nmult
               ftmp_int = getInterpolation(last2,r,ftmp(1:),Core(id)%r_mesh(j),err)
               Core(id)%semden(j,is,ia)=Core(id)%semden(j,is,ia) + fac1*ftmp_int
            enddo
            do j=Core(id)%Grid%jmt+1,Core(id)%Grid%jmt+(last2-Core(id)%Grid%jmt)*Core(id)%Grid%nmult
               dftmp_int = getInterpolation(last2,r,dftmp,Core(id)%r_mesh(j),err)
               Core(id)%dsemden(j,is,ia)=Core(id)%dsemden(j,is,ia) + fac1*dftmp_int
            enddo
         else
            do j=Core(id)%Grid%jmt+1,last2
               ftmp_int = ftmp(j)
               Core(id)%semden(j,is,ia)=Core(id)%semden(j,is,ia) + fac1*ftmp_int
            enddo
            do j=Core(id)%Grid%jmt+1,last2
               dftmp_int = dftmp(j)
               Core(id)%dsemden(j,is,ia)=Core(id)%dsemden(j,is,ia) + fac1*dftmp_int
            enddo
         endif
         Core(id)%esemv(is,ia)=Core(id)%esemv(is,ia)+Core(id)%ec(i,is,ia)*fac1
      else
         do j=1,Core(id)%Grid%jmt
            Core(id)%corden(j,is,ia)= Core(id)%corden(j,is,ia) + fac1*ftmp(j)
         enddo
         do j=1,Core(id)%Grid%jmt
            Core(id)%dcorden(j,is,ia)= Core(id)%dcorden(j,is,ia) + fac1*dftmp(j)
         enddo
         if ( Core(id)%Grid%nmult>1 ) then
            do j=Core(id)%Grid%jmt+1,Core(id)%Grid%jmt+(last2-Core(id)%Grid%jmt)*Core(id)%Grid%nmult
               ftmp_int = getInterpolation(last2,r,ftmp(1:),Core(id)%r_mesh(j),err)
               Core(id)%corden(j,is,ia)=Core(id)%corden(j,is,ia) + fac1*ftmp_int
            enddo
            do j=Core(id)%Grid%jmt+1,Core(id)%Grid%jmt+(last2-Core(id)%Grid%jmt)*Core(id)%Grid%nmult
               dftmp_int = getInterpolation(last2,r,dftmp,Core(id)%r_mesh(j),err)
               Core(id)%dcorden(j,is,ia)=Core(id)%dcorden(j,is,ia) + fac1*dftmp_int
            enddo
         else
            do j=Core(id)%Grid%jmt+1,last2
               ftmp_int = ftmp(j)
               Core(id)%corden(j,is,ia)=Core(id)%corden(j,is,ia) + fac1*ftmp_int
            enddo
            do j=Core(id)%Grid%jmt+1,last2
               dftmp_int = dftmp(j)
               Core(id)%dcorden(j,is,ia)=Core(id)%dcorden(j,is,ia) + fac1*dftmp_int
            enddo
         endif
         Core(id)%ecorv(is,ia)=Core(id)%ecorv(is,ia)+Core(id)%ec(i,is,ia)*fac1
      endif
   enddo LOOP_i
!
   Core(id)%numc_below(ia) = num_cores
   if (num_cores < Core(id)%numc(ia)) then
      an = getLocalAtomNickName(id,ia)
      if (ndeep >= ndeepz) then
         Core(id)%zsemss(ia)=n_spin_pola*(num_electrons-ndeepz)
!        -------------------------------------------------------------
         call setConfiguration(an, Zs = n_spin_pola*(num_electrons-ndeepz))
!        -------------------------------------------------------------
      else
         Core(id)%zcorss(ia)=n_spin_pola*num_electrons
         Core(id)%zsemss(ia)=0
!        -------------------------------------------------------------
         call setConfiguration(an, Zc = n_spin_pola*num_electrons, Zs = 0)
!        -------------------------------------------------------------
      endif
!     Core(id)%numc(ia) = nc
      if (print_level >= 0) then
         write(6,'(a)')   'The electronic configuration has been reset to the follows:'
         write(6,'(a,i3)')'    Number of Deep Core Electrons: ',getZcor(an)
         write(6,'(a,i3)')'    Number of Semi Core Electrons: ',getZsem(an)
!        write(6,'(a,i3)')'    Number of Core Electron States:',getNumCoreStates(an)
      endif
   endif
!
   if(stop_routine.eq.sname) then
      call StopHandler(sname)
   endif
!
   end subroutine corslv
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deepst(nqn,lqn,kqn,en,rv,r,sqrt_r,rf,der_rf,h2nrm,h,z,nws,last)
!  ===================================================================
!
!     core states solver ...............................bg...june 1990
!     energy eigenvalues are found by newton-raphson method matching
!     at the classical inversion point the f/g ratio. the energy
!     derivative of f/g is evaluated analitically: see the book by Rose,
!     chapter 4, pages 163-169, and also jp desclaux code.............
!
!     calls: outws(invals) and inws
!
!     variables explanation:
!
!               nqn:     principal quantum number;
!               lqn:     orbital quantum number;
!               kqn:     kappa quantum number;
!               en:      energy;
!               rv:      potential in rydbergs times r;
!               r:       radial log. grid;
!               rg:      big component times r;
!               rf:      small component times r; 
!                        In return, it is density times r**2
!               der_rf:  derivative of density times r**2
!               h:       exp. step;
!               z:       atomic number;
!               nws:     bounding sphere radius index;
!               c:       speed of light in rydbergs;
!               drg,drf: wavefunctions derivatives times r;
!               gam:     first power for small r expansion;
!               slp:     slope at the origin;
!               dm:      h/720;
!  ===================================================================
   implicit   none
!
   character (len=6), parameter :: sname='deepst'
!
   integer (kind=IntKind), intent(in) :: nqn
   integer (kind=IntKind), intent(in) :: lqn
   integer (kind=IntKind), intent(in) :: kqn
   integer (kind=IntKind), intent(in) :: nws
   integer (kind=IntKind), intent(in) :: last
!
   integer (kind=IntKind) :: iter
   integer (kind=IntKind) :: imm
   integer (kind=IntKind) :: nodes
   integer (kind=IntKind) :: lll
   integer (kind=IntKind) :: invp
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: nmax
!
   real (kind=RealKind), intent(in) :: z
   real (kind=RealKind), intent(out) :: rf(last)
   real (kind=RealKind), intent(out) :: der_rf(last)
   real (kind=RealKind), intent(out) :: h2nrm
   real (kind=RealKind), intent(in) :: rv(last)
   real (kind=RealKind), intent(in) :: r(last)
   real (kind=RealKind), intent(in) :: sqrt_r(0:last)
   real (kind=RealKind), intent(inout) :: en
   real (kind=RealKind), intent(in) :: h
!
   real (kind=RealKind) :: drg(ipdeq2)
   real (kind=RealKind) :: drf(ipdeq2)
   real (kind=RealKind) :: rg(last), der_rg(last)
   real (kind=RealKind) :: dk
   real (kind=RealKind) :: dm
   real (kind=RealKind) :: gam
   real (kind=RealKind) :: sign
   real (kind=RealKind) :: slp
   real (kind=RealKind) :: elim
   real (kind=RealKind) :: rfm
   real (kind=RealKind) :: rgm
   real (kind=RealKind) :: rose
   real (kind=RealKind) :: fnrm
   real (kind=RealKind) :: val
   real (kind=RealKind) :: de
   real (kind=RealKind) :: enew, en_save
!
!  ===================================================================
!  check grid and see if core routine will work on this grid
!  ===================================================================
   if(nws+ipdeq2 > last) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'core cannot work: nws+ipdeq2 > last', &
                        nws+ipdeq2,last)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  initialize quantities
!  ===================================================================
   imm=0
   iter=0
   dk=kqn
   dm=h/720.0d+00
!
!  ===================================================================
!  no. of nodes for the current states big component
!  ===================================================================
   nodes=nqn-lqn
!
!  ===================================================================
!  first power of small r expansion
!  ===================================================================
!  gam=sqrt(dk*dk-FOUR*z*z/(MyLightSpeed*MyLightSpeed))
   gam=sqrt(dk*dk-FOUR*z*z*c2inv)      ! 02/11/18 -Yang
!
!  ===================================================================
!  slope of the wavefcns at the origine
!  ===================================================================
   if(nodes-2*(nodes/2) /= 0) then
      sign= ONE
   else
      sign=-ONE
   endif
   slp=sign*kqn/abs(kqn)
!
!  ===================================================================
!  find the lowest energy limit
!  ===================================================================
   lll=(lqn*(lqn+1))/2
   if (lll.eq.0) then
      elim=-z*z/(0.75d0*nqn*nqn)
   else
      elim=(rv(1)+lll/r(1))/r(1)
!ywg  do j=2,last
      do j=2,nws
         elim=min((rv(j)+lll/r(j))/r(j),elim)
      enddo
   endif
!
!  ===================================================================
!  check potential
!  ===================================================================
   if(elim.ge.0) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'v+l*(l+1)/r**2 always positive',elim)
!     ----------------------------------------------------------------
   endif
!
   if(en <= elim .or. abs(en) < TEN2m8) en=elim*HALF
!
!  ===================================================================
!  routine outws performs the outward integration
!  ===================================================================
   en_save = en  ! Inserted here by Yang Wang on 12/17/2015
   LOOP_iter: do while(iter <= nitmax)
!     ----------------------------------------------------------------
      call outws(invp,rg,rf,der_rg,der_rf,rv,r,en,drg,drf,elim,z,          &
                 gam,slp,imm,lll,dk,dm,nodes,last)
!     ----------------------------------------------------------------
!
!     ================================================================
!     store the inversion point values
!     ================================================================
      rfm=rf(invp)
      rgm=rg(invp)
!
!     ================================================================
!     routine inws performs the inward integration
!     ----------------------------------------------------------------
      call inws(invp,nmax,rg,rf,der_rg,der_rf,rv,r,en,drg,drf,dk,dm,last,imm)
!     ----------------------------------------------------------------
!
!     ================================================================
!     components match
!     ================================================================
      fnrm=rgm/rg(invp)
!
!     ================================================================
!     check first the sign of the big component
!     ================================================================
      if (fnrm.lt.ZERO) then
!        -------------------------------------------------------------
         call ErrorHandler(sname,'wrong big component change, fnrm < 0',fnrm)
!        -------------------------------------------------------------
      endif
!
      do j=invp,nmax
         rg(j)=rg(j)*fnrm
         rf(j)=rf(j)*fnrm
      enddo
!
!     ================================================================
!     Inserted by Yang Wang @12/21/2018
!     ----------------------------------------------------------------
      do j=invp,nmax
         der_rg(j) = der_rg(j)*fnrm
         der_rf(j) = der_rf(j)*fnrm
      enddo
!
!     ================================================================
!     energy derivative of the wvfcns "log. derivative"
!     ================================================================
      tmp(0) = ZERO
      do j=1,nmax
         tmp(j)=rg(j)*rg(j)+rf(j)*rf(j)
      enddo
!     ----------------------------------------------------------------
      call calIntegration(nmax+1,sqrt_r(0:nmax),tmp(0:nmax),qmp(0:nmax),1)
!     ----------------------------------------------------------------
      rose=TWO*qmp(nmax)+h*r(invp)*(rfm*rfm-rf(invp)*rf(invp))*THIRD
!     ----------------------------------------------------------------
!     call calIntegration(nmax,r(1:nmax),tmp(1:nmax),qmp(1:nmax),0)
!     ----------------------------------------------------------------
!     rose=qmp(nmax)+HALF*tmp(1)*r(1)+                                &
!          h*r(invp)*(rfm*rfm-rf(invp)*rf(invp))*THIRD
!
!     ================================================================
!     energy update
!     ================================================================
      de=rg(invp)*(rfm-rf(invp))*MyLightSpeed/rose
      imm=0
      val=abs(de/en)
      if (val > tol) then
         LOOP_val: do
            enew=en+de
!           if(enew >= ZERO) then
            if(enew >= evbot) then ! Modified by Yang Wang on 12/17/2015
!              de=de+HALF   ! It should be de=de*HALF
               de=de*HALF   ! Modified by Yang Wang on 12/17/2015
               val=val*HALF
               if (val > tol) then
                  cycle LOOP_val
               else
!                 ====================================================
!                 just in case the energy becomes zero ...
!                 Instead of terminating the code, I force the code to
!                 exit the loop LOOP_iter. Modified by Yang Wang on 12/17/2015
!                 ----------------------------------------------------
!                 if (print_level >= 0) then
!                    write(6,'(a,3d15.6)')'Forced to exit the iteration loop. de,en,val = ',de,en,val
!                 endif
                  en = enew
                  iter = nitmax + 1
                  exit LOOP_iter
!                 call ErrorHandler(sname,'zero energy.',enew)
!                 ----------------------------------------------------
               endif
            else
               exit LOOP_val
            endif
         enddo LOOP_val
         en=enew
         if (val <= 0.2d+00) imm=1
         iter=iter+1
      else
         exit LOOP_iter
      endif
   enddo LOOP_iter
!  if (iter > nitmax) then
!     ================================================================
!     not converged, too small tolerance or too small nitmax
!        The state is possibly entering the valence band. Modified
!        by Yang Wang on 12/17/2015
!     ================================================================
!     if (print_level >= 0) then
!        -------------------------------------------------------------
!        call WarningHandler(sname,'too many iterations to converge energy',en)
!        -------------------------------------------------------------
!     endif
!  endif
!
!  ===================================================================
!  Inserted by Yang Wang @12/21/2018
!  -------------------------------------------------------------------
   do j=1,nmax
      der_rf(j) = TWO*(rf(j)*(der_rf(j)-rf(j)) + rg(j)*(der_rg(j)-rg(j)))/r(j)
   enddo
   if(nmax.lt.last) then
      do j=nmax+1,last
         der_rf(j)=ZERO
      enddo
   endif
!  ===================================================================
!
   do j=1,nmax
      rf(j)=rf(j)*rf(j)+rg(j)*rg(j)
   enddo    
   if(nmax.lt.last) then
      do j=nmax+1,last
         rf(j)=ZERO
      enddo
   endif
!
   h2nrm = fnrm**2
!
   if (sname == stop_routine) then
      call StopHandler(sname,'Forced to stop')
   endif
!
   end subroutine deepst
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine invals(rg,rf,r,slp,gam,v,z,drg,drf,en,dk)
!  ===================================================================
!
!     initial values for outward integration..........................
!     modified desclaux code.................bg.....june 1990.........
!
!     variables explanation: .........................................
!
!               rg:        big component times r;
!               rf:        small component times r;
!               r:         radial log. grid;
!               slp:       slope at the origine;
!               gam:       1st term power;
!               v:         potential at the first point;
!               z:         atomic number;
!               drg,drf:   derivatives times r;
!               c:         speed of light in rydbergs;
!               en:        energy in rydbergs;
!               dk:        spin angular quantum number;
!  ===================================================================
   implicit   none
!
   character (len=6), parameter :: sname = 'invals'
!
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: k
!
   real (kind=RealKind), intent(in) :: z
   real (kind=RealKind), intent(in) :: v
   real (kind=RealKind), intent(in) :: en
   real (kind=RealKind), intent(out) :: rg(:)
   real (kind=RealKind), intent(out) :: rf(:)
   real (kind=RealKind), intent(out) :: drg(ipdeq2)
   real (kind=RealKind), intent(out) :: drf(ipdeq2)
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: dk
   real (kind=RealKind), intent(in) :: gam
   real (kind=RealKind), intent(in) :: slp
   real (kind=RealKind) :: zalfa
   real (kind=RealKind) :: term1
   real (kind=RealKind) :: term2
   real (kind=RealKind) :: pow1
   real (kind=RealKind) :: pow2
   real (kind=RealKind) :: dm
   real (kind=RealKind) :: ratfg
   real (kind=RealKind) :: sum0
   real (kind=RealKind) :: rfr
   real (kind=RealKind) :: rgr
!
!  ==================================================================
!  zalfa=TWO*z/MyLightSpeed
   zalfa=TWO*z*cinv  ! 02/11/18 -Yang
   term1=-zalfa
!  term2=zalfa/r(1)+(v-en)/MyLightSpeed
   term2=zalfa/r(1)+(v-en)*cinv  ! 02/11/18 -Yang
!
   do j=1,ipdeq2
      rg(j)=ZERO
      rf(j)=ZERO
   enddo
!
   if( dk .le. ZERO ) then
!     ratfg=(dk-gam)/zalfa
      ratfg=(dk-gam)*MyLightSpeed*HALF/z ! 02/11/18 -Yang  Careful!!!
   else
      ratfg=zalfa/(dk+gam)
   endif
!
   rf(ipdeq2)=slp
   rg(ipdeq2)=ratfg*slp
!
   do j=1,ipdeq
      rg(j)=rg(ipdeq2)
      rf(j)=rf(ipdeq2)
      drg(j)=rg(j)*gam
      drf(j)=rf(j)*gam
   enddo
!
   m=1
   do while(m <= 20)
      dm=m+gam
      sum0=dm*dm-dk*dk+term1*term1
      rfr=(MyLightSpeed-term2)*rf(m+ipdeq2-1)
      rgr=term2*rg(m+ipdeq2-1)
      pow1=((dm-dk)*rfr-term1*rgr)/sum0
      pow2=((dm+dk)*rgr+term1*rfr)/sum0
      k=-1
!
      do j=1,ipdeq
         rgr=r(j)**m
         rfr=pow2*rgr
         rgr=pow1*rgr
         if (m.ne.1 .and. abs(rgr/rg(j)).le.tol .and.                 &
             abs(rfr/rf(j)).le.tol) then
            k=1
         endif
         rg(j)=rg(j)+rgr
         rf(j)=rf(j)+rfr
         drg(j)=drg(j)+rgr*dm
         drf(j)=drf(j)+rfr*dm
      enddo
!
      if( k .eq. 1 ) then
         return
      endif
!
      rg(m+ipdeq2)=pow1
      rf(m+ipdeq2)=pow2
      m=m+1
   enddo
!
!  -------------------------------------------------------------------
   call ErrorHandler(sname,'no convergence in small r expansion')
!  -------------------------------------------------------------------
!
   end subroutine invals
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine semcst(nqn,lqn,kqn,en,rv,r,sqrt_r,rf,der_rf,h2nrm,      &
                     h,z,jmt,nws,last2)
!  ===================================================================
!
!     semi-core states solver ..........................bg...june 1990
!     energy eigenvalues are found by newton-raphson method matching
!     at the classical inversion point the f/g ratio. the energy
!     derivative of f/g is evaluated analitically:
!     see the book by Rose, chapter 4, pages 163-169,
!     and also jp desclaux code.
!
!     calls: outws(invals) and inws
!
!     variables explanation:
!
!               nqn:      principal quantum number;
!               lqn:      orbital quantum number;
!               kqn:      kappa quantum number;
!               en:       energy;
!               rv:       potential in rydbergs times r;
!               r:        radial log. grid;
!               rg:       big component;
!               rf:       small component;
!                         In return, it is density times r**2
!               der_rf:   derivative of density times r**2
!               h:        exp. step;
!               z:        atomic number;
!               jmt:      muffin-tin radius index;
!               nws:      bopunding sphere radius index;
!               drg,drf:  wavefunctions derivatives times r;
!               gam:      first power for small r expansion;
!               slp:      slope at the origin;
!               dm:       h/720;
!  ===================================================================
   implicit   none
!
   character (len=6), parameter :: sname = 'semcst'
!
   integer (kind=IntKind), intent(in) :: nqn
   integer (kind=IntKind), intent(in) :: lqn
   integer (kind=IntKind), intent(in) :: kqn
   integer (kind=IntKind), intent(in) :: jmt
   integer (kind=IntKind), intent(in) :: nws
   integer (kind=IntKind), intent(in) :: last2
   integer (kind=IntKind) :: iter
   integer (kind=IntKind) :: nmax
   integer (kind=IntKind) :: invp
   integer (kind=IntKind) :: imm
   integer (kind=IntKind) :: nodes
   integer (kind=IntKind) :: lll
   integer (kind=IntKind) :: j
!
   real (kind=RealKind), intent(inout) :: en
   real (kind=RealKind), intent(in) :: h
   real (kind=RealKind), intent(in) :: z
   real (kind=RealKind), intent(out) :: rf(last2)
   real (kind=RealKind), intent(out) :: der_rf(last2)
   real (kind=RealKind), intent(out) :: h2nrm
   real (kind=RealKind), intent(in) :: rv(last2)
   real (kind=RealKind), intent(in) :: r(last2)
   real (kind=RealKind), intent(in) :: sqrt_r(0:last2)
   real (kind=RealKind) :: rg(last2), der_rg(last2)
   real (kind=RealKind) :: drg(ipdeq*2)
   real (kind=RealKind) :: drf(ipdeq*2)
   real (kind=RealKind) :: fnrm
   real (kind=RealKind) :: dk
   real (kind=RealKind) :: dm
   real (kind=RealKind) :: gam
   real (kind=RealKind) :: sign
   real (kind=RealKind) :: slp
   real (kind=RealKind) :: elim
   real (kind=RealKind) :: rfm
   real (kind=RealKind) :: rgm
   real (kind=RealKind) :: rose
   real (kind=RealKind) :: de
   real (kind=RealKind) :: val
   real (kind=RealKind) :: enew
!
!  ===================================================================
!  initialize quantities
!  ===================================================================
   iter=0
   dk=kqn
   dm=h/720.0d+00
!
!  ===================================================================
!  for semicore states you want to fix the inversion point at
!  the muffin-in, or ASA radius. set imm=1.........................
!  ===================================================================
   invp=jmt
   imm=1
!
!  ===================================================================
!  for MT case nmax=last, but for ASA case nmax=jws.
!  this sets the charge within the proper regions.
!  ===================================================================
   nmax=last2
!
!  ===================================================================
!  no. of nodes for the current states big component
!  ===================================================================
   nodes=nqn-lqn
!
!  ===================================================================
!  first power of small r expansion
!  ===================================================================
!  gam=sqrt(dk*dk-FOUR*z*z/(MyLightSpeed*MyLightSpeed))
   gam=sqrt(dk*dk-FOUR*z*z*c2inv) ! 02/11/18 -Yang
!
!  ===================================================================
!  slope of the wavefcns at the origine
!  ===================================================================
   if(nodes-2*(nodes/2) /= 0) then
      sign= ONE
   else
      sign=-ONE
   endif
   slp=sign*kqn/abs(kqn)
!
!  ===================================================================
!  find the lowest energy limit
!  ===================================================================
   lll=(lqn*(lqn+1))/2
   if (lll /= 0) then
      elim=-z*z/(0.75d+00*nqn*nqn)
   else
      elim=(rv(1)+lll/r(1))/r(1)
!ywg  do j=2,last
      do j=2,nws
         elim=min((rv(j)+lll/r(j))/r(j),elim)
      enddo
   endif
!
!  ===================================================================
!  check potential
!  ===================================================================
   if(elim >= ZERO) then
      call ErrorHandler(sname,'v+l*(l+1)/r**2 always positive')
   endif
!
   if(en <= elim) en=elim*HALF
!
!  ===================================================================
!  routine outws performs the outward integration
   LOOP_while: do while(iter<=nitmax)
!     ----------------------------------------------------------------
      call outws(invp,rg,rf,der_rg,der_rf,rv,r,en,drg,drf,elim,z,     &
                 gam,slp,imm,lll,dk,dm,nodes,last2)
!     ----------------------------------------------------------------
!
!     ================================================================
!     store the inversion point values
!     ================================================================
      rfm=rf(invp)
      rgm=rg(invp)
!
!     ================================================================
!     sets free solution directly to Riccati-Hankel for r*g and r*f.
!     ----------------------------------------------------------------
      call inwhnk(invp,rg,rf,der_rg,der_rf,r,en,drg,drf,dk,last2)
!     ----------------------------------------------------------------
!
!     ================================================================
!     routine inws performs the inward integration
!     ----------------------------------------------------------------
!     call inws(invp,nmax,rg,rf,der_rg,der_rf,rv,r,en,drg,drf,dk,dm,nws,imm)
!     ----------------------------------------------------------------
!
!     ================================================================
!     components match at r(invp)
!     ================================================================
      fnrm=rgm/rg(invp)
!
!     ================================================================
!     check first the sign of the big component
!     ================================================================
      if (fnrm.lt.ZERO) then
!        -------------------------------------------------------------
         call ErrorHandler(sname,'wrong big component change')
!        -------------------------------------------------------------
      endif
!
      do j=invp,last2
         rg(j)=rg(j)*fnrm
         rf(j)=rf(j)*fnrm
      enddo
!
!     ================================================================
!     Inserted by Yang Wang @12/21/2018
!     ----------------------------------------------------------------
      do j=invp,last2
         der_rg(j) = der_rg(j)*fnrm
         der_rf(j) = der_rf(j)*fnrm
      enddo
!     ================================================================
!
!     ================================================================
!     energy derivative of the wvfcns "log. derivative"
!     ================================================================
      tmp(0) = ZERO
      do j=1,last2
         tmp(j)=rg(j)**2+rf(j)**2
      enddo
!     ----------------------------------------------------------------
      call calIntegration(last2+1,sqrt_r(0:last2),tmp(0:last2),qmp(0:last2),1)
!     ----------------------------------------------------------------
      rose=TWO*qmp(last2)+h*r(invp)*(rfm*rfm-rf(invp)*rf(invp))*THIRD
!     ----------------------------------------------------------------
!     call calIntegration(last,r(1:last),tmp(1:last),qmp(1:last),0)
!     ----------------------------------------------------------------
!     rose=qmp(last)+HALF*tmp(1)*r(1)+                                &
!          h*r(invp)*(rfm*rfm-rf(invp)*rf(invp))*THIRD
!
!     ================================================================
!     energy update
!     ================================================================
      de=rg(invp)*(rfm-rf(invp))*MyLightSpeed/rose
      val=abs(de/en)
      if (val <= tol) then
         exit LOOP_while
      endif
      LOOP_while2: do while (val > tol)
         enew=en+de
         if(enew < ZERO) then
            exit LOOP_while2
         endif
!        de=de+HALF  ! It should be de = de*HALF
         de=de*HALF  ! Modified by Yang Wang on 12/17/2015
         val=val*HALF
      enddo LOOP_while2
!
!     if (enew >= ZERO) then   ! modified by Yang Wang on 12/17/2015
      if (enew >= -0.0001d0) then  ! modified by Yang Wang on 12/17/2015
!        =============================================================
!        just in case the energy becomes zero
!        -------------------------------------------------------------
!        call ErrorHandler(sname,'zero energy')
!        -------------------------------------------------------------
         iter = nitmax
      endif
!
!     ================================================================
!     not yet convergence: try again
!     ================================================================
      en=enew
      iter=iter+1
   enddo LOOP_while
!
   if (en >= evbot .and. print_level >= 0) then
!     ================================================================
!     The core state is trated as a valence state
!     ----------------------------------------------------------------
      call WarningHandler(sname,'The energy is entering the valence band',en)
!     ----------------------------------------------------------------
   else if (iter > nitmax .and. print_level >= 0) then
!     ================================================================
!     not converged, too small tolerance or too small nitmax
!     ----------------------------------------------------------------
      call WarningHandler(sname,'The energy is not converged',en)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  Inserted by Yang Wang @12/21/2018
!  -------------------------------------------------------------------
   do j=1,last2
      der_rf(j) = TWO*(rf(j)*(der_rf(j)-rf(j)) + rg(j)*(der_rg(j)-rg(j)))/r(j)
   enddo
!  ===================================================================
!
   do j=1,last2
      rf(j)=rf(j)*rf(j)+rg(j)*rg(j)
   enddo
!
   h2nrm = fnrm**2
!
   end subroutine semcst
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine inwhnk(invp,rg,rf,der_rg,der_rf,r,en,drg,drf,dk,last2)
!  ===================================================================
!
!     drg,drf:   derivative of dp and dq times r;
!     c      :   speed of light.
!
!     free solutions ....... riccati-hankel functions
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: last2
   integer (kind=IntKind) :: invp
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: ll
   integer (kind=IntKind) :: lp
   integer (kind=IntKind) :: n
!
   real (kind=RealKind) :: kappa
   real (kind=RealKind) :: rg(:)
   real (kind=RealKind) :: rf(:)
   real (kind=RealKind) :: der_rg(:)
   real (kind=RealKind) :: der_rf(:)
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind) :: en
   real (kind=RealKind) :: drf(ipdeq2)
   real (kind=RealKind) :: drg(ipdeq2)
   real (kind=RealKind), intent(in) :: dk
!
   real (kind=RealKind) :: mass
   real (kind=RealKind) :: ak
   real (kind=RealKind) :: sgnk
   real (kind=RealKind) :: factor
!
   complex (kind=CmplxKind) :: bh(10)
   complex (kind=CmplxKind) :: dh(10)
   complex (kind=CmplxKind) :: x
   complex (kind=CmplxKind) :: p
   complex (kind=CmplxKind) :: mil
!
   mass=(ONE + en*c2inv)
   p=sqrt(cmplx(en*mass,ZERO,CmplxKind) )
!
   kappa=dk
   ak=abs(dk)
   sgnk=dk/ak
   factor= sgnk*sqrt( -en/mass )*cinv
!
!  ===================================================================
!  officially this is L+1
!  ===================================================================
   l = ak + (1+sgnk)/2
   lp=l - sgnk
   mil=sqrtm1**(l-1)
!
!  ===================================================================
!  maximum l needed for this kappa
!  ===================================================================
   ll=max(l,lp)
!
   do n=last2,invp,-1
      x=p*r(n)
!     ----------------------------------------------------------------
      call richnk(ll,x,bh,dh)
!     ----------------------------------------------------------------
!
!     ================================================================
!     NOTE: working on log-mesh so need extra factor of r for rdg
!           and drf.
!     ================================================================
      rg(n)=SQRTm1*mil*bh(l)
      rf(n)=-factor*mil*bh(lp)
!
!     ================================================================
!     Inserted by Yang Wang @12/21/2018
!     ----------------------------------------------------------------
      der_rg(n) = x*(SQRTm1*mil*dh(l))
      der_rf(n) = -x*factor*mil*dh(lp)
!     ================================================================
   enddo
!
   drg(ipdeq)=x*(SQRTm1*mil*dh(l))
   drf(ipdeq)=-x*factor*mil*dh(lp)
!
   end subroutine inwhnk
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine inws(invp,nmax,rg,rf,der_rg,der_rf,rv,r,en,drg,drf,dk,dm,nws,imm)
!  ===================================================================
!
   implicit   none
!
   character (len=4), parameter :: sname = 'inws'
!
   integer (kind=IntKind), intent(in) :: invp
   integer (kind=IntKind), intent(out) :: nmax
   integer (kind=IntKind), intent(in) :: nws
   integer (kind=IntKind), intent(in) :: imm
   integer (kind=IntKind) :: i,j,itop,l
!
   real (kind=RealKind), intent(in) :: en,dk,dm
   real (kind=RealKind), intent(in) :: r(:),rv(:)
   real (kind=RealKind), intent(out) :: rg(:),rf(:)
   real (kind=RealKind), intent(out) :: der_rg(:),der_rf(:)
   real (kind=RealKind), intent(out) :: drf(ipdeq2),drg(ipdeq2)
   real (kind=RealKind) :: p,pq,dpr,dqr,er,emvoc,dq,dp
!
   real (kind=RealKind), parameter :: twtnth=TWO/TEN
   real (kind=RealKind), parameter :: thrqrt=THREE/FOUR
   real (kind=RealKind), parameter :: twntsvn=27.d0
   real (kind=RealKind), parameter :: fhnd=475.d0
   real (kind=RealKind), parameter :: fhd=502.d0
   real (kind=RealKind), parameter :: sixhnd=600.d0
   real (kind=RealKind), parameter :: coef1=.9462151394422310d+00
   real (kind=RealKind), parameter :: coef2=.0537848605577689d+00
!
!  ===================================================================
!  adams 5 points  diff. eq. solver for dirac equations
!
!     drg,drf enrivatives of rg and rf times r;  c speed of light
!     rv potential times r; r: radial grid; en energy guess
!     coef1=475./502., coef2=27./502.
!
!  at first find the last point where the wavefcns are not zero
!  ===================================================================
   if( imm .ne. 1 ) then
      do j=1,nws,3
         nmax=nws+1-j
         if( (rv(nmax)-en*r(nmax))*r(nmax) .le. sixhnd ) then
            exit
         endif
      enddo
   endif
!
!  ===================================================================
!  initial values for inward integration
!  ===================================================================
!  p=sqrt(-en*(one+en/(MyLightSpeed*MyLightSpeed)))
   p=sqrt(-en*(one+en*c2inv)) ! 02/11/18 -Yang
!  pq=-p/(MyLightSpeed+en/MyLightSpeed)
   pq=-p*cinv/(ONE+en*c2inv) ! 02/11/18 -Yang
   do j = nmax-ipdeq+1, nmax
      rg(j)=exp(-p*r(j))
      rf(j)=pq*rg(j)
!     ================================================================
!     The following if statement has to be used to fool the STUPID(!!!)
!     pgf90 compiler.
!     ================================================================
      if (j == 1) then
         write(6,'(a,i5,2x,d15.8,2x,d15.8)')'  j, rg, rf = ',j,rg(j),rf(j)
      endif
   enddo
   do i=1,ipdeq
      drg(i)=-p*rg(nmax+1-i)*r(nmax+1-i)
      drf(i)=pq*drg(i)
!
!     ================================================================
!     Inserted by Yang Wang @12/21/2018
!     ----------------------------------------------------------------
      der_rf(nmax+1-i) = drf(i)
      der_rg(nmax+1-i) = drg(i)
!     ================================================================
   enddo
   itop=nmax-ipdeq
!
!  ===================================================================
!  solve dirac equations now
!  ===================================================================
   do i=invp,itop
      j=itop-i+invp
!
!     ================================================================
!     5 points predictor
!     ================================================================
      dpr= rg(j+1) - dm*( 2.51d+02*drg(1)-1.274d+03*drg(2)            &
                         +2.616d+03*drg(3)-2.774d+03*drg(4)+1.901d+03*drg(5) )
      dqr= rf(j+1) - dm*( 2.51d+02*drf(1)-1.274d+03*drf(2)            &
                         +2.616d+03*drf(3)-2.774d+03*drf(4)+1.901d+03*drf(5) )
!
!     ================================================================
!     shift derivatives
!     ================================================================
      do l=2,ipdeq
         drg(l-1)=drg(l)
         drf(l-1)=drf(l)
      enddo
!
!     ================================================================
!     dirac equations (log. mesh)
!     ================================================================
      er=en*r(j)
      emvoc=(er-rv(j))/MyLightSpeed
!     emvoc=(er-rv(j))*cinv ! 02/11/18 -Yang  Careful!!!
      drg(ipdeq)=-dk*dpr+(MyLightSpeed*r(j)+emvoc)*dqr
      drf(ipdeq)=dk*dqr-emvoc*dpr
!
!     ================================================================
!     5 points corrector
!     ================================================================
      dp= rg(j+1) - dm*(-1.9d+01*drg(1)+1.06d+02*drg(2)               &
                 -2.64d+02*drg(3)+ 6.46d+02*drg(4)+2.51d+02*drg(5))
      dq= rf(j+1) - dm*(-1.9d+01*drf(1)+1.06d+02*drf(2)               &
                 -2.64d+02*drf(3)+ 6.46d+02*drf(4)+2.51d+02*drf(5))
!
!     ================================================================
!     mixing
!     ================================================================
      dp=coef1*dp+coef2*dpr
      dq=coef1*dq+coef2*dqr
      rg(j)=dp
      rf(j)=dq
!
!     ================================================================
!     Inserted by Yang Wang @12/21/2018
!     ----------------------------------------------------------------
      der_rf(j) = drf(5)
      der_rg(j) = drg(5)
!     ================================================================
!
!     ================================================================
!     update derivative
!     ================================================================
      drg(ipdeq)=-dk*dp+(MyLightSpeed*r(j)+emvoc)*dq
      drf(ipdeq)=dk*dq-emvoc*dp
   enddo
!
   end subroutine inws
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine outws(invp,rg,rf,der_rg,der_rf,rv,r,en,drg,drf,elim,z,  &
                    gam,slp,imm,lll,dk,dm,nodes,nws)
!  ===================================================================
!
   implicit   none
!
   character (len=5), parameter :: sname='outws'
!
   integer (kind=IntKind), intent(out) :: invp
   integer (kind=IntKind), intent(in) :: imm
   integer (kind=IntKind), intent(in) :: lll
   integer (kind=IntKind), intent(in) :: nodes
   integer (kind=IntKind), intent(in) :: nws
   integer (kind=IntKind) :: nd
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: j
!
   real (kind=RealKind), intent(in) :: z
   real (kind=RealKind), intent(inout) :: en
   real (kind=RealKind), intent(in) :: elim
   real (kind=RealKind), intent(in) :: gam
   real (kind=RealKind), intent(in) :: slp
   real (kind=RealKind), intent(in) :: dk
   real (kind=RealKind), intent(in) :: dm
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: rv(:)
   real (kind=RealKind), intent(out) :: rg(:)
   real (kind=RealKind), intent(out) :: rf(:)
   real (kind=RealKind), intent(out) :: der_rg(:)
   real (kind=RealKind), intent(out) :: der_rf(:)
   real (kind=RealKind), intent(out) :: drf(ipdeq2)
   real (kind=RealKind), intent(out) :: drg(ipdeq2)
   real (kind=RealKind) :: vor
   real (kind=RealKind) :: rpow
   real (kind=RealKind) :: emvoc
   real (kind=RealKind) :: dp
   real (kind=RealKind) :: dq
   real (kind=RealKind) :: dqr
   real (kind=RealKind) :: er
   real (kind=RealKind) :: dpr
!  real (kind=RealKind) :: drfdr(nws)
   real (kind=RealKind), parameter :: twnsvn=27.d0
   real (kind=RealKind), parameter :: fhnd=475.d0
   real (kind=RealKind), parameter :: fnd=502.d0
   real (kind=RealKind), parameter :: coef1=.9462151394422310d+00
   real (kind=RealKind), parameter :: coef2=.0537848605577689d+00
!
!  ===================================================================
!     outward solution
!     adams 5 points  diff. eq. solver for dirac equations
!     drg,drf derivatives of rg and rf times r;  c speed of light
!     rv potential times r; r: radial grid; en energy guess
!     coef1=475./502., coef2=27./502.
!  ===================================================================
!
   nd=0
   vor=rv(1)/r(1)
!
!  ===================================================================
!  find classical inversion point
!  ===================================================================
   LOOP_do : do
      if( imm /= 1 ) then
         LOOP_j1 : do j=ipdeq+2,nws,2
            invp=nws+1-j
            if( ((rv(invp)+lll/r(invp))/r(invp)-en) <= ZERO ) then
               exit LOOP_j1
            endif
         enddo LOOP_j1
         if( invp <= ipdeq ) then
            en=en*HALF
            if( en < -tol .and. nd <= nodes ) then
               cycle LOOP_do
            else
!              -------------------------------------------------------
               call ErrorHandler(sname,'screwed potential')
!              -------------------------------------------------------
            endif
         endif
      endif
!
!     ================================================================
!     initial values for outward integration
!     ----------------------------------------------------------------
      call invals(rg,rf,r,slp,gam,vor,z,drg,drf,en,dk)
!     ----------------------------------------------------------------
!
      nd=1
      do j=1,ipdeq
         rpow=r(j)**gam
         if( j .ne. 1 ) then
            if( rg(j-1) .ne. ZERO ) then
               if( (rg(j)/rg(j-1)) .le. ZERO ) then
                  nd=nd+1
               endif
            endif
         endif
         rg(j)=rg(j)*rpow
         rf(j)=rf(j)*rpow
         drg(j)=drg(j)*rpow
         drf(j)=drf(j)*rpow
!
!        =============================================================
!        Inserted by Yang Wang @12/21/2018
!        -------------------------------------------------------------
         der_rf(j) = drf(j)
         der_rg(j) = drg(j)
!        =============================================================
      enddo
!
!     ================================================================
!     check consistence of signs of big and small component
!     ================================================================
      k=-1+2*(nodes-2*(nodes/2))
!
      if( (rg(1)*k) <= ZERO .or. (k*dk*rf(1)) < ZERO ) then
!        -------------------------------------------------------------
         call ErrorHandler(sname,'errors in small r expansion')
!        -------------------------------------------------------------
      endif
!
!     ================================================================
!     solve dirac eqs. now
!     ================================================================
      LOOP_j2 : do j=ipdeq+1,invp
!        =============================================================
!        5 points predictor
!        =============================================================
         dpr=rg(j-1)+dm*( 2.51d+02*drg(1)-1.274d+03*drg(2)            &
                         +2.616d+03*drg(3)-2.774d+03*drg(4)+1.901d+03*drg(5) )
         dqr=rf(j-1)+dm*( 2.51d+02*drf(1)-1.274d+03*drf(2)            &
                         +2.616d+03*drf(3)-2.774d+03*drf(4)+1.901d+03*drf(5) )
!
!        =============================================================
!        shift derivatives
!        =============================================================
         do l=2,ipdeq
            drg(l-1)=drg(l)
            drf(l-1)=drf(l)
         enddo
!
!        =============================================================
!        dirac equations (log. mesh)
!        =============================================================
         er=en*r(j)
         emvoc=(er-rv(j))/MyLightSpeed
!        emvoc=(er-rv(j))*cinv ! 02/11/18 -Yang  Careful!!!
         drg(ipdeq)=-dk*dpr+(MyLightSpeed*r(j)+emvoc)*dqr
         drf(ipdeq)=dk*dqr-emvoc*dpr
!
!        =============================================================
!        5 points corrector
!        =============================================================
         dp=rg(j-1)+dm*(-1.9d+01*drg(1)+1.06d+02*drg(2)               &
                        -2.64d+02*drg(3)+6.46d+02*drg(4)+2.51d+02*drg(5) )
         dq=rf(j-1)+dm*(-1.9d+01*drf(1)+1.06d+02*drf(2)               &
                        -2.64d+02*drf(3)+6.46d+02*drf(4)+2.51d+02*drf(5) )
!
!        =============================================================
!        mixing
!        =============================================================
         dp=coef1*dp+coef2*dpr
         dq=coef1*dq+coef2*dqr
         rg(j)=dp
         rf(j)=dq
!
!        =============================================================
!        Inserted by Yang Wang @12/21/2018
!        -------------------------------------------------------------
         der_rf(j) = drf(5)
         der_rg(j) = drg(5)
!        =============================================================
!
!        =============================================================
!        update derivative
!        =============================================================
         drg(ipdeq)=-dk*dp+(MyLightSpeed*r(j)+emvoc)*dq
         drf(ipdeq)=dk*dq-emvoc*dp
!
!        =============================================================
!        check number of nodes
!        =============================================================
         if( rg(j-1) /= ZERO .and. rg(j)/rg(j-1) <= ZERO ) then
            nd=nd+1
         endif
         if( nd .gt. nodes ) then
!            =========================================================
!            if no. of nodes is too big decrease the energy and start again
!            =========================================================
             en=1.2d+00*en
             if( en <= elim ) then
                 exit LOOP_do
             else
                 cycle LOOP_do
             endif
         endif
      enddo LOOP_j2
!
!     ================================================================
!     if no. of nodes is too small increase the energy and start again
!     ================================================================
      if( nd == nodes ) then
!        if (print_level >= 0) then
!           call derv5(rf,drfdr,r,invp)
!           write(6,'(/,a)')'================================================='
!           write(6,'(a,2i5)')'Core state wave function and derivatives with nodes, invp = ',nodes,invp
!           write(6,'(  a)')'-------------------------------------------------'
!           do j=1,invp
!              write(6,'(f13.8,2x,3d16.8)')r(j),rf(j),der_rf(j),drfdr(j)*r(j)
!           enddo
!           write(6,'(a,/)')'================================================='
!        endif
         return
      endif
!
      en=0.8d+00*en
!
      if( en .lt. -tol ) then
         cycle LOOP_do
      endif
!
!     ================================================================
!     this energy guess has become too high
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'not enough nodes')
!     ----------------------------------------------------------------
   enddo LOOP_do
!
!  ===================================================================
!  this energy guess has become too low
!  -------------------------------------------------------------------
   call ErrorHandler(sname,'too many nodes')
!  -------------------------------------------------------------------
!
   end subroutine outws
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine richnk (nn,y,bh,dh)
!  ===================================================================
!
!     calculates the riccati-bessel functions hl.
!     this version calculates the derivatives.
!  ===================================================================
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: nn
   integer (kind=IntKind) :: l
!
   complex (kind=CmplxKind), intent(out) :: bh(nn)
   complex (kind=CmplxKind), intent(out) :: dh(nn)
   complex (kind=CmplxKind), intent(in) :: y
!
   complex (kind=CmplxKind) :: x
   complex (kind=CmplxKind) :: xinv
   complex (kind=CmplxKind) :: xpon
!
   x=y
   if (abs(x)==ZERO) then
!     ----------------------------------------------------------------
      call ErrorHandler('richnk','x = 0',x)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  recursion relations
!  ===================================================================
   xpon=exp( sqrtm1*x )
   xinv=one/x
   bh(1)=-sqrtm1
   bh(2)= - one - sqrtm1*xinv
   dh(1)= one
   dh(2)= sqrtm1 * ( xinv - ( sqrtm1 + x ) ) * xinv
!
!  ===================================================================
!  flm=2l+1 for real l=0,1,2,3, ...  quantum number
!  ===================================================================
   do l=3,nn
      bh(l) = (2*l-3)*bh(l-1)*xinv - bh(l-2)
      dh(l) = bh(l-1) - (l-1)*bh(l)*xinv
   enddo
   do l=1,nn
      bh(l)=bh(l)*xpon
      dh(l)=dh(l)*xpon
   enddo
!
   if( abs(x) <= 0.01d0 ) then
!     ================================================================
!     power-series for j if abs(x) is smaller than 0.9. trouble!!!
!     ----------------------------------------------------------------
!     call ErrorHandler('richnk','small argument',x)  ! commented out by Yang
!     Wang on 12/17/2015
!     ----------------------------------------------------------------
   endif
!
   end subroutine richnk
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDCD0(id,ia,is) result(dcdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
!
   real (kind=RealKind), pointer :: dcdp(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreDensity','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getDeepCoreDensity','Invalid local species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getDeepCoreDensity','Invalid spin index',is)
   endif
!
!  n = Core(id)%rsize
!  dcdp => Core(id)%corden(1:n,is,ia)
   dcdp => Core(id)%corden(:,is,ia)
!
   end function getDCD0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDCD1(id,ia) result(dcdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), pointer :: dcdp(:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreDensity','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getDeepCoreDensity','Invalid local species index',ia)
   endif
!
!  n = Core(id)%rsize
   dcdp => Core(id)%corden(:,:,ia)
!
   end function getDCD1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDCD2(id) result(dcdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), pointer :: dcdp(:,:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreDensity','Invalid local atom index',id)
   endif
!
!  n = Core(id)%rsize
   dcdp => Core(id)%corden(:,:,:)
!
   end function getDCD2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDCDDer0(id,ia,is) result(dcdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
!
   real (kind=RealKind), pointer :: dcdp(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreDenDeri','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getDeepCoreDenDeri','Invalid local species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getDeepCoreDenDeri','Invalid spin index',is)
   endif
!
!  n = Core(id)%rsize
   dcdp => Core(id)%dcorden(:,is,ia)
!
   end function getDCDDer0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDCDDer1(id,ia) result(dcdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), pointer :: dcdp(:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreDenDeri','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getDeepCoreDenDeri','Invalid local species index',ia)
   endif
!
!  n = Core(id)%rsize
   dcdp => Core(id)%dcorden(:,:,ia)
!
   end function getDCDDer1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDCDDer2(id) result(dcdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), pointer :: dcdp(:,:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreDenDeri','Invalid local atom index',id)
   endif
!
!  n = Core(id)%rsize
   dcdp => Core(id)%dcorden(:,:,:)
!
   end function getDCDDer2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSCD0(id,ia,is) result(scdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
!
   real (kind=RealKind), pointer :: scdp(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreDensity','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getSemiCoreDensity','Invalid local species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSemiCoreDensity','Invalid spin index',is)
   endif
!
!  n = Core(id)%rsize
!  scdp => Core(id)%semden(1:n,is,ia)
   scdp => Core(id)%semden(:,is,ia)
!
   end function getSCD0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSCD1(id,ia) result(scdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), pointer :: scdp(:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreDensity','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getSemiCoreDensity','Invalid local species index',ia)
   endif
!
!  n = Core(id)%rsize
   scdp => Core(id)%semden(:,:,ia)
!
   end function getSCD1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSCD2(id) result(scdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), pointer :: scdp(:,:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreDensity','Invalid local atom index',id)
   endif
!
!  n = Core(id)%rsize
   scdp => Core(id)%semden(:,:,:)
!
   end function getSCD2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSCDDer0(id,ia,is) result(scdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
!
   real (kind=RealKind), pointer :: scdp(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreDenDeri','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getSemiCoreDenDeri','Invalid local species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSemiCoreDenDeri','Invalid spin index',is)
   endif
!
!  n = Core(id)%rsize
   scdp => Core(id)%dsemden(:,is,ia)
!
   end function getSCDDer0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSCDDer1(id,ia) result(scdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), pointer :: scdp(:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreDenDeri','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getSemiCoreDenDeri','Invalid local species index',ia)
   endif
!
!  n = Core(id)%rsize
   scdp => Core(id)%dsemden(:,:,ia)
!
   end function getSCDDer1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSCDDer2(id) result(scdp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), pointer :: scdp(:,:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreDenDeri','Invalid local atom index',id)
   endif
!
!  n = Core(id)%rsize
   scdp => Core(id)%dsemden(:,:,:)
!
   end function getSCDDer2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDeepCoreE0(id,is) result(dce)
!  ===================================================================
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is
!
   integer (kind=IntKind) :: ia
!
   real (kind=RealKind) :: dce
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreEnergy','Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getDeepCoreEnergy','Invalid spin index',is)
   endif
!
   dce = ZERO
   do ia = 1, Core(id)%NumSpecies
      dce = dce + Core(id)%ecorv(is,ia)*getLocalSpeciesContent(id,ia)
   enddo
!
   end function getDeepCoreE0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDeepCoreE1(id,ia,is) result(dce)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
!
   real (kind=RealKind) :: dce
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreEnergy','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getDeepCoreEnergy','Invalid local species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getDeepCoreEnergy','Invalid spin index',is)
   endif
!
   dce = Core(id)%ecorv(is,ia)
!
   end function getDeepCoreE1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSemiCoreE0(id,is) result(sce)
!  ===================================================================
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is
!
   integer (kind=IntKind) :: ia
!
   real (kind=RealKind) :: sce
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreEnergy','Invalid local atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSemiCoreEnergy','Invalid spin index',is)
   endif
!
   sce = ZERO
   do ia = 1, Core(id)%NumSpecies
      sce = sce + Core(id)%esemv(is,ia)*getLocalSpeciesContent(id,ia)
   enddo
!
   end function getSemiCoreE0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSemiCoreE1(id,ia,is) result(sce)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
!
   real (kind=RealKind) :: sce
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreEnergy','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getSemiCoreEnergy','Invalid local species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSemiCoreEnergy','Invalid spin index',is)
   endif
!
   sce = Core(id)%esemv(is,ia)
!
   end function getSemiCoreE1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreVPC0(id) result(wsc)
!  ===================================================================
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: ia
!
   real (kind=RealKind) :: wsc
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCoreVPCharge','Invalid local atom index',id)
   endif
!
   wsc = ZERO
   do ia = 1, Core(id)%NumSpecies
      wsc = wsc + Core(id)%qcpsc_ws(ia)*getLocalSpeciesContent(id,ia)
   enddo
!
   end function getCoreVPC0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreVPC1(id,ia) result(wsc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind) :: wsc
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCoreVPCharge','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getCoreVPCharge','Invalid local species index',ia)
   endif
!
   wsc = Core(id)%qcpsc_ws(ia)
!
   end function getCoreVPC1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreMTC0(id) result(mtc)
!  ===================================================================
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: ia
!
   real (kind=RealKind) :: mtc
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCoreMTCharge','Invalid local atom index',id)
   endif
!
   mtc = ZERO
   do ia = 1, Core(id)%NumSpecies
      mtc = mtc + Core(id)%qcpsc_mt(ia)*getLocalSpeciesContent(id,ia)
   enddo
!
   end function getCoreMTC0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreMTC1(id,ia) result(mtc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind) :: mtc
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCoreMTCharge','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getCoreMTCharge','Invalid local species index',ia)
   endif
!
   mtc = Core(id)%qcpsc_mt(ia)
!
   end function getCoreMTC1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreVPM0(id) result(wsm)
!  ===================================================================
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: ia
!
   real (kind=RealKind) :: wsm
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCoreVPMoment','Invalid local atom index',id)
   endif
!
   wsm = ZERO
   do ia = 1, Core(id)%NumSpecies
      wsm = wsm + Core(id)%mcpsc_ws(ia)*getLocalSpeciesContent(id,ia)
   enddo
!
   end function getCoreVPM0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreVPM1(id,ia) result(wsm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind) :: wsm
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCoreVPMoment','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getCoreVPMoment','Invalid local species index',ia)
   endif
!
   wsm = Core(id)%mcpsc_ws(ia)
!
   end function getCoreVPM1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreMTM0(id) result(mtm)
!  ===================================================================
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: ia
!
   real (kind=RealKind) :: mtm
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCoreMTMoment','Invalid local atom index',id)
   endif
!
   mtm = ZERO
   do ia = 1, Core(id)%NumSpecies
      mtm = mtm + Core(id)%mcpsc_mt(ia)*getLocalSpeciesContent(id,ia)
   enddo
!
   end function getCoreMTM0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreMTM1(id,ia) result(mtm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind) :: mtm
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCoreMTMoment','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getCoreMTMoment','Invalid local species index',ia)
   endif
!
   mtm = Core(id)%mcpsc_mt(ia)
!
   end function getCoreMTM1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterstitialSemiCoreDensity() result(r0)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: r0
!
   r0 = rhoint_semi
!
   end function getInterstitialSemiCoreDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterstitialDeepCoreDensity() result(r0)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: r0
!
   r0 = rhoint_deep
!
   end function getInterstitialDeepCoreDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printCoreDensity(id,derivative)
!  ===================================================================
   use MathParamModule, only: Ten2m8
   use InterpolationModule, only : FitInterp
   use Atom2ProcModule, only : getGlobalIndex
   use WriteFunctionModule, only : writeFunction
!
   implicit none
!
   character (len=6) :: denFlag
   character (len=2) :: specFlag
   character(len=20) :: file_den
   character(len=20) :: sname = 'printCoreDensity'
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: ia
!
   logical, intent(in), optional :: derivative
   logical :: deriv
!
   real (kind=RealKind), pointer :: den(:,:), der_den(:,:), r_mesh(:)
   real (kind=RealKind), allocatable :: total_den(:,:)
   real (kind=RealKind), allocatable :: copy_der_den(:,:)
   real (kind=RealKind), allocatable :: int_der_den(:,:), sqrt_r(:)
   real (kind=RealKind) :: dps, den0
!
   integer (kind=IntKind) :: funit, nrs, is, ir
!
   if (present(derivative)) then
      deriv = derivative
   else
      deriv = .false.
   endif
!
   funit = id*100
   write(denFlag,'(i6)')100000+getGlobalIndex(id)
   denFlag(1:1) = 'n'
   nrs = Core(id)%rsize
   r_mesh => Core(id)%r_mesh
   allocate(total_den(nrs,n_spin_pola))
   if (deriv) then
      allocate(copy_der_den(0:nrs,n_spin_pola))
      allocate(int_der_den(0:nrs,n_spin_pola))
      allocate(sqrt_r(0:nrs))
      sqrt_r(0) = ZERO
      do ir = 1,nrs
         sqrt_r(ir) = sqrt(r_mesh(ir))
      enddo
   endif
   do ia = 1, Core(id)%NumSpecies
      write(specFlag,'(i2)')10+ia
      specFlag(1:1) = 'c'
!
      den => Core(id)%corden(:,:,ia)
      funit = funit + 1
      file_den = 'DeepCore'//'_'//denFlag//specFlag
      if (deriv) then
         der_den => Core(id)%dcorden(:,:,ia)
         do is = 1, n_spin_pola
!           copy_der_den(0,is) = ZERO
!           do ir = 1, nrs
!              copy_der_den(ir,is) = der_den(ir,is)
!           enddo
!           ----------------------------------------------------------
!           call FitInterp( 4, sqrt_r(1:4), den(1:4,is), ZERO, den0, dps )
!           call FitInterp( 4, sqrt_r(1:4), copy_der_den(1:4,is),     &
!                           ZERO, copy_der_den(0,is), dps )
!           ----------------------------------------------------------
!           call calIntegration(nrs+1,sqrt_r(0:nrs),copy_der_den(0:nrs,is),&
!                               int_der_den(0:nrs,is),1)
!           ----------------------------------------------------------
!           int_der_den(:,is) = TWO*int_der_den(:,is) + den0
!           ----------------------------------------------------------
            call derv5(den(1:,is),copy_der_den(1:,is),r_mesh,nrs)
!           ----------------------------------------------------------
         enddo
!        -------------------------------------------------------------
         call writeFunction(file_den,nrs,n_spin_pola,r_mesh,den,      &
                            der_den,copy_der_den(1:,:))
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call writeFunction(file_den,nrs,n_spin_pola,r_mesh,den)
!        -------------------------------------------------------------
      endif
!
      den => Core(id)%semden(:,:,ia)
      funit = funit + 1
      file_den = 'SemiCore'//'_'//denFlag//specFlag
      if (deriv) then
         der_den => Core(id)%dsemden(:,:,ia)
         do is = 1, n_spin_pola
!           copy_der_den(0,is) = ZERO
!           do ir = 1, nrs
!              copy_der_den(ir,is) = der_den(ir,is)
!           enddo
!           ----------------------------------------------------------
!           call FitInterp( 4, sqrt_r(1:4), den(1:4,is), ZERO, den0, dps )
!           call FitInterp( 4, sqrt_r(1:4), copy_der_den(1:4,is),     &
!                           ZERO, copy_der_den(0,is), dps )
!           ----------------------------------------------------------
!           call calIntegration(nrs+1,sqrt_r(0:nrs),copy_der_den(0:nrs,is),&
!                               int_der_den(0:nrs,is),1)
!           ----------------------------------------------------------
!           int_der_den(:,is) = TWO*int_der_den(:,is) + den0
!           ----------------------------------------------------------
            call derv5(den(1:,is),copy_der_den(1:,is),r_mesh,nrs)
!           ----------------------------------------------------------
         enddo
!        -------------------------------------------------------------
         call writeFunction(file_den,nrs,n_spin_pola,r_mesh,den,      &
                            der_den,copy_der_den(1:,:))
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call writeFunction(file_den,nrs,n_spin_pola,r_mesh,den)
!        -------------------------------------------------------------
      endif
!
      do is = 1, n_spin_pola
         do ir = 1, nrs
            total_den(ir,is) = Core(id)%semden(ir,is,ia) + Core(id)%corden(ir,is,ia)
         enddo
      enddo
      funit = funit + 1
      file_den = 'TotalCore'//'_'//denFlag//specFlag
      if (deriv) then
         do is = 1, n_spin_pola
            int_der_den(0,is) = ZERO
            do ir = 1, nrs
               int_der_den(ir,is) = Core(id)%dsemden(ir,is,ia) + Core(id)%dcorden(ir,is,ia)
            enddo
!           ----------------------------------------------------------
!           call FitInterp( 4, sqrt_r(1:4), total_den(1:4,is), ZERO, den0, dps )
!           call FitInterp( 4, sqrt_r(1:4), copy_der_den(1:4,is),    &
!                           ZERO, copy_der_den(0,is), dps )
!           ----------------------------------------------------------
!           call calIntegration(nrs+1,sqrt_r(0:nrs),copy_der_den(0:nrs,is),&
!                               int_der_den(0:nrs,is),1)
!           ----------------------------------------------------------
!           int_der_den(:,is) = TWO*int_der_den(:,is) + den0
!           ----------------------------------------------------------
            call derv5(total_den(1:,is),copy_der_den(1:,is),r_mesh,nrs)
!           ----------------------------------------------------------
         enddo
!
!        -------------------------------------------------------------
         call writeFunction(file_den,nrs,n_spin_pola,r_mesh,total_den,&
                            int_der_den(1:,:),copy_der_den(1:,:))
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call writeFunction(file_den,nrs,n_spin_pola,r_mesh,total_den)
!        -------------------------------------------------------------
      endif
   enddo
   deallocate(total_den)
   if (deriv) then
      deallocate(sqrt_r,copy_der_den,int_der_den)
   endif
   nullify(den, der_den)
!
   if (sname == stop_routine) then
      call StopHandler(sname)
   endif
!
   end subroutine printCoreDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDeepCoreKineticEnergy(id,ia,is) result(ke)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
!
   real (kind=RealKind) :: ke
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDeepCoreKineticEnergy','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getDeepCoreKineticEnergy','Invalid local species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getDeepCoreKineticEnergy','Invalid spin index',is)
   endif
!
   ke = Core(id)%core_ke(is,ia)
!
   end function getDeepCoreKineticEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSemiCoreKineticEnergy(id,ia,is) result(ke)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
!
   real (kind=RealKind) :: ke
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSemiCoreKineticEnergy','Invalid local atom index',id)
   else if (ia < 1 .or. ia > Core(id)%NumSpecies) then
      call ErrorHandler('getSemiCoreKineticEnergy','Invalid local species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSemiCoreKineticEnergy','Invalid spin index',is)
   endif
!
   ke = Core(id)%semi_ke(is,ia)
!
   end function getSemiCoreKineticEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateGlobalCoreStatesTable()
!  ===================================================================
   use GroupCommModule, only : GlobalMaxInGroup, GlobalSumInGroup,    &
                               GlobalCollectInGroup
   use Atom2ProcModule, only : getGlobalIndex
!
   implicit none
!
   character (len=2) :: symb
!
   integer (kind=IntKind) :: maxnc, ia, id, ig, ic, nc, lc
!
   maxnc = 0
   do id = 1,LocalNumAtoms
      ia = 1 ! We only consider 1 species per site
      maxnc = max(maxnc, Core(id)%numc(ia))
   enddo
!  -------------------------------------------------------------------
   call GlobalMaxInGroup(GroupID,maxnc)
!  -------------------------------------------------------------------
!
   if (.not.GlobalTable_Allocated) then
      allocate(NumStatesTable(GlobalNumAtoms))
      allocate(TmpDesTable(maxnc,GlobalNumAtoms))
      allocate(DescriptionTable(maxnc,GlobalNumAtoms))
      allocate(SplitTable(maxnc,GlobalNumAtoms))
      GlobalTable_Allocated = .true.
      maxnc_save = maxnc
   else if (maxnc_save < maxnc) then
      deallocate(NumStatesTable,TmpDesTable,DescriptionTable,SplitTable)
      allocate(NumStatesTable(GlobalNumAtoms))
      allocate(TmpDesTable(maxnc,GlobalNumAtoms))
      allocate(DescriptionTable(maxnc,GlobalNumAtoms))
      allocate(SplitTable(maxnc,GlobalNumAtoms))
      maxnc_save = maxnc
   endif
!
   NumStatesTable = 0
   TmpDesTable = 0
   SplitTable = ZERO
   do id = 1,LocalNumAtoms
      ig = getGlobalIndex(id)
      ia = 1 ! We only consider 1 species per site
      NumStatesTable(ig) = Core(id)%numc(ia)
      do ic = 1, Core(id)%numc(ia)
         SplitTable(ic,ig) = Core(id)%ec(ic,n_spin_pola,ia) -         &
                             Core(id)%ec(ic,1,ia)
         TmpDesTable(ic,ig) = Core(id)%nc(ic,ia)*10+Core(id)%lc(ic,ia)
      enddo
   enddo  
   call GlobalSumInGroup(GroupID,NumStatesTable,GlobalNumAtoms)
   call GlobalSumInGroup(GroupID,SplitTable,maxnc,GlobalNumAtoms)
   call GlobalSumInGroup(GroupID,TmpDesTable,maxnc,GlobalNumAtoms)
!
   DescriptionTable = '  '
   do ig = 1, GlobalNumAtoms
      do ic = 1, maxnc
         if (TmpDesTable(ic,ig) > 0) then
            nc = floor(TmpDesTable(ic,ig)/10.0)
            lc = TmpDesTable(ic,ig)-nc*10
            if (lc == 0) then
               write(symb,'(i1,''s'')')nc
            else if (lc == 1) then
               write(symb,'(i1,''p'')')nc
            else if (lc == 2) then
               write(symb,'(i1,''d'')')nc
            else if (lc == 3) then
               write(symb,'(i1,''f'')')nc
            else
               call ErrorHandler('updateGlobalCoreStateTable','Invalid lc',lc)
            endif
            DescriptionTable(ic,ig) = symb
         endif
      enddo
   enddo
!
   end subroutine updateGlobalCoreStatesTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreSplitTable() result(tbl)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: tbl(:,:)
!
   tbl => SplitTable
!
   end function getCoreSplitTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreNumStatesTable() result(tbl)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: tbl(:)
!
   tbl => NumStatesTable
!
   end function getCoreNumStatesTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreDescriptionTable() result(tbl)
!  ===================================================================
   implicit none
!
   character (len=2), pointer :: tbl(:,:)
!
   tbl => DescriptionTable
!
   end function getCoreDescriptionTable
!  ===================================================================
end module CoreStatesModule
