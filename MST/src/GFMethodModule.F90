module GFMethodModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, THREE, TEN, SQRT_PI, PI, PI2, PI4, &
                               CZERO, CONE, TWO, HALF, SQRTm1, Y0, FIVE, THIRD
   use MathParamModule, only : Ten2m2, Ten2m3, Ten2m6, Ten2m4, Ten2m5,&
                       Ten2m7, TEN2m12, TEN2m10, TEN2m8, TEN2m14, TEN2m16
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
   use IntegerFactorsModule, only : lofk, lofj, mofj, m1m, mofk, jofk
   use PublicTypeDefinitionsModule, only : GridStruct, LloydStruct
   use TimerModule, only : getTime
   use MPPModule, only : MyPE, syncAllPEs
!
public :: initGFMethod,              &
          endGFMethod,               &
          printValenceStates,        &
          calValenceStates,          &
          printExchangeParam,        &
          calValenceDOS,             &
          calPartialDOS,             &
          writeSS_DOS,               &
          writeDOS
!
private
!
   interface adjustEnergy
      module procedure adjustEnergy_r, adjustEnergy_c
   end interface
!
   interface addElectroStruct
      module procedure addElectroStruct_r, addElectroStruct_c
   end interface
!
   interface addSpinCantStruct
      module procedure addSpinCantStruct_r, addSpinCantStruct_c
   end interface
!
   type SpinCantStruct
      integer (kind=IntKind) :: NumJs
      integer (kind=IntKind), pointer :: jid(:)
      real (kind=RealKind), pointer :: r0j(:,:)
      real (kind=RealKind), pointer :: pair_exchab(:,:)
!
      real (kind=RealKind) :: stoner(1:3)
      real (kind=RealKind) :: torque(1:3)
      real (kind=RealKind) :: onesite_stoner
      real (kind=RealKind) :: onesite_susab(0:9)
      real (kind=RealKind) :: onesite_exchange
      real (kind=RealKind) :: onesite_exchab(1:9)
   end type SpinCantStruct
!
   type ElectroStruct
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind) :: NumRs
      integer (kind=IntKind) :: jmax
      integer (kind=IntKind) :: kmax
      integer (kind=IntKind) :: indx  ! the starting index of wk_green that is associated with dos_r_jl
      integer (kind=IntKind) :: size  ! size of dos_r_jl
!
      complex (kind=CmplxKind), pointer :: dos(:,:)
      complex (kind=CmplxKind), pointer :: dos_mt(:,:)
      complex (kind=CmplxKind), pointer :: evalsum(:,:)
      complex (kind=CmplxKind), pointer :: dos_r_jl(:,:,:,:)
      complex (kind=CmplxKind), pointer :: der_dos_r_jl(:,:,:,:)
!
      type (SpinCantStruct), pointer :: pSC(:)
      complex (kind=CmplxKind), pointer :: density_matrix(:,:,:,:)
   end type ElectroStruct
!
   type (ElectroStruct), allocatable, target :: LastValue(:)
   type (ElectroStruct), allocatable, target :: IntegrValue(:)
   type (ElectroStruct), allocatable, target :: ssLastValue(:)
   type (ElectroStruct), allocatable, target :: ssIntegrValue(:)
!
   logical :: Initialized = .false.
   logical :: isIterateEfOn = .true.
   logical :: isEfPinning = .false.
!
!  Flags to stop the code after the specific task is completed( input(external) flags )
!
   character (len = 50) :: stop_routine
!
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: NumVacancies
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: jend_max
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: iharris
   integer (kind=IntKind) :: max_print_level
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind), allocatable :: print_level(:)
   integer (kind=IntKind), allocatable :: AtomIndex(:)
!
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: lmax_green(:)
!
   integer (kind=IntKind) :: RECORD_LENGTH
!
   integer (kind=IntKind), parameter :: n_inter = 5 ! order of polynomial
                                                    ! interpolation
!
   integer (kind=IntKind) :: RelativisticFlag
!
   integer (kind=IntKind) :: Fit_Method = 1
   integer (kind=IntKind) :: ie_count, ews_size
   integer (kind=IntKind) :: NumPEsInAGroup, MyPEinAGroup, aGID
   integer (kind=IntKind) :: NumPEsInEGroup, MyPEinEGroup, eGID
!
   real (kind=RealKind) :: chempot
   real (kind=RealKind) :: zvaltss
   real (kind=RealKind), allocatable :: posi(:,:)
   real (kind=RealKind), allocatable :: evec(:,:)
   real (kind=RealKind), allocatable :: exc(:,:) ! exchange splitting energy
                                               ! between spin-up and spin-down states
!
   type RealArray2Struct
      integer (kind=IntKind) :: n1, n2
      real (kind=RealKind), pointer :: rarray2(:,:)
   end type RealArray2Struct
!
   type RealArray3Struct
      integer (kind=IntKind) :: n1, n2, n3
      real (kind=RealKind), pointer :: rarray3(:,:,:)
   end type RealArray3Struct
!
   type (RealArray2Struct), allocatable ::  ssIDOS_out(:)
   type (RealArray3Struct), allocatable ::  SS_dosdata(:)
   type (RealArray3Struct), allocatable ::  dosdata(:)
   real (kind=RealKind), allocatable, target ::  wk_ssIDOS_out(:)
   real (kind=RealKind), allocatable, target ::  wk_SS_dosdata(:)
   real (kind=RealKind), allocatable, target ::  wk_dosdata(:)
!
   complex (kind=CmplxKind), allocatable, target :: wk_green(:)
   complex (kind=CmplxKind), allocatable, target :: wk_dos(:)
   complex (kind=CmplxKind), allocatable, target :: wk_dgreen(:)
!
   logical :: isDensityMatrixNeeded = .false.
   logical :: isSingleSiteCluster = .false.
!
   complex(kind=CmplxKind) :: Lloyd_lastE(2)
!
   logical :: is_Bxyz=.true., rel_B=.false. !xianglin 
   integer (kind=IntKind) :: lmax_kkr_max
!
   logical :: rad_derivative = .false.
!
!  findResonance variables
!================================================================
   type PolesStruct
      integer (kind=IntKind) :: NumElements
      integer (kind=IntKind), allocatable :: NumPoles(:,:)
      real (kind=RealKind), allocatable :: Poles(:,:,:)
      integer (kind=IntKind), allocatable :: NumPoles_plus(:,:)
      complex (kind=CmplxKind), allocatable :: Poles_plus(:,:,:)
   end type PolesStruct
   type (PolesStruct), allocatable :: MatrixPoles(:)
   logical :: isRel, isPole_plus=.true., isPole=.true. !isRel is just a variable to be used. Relativity is controlled by RelativisticFlag
!================================================================
   integer (kind=IntKind) :: nw_Pole = 5  !global variables
   integer (kind=IntKind) :: nw_Pole_new, nwPG
   integer (kind=IntKind) :: nw_Pole_plus = 5  
   integer (kind=IntKind) :: nw_Pole_new_plus, nwPG_plus
!
   integer (kind=IntKind) :: NumCalls_SS = 0
   real (kind=RealKind) :: Timing_SS = ZERO
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initGFMethod(na,gindex,posi_in,                         &
                           lmax_kkr_in,lmax_phi_in,lmax_pot_in,       &
                           lmax_step_in,lmax_green_in,pola,cant,      &
                           istop,iprint,isGGA)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   use ProcMappingModule, only : getNumEsOnMyProc
!
   use OutputModule, only : getStandardOutputLevel, isOutputAtomBased
!
   use MathParamModule, only : ZERO, CZERO
!
   use RadialGridModule, only : getGrid, getNumRmesh
!
   use SystemModule, only : getNumAtoms, getNumVacancies, getLmaxMax
!
   use AtomModule, only : getLocalEvecOld, getLocalNumSpecies, &
                          getLocalAtomicNumber
!
   use KreinModule, only : setLloydStatus, setLloydQ
!
   use SpinRotationModule, only :  initSpinRotation
!
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       isDataStorageExisting, &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark
!
   use ScfDataModule, only : isLSMS
   use ScfDataModule, only : isNonRelativisticValence
   use ScfDataModule, only : isScalarRelativisticValence
   use ScfDataModule, only : isEfIterateOn
   use ScfDataModule, only : isLloyd
   use ScfDataModule, only : isReadEmesh, getEmeshFileName, Harris
   use ScfDataModule, only : getAdaptiveIntegrationMethod
!
   use NeighborModule, only  : getNumNeighbors
!
   use ValenceDensityModule, only : getFermiEnergy
!
   use AdaptIntegrationModule, only : initAdaptIntegration
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: gindex(na)
   integer (kind=IntKind), intent(in) :: lmax_kkr_in(na)
   integer (kind=IntKind), intent(in) :: lmax_phi_in(na)
   integer (kind=IntKind), intent(in) :: lmax_pot_in(na)
   integer (kind=IntKind), intent(in) :: lmax_step_in(na)
   integer (kind=IntKind), intent(in) :: lmax_green_in(na)
   integer (kind=IntKind), intent(in) :: pola, cant
   integer (kind=IntKind), intent(in) :: iprint(na)
   integer (kind=IntKind) :: id, nsize, n, ia, MaxNumSpecies
   integer (kind=IntKind) :: lmax_max, jmax, green_size, green_ind
   integer (kind=IntKind), allocatable :: NumRs(:)
!
   real (kind=RealKind), intent(in) :: posi_in(3,na)
   real (kind=RealKind), pointer :: p1(:)
!
   logical, intent(in), optional :: isGGA
!
   type (GridStruct), pointer :: Grid
!
   if (pola == 1 .or. pola == 2) then
      n_spin_pola = pola
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initGFMethod','Invalid spin polarizing parameter',pola)
!     ----------------------------------------------------------------
   endif
!
   if (cant == 1 .or. cant == 2) then
      n_spin_cant = cant
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initGFMethod','Invalid spin canting parameter',cant)
!     ----------------------------------------------------------------
   endif
!
   if (n_spin_pola < n_spin_cant) then
!     ----------------------------------------------------------------
      call ErrorHandler('initGFMethod','n_spin_cant > n_spin_pola')
!     ----------------------------------------------------------------
   endif
!
   stop_routine = istop
   isIterateEfOn = isEfIterateOn(isEfPinning)
   iharris = Harris
!
   if (present(isGGA)) then
      rad_derivative = isGGA
   else
      rad_derivative = .false.
   endif
!
   if ( isLloyd() ) then
!     ----------------------------------------------------------------
      call ErrorHandler('initGFMethod','isLloyd is true: It is not verified')
!     ----------------------------------------------------------------
      call setLloydStatus(.true.,.false.)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call setLloydStatus(.false.,.false.)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call setLloydQ(CZERO)
!  -------------------------------------------------------------------
!
   if ( cant >1 .or. pola >1) then !xianglin
      is_Bxyz=.true.
   else
      is_Bxyz=.false.
   endif 
   if (isNonRelativisticValence()) then
      RelativisticFlag = 0
   else if (isScalarRelativisticValence()) then
      RelativisticFlag = 1
   else
      RelativisticFlag = 2
   endif
!
!   xianglin
   if (RelativisticFlag == 2 .and. is_Bxyz) then
      n_spin_cant=2
      n_spin_pola=2
      rel_B=.true.
   endif   
!
   aGID = getGroupID('Unit Cell')
   NumPEsInAGroup = getNumPEsInGroup(aGID)
   MyPEinAGroup = getMyPEinGroup(aGID)
   eGID = getGroupID('Energy Mesh')
   NumPEsInEGroup = getNumPEsInGroup(eGID)
   MyPEinEGroup = getMyPEinGroup(eGID)
!
   LocalNumAtoms = na
   GlobalNumAtoms = getNumAtoms()
   NumVacancies = getNumVacancies()
   if (GlobalNumAtoms == NumVacancies) then
      NumVacancies = 0
   endif
!
   allocate( print_level(LocalNumAtoms) )
   allocate( lmax_green(LocalNumAtoms) )
!
   allocate( AtomIndex(LocalNumAtoms) )
   allocate( posi(3,LocalNumAtoms), evec(3,LocalNumAtoms) )
!
   nsize = 0
   MaxNumSpecies = 0
   do id = 1,LocalNumAtoms
      nsize = nsize + n_spin_pola*n_spin_cant*getLocalNumSpecies(id)
      MaxNumSpecies = max(MaxNumSpecies,getLocalNumSpecies(id))
   enddo
   allocate( wk_ssIDOS_out(nsize), ssIDOS_out(LocalNumAtoms) )
   nsize = 0
   do id = 1,LocalNumAtoms
      ssIDOS_out(id)%n1 = n_spin_pola*n_spin_cant
      ssIDOS_out(id)%n2 = getLocalNumSpecies(id)
      n = ssIDOS_out(id)%n1*ssIDOS_out(id)%n2
      p1 => wk_ssIDOS_out(nsize+1:nsize+n)
      ssIDOS_out(id)%rarray2 => aliasArray2_r(p1,ssIDOS_out(id)%n1,ssIDOS_out(id)%n2)
      nsize = nsize + n
   enddo
!
   allocate(lmax_pot(LocalNumAtoms))
   allocate(lmax_phi(LocalNumAtoms))
   allocate(lmax_kkr(LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms))
   allocate(exc(MaxNumSpecies,LocalNumAtoms)); exc = ZERO
   do id = 1,LocalNumAtoms
      lmax_kkr(id) = lmax_kkr_in(id)
      lmax_pot(id) = lmax_pot_in(id)
      lmax_phi(id) = lmax_phi_in(id)
      lmax_step(id) = lmax_step_in(id)
   enddo
!
   jend_max = 0
   do id = 1,LocalNumAtoms
      Grid => getGrid(id)
      jend_max = max(jend_max,Grid%jend)
   enddo
!
   lmax_max = 0
   lmax_kkr_max = 0
   node_print_level = getStandardOutputLevel()
   if (isOutputAtomBased()) then
      if (getMyPEinGroup(aGID) > 0 .or. getMyPEinGroup(eGID) > 0) then
         node_print_level = -1
      else if (.not.isLSMS()) then
         if (getMyPEinGroup(getGroupID('K-Mesh')) > 0) then
            node_print_level = -1
         endif
      endif
   endif
   do id = 1, LocalNumAtoms
      lmax_max=max(lmax_max, lmax_green_in(id))
      lmax_kkr_max=max(lmax_kkr_max,lmax_kkr_in(id))
      print_level(id) = iprint(id)
      max_print_level = max(iprint(id), max_print_level)
      lmax_green(id)  = lmax_green_in(id)
   enddo
   lmax_max = max(lmax_max,getLmaxMax())
!  ------------------------------------------------------------------
   call initIntegerFactors(lmax_max)
!  ------------------------------------------------------------------
!
   do id = 1, LocalNumAtoms
      Grid => getGrid(id)
      AtomIndex(id) = gindex(id)
      posi(1:3,id)=posi_in(1:3,id)
   enddo
!
   Initialized = .true.
   chempot = getFermiEnergy()
!
   do id = 1,LocalNumAtoms
      evec(1:3,id) = getLocalEvecOld(id)
   enddo
!
   call initSpinRotation(LocalNumAtoms,evec)
!
!  ===================================================================
!  setup the data structure
!  ===================================================================
   allocate( NumRs(LocalNumAtoms) )
!
   green_size = 0
   nsize = 0
   do id = 1,LocalNumAtoms
      Grid => getGrid(id)
      NumRs(id) = Grid%jend
      jmax = (lmax_green(id)+1)*(lmax_green(id)+2)/2
      green_size = green_size +                                       &
                   NumRs(id)*jmax*n_spin_cant*n_spin_pola*getLocalNumSpecies(id)
      if (rad_derivative) then
         nsize = max(nsize,2*NumRs(id)*jmax*getLocalNumSpecies(id))
      else
         nsize = max(nsize,NumRs(id)*jmax*getLocalNumSpecies(id))
      endif
   enddo
!
   if (getAdaptiveIntegrationMethod() == 1) then
      isPole_plus = .true.
   else
      isPole_plus = .false.
   endif
!
   if (RelativisticFlag == 2 .and. is_Bxyz) then
      call initAdaptIntegration(nsize+4,eGID,rel_B)
   else
!     -------------------------------------------------------------------
      call initAdaptIntegration(nsize+4,eGID)
!     -------------------------------------------------------------------
   endif
!
   allocate(wk_green(green_size*4),wk_dos((nsize+4*MaxNumSpecies)*n_spin_cant*n_spin_cant+12))
   wk_green = CZERO; wk_dos = CZERO
   if (rad_derivative) then
      allocate(wk_dgreen(green_size*4))
      wk_dgreen = CZERO
   endif
!
   allocate( LastValue(LocalNumAtoms) )
   allocate( IntegrValue(LocalNumAtoms) )
   allocate( ssLastValue(LocalNumAtoms) )
   allocate( ssIntegrValue(LocalNumAtoms) )
!
   green_ind = 1
   do id = 1, LocalNumAtoms
!     ----------------------------------------------------------------
      call setupElectroStruct(id,NumRs(id),lmax_green(id),IntegrValue(id),green_ind)
!     ----------------------------------------------------------------
      call setupElectroStruct(id,NumRs(id),lmax_green(id),LastValue(id),green_ind)
!     ----------------------------------------------------------------
      call setupElectroStruct(id,NumRs(id),lmax_green(id),ssIntegrValue(id),green_ind)
!     ----------------------------------------------------------------
      call setupElectroStruct(id,NumRs(id),lmax_green(id),ssLastValue(id),green_ind)
!     ----------------------------------------------------------------
   enddo
!
   deallocate(NumRs)
!
   if (isLSMS()) then
      isSingleSiteCluster = .true.
      do id = 1, LocalNumAtoms
         if (getNumNeighbors(id) > 0) then
            isSingleSiteCluster = .false.
            exit
         endif
      enddo
   endif
!
   end subroutine initGFMethod
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endGFMethod()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
!
   use SpinRotationModule, only :  endSpinRotation
!
   use AdaptIntegrationModule, only : endAdaptIntegration
!
   implicit none
   integer (kind=IntKind) :: id
!
   do id = 1, LocalNumAtoms
!     ----------------------------------------------------------------
      call cleanElectroStruct(LastValue(id))
      call cleanElectroStruct(IntegrValue(id))
      call cleanElectroStruct(ssLastValue(id))
      call cleanElectroStruct(ssIntegrValue(id))
!     ----------------------------------------------------------------
   enddo
!
   deallocate(LastValue, IntegrValue, ssLastValue, ssIntegrValue)
   deallocate( wk_green, wk_dos )
   if (rad_derivative) then
      deallocate( wk_dgreen )
      rad_derivative = .false.
   endif
!  -------------------------------------------------------------------
   call endAdaptIntegration()
   call endSpinRotation()
   call endIntegerFactors()
!  -------------------------------------------------------------------
!
   deallocate(print_level)
   deallocate(AtomIndex, posi, evec)
   do id = 1, LocalNumAtoms
      nullify(ssIDOS_out(id)%rarray2)
   enddo
   deallocate(ssIDOS_out, wk_ssIDOS_out)
   deallocate(lmax_green, lmax_pot, lmax_phi, lmax_kkr, lmax_step)
   deallocate(exc)
!
   if (allocated(SS_dosdata)) then
      do id = 1, LocalNumAtoms
         nullify(SS_dosdata(id)%rarray3)
      enddo
      deallocate(SS_dosdata,wk_SS_dosdata)
   endif
   if (allocated(dosdata)) then
      do id = 1, LocalNumAtoms
         nullify(dosdata(id)%rarray3)
      enddo
      deallocate(dosdata,wk_dosdata)
   endif
!
   Initialized = .false.
!
   end subroutine endGFMethod
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupElectroStruct(id,NumRs,lmax,ES,green_ind)
!  ===================================================================
   use ScfDataModule, only : isExchangeParamNeeded
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor
   use AtomModule, only : getLocalNumSpecies
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: id, NumRs, lmax
   integer(kind=IntKind), intent(inout), optional :: green_ind
!
   complex (kind=CmplxKind), pointer :: p1(:)
!
   type (ElectroStruct), intent(inout) :: ES
!
   integer(kind=IntKind) :: iend, jmax, kmax, green_size, n, ia
!
   type (NeighborStruct), pointer :: Neighbor
!
   ES%NumSpecies = getLocalNumSpecies(id)
   ES%NumRs = NumRs
   ES%jmax = (lmax+1)*(lmax+2)/2
   ES%kmax = (lmax+1)**2
!
   n = ES%NumSpecies
   iend = ES%NumRs
   jmax = ES%jmax
   green_size = iend*jmax*n_spin_cant*n_spin_pola*n
   ES%indx = green_ind
   ES%size = green_size
   allocate(ES%dos(4,n))
   allocate(ES%dos_mt(4,n))
   allocate(ES%evalsum(4,n))
!
   p1 => wk_green(green_ind:green_ind+green_size-1)
   ES%dos_r_jl => aliasArray4_c( p1, iend, jmax, n_spin_cant*n_spin_pola, n )
   if (rad_derivative) then
      p1 => wk_dgreen(green_ind:green_ind+green_size-1)
      ES%der_dos_r_jl => aliasArray4_c( p1, iend, jmax, n_spin_cant*n_spin_pola, n )
   endif
   green_ind = green_ind + green_size
!
   if (isDensityMatrixNeeded) then
      kmax = ES%kmax
      allocate(ES%density_matrix(1:kmax,1:kmax,1:n_spin_cant*n_spin_pola,1:n))
   else
      nullify(ES%density_matrix)
   endif
!
   if (n_spin_cant == 2) then
      allocate(ES%pSC(n))
      Neighbor => getNeighbor(id)
      do ia = 1, n
         if ( isExchangeParamNeeded() ) then
            ES%pSC(ia)%NumJs = Neighbor%NumAtoms
         else
            ES%pSC(ia)%NumJs = 0
         endif
         if (ES%pSC(ia)%NumJs > 0) then
            allocate(ES%pSC(ia)%pair_exchab(1:9,1:ES%pSC(ia)%NumJs))
            allocate(ES%pSC(ia)%jid(1:ES%pSC(ia)%NumJs))
            allocate(ES%pSC(ia)%r0j(1:3,1:ES%pSC(ia)%NumJs))
         endif
      enddo
   endif
!  -------------------------------------------------------------------
   call zeroElectroStruct(ES)
!  -------------------------------------------------------------------
!
   end subroutine setupElectroStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cleanElectroStruct(ES)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: ia
!
   type(ElectroStruct), intent(inout) :: ES
!
   deallocate(ES%dos, ES%dos_mt, ES%evalsum)
   if ( associated(ES%dos_r_jl) ) then
      ES%indx = 0; ES%size = 0;
      nullify( ES%dos_r_jl )
      if (rad_derivative) then
         nullify( ES%der_dos_r_jl )
      endif
   endif
   if ( associated(ES%density_matrix) ) then
      deallocate( ES%density_matrix ); nullify( ES%density_matrix )
   endif
   if (n_spin_cant == 2) then
      do ia = 1, ES%NumSpecies
         if ( ES%pSC(ia)%NumJs > 0 ) then
            deallocate( ES%pSC(ia)%pair_exchab )
            deallocate( ES%pSC(ia)%jid )
            deallocate( ES%pSC(ia)%r0j )
         endif
         ES%pSC(ia)%NumJs = 0
      enddo
      deallocate( ES%pSC ); nullify( ES%pSC )
   endif
!
   end subroutine cleanElectroStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine zeroElectroStruct(ES)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: ia
!
   type(ElectroStruct), intent(inout) :: ES
!
   if(associated(ES%dos_r_jl)) then
      ES%dos_r_jl = CZERO
      if (rad_derivative) then
         ES%der_dos_r_jl = CZERO
      endif
   endif
!
   if (isDensityMatrixNeeded .and. associated(ES%density_matrix)) then
      ES%density_matrix = CZERO
   endif
!
   ES%dos = CZERO
   ES%dos_mt = CZERO
   ES%evalsum = CZERO
!
   if (n_spin_cant == 2) then
      do ia = 1, ES%NumSpecies
         ES%pSC(ia)%torque(1:3) = ZERO
         ES%pSC(ia)%stoner(1:3) = ZERO
         ES%pSC(ia)%onesite_stoner = ZERO
         ES%pSC(ia)%onesite_susab(0:9) = ZERO
         ES%pSC(ia)%onesite_exchange = ZERO
         if (ES%pSC(ia)%NumJs > 0) then
            ES%pSC(ia)%pair_exchab = ZERO
            ES%pSC(ia)%jid = 0
            ES%pSC(ia)%r0j = ZERO
            ES%pSC(ia)%pair_exchab = ZERO
         endif
      enddo
   endif
!
   end subroutine zeroElectroStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printValenceStates()
!  ===================================================================
   use AtomModule, only : getLocalEvecNew
!
   use ScfDataModule, only : isExchangeParamNeeded
!
   implicit none
!
   integer (kind=Intkind) :: id, j, ia
!
   real (kind=RealKind) :: evec_new(3)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('printValenceStates','valence module not initialized')
!     ----------------------------------------------------------------
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,23x,a)')'**********************************'
   write(6,'( 23x,a )')'* Output from printValenceStates *'
   write(6,'(23x,a,/)')'**********************************'
!
   write(6,'(2x,a,t40,a,f18.11)')'Fermi Energy (Ryd)','=',chempot
   if (n_spin_cant == 2) then
      do id =1, LocalNumAtoms
         do ia = 1, IntegrValue(id)%NumSpecies
            evec_new(1:3) = getLocalEvecNew(id)
            write(6,'(10x,a,t40,''='',2x,3f10.5)') &
                    'LOCAL Moment orientation',evec_new(1:3)
            write(6,'(10x,a,t40,''='',2x,3f10.5)') &
                    'LOCAL torque on the moment',IntegrValue(id)%pSC(ia)%torque(1:3)
            if ( isExchangeParamNeeded() ) then
               write(6,'(10x,a,t40,''='',2x,3f10.5)') &
                    'LOCAL Stoner parameter vec.',IntegrValue(id)%pSC(ia)%stoner(1:3)
               write(6,'(10x,a,t40,''='',2x,1f10.5)') &
                    'LOCAL On-site Stoner param.',IntegrValue(id)%pSC(ia)%onesite_stoner
               write(6,'(10x,a,t40,''='',2x,1f10.5)') &
                    'LOCAL On-site magn. sus. 00',IntegrValue(id)%pSC(ia)%onesite_susab(0)
               write(6,'(10x,a,t40,''='',2x,3f10.5)') &
                    'LOCAL On-site magn. sus. ab',IntegrValue(id)%pSC(ia)%onesite_susab(1:7:3)
               write(6,'(t40,3x,3f10.5)')IntegrValue(id)%pSC(ia)%onesite_susab(2:8:3)
               write(6,'(t40,3x,3f10.5)')IntegrValue(id)%pSC(ia)%onesite_susab(3:9:3)
               write(6,'(10x,a,t40,''='',2x,3f10.5)') &
                    'LOCAL On-site exch. tensor',IntegrValue(id)%pSC(ia)%onesite_exchab(1:7:3)
               write(6,'(t40,3x,3f10.5)')IntegrValue(id)%pSC(ia)%onesite_exchab(2:8:3)
               write(6,'(t40,3x,3f10.5)')IntegrValue(id)%pSC(ia)%onesite_exchab(3:9:3)
               write(6,'(10x,a,t40,''='',2x,1f10.5)') &
                    'LOCAL On-site exchange param',IntegrValue(id)%pSC(ia)%onesite_exchange
               do j = 1, IntegrValue(id)%pSC(ia)%NumJs
                  write(6,'(10x,a,t40,''='',2x,i3,2x,3f10.5)') &
                   'LOCAL Pair exchange param',j,IntegrValue(id)%pSC(ia)%pair_exchab(1:7:3,j)
                  write(6,'(t40,8x,3f10.5)')IntegrValue(id)%pSC(ia)%pair_exchab(2:8:3,j)
                  write(6,'(t40,8x,3f10.5)')IntegrValue(id)%pSC(ia)%pair_exchab(3:9:3,j)
               enddo
            endif
         enddo
      enddo
   endif
   write(6,'(/,80(''-''))')
!
   end subroutine printValenceStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calValenceStates(PartialDensity_Ef)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use ChemElementModule, only : getZval
!
   use AtomModule, only : getLocalEvecOld, getLocalNumSpecies, &
                          getLocalSpeciesContent, getLocalAtomicNumber
!
   use SystemModule, only : getAdditionalElectrons
!
   use SSSolverModule, only : initSSSolver, endSSSolver
   use RelSSSolverModule, only : initRelSSSolver, endRelSSSolver !xianglin
   use RelSSSolverModule, only : SingleDiracScattering, computeRelSingleSiteDOS !xianglin
!
   use MSSolverModule, only : initMSSolver, endMSSolver
   use RelMSSolverModule, only : initRelMSSolver, endRelMSSolver, computeRelMST, getRelMSDOS !xianglin
!
!   use RelScattererModule, only : initRelScatterer, endRelScatterer !xianglin
!
!  use RelGreenFunctionModule, only : initRelGreenFunction, &
!                                     endRelGreenFunction
!
   use LdaCorrectionModule, only : checkLdaCorrection, computeLdaPlusU
!
   use SpinRotationModule, only :  resetSpinRotation
!
   use KreinModule,only : calJost_E0, isLloydOn, setLloydStatus, setLloydQ
!
   use ValenceDensityModule, only : getFermiEnergy
!
   use ScfDataModule, only : ErBottom, ErTop !xianglin
!
   implicit none
!
   integer (kind=IntKind) :: id, ia, kmax, num_species
!
   real (kind=RealKind) :: t1, t2, NLloyd, efermi
   real (kind=RealKind), allocatable :: local_density_matrix(:,:,:)
!
   logical, optional, intent(in) :: PartialDensity_Ef
   logical :: isBxyz(LocalNumAtoms) !xianglin
   complex (kind=CmplxKind) :: energy !xianglin
   real (kind=RealKind) :: etop, ebot !for findResonance
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
#ifdef DEBUG
   write(6,*) "calValenceStates: Begin"
   call FlushFile(6)
#endif
!
   zvaltss = ZERO
   do id = 1,LocalNumAtoms
      do ia = 1, getLocalNumSpecies(id)
         zvaltss = zvaltss + &
                   getLocalSpeciesContent(id,ia)*getZval(getLocalAtomicNumber(id,ia))
      enddo
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(aGID,zvaltss)
!  -------------------------------------------------------------------
   zvaltss = zvaltss + getAdditionalElectrons()
!
   if (node_print_level >= 0) then
      write(6,'(/,a,f10.2,/)')'Number of electrons used for finding Fermi energy = ',zvaltss
   endif
!
   if (.not.Initialized) then 
!     ----------------------------------------------------------------
      call ErrorHandler('calValenceStates',                           &
                        'module needs to be initialized first')
!     ----------------------------------------------------------------
   endif
! 
   do id = 1,LocalNumAtoms
      evec(1:3,id) = getLocalEvecOld(id)
      if (RelativisticFlag .ne. 2) then
!        ----------------------------------------------------------------
         call resetSpinRotation(id,evec(1:3,id))
!        ----------------------------------------------------------------
      endif
   enddo
!
   NumCalls_SS = 0
   Timing_SS = ZERO
!
!  ===================================================================
!  initialize Single Site Scatterer
!  -------------------------------------------------------------------
   if (RelativisticFlag == 2) then
!     ----------------------------------------------------------------
      isBxyz=is_Bxyz !xianglin
!YW   call initRelSSSolver(LocalNumAtoms, lmax_kkr, 2*(lmax_kkr), 2*(lmax_kkr), &
      call initRelSSSolver(LocalNumAtoms, lmax_kkr, lmax_pot, lmax_green, &
                           getLocalNumSpecies, getLocalAtomicNumber,      &
                           isBxyz, print_level, stop_routine, evec)
!YW   call initRelMSSolver(LocalNumAtoms, AtomIndex, lmax_kkr,lmax_kkr,2*(lmax_kkr),&
      call initRelMSSolver(LocalNumAtoms, AtomIndex, lmax_kkr,lmax_kkr,lmax_green, &
                           posi,is_Bxyz,stop_routine, print_level )
   else
!     ----------------------------------------------------------------
      call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, &
                        lmax_kkr, lmax_phi, lmax_pot, lmax_step, lmax_green, &
                        n_spin_pola, n_spin_cant, RelativisticFlag,   &
                        stop_routine, print_level, derivative=rad_derivative)
!     ----------------------------------------------------------------
!     initialize Multiple scattering solver
!     ----------------------------------------------------------------
      call initMSSolver(LocalNumAtoms, AtomIndex,                        &
                        lmax_kkr, lmax_phi, lmax_green, posi,            &
                        n_spin_pola, n_spin_cant, RelativisticFlag,      &
                        stop_routine, print_level, derivative=rad_derivative)
!  -------------------------------------------------------------------
   endif
!
   chempot = getFermiEnergy()
   efermi = chempot
!  ===================================================================
!  find resonance part for negative energy
   if ( (ErBottom < ZERO .and. isPole) .or. isPole_plus) then
      allocate(MatrixPoles(LocalNumAtoms))
      do id = 1,LocalNumAtoms
         MatrixPoles(id)%NumElements = getLocalNumSpecies(id)
      enddo
   endif
!
   if ( ErBottom < ZERO .and. isPole) then !xianglin
      etop = min(ZERO,efermi)
      ebot = ErBottom
      nwPG = (nw_Pole/NumPEsInEGroup)+1
      nw_Pole_new = nwPG*NumPEsInEGroup ! to ensure exact division of NW
      do id = 1,LocalNumAtoms
         num_species = getLocalNumSpecies(id)
         if ( RelativisticFlag == 2 ) then
            allocate( MatrixPoles(id)%Poles(2*(lmax_kkr_max+1)**2,num_species,1) ) !"1" here is n_spin
            allocate( MatrixPoles(id)%NumPoles(num_species,1) )
            MatrixPoles(id)%Poles= ZERO
            call calQuadraticPoles(id,num_species,lmax_kkr(id),lmax_kkr(id), &
                                   1,ebot,etop,                       &
                                   2*(lmax_kkr_max+1)**2,             &
                                   MatrixPoles(id)%NumPoles,MatrixPoles(id)%Poles,isRel=.true.)
         else
            allocate( MatrixPoles(id)%Poles((lmax_kkr_max+1)**2,num_species,n_spin_pola) )
            allocate( MatrixPoles(id)%NumPoles(num_species,n_spin_pola) )
            MatrixPoles(id)%Poles= ZERO
!           find the poles for e<0 in nonrelativistic case. Commented because haven't been fully tested
!           call calQuadraticPoles(id,num_species,lmax_kkr(id),lmax_phi(id),     &
!                                  n_spin_pola,ebot,etop,             &
!                                  (lmax_kkr_max+1)**2,MatrixPoles(id)%NumPoles,MatrixPoles(id)%Poles)
         endif
      enddo
   endif
!
   if ( isPole_plus ) then
      etop=efermi
      ebot=max(ErBottom,1.0d-20)
      nwPG_plus = (nw_Pole_plus/NumPEsInEGroup)+1
      nw_Pole_new_plus = nwPG_plus*NumPEsInEGroup ! to ensure exact division of NW
      do id = 1,LocalNumAtoms
         num_species = getLocalNumSpecies(id)
         if ( RelativisticFlag == 2 ) then
            allocate( MatrixPoles(id)%Poles_plus(2*(lmax_kkr_max+1)**2,num_species,1) ) !"1" here is n_spin
            allocate( MatrixPoles(id)%NumPoles_plus(num_species,1) )
            MatrixPoles(id)%Poles_plus= ZERO
            call calQuadraticPoles_cmplx(id,num_species,lmax_kkr(id),lmax_kkr(id), &
                                         1,ebot,etop,                   &
                                         2*(lmax_kkr_max+1)**2,         &
                                         MatrixPoles(id)%NumPoles_plus,MatrixPoles(id)%Poles_plus,isRel=.true.)
!           print*,"NumPoles_plus= ",MatrixPoles(id)%NumPoles_plus
!           print*,"2*(lmax_kkr_max+1)**2= ",2*(lmax_kkr_max+1)**2
!           print*,"Poles_plus= ",MatrixPoles(id)%Poles_plus,SIZE(MatrixPoles(id)%Poles_plus)
         else
            allocate( MatrixPoles(id)%Poles_plus((lmax_kkr_max+1)**2,num_species,n_spin_pola) )
            allocate( MatrixPoles(id)%NumPoles_plus(num_species,n_spin_pola) )
            MatrixPoles(id)%Poles_plus= ZERO
!           find the poles for e<0 in nonrelativistic case. Commented because haven't been fully tested
!           call calQuadraticPoles_cmplx(id,num_species,lmax_kkr(id),lmax_phi(id), &
!                                        n_spin_pola,ebot,etop,         &
!                                        (lmax_kkr_max+1)**2,MatrixPoles(id)%NumPoles_plus,MatrixPoles(id)%Poles_plus)
         endif
      enddo
   endif
!  ===================================================================

   if ( isLloydOn() ) then
!     ----------------------------------------------------------------
      call calJost_E0()
!     ----------------------------------------------------------------
      Lloyd_lastE(1:2)=CZERO
!     ----------------------------------------------------------------
      call setLloydStatus(.true.,.false.)
      call setLloydQ(CZERO)
!     ----------------------------------------------------------------
      if ( isIterateEfOn ) then  ! I modifed the original code and deleted isLloydEf -YW
!        -------------------------------------------------------------
         call findLloydEf(efermi)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         NLloyd = getNLloyd(efermi)
!        -------------------------------------------------------------
      endif
   endif
!
!  ===================================================================
!  calulate the integrated DOS and the Fermi energy
!  -------------------------------------------------------------------
   if (RelativisticFlag == 2) then !xianglin
      call calRelIntegratedDOS(efermi)
   else
     call calIntegratedDOS(efermi)
   endif
!
!  -------------------------------------------------------------------
!  call averageElectroStruct(IntegrValue)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Update the fermi energy: chempot
!  ===================================================================
   chempot = efermi
!
!  ===================================================================
!  call calDensity to compute the electron and moment densities
!  -------------------------------------------------------------------
   call calDensity(efermi)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  If LDA+U is needed, call computeLdaPlusU to compute the potential
!  and energy corrections. The following code needs some work
!  ===================================================================
   if (isDensityMatrixNeeded) then
      kmax = 0
      do id = 1,LocalNumAtoms
         kmax=max(kmax,IntegrValue(id)%kmax)
      enddo
      allocate( local_density_matrix(1:kmax,1:kmax,1:2) )
!
      do id = 1,LocalNumAtoms
         kmax = IntegrValue(id)%kmax
         do ia = 1, IntegrValue(id)%NumSpecies
            if (n_spin_cant == 2) then
               local_density_matrix(1:kmax,1:kmax,1) =                  &
                            real(IntegrValue(id)%density_matrix(1:kmax,1:kmax,1,ia),kind=RealKind)
               local_density_matrix(1:kmax,1:kmax,2) =                  &
                 real(IntegrValue(id)%density_matrix(1:kmax,1:kmax,2,ia)*evec(1,id) &
                    + IntegrValue(id)%density_matrix(1:kmax,1:kmax,3,ia)*evec(2,id) &
                    + IntegrValue(id)%density_matrix(1:kmax,1:kmax,4,ia)*evec(3,id),kind=RealKind)
            else
               local_density_matrix(1:kmax,1:kmax,1:2) =                  &
                            real(IntegrValue(id)%density_matrix(1:kmax,1:kmax,1:2,ia),kind=RealKind)
            endif
!
            kmax = size(local_density_matrix,1)
            if (checkLdaCorrection(id,ia)) then  ! temperary dix
!              -------------------------------------------------------
               call computeLdaPlusU(id,ia,kmax,local_density_matrix)
!              -------------------------------------------------------
            endif
         enddo
      enddo
!
      deallocate( local_density_matrix )
   endif
!
!  ===================================================================
   if (node_print_level >= 0) then
      write(6,'(/,a,i5)')'Number of single site scattering solver calls:', &
                         NumCalls_SS
      write(6,'(a,f12.3,a,/)')'Total of single site scattering solver timing:', &
                         Timing_SS,' (sec)'
   endif
!
   if (RelativisticFlag == 2) then
      call endRelMSSolver() 
      call endRelSSSolver()
   else
      call endMSSolver()
      call endSSSolver()
   endif
!
   if (allocated(MatrixPoles)) then
      if (ErBottom < ZERO .and. isPole) then !xianglin
         do id = 1, LocalNumAtoms
            deallocate(MatrixPoles(id)%NumPoles,MatrixPoles(id)%Poles)
         enddo
      endif
      if (isPole_plus) then
         do id = 1, LocalNumAtoms
            deallocate(MatrixPoles(id)%NumPoles_plus,MatrixPoles(id)%Poles_plus)
         enddo
      endif
      deallocate(MatrixPoles)
   endif
!
   if (nocaseCompare(stop_routine,'calValenceStates')) then
      call StopHandler('calValenceStates','Stop at the end of this routine')
   endif
!
   end subroutine calValenceStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calValenceDOS()
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use AtomModule, only : getLocalEvecOld, getLocalNumSpecies,        &
                          getLocalAtomicNumber
!
   use SSSolverModule, only : initSSSolver, endSSSolver, computeDOS
!
   use MSSolverModule, only : initMSSolver, endMSSolver,              &
                              computeMSGreenFunction
!
!  use RelScattererModule, only : initRelScatterer, endRelScatterer
!
!  use RelGreenFunctionModule, only : initRelGreenFunction, &
!                                      endRelGreenFunction
   use RelSSSolverModule, only : initRelSSSolver, endRelSSSolver !xianglin
   use RelMSSolverModule, only : initRelMSSolver, endRelMSSolver
   use RelMSSolverModule, only : computeRelMST
!
   use SpinRotationModule, only :  resetSpinRotation
!
   use ScfDataModule, only : ErBottom, ErTop, EiBottom, NumSS_IntEs, isLSMS
!
   use DerivativeModule, only : derv5
!
   implicit none
!
   logical :: redundant
   logical :: isBxyz(LocalNumAtoms) !xianglin
!
   integer (kind=IntKind) :: id, is, ne, ie, me, info(4), nvals, js, ns, ia
   integer (kind=IntKind) :: num_species, nsize
!
   real (kind=RealKind) :: e_imag, de, msDOS(2), ssDOS, eb
   real (kind=RealKind), allocatable :: e_real(:)
   real (kind=RealKind), pointer :: p1(:)
!
   complex (kind=CmplxKind) :: energy, ec
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('calValenceDOS',                              &
                        'module needs to be initialized first')
!     ----------------------------------------------------------------
   endif
!
   do id = 1,LocalNumAtoms
      evec(1:3,id) = getLocalEvecOld(id)
      if (RelativisticFlag .ne. 2) then
!        ----------------------------------------------------------------
         call resetSpinRotation(id,evec(1:3,id))
!        ----------------------------------------------------------------
      endif
   enddo
!
   NumCalls_SS = 0
   Timing_SS = ZERO
!
!  ===================================================================
!  initialize Single Site Scatterer
!  -------------------------------------------------------------------
   if (RelativisticFlag == 2) then
!     ----------------------------------------------------------------
      isBxyz=is_Bxyz !xianglin
      call initRelSSSolver(LocalNumAtoms, lmax_kkr, 2*(lmax_kkr), 2*(lmax_kkr), &
                           getLocalNumSpecies, getLocalAtomicNumber,            &
                           isBxyz, print_level, stop_routine, evec)
      call initRelMSSolver(LocalNumAtoms, AtomIndex, lmax_kkr,lmax_kkr,2*(lmax_kkr),&
                           posi,is_Bxyz,stop_routine, print_level )
   else
!     ----------------------------------------------------------------
      call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, &
                        lmax_kkr, lmax_phi, lmax_pot, lmax_step, lmax_green, &
                        n_spin_pola, n_spin_cant, RelativisticFlag,   &
                        stop_routine, print_level)
!     ----------------------------------------------------------------
!     initialize Multiple scattering solver
!     ----------------------------------------------------------------
      call initMSSolver(LocalNumAtoms, AtomIndex,                        &
                        lmax_kkr, lmax_phi, lmax_green, posi,            &
                        n_spin_pola, n_spin_cant, RelativisticFlag,      &
                        stop_routine, print_level)
!  -------------------------------------------------------------------
   endif
!
!  ===================================================================
!  calulate DOS along the real energy axis............................
!  ===================================================================
   eb = ErBottom
!
!  ne = int((chempot+0.2d0-eb)/chempot*NumSS_IntEs)+1
   ne = NumSS_IntEs
!
   allocate(e_real(ne))
!
   if (allocated(SS_dosdata)) then
      do id = 1, LocalNumAtoms
         nullify(SS_dosdata(id)%rarray3)
      enddo
      deallocate(SS_dosdata,wk_SS_dosdata)
   endif
   if (allocated(dosdata)) then
      do id = 1, LocalNumAtoms
         nullify(dosdata(id)%rarray3)
      enddo
      deallocate(dosdata,wk_dosdata)
   endif
!
   if (n_spin_cant == 2) then
      nvals = 7
   else
      nvals = n_spin_pola+1
   endif
!
   nsize = 0
   do id = 1, LocalNumAtoms
      num_species = getLocalNumSpecies(id)
      nsize = nsize+nvals*ne*num_species
   enddo
   allocate(wk_SS_dosdata(nsize))
   allocate(wk_dosdata(nsize))
   allocate(SS_dosdata(LocalNumAtoms))
   allocate(dosdata(LocalNumAtoms))
   nsize = 0
   do id = 1, LocalNumAtoms
      num_species = getLocalNumSpecies(id)
      SS_dosdata(id)%n1=nvals
      SS_dosdata(id)%n2=ne
      SS_dosdata(id)%n3=num_species
      dosdata(id)%n1=nvals
      dosdata(id)%n2=ne
      dosdata(id)%n3=num_species
      p1 => wk_SS_dosdata(nsize+1:nsize+nvals*ne*num_species)
      SS_dosdata(id)%rarray3=>aliasArray3_r(p1,nvals,ne,num_species)
      p1 => wk_dosdata(nsize+1:nsize+nvals*ne*num_species)
      dosdata(id)%rarray3=>aliasArray3_r(p1,nvals,ne,num_species)
      nsize = nsize+nvals*ne*num_species 
   enddo
   wk_SS_dosdata = ZERO
   wk_dosdata = ZERO
!
   me = ne-mod(ne,NumPEsInEGroup)
!  de = (chempot+0.20d0-eb)/real(ne-1,RealKind)
   if (abs(ErTop) < TEN2m8) then
      de = (chempot-eb)/real(ne-1,RealKind)
   else
      de = (ErTop - eb)/real(ne-1,RealKind)
   endif
   do ie = 1, ne
      e_real(ie) = eb + (ie-1)*de
   enddo
   do is = 1, n_spin_pola/n_spin_cant
      ie = 0
      do while (ie < ne)
         ie = ie + 1
         if (abs(e_real(ie)) < 0.001d0) then  ! In case of the energy is too close to the origin
            e_imag = max(0.01d0,EiBottom) ! we take a larger imaginary part
         else if (e_real(ie) < 0.1d0) then     ! In case of existing shallow bound states or low lying resonance states
            e_imag = max(0.005d0,EiBottom) ! we take a larger imaginary part
         else
            e_imag = max(0.001d0,EiBottom) ! in case EiBottom is zero or too small
         endif
         if (ie > me) then
            redundant = .true.
         else
            redundant = .false.
         endif
         if (RelativisticFlag == 2) then !relativistic or not, note n_spin_pola/n_spin_cant = 1 in our formalism
            if (redundant .or. mod(ie,NumPEsInEGroup)==MyPEinEGroup) then
               if (isSingleSiteCluster) then
                  energy = cmplx(e_real(ie),e_imag)
                  do js = 1, 1
                     do id = 1, LocalNumAtoms
                        info(1) = max(is,js); info(2) = id; info(3) = -1; info(4) = -1
!                       -------------------------------------------------
                        ssDOS = returnRelSingleSiteDOSinCP(info,energy,wk_dos)
                        call calElectroStruct(info,4,wk_dos,LastValue(id),ss_int=.true.)
!                       -------------------------------------------------
                     enddo
                  enddo
                  if (node_print_level >= 0) then
                     do js = 1, n_spin_cant*n_spin_cant
                        do id = 1, LocalNumAtoms
                           do ia = 1, LastValue(id)%NumSpecies
                              write(6,'(/,''returnRelSingleSiteDOS:   energy ='',2f18.12,'', id ='',i4)')energy,id,ia
                              write(6,'(''                       Single Site DOS_MT   ='',f18.12)')real(LastValue(id)%dos_mt(js,ia))
                              write(6,'(''                       Single Site DOS_VP   ='',f18.12)')real(LastValue(id)%dos(js,ia))
                           enddo
                        enddo
                     enddo
                  endif
               else !   single-site cluster or not
                  energy = cmplx(e_real(ie),e_imag)
!                 ----------------------------------------------------
                  call computeRelMST(adjustEnergy(is,energy))
!                 ----------------------------------------------------
                  do id =  1,LocalNumAtoms
                     info(1) = is; info(2) = id; info(3) = -1; info(4) = 0
!                    ----------------------------------------------------
                     msDOS = returnRelMultipleSiteDOS(info,energy,wk_dos)
                     call calElectroStruct(info,4,wk_dos,LastValue(id))
!                    ----------------------------------------------------
                     if (e_real(ie) > ZERO) then
                        do js = 1, 1
                           info(1) = max(is,js); info(2) = id; info(3) = -1
                           info (4) = -1 ! do not set internal print to 1, otherwise, it
                                         ! could possiblly hang the MPI process
!                          ----------------------------------------------
                           ssDOS = returnRelSingleSiteDOS(info,e_real(ie),wk_dos)
!                          ----------------------------------------------
                           call calElectroStruct(info,4,wk_dos,ssLastValue(id),ss_int=.true.)
!                          ----------------------------------------------
                        enddo
                        if (node_print_level >= 0) then
                           do ia = 1, LastValue(id)%NumSpecies
                              do js = 1, n_spin_cant*n_spin_cant
                                 write(6,'(''                       Single Site DOS_MT   ='',f18.12)') &
                                       real(ssLastValue(id)%dos_mt(js,ia))
                                 write(6,'(''                       Single Site DOS_VP   ='',f18.12)') &
                                       real(ssLastValue(id)%dos(js,ia))
                              enddo
                           enddo
                        endif             
!                       ----------------------------------------------
                        call addElectroStruct(ONE,ssLastValue(id),LastValue(id))
!                       ----------------------------------------------
                        do ia = 1, LastValue(id)%NumSpecies
                           do js = 1, 1 !n_spin_cant*n_spin_cant
                              ns = max(js,is) !ns = 1:4
                              if (real(LastValue(id)%dos_mt(ns,ia)) < ZERO) then
                                 if (node_print_level >= 0) then
                                    call WarningHandler('calValenceDOS',       &
                                              'dos_mt < 0 due to the energy in MS term having imaginary part', &
                                              real(LastValue(id)%dos_mt(ns,ia)))
                                    call WarningHandler('calValenceDOS','dos_mt is set to 0.0')
                                 endif
                                 LastValue(id)%dos_mt(ns,ia) = CZERO
                              endif
                              if (real(LastValue(id)%dos(ns,ia)) < ZERO) then
                                 if (node_print_level >= 0) then
                                    call WarningHandler('calValenceDOS',       &
                                              'dos < 0 due to the energy in MS term having imaginary part', &
                                              real(LastValue(id)%dos(ns,ia)))
                                    call WarningHandler('calValenceDOS','dos is set to 0.0')
                                 endif
                                 LastValue(id)%dos(ns,ia) = CZERO
                              endif
                           enddo
                        enddo
                     endif
                  enddo
               endif
!              ----------------------------------------------------------
               call calSS_dosdata(is,ie,e_real(ie),ssLastValue,redundant)
               call calDOSDATA(is,ie,e_real(ie),LastValue,redundant)
!              ----------------------------------------------------------
            endif
         else ! relativistic or not
            if (redundant .or. mod(ie,NumPEsInEGroup)==MyPEinEGroup) then
               if (isSingleSiteCluster) then
                  energy = cmplx(e_real(ie),e_imag)
                  do js = 1, n_spin_cant
                     do id = 1, LocalNumAtoms
                        info(1) = max(is,js); info(2) = id; info(3) = -1; info(4) = -1
!                       -------------------------------------------------
                        ssDOS = returnSingleSiteDOSinCP(info,energy,wk_dos)
                        call calElectroStruct(info,1,wk_dos,LastValue(id),ss_int=.true.)
!                       -------------------------------------------------
                     enddo
                  enddo
                  if (node_print_level >= 0) then
                     do js = 1, n_spin_cant*n_spin_cant
                        do id = 1, LocalNumAtoms
                           do ia = 1, LastValue(id)%NumSpecies
                              write(6,'(/,''returnSingleSiteDOS:   energy ='',2f18.12,'', id ='',i4)')energy,id
                              write(6,'(''                       Single Site DOS_MT   ='',f18.12)') &
                                    real(LastValue(id)%dos_mt(js,ia))
                              write(6,'(''                       Single Site DOS_VP   ='',f18.12)') &
                                    real(LastValue(id)%dos(js,ia))
                           enddo
                        enddo
                     enddo
                  endif
               else
                  energy = cmplx(e_real(ie),e_imag)
                  if (e_real(ie) <= ZERO) then
!                    ----------------------------------------------------
                     call computeMSGreenFunction(is, adjustEnergy(is,energy), add_Gs=.true., isSphSolver=.true.)
!                    ----------------------------------------------------
                  else
!                    ----------------------------------------------------
                     call computeMSGreenFunction(is, adjustEnergy(is,energy))
!                    ----------------------------------------------------
                  endif
                  do id =  1,LocalNumAtoms
                     info(1) = is; info(2) = id; info(3) = -1; info(4) = 0
!                    ----------------------------------------------------
                     msDOS = returnMultipleSiteDOS(info,energy,wk_dos)
                     call calElectroStruct(info,n_spin_cant,wk_dos,LastValue(id))
!                    ----------------------------------------------------
                     if (e_real(ie) > ZERO) then
                        do js = 1, n_spin_cant
                           info(1) = max(is,js); info(2) = id; info(3) = -1 
                           info(4) = -1 ! do not set internal print to 1, otherwise, it
                                        ! could possiblly hang the MPI process
!                          ----------------------------------------------
                           ssDOS = returnSingleSiteDOS(info,e_real(ie),wk_dos)
!                          ----------------------------------------------
                           call calElectroStruct(info,1,wk_dos,ssLastValue(id),ss_int=.true.)
!                          ----------------------------------------------
                        enddo
                        if (node_print_level >= 0) then
                           do ia = 1, LastValue(id)%NumSpecies
                              do js = 1, n_spin_cant*n_spin_cant
                                 write(6,'(''                       Single Site DOS_MT   ='',f18.12)') &
                                       real(ssLastValue(id)%dos_mt(js,ia))
                                 write(6,'(''                       Single Site DOS_VP   ='',f18.12)') &
                                       real(ssLastValue(id)%dos(js,ia))
                              enddo
                           enddo
                        endif
                        if (n_spin_cant == 2) then
!                          ----------------------------------------------
                           call addElectroStruct(ONE,ssLastValue(id),LastValue(id))
!                          ----------------------------------------------
                        else
!                          ----------------------------------------------
                           call addElectroStruct(ONE,ssLastValue(id),LastValue(id),is)
!                          ----------------------------------------------
                        endif
                        do ia = 1, LastValue(id)%NumSpecies
                           do js = 1, 1   !n_spin_cant*n_spin_cant
                              ns = max(js,is)
                              if (real(LastValue(id)%dos_mt(ns,ia)) < ZERO) then
                                 if (node_print_level >= 0) then
                                    call WarningHandler('calValenceDOS',       &
                                             'dos_mt < 0 due to the energy in MS term having imaginary part', &
                                             real(LastValue(id)%dos_mt(ns,ia)))
                                    call WarningHandler('calValenceDOS','dos_mt is set to 0.0')
                                 endif
                                 LastValue(id)%dos_mt(ns,ia) = CZERO
                              endif
                              if (real(LastValue(id)%dos(ns,ia)) < ZERO) then
                                 if (node_print_level >= 0) then
                                    call WarningHandler('calValenceDOS',       &
                                            'dos < 0 due to the energy in MS term having imaginary part', &
                                            real(LastValue(id)%dos(ns,ia)))
                                    call WarningHandler('calValenceDOS','dos is set to 0.0')
                                 endif
                                 LastValue(id)%dos(ns,ia) = CZERO
                              endif
                           enddo
                        enddo
                     endif
                  enddo
               endif
!              ----------------------------------------------------------
               call calSS_dosdata(is,ie,e_real(ie),ssLastValue,redundant)
               call calDOSDATA(is,ie,e_real(ie),LastValue,redundant)
!              ----------------------------------------------------------
            endif
         endif  ! relativistic or not
      enddo
   enddo
!
   if ( NumPEsInEGroup > 1 ) then
!     ----------------------------------------------------------------
      call GlobalSumInGroup(eGID, wk_SS_dosdata, nsize)
      call GlobalSumInGroup(eGID, wk_dosdata, nsize)
!     ----------------------------------------------------------------
   endif
!
   if (node_print_level >= 0) then
      write(6,'(/,a,i5)')'Number of single site scattering solver calls:', &
                         NumCalls_SS
      write(6,'(a,f12.3,a,/)')'Total of single site scattering solver timing:', &
                         Timing_SS,' (sec)'
   endif
!
   deallocate(e_real)
!
   if (RelativisticFlag == 2) then
!     ----------------------------------------------------------------
      call endRelSSSolver()
      call endRelMSSolver()
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call endMSSolver()
      call endSSSolver()
!     ----------------------------------------------------------------
   endif
!
   if (nocaseCompare(stop_routine,'calValenceDOS')) then
      call StopHandler('calValenceDOS','Stop at the end of this routine')
   endif
!
   end subroutine calValenceDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntegratedDOS(efermi)
!  ===================================================================
   use PublicParamDefinitionsModule, only : ButterFly
!
   use PhysParamModule, only : Boltzmann
!
   use MPPModule, only : endMPP
!
   use GroupCommModule, only : GlobalMaxInGroup
!
   use ScfDataModule, only : ErBottom, ErTop, EiBottom, EiTop, Temperature
   use ScfDataModule, only : NumEs, ContourType, eGridType, isReadEmesh, getEmeshFileName
   use ScfDataModule, only : isKKR, isKKRCPA, isSSIrregularSolOn
   use ScfDataModule, only : NumSS_IntEs
!
   use AtomModule, only : getLocalSpeciesContent
!
   use PotentialTypeModule, only : isFullPotential
!
   use KreinModule, only : isLloydOn
!
   use ContourModule, only : initContour, endContour, setupContour
   use ContourModule, only : printContour
!
   use MSSolverModule, only : computeMSGreenFunction
!
   use ValenceDensityModule, only : printValenceDensity
!
   implicit none
!
   logical :: useIrregularSolution = .true.
   logical :: isContourInitialized = .false.
!
   integer (kind=IntKind) :: id, BadFermiEnergy, is, info(4), ns, ia
   integer (kind=IntKind), parameter :: MaxIterations = 20
!
   real (kind=RealKind), intent(inout) :: efermi
   real (kind=RealKind) :: efermi_old, dosmt, dosws, kBT, rfac
   real (kind=RealKind) :: Lloyd_factor(n_spin_pola), ssDOS, msDOS(2)
   real (kind=RealKind) :: int_dos(n_spin_cant*n_spin_pola,LocalNumAtoms)
   real (kind=RealKind) :: last_dos(n_spin_cant*n_spin_pola,LocalNumAtoms)
!
   complex (kind=CmplxKind) :: eLast
!
!  ===================================================================
!  Determine if using Z*Tau*Z - Z*J formula for the Green function
!  ===================================================================
   useIrregularSolution = isSSIrregularSolOn() .and. .not.(isFullPotential()) 
!
   kBT = Temperature*Boltzmann
!
!  ===================================================================
!  Initialize the enegy contour module
!  ===================================================================
   isContourInitialized = .false.
   if (isReadEmesh()) then
!     ----------------------------------------------------------------
      call initContour(getEmeshFileName(), stop_routine, max_print_level)
!     ----------------------------------------------------------------
      isContourInitialized = .true.
      if (node_print_level >= 0) then
!        -------------------------------------------------------------
         call printContour()
!        -------------------------------------------------------------
      endif
   else if (.not.isSingleSiteCluster .or. useIrregularSolution) then
!     ----------------------------------------------------------------
      call initContour( ContourType, eGridType, NumEs, Temperature,   &
                        stop_routine, maxval(print_level(1:LocalNumAtoms)), .true. )
!     ----------------------------------------------------------------
      isContourInitialized = .true.
!
!     ================================================================
!     set up the enegy contour that ends at (efermi,0) for Gaussian grids
!     on a semi-circle contour or ends at (efermi,Eibottom) for uniform
!     grids on arectangular contour.
!     New feature: No extra last energy point is added!
!     ----------------------------------------------------------------
      call setupContour( ErBottom, efermi, EiBottom, EiTop )
!     ----------------------------------------------------------------
      if (node_print_level >= 0) then
!        -------------------------------------------------------------
         call printContour()
!        -------------------------------------------------------------
      endif
   endif
!
   do id =  1, LocalNumAtoms
!     ----------------------------------------------------------------
      call zeroElectroStruct(IntegrValue(id))
!     ----------------------------------------------------------------
   enddo
!
!  ===================================================================
!  compute the DOS arising from the multiple scattering term along
!  the energy contour set up by setupContour
!  ===================================================================
 ! if (isSingleSiteCluster .and.  ContourType == ButterFly) then
!!    ----------------------------------------------------------------
 !    call calSingleScatteringIDOS(LowerContour=.true.)
!!    ----------------------------------------------------------------
 ! else
   if (.not.isSingleSiteCluster .or. useIrregularSolution) then
!     ----------------------------------------------------------------
      call calMultipleScatteringIDOS(useIrregularSolution)
!     ----------------------------------------------------------------
   endif
!
   if ( node_print_level >= 0) then
      write(6,'(/,a)')'The IDOS of the MS term'
      do id =  1,LocalNumAtoms
         do ia = 1, IntegrValue(id)%NumSpecies
            do is = 1, n_spin_pola*n_spin_cant
               write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                       ', MS    IDOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)),&
                       ', MS    IDOS_ws = ',real(IntegrValue(id)%dos(is,ia))
            enddo
         enddo
      enddo
   !  else
   !     do id =  1,LocalNumAtoms
   !        do ia = 1, IntegrValue(id)%NumSpecies
   !           do is = 1, n_spin_pola
   !              dosmt=HALF*( real(IntegrValue(id)%dos_mt(1,ia)) +                      &
   !                           (3-2*is)*(real(IntegrValue(id)%dos_mt(2,ia))*evec(1,id) + &
   !                                     real(IntegrValue(id)%dos_mt(3,ia))*evec(2,id) + &
   !                                     real(IntegrValue(id)%dos_mt(4,ia))*evec(3,id)) )
   !              dosws=HALF*( real(IntegrValue(id)%dos(1)) +                      &
   !                           (3-2*is)*(real(IntegrValue(id)%dos(2,ia))*evec(1,id) + &
   !                                     real(IntegrValue(id)%dos(3,ia))*evec(2,id) + &
   !                                     real(IntegrValue(id)%dos(4,ia))*evec(3,id)) )
   !              write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is =',is, &
   !                      ', MS    IDOS_mt = ',dosmt,', MS    IDOS_ws = ',dosws
   !           enddo
   !        enddo
   !     enddo
   !  endif
   endif
!
   if (.not.useIrregularSolution) then
!     print *,'Not use irregular solution'
!!    if (ContourType == ButterFly) then
         do id =  1,LocalNumAtoms
            call zeroElectroStruct(ssIntegrValue(id))
         enddo
!        =============================================================
!        compute the DOS arising from the single site scattering term along
!        the real energy axis from Eb < 0 upto E = 0, for taking care
!        if shallow states in single scattering potential
!        =============================================================
!        -------------------------------------------------------------
         call calSingleScatteringIDOS(LowerContour=.true.)
!        -------------------------------------------------------------
         do id =  1,LocalNumAtoms
           if ( node_print_level >= 0) then
               do ia = 1, IntegrValue(id)%NumSpecies
                  do is = 1, n_spin_pola*n_spin_cant
                    write(6,'(3(a,i2),2(a,d15.8))')'Before id = ',id,  &
                          ', ia = ',ia,', is = ',is,                   &
                          ', MS    IDOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)),&
                          ', MS    IDOS_ws = ',real(IntegrValue(id)%dos(is,ia))
                  enddo
               enddo
            endif
!           ----------------------------------------------------------
            call addElectroStruct(CONE,ssIntegrValue(id),IntegrValue(id))
!           ----------------------------------------------------------
            if ( node_print_level >= 0) then
               do ia = 1, IntegrValue(id)%NumSpecies
                  do is = 1, n_spin_pola*n_spin_cant
                    write(6,'(3(a,i2),2(a,d15.8))')'After  id = ',id, &
                          ', ia = ',ia,', is = ',is,                  &
                          ', MS    IDOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)),&
                          ', MS    IDOS_ws = ',real(IntegrValue(id)%dos(is,ia))
                  enddo
               enddo
            endif
         enddo
!!    endif
!
!     ================================================================
!     compute the DOS arising from the single site scattering term along
!     the real energy axis from E = 0 upto E = efermi
!     ================================================================
      do id =  1,LocalNumAtoms
!        -------------------------------------------------------------
         call zeroElectroStruct(ssIntegrValue(id))
!        -------------------------------------------------------------
      enddo
!!    if (NumSS_IntEs > 1 .and. efermi > ZERO) then
!        =============================================================
!        For finite temperature, extend the integration to efermi+6*log(10)*kB*T
!        where 6*log(10)*kB*T is a result of the fact that at energy = efermi+13.82*kB*T
!        FD(energy,T) = 1.0d-6.
!        -------------------------------------------------------------
         call calSingleScatteringIDOS(Ebegin=0.00001d0,               &
                                      Eend=efermi+8.0d0*log(10.0d0)*kBT)
!        -------------------------------------------------------------
!!    else
!        -------------------------------------------------------------
!!       call calSingleScatteringIDOS(UpperContour=.true.)
!        -------------------------------------------------------------
!!    endif
!
!     ================================================================
!     Note: IntegrValue(id)%dos = is the MST part of integrated DOS up to (chempot,0) on a 
!                                 Gaussian grid contour for atom id
!         ssIntegrValue(id)%dos = is the single site part of integrated DOS up to efermi
!                                 on the real energy axis for atom id
!
!     Now, add the single site IDOS to the multiple scattering IDOS term
!     ================================================================
      if ( node_print_level >= 0) then
         write(6,'(/,a)')'Before adding the IDOS of the SS term to the IDOS of the MS term'
         do id =  1,LocalNumAtoms
            do ia = 1, ssIntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is,  &
                          ', SS    IDOS_mt = ',real(ssIntegrValue(id)%dos_mt(is,ia)), &
                          ', SS    IDOS_ws = ',real(ssIntegrValue(id)%dos(is,ia))
               enddo
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is,  &
                          ', MS    IDOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)),&
                          ', MS    IDOS_ws = ',real(IntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
!
      do id =  1,LocalNumAtoms
!        -------------------------------------------------------------
         call addElectroStruct(CONE,ssIntegrValue(id),IntegrValue(id))
!        -------------------------------------------------------------
      enddo
!
      if ( node_print_level >= 0) then
         write(6,'(/,a)')'After adding the IDOS of the SS term to the IDOS of the MS term'
         do id =  1,LocalNumAtoms
            do ia = 1, IntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', MS+SS IDOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)),&
                          ', MS+SS IDOS_ws = ',real(IntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
   else
!     print *,'Use irregular solution'
      if ( node_print_level >= 0) then
         write(6,'(/,a)')'The IDOS associated with the atom'
         do id =  1,LocalNumAtoms
            do ia = 1, IntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', IDOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)),     &
                          ', IDOS_ws = ',real(IntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
   endif
!
!  ===================================================================
!  Transform the data structure into moment density representation
!  so that IntegrValue%dos(1)   = charge density 
!          IntegrValue%dos(2:4) = moment density
!  ===================================================================
   if (n_spin_cant == 2) then
      do id =  1,LocalNumAtoms
!        -------------------------------------------------------------
         call transformElectroStruct(id,IntegrValue(id))
!        -------------------------------------------------------------
         if ( node_print_level >= 0) then
            write(6,'(/,a)')'After transform, the IDOS associated with the atom'
            do ia = 1, IntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola
                  dosmt=HALF*( real(IntegrValue(id)%dos_mt(1,ia)) +                      &
                               (3-2*is)*(real(IntegrValue(id)%dos_mt(2,ia))*evec(1,id) + &
                                         real(IntegrValue(id)%dos_mt(3,ia))*evec(2,id) + &
                                         real(IntegrValue(id)%dos_mt(4,ia))*evec(3,id)) )
                  dosws=HALF*( real(IntegrValue(id)%dos(1,ia)) +                      &
                               (3-2*is)*(real(IntegrValue(id)%dos(2,ia))*evec(1,id) + &
                                         real(IntegrValue(id)%dos(3,ia))*evec(2,id) + &
                                         real(IntegrValue(id)%dos(4,ia))*evec(3,id)) )
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', MS    IDOS_mt = ',dosmt,', MS    IDOS_ws = ',dosws
               enddo
            enddo
         endif
      enddo
   endif
!
!  ===================================================================
!  This is an important stage, may need to stop for debugging purposes
!  ===================================================================
   if (stop_routine == 'BeforeMufind') then
      call syncAllPEs()
      call endMPP()
      call StopHandler('BeforeMufind')
   endif
!  ===================================================================
!
!  ===================================================================
!  Iterate the ending point of the energy contour to find the Fermi energy
!  ===================================================================
   BadFermiEnergy = 1
   LOOP_LastE: do while (BadFermiEnergy > 0 .and. BadFermiEnergy <= MaxIterations)
!     ===============================================================
!     Solve the multiple scattering problem for e = eLast, which is
!     set to be (efermi,0.001) in KKR case.
!     ===============================================================
      if (isKKR() .or. isKKRCPA()) then
         eLast = cmplx(efermi,0.001d0,kind=CmplxKind)
      else
         eLast = cmplx(efermi,0.000d0,kind=CmplxKind)
      endif
!
!     ===============================================================
!     Calculate the DOS of the multiple scattering term for e = eLast,
!     store the data in LastValue structure, and add the single site 
!     results (at efermi, not eLast) to it.
!          LastValue(id)%dos = is the MST part of the DOS at (chempot,eib) for atom id
!     ===============================================================
      if ( node_print_level >= 0) then
         write(6,'(/,a,2d15.8)')'M.S. Term DOS at the last energy: ',eLast
      endif
      do id =  1,LocalNumAtoms
!        -------------------------------------------------------------
         call zeroElectroStruct(LastValue(id))
!        -------------------------------------------------------------
      enddo
      do is = 1, n_spin_pola/n_spin_cant
         if (isFullPotential() .or. .not.useIrregularSolution) then
!           =========================================================
!           In this case, the multiple scattering module returns DOS of
!           Z*(tau-t)*Z, instead of Z*tau*Z-Z*J. One needs to add the single
!           site DOS to the result.
!           ----------------------------------------------------------
            call computeMSGreenFunction(is,adjustEnergy(is,eLast))
!           ----------------------------------------------------------
         else
!           =========================================================
!           In this case, the multiple scattering module returns DOS of
!           Z*tau*Z-Z*J.
!           ----------------------------------------------------------
            call computeMSGreenFunction(is,adjustEnergy(is,eLast),    &
                                        add_Gs=.true.,isSphSolver=.true.)
!           ----------------------------------------------------------
         endif
         do id =  1,LocalNumAtoms
            info(1) = is; info(2) = id; info(3) = -1; info(4) = 1
!           ----------------------------------------------------------
            msDOS = returnMultipleSiteDOS(info,eLast,wk_dos)
            call calElectroStruct(info,n_spin_cant,wk_dos,LastValue(id))
!           ----------------------------------------------------------
            if ( node_print_level >= 0) then
               do ia = 1, LastValue(id)%NumSpecies
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', M.Site DOS_mt = ',real(LastValue(id)%dos_mt(is,ia)), &
                          ', M.Site DOS_ws = ',real(LastValue(id)%dos(is,ia))
               enddo
            endif
         enddo
      enddo
      if (isFullPotential() .or. .not.useIrregularSolution) then
!        ============================================================
!        In this case, the multiple scattering module returns DOS of
!        Z*(tau-t)*Z, instead of Z*tau*Z-Z*J. One needs to add the single
!        site DOS to the result.
!        Solve the single scattering problem for e = efermi.
!        ssLastValue(id)%dos = is the single site part of the DOS at efermi on the
!                              real energy axis for atom id
!        ============================================================
         if ( node_print_level >= 0) then
            write(6,'(/,a,d15.8)')'S.S. Term DOS at the last energy: ',efermi
         endif
         do id =  1,LocalNumAtoms
!           ----------------------------------------------------------
            call zeroElectroStruct(ssLastValue(id))
!           ----------------------------------------------------------
            do is = 1, n_spin_pola
               info(1) = is; info(2) = id; info(3) = -1; info(4) = 0
!              -------------------------------------------------------
               ssDOS = returnSingleSiteDOS(info,efermi,wk_dos)
               call calElectroStruct(info,1,wk_dos,ssLastValue(id),ss_int=.true.)
!              -------------------------------------------------------
            enddo
            if ( node_print_level >= 0) then
               do ia = 1, ssLastValue(id)%NumSpecies
                  do is = 1, n_spin_pola
                     write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                             ', S.Site DOS_mt = ',real(ssLastValue(id)%dos_mt(is,ia)), &
                             ', S.Site DOS_ws = ',real(ssLastValue(id)%dos(is,ia))
                  enddo
               enddo
            endif
!
!           =========================================================
!           Now, add the single site DOS to the multiple scattering DOS term at last
!           energy point
!           ----------------------------------------------------------
            call addElectroStruct(ONE,ssLastValue(id),LastValue(id))
!           ----------------------------------------------------------
         enddo
      endif
!
      if ( node_print_level >= 0) then
         write(6,'(/,a,2d15.8)')'M.S.+S.S. DOS at the last energy: ',eLast
         do id =  1,LocalNumAtoms
            do ia = 1, LastValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', MS+SS  DOS_mt = ',real(LastValue(id)%dos_mt(is,ia)), &
                          ', MS+SS  DOS_ws = ',real(LastValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
!
      if (n_spin_cant == 2) then
         do id =  1,LocalNumAtoms
!           ----------------------------------------------------------
            call transformElectroStruct(id,LastValue(id))
!           ----------------------------------------------------------
         enddo
      endif
!
      last_dos = ZERO; int_dos = ZERO
      do id = 1, LocalNumAtoms
         do ia = 1, LastValue(id)%NumSpecies
            rfac = getLocalSpeciesContent(id,ia)
            do is = 1, n_spin_cant*n_spin_pola
               last_dos(is,id) = last_dos(is,id) + rfac*real(LastValue(id)%dos(is,ia),kind=RealKind)
               int_dos(is,id) =  int_dos(is,id) + rfac*real(IntegrValue(id)%dos(is,ia),kind=RealKind)
            enddo
         enddo
      enddo
!
      efermi_old = efermi
!     ================================================================
!     Compute the Fermi energy and efermi will be renewed.
!     ----------------------------------------------------------------
      call mufind(efermi,efermi_old,int_dos,last_dos,BadFermiEnergy,Lloyd_factor)
!     ----------------------------------------------------------------
      do id = 1, LocalNumAtoms
!        -------------------------------------------------------------
         call addElectroStruct(efermi-efermi_old,LastValue(id),IntegrValue(id))
!        -------------------------------------------------------------
      enddo
      if ( node_print_level >= 0) then
         write(6,'(/,a)')'After mufind, the integrated DOS are:'
         do id =  1,LocalNumAtoms
            do ia = 1, IntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', MS+SS  DOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)), &
                          ', MS+SS  DOS_ws = ',real(IntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
!
      if (isEfPinning) then
         efermi = efermi_old
      endif
!
      if ( .not.isIterateEfOn .or. isLloydOn() ) then
         exit LOOP_LastE
      else
!        -------------------------------------------------------------
         call GlobalMaxInGroup(aGID,BadFermiEnergy)
!        -------------------------------------------------------------
      endif
   enddo Loop_LastE
!
!  -------------------------------------------------------------------
!  call updateValenceDOS(efermi,Lloyd_factor)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  end Energy contour modules
!  ===================================================================
   if (isContourInitialized) then
!     ----------------------------------------------------------------
      call endContour()
!     ----------------------------------------------------------------
   endif
!
   if (stop_routine == 'calIntegratedDOS') then
      call syncAllPEs()
      call endMPP()
      call StopHandler('calIntegratedDOS')
   endif
!
   end subroutine calIntegratedDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calPartialDOS(er)
!  ===================================================================
   use Atom2ProcModule, only : getGlobalIndex
!
   use SSSolverModule, only : initSSSolver, endSSSolver
   use SSSolverModule, only : solveSingleScattering, computePDOS,     &
                              computePhaseShift, getPhaseShift,       &
                              getCellDOS, getMTSphereDOS,             &
                              getCellPDOS, getMTSpherePDOS
!
   use MSSolverModule, only : initMSSolver, endMSSolver
   use MSSolverModule, only : computeMSPDOS, getMSCellDOS, getMSMTSphereDOS, &
                              getMSCellPDOS, getMSMTSpherePDOS
!
   use RelScattererModule, only : initRelScatterer, endRelScatterer
!
   use RelGreenFunctionModule, only : initRelGreenFunction, &
                                      endRelGreenFunction
!
   use RadialGridModule, only : getGrid
!
   use ScfDataModule, only : isKKR, isKKRCPA
!
   use PotentialTypeModule, only : isASAPotential
!
   use AtomModule, only : getLocalEvecNew, getLocalEvecOld, getLocalAtomName
   use AtomModule, only : getLocalNumSpecies, getLocalAtomicNumber
!
   use SpinRotationModule, only : transformDensityMatrix
   use SpinRotationModule, only : resetSpinRotation
!
   implicit none
!
   character (len=10) :: fname
!
   integer (kind=IntKind) :: kmax_phi(LocalNumAtoms)
   integer (kind=IntKind) :: ia, id, is, js, iend, kl, n, ns_sqr, ks, l, m, klc
   integer (kind=IntKind) :: kmax_phi_max, iprint, num_species
!
   real (kind=RealKind), intent(in) :: er
   real (kind=RealKind) :: sfac, evec(3), check_ws, check_mt, t0
   real (kind=RealKind), pointer :: ps(:)
   real (kind=RealKind), pointer :: pdos_ws(:), pdos_mt(:)
   real (kind=RealKind) :: ms_dos_ws(n_spin_cant*n_spin_pola), ms_dos_mt(n_spin_cant*n_spin_pola)
!
   complex (kind=CmplxKind) :: ec
!
   complex (kind=CmplxKind), pointer :: green(:,:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:)
!
   type (GridStruct), pointer :: Grid
!
   type PDOSStruct
      real (kind=RealKind), pointer :: dos_ws(:,:), dos_mt(:,:)
      real (kind=RealKind), pointer :: ss_dos_ws(:,:), ss_dos_mt(:,:)
      real (kind=RealKind), pointer :: phase_shift(:,:,:)
      real (kind=RealKind), pointer :: ss_pdos_ws(:,:,:), ss_pdos_mt(:,:,:)
      real (kind=RealKind), pointer :: ms_pdos_ws(:,:,:), ms_pdos_mt(:,:,:)
      real (kind=RealKind), pointer :: partial_dos_ws(:,:,:), partial_dos_mt(:,:,:)
   end type PDOSStruct
!
   type (PDOSStruct), allocatable :: PartialDOS(:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('calpartialDOS',                              &
                        'module needs to be initialized first')
!     ----------------------------------------------------------------
   endif
!
   do id = 1,LocalNumAtoms
      evec(1:3) = getLocalEvecOld(id)
!     ----------------------------------------------------------------
      call resetSpinRotation(id,evec(1:3))
!     ----------------------------------------------------------------
   enddo
!
   NumCalls_SS = 0
   Timing_SS = ZERO
!
!  ===================================================================
!  initialize Single Site Scatterer
!  -------------------------------------------------------------------
   if (RelativisticFlag == 2) then
!     ----------------------------------------------------------------
      call initRelScatterer(LocalNumAtoms, lmax_kkr, getLocalAtomicNumber,    &
                            print_level, stop_routine)
!     ----------------------------------------------------------------
      call initRelGreenFunction(LocalNumAtoms, AtomIndex, lmax_kkr,    &
                                lmax_phi, lmax_pot, posi, n_spin_pola, &
                                n_spin_cant, stop_routine, print_level)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call initSSSolver(LocalNumAtoms, getLocalNumSpecies, getLocalAtomicNumber, &
                        lmax_kkr, lmax_phi, lmax_pot, lmax_step, lmax_green, &
                        n_spin_pola, n_spin_cant, RelativisticFlag,   &
                        stop_routine, print_level)
!     ----------------------------------------------------------------
   endif
!
   kmax_phi_max = 1
   do id = 1, LocalNumAtoms
      kmax_phi(id) = (lmax_phi(id)+1)**2
      kmax_phi_max = max(kmax_phi(id),kmax_phi_max)
   enddo
!  
   allocate( PartialDOS(LocalNumAtoms) )
   do id = 1, LocalNumAtoms
      num_species = getLocalNumSpecies(id)
      allocate( PartialDOS(id)%dos_ws(n_spin_pola,num_species) )      
      allocate( PartialDOS(id)%dos_mt(n_spin_pola,num_species) )
      allocate( PartialDOS(id)%ss_dos_ws(n_spin_pola,num_species) )      
      allocate( PartialDOS(id)%ss_dos_mt(n_spin_pola,num_species) )
      allocate( PartialDOS(id)%phase_shift(kmax_phi_max,n_spin_pola,num_species) )
      allocate( PartialDOS(id)%ss_pdos_mt(kmax_phi_max,n_spin_pola,num_species) )
      allocate( PartialDOS(id)%ss_pdos_ws(kmax_phi_max,n_spin_pola,num_species) )
      allocate( PartialDOS(id)%ms_pdos_mt(kmax_phi_max,n_spin_cant*n_spin_pola,num_species) )
      allocate( PartialDOS(id)%ms_pdos_ws(kmax_phi_max,n_spin_cant*n_spin_pola,num_species) )
      allocate( PartialDOS(id)%partial_dos_mt(kmax_phi_max,n_spin_pola,num_species) )
      allocate( PartialDOS(id)%partial_dos_ws(kmax_phi_max,n_spin_pola,num_species) )
   enddo
!  
   sfac= TWO/real(n_spin_pola,kind=RealKind)
!  
   if ( (MyPEinEGroup == 0 .and. GlobalNumAtoms < 10) .or. node_print_level >= 0 ) then
      do id = 1, LocalNumAtoms
         if (getGlobalIndex(id) < 10) then
            write(fname,'(a,i1)')trim(getLocalAtomName(id)),getGlobalIndex(id)
         else if (getGlobalIndex(id) < 100) then
            write(fname,'(a,i2)')trim(getLocalAtomName(id)),getGlobalIndex(id)
         else
            write(fname,'(a,i3)')trim(getLocalAtomName(id)),getGlobalIndex(id)
         endif
         open(unit=1001+id,file=trim(fname)//'_AtomicPartialDOS.dat',status='unknown',form='formatted')
         open(unit=2001+id,file=trim(fname)//'_CrystalPartialDOS.dat',status='unknown',form='formatted')
      enddo
   endif
!
   do id = 1, LocalNumAtoms
      Grid => getGrid(id)
      do is = 1, n_spin_pola
         ec = adjustEnergy(is,cmplx(er,0.0d0,kind=CmplxKind))
!
         t0 = getTime()
!        -------------------------------------------------------------
         call solveSingleScattering(is,id,ec,CZERO)
!        -------------------------------------------------------------
         Timing_SS = Timing_SS + (getTime() - t0)
         NumCalls_SS = NumCalls_SS + 1
!
         call computePDOS(spin=is,site=id)
         do ia = 1, getLocalNumSpecies(id)
            PartialDOS(id)%ss_dos_ws(is,ia) = sfac*getCellDOS(spin=is,site=id,atom=ia)
            PartialDOS(id)%ss_dos_mt(is,ia) = sfac*getMTSphereDOS(spin=is,site=id,atom=ia)
            pdos_ws => getCellPDOS(spin=is,site=id,atom=ia)
            pdos_mt => getMTSpherePDOS(spin=is,site=id,atom=ia)
            do kl = 1, kmax_phi(id) 
               PartialDOS(id)%ss_pdos_ws(kl,is,ia) = sfac*pdos_ws(kl)
               PartialDOS(id)%ss_pdos_mt(kl,is,ia) = sfac*pdos_mt(kl)
            enddo
         enddo
!        -------------------------------------------------------------
         call computePhaseShift(spin=is,site=id)
!        -------------------------------------------------------------
         do ia = 1, getLocalNumSpecies(id)
            ps => getPhaseShift(spin=is,site=id,atom=ia)
            do kl = 1, kmax_phi(id)
               PartialDOS(id)%phase_shift(kl,is,ia) = ps(kl)
            enddo
!        
            if ( (MyPEinEGroup == 0 .and. GlobalNumAtoms < 10) .or. node_print_level >= 0 ) then
               write(1001+id,'(3a,i1,a,f13.8)')'Atom: ',getLocalAtomName(id,ia),', Spin: ',is, &
                                               ', Partial Phase Shift and DOS at Energy = ',er
               write(1001+id,'(a,2f13.8)')'Total DOS in Cell and MT Volumes = ',PartialDOS(id)%ss_dos_ws(is,ia), &
                                          PartialDOS(id)%ss_dos_mt(is,ia)
               write(1001+id,'(a)')'=========================================================================='
               write(1001+id,'(a)')'  l     m         Partial_Phase         Partial_DOS         Partial_DOS_MT'
               do kl = 1, kmax_phi(id)
                  write(1001+id,'(i3,3x,i3,3(9x,f12.8))')lofk(kl),mofk(kl),PartialDOS(id)%phase_shift(kl,is,ia), &
                                  PartialDOS(id)%ss_pdos_ws(kl,is,ia),PartialDOS(id)%ss_pdos_mt(kl,is,ia)
               enddo
            endif
            if ( node_print_level >= 0) then
               check_ws = ZERO; check_mt = ZERO
               do kl = kmax_phi(id), 1, -1
                  check_ws = check_ws + PartialDOS(id)%ss_pdos_ws(kl,is,ia)
                  check_mt = check_mt + PartialDOS(id)%ss_pdos_mt(kl,is,ia)
               enddo
               write(6,'(/,3a,i1,a,f11.8,a)')'Atom: ',getLocalAtomName(id,ia),', Spin: ',is, &
                                             ', Total Single Atom DOS at Energy = ',er,':'
               write(6,'(a,2(f13.8,a))')'SS_DOS          = ',check_ws,'(ws)',check_mt,'(mt)'
            endif
         enddo
      enddo
   enddo
!
!  ==================================================================
!  Solve the multiple scattering problem for e = er, which is
!  set to be (er,0.001) in KKR case.
!
!  initialize Multiple scattering solver
!  -------------------------------------------------------------------
   call initMSSolver(LocalNumAtoms, AtomIndex,                        &
                     lmax_kkr, lmax_phi, lmax_green, posi,            &
                     n_spin_pola, n_spin_cant, RelativisticFlag,      &
                     stop_routine, print_level)
!  -------------------------------------------------------------------
   ns_sqr = n_spin_cant*n_spin_cant
   do is = 1, n_spin_pola/n_spin_cant
      if (isKKR() .or. isKKRCPA()) then
         ec = adjustEnergy(is,cmplx(er,0.001d0,kind=CmplxKind))
      else
         ec = adjustEnergy(is,cmplx(er,0.000d0,kind=CmplxKind))
      endif
!     
!     ===============================================================
!     Calculate the DOS of the multiple scattering term for e = ec,
!     and add the single site results (at er, not ec) to it.
!     ----------------------------------------------------------------
      call computeMSPDOS(is,ec)
!     ----------------------------------------------------------------
      do id = 1, LocalNumAtoms
         do ia = 1, getLocalNumSpecies(id)
            do js = 1, ns_sqr  ! Note: If n_spin_cant = 1, ks = 1
               ks = max(js,is)
               ms_dos_ws(ks) = sfac*getMSCellDOS(ks,id,ia)
               ms_dos_mt(ks) = sfac*getMSMTSphereDOS(ks,id,ia)
               pdos_ws => getMSCellPDOS(ks,id,ia)
               pdos_mt => getMSMTSpherePDOS(ks,id,ia)
               do kl = 1, kmax_phi(id)
                  PartialDOS(id)%ms_pdos_ws(kl,ks,ia) = sfac*pdos_ws(kl)
                  PartialDOS(id)%ms_pdos_mt(kl,ks,ia) = sfac*pdos_mt(kl)
               enddo
            enddo
!
            if (isASAPotential()) then
               ms_dos_ws = ms_dos_mt
               PartialDOS(id)%ms_pdos_ws(:,:,ia) = PartialDOS(id)%ms_pdos_mt(:,:,ia)
            endif
            if (n_spin_cant == 2) then
               PartialDOS(id)%dos_mt(1,ia) = ms_dos_mt(1) + PartialDOS(id)%ss_dos_mt(1,ia)
               PartialDOS(id)%dos_mt(2,ia) = ms_dos_mt(4) + PartialDOS(id)%ss_dos_mt(2,ia)
               PartialDOS(id)%dos_ws(1,ia) = ms_dos_ws(1) + PartialDOS(id)%ss_dos_ws(1,ia)
               PartialDOS(id)%dos_ws(2,ia) = ms_dos_ws(4) + PartialDOS(id)%ss_dos_ws(2,ia)
               do kl = 1, kmax_phi(id)
                  PartialDOS(id)%partial_dos_mt(kl,2,ia) = PartialDOS(id)%ms_pdos_mt(kl,4,ia) + &
                                                           PartialDOS(id)%ss_pdos_mt(kl,2,ia)
                  PartialDOS(id)%partial_dos_ws(kl,1,ia) = PartialDOS(id)%ms_pdos_ws(kl,1,ia) + &
                                                           PartialDOS(id)%ss_pdos_ws(kl,1,ia)
                  PartialDOS(id)%partial_dos_ws(kl,2,ia) = PartialDOS(id)%ms_pdos_ws(kl,4,ia) + &
                                                           PartialDOS(id)%ss_pdos_ws(kl,2,ia)
               enddo
            else
               PartialDOS(id)%dos_mt(is,ia) = ms_dos_mt(is) + PartialDOS(id)%ss_dos_mt(is,ia)
               PartialDOS(id)%dos_ws(is,ia) = ms_dos_ws(is) + PartialDOS(id)%ss_dos_ws(is,ia)
               do kl = 1, kmax_phi(id)
                  PartialDOS(id)%partial_dos_mt(kl,is,ia) = PartialDOS(id)%ms_pdos_mt(kl,is,ia) + &
                                                            PartialDOS(id)%ss_pdos_mt(kl,is,ia)
                  PartialDOS(id)%partial_dos_ws(kl,is,ia) = PartialDOS(id)%ms_pdos_ws(kl,is,ia) + &
                                                            PartialDOS(id)%ss_pdos_ws(kl,is,ia)
               enddo
            endif
         enddo
      enddo
   enddo
!  -------------------------------------------------------------------
   call endMSSolver()
!  -------------------------------------------------------------------
!
   do is = 1, n_spin_pola
      ks = (2*n_spin_cant-1)*is - 2*(n_spin_cant-1)
      do id = 1, LocalNumAtoms
         do ia = 1, getLocalNumSpecies(id)
            if ( (MyPEinEGroup == 0 .and. GlobalNumAtoms < 10) .or. node_print_level >= 0 ) then
               write(2001+id,'(3a,i1,a,f13.8)')'Atom: ',getLocalAtomName(id,ia),', Spin: ',is, &
                                               ', Partial Phase Shift and DOS at Energy = ',er
               write(2001+id,'(a,2f13.8)')'Total DOS in Cell and MT Volumes = ',PartialDOS(id)%dos_ws(is,ia), &
                                          PartialDOS(id)%dos_mt(is,ia)
               write(2001+id,'(a)')'=========================================================================='
               write(2001+id,'(a)')'  l     m         Partial_Phase         Partial_DOS         Partial_DOS_MT'
               do kl = 1, kmax_phi(id)
                  write(2001+id,'(i3,3x,i3,3(9x,f12.8))')lofk(kl),mofk(kl),PartialDOS(id)%phase_shift(kl,is,ia), &
                                  PartialDOS(id)%partial_dos_ws(kl,is,ia),PartialDOS(id)%partial_dos_mt(kl,is,ia)
               enddo
            endif
            if ( node_print_level >= 0) then
               write(6,'(/,3a,i1,a,f11.8,a)')'Atom: ',getLocalAtomName(id,ia),', Spin: ',is, &
                                             ', Total Atom in Crystal DOS at Energy = ',er,':'
               check_ws = ZERO; check_mt = ZERO
               do kl = kmax_phi(id), 1, -1
                  check_ws = check_ws + PartialDOS(id)%ms_pdos_ws(kl,ks,ia)
                  check_mt = check_mt + PartialDOS(id)%ms_pdos_mt(kl,ks,ia)
               enddo
               write(6,'(a,2(f13.8,a))')'MS_DOS          = ',check_ws,'(ws)',check_mt,'(mt)'
               check_ws = ZERO; check_mt = ZERO
               do kl = kmax_phi(id), 1, -1
                  check_ws = check_ws + PartialDOS(id)%partial_dos_ws(kl,is,ia)
                  check_mt = check_mt + PartialDOS(id)%partial_dos_mt(kl,is,ia)
               enddo
               write(6,'(a,2(f13.8,a))')'MS_DOS + SS_DOS = ',check_ws,'(ws)',check_mt,'(mt)'
            endif
         enddo
      enddo
   enddo
!
   if ( (MyPEinEGroup == 0 .and. GlobalNumAtoms < 10) .or. node_print_level >= 0 ) then
      iprint = 1
   else
      iprint = 0
   endif
!
   do id = 1, LocalNumAtoms
      do ia = 1, getLocalNumSpecies(id)
!        -------------------------------------------------------------
         call gaspari_gyorffy_formula(kmax_phi_max,kmax_phi(id),n_spin_pola, &
                                      getLocalAtomicNumber(id,ia),           &
                                      PartialDOS(id)%ss_dos_mt(1,ia),        &
                                      PartialDOS(id)%dos_mt(1,ia),           &
                                      PartialDOS(id)%phase_shift(1,1,ia),    &
                                      PartialDOS(id)%ss_pdos_mt(1,1,ia),     &
                                      PartialDOS(id)%partial_dos_mt(1,1,ia),iprint)
!        -------------------------------------------------------------
      enddo
   enddo
!
   if ( (MyPEinEGroup == 0 .and. GlobalNumAtoms < 10) .or. node_print_level >= 0 ) then
      do id = 1, LocalNumAtoms
         close(unit=1001+id)
         close(unit=2001+id)
      enddo
   endif
!
   if (node_print_level >= 0) then
      write(6,'(/,a,i5)')'Number of single site scattering solver calls:', &
                         NumCalls_SS
      write(6,'(a,f12.3,a,/)')'Total of single site scattering solver timing:', &
                         Timing_SS,' (sec)'
   endif
!
   do id = 1, LocalNumAtoms
      deallocate( PartialDOS(id)%dos_ws )
      deallocate( PartialDOS(id)%dos_mt )
      deallocate( PartialDOS(id)%ss_dos_ws )
      deallocate( PartialDOS(id)%ss_dos_mt )
      deallocate( PartialDOS(id)%phase_shift )
      deallocate( PartialDOS(id)%ss_pdos_mt )
      deallocate( PartialDOS(id)%ss_pdos_ws )
      deallocate( PartialDOS(id)%ms_pdos_mt )
      deallocate( PartialDOS(id)%ms_pdos_ws )
      deallocate( PartialDOS(id)%partial_dos_mt )
      deallocate( PartialDOS(id)%partial_dos_ws )
   enddo
   deallocate( PartialDOS )
!
   if (RelativisticFlag == 2) then
!     ----------------------------------------------------------------
      call endRelScatterer()
      call endRelGreenFunction()
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call endSSSolver()
!     ----------------------------------------------------------------
   endif
!
   end subroutine calPartialDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function adjustEnergy_r(is,e) result(energy)
!  ===================================================================
   use PotentialModule, only : getVdif
   use ScfDataModule, only : isInterstitialElectronPolarized
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind), pointer :: p_vdif(:)
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind) :: energy, vdif
!
   if (isInterstitialElectronPolarized()) then
      p_vdif => getVdif()
      vdif = p_vdif(1)
   else
      vdif = ZERO
   endif
!
   energy = e - (is-1)*vdif
!
   end function adjustEnergy_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function adjustEnergy_c(is,e) result(energy)
!  ===================================================================
   use PotentialModule, only : getVdif
   use ScfDataModule, only : isInterstitialElectronPolarized
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind), pointer :: p_vdif(:)
   real (kind=RealKind) :: vdif
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: energy
!
   if (isInterstitialElectronPolarized()) then
      p_vdif => getVdif()
      vdif = p_vdif(1)
   else
      vdif = ZERO
   endif
!
   energy = e - (is-1)*vdif
!
   end function adjustEnergy_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function returnSingleSitePS(info,e) result(tps)
!  ===================================================================
   use SSSolverModule, only : solveSingleScattering, computePhaseShift, getPhaseShift
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: kl, kmax_phi, is, id, atom
!
   real (kind=RealKind), intent(in) :: e
!
   complex (kind=CmplxKind) :: energy
!
   real (kind=RealKind) :: tps, t0
   real (kind=RealKind), pointer :: ps(:)
!
   is = info(1); id = info(2); atom = info(3)
   energy = adjustEnergy(is,e)
!
   if (atom < 1) then
      call ErrorHandler('returnSingleSitePS','atom < 1',atom)
   endif
!
   t0 = getTime()
!  -------------------------------------------------------------------
   call solveSingleScattering(spin=is,site=id,atom=atom,e=energy,vshift=CZERO)
!  -------------------------------------------------------------------
   Timing_SS = Timing_SS + (getTime() - t0)
   NumCalls_SS = NumCalls_SS + 1
!
   call computePhaseShift(spin=is,site=id,atom=atom)
   ps => getPhaseShift(spin=is,site=id,atom=atom)
!  -------------------------------------------------------------------
   kmax_phi = (lmax_phi(id)+1)**2
   tps = ZERO
   do kl = 1, kmax_phi
      tps = tps + ps(kl)
   enddo
!
   end function returnSingleSitePS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function returnSingleSiteDOS(info,e,aux,rfac,redundant) result(ssdos)
!  ===================================================================
!
!  This function returns single site DOS for a given energy on the real 
!  energy axis.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!
!  In the no-spin-polarized case, a factor of 2 is included in dos.
!  ===================================================================
   use PhysParamModule, only : Boltzmann
!
   use GroupCommModule, only : GlobalSumInGroup
!
   use ScfDataModule, only : Temperature
!
   use RadialGridModule, only : getGrid
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use SSSolverModule, only : solveSingleScattering, computeDOS, getDOS
   use SSSolverModule, only : getOutsideDOS, computePhaseShift, getPhaseShift
   use SSSolverModule, only : getDOSDerivative
!
   use AtomModule, only : getLocalAtomicNumber, getLocalAtomName
   use AtomModule, only : getLocalNumSpecies
!
   implicit none
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind), intent(in), optional :: rfac
!
   logical, intent(in), optional :: redundant
   logical :: red
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: jmax_dos, kmax_phi, iend, n, kl, jl, ir
   integer (kind=IntKind) :: is, id, print_dos, ia, i, atom
!
   real (kind=RealKind), pointer :: dos(:), dos_mt(:), dos_out(:), tps(:)
   real (kind=RealKind) :: sfac, rmul, t0, ssdos
   real (kind=RealKind), pointer :: ps(:), p1(:)
   real (kind=RealKind), pointer :: msgbuf(:,:)
   real (kind=RealKind), allocatable, target :: wks_loc(:)
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: der_dos_r_jl(:,:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3);  print_dos = info(4)
   sfac= TWO/real(n_spin_pola,kind=RealKind)
!  ===================================================================
!  The following piece of code for multiplying the Fermi-Dirac function
!  needs to be carefully thought since for T = 0, if e > chempot, 
!  getFermiDiracFunc = 0, which causes sfac = 0 and has effects on
!  finding the Fermi energy properly.
!  ===================================================================
   if (Temperature > TEN2m6) then
!     ----------------------------------------------------------------
      sfac = sfac*getFermiDiracFunc(adjustEnergy(is,e),chempot,       &
                                    Boltzmann*Temperature)
!     ----------------------------------------------------------------
   endif
!
   n = getLocalNumSpecies(id)
   allocate(wks_loc(4*n+(4*n+1)*NumPEsInEGroup))
   dos => wks_loc(1:n)
   dos_mt => wks_loc(n+1:2*n)
   dos_out => wks_loc(2*n+1:3*n)
   tps => wks_loc(3*n+1:4*n)
   p1 => wks_loc(4*n+1:4*n+(4*n+1)*NumPEsInEGroup)
   msgbuf => aliasArray2_r(p1,4*n+1,NumPEsInEGroup)
!
   if (present(rfac)) then
      rmul = rfac
   else
      rmul = ONE
   endif
!
   if (present(redundant)) then
      red = redundant
   else
      red = .false.
   endif
!
   energy = adjustEnergy(is,e)
!
   t0 = getTime()
!  -------------------------------------------------------------------
   call solveSingleScattering(spin=is,site=id,atom=atom,e=energy,vshift=CZERO)
!  -------------------------------------------------------------------
   Timing_SS = Timing_SS + (getTime() - t0)
   NumCalls_SS = NumCalls_SS + 1
!
   if (getLocalAtomicNumber(id,1) == 0) then
!     call computeDOS(add_highl_fec=.true.)
      call computeDOS(spin=is,site=id,atom=atom)
   else
      call computeDOS(spin=is,site=id,atom=atom)
   endif
!
   Grid => getGrid(id)
   n = 0
   do ia = 1, getLocalNumSpecies(id)
      if (atom < 0 .or. ia == atom) then
         dos_r_jl => getDOS(spin=is,site=id,atom=ia)
         if (rad_derivative) then
            der_dos_r_jl => getDOSDerivative(spin=is,site=id,atom=ia)
         endif
!        -------------------------------------------------------------
         iend = size(dos_r_jl,1); jmax_dos = size(dos_r_jl,2)
!        -------------------------------------------------------------
         dos(ia) = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                             jmax_dos, 2, dos_r_jl, dos_mt(ia) )
         dos_mt(ia) = sfac*dos_mt(ia)
         dos_out(ia) = sfac*getOutsideDOS(spin=is,site=id,atom=ia)
!        -------------------------------------------------------------
!
         do jl = 1, jmax_dos
            do ir = 1, LastValue(id)%NumRs
               aux(n+ir) = rmul*sfac*dos_r_jl(ir,jl)
            enddo
!           ----------------------------------------------------------
!           call zcopy(LastValue(id)%NumRs,dos_r_jl(1,jl),1,aux(n+1),1)
!           ----------------------------------------------------------
            n = n + LastValue(id)%NumRs
         enddo
         if (rad_derivative) then
            do jl = 1, jmax_dos
               do ir = 1, LastValue(id)%NumRs
                  aux(n+ir) = rmul*sfac*der_dos_r_jl(ir,jl)
               enddo
               n = n + LastValue(id)%NumRs
            enddo
         endif
         aux(n+1) = rmul*dos(ia)
         aux(n+2) = rmul*dos_mt(ia)
         aux(n+3) = rmul*dos(ia)*energy
         aux(n+4) = rmul*dos_out(ia)
         n = n + 4
      endif
   enddo
!
!  -------------------------------------------------------------------
   call computePhaseShift(spin=is,site=id,atom=atom)
!  -------------------------------------------------------------------
   kmax_phi = (lmax_phi(id)+1)**2
   tps = ZERO
   do ia = 1, getLocalNumSpecies(id)
      if (atom < 0 .or. ia == atom) then
         ps => getPhaseShift(spin=is,site=id,atom=ia)
         do kl = 1, kmax_phi
            tps(ia) = tps(ia) + ps(kl)
         enddo
      endif
   enddo
!
   if (print_dos > 0) then
      if (.not.red) then
         msgbuf = ZERO
         msgbuf(1,MyPEinEGroup+1) = real(energy)
         n = 1
         do ia = 1, getLocalNumSpecies(id)
            if (atom < 0 .or. ia == atom) then
               msgbuf(n+1,MyPEinEGroup+1) = dos(ia)
               msgbuf(n+2,MyPEinEGroup+1) = dos_mt(ia)
               msgbuf(n+3,MyPEinEGroup+1) = dos_out(ia)
               msgbuf(n+4,MyPEinEGroup+1) = tps(ia)
               n = n + 4
            endif
         enddo
!        -------------------------------------------------------------
         call GlobalSumInGroup(eGID,msgbuf,n,NumPEsInEGroup)
!        -------------------------------------------------------------
         if ( node_print_level >= 0) then
            if (getLocalNumSpecies(id) == 1) then
               do i = 1, NumPEsInEGroup
!                 write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy),dos,dos_mt,dos_out,tps
                  write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')msgbuf(1:5,i)
               enddo
            else
               do i = 1, NumPEsInEGroup
!                 write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy),dos,dos_mt,dos_out,tps
                  n = 1
                  do ia = 1, getLocalNumSpecies(id)
                     if (atom < 0 .or. ia == atom) then
                        write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8),4x,a)')msgbuf(1,i),msgbuf(n+1,i), &
                              msgbuf(n+2,i), msgbuf(n+3,i), msgbuf(n+4,i), getLocalAtomName(id,ia)
                        n = n + 4
                     endif
                  enddo
               enddo
            endif
         endif
      else if ( node_print_level >= 0) then
         if (getLocalNumSpecies(id) == 1) then
            write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy),dos(1),dos_mt(1), &
                  dos_out(1),tps(1)
         else
            do ia = 1, getLocalNumSpecies(id)
               if (atom < 0 .or. ia == atom) then
                  write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8),4x,a)')real(energy),dos(ia),dos_mt(ia), &
                        dos_out(ia),tps(ia), getLocalAtomName(id,ia)
               endif
            enddo
         endif
      endif
   endif
!
   if (atom < 1) then
      ssdos = dos(1)
   else
      ssdos = dos(atom)
   endif
!
   nullify(msgbuf, dos, dos_mt, dos_out, tps)
   deallocate(wks_loc)
!
   end function returnSingleSiteDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function returnSingleSiteDOSinCP(info,e,aux,cfac) result(dos)
!  ===================================================================
!
!  This function returns single site DOS for a given energy on the real 
!  energy axis.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!
!  In the no-spin-polarized case, a factor of 2 is included in dos.
!  ===================================================================
   use PotentialTypeModule, only : isASAPotential
!
   use RadialGridModule, only : getGrid
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : computeGreenFunction, getGreenFunction
   use SSSolverModule, only : getGreenFunctionDerivative
!
   use AtomModule, only : getLocalNumSpecies
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in), optional :: cfac
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: is, id, atom, ia
   integer (kind=IntKind) :: kmax, jmax, iend, n, ir, jl, l, m, kl, klc
!
   real (kind=RealKind) :: dos, sfac, t0
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy, cmul, greenint, greenint_mt, ede
   complex (kind=CmplxKind), pointer :: green(:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: der_green(:,:), der_dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:), p1(:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3)
   sfac= TWO/real(n_spin_pola,kind=RealKind)
!
   if (present(cfac)) then
      cmul = cfac
   else
      cmul = CONE
   endif
!
!  ===================================================================
!  For energy in the upper complex plane, use spherical solver
!  together with irregular solution
!  ===================================================================
   energy = adjustEnergy(is,e)
!
   t0 = getTime()
!  -------------------------------------------------------------------
   call solveSingleScattering(spin=is,site=id,atom=atom,e=energy,vshift=CZERO, &
                              isSphSolver=.true., useIrrSol='h')
!  -------------------------------------------------------------------
   Timing_SS = Timing_SS + (getTime() - t0)
   NumCalls_SS = NumCalls_SS + 1
!
!  -------------------------------------------------------------------
   call computeGreenFunction(spin=is,site=id,atom=atom)
!  -------------------------------------------------------------------
!
   Grid => getGrid(id)
   iend = LastValue(id)%NumRs
   jmax = LastValue(id)%jmax
   n = 0
   do ia = 1, getLocalNumSpecies(id)
      if (atom < 0 .or. ia == atom) then
!        -------------------------------------------------------------
         green=>getGreenFunction(spin=is,site=id,atom=ia)
!        -------------------------------------------------------------
         iend = size(green,1); kmax = size(green,2)
!        -------------------------------------------------------------
         greenint = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,    &
                                               kmax, 2, green, greenint_mt )
!        -------------------------------------------------------------
         greenint_mt = sfac*greenint_mt
!
         if (atom > 0 .or. ia == 1) then
            dos = real(SQRTm1*greenint/PI,kind=RealKind)
         endif
!
         if (isASAPotential()) then
            greenint = greenint_mt
         endif
!
         if (iharris <= 1) then
            ede = energy
         else
            ede = (energy-chempot)
         endif
!
         p1 => aux(n+1:n+iend*jmax)
         dos_r_jl => aliasArray2_c(p1, iend, jmax)
         do jl = 1, jmax
            l = lofj(jl); m = mofj(jl)
            kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
            pcv_x => green(1:iend,kl)
            pcv_y => green(1:iend,klc)
            pcv_z => dos_r_jl(1:iend,jl)
            pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
         enddo
         n = n + iend*jmax
         if (rad_derivative) then
            der_green=>getGreenFunctionDerivative(spin=is,site=id,atom=ia)
            p1 => aux(n+1:n+iend*jmax)
            der_dos_r_jl => aliasArray2_c(p1, iend, jmax)
            do jl = 1, jmax
               l = lofj(jl); m = mofj(jl)
               kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
               pcv_x => der_green(1:iend,kl)
               pcv_y => der_green(1:iend,klc)
               pcv_z => der_dos_r_jl(1:iend,jl)
               pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
            enddo
            n = n + iend*jmax
         endif
         aux(n+1) = SQRTm1*cmul*greenint/PI
         aux(n+2) = SQRTm1*cmul*greenint_mt/PI
         aux(n+3) = SQRTm1*cmul*ede*greenint/PI
!        -------------------------------------------------------------
         greenint = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,       &
                                               kmax, 2, green, greenint_mt, &
                                               truncated=.false. ) - greenint
!        -------------------------------------------------------------
         aux(n+4) = SQRTm1*greenint/PI
         n = n + 4
      endif
   enddo
!
   end function returnSingleSiteDOSinCP
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNLloyd(e_in)                                result(qlloyd)
!  ===================================================================
   use PotentialModule, only : getVdif
   use KreinModule, only : aimag_tol, setLloydStatus, getLloydQ
!
   use RelScattererModule, only : solveRelSST
   use MSSolverModule, only : computeMSGreenFunction
!
   implicit none
!
   real(kind=RealKind), intent(in):: e_in
!
   integer (kind=IntKind) :: is
   real (kind=RealKind) :: qlloyd, e
   complex (kind=CmplxKind) ::e_lloyd
!
   do is = 1, n_spin_pola/n_spin_cant
      e  = adjustEnergy(is,e_in)
!     ----------------------------------------------------------------
!     finds the Lloyd factors for the specified Fermi level
!     ----------------------------------------------------------------
      call setLloydStatus(.true.,.true.)
!     ----------------------------------------------------------------
      e_lloyd = cmplx(e, aimag_tol, kind=CmplxKind )
!
!     ===============================================================
!     Solve the single site scattering problem for a given
!     energy and spin
!     In the future, the relativity code will be integrated into the MSSolver module.
!     ===============================================================
      if (RelativisticFlag == 2) then
!        -------------------------------------------------------------
         call solveRelSST(e_lloyd)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call computeMSGreenFunction(is, e_lloyd)
!        -------------------------------------------------------------
      endif
#ifdef ONESIDED
#ifdef Sync_EiLoop
!     ----------------------------------------------------------------
      call syncAllPEs()
!     ----------------------------------------------------------------
#endif
#endif
!     ----------------------------------------------------------------
      Lloyd_lastE(1:n_spin_pola)=getLloydQ()
!     ----------------------------------------------------------------
      call setLloydStatus(.false.,.false.)
!     ----------------------------------------------------------------
   enddo
!
   if ( n_spin_pola==1 ) then
      qlloyd = TWO*real(-sqrtm1*Lloyd_lastE(1),kind=RealKind)
   else
      qlloyd = real(-sqrtm1*(Lloyd_lastE(1) + Lloyd_lastE(2)),kind=RealKind)
   endif
!
   end function getNLloyd
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine findLloydEf(efermi)
!  ===================================================================
   implicit none
!
   logical :: mflag
!
   integer (kind=IntKind) :: n
!
   real (kind=RealKind), intent(inout) :: efermi
   real (kind=RealKind) :: time_ie, rtmp, rtmp2, qtol, etol
   real (kind=RealKind) :: ea, eb, ec, ed, es, qa, qb, qc, qs
   real (kind=RealKind) :: qdif_a, qdif_b, qdif_c, qdif_s
!
   mflag=.true.
!
   ed = efermi+ONE
   qtol = TEN2m10
   etol = FIVE*TEN2m16
   qa = getNLloyd(efermi)
!
   qdif_a = (qa-zvaltss)/GlobalNumAtoms
   if ( abs(qa) < qtol  ) then
      return
   endif
   ea = efermi
!
   if (qdif_a<ZERO) then
      eb = efermi + max(.1d0*abs(qdif_a),0.10d0)
   else
      eb = efermi - max(.1d0*abs(qdif_a),0.10d0)
   endif
   qb = getNLloyd(eb)
   qdif_b = (qb-zvaltss)/GlobalNumAtoms
!
   if ( abs(qdif_b) < qtol  ) then
      efermi=eb
      return ! early exit
   endif
!
   do while ( qdif_a*qdif_b > ZERO  )
      ea = eb
      qa = qb
      qdif_a = qdif_b
      if ( qdif_b > ZERO ) then
         eb = ea - min(.1d0*abs(qdif_b),0.10d0)
      else
         eb = ea + min(.1d0*abs(qdif_b),0.10d0)
      endif
      qb = getNLloyd(eb)
      qdif_b = (qb-zvaltss)/GlobalNumAtoms
   enddo
!
   if ( ea>eb ) then
!      swap
       rtmp = ea; ea = eb; eb = rtmp
       rtmp = qa; qa = qb; qb = rtmp
       qdif_a = (qa-zvaltss)/GlobalNumAtoms
       qdif_b = (qb-zvaltss)/GlobalNumAtoms
   endif
!
#ifdef DEBUG_Lloyd
   write(6,'(a,2f16.8)') ' LloydSum: ',qb/GlobalNumAtoms, eb
#endif
!
   ec = ea 
   qc = qa
   qdif_c = qdif_a
   n = 0
   do while ( .not.(abs(eb-ea) < etol .or. abs(qdif_b)<qtol .or. n>100) )
      if ( qa/=qc .and. qb/=qc) then
         es = ea*qdif_b*qdif_c/(qdif_a-qdif_b)/(qdif_a-qdif_c) + &
              eb*qdif_a*qdif_c/(qdif_b-qdif_a)/(qdif_b-qdif_c) + &
              ec*qdif_a*qdif_b/(qdif_c-qdif_a)/(qdif_c-qdif_b)
      else
         es = eb-qdif_b*(eb-ea)/(qdif_b-qdif_a)
      endif
!
      rtmp2 = (3.0d0*ea +eb)/4.0d0
!
      if ( (.not.(((es > rtmp2) .and. (es < eb)) .or. ((es < rtmp2) .and. (es > eb)))) &
           .or. (mflag .and. (abs(es - eb) >= (abs(eb - ec)*HALF))) &
           .or. ( .not.mflag .and. (abs(es - eb) >= (abs(ec - ed)*HALF))) ) then
         es = (ea+eb)*HALF
         mflag=.true.
      else
         if ((mflag .and. (abs(eb - ec) < etol)) .or. &
             (.not.mflag .and. (abs(ec - ed) < etol))) then
            es = (ea + eb)*HALF
            mflag = .true.
         else
            mflag = .false.
         endif
      endif
!     ------------------
      qs = getNLloyd(es)
!     ------------------
#ifdef DEBUG_Lloyd
      write(6,'(a,2f22.14)') ' LloydSum: ',qs/GlobalNumAtoms, es
#endif
      ed = ec; ec = eb; qc=qb
      qdif_a = (qa-zvaltss)/GlobalNumAtoms
      qdif_s = (qs-zvaltss)/GlobalNumAtoms
      qdif_c = (qc-zvaltss)/GlobalNumAtoms
      if ( qdif_a*qdif_s<0 ) then
         eb = es; qb = qs; qdif_b = qdif_s
      else
         ea = es; qa = qs; qdif_a = qdif_s
      endif
!
      if ( abs(qdif_a)<abs(qdif_b) ) then
!        swap
         rtmp = ea; ea = eb; eb = rtmp
         rtmp = qa; qa = qb; qb = rtmp
         qdif_a = (qa-zvaltss)/GlobalNumAtoms
         qdif_b = (qb-zvaltss)/GlobalNumAtoms
      endif
!
      n = n+1
   enddo
!
   efermi = eb
!  ---------------------
   qb = getNLloyd(eb)
!  ---------------------
   if ( node_print_level>=0 ) then
      write(6,'(a,2d24.16,i4)')"Lloyd: New Fermi Level::",efermi,qb/GlobalNumAtoms,n
   endif
!
   end subroutine findLloydEf
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSingleScatteringIDOS(LowerContour,UpperContour,Ebegin,Eend,relativity)
!  ===================================================================
   use PublicParamDefinitionsModule, only : ButterFly
!
   use GroupCommModule, only : GlobalSumInGroup
!
   use ScfDataModule, only : NumSS_IntEs, getAdaptiveIntegrationMethod, &
                             Contourtype, ErBottom, NumExtraEs
!
   use AdaptIntegrationModule, only : setupAdaptMesh, getAdaptIntegration
   use AdaptIntegrationModule, only : getAuxDataAdaptIntegration
   use AdaptIntegrationModule, only : getUniFormIntegration, getPeakPos,getWeightedIntegration !xianglin
   use AdaptIntegrationModule, only : getUniMeshIntegration
!
   use PolyhedraModule, only : getVolume
!
   use ContourModule, only : getEPoint, getEWeight, getNumEs
!
   implicit none
!
   logical, intent(in), optional :: LowerContour, UpperContour 
   logical :: LC, UC, REL !xianglin
!
   integer (kind=IntKind) :: id, is, ns, info(4), nm, ie, NumEs, ilc !xianglin
   integer (kind=IntKind) :: ia
!
   real (kind=RealKind), intent(in), optional :: Ebegin, Eend
   logical, intent(in), optional :: relativity !xianglin
!
   real (kind=RealKind) :: ssdos_int, IDOS_cell, ps, peak_pos, e0, ps0, ssDOS
   real (kind=RealKind) :: ebot, etop, er, ei
   real (kind=RealKind), allocatable :: xg(:), wg(:)
!
   complex (kind=CmplxKind) :: int_test, ec, es
   complex (kind=CmplxKind), pointer :: p_aux(:)
   complex (kind=CmplxKind), pointer :: EPoint(:)
   complex (kind=CmplxKind), pointer :: EWght(:)
!
   type (ElectroStruct), pointer :: pCurrentValue
!
!  Lower contour with small semi-circle encircle the poles
   integer (kind=IntKind) :: NumLC !xinglin
   integer (kind=IntKind) :: LC_list(2)
   real (kind=RealKind) :: Pole_loc, LC_Width
   LC_list(1)=1
   LC_list(2)=4
   LC_Width=2.d0*Ten2m4
!
   if (present(relativity)) then
      REL = relativity
   else
      REL = .false.
   endif
!
   if (present(LowerContour)) then
      LC = LowerContour
   else
      LC = .false.
   endif
!
   if (present(UpperContour)) then
      UC = UpperContour
   else
      UC = .false.
   endif
!
   if (present(Ebegin)) then
      ebot = Ebegin
   else
      ebot = 1.0d+20
   endif
!
   if (present(Eend)) then
      etop = Eend
   else
      etop = -1.0d+20
   endif
!
   if (.not.LC .and. .not.UC .and. etop < ebot) then
      call ErrorHandler('calSingleScatteringIDOS','ebot > etop',ebot,etop)
   endif
!
!  ===================================================================
!  Calculate the single site Green function and DOS along the real
!  energy axis between (0.0d0,etop)
!  ===================================================================
   if ( node_print_level >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,20x,a)')'***************************************'
      write(6,'(  20x,a)')'* Output from calSingleScatteringIDOS *'
      write(6,'(  20x,a)')'***************************************'
   endif
!
   do id =  1, LocalNumAtoms
!     ----------------------------------------------------------------
      call zeroElectroStruct(ssLastValue(id)) ! Use ssLastValue as temporary space
!     ----------------------------------------------------------------
   enddo
!
   if (REL) then !added by xianglin
      if (LC) then
         if (ContourType == ButterFly) then
            call ErrorHandler('calSingleScatteringIDOS','ButterFly ContourType not implemented')
         endif
         NumEs = NumExtraEs
         if (NumEs < 2) then
            call ErrorHandler('calSingleScatteringIDOS','On lower contour, NumEs < 2',NumEs)
         endif
         allocate(EPoint(NumEs), EWght(NumEs), xg(NumEs), wg(NumEs))
!        ----------------------------------------------------------
         call gauleg(-ONE, ONE, xg, wg, NumEs)
!        ----------------------------------------------------------
         do is = 1, 1
            do id =  1, LocalNumAtoms !this part only works for less than 2+2 poles at negative energy
               ia = 1 ! This is temporarily set to 1. It should be a loop over number of species
               if (MatrixPoles(id)%NumPoles(ia,1)==2) then !A more general method is needed
                  NumLC=1
               else if (MatrixPoles(id)%NumPoles(ia,1)==4) then
                  NumLC=2
               else if (MatrixPoles(id)%NumPoles(ia,1)==0) then
                  NumLC=0
               else 
                  call ErrorHandler('calSingleScatteringIDOS','NumPoles .ne. 2 or 4')
               endif
!
               if (NumLC .ne. 0) then
                  do ilc = 1, NumLC
                     er = abs(LC_Width)*HALF
                     ec = PI*HALF*SQRTm1
                     Pole_loc = MatrixPoles(id)%Poles(LC_list(ilc),ia,1)
                     do ie = 1, NumEs
                        es = er*exp(ec*(ONE-xg(ie)))
                        EPoint(ie)=Pole_loc + es
                        EWght(ie)=-ec*es*wg(ie)
                     enddo
!           print*,'Pole_loc=',Pole_loc,"Epoint for LC: ", EPoint
                     pCurrentValue => ssLastValue(id) ! Use ssLastValue space for temporary working space.
!                    int_test = CZERO
                     do ie = MyPEinEGroup+1, NumEs, NumPEsInEGroup
!                    int_test = int_test + EPoint(ie)*EWght(ie)
                        er = real(EPoint(ie),kind=RealKind)
                        ei = aimag(EPoint(ie))
!                       =======================================================
!                       Solve the single scattering problem e = (er,ei).
!                       In this case, the calculated Green function is Z*tau*Z - Z*J
!                       =======================================================
                        info(1) = is; info(2) = id; info(3) = -1; info(4) = 0
!                       ----------------------------------------------
                        ssDOS = returnRelSingleSiteDOSinCP(info,EPoint(ie),wk_dos,EWght(ie))
                        if (is_Bxyz) then
                           call calElectroStruct(info,4,wk_dos,pCurrentValue,.true.)
                        else
                           call calElectroStruct(info,1,wk_dos,pCurrentValue,.true.)
                        endif
!                       call calElectroStruct(info,n_spin_cant,wk_dos,pCurrentValue)
!                       ----------------------------------------------
                        call addElectroStruct(CONE,pCurrentValue,ssIntegrValue(id))
!                       ----------------------------------------------
                     enddo
!                    ----------------------------------------------------
!                    call GlobalSumInGroup(eGID,int_test)
!                    ----------------------------------------------------
!                    if ( node_print_level >= 0) then
!                       write(6,'(a,2d15.8,a,d15.8)')'Int[e] = ',int_test,', versus ',HALF*(Etop**2-ErBottom**2)
!                    endif
                  enddo
               endif !if NumLC 
            enddo !LocalNumAtom loop
!           ----------------------------------------------------------
            call sumESVAL(is,ssIntegrValue,eGID) ! Sum over the processors to get the integrated DOS
                                                 ! up to e_loc = NumEsOnMyProc - NumRedunEs
                                                 ! For e_loc > NumEsOnMyProc - NumRedunEs, the
                                                 ! energy integration is only carried in calElectroStruct
                                                 ! on local processor
!           ----------------------------------------------------------
!           if ( node_print_level >= 0) then
!              if ( MyPE == 0 ) then
!                 do id =  1, LocalNumAtoms
!                    do ia = 1, ssIntegrValue(id)%NumSpecies
!                    write(6,'(/,a)')'================================================================================'
!                       write(6,'(a,2i4,2x,a,2i4,f14.8)')'The integrated DOS in atomic cell over the lower contour = ', &
!                             id,ia,'ilc=',ilc,real(ssIntegrValue(1)%dos(1,ia),RealKind)
!                       write(6,'(a,2i4,2x,a,2i4,f14.8)')'The integrated DOS in muffin-tin  over the lower contour = ', &
!                             id,ia,'ilc=',ilc,real(ssIntegrValue(1)%dos_mt(1,ia),RealKind)
!                       write(6,'(a,/)')'================================================================================'
!                    enddo
!                 enddo
!                 do id=1,LocalNumAtoms
!                    do ia = 1, ssIntegrValue(id)%NumSpecies
!                       print*,'The integrated DOS in atomic cell over the lower contour = ', &
!                                id,ia,'NumLC=',NumLC,real(ssIntegrValue(id)%dos(is),RealKind)
!                       print*,'The integrated DOS in muffin-tin  over the lower contour = ', &
!                                id,ia,'NumLC=',NumLC,real(ssIntegrValue(id)%dos_mt(is),RealKind)
!                    enddo
!                 enddo
!              endif
!           endif
         enddo
         if (ContourType /= ButterFly) then
            deallocate(EPoint, EWght, xg, wg)
         endif

!        print*,"findResonce used"
!        do id = 1, LocalNumAtoms
!           do is = 1, 1
!              do ie = 1, MatrixPoles(id)%NumPoles(ia,1)
!                 call zeroElectroStruct(ssLastValue(id))
!                 do ia = 1, ssLastValue(id)%NumSpecies
!                    info(1) = is; info(2) = id; info(3) = ia; info(4) = 1
!                    ssDOS = returnRelSingleSiteDOS_pole(info,real(MatrixPoles(id)%Poles(ie,ia,1),RealKind),wk_dos)
!                    if (is_Bxyz) then
!                       call calElectroStruct(info,4,wk_dos,ssLastValue(id),.true.)
!                    else
!                       call calElectroStruct(info,1,wk_dos,ssLastValue(id),.true.)
!                    endif
!                 enddo
!                 call addElectroStruct(CONE,ssLastValue(id),ssIntegrValue(id))
!              enddo
!           enddo
!        enddo
      else !(if LC)
         exc = ZERO
         e0 = Ten2m8 ! initial energy slightly above e = 0.
         e0 = min(e0,abs(ebot))
         do id = 1, LocalNumAtoms
            do ia = 1, ssLastValue(id)%NumSpecies
               do is = 1, 1!n_spin_cant*n_spin_cant
                  info(1) = is; info(2) = id; info(3) = ia; info(4) = 1
                  if (getAdaptiveIntegrationMethod() == 1) then
                     call ErrorHandler('calSingleScatteringIDOS',&
                                       'relativistic AdaptiveIntegration not implemented')
                  else if (isPole_plus) then
                     ssdos_int = getWeightedIntegration(NumSS_IntEs,ebot,etop,info, &
                                                        returnRelSingleSiteDOS,     &
                                                        MatrixPoles(id)%NumPoles_plus(ia,1),MatrixPoles(id)%Poles_plus(:,ia,1))
                     nm=NumSS_IntEs
                  else
                     ssdos_int = getUniFormIntegration(NumSS_IntEs,ebot,etop,info,returnRelSingleSiteDOS,nm)
 !                   peak_pos = peak_pos + (3.d0-2.d0*is)*getPeakPos()
!                    exc(ia,id) = exc(ia,id) + (3.d0-2.d0*is)*getPeakPos()
                  endif
                  if ( node_print_level >= 0) then
                     write(6,'(a)')   '================================================================================'
                     write(6,'(a,i4)')'Number of mesh points for the integration: ',nm
                  endif
!                 ----------------------------------------------------
                  p_aux => getAuxDataAdaptIntegration()
!                 ----------------------------------------------------
                  if (is_Bxyz) then
                     call calElectroStruct(info,4,p_aux,ssLastValue(id),.true.)
                  else
                     call calElectroStruct(info,1,p_aux,ssLastValue(id),.true.)
                  endif
               enddo
            enddo
!           -------------------------------------------------------
            call addElectroStruct(CONE,ssLastValue(id),ssIntegrValue(id))
!           -------------------------------------------------------
         enddo
      endif !end if REL
   else if (LC) then !if nonrelativistic and lower contour needed
      if (ContourType == ButterFly) then
         EPoint => getEPoint()
         EWght => getEWeight()
         NumEs = getNumEs(LowContour=.true.,HighContour=.false.)
         if (real(EPoint(1),RealKind) > ZERO) then
            call ErrorHandler('calSingleScatteringIDOS','On lower contour, e(1) > 0.0',EPoint(1))
         endif
      else
         NumEs = NumExtraEs
         if (NumEs < 2) then
            call ErrorHandler('calSingleScatteringIDOS','On lower contour, NumEs < 2',NumEs)
         endif
         allocate(EPoint(NumEs), EWght(NumEs), xg(NumEs), wg(NumEs))
!        -------------------------------------------------------------
         call gauleg(-ONE, ONE, xg, wg, NumEs)
!        -------------------------------------------------------------
         er = abs(ErBottom)*HALF
         ec = PI*HALF*SQRTm1
         do ie = 1, NumEs
            es = er*exp(ec*(ONE-xg(ie)))
            EPoint(ie)=HALF*ErBottom + es
            EWght(ie)=-ec*es*wg(ie)
         enddo
      endif
      do is = 1, n_spin_pola
         ns = (2*n_spin_cant-1)*is - (n_spin_cant-1)*2  ! ns = 1 or 2, if n_spin_cant = 1
                                                        ! ns = 1 or 4, if n_spin_cant = 2
         do id =  1, LocalNumAtoms
            pCurrentValue => ssLastValue(id) ! Use ssLastValue space for temporary working space.
!           int_test = CZERO
            do ie = MyPEinEGroup+1, NumEs, NumPEsInEGroup
!              int_test = int_test + EPoint(ie)*EWght(ie)
               er = real(EPoint(ie),kind=RealKind)
               ei = aimag(EPoint(ie))
!              =======================================================
!              Solve the single scattering problem e = (er,ei).
!              In this case, the calculated Green function is Z*tau*Z - Z*J
!              =======================================================
               info(1) = is; info(2) = id; info(3) = -1; info(4) = 0
!              -------------------------------------------------------
               ssDOS = returnSingleSiteDOSinCP(info,EPoint(ie),wk_dos,EWght(ie))
               call calElectroStruct(info,1,wk_dos,pCurrentValue,ss_int=.true.)
!              -------------------------------------------------------
               call addElectroStruct(CONE,pCurrentValue,ssIntegrValue(id),ns)
!              -------------------------------------------------------
            enddo
!           ----------------------------------------------------------
!           call GlobalSumInGroup(eGID,int_test)
!           ----------------------------------------------------------
!           if ( node_print_level >= 0) then
!              write(6,'(a,2d15.8,a,d15.8)')'Int[e] = ',int_test,', versus ',HALF*(Etop**2-ErBottom**2)
!           endif
         enddo
      enddo
!     ----------------------------------------------------------------
      call sumESVAL_r(ssIntegrValue,eGID) ! Sum over the processors to get the integrated DOS
                                          ! up to e_loc = NumEsOnMyProc - NumRedunEs
                                          ! For e_loc > NumEsOnMyProc - NumRedunEs, the
                                          ! energy integration is only carried in calElectroStruct
                                          ! on local processor
!     ----------------------------------------------------------------
      if ( node_print_level >= 0) then
         do is = 1, n_spin_pola
            ns = (2*n_spin_cant-1)*is - (n_spin_cant-1)*2  ! ns = 1 or 2, if n_spin_cant = 1
                                                           ! ns = 1 or 4, if n_spin_cant = 2
            do id =  1, LocalNumAtoms
               do ia = 1, ssIntegrValue(id)%NumSpecies
                  write(6,'(/,a)')'================================================================================'
                  write(6,'(a,3i4,2x,f14.8)')'The integrated DOS in atomic cell over the lower contour = ', &
                        id,ia,is,real(ssIntegrValue(id)%dos(ns,ia),RealKind)
                  write(6,'(a,3i4,2x,f14.8)')'The integrated DOS in muffin-tin  over the lower contour = ', &
                        id,ia,is,real(ssIntegrValue(id)%dos_mt(ns,ia),RealKind)
                  write(6,'(a,/)')'================================================================================'
               enddo
            enddo
         enddo
      endif
      if (ContourType /= ButterFly) then
         deallocate(EPoint, EWght, xg, wg)
      endif
   else if (UC) then !if nonrelativistic and upper contour needed
      EPoint => getEPoint()
      EWght => getEWeight()
      do is = 1, n_spin_pola
         ns = (2*n_spin_cant-1)*is - (n_spin_cant-1)*2  ! ns = 1 or 2, if n_spin_cant = 1
                                                        ! ns = 1 or 4, if n_spin_cant = 2
         do id =  1, LocalNumAtoms
            pCurrentValue => ssLastValue(id) ! Use ssLastValue space for temporary working space.
!           int_test = CZERO
            do ie = getNumEs(LowContour=.true.,HighContour=.false.)+MyPEinEGroup+1, getNumEs(), NumPEsInEGroup
!              int_test = int_test + EPoint(ie)*EWght(ie)
               er = real(EPoint(ie),kind=RealKind)
               ei = aimag(EPoint(ie))
!              =======================================================
!              Solve the single scattering problem e = (er,ei).
!              In this case, the calculated Green function is Z*tau*Z - Z*J
!              =======================================================
               info(1) = is; info(2) = id; info(3) = -1; info(4) = 0
!              -------------------------------------------------------
               ssDOS = returnSingleSiteDOSinCP(info,EPoint(ie),wk_dos,EWght(ie))
               call calElectroStruct(info,1,wk_dos,pCurrentValue,ss_int=.true.)
               call addElectroStruct(CONE,pCurrentValue,ssIntegrValue(id),ns)
!              -------------------------------------------------------
            enddo
!           ----------------------------------------------------------
!           call GlobalSumInGroup(eGID,int_test)
!           ----------------------------------------------------------
!           if ( node_print_level >= 0) then
!              write(6,'(a,2d15.8,a,d15.8)')'Int[e] = ',int_test,', versus ',HALF*Etop**2
!           endif
         enddo
      enddo
!     ----------------------------------------------------------------
      call sumESVAL_r(ssIntegrValue,eGID) ! Sum over the processors to get the integrated DOS
                                          ! up to e_loc = NumEsOnMyProc - NumRedunEs
                                          ! For e_loc > NumEsOnMyProc - NumRedunEs, the
                                          ! energy integration is only carried in calElectroStruct
                                          ! on local processor
!     ----------------------------------------------------------------
      if (node_print_level >= 0) then
         do is = 1, n_spin_pola
            ns = (2*n_spin_cant-1)*is - (n_spin_cant-1)*2  ! ns = 1 or 2, if n_spin_cant = 1
                                                           ! ns = 1 or 4, if n_spin_cant = 2
            do id =  1, LocalNumAtoms
               write(6,'(a,/)')'================================================================================'
               do ia = 1, ssIntegrValue(id)%NumSpecies
                  write(6,'(a,3i4,2x,f14.8)')'The integrated DOS in atomic cell over the higher contour = ', &
                        id,ia,is,real(ssIntegrValue(id)%dos(ns,ia),RealKind)
                  write(6,'(a,3i4,2x,f14.8)')'The integrated DOS in muffin-tin  over the higher contour = ', &
                        id,ia,is,real(ssIntegrValue(id)%dos_mt(ns,ia),RealKind)
               enddo
               write(6,'(a,/)')'================================================================================'
            enddo
         enddo
      endif
   else !if nonrelativistic and integration of G_ss on real axis is needed
      exc = ZERO
      e0 = 0.00001d0 ! initial energy slightly above e = 0.
      e0 = min(e0,abs(ebot))
!
      exc = ZERO
      do id =  1, LocalNumAtoms
         do is = 1, n_spin_pola
            ns = (2*n_spin_cant-1)*is - (n_spin_cant-1)*2  ! ns = 1 or 2, if n_spin_cant = 1
                                                           ! ns = 1 or 4, if n_spin_cant = 2
            do ia = 1, ssLastValue(id)%NumSpecies
               info(1) = is; info(2) = id; info(3) = ia; info(4) = 1
!              -------------------------------------------------------
               ps0 = returnSingleSitePS(info,e0)
!              -------------------------------------------------------
               if ( node_print_level >= 0) then
                  write(6,'(/,3(a,i2))')'is = ',is,', id = ',id,', ia = ',ia
                  write(6,'(a,f11.8,a)')'Integration over the real energy interval:  [0.001,',etop,']'
                  write(6,'(a,i5)')'The number of processors employed for parallelizing the DOS calculation: ',NumPEsInEGroup
                  write(6,'(a)')   '=========================================================================================='
                  write(6,'(4x,a)')'Energy    Single_Site_DOS_ws   Single_Site_DOS_mt    DOS_Outside     Total_Phase_Shift'
                  write(6,'(a)')   '------------------------------------------------------------------------------------------'
               endif
               if (getAdaptiveIntegrationMethod() == 0) then
!                 ----------------------------------------------------
                  ssdos_int = getUniMeshIntegration(NumSS_IntEs,ebot,etop,info,returnSingleSiteDOS,nm)
!                 ----------------------------------------------------
               else if (getAdaptiveIntegrationMethod() == 1) then
!                 ----------------------------------------------------
                  call setupAdaptMesh(ebot,etop,NumSS_IntEs,info,returnSingleSitePS)
                  ssdos_int = getAdaptIntegration(info,returnSingleSiteDOS,nm)
!                 ----------------------------------------------------
               else if (getAdaptiveIntegrationMethod() == 2) then
!                 ----------------------------------------------------
                  ssdos_int = getUniformIntegration(NumSS_IntEs,ebot,etop,info,returnSingleSiteDOS,nm)
!                 ----------------------------------------------------
               else
!                 ----------------------------------------------------
                  call ErrorHandler('calSingleScatteringIDOS','Unknown integration scheme',&
                                    getAdaptiveIntegrationMethod())
!                 ----------------------------------------------------
               endif
               exc(ia,id) = exc(ia,id) + (3.d0-2.d0*is)*getPeakPos()
!
               if ( node_print_level >= 0) then
                  write(6,'(a)')   '=========================================================================================='
                  write(6,'(a,i4)')'Number of mesh points for the integration: ',nm
               endif
!              -------------------------------------------------------
               p_aux => getAuxDataAdaptIntegration()
!              -------------------------------------------------------
!              write(6,'(a,i5)')'size of p_aux = ',size(p_aux)
!              ----------------------------------------------------------
               call calElectroStruct(info,1,p_aux,ssLastValue(id),ss_int=.true.)
               ps = returnSingleSitePS(info,etop) - ps0  ! Relative to the phase shift at energy = e0.
!              -------------------------------------------------------
               if ( node_print_level >= 0) then
                  write(6,'(a,f12.8,a,d15.8)')'At energy = ',etop,     &
                                              ': Single site IDOS given by Green function  =',ssdos_int
                  write(6,'(26x,a,d15.8)')      'Integrated DOS outside the atomic cell    =', &
                                              ssIDOS_out(id)%rarray2(ns,ia)
!                 ====================================================
                  write(6,'(26x,a,d15.8)')      'Single site sum. of partial phase shifts  =',ps
                  write(6,'(26x,a,d15.8)')      'Integrated DOS in space due to single site=',(2/n_spin_pola)*ps/PI
                  write(6,'(26x,a,d15.8)')      'Free electron DOS in the atomic cell      =', &
                                          (2/n_spin_pola)*getVolume(id)*sqrt(etop**3)/(6.0d0*PI**2)
                  IDOS_cell = (2/n_spin_pola)*(ps+getVolume(id)*sqrt(etop**3)/(6.0d0*PI))/PI - &
                              ssIDOS_out(id)%rarray2(ns,ia)
                  write(6,'(26x,a,d15.8)')      'Single site {phase shift sum-OutsideIDOS} =',IDOS_cell
!                 ====================================================
               endif
            enddo
!           ----------------------------------------------------------
            call addElectroStruct(CONE,ssLastValue(id),ssIntegrValue(id),ns)
!           ----------------------------------------------------------
         enddo
      enddo
   endif
!
   end subroutine calSingleScatteringIDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calMultipleScatteringIDOS(useIrregularSolution,relativity)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use ProcMappingModule, only : getNumEsOnMyProc, getEnergyIndex,    &
                                 getNumRedundantEsOnMyProc
!
   use PublicParamDefinitionsModule, only : ButterFly
!
   use ContourModule, only : getEPoint, getEWeight, getNumEs
!
   use ScfDataModule, only : isLloyd, getLloydMode, ContourType
!
   use PotentialTypeModule, only : isFullPotential
!
   use MSSolverModule, only : computeMSGreenFunction
   use RelMSSolverModule, only : computeRelMST
!
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   logical, intent(in) :: useIrregularSolution
   logical, intent(in), optional :: relativity
   logical :: REL
!
   integer (kind=IntKind) :: ie, id, is, js, info(4), ia
   integer (kind=IntKind) :: e_loc
   integer (kind=IntKind) :: NumEsOnMyProc, NumRedunEs
!
   real (kind=RealKind) :: time_ie, ssDOS, msDOS(2), ChargeInLowContour(LocalNumAtoms)
!
   complex (kind=CmplxKind), pointer :: EPoint(:)
   complex (kind=CmplxKind), pointer :: EWght(:)
!  complex (kind=CmplxKind) :: test_function, test_result(LocalNumAtoms)
!  complex (kind=CmplxKind) :: test_result_c1(LocalNumAtoms)
!
   type (ElectroStruct), pointer :: pCurrentValue
!
!  if (isSingleSiteCluster) then
!     return
!  endif
!
!  subroutine changed by xianglin
   if (present(relativity)) then
      REL = relativity
   else
      REL = .false.
   endif 

   EPoint => getEPoint()
   EWght => getEWeight()
!
   NumEsOnMyProc = getNumEsOnMyProc()
   NumRedunEs = getNumRedundantEsOnMyProc()
!
   if (REL) then
      do is = 1, n_spin_pola/n_spin_cant
!     test_result = CZERO
         do e_loc = 1,NumEsOnMyProc
            ie = getEnergyIndex(e_loc)
            time_ie = getTime()
            call computeRelMST( adjustEnergy(is,EPoint(ie)) )
            do id =1, LocalNumAtoms
               pCurrentValue => LastValue(id) ! Use LastValue space for temporary working space.
               info(1) = is; info(2) = id; info(3) = -1; info(4) = 0
               msDOS = returnRelMultipleSiteDOS(info,EPoint(ie),wk_dos,EWght(ie))
               call calElectroStruct(info,n_spin_cant,wk_dos,pCurrentValue)
               if (e_loc <= NumEsOnMyProc - NumRedunEs .or. MyPEinEGroup == 0) then
!                 test_result(id) = test_result(id) + test_function*EWght(ie)
                  if (n_spin_cant == 2) then
!                    -------------------------------------------------
                     call addElectroStruct(CONE,pCurrentValue,IntegrValue(id))
!                    -------------------------------------------------
                  else
!                    -------------------------------------------------
                     call addElectroStruct(CONE,pCurrentValue,IntegrValue(id),is)
!                    -------------------------------------------------
                  endif
               else
                  if (n_spin_cant == 2) then
!                    -------------------------------------------------
                     call addElectroStruct(CZERO,pCurrentValue,IntegrValue(id))
!                    -------------------------------------------------
                  else
!                    -------------------------------------------------
                     call addElectroStruct(CZERO,pCurrentValue,IntegrValue(id),is)
!                    -------------------------------------------------
                  endif
               endif
            enddo
            time_ie = getTime()-time_ie
            if (node_print_level >= 0) then
               write(6,'(a,i3,a,f12.5)') 'ie = ', ie,' , Time: ', time_ie
               write(6,'(92(''-''),/)')
            endif
            call FlushFile(6)
         enddo   ! Loop over energy parameter
!        ================================================================
!        Sum over the processors to get the integrated DOS
!        ----------------------------------------------------------------
         call sumESVAL(is,IntegrValue,eGID)
      enddo ! Loop over spin index
   else
      do is = 1, n_spin_pola/n_spin_cant
!        test_result = CZERO
         do e_loc = 1,NumEsOnMyProc
            ie = getEnergyIndex(e_loc)
            time_ie = getTime()
!
!           =============================================================
!           Solve the multiple scattering problem e = energy.
!           In the future, the relativity code will be implemented here
!           =============================================================
            if (useIrregularSolution) then
!              ==========================================================
!              In this case, the calculated Green function is Z*tau*Z - Z*J
!              ----------------------------------------------------------
               call computeMSGreenFunction(is, adjustEnergy(is,EPoint(ie)), &
                                           add_Gs=.true., isSphSolver=.true.)
!              ----------------------------------------------------------
            else if (isFullPotential()) then
!              ==========================================================
!              In this case, the calculated Green function is Z*(tau-t)*Z
!              ----------------------------------------------------------
               call computeMSGreenFunction(is, adjustEnergy(is,EPoint(ie)))
!              ----------------------------------------------------------
            else
!              ==========================================================
!              In this case, the calculated Green function is Z*(tau-t)*Z
!              ----------------------------------------------------------
               call computeMSGreenFunction(is, adjustEnergy(is,EPoint(ie)), &
                                           isSphSolver=.true.)
!              ----------------------------------------------------------
            endif
            do id = 1, LocalNumAtoms
               pCurrentValue => LastValue(id) ! Use LastValue space for temporary working space.
               info(1) = is; info(2) = id; info(3) = -1; info(4) = 0
!              ----------------------------------------------------------
               msDOS = returnMultipleSiteDOS(info,EPoint(ie),wk_dos,EWght(ie))
               call calElectroStruct(info,n_spin_cant,wk_dos,pCurrentValue)
!              ==========================================================
!              Sum over the energy points up to e_loc = NumEsOnMyProc - NumRedunEs
!              on local processor.
!              For e_loc > NumEsOnMyProc - NumRedunEs, the energy integration is
!              carried out with a prefactor 0.0.
!              ==========================================================
!              test_function = CONE + TWO*EPoint(ie) + THREE*EPoint(ie)**2
               if (e_loc <= NumEsOnMyProc - NumRedunEs .or. MyPEinEGroup == 0) then
                  if (n_spin_cant == 1) then
!                    ----------------------------------------------------
                     call addElectroStruct(CONE,pCurrentValue,IntegrValue(id),is)
!                    ----------------------------------------------------
                  else
!                    ----------------------------------------------------
                     call addElectroStruct(CONE,pCurrentValue,IntegrValue(id))
!                    ----------------------------------------------------
                  endif
               else
                  if (n_spin_cant == 1) then
!                    ----------------------------------------------------
                     call addElectroStruct(CZERO,pCurrentValue,IntegrValue(id),is)
!                    ----------------------------------------------------
                  else
!                    ----------------------------------------------------
                     call addElectroStruct(CZERO,pCurrentValue,IntegrValue(id))
!                    ----------------------------------------------------
                  endif
               endif
               if (ContourType == ButterFly .and.                        &
                   ie <= getNumEs(LowContour=.true.,HighContour=.false.)) then
                  ChargeInLowContour(id) = ZERO
                  do ia = 1, IntegrValue(id)%NumSpecies
                     ChargeInLowContour(id) = ChargeInLowContour(id) +   &
                         getLocalSpeciesContent(id,ia)*real(IntegrValue(id)%dos(is,ia),RealKind)
                  enddo
               endif
!
!              =======================================================
!              For testing purposes:
!              Check the single site DOS for the real part of the energy
!              =======================================================
!!             do js = 1, n_spin_cant
!!                info(1) = js; info(2) = id; info(3) = -1; info(4) = 1
!!                ssDOS = returnSingleSiteDOS(info,real(EPoint(ie),kind=RealKind),wk_dos)
!!                write(6,'(a,i2,a,i2,a,d15.8,a,d15.8)')'is = ',max(is,js),', id = ',id, &
!!                     ', e = ',real(EPoint(ie),kind=RealKind),' 0.00000000D+00, Single Site DOS = ',ssDOS
!!                write(6,'(a,i2,a,i2,a,d15.8,a,d15.8)')'is = ',max(is,js),', id = ',id, &
!!                     ', e = ',real(EPoint(ie),kind=RealKind),' 0.00000000D+00, SS DOS + MS DOS = ',ssDOS + msDOS(js)
!!             enddo
            enddo !e_loc loop
!
            time_ie = getTime()-time_ie
            if (node_print_level >= 0) then
               write(6,'(a,i3,a,f12.5)') 'ie = ', ie,' , Time: ', time_ie
               write(6,'(92(''-''),/)')
            endif
            call FlushFile(6)
         enddo   ! Loop over energy parameter
!        -------------------------------------------------------------
!        call GlobalSumInGroup(eGID,test_result,LocalNumAtoms)
!        -------------------------------------------------------------
!
         if (ContourType == ButterFly) then
!           ----------------------------------------------------------
            call GlobalSumInGroup(eGID,ChargeInLowContour,LocalNumAtoms)
!           call GlobalSumInGroup(eGID,test_result_c1,LocalNumAtoms)
!           ----------------------------------------------------------
            if (node_print_level >= 0) then
               write(6,'(/,a)')'========================================================='
               write(6,'(a)')'The integrated DOS of the MS term over the lower contour:'
               do id = 1, LocalNumAtoms
                  write(6,'(4x,2(a,i4),2x,a,f14.8)')'atom = ',id,', spin = ',is,', IDOS = ', &
                                                 ChargeInLowContour(id)
!                 write(6,'(16x,a,2f14.8)')'Integration test = ', test_result_c1(id)
               enddo
               write(6,'(a,/)')'========================================================='
            endif
         endif
!        do id = 1, LocalNumAtoms
!           write(6,'(/,a,2f14.8)')'Full integration test = ', test_result(id)
!        enddo
      enddo ! Loop over spin index
!     ================================================================
!     Sum over the processors to get the integrated DOS
!     ----------------------------------------------------------------
      call sumESVAL_r(IntegrValue,eGID)
!     ----------------------------------------------------------------
   endif
!
   end subroutine calMultipleScatteringIDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calElectroStruct(info,ns,dos_array,CurrentValue,ss_int)
!  ===================================================================
   implicit none !subroutine added by xianglin
!
   integer (kind=IntKind), intent(in) :: info(*), ns
!
   logical, optional, intent(in) :: ss_int
!
   logical :: IntegratedSingleSite, MomRepresentation
!
   integer (kind=IntKind) :: is, id, ia, atom
   integer (kind=IntKind) :: js, n, p0, iend, jmax, kmax, n0
! 
   complex (kind=CmplxKind), intent(in), target :: dos_array(:)
!
   type (ElectroStruct), intent(out) :: CurrentValue
!
   complex (kind=CmplxKind), pointer :: pca_x(:,:), pca_y(:,:), p1(:)
!
   if (present(ss_int)) then
      IntegratedSingleSite = ss_int
   else
      IntegratedSingleSite = .false.
   endif
!
   is = info(1); id = info(2); atom = info(3)
!
   if (ns == 4) then !when ns=4, is=1:4,for relativistic Bxyz calculation. Added by xianglin
      if (is_Bxyz) then
         n0 = 0
         n = CurrentValue%NumRs*CurrentValue%jmax
         do ia = 1, CurrentValue%NumSpecies
            if (atom < 0 .or. atom == ia) then
               do js = 1,4
!                 n0 = (js-1)*(n+4)
!                 ----------------------------------------------------------
                  call zcopy(n,dos_array(n0+1:n0+n),1,CurrentValue%dos_r_jl(1,1,js,ia),1)
!                 ----------------------------------------------------------
                  n0 = n0 + n
                  if (rad_derivative) then
!                 ----------------------------------------------------------
                     call zcopy(n,dos_array(n0+1:n0+n),1,CurrentValue%der_dos_r_jl(1,1,js,ia),1)
!                 ----------------------------------------------------------
                     n0 = n0 + n
                  endif
!                 write(6,'(a,i5)')'Size of dos_array = ',n+4
                  CurrentValue%dos(js,ia) = dos_array(n0+1)
                  CurrentValue%dos_mt(js,ia) = dos_array(n0+2)
                  CurrentValue%evalsum(js,ia) = dos_array(n0+3)
                  if (IntegratedSingleSite) then
                     ssIDOS_out(id)%rarray2(js,ia) = real(dos_array(n0+4),kind=RealKind)
                  endif
                  n0 = n0 + 4
               enddo
            endif
         enddo
      else
         js=is
         n0 = 0
         n = CurrentValue%NumRs*CurrentValue%jmax
         do ia = 1, CurrentValue%NumSpecies
            if (atom < 0 .or. atom == ia) then
!              ----------------------------------------------------------
               call zcopy(n,dos_array(n0+1:n0+n),1,CurrentValue%dos_r_jl(1,1,js,ia),1)
!              ----------------------------------------------------------
               n0 = n0 + n
               if (rad_derivative) then
!                 -------------------------------------------------------
                  call zcopy(n,dos_array(n0+1:n0+n),1,CurrentValue%der_dos_r_jl(1,1,js,ia),1)
!                 -------------------------------------------------------
                  n0 = n0 + n
               endif
!              write(6,'(a,i5)')'Size of dos_array = ',n+4
               CurrentValue%dos(js,ia) = dos_array(n0+1)
               CurrentValue%dos_mt(js,ia) = dos_array(n0+2)
               CurrentValue%evalsum(js,ia) = dos_array(n0+3)
               if (IntegratedSingleSite) then
                  ssIDOS_out(id)%rarray2(js,ia) = real(dos_array(n0+4),kind=RealKind)
               endif
               n0 = n0 + 4
            endif
         enddo
      endif
   else if (ns == 1) then
      if (n_spin_cant == 1) then
         js = is      ! js = 1 or 2
      else
         js = is*is   ! js = 1 or 4
      endif
      n0 = 0
      n = CurrentValue%NumRs*CurrentValue%jmax
      do ia = 1, CurrentValue%NumSpecies
         if (atom < 0 .or. atom == ia) then
!           ----------------------------------------------------------
            call zcopy(n,dos_array(n0+1:n0+n),1,CurrentValue%dos_r_jl(1,1,js,ia),1)
!           ----------------------------------------------------------
            n0 = n0 + n
            if (rad_derivative) then
!              -------------------------------------------------------
               call zcopy(n,dos_array(n0+1:n0+n),1,CurrentValue%der_dos_r_jl(1,1,js,ia),1)
!              -------------------------------------------------------
               n0 = n0 + n
            endif
!           write(6,'(a,i5)')'Size of dos_array = ',n+4
            CurrentValue%dos(js,ia) = dos_array(n0+1)
            CurrentValue%dos_mt(js,ia) = dos_array(n0+2)
            CurrentValue%evalsum(js,ia) = dos_array(n0+3)
            if (IntegratedSingleSite) then
               ssIDOS_out(id)%rarray2(js,ia) = real(dos_array(n0+4),kind=RealKind)
            endif
            n0 = n0 + 4
!           write(6,'(a,3d15.8)')'dos_array(n+1:n+3) = ',real(dos_array(n+1:n+3))*PI
            if (isDensityMatrixNeeded) then
               kmax = CurrentValue%kmax
               p1 => dos_array(n0+1:n0+kmax*kmax)
               pca_x => aliasArray2_c(p1,kmax,kmax)
               pca_y => CurrentValue%density_matrix(:,:,js,ia)
               pca_y = pca_x
               n0 = n0 + kmax*kmax
            endif
         endif
      enddo
   else if (is == 1 .and. ns == 2 .and. n_spin_cant == 2) then
      iend = CurrentValue%NumRs
      jmax = CurrentValue%jmax
      p0 = 0
      n = iend*jmax
      do ia = 1, CurrentValue%NumSpecies
         if (atom < 0 .or. atom == ia) then
            do js = 1, ns*ns
               p1 => dos_array(p0+1:p0+n)
               pca_x => aliasArray2_c(p1,iend,jmax)
               pca_y => CurrentValue%dos_r_jl(:,:,js,ia)
               pca_y = pca_x
               if (rad_derivative) then
                  p0 = p0 + n
                  p1 => dos_array(p0+1:p0+n)
                  pca_x => aliasArray2_c(p1,iend,jmax)
                  pca_y => CurrentValue%der_dos_r_jl(:,:,js,ia)
                  pca_y = pca_x
               endif
               p0 = p0 + n
               CurrentValue%dos(js,ia) = dos_array(p0+1)
               CurrentValue%dos_mt(js,ia) = dos_array(p0+2)
               CurrentValue%evalsum(js,ia) = dos_array(p0+3)
               if (IntegratedSingleSite) then
                  ssIDOS_out(id)%rarray2(js,ia) = real(dos_array(p0+4),kind=RealKind)
               endif
               p0 = p0 + 4
               if (isDensityMatrixNeeded) then
                  kmax = CurrentValue%kmax
                  p1 => dos_array(p0+1:p0+kmax*kmax)
                  pca_x => aliasArray2_c(p1,kmax,kmax)
                  pca_y => CurrentValue%density_matrix(:,:,js,ia)
                  pca_y = pca_x
                  p0 = p0 + kmax*kmax
               endif
            enddo
         endif
      enddo
   else
!     ----------------------------------------------------------------
      call ErrorHandler('calElectroStruct','Invalid is and ns combination',is,ns)
!     ----------------------------------------------------------------
   endif
!
   end subroutine calElectroStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine transformElectroStruct(id,ESV)
!  ===================================================================
   use SpinRotationModule, only : transformDensityMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: iend, jmax, kmax, ia
!
   complex (kind=CmplxKind), pointer :: p1(:), p3(:,:,:)
!
   type (ElectroStruct), intent(inout) :: ESV
!
   if ( n_spin_cant == 1) then
      return
   endif
!
!  ===================================================================
!  If spin canted, transform the dos arrays in the ESV structure by
!  taking the following transform so the the quantity is in moment
!  representation
!  
!        dos(1) = Tr[dos(1:4)]     <- density
!                    --------
!
!        dos(2) = Tr[dos(1:4)*wx]  <- moment density along x
!                    -------- --
!  
!        dos(3) = Tr[dos(1:4)*wy]  <- moment density along y
!                    -------- --
!
!        dos(4) = Tr[dos(1:4)*wz]  <- moment density along z
!                    -------- --
!  ===================================================================
   iend = ESV%NumRs; jmax = ESV%jmax
   do ia = 1, ESV%NumSpecies
!     ----------------------------------------------------------------
      p3 => ESV%dos_r_jl(:,:,:,ia)
      call transformDensityMatrix(id,iend,jmax,p3)
!     ----------------------------------------------------------------
      if (rad_derivative) then
!        -------------------------------------------------------------
         p3 => ESV%der_dos_r_jl(:,:,:,ia)
         call transformDensityMatrix(id,iend,jmax,p3)
!        -------------------------------------------------------------
      endif
!     ----------------------------------------------------------------
      p1 => ESV%dos(:,ia)
      call transformDensityMatrix(id,p1)
!     ----------------------------------------------------------------
      p1 => ESV%dos_mt(:,ia)
      call transformDensityMatrix(id,p1)
!     ----------------------------------------------------------------
      p1 => ESV%evalsum(:,ia)
      call transformDensityMatrix(id,p1)
!     ----------------------------------------------------------------
      if (isDensityMatrixNeeded) then
         kmax = ESV%kmax
!        -------------------------------------------------------------
         p3 => ESV%density_matrix(:,:,:,ia)
         call transformDensityMatrix(id,kmax,kmax,p3)
!        -------------------------------------------------------------
      endif
   enddo
!
   end subroutine transformElectroStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function returnMultipleSiteDOS(info,e,dos_array,cfac) result(dos)
!  ===================================================================
   use MathParamModule, only : CONE, ZERO, PI
!
   use MSSolverModule, only : getMSGreenFunction, getMSGreenMatrix
   use MSSolverModule, only : getMSGreenFunctionDerivative
!
   use PotentialTypeModule, only : isASAPotential
!
   use RadialGridModule, only : getGrid
!
   use AtomModule, only : getLocalEvecNew, getLocalNumSpecies
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use SpinRotationModule, only : transformDensityMatrix
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in), optional :: cfac
!
   complex (kind=CmplxKind), intent(inout), target :: dos_array(:)
!
   integer (kind=IntKind), intent(in) :: info(:)
   integer (kind=IntKind) :: is, id, atom, ia
   integer (kind=IntKind) :: ks, jid, j, l, m, jl, kl, klc, p0, n, i
   integer (kind=IntKind) :: jmax, kmax, iend, ns_sqr, gform
!
   real (kind=RealKind) :: r0j(3), evec_new(3), sfac, dos(2)
   real (kind=RealKind), pointer :: r_mesh(:), pra_y(:,:)
!
   complex (kind=CmplxKind) :: dosmt_pola(2)
   complex (kind=CmplxKind) :: dosws_pola(2), ede, energy, cmul
!
   complex (kind=CmplxKind) :: greenint(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: greenint_mt(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: gm(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: gm_mt(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind), pointer :: green_matrix(:,:,:,:)
   complex (kind=CmplxKind), pointer :: green(:,:,:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:), pca_x(:,:), pca_y(:,:)
   complex (kind=CmplxKind), pointer :: der_green(:,:,:,:), der_dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:), p1(:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3)
   energy = adjustEnergy(is,e)
!
   if (present(cfac)) then
      cmul = cfac
   else
      cmul = CONE
   endif
!
!  =================================================================
!  Note: green is the multiple scattering term of the Green function
!        multiplied by r^2 and is expanded in spherical harmonics.
!
!        For nonspin-polarized case, green, greenint, and
!        greenint_mt also contain a factor 2.0.
!
!        green and greenint are in local spin reference frame
!  =================================================================
   if (RelativisticFlag == 2) then
      stop 'Not implemented yet!'
!     --------------------------------------------------------------
  !   green  => getRelGreenFunction(id)
!     --------------------------------------------------------------
   else
!     --------------------------------------------------------------
      green  => getMSGreenFunction(id,gform)
!     --------------------------------------------------------------
      if (rad_derivative) then
!        -----------------------------------------------------------
         der_green => getMSGreenFunctionDerivative(id)
!        -----------------------------------------------------------
      endif
   endif
!
   sfac= TWO/real(n_spin_pola,kind=RealKind)
   ns_sqr = n_spin_cant*n_spin_cant
   iend = size(green,1); kmax = size(green,2); jmax = jofk(kmax)
!
   if (iend /= LastValue(id)%NumRs) then
      call ErrorHandler('returnMultipleSiteDOS','Inconsistent r-mesh size of green and dos_r_jl', &
                        iend,LastValue(id)%NumRs)
   else if (jmax /= LastValue(id)%jmax) then
      call ErrorHandler('returnMultipleSiteDOS','Inconsistent jmax size of green and dos_r_jl', &
                        jmax,LastValue(id)%jmax)
   endif
!
   Grid => getGrid(id)
!
   p0 = 0
   do ia = 1, getLocalNumSpecies(id)
      if (atom < 0 .or. ia == atom) then
         do ks = 1, ns_sqr
!           --------------------------------------------------------
            greenint(ks) = cmul*sfac*getVolumeIntegration( id, iend, Grid%r_mesh, &
                                                        kmax, 2, green(:,:,ks,ia), greenint_mt(ks) )
!           --------------------------------------------------------
            greenint_mt(ks) = cmul*sfac*greenint_mt(ks)
         enddo
!
         if (isASAPotential()) then
            gm_mt = greenint_mt
            if (n_spin_cant == 2) then
!              ------------------------------------------------------
               call transformDensityMatrix(id,gm_mt)
!              ------------------------------------------------------
            endif
            greenint = greenint_mt
            gm = gm_mt
         else
            gm = greenint
            gm_mt = greenint_mt
            if (n_spin_cant == 2) then
!              ------------------------------------------------------
               call transformDensityMatrix(id,gm)
               call transformDensityMatrix(id,gm_mt)
!              ------------------------------------------------------
            endif
         endif
!
!        ===========================================================
!        Note: dosmt, dosws, and dos_array(p0+1:p0+3) are complex type, and
!              their real part are MT-DOS, WS-DOS, and e*DOS, respectively.
!        ===========================================================
         dosmt_pola = CZERO; dosws_pola = CZERO
         if(n_spin_cant == 2) then
!           ---------------------------------------------------------
            evec_new(1:3) = getLocalEvecNew(id)
            dosmt_pola(1) = SQRTm1*HALF*(gm_mt(1)+gm_mt(2)*evec_new(1)+gm_mt(3)*evec_new(2)+gm_mt(4)*evec_new(3))/PI
            dosmt_pola(2) = SQRTm1*HALF*(gm_mt(1)-gm_mt(2)*evec_new(1)-gm_mt(3)*evec_new(2)-gm_mt(4)*evec_new(3))/PI
            dosws_pola(1) = SQRTm1*HALF*(gm(1)+gm(2)*evec_new(1)+gm(3)*evec_new(2)+gm(4)*evec_new(3))/PI
            dosws_pola(2) = SQRTm1*HALF*(gm(1)-gm(2)*evec_new(1)-gm(3)*evec_new(2)-gm(4)*evec_new(3))/PI
         else
            dosmt_pola(1) = SQRTm1*greenint_mt(1)/PI
            dosws_pola(1) = SQRTm1*greenint(1)/PI
         endif
!
         if (node_print_level >= 0) then
            write(6,'(/,''returnMultipleSiteDOS: energy ='',2f18.12,'', id ='',i4,'', ia ='',i4)')energy,id,ia
            if (gform == 0) then
               write(6,'(  ''                       Int [Z*(Tau-t)*Z] on MT         ='',2f18.12)') &
                     dosmt_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Int [Z*(Tau-t)*Z] on VP         ='',2f18.12)') &
                     dosws_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Im{-Int [Z*(Tau-t)*Z]/pi on MT} ='',f18.12)')  &
                     real(dosmt_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
               write(6,'(  ''                       Im{-Int [Z*(Tau-t)*Z]/pi on VP} ='',f18.12)')  &
                     real(dosws_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
            else if (gform == 1) then
               write(6,'(  ''                       Int [Z*Tau*Z-Z*J] on MT         ='',2f18.12)') &
                     dosmt_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Int [Z*Tau*Z-Z*J] on VP         ='',2f18.12)') &
                     dosws_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Im{-Int [Z*Tau*Z-Z*J]/pi on MT} ='',f18.12)')  &
                     real(dosmt_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
               write(6,'(  ''                       Im{-Int [Z*Tau*Z-Z*J]/pi on VP} ='',f18.12)')  &
                     real(dosws_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
            else
               write(6,'(  ''                       Int [Z*Tau*Z] on MT         ='',2f18.12)')     &
                     dosmt_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Int [Z*Tau*Z] on VP         ='',2f18.12)')     &
                     dosws_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Im{-Int [Z*Tau*Z]/pi on MT} ='',f18.12)')      &
                     real(dosmt_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
               write(6,'(  ''                       Im{-Int [Z*Tau*Z]/pi on VP} ='',f18.12)')      &
                     real(dosws_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
            endif
         endif ! print_level
!
         if (atom > 0 .or. ia == 1) then
            dos = real(dosws_pola,kind=RealKind)
         endif
!
         if (isDensityMatrixNeeded) then
!           ---------------------------------------------------------
            green_matrix => getMSGreenMatrix(id)
!           ---------------------------------------------------------
         endif
!
         if (iharris <= 1) then
            ede = energy
         else
            ede = (energy-chempot)
         endif
!
         n = iend*jmax
         do ks = 1, ns_sqr
            p1 => dos_array(p0+1:p0+n)
            dos_r_jl => aliasArray2_c(p1, iend, jmax)
            do jl = 1, jmax
               l = lofj(jl); m = mofj(jl)
               kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
               pcv_x => green(:,kl,ks,ia)
               pcv_y => green(:,klc,ks,ia)
               pcv_z => dos_r_jl(:,jl)
               pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
            enddo
            p0 = p0 + n
            if (rad_derivative) then
               p1 => dos_array(p0+1:p0+n)
               der_dos_r_jl => aliasArray2_c(p1, iend, jmax)
               do jl = 1, jmax
                  l = lofj(jl); m = mofj(jl)
                  kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
                  pcv_x => der_green(:,kl,ks,ia)
                  pcv_y => der_green(:,klc,ks,ia)
                  pcv_z => der_dos_r_jl(:,jl)
                  pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
               enddo
               p0 = p0 + n
            endif
            dos_array(p0+1)  = SQRTm1*greenint(ks)/PI
            dos_array(p0+2) = SQRTm1*greenint_mt(ks)/PI
            dos_array(p0+3) = SQRTm1*ede*greenint(ks)/PI
            dos_array(p0+4) = CZERO
            p0 = p0 + 4
            if (isDensityMatrixNeeded) then
               kmax = size(green_matrix,1)
               pca_x => green_matrix(:,:,ks,ia)
               p1 => dos_array(p0+1:p0+kmax*kmax)
               pca_y => aliasArray2_c(p1,kmax,kmax)
               pca_y = cmul*SQRTm1*pca_x
               p0 = p0 + kmax*kmax
            endif
         enddo
      endif
   enddo
!
!! The following codes will be uncommented and implemented later....................................
!! if (n_spin_cant == 2) then
!!    if (RelativisticFlag == 2) then
!!       stop 'Not implemented yet!'
!!    else   ! The following codes are temporary, and the get-functions need to be implemented................
!!       ESVAL%pSC%torque(1:3) = ESVAL%pSC%torque(1:3) + 0 ! getSpinTorqueMoment(id,beta)
!!       ESVAL%pSC%stoner(1:3) = ESVAL%pSC%stoner(1:3) + 0 ! getStonerParamMoment(id,beta)
!!       ESVAL%pSC%onesite_stoner     = ESVAL%pSC%onesite_stoner + 0 ! getOneSiteStonerParamMoment(id,beta)
!!       ESVAL%pSC%onesite_susab(0:9) = ESVAL%pSC%onesite_susab(0:9) + 0 ! getOneSiteSusceptibilityMoment(id,beta)
!!       ESVAL%pSC%onesite_exchange   =  ESVAL%pSC%onesite_exchange  + 0 ! getOneSiteExchangeParamMoment(id,beta)
!!       do j = 1, ESVAL%pSC%NumJs
!!          ESVAL%pSC%pair_exchab(1:9,j) = ESVAL%pSC%pair_exchab(1:9,j) + 0 ! getPairExchangeParamMoment(id,j,beta,jid,r0j)
!!          ESVAL%pSC%jid(j) = jid
!!          ESVAL%pSC%r0j(1:3,j) = r0j(1:3)
!!       enddo
!!    endif
!! endif
!
   end function returnMultipleSiteDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sumESVAL_r(ESVAL,eGID)
!  ===================================================================
   use GroupCommModule, only : getNumPEsInGroup, GlobalSumInGroup
!
   use RadialGridModule, only : getGrid
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: eGID
   integer (kind=IntKind) :: id, iend, jmax, kmax, j, n, buff_size, ns, ia
!
   complex (kind=CmplxKind), allocatable :: buff(:)
!
   type(ElectroStruct), target, intent(inout) :: ESVAL(LocalNumAtoms)
!
   type(ElectroStruct), pointer :: pESVAL
!
   if (getNumPEsInGroup(eGID) < 2) then
      return
   endif
!
   if (n_spin_cant == 1) then
      buff_size = 0
      do id = 1, LocalNumAtoms
         buff_size = buff_size + 3*n_spin_pola*ESVAL(id)%NumSpecies
      enddo
   else
      buff_size = 0
      do id = 1, LocalNumAtoms
         buff_size = buff_size + 30*ESVAL(id)%NumSpecies
         do ia = 1, ESVAL(id)%NumSpecies
            buff_size = buff_size + 9*ESVAL(id)%pSC(ia)%NumJs
         enddo
      enddo
   endif
   allocate( buff(1:buff_size) )
!
   ns = n_spin_pola*n_spin_cant
!
   n = 0
   do id = 1, LocalNumAtoms
      iend = ESVAL(id)%NumRs
      jmax = ESVAL(id)%jmax
      kmax = ESVAL(id)%kmax
      pESVAL => ESVAL(id)
!     ----------------------------------------------------------------
      call GlobalSumInGroup(eGID,pESVAL%dos_r_jl,iend,jmax,ns,ESVAL(id)%NumSpecies)
!     ----------------------------------------------------------------
      if (rad_derivative) then
!        -------------------------------------------------------------
         call GlobalSumInGroup(eGID,pESVAL%der_dos_r_jl,iend,jmax,ns,ESVAL(id)%NumSpecies)
!        -------------------------------------------------------------
      endif
      if (isDensityMatrixNeeded) then
!        -------------------------------------------------------------
         call GlobalSumInGroup(eGID,pESVAL%density_matrix,kmax,kmax,ns,ESVAL(id)%NumSpecies)
!        -------------------------------------------------------------
      endif
      do ia = 1, ESVAL(id)%NumSpecies
         buff(n+1:n+ns)  = pESVAL%dos(1:ns,ia)
         n = n + ns
         buff(n+1:n+ns)  = pESVAL%dos_mt(1:ns,ia)
         n = n + ns
         buff(n+1:n+ns)  = pESVAL%evalsum(1:ns,ia)
         n = n + ns
         if (n_spin_cant == 2) then
            buff(n+1:n+3)  = pESVAL%pSC(ia)%torque(1:3)
            n = n + 3
            buff(n+1:n+3)  = pESVAL%pSC(ia)%stoner(1:3)
            n = n + 3
            buff(n+1)      = pESVAL%pSC(ia)%onesite_stoner
            n = n + 1
            buff(n+1:n+10) = pESVAL%pSC(ia)%onesite_susab(0:9)
            n = n + 10
            buff(n+1)      = pESVAL%pSC(ia)%onesite_exchange
            n = n + 1
            do j = 1, pESVAL%pSC(ia)%NumJs
               buff(n+1:n+9) = pESVAL%pSC(ia)%pair_exchab(1:9,j)
               n = n + 9
            enddo
         endif
      enddo
   enddo
   if (n /= buff_size) then
!     ----------------------------------------------------------------
      call ErrorHandler('sumESVAL','Inconsistent buffer size',buff_size,n)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call GlobalSumInGroup(eGID,buff(1:n),n)
!  -------------------------------------------------------------------
!
   n = 0
   do id = 1, LocalNumAtoms
      pESVAL => ESVAL(id)
      do ia = 1, ESVAL(id)%NumSpecies
         pESVAL%dos(1:ns,ia)     = buff(n+1:n+ns)
         n = n + ns
         pESVAL%dos_mt(1:ns,ia)  = buff(n+1:n+ns)
         n = n + ns
         pESVAL%evalsum(1:ns,ia) = buff(n+1:n+ns)
         n = n + ns
         if (n_spin_cant == 2) then
            pESVAL%pSC(ia)%torque(1:3)        = buff(n+1:n+3)
            n = n + 3
            pESVAL%pSC(ia)%stoner(1:3)        = buff(n+1:n+3)
            n = n + 3
            pESVAL%pSC(ia)%onesite_stoner     = buff(n+1)
            n = n + 1
            pESVAL%pSC(ia)%onesite_susab(0:9) = buff(n+1:n+10)
            n = n + 10
            pESVAL%pSC(ia)%onesite_exchange   = buff(n+1)
            n = n + 1
            do j = 1, pESVAL%pSC(ia)%NumJs
               pESVAL%pSC(ia)%pair_exchab(1:9,j) = buff(n+1:n+9)
               n = n + 9
            enddo
         endif
      enddo
   enddo
!
   deallocate( buff )
!
   end subroutine sumESVAL_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sumESVAL(is,ESVAL,eGID)
!  ===================================================================
   use GroupCommModule, only : getNumPEsInGroup, GlobalSumInGroup
!
   use RadialGridModule, only : getGrid
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, eGID
   integer (kind=IntKind) :: id, iend, jmax, kmax, j, n, buff_size, ia
!
   complex (kind=CmplxKind), allocatable :: buff(:)
!
   type(ElectroStruct), target, intent(inout) :: ESVAL(LocalNumAtoms)
!
   type(ElectroStruct), pointer :: pESVAL
!
   if (getNumPEsInGroup(eGID) < 2) then
      return
   endif
!
   if (n_spin_cant == 1) then
      buff_size = 0
      do id = 1, LocalNumAtoms
         buff_size = buff_size + 3*ESVAL(id)%NumSpecies
      enddo
   else
      buff_size = 0
      do id = 1, LocalNumAtoms
         buff_size = buff_size + 30*ESVAL(id)%NumSpecies
         do ia = 1, ESVAL(id)%NumSpecies
            buff_size = buff_size + 9*ESVAL(id)%pSC(ia)%NumJs
         enddo
      enddo
   endif
   allocate( buff(1:buff_size) )
!
   n = 0
   do id = 1, LocalNumAtoms
      iend = ESVAL(id)%NumRs
      jmax = ESVAL(id)%jmax
      kmax = ESVAL(id)%kmax
      pESVAL => ESVAL(id)
      if (n_spin_cant == 1) then
         do ia = 1, ESVAL(id)%NumSpecies
!           ----------------------------------------------------------
            call GlobalSumInGroup(eGID,pESVAL%dos_r_jl(:,:,is,ia),iend,jmax)
!           ----------------------------------------------------------
            if (rad_derivative) then
!              -------------------------------------------------------
               call GlobalSumInGroup(eGID,pESVAL%der_dos_r_jl(:,:,is,ia),iend,jmax)
!              -------------------------------------------------------
            endif
            if (isDensityMatrixNeeded) then
!              -------------------------------------------------------
               call GlobalSumInGroup(eGID,pESVAL%density_matrix(:,:,is,ia),kmax,kmax)
!              -------------------------------------------------------
            endif
            buff(n+1) = pESVAL%dos(is,ia)
            buff(n+2) = pESVAL%dos_mt(is,ia)
            buff(n+3) = pESVAL%evalsum(is,ia)
            n = n + 3
         enddo
      else
!        -------------------------------------------------------------
         call GlobalSumInGroup(eGID,pESVAL%dos_r_jl,iend,jmax,4,ESVAL(id)%NumSpecies)
!        -------------------------------------------------------------
         if (rad_derivative) then
!           ----------------------------------------------------------
            call GlobalSumInGroup(eGID,pESVAL%der_dos_r_jl,iend,jmax,4,ESVAL(id)%NumSpecies)
!           ----------------------------------------------------------
         endif
         if (isDensityMatrixNeeded) then
!           ----------------------------------------------------------
            call GlobalSumInGroup(eGID,pESVAL%density_matrix,kmax,kmax,4,ESVAL(id)%NumSpecies)
!           ----------------------------------------------------------
         endif
         do ia = 1, ESVAL(id)%NumSpecies
            buff(n+1:n+4)   = pESVAL%dos(1:4,ia)
            buff(n+5:n+8)   = pESVAL%dos_mt(1:4,ia)
            buff(n+9:n+12)  = pESVAL%evalsum(1:4,ia)
            buff(n+13:n+15) = pESVAL%pSC(ia)%torque(1:3)
            buff(n+16:n+18) = pESVAL%pSC(ia)%stoner(1:3)
            buff(n+19)      = pESVAL%pSC(ia)%onesite_stoner
            buff(n+20:n+29) = pESVAL%pSC(ia)%onesite_susab(0:9)
            buff(n+30)      = pESVAL%pSC(ia)%onesite_exchange
            n = n + 30
            do j = 1, pESVAL%pSC(ia)%NumJs
               buff(n+1:n+9) = pESVAL%pSC(ia)%pair_exchab(1:9,j)
               n = n + 9
            enddo
         enddo
      endif
   enddo
   if (n /= buff_size) then
!     ----------------------------------------------------------------
      call ErrorHandler('sumESVAL','Inconsistent buffer size',buff_size,n)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call GlobalSumInGroup(eGID,buff(1:n),n)
!  -------------------------------------------------------------------
!
   n = 0
   do id = 1, LocalNumAtoms
      pESVAL => ESVAL(id)
      if (n_spin_cant == 1) then
         do ia = 1, ESVAL(id)%NumSpecies
            pESVAL%dos(is,ia) = buff(n+1)
            pESVAL%dos_mt(is,ia) = buff(n+2)
            pESVAL%evalsum(is,ia) = buff(n+3)
            n = n + 3
         enddo
      else
         do ia = 1, ESVAL(id)%NumSpecies
            pESVAL%dos(1:4,ia)                = buff(n+1:n+4)
            pESVAL%dos_mt(1:4,ia)             = buff(n+5:n+8)
            pESVAL%evalsum(1:4,ia)            = buff(n+9:n+12)
            pESVAL%pSC(ia)%torque(1:3)        = buff(n+13:n+15)
            pESVAL%pSC(ia)%stoner(1:3)        = buff(n+16:n+18)
            pESVAL%pSC(ia)%onesite_stoner     = buff(n+19)
            pESVAL%pSC(ia)%onesite_susab(0:9) = buff(n+20:n+29)
            pESVAL%pSC(ia)%onesite_exchange   = buff(n+30)
            n = n + 30
            do j = 1, pESVAL%pSC(ia)%NumJs
               pESVAL%pSC(ia)%pair_exchab(1:9,j) = buff(n+1:n+9)
               n = n + 9
            enddo
         enddo
      endif
   enddo
!
   deallocate( buff )
!
   end subroutine sumESVAL
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine addElectroStruct_r(alp,ESV0,ESV1,is)
!  ===================================================================
!
!  Perform ESV1 = alp*ESV0 + ESV1
!
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: is
   integer (kind=IntKind) :: ia
!
   type(ElectroStruct), intent(in) :: ESV0
   type(ElectroStruct), intent(inout) :: ESV1
!
   real (kind=RealKind), intent(in) :: alp
!
   complex (kind=CmplxKind), pointer :: pca_x(:,:), pca_y(:,:)
!
   if (.not.present(is)) then
      ESV1%dos_r_jl = ESV1%dos_r_jl + alp*ESV0%dos_r_jl
      if (rad_derivative) then
         ESV1%der_dos_r_jl = ESV1%der_dos_r_jl + alp*ESV0%der_dos_r_jl
      endif
      ESV1%dos = ESV1%dos + alp*ESV0%dos
      ESV1%dos_mt = ESV1%dos_mt + alp*ESV0%dos_mt
      ESV1%evalsum = ESV1%evalsum + alp*ESV0%evalsum
      if (isDensityMatrixNeeded) then
         ESV1%density_matrix = ESV1%density_matrix + alp*ESV0%density_matrix
      endif
      if (n_spin_cant == 2) then
         do ia = 1, ESV1%NumSpecies
!           ----------------------------------------------------------
            call addSpinCantStruct(alp,ESV0%pSC(ia),ESV1%pSC(ia))
!           ----------------------------------------------------------
         enddo
      endif
   else
      do ia = 1, ESV1%NumSpecies
         pca_x => ESV0%dos_r_jl(:,:,is,ia)
         pca_y => ESV1%dos_r_jl(:,:,is,ia)
         pca_y = pca_y + alp*pca_x
         if (rad_derivative) then
            pca_x => ESV0%der_dos_r_jl(:,:,is,ia)
            pca_y => ESV1%der_dos_r_jl(:,:,is,ia)
            pca_y = pca_y + alp*pca_x
         endif
         ESV1%dos(is,ia) = ESV1%dos(is,ia) + alp*ESV0%dos(is,ia)
         ESV1%dos_mt(is,ia) = ESV1%dos_mt(is,ia) + alp*ESV0%dos_mt(is,ia)
         ESV1%evalsum(is,ia) = ESV1%evalsum(is,ia) + alp*ESV0%evalsum(is,ia)
         if (isDensityMatrixNeeded) then
            pca_x => ESV0%density_matrix(:,:,is,ia)
            pca_y => ESV1%density_matrix(:,:,is,ia)
            pca_y = pca_y + alp*pca_x
         endif
         if (n_spin_cant == 2 .and. is == 4) then
!           ----------------------------------------------------------
            call addSpinCantStruct(alp,ESV0%pSC(ia),ESV1%pSC(ia))
!           ----------------------------------------------------------
         endif
      enddo
   endif
!
   end subroutine addElectroStruct_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine addElectroStruct_c(alp,ESV0,ESV1,is)
!  ===================================================================
!
!  Perform ESV1 = alp*ESV0 + ESV1
!
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: is
   integer (kind=IntKind) :: ia
!
   type(ElectroStruct), intent(in) :: ESV0
   type(ElectroStruct), intent(inout) :: ESV1
!
   complex (kind=CmplxKind), intent(in) :: alp
!
   complex (kind=CmplxKind), pointer :: pca_x(:,:), pca_y(:,:)
!
   if (.not.present(is)) then
      ESV1%dos_r_jl = ESV1%dos_r_jl + alp*ESV0%dos_r_jl
      if (rad_derivative) then
         ESV1%der_dos_r_jl = ESV1%der_dos_r_jl + alp*ESV0%der_dos_r_jl
      endif
      ESV1%dos = ESV1%dos + alp*ESV0%dos
      ESV1%dos_mt = ESV1%dos_mt + alp*ESV0%dos_mt
      ESV1%evalsum = ESV1%evalsum + alp*ESV0%evalsum
      if (isDensityMatrixNeeded) then
         ESV1%density_matrix = ESV1%density_matrix + alp*ESV0%density_matrix
      endif
      if (n_spin_cant == 2) then
         do ia = 1, ESV1%NumSpecies
!           ----------------------------------------------------------
            call addSpinCantStruct(alp,ESV0%pSC(ia),ESV1%pSC(ia))
!           ----------------------------------------------------------
         enddo
      endif
   else
      do ia = 1, ESV1%NumSpecies
         pca_x => ESV0%dos_r_jl(:,:,is,ia)
         pca_y => ESV1%dos_r_jl(:,:,is,ia)
         pca_y = pca_y + alp*pca_x
         if (rad_derivative) then
            pca_x => ESV0%der_dos_r_jl(:,:,is,ia)
            pca_y => ESV1%der_dos_r_jl(:,:,is,ia)
            pca_y = pca_y + alp*pca_x
         endif
         ESV1%dos(is,ia) = ESV1%dos(is,ia) + alp*ESV0%dos(is,ia)
         ESV1%dos_mt(is,ia) = ESV1%dos_mt(is,ia) + alp*ESV0%dos_mt(is,ia)
         ESV1%evalsum(is,ia) = ESV1%evalsum(is,ia) + alp*ESV0%evalsum(is,ia)
         if (isDensityMatrixNeeded) then
            pca_x => ESV0%density_matrix(:,:,is,ia)
            pca_y => ESV1%density_matrix(:,:,is,ia)
            pca_y = pca_y + alp*pca_x
         endif
         if (n_spin_cant == 2 .and. is == 4) then
!           ----------------------------------------------------------
            call addSpinCantStruct(alp,ESV0%pSC(ia),ESV1%pSC(ia))
!           ----------------------------------------------------------
         endif
      enddo
   endif
!
   end subroutine addElectroStruct_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine addSpinCantStruct_r(alp,SCS0,SCS1)
!  ===================================================================
!
!  Perform SCS1 = alp*SCS0 + SCS1
!
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: alp
!
   type(SpinCantStruct), intent(in) :: SCS0
   type(SpinCantStruct), intent(inout) :: SCS1
!
   integer (kind=IntKind) :: j
!
   SCS1%torque(1:3) = SCS1%torque(1:3) + alp*SCS0%torque(1:3)
   SCS1%stoner(1:3) = SCS1%stoner(1:3) + alp*SCS0%stoner(1:3)
   SCS1%onesite_stoner     = SCS1%onesite_stoner + alp*SCS0%onesite_stoner
   SCS1%onesite_susab(0:9) = SCS1%onesite_susab(0:9) + alp*SCS0%onesite_susab(0:9)
   SCS1%onesite_exchange   =  SCS1%onesite_exchange  + alp*SCS0%onesite_exchange
   do j = 1, SCS1%NumJs
      SCS1%pair_exchab(1:9,j) = SCS1%pair_exchab(1:9,j) + alp*SCS0%pair_exchab(1:9,j)
!     SCS1%jid(j) = SCS1%jid(j) + alp*SCS0%jid(j)
!     SCS1%r0j(1:3,j) = SCS1%r0j(1:3,j) + alp*SCS0%r0j(1:3,j)
   enddo
!
   end subroutine addSpinCantStruct_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine addSpinCantStruct_c(alp,SCS0,SCS1)
!  ===================================================================
!
!  Perform SCS1 = alp*SCS0 + SCS1
!
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: alp
!
   type(SpinCantStruct), intent(in) :: SCS0
   type(SpinCantStruct), intent(inout) :: SCS1
!
   integer (kind=IntKind) :: j
!
   SCS1%torque(1:3) = SCS1%torque(1:3) + alp*SCS0%torque(1:3)
   SCS1%stoner(1:3) = SCS1%stoner(1:3) + alp*SCS0%stoner(1:3)
   SCS1%onesite_stoner     = SCS1%onesite_stoner + alp*SCS0%onesite_stoner
   SCS1%onesite_susab(0:9) = SCS1%onesite_susab(0:9) + alp*SCS0%onesite_susab(0:9)
   SCS1%onesite_exchange   =  SCS1%onesite_exchange  + alp*SCS0%onesite_exchange
   do j = 1, SCS1%NumJs
      SCS1%pair_exchab(1:9,j) = SCS1%pair_exchab(1:9,j) + alp*SCS0%pair_exchab(1:9,j)
!     SCS1%jid(j) = SCS1%jid(j) + alp*SCS0%jid(j)
!     SCS1%r0j(1:3,j) = SCS1%r0j(1:3,j) + alp*SCS0%r0j(1:3,j)
   enddo
!
   end subroutine addSpinCantStruct_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine scaleElectroStruct(alp,ESV)
!  ===================================================================
!
!  Perform ESV = alp*ESV
!
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: alp
!
   type(ElectroStruct), intent(inout) :: ESV
!
   complex (kind=CmplxKind) :: alpc
!
   integer (kind=IntKind) :: ia
!
   alpc = alp
!
   ESV%dos_r_jl = alpc*ESV%dos_r_jl
   if (rad_derivative) then
      ESV%der_dos_r_jl = alpc*ESV%der_dos_r_jl
   endif
   ESV%dos = alpc*ESV%dos
   ESV%dos_mt = alpc*ESV%dos_mt
   ESV%evalsum = alpc*ESV%evalsum
   if (isDensityMatrixNeeded) then
      ESV%density_matrix = alpc*ESV%density_matrix
   endif
!
   if (n_spin_cant == 2) then
      do ia = 1, ESV%NumSpecies
!        -------------------------------------------------------------
         call scaleSpinCantStruct(alp,ESV%pSC(ia))
!        -------------------------------------------------------------
      enddo
   endif
!
   end subroutine scaleElectroStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine scaleSpinCantStruct(alp,SCS)
!  ===================================================================
!
!  Perform SCS = alp*SCS
!
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: alp
!
   type(SpinCantStruct), intent(inout) :: SCS
!
   integer (kind=IntKind) :: j
!
   SCS%torque(1:3) = alp*SCS%torque(1:3)
   SCS%stoner(1:3) = alp*SCS%stoner(1:3)
   SCS%onesite_stoner     = alp*SCS%onesite_stoner
   SCS%onesite_susab(0:9) = alp*SCS%onesite_susab(0:9)
   SCS%onesite_exchange   = alp*SCS%onesite_exchange
   do j = 1, SCS%NumJs
      SCS%pair_exchab(1:9,j) = alp*SCS%pair_exchab(1:9,j)
!     SCS%jid(j) = alp*SCS%jid(j)
!     SCS%r0j(1:3,j) = alp*SCS%r0j(1:3,j)
   enddo
!
   end subroutine scaleSpinCantStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calOnsiteExchParam(SCS)
!  ===================================================================
   implicit none
!
   type(SpinCantStruct), intent(inout) :: SCS
!
!  ===================================================================
!  calculate onesite exchange parameters:
!
!            ab              Ef        ->  ->
!           J  = + 1/pi * Im int Tr [ -p * Tau * d  -
!            i              -inf   L              ab
!
!                                     -p * Tau * p * Tau * d
!                                             0         0   ab
!
!                                     +p * Tau * p * Tau  ]
!                                             a         b
!
!                    ->         ^  00  ->
!           where,   Tau = Tr [ Tau  * Sigma ] / 2  = (Tau , Tau , Tau )
!                            s                            1     2     3
!
!                               ^  00
!                    Tau = Tr [ Tau   ] / 2
!                       0    s
!
!                     ->        ^      ->
!                     p  = Tr [ Pmat * Sigma ] / 2
!                            s
!
!                          ->   ->
!                     p  = p  * e  = pmat_m / 2
!
!
!                    a,b = 1, 2, 3, or x, y, z
!
!                    d   = 1, if a =  b
!                     ab
!                        = 0, if a <> b
!  ===================================================================
   SCS%onesite_exchab(1:9)=SCS%onesite_susab(1:9)
   SCS%onesite_exchab(1)=SCS%onesite_exchab(1)-SCS%onesite_stoner-SCS%onesite_susab(0)
   SCS%onesite_exchab(5)=SCS%onesite_exchab(5)-SCS%onesite_stoner-SCS%onesite_susab(0)
   SCS%onesite_exchab(9)=SCS%onesite_exchab(9)-SCS%onesite_stoner-SCS%onesite_susab(0)
!
   end subroutine calOnsiteExchParam
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine mufind(efermi,efermi_old,int_dos,last_dos,BadFermiEnergy,Lloyd_factor)
!  ===================================================================
   use PhysParamModule, only : Boltzmann, BohrMagneton
!
   use GroupCommModule, only : GlobalSumInGroup
!
   use ScfDataModule, only : isLloyd, getLloydMode
!
   use KreinModule, only : isLloydOn
!
   implicit   none
!
   character (len=6) :: sname='mufind'
!
   integer (kind=IntKind), intent(inout) :: BadFermiEnergy
!
   integer (kind=IntKind) :: id, js
!
   real (kind=RealKind), intent(out) :: efermi
   real (kind=RealKind), intent(in) :: efermi_old
   real (kind=RealKind), intent(in) :: int_dos(n_spin_cant*n_spin_pola,LocalNumAtoms)
   real (kind=RealKind), intent(in) :: last_dos(n_spin_cant*n_spin_pola,LocalNumAtoms)
!
   real (kind=RealKind) :: tnen, xtws
   real (kind=RealKind) :: zvaltss_avr
   real (kind=RealKind) :: efdif
   real (kind=RealKind) :: dosefpa
   real (kind=RealKind) :: ztdif
   real (kind=RealKind) :: N_Green_Contour(2)
   real (kind=RealKind) :: wspace(6)
!
   real (kind=RealKind), intent(out) :: Lloyd_factor(n_spin_pola)
   integer (kind=IntKind) :: Lloyd_nspin
!
!  *******************************************************************
!  calculates the chemical potential and eigenvalue sum
!  *******************************************************************
!
!  ===================================================================
!  calculates ::: integrated n(e) .................................
!  calculate VP valence charges up to chem.pot.............
!  xtws is obtained from a global sum..............................
!  Note: IntegrValue(id)%dos = is the MST part of integrated DOS up to (efermi_old,0) on a 
!                              Gaussian grid contour for atom id
!      ssIntegrValue(id)%dos = is the single site part of integrated DOS up to efermi
!                              on the real energy axis for atom id
!          LastValue(id)%dos = is the MST part of the DOS at (efermi_old,eib) for atom id
!        ssLastValue(id)%dos = is the single site part of the DOS at efermi on the
!                              real energy axis for atom id
!  ===================================================================
!  zvaltss_avr = zvaltss/real(GlobalNumAtoms-NumVacancies,RealKind)
   zvaltss_avr = zvaltss/real(GlobalNumAtoms,RealKind)
   Lloyd_factor=ONE
   wspace = ZERO
   do id = 1, LocalNumAtoms
!     ================================================================
!     Note: wspace(1) here is the difference between the IDOS and the
!           average number of valence electrons per atom
!     ================================================================
      if(n_spin_cant == 2 .or. n_spin_pola == 1) then
         wspace(1) = wspace(1) + int_dos(1,id) - zvaltss_avr
         wspace(2) = wspace(2) + last_dos(1,id)
         wspace(3) = wspace(3) + int_dos(1,id) - zvaltss_avr
         wspace(5) = wspace(5) + last_dos(1,id)
         if ( n_spin_cant==2 ) then
            wspace(4) = wspace(4) + sqrt(int_dos(2,id)**2+int_dos(3,id)**2+int_dos(4,id)**2)
         endif
      else   ! spin polarized case
         wspace(1) = wspace(1) + (int_dos(1,id)+int_dos(2,id)) - zvaltss_avr
         wspace(2) = wspace(2) + (last_dos(1,id)+last_dos(2,id))
         wspace(3) = wspace(3) + int_dos(1,id) - zvaltss_avr*HALF
         wspace(4) = wspace(4) + int_dos(2,id) - zvaltss_avr*HALF
         wspace(5) = wspace(5) + last_dos(1,id)
         wspace(6) = wspace(6) + last_dos(2,id)
      endif
   enddo ! id
!  -------------------------------------------------------------------
   call GlobalSumInGroup(aGID,wspace,6)
!  -------------------------------------------------------------------
   xtws = wspace(1)            ! xtws = the difference between the total 
                               ! integrated DOS of the system up to (efermi_old,0)
                               ! and the total number of valence electrons 
   tnen = wspace(2)            ! tnen = the DOS of the system at (efermi_old,0)
!
!  ===================================================================
!  find the fermi energy
!  tnen is obtained from a global sum
!  ===================================================================
!  if ( iharris < 0 ) then  ! For semiconductors  -- this is a test
!     if (abs(xtws) < 0.001 .or. abs(tnen) < 0.01d0) then
!        efermi = efermi_old
!     else if (xtws/tnen > 0.05d0) then
!        efermi = efermi_old - 0.05d0
!     else if (xtws/tnen < -0.05d0) then
!        efermi = efermi_old + 0.05d0
!     else
!        efermi = efermi_old - xtws/tnen
!     endif
!  else if ( iharris <= 1 ) then
!     if (abs(tnen) > 0.0001d0) then
!        efermi = efermi_old - xtws/tnen
!     else
!        call WarningHandler(sname,'The DOS at the last energy is too small',tnen)
!        efermi = efermi_old
!     endif
!  else
!     efermi = efermi_old
!  endif
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  In code revision 788, the above lines are commented out and are replaced
!  with the following lines, so that in the case when Fermi energy
!  falls into a band gap, where the DOS is very small, the new Fermi
!  energy is adjusted by an amount <= 0.05.
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if (iharris > 1) then
      efermi = efermi_old
   else if (abs(xtws/tnen) > 0.05d0) then
      if (xtws > ZERO) then
         efermi = efermi_old - min(0.01d0,xtws)
      else
         efermi = efermi_old + min(0.01d0,-xtws)
      endif
   else
      efermi = efermi_old - xtws/tnen
   endif
!  ===================================================================
!
   efdif = efermi-efermi_old
!  dosefpa = tnen/real(GlobalNumAtoms-NumVacancies,kind=RealKind)
!  ztdif   = -xtws/real(GlobalNumAtoms-NumVacancies,kind=RealKind)
   dosefpa = tnen/real(GlobalNumAtoms,kind=RealKind)
   ztdif   = -xtws/real(GlobalNumAtoms,kind=RealKind)
!
   if (isLloyd().and.(getLloydMode().eq.1) ) then
      if ( n_spin_pola == 2 ) then
         wspace(3) = wspace(3) + HALF*zvaltss
         wspace(4) = wspace(4) + HALF*zvaltss
      else
         wspace(3) = wspace(3) + zvaltss
      endif
!
      if (n_spin_cant == 1) then
         N_Green_Contour(1) = wspace(3)+wspace(5)*efdif
         N_Green_Contour(2) = wspace(4)+wspace(6)*efdif
      else
         N_Green_Contour(1) = (wspace(3)+wspace(4)+(wspace(5)+wspace(6))*efdif)*HALF
         N_Green_Contour(2) = (N_Green_Contour(1)-(wspace(4)+wspace(6)*efdif))
      endif
!     ----------------------------------------------------------------
      if ( n_spin_pola == 1) then
         Lloyd_nspin = 1
         Lloyd_lastE(1) =two*Lloyd_lastE(1)
         Lloyd_factor(1) = real(-sqrtm1*Lloyd_lastE(1),kind=RealKind)/N_Green_Contour(1)
      else if ( n_spin_cant == 2 ) then
         Lloyd_nspin=2
         Lloyd_factor(1) = real(-sqrtm1*Lloyd_lastE(1),kind=RealKind)/N_Green_Contour(1)
         Lloyd_factor(2) = real(-sqrtm1*Lloyd_lastE(2),kind=RealKind)/N_Green_Contour(2)
      else   ! spin polarized case
         Lloyd_nspin=2
         Lloyd_factor(1) = real(-sqrtm1*Lloyd_lastE(1),kind=RealKind)/N_Green_Contour(1)
         Lloyd_factor(2) = real(-sqrtm1*Lloyd_lastE(2),kind=RealKind)/N_Green_Contour(2)
      endif ! spin
   endif ! Lloyd on and mode=1
!  ----------------------------------------------------------------
!
!  ===================================================================
!  Major printout....... Integrated n(e) at end of contour.........
!  ===================================================================
   if ( node_print_level >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,29x,a)')'**********************'
      write(6,'( 29x,a )')'* Output from mufind *'
      write(6,'(29x,a,/)')'**********************'
      write(6,'(/,''At Old Fermi Energy'')')
      write(6,'(4x,''GLOBAL VP Int[n(e)] per Atom'',t40,''='',1f18.12)') &
               (zvaltss+xtws)/real(GlobalNumAtoms-NumVacancies,kind=RealKind)
      write(6,'(4x,''GLOBAL (Val-Int[n(e)]) per Site'',t40,''='',1f18.12)') &
                ztdif
      write(6,'(4x,''Fermi Energy DOS per Site'',t40,''='',1f18.12)')dosefpa
!
      write(6,'(4x,a,t40,a,1f18.12)')'Specific Heat Coefficient per Site', &
                '=',PI*PI*THIRD*Boltzmann*Boltzmann*dosefpa
      write(6,'(4x,a,t40,a,1f18.12)')'Magnetic Susceptibility per Site',   &
                '=',BohrMagneton*BohrMagneton*dosefpa
!
      write(6,'(4x,a,t40,a,1f18.12)')                                 &
           'Num Valence Electrons per Atom','=',                      &
            zvaltss/real(GlobalNumAtoms-NumVacancies,kind=RealKind)
!
      if (isLloyd().and.(getLloydMode().eq.1) ) then
         write (6,'(4x,a,t40,a,2f18.12)') 'Lloyd values per Atom','=', &
              & aimag(Lloyd_lastE(1:Lloyd_nspin))/real(GlobalNumAtoms-NumVacancies,kind=RealKind)
         write (6,'(4x,a,t40,a,2f18.12)') 'N Green per Atom','=', &
              & N_Green_Contour(1:Lloyd_nspin)/real(GlobalNumAtoms-NumVacancies,kind=RealKind)
         write (6,'(4x,a,t40,a,3f18.12)') 'Lloyd factor','=', Lloyd_factor(1:Lloyd_nspin)
      end if ! isLloyd, mode=1
   endif
!
   if ( isIterateEfOn ) then
      ztdif   = -xtws/real(GlobalNumAtoms,kind=RealKind)
      if (abs(efdif) > 0.10d0 .or. abs(ztdif)> 0.010 ) then
!        if ( efdif > 0.1d0 ) then
!           efdif = 0.1d0   ! This is dangerous!
!        else if (efdif < -0.1d0) then
!           efdif = -0.10d0   ! This is dangerous!
!        endif
!        if ( ztdif*efdif < 0.0 ) then
!           efdif = -efdif
!        endif
!        efermi = efermi_old + efdif
         if ( isLloydOn() ) then
            BadFermiEnergy=0
         else
            BadFermiEnergy = BadFermiEnergy + 1
         endif
      else
         BadFermiEnergy = 0
      endif
   endif
!
   if ( node_print_level >= 0) then
      write(6,'(/,''New Fermi energy'',t40,''='',f18.12)') efermi
      write(6,'(''Change in Fermi energy'',t40,''='',f18.12)') efermi-efermi_old
   endif
!
   if (stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine mufind
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkDensity(id,ia,CurrentValue,message)
!  ===================================================================
   use IntegrationModule, only : calIntegration
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use RadialGridModule, only : getGrid
!
   implicit none
!
   character (len=*), intent(in) :: message
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: js, ir
!
   real (kind=RealKind) :: rho_0(jend_max), intr(jend_max)
   real (kind=RealKind) :: rfac, q_VP, q_MT
!
   type (ElectroStruct), intent(in) :: CurrentValue
   type (GridStruct), pointer :: Grid
!
   if ( node_print_level >= 0) then
      write(6,'(/,2x,60(''-''))')
      write(6,'(2x,a,t40,a,i4,a,i4)')'Local Atom Index','=',id,', Atomic Species = ',ia
      write(6,'(a)')message
   endif
!
   Grid => getGrid(id)
   rfac = ONE/(sqrt(PI4))
   if (n_spin_pola == 1 .or. n_spin_cant == 2) then
      do ir = 1, CurrentValue%NumRs
!        rho_0(ir) = real(CurrentValue%dos_r_jl(ir,1,1,ia),kind=RealKind)*rfac &
!                    /(Grid%r_mesh(ir)*Grid%r_mesh(ir))
         rho_0(ir) = real(CurrentValue%dos_r_jl(ir,1,1,ia),kind=RealKind)*rfac
      enddo
   else
      do ir = 1, CurrentValue%NumRs
         rho_0(ir) = real(CurrentValue%dos_r_jl(ir,1,1,ia)+CurrentValue%dos_r_jl(ir,1,2,ia),kind=RealKind)*rfac
!        rho_0(ir) = real(CurrentValue%dos_r_jl(ir,1,1)+CurrentValue%dos_r_jl(ir,1,2),kind=RealKind)*rfac &
!                    /(Grid%r_mesh(ir)*Grid%r_mesh(ir))
      enddo
   endif
!  -------------------------------------------------------------------
   call calIntegration(Grid%jmt,Grid%r_mesh,rho_0,intr,0)
!  call calIntegration(Grid%jmt,Grid%r_mesh,rho_0,intr,2)
   do js=1,n_spin_pola*n_spin_cant
!     ----------------------------------------------------------------
      q_VP=getVolumeIntegration( id, CurrentValue%NumRs, Grid%r_mesh,    &
                                 CurrentValue%kmax, CurrentValue%jmax, 2, CurrentValue%dos_r_jl(:,:,js,ia), q_MT )
!     ----------------------------------------------------------------
      if ( node_print_level >= 0) then
         write(6,'(4x,a,t40,a,f18.11)')'Integrated spherical Density in MT','=', intr(Grid%jmt)*PI4
         write(6,'(4x,a,t40,a,2f18.11)')'Integrated Electron Density in VP, MT','=',q_VP, q_MT
         write(6,'(4x,a,t40,a,3f18.11)')'Integrated DOS in MT','=',real(CurrentValue%dos_mt(js,ia),RealKind)
         write(6,'(4x,a,t40,a,3f18.11)')'Integrated DOS in VP','=',real(CurrentValue%dos(js,ia),RealKind)
      endif
   enddo
!
   end subroutine checkDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateValenceDOS(efermi,Lloyd_factor)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use PotentialModule, only : getV0, getVdif
!
   use SystemVolumeModule, only : getSystemVolume
!
   use ScfDataModule, only : isExchangeParamNeeded
   use ScfDataModule, only : isLloyd, getLloydMode
!
   use AtomModule, only : getLocalEvecNew
!
   implicit   none
!
   character (len=17) :: sname='updateValenceDOS'
!
   integer (kind=IntKind) :: id, is, js, ia
!
   real (kind=RealKind), intent(in):: efermi
   real (kind=RealKind) :: vtot, fef
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), intent(in) :: Lloyd_factor(n_spin_pola)
!
   if ( node_print_level >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,25x,a)')'********************************'
      write(6,'( 25x,a )')'* Output from updateValenceDOS *'
      write(6,'(25x,a,/)')'********************************'
!
      write(6,'(/,''New Fermi energy'',t40,''='',f18.12)') efermi
      write(6,'(''Change in Fermi energy'',t40,''='',f18.12)') efermi-chempot
!
      vtot = getSystemVolume()
      fef = PI*PI*(THREE*zvaltss/(PI*vtot))**(TWO*THIRD)-getV0(1) ! getV0 usually returns 0 since after the potential
                                                                  ! has been shifted accordingly, the v0 is set to be zero.
      write(6,'(/,a)')'For free electron model:'
      write(6,'(4x,a,t40,a,f18.12)')'Fermi energy','=',fef
      write(6,'(4x,a,t40,a,f18.12,/)')'DOS per atom at the Fermi energy','=', &
               vtot*sqrt(fef)/(TWO*PI*(GlobalNumAtoms-NumVacancies))
   endif
!
!  ===================================================================
!  update integrated densities of states and Green function up to fermi
!  energy
!  ===================================================================
   do id = 1, LocalNumAtoms
      do is = 1, n_spin_pola/n_spin_cant
         if (isLloyd().and.(getLloydMode().eq.1) ) then
!           ----------------------------------------------------------
            call scaleElectroStruct(Lloyd_factor(is),IntegrValue(id))
!           ----------------------------------------------------------
         endif
      enddo
!
      do js=1,n_spin_pola*n_spin_cant
!        -------------------------------------------------------------
!        call checkDOSsign()   ! Will build this code later
!        -------------------------------------------------------------
         if (print_level(id) >= 0) then
            do ia = 1, IntegrValue(id)%NumSpecies
               write(6,'(3(a,i3),t30,a,d15.8)')'id =',id ,', js =',js,', ia =',ia, &
                  'IDOS in M.T. =', real(IntegrValue(id)%dos_mt(js,ia),kind=RealKind)
               write(6,'(t30,a,d15.8)')                                  &
                  'IDOS in Cell =', real(IntegrValue(id)%dos(js,ia),kind=RealKind)
               write(6,'(t30,a,d15.8)')                                  &
                  'Band Energy  =', real(IntegrValue(id)%evalsum(js,ia),kind=RealKind)
            enddo
         endif
      enddo
      if (n_spin_cant == 2) then
         do ia = 1, IntegrValue(id)%NumSpecies
!           ----------------------------------------------------------
            call calOnsiteExchParam(IntegrValue(id)%pSC(ia))
!           ----------------------------------------------------------
         enddo
      endif
   enddo
!
   if(stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine updateValenceDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageElectroStruct(eValue)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   integer (kind=IntKind) :: id, n, ia
!
   type (ElectroStruct), intent(inout) :: eValue(LocalNumAtoms)
!
   complex (kind=CmplxKind) :: cfac
!
   cfac = CONE/real(NumPEsInEGroup,kind=RealKind)
!
   n = 0
   do id = 1, LocalNumAtoms
      do ia = 1, eValue(id)%NumSpecies
         wk_dos(n+1:n+4) = eValue(id)%dos(1:4,ia)
         wk_dos(n+5:n+8) = eValue(id)%dos_mt(1:4,ia)
         wk_dos(n+9:n+12) = eValue(id)%evalsum(1:4,ia)
         n = n + 12
!        -------------------------------------------------------------
         call zcopy(eValue(id)%size,eValue(id)%dos_r_jl(1,1,1,ia),1,wk_dos(n+1),1)
!        -------------------------------------------------------------
         n = n + eValue(id)%size
         if (rad_derivative) then
!           ----------------------------------------------------------
            call zcopy(eValue(id)%size,eValue(id)%der_dos_r_jl(1,1,1,ia),1,wk_dos(n+1),1)
!           ----------------------------------------------------------
            n = n + eValue(id)%size
         endif
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call GlobalSumInGroup(eGID,wk_dos,n)
!  -------------------------------------------------------------------
   wk_dos = cfac*wk_dos
!
   n = 0
   do id = 1, LocalNumAtoms
      do ia = 1, eValue(id)%NumSpecies
         eValue(id)%dos(1:4,ia) = wk_dos(n+1:n+4)
         eValue(id)%dos_mt(1:4,ia) = wk_dos(n+5:n+8)
         eValue(id)%evalsum(1:4,ia) = wk_dos(n+9:n+12)
         n = n + 12
!        -------------------------------------------------------------
         call zcopy(eValue(id)%size,wk_dos(n+1),1,eValue(id)%dos_r_jl(1,1,1,ia),1)
!        -------------------------------------------------------------
         n = n + eValue(id)%size
         if (rad_derivative) then
!           ----------------------------------------------------------
            call zcopy(eValue(id)%size,wk_dos(n+1),1,eValue(id)%der_dos_r_jl(1,1,1,ia),1)
!           ----------------------------------------------------------
            n = n + eValue(id)%size
         endif
      enddo
   enddo
!
   end subroutine averageElectroStruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calDensity(efermi) ! This code needs to be modified for random alloys case
!  ===================================================================
   use ValenceDensityModule, only : updateValenceEvec,                &
                                    updateValenceCharge,              &
                                    updateValenceEnergy,              &
                                    updateValenceDensity, updateFermiEnergy
   use ValenceDensityModule, only : printValenceDensity
!
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   character (len=17) :: sname='calDensity'
!
   integer (kind=IntKind) :: id, ir, ia
!
   real (kind=RealKind), intent(in):: efermi
   real (kind=RealKind) :: dos(4)
   real (kind=RealKind) :: dos_mt(4)
   real (kind=RealKind) :: evalsum(4)
   real (kind=RealKind) :: torque(3)
!
   do id = 1, LocalNumAtoms
      torque = ZERO
      do ia = 1, IntegrValue(id)%NumSpecies
         dos = real(IntegrValue(id)%dos(:,ia),kind=RealKind)
         dos_mt = real(IntegrValue(id)%dos_mt(:,ia),kind=RealKind)
         evalsum = real(IntegrValue(id)%evalsum(:,ia),kind=RealKind)
!        -------------------------------------------------------------
         call updateValenceCharge(id,ia,dos,dos_mt)
         call updateValenceEnergy(id,ia,evalsum,exc(ia,id))
         call updateValenceDensity(id,ia,IntegrValue(id)%dos_r_jl(:,:,:,ia))
!        -------------------------------------------------------------
         if (rad_derivative) then
!           ----------------------------------------------------------
            call updateValenceDensity(id,ia,IntegrValue(id)%der_dos_r_jl(:,:,:,ia),   &
                                      isDerivative=.true.)
!           ----------------------------------------------------------
         endif
         if (n_spin_cant == 2 .and. RelativisticFlag .ne. 2) then
            torque = torque + getLocalSpeciesContent(id,ia)*IntegrValue(id)%pSC(ia)%torque
         endif
      enddo
      if (n_spin_cant == 2 .and. RelativisticFlag .ne. 2) then
         call updateValenceEvec(id,torque)
      endif
   enddo
!
!  ===================================================================
!  Finally, update the Fermi energy and other global quantities...
!  -------------------------------------------------------------------
   call updateFermiEnergy(efermi,zvaltss)
!  -------------------------------------------------------------------
   if (node_print_level >= 0) then
!     ----------------------------------------------------------------
      call printValenceDensity()
!     ----------------------------------------------------------------
   endif
!
   if(stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine calDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printExchangeParam(fname18,fname19)
!  ===================================================================
   use MPPModule, only : AnyPE
   use MPPModule, only : sendPackage, recvPackage, packMessage, unpackMessage
!
   use GroupCommModule, only : GlobalMaxInGroup, isInSamegroup
!
   use SystemModule, only : getAtomPosition
!
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   character (len=*), intent(in) :: fname18
   character (len=*), intent(in) :: fname19
!
   integer (kind=IntKind) :: id, ia, ja, ic, nla, MaxJs, MaxNla, j, itmp(2), njs_l
   integer (kind=IntKind), allocatable :: jid_g(:,:), njs_g(:)
!
   real (kind=RealKind), allocatable :: onesite_exchab(:,:)
   real (kind=RealKind), allocatable :: pair_exchab(:,:,:)
   real (kind=RealKind), allocatable :: r0j_g(:,:,:)
   real (kind=RealKind) :: dis, cfac
!
   itmp(1) = 0
   do id = 1, LocalNumAtoms
      do ia = 1, IntegrValue(id)%NumSpecies
         itmp(1) = max(itmp(1), IntegrValue(id)%pSC(ia)%NumJs)
      enddo
   enddo
   itmp(2) = LocalNumAtoms
   call GlobalMaxInGroup(aGID,itmp,2)
   MaxJs = itmp(1)
   MaxNla = itmp(2)
!
   if (MyPE == 0) then
      allocate( onesite_exchab(1:9,1:GlobalNumAtoms) )
      onesite_exchab(1:9,1:GlobalNumAtoms) = ZERO
!
      open(unit=18,file=trim(fname18),status='unknown',form='formatted')
!
      write(18,'(112(''=''))')
      write(18,'(a,a,a)')'  Index                J(1:3,1)',           &
                         '                           J(1:3,2)',       &
                         '                           J(1:3,3)'
      write(18,'(112(''-''))')
!
      if (MaxJs > 0) then
         allocate( pair_exchab(1:9,1:MaxJs,1:GlobalNumAtoms) )
         pair_exchab(1:9,1:MaxJs,1:GlobalNumAtoms) = ZERO
         allocate( jid_g(1:MaxJs,1:GlobalNumAtoms) )
         allocate( r0j_g(1:3,1:MaxJs,1:GlobalNumAtoms) )
         allocate( njs_g(1:GlobalNumAtoms) )
         jid_g(:,:) = -1
         njs_g(:) = -1
!
         open(unit=19,file=trim(fname19),status='unknown',form='formatted')
         write(19,'(134(''=''))')
         write(19,'(a,a,a)')'  Index   j     id              J_0j(1:3,1)', &
                            '                        J_0j(1:3,2)',    &
                       '                        J_0j(1:3,3)               r'
         write(19,'(134(''-''))')
      endif
!
      do ic = 2, NumPEsInAGroup
!        -------------------------------------------------------------
         call recvPackage(110110,AnyPE)
         call unpackMessage(nla)
!        -------------------------------------------------------------
         do id = 1, nla
!           ----------------------------------------------------------
            call unpackMessage(ja)
!           ----------------------------------------------------------
            if (ja < 1 .or. ja > GlobalNumAtoms) then
!              -------------------------------------------------------
               call ErrorHandler('printExchangeParam','Invalid ja',ja)
!              -------------------------------------------------------
            endif
!           ----------------------------------------------------------
            call unpackMessage(njs_g(ja))
!           ----------------------------------------------------------
            if (njs_g(ja) < 0 .or. njs_g(ja) > MaxJs) then
!              -------------------------------------------------------
               call ErrorHandler('printExchangeParam','Invalid njs',  &
                                 njs_g(ja))
!              -------------------------------------------------------
            else if (njs_g(ja) > 0) then
!              -------------------------------------------------------
               call unpackMessage(jid_g(1:njs_g(ja),ja),njs_g(ja))
               call unpackMessage(r0j_g(1:3,1:njs_g(ja),ja),3,njs_g(ja))
               call unpackMessage(pair_exchab(1:9,1:njs_g(ja),ja),9,njs_g(ja))
!              -------------------------------------------------------
            endif
!           ----------------------------------------------------------
            call unpackMessage(onesite_exchab(1:9,ja),9)
!           ----------------------------------------------------------
         enddo
      enddo
      do id = 1, LocalNumAtoms
         ja = AtomIndex(id)
         onesite_exchab(1:9,ja) = ZERO
         njs_g(ja) = IntegrValue(id)%pSC(1)%NumJs
         do ia = 1, IntegrValue(id)%NumSpecies
            onesite_exchab(1:9,ja) = onesite_exchab(1:9,ja) +         &
                   getLocalSpeciesContent(id,ia)*IntegrValue(id)%pSC(ia)%onesite_exchab(1:9)
         enddo
         pair_exchab(:,:,ja) = ZERO
         do j = 1, njs_g(ja)
            jid_g(j,ja) = IntegrValue(id)%pSC(1)%jid(j)
            r0j_g(1:3,j,ja) = IntegrValue(id)%pSC(1)%r0j(1:3,j)
            do ia = 1, IntegrValue(id)%NumSpecies
               pair_exchab(1:9,j,ja) = pair_exchab(1:9,j,ja) +  &
                   getLocalSpeciesContent(id,ia)*IntegrValue(id)%pSC(ia)%pair_exchab(1:9,j)
            enddo
         enddo
      enddo
!
      do ja = 1, GlobalNumAtoms
         write(18,'(1i7,3(2x,3f11.5))')ja,onesite_exchab(1:9,ja)
      enddo
!
      close(unit=18)
!
      deallocate( onesite_exchab )
!
      if (MaxJs > 0) then
         do ja = 1, GlobalNumAtoms
            j = 1
            dis = sqrt( r0j_g(1,j,ja)*r0j_g(1,j,ja) +                &
                        r0j_g(2,j,ja)*r0j_g(2,j,ja) +                &
                        r0j_g(3,j,ja)*r0j_g(3,j,ja) )
            write(19,'(1i7,i4,i7,3(2x,3f11.5),1x,f10.5)')ja,j,jid_g(j,ja),  &
                                              pair_exchab(1:9,j,ja),dis
            do j = 2, njs_g(ja)
               dis = sqrt( r0j_g(1,j,ja)*r0j_g(1,j,ja) +             &
                           r0j_g(2,j,ja)*r0j_g(2,j,ja) +             &
                           r0j_g(3,j,ja)*r0j_g(3,j,ja) )
               write(19,'(7x,i4,i7,3(2x,3f11.5),1x,f10.5)')j,jid_g(j,ja),   &
                                              pair_exchab(1:9,j,ja),dis
            enddo
         enddo
!
         close(unit=19)
!
         deallocate( pair_exchab, jid_g, njs_g, r0j_g )
      endif
   else if (isInSameGroup(aGID,0)) then
      allocate( onesite_exchab(1:9,LocalNumAtoms) )
      allocate( pair_exchab(1:9,1:MaxJs,LocalNumAtoms) )
      onesite_exchab = ZERO; pair_exchab = ZERO
      do id = 1, LocalNumAtoms
         njs_l = IntegrValue(id)%pSC(1)%NumJs
!        =============================================================
!        Average pair_exchab and onsite_exchab over species.
!        =============================================================
         do ia = 1, IntegrValue(id)%NumSpecies
            cfac = getLocalSpeciesContent(id,ia)
            if (njs_l > 0) then
               pair_exchab(1:9,1:njs_l,id) = pair_exchab(1:9,1:njs_l,id) + &
                           cfac*IntegrValue(id)%pSC(ia)%pair_exchab(1:9,1:njs_l)
            endif
            onesite_exchab(1:9,id) = onesite_exchab(1:9,id) +           &
                              cfac*IntegrValue(id)%pSC(ia)%onesite_exchab(1:9)
         enddo
      enddo
!     ----------------------------------------------------------------
      call packMessage(LocalNumAtoms)
!     ----------------------------------------------------------------
      do id = 1, LocalNumAtoms
         njs_l = IntegrValue(id)%pSC(1)%NumJs
!        -------------------------------------------------------------
         call packMessage(AtomIndex(id))
         call packMessage(njs_l)
!        -------------------------------------------------------------
         if (njs_l > 0) then
!           ----------------------------------------------------------
            call packMessage(IntegrValue(id)%pSC(1)%jid(1:njs_l),njs_l)
            call packMessage(IntegrValue(id)%pSC(1)%r0j(1:3,1:njs_l),3,njs_l)
            call packMessage(pair_exchab(1:9,1:njs_l,id),9,njs_l)
!           ----------------------------------------------------------
         endif
!        -------------------------------------------------------------
         call packMessage(onesite_exchab(1:9,id),9)
!        -------------------------------------------------------------
      enddo
!     ----------------------------------------------------------------
      call sendPackage(110110,0)
!     ----------------------------------------------------------------
   endif
!
   call syncAllPEs()
!
   end subroutine printExchangeParam
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSS_dosdata(is,ie,e_real,CurrentValue,redundant)
!  ===================================================================
   use AtomModule, only : getLocalEvecNew
!
   implicit none
!
   logical, intent(in) :: redundant
!
   integer (kind=IntKind), intent(in) :: is, ie
   integer (kind=IntKind) :: id, ks, ns_sqr, ia
!
   real (kind=RealKind), intent(in) :: e_real
   real (kind=RealKind) :: evec_new(3)
!
   type (ElectroStruct), intent(in) :: CurrentValue(LocalNumAtoms)
!
   ns_sqr = n_spin_cant*n_spin_cant
!
   do id = 1, LocalNumAtoms
      do ia = 1, CurrentValue(id)%NumSpecies
         if (.not.redundant .or. MyPEinEGroup == 0) then
            SS_dosdata(id)%rarray3(1,ie,ia) = e_real
            if (n_spin_cant == 2) then
               do ks = 1, ns_sqr
                  SS_dosdata(id)%rarray3(ks+1,ie,ia) = real(CurrentValue(id)%dos(ks,ia),kind=RealKind)
               enddo
               evec_new(1:3) = getLocalEvecNew(id)
!              =======================================================
!              dos(1) = electron density of states
!              dos(2:4) = moment density of states
!              =======================================================
               SS_dosdata(id)%rarray3(6,ie,ia) = HALF*( real(CurrentValue(id)%dos(1,ia),RealKind)      &
                                                +real(CurrentValue(id)%dos(2,ia),RealKind)*evec_new(1) &
                                                +real(CurrentValue(id)%dos(3,ia),RealKind)*evec_new(2) &
                                                +real(CurrentValue(id)%dos(4,ia),RealKind)*evec_new(3) )
               SS_dosdata(id)%rarray3(7,ie,ia) = HALF*( real(CurrentValue(id)%dos(1,ia),RealKind)      &
                                                -real(CurrentValue(id)%dos(2,ia),RealKind)*evec_new(1) &
                                                -real(CurrentValue(id)%dos(3,ia),RealKind)*evec_new(2) &
                                                -real(CurrentValue(id)%dos(4,ia),RealKind)*evec_new(3) )
            else
               SS_dosdata(id)%rarray3(is+1,ie,ia) = real(CurrentValue(id)%dos(is,ia),kind=RealKind)
            endif
         else
            if (n_spin_cant == 2) then
               SS_dosdata(id)%rarray3(1:7,ie,ia) = ZERO
            else
               SS_dosdata(id)%rarray3(1,ie,ia) = ZERO
               SS_dosdata(id)%rarray3(is+1,ie,ia) = ZERO
            endif
         endif
      enddo
   enddo
!
   end subroutine calSS_dosdata
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calDOSDATA(is,ie,e_real,CurrentValue,redundant)
!  ===================================================================
   use AtomModule, only : getLocalEvecNew
!
   implicit none
!
   logical, intent(in) :: redundant
!
   integer (kind=IntKind), intent(in) :: is, ie
   integer (kind=IntKind) :: id, ks, ns_sqr, ia
!
   real (kind=RealKind), intent(in) :: e_real
   real (kind=RealKind) :: evec_new(3)
!
   type (ElectroStruct), intent(in) :: CurrentValue(LocalNumAtoms)
!
   ns_sqr = n_spin_cant*n_spin_cant
!
   do id = 1, LocalNumAtoms
      do ia = 1, CurrentValue(id)%NumSpecies
         if (.not.redundant .or. MyPEinEGroup == 0) then
            dosdata(id)%rarray3(1,ie,ia) = e_real
            if (n_spin_cant == 2) then
               do ks = 1, ns_sqr
                  dosdata(id)%rarray3(ks+1,ie,ia) = real(CurrentValue(id)%dos(ks,ia),kind=RealKind)
               enddo
               evec_new(1:3) = getLocalEvecNew(id)
!              =======================================================
!              dos(1) = electron density of states
!              dos(2:4) = moment density of states
!              =======================================================
               dosdata(id)%rarray3(6,ie,ia) = HALF*( real(CurrentValue(id)%dos(1,ia),RealKind)      &
                                             +real(CurrentValue(id)%dos(2,ia),RealKind)*evec_new(1) &
                                             +real(CurrentValue(id)%dos(3,ia),RealKind)*evec_new(2) &
                                             +real(CurrentValue(id)%dos(4,ia),RealKind)*evec_new(3) )
               dosdata(id)%rarray3(7,ie,ia) = HALF*( real(CurrentValue(id)%dos(1,ia),RealKind)      &
                                             -real(CurrentValue(id)%dos(2,ia),RealKind)*evec_new(1) &
                                             -real(CurrentValue(id)%dos(3,ia),RealKind)*evec_new(2) &
                                             -real(CurrentValue(id)%dos(4,ia),RealKind)*evec_new(3) )
            else
               dosdata(id)%rarray3(is+1,ie,ia) = real(CurrentValue(id)%dos(is,ia),kind=RealKind)
            endif
         else
            if (n_spin_cant == 2) then
               dosdata(id)%rarray3(1:7,ie,ia) = ZERO
            else
               dosdata(id)%rarray3(1,ie,ia) = ZERO
               dosdata(id)%rarray3(is+1,ie,ia) = ZERO
            endif
         endif
      enddo
   enddo
!
   end subroutine calDOSDATA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine correctDOSDATA(js,id,ne,dimdos,e_imag)
!  ===================================================================
   use AtomModule, only : getLocalEvecNew, getLocalNumSpecies
   implicit none
!
   integer (kind=IntKind), intent(in) :: js, id, ne
   integer (kind=IntKind) :: ie, ia
!
   real (kind=RealKind), intent(in) :: dimdos(ne), e_imag
   real (kind=RealKind) :: evec_new(3)
!
   do ia = 1, getLocalNumSpecies(id)
      do ie = 1, ne
!        write(6,'(a,2d15.8)') 'dos, dimdos(ie)*e = ',dosdata(id)%rarray3(js+1,ie,ia), dimdos(ie)*e_imag
         dosdata(id)%rarray3(js+1,ie,ia) = dosdata(id)%rarray3(js+1,ie,ia) + dimdos(ie)*e_imag
      enddo
!
      if (n_spin_cant == 2 .and. js == 4) then
         evec_new(1:3) = getLocalEvecNew(id)
         dosdata(id)%rarray3(6,ie,ia) = HALF*( dosdata(id)%rarray3(2,ie,ia) &
                                             + dosdata(id)%rarray3(3,ie,ia)*evec_new(1) &
                                             + dosdata(id)%rarray3(4,ie,ia)*evec_new(2) &
                                             + dosdata(id)%rarray3(5,ie,ia)*evec_new(3) )
         dosdata(id)%rarray3(7,ie,ia) = HALF*( dosdata(id)%rarray3(2,ie,ia) &
                                             - dosdata(id)%rarray3(3,ie,ia)*evec_new(1) &
                                             - dosdata(id)%rarray3(4,ie,ia)*evec_new(2) &
                                             - dosdata(id)%rarray3(5,ie,ia)*evec_new(3) )
      endif
   enddo
!
   end subroutine correctDOSDATA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeSS_DOS(cpath,systemid,ig)
!  ===================================================================
   use MPPModule, only : getMyPE
!
   use PublicParamDefinitionsModule, only : PrintEachAtomDOS, PrintDOSswitchOff
!
   use PublicTypeDefinitionsModule, only : NeighborStruct
!
   use GroupCommModule, only : getGroupID, getMyPEinGroup
   use GroupCommModule, only : GlobalSumInGroup
!
   use ProcMappingModule, only : isAtomOnMyProc
!
   use Atom2ProcModule, only : getAtom2ProcInGroup, getLocalIndex
!
   use SystemModule, only : getNumAtomTypes, getAtomType, getNumAtomsOfType, &
                            getAtomTypeName
!
   use NeighborModule, only : getNeighbor
!
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   character (len=*), intent(in) :: cpath
   character (len=*), intent(in) :: systemid
   character (len=150) :: fname
   character (len=60) :: dos_title(3)
   character (len=15) :: anm
!
   integer (kind=IntKind), intent(in) :: ig ! = -1: print DOS for each atom
                                            ! =  0: will not print DOS
                                            ! >  0: print the DOS of the atom with global index ig
!
   integer (kind=IntKind) :: nvals
   integer (kind=IntKind) :: DOSpid
   integer (kind=IntKind) :: s, ie, id, iv, n, gid, kGID, MyPEinKGroup
   integer (kind=IntKind) :: nume_dos, it, ntypes, ia
   integer (kind=IntKind), parameter :: fu = 250, maxatoms = 20
!
   real (kind=RealKind), allocatable, target :: aver_dos(:,:,:)
   real (kind=RealKind), pointer :: p_dos(:,:), p_tdos(:,:)
   real (kind=RealKind) :: rfac, cfac
!
   type (NeighborStruct), pointer :: Neighbor
!
   if (.not.allocated(SS_dosdata)) then
      call ErrorHandler('writeSS_DOS', 'Need to call calValenceDOS first')
   endif
!
!
!  nume_dos = size(SS_dosdata,dim=2)
!  nvals = size(SS_dosdata,dim=1)
!
   kGID = getGroupID('K-Mesh')
   MyPEinKGroup = getMyPEinGroup(kGID)
!
   if (ig == PrintDOSswitchOff) then
      return
   else if (ig > 0 .and. ig <= GlobalNumAtoms) then
      if (isAtomOnMyProc(ig) .and. MyPEinEGroup == 0 .and. MyPEinKGroup == 0) then
         id = getLocalIndex(ig)
         nume_dos = SS_dosdata(id)%n2
         nvals = SS_dosdata(id)%n1
         Neighbor => getNeighbor(id)
         s = Neighbor%NumShells 
         do ia = 1, SS_dosdata(id)%n3
            if (SS_dosdata(id)%n3 == 1) then
               write(anm,'(a,i7)')'_ato',1000000+ig
            else
               if (ia < 10) then
                  write(anm,'(a,i7,a,i1)')'_ato',1000000+ig,'_c',ia
               else
                  write(anm,'(a,i7,a,i2)')'_ato',1000000+ig,'_c',ia
               endif
            endif
            anm(5:5) = 'm'
            fname = trim(cpath)//'SS_DOS_'//trim(systemid)//anm
            open(unit=fu, file=fname, form='formatted', status='unknown')
            p_dos => SS_dosdata(id)%rarray3(:,:,ia)
            write(fu,'(a,i4)')'SS_DOS Data: Atom #',ig
            write(dos_title(1),'(a,f12.8)')'Fermi energy:',chempot
            write(dos_title(2),'(a,f12.8)')'LIZ Radius:',Neighbor%ShellRad(s)
            write(dos_title(3),'(a,i6)')'Num Energy Points:', nume_dos
!           ----------------------------------------------------------
            call printDOSDATA(fu,nvals,nume_dos,3,p_dos,dos_title)
!           ----------------------------------------------------------
            close(fu)
         enddo
      endif
   else if (ig == PrintEachAtomDOS) then
      n = min(GlobalNumAtoms,maxatoms)
      do gid = 1, n
         if (isAtomOnMyProc(gid) .and. MyPEinEGroup == 0 .and. MyPEinKGroup == 0) then
            id = getLocalIndex(gid)
            nume_dos = SS_dosdata(id)%n2
            nvals = SS_dosdata(id)%n1
            Neighbor => getNeighbor(id)
            s = Neighbor%NumShells 
            do ia = 1, SS_dosdata(id)%n3
               if (SS_dosdata(id)%n3 == 1) then
                  write(anm,'(a,i7)')'_ato',1000000+gid
               else
                  if (ia < 10) then
                     write(anm,'(a,i7,a,i1)')'_ato',1000000+gid,'_c',ia
                  else
                     write(anm,'(a,i7,a,i2)')'_ato',1000000+gid,'_c',ia
                  endif
               endif
               anm(5:5) = 'm'
               fname = trim(cpath)//'SS_DOS_'//trim(systemid)//anm
               open(unit=fu+gid, file=fname, form='formatted', status='unknown')
               p_dos => SS_dosdata(id)%rarray3(:,:,ia)
               write(fu+gid,'(a,i4)')'SS_DOS Data: Atom #',gid
               write(dos_title(1),'(a,f12.8)')'Fermi energy:',chempot
               write(dos_title(2),'(a,f12.8)')'LIZ Radius:',Neighbor%ShellRad(s)
               write(dos_title(3),'(a,i6)')'Num Energy Points:', nume_dos
!              -------------------------------------------------------
               call printDOSDATA(fu+gid,nvals,nume_dos,3,p_dos,dos_title)
!              -------------------------------------------------------
               close(fu+gid)
            enddo
         endif
      enddo
   else
      call ErrorHandler('writeSS_DOS','Invalid global atom index',ig)
   endif
!
   if (GlobalNumAtoms == 1) then
      return
   endif
!
!  ===================================================================
!  Always print average DOS ..........................................
!  if (ig == PrintAverageDOS) then
!  ===================================================================
   nume_dos = SS_dosdata(1)%n2
   nvals = SS_dosdata(1)%n1
   ntypes = getNumAtomTypes()
   if (ntypes < 1) then
      call ErrorHandler('writeSS_DOS','Invalid number of atom types',ntypes)
   endif
   allocate(aver_dos(nvals,nume_dos,ntypes))
!
   aver_dos = ZERO
   do id = 1, LocalNumAtoms
      it = getAtomType(AtomIndex(id))
      do ia = 1, SS_dosdata(id)%n3
         cfac = getLocalSpeciesContent(id,ia)
         do ie = 1, nume_dos
            do iv = 2, nvals
               aver_dos(iv,ie,it) = aver_dos(iv,ie,it) + cfac*SS_dosdata(id)%rarray3(iv,ie,ia)
            enddo
         enddo
      enddo
   enddo
   if (NumPEsInAGroup > 1) then
!     ----------------------------------------------------------------
      call GlobalSumInGroup(aGID, aver_dos, nvals, nume_dos, ntypes)
!     ----------------------------------------------------------------
   endif
!
   do it = 1, ntypes
      n = getNumAtomsOfType(it)
      write(anm,'(a,a1,i3)')trim(getAtomTypeName(it)),'_',100+it
      if (len(trim(getAtomTypeName(it))) == 2) then
         anm(4:4) = 't'
      else
         anm(3:3) = 't'
      endif
      p_dos => aver_dos(:,:,it)
      rfac = ONE/real(n,kind=RealKind)
      p_dos = p_dos*rfac
      do ie = 1, nume_dos
         p_dos(1,ie) = SS_dosdata(1)%rarray3(1,ie,1)
      enddo
!     Print the average DOS of each type .............................
      if (getMyPE() == 0) then
         fname = trim(cpath)//'SS_DOS_'//trim(systemid)//'_average_'//trim(anm)
         open(unit=fu-it, file=fname, form='formatted', status='unknown')
         write(fu-it,'(a,i6,a,a)')'SS_DOS Data: Averaged over ',n,trim(anm),' atoms'
         write(dos_title(1),'(a,f12.8)')'Fermi energy:',chempot
         write(dos_title(2),'(a,f12.8)')'LIZ Radius:',Neighbor%ShellRad(s)
         write(dos_title(3),'(a,i6)')'Num Energy Points:', nume_dos
!        -------------------------------------------------------------
         call printDOSDATA(fu-it,nvals,nume_dos,3,p_dos,dos_title)
!        -------------------------------------------------------------
         close(fu-it)
      endif
   enddo
!
!  Print the averaged total DOS ......................................
   if (ntypes > 1 .and. getMyPE() == 0) then
      do it = 1, ntypes
         n = getNumAtomsOfType(it)
         p_dos => aver_dos(:,:,it)
         rfac = real(n,kind=RealKind)
         p_dos = p_dos*rfac
      enddo
      p_tdos => aver_dos(:,:,1)
      do it = 2, ntypes
         p_dos => aver_dos(:,:,it)
         p_tdos = p_tdos + p_dos
      enddo
      rfac = ONE/real(GlobalNumAtoms,kind=RealKind)
      p_tdos = p_tdos*rfac
      do ie = 1, nume_dos
         p_tdos(1,ie) = SS_dosdata(1)%rarray3(1,ie,1)
      enddo
      fname = trim(cpath)//'SS_DOS_'//trim(systemid)//'_average'
      open(unit=fu-ntypes-1, file=fname, form='formatted', status='unknown')
      write(fu-ntypes-1,'(a,i6,a)')'SS_DOS Data: Averaged over ',GlobalNumAtoms,' atoms'
      write(dos_title(1),'(a,f12.8)')'Fermi energy:',chempot
      write(dos_title(2),'(a,f12.8)')'LIZ Radius:',Neighbor%ShellRad(s)
      write(dos_title(3),'(a,i6)')'Num Energy Points:', nume_dos
!     ----------------------------------------------------------------
      call printDOSDATA(fu-ntypes-1,nvals,nume_dos,3,p_tdos,dos_title)
!     ----------------------------------------------------------------
      close(fu-ntypes-1)
   endif
!
   deallocate(aver_dos)
!
   end subroutine writeSS_DOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeDOS(cpath,systemid,ig)
!  ===================================================================
   use MPPModule, only : getMyPE
!
   use PublicParamDefinitionsModule, only : PrintEachAtomDOS, PrintDOSswitchOff
!
   use PublicTypeDefinitionsModule, only : NeighborStruct
!
   use GroupCommModule, only : getGroupID, getMyPEinGroup
   use GroupCommModule, only : GlobalSumInGroup
!
   use ProcMappingModule, only : isAtomOnMyProc
!
   use Atom2ProcModule, only : getAtom2ProcInGroup, getLocalIndex
!
   use SystemModule, only : getNumAtomTypes, getAtomType, getNumAtomsOfType, &
                            getAtomTypeName
!
   use NeighborModule, only : getNeighbor
!
   use AtomModule, only : getLocalSpeciesContent
!
   implicit none
!
   character (len=*), intent(in) :: cpath
   character (len=*), intent(in) :: systemid
   character (len=150) :: fname
   character (len=60) :: dos_title(3)
   character (len=15) :: anm
!
   integer (kind=IntKind), intent(in) :: ig ! = -1: print DOS for each atom
                                            ! =  0: will not print DOS
                                            ! >  0: print the DOS of the atom with global index ig
!
   integer (kind=IntKind) :: nvals
   integer (kind=IntKind) :: DOSpid
   integer (kind=IntKind) :: s, ie, id, iv, n, gid, kGID, MyPEinKGroup
   integer (kind=IntKind) :: nume_dos, it, ntypes, ia
   integer (kind=IntKind), parameter :: fu = 250, maxatoms = 20
!
   real (kind=RealKind), allocatable, target :: aver_dos(:,:,:)
   real (kind=RealKind), pointer :: p_dos(:,:), p_tdos(:,:)
   real (kind=RealKind) :: rfac, cfac
!
   type (NeighborStruct), pointer :: Neighbor
!
   if (.not.allocated(dosdata)) then
      call ErrorHandler('writeDOS', 'Need to call calValenceDOS first')
   endif
!
!  nume_dos = size(dosdata,dim=2)
!  nvals = size(dosdata,dim=1)
!
   kGID = getGroupID('K-Mesh')
   MyPEinKGroup = getMyPEinGroup(kGID)
!
   if (ig == PrintDOSswitchOff) then
      return
   else if (ig > 0 .and. ig <= GlobalNumAtoms) then
      if (isAtomOnMyProc(ig) .and. MyPEinEGroup == 0 .and. MyPEinKGroup == 0) then
         id = getLocalIndex(ig)
         nume_dos = dosdata(id)%n2
         nvals = dosdata(id)%n1
         Neighbor => getNeighbor(id)
         s = Neighbor%NumShells 
         do ia = 1, dosdata(id)%n3
            if (dosdata(id)%n3 == 1) then
               write(anm,'(a,i7)')'_ato',1000000+ig
            else
               if (ia < 10) then
                  write(anm,'(a,i7,a,i1)')'_ato',1000000+ig,'_c',ia
               else
                  write(anm,'(a,i7,a,i2)')'_ato',1000000+ig,'_c',ia
               endif
            endif
            anm(5:5) = 'm'
            fname = trim(cpath)//'DOS_'//trim(systemid)//anm
            open(unit=fu, file=fname, form='formatted', status='unknown')
            p_dos => dosdata(id)%rarray3(:,:,ia)
            write(fu,'(a,i4)')'DOS Data: Atom #',ig
            write(dos_title(1),'(a,f12.8)')'Fermi energy:',chempot
            write(dos_title(2),'(a,f12.8)')'LIZ Radius:',Neighbor%ShellRad(s)
            write(dos_title(3),'(a,i6)')'Num Energy Points:', nume_dos
!           ----------------------------------------------------------
            call printDOSDATA(fu,nvals,nume_dos,3,p_dos,dos_title)
!           ----------------------------------------------------------
            close(fu)
         enddo
      endif
   else if (ig == PrintEachAtomDOS) then
      n = min(GlobalNumAtoms,maxatoms)
      do gid = 1, n
         if (isAtomOnMyProc(gid) .and. MyPEinEGroup == 0 .and. MyPEinKGroup == 0) then
            id = getLocalIndex(gid)
            nume_dos = dosdata(id)%n2
            nvals = dosdata(id)%n1
            Neighbor => getNeighbor(id)
            s = Neighbor%NumShells 
            do ia = 1, dosdata(id)%n3
               if (dosdata(id)%n3 == 1) then
                  write(anm,'(a,i7)')'_ato',1000000+gid
               else
                  if (ia < 10) then
                     write(anm,'(a,i7,a,i1)')'_ato',1000000+gid,'_c',ia
                  else
                     write(anm,'(a,i7,a,i2)')'_ato',1000000+gid,'_c',ia
                  endif
               endif
               anm(5:5) = 'm'
               fname = trim(cpath)//'DOS_'//trim(systemid)//anm
               open(unit=fu+gid, file=fname, form='formatted', status='unknown')
               p_dos => dosdata(id)%rarray3(:,:,ia)
               write(fu+gid,'(a,i4)')'DOS Data: Atom #',gid
               write(dos_title(1),'(a,f12.8)')'Fermi energy:',chempot
               write(dos_title(2),'(a,f12.8)')'LIZ Radius:',Neighbor%ShellRad(s)
               write(dos_title(3),'(a,i6)')'Num Energy Points:', nume_dos
!              -------------------------------------------------------
               call printDOSDATA(fu+gid,nvals,nume_dos,3,p_dos,dos_title)
!              -------------------------------------------------------
               close(fu+gid)
            enddo
         endif
      enddo
   else
      call ErrorHandler('writeDOS','Invalid global atom index',ig)
   endif
!
   if (GlobalNumAtoms == 1) then
      return
   endif
!
!  ===================================================================
!  Always print average DOS ..........................................
!  if (ig == PrintAverageDOS) then
!  ===================================================================
   nume_dos = dosdata(1)%n2
   nvals = dosdata(1)%n1
   ntypes = getNumAtomTypes()
   if (ntypes < 1) then
      call ErrorHandler('writeDOS','Invalid number of atom types',ntypes)
   endif
   allocate(aver_dos(nvals,nume_dos,ntypes))
!
   aver_dos = ZERO
   do id = 1, LocalNumAtoms
      it = getAtomType(AtomIndex(id))
      do ia = 1, SS_dosdata(id)%n3
         cfac = getLocalSpeciesContent(id,ia)
         do ie = 1, nume_dos
            do iv = 2, nvals
               aver_dos(iv,ie,it) = aver_dos(iv,ie,it) + cfac*dosdata(id)%rarray3(iv,ie,ia)
            enddo
         enddo
      enddo
   enddo
   if (NumPEsInAGroup > 1) then
!     ----------------------------------------------------------------
      call GlobalSumInGroup(aGID, aver_dos, nvals, nume_dos, ntypes)
!     ----------------------------------------------------------------
   endif
   do it = 1, ntypes
      n = getNumAtomsOfType(it)
      write(anm,'(a,a1,i3)')trim(getAtomTypeName(it)),'_',100+it
      if (len(trim(getAtomTypeName(it))) == 2) then
         anm(4:4) = 't'
      else
         anm(3:3) = 't'
      endif
      p_dos => aver_dos(:,:,it)
      rfac = ONE/real(n,kind=RealKind)
      p_dos = p_dos*rfac
      do ie = 1, nume_dos
         p_dos(1,ie) = dosdata(1)%rarray3(1,ie,1)
      enddo
!     Print the average DOS of each type .............................
      if (getMyPE() == 0) then
         fname = trim(cpath)//'DOS_'//trim(systemid)//'_average_'//trim(anm)
         open(unit=fu-it, file=fname, form='formatted', status='unknown')
         write(fu-it,'(a,i6,a,a)')'DOS Data: Averaged over ',n,trim(anm),' atoms'
         write(dos_title(1),'(a,f12.8)')'Fermi energy:',chempot
         write(dos_title(2),'(a,f12.8)')'LIZ Radius:',Neighbor%ShellRad(s)
         write(dos_title(3),'(a,i6)')'Num Energy Points:', nume_dos
!        -------------------------------------------------------------
         call printDOSDATA(fu-it,nvals,nume_dos,3,p_dos,dos_title)
!        -------------------------------------------------------------
         close(fu-it)
      endif
   enddo
!
!  Print the averaged total DOS ......................................
   if (ntypes > 1 .and. getMyPE() == 0) then
      do it = 1, ntypes
         n = getNumAtomsOfType(it)
         p_dos => aver_dos(:,:,it)
         rfac = real(n,kind=RealKind)
         p_dos = p_dos*rfac
      enddo
      p_tdos => aver_dos(:,:,1)
      do it = 2, ntypes
         p_dos => aver_dos(:,:,it)
         p_tdos = p_tdos + p_dos
      enddo
      rfac = ONE/real(GlobalNumAtoms,kind=RealKind)
      p_tdos = p_tdos*rfac
      do ie = 1, nume_dos
         p_tdos(1,ie) = dosdata(1)%rarray3(1,ie,1)
      enddo
      fname = trim(cpath)//'DOS_'//trim(systemid)//'_average'
      open(unit=fu-ntypes-1, file=fname, form='formatted', status='unknown')
      write(fu-ntypes-1,'(a,i6,a)')'DOS Data: Averaged over ',GlobalNumAtoms,' atoms'
      write(dos_title(1),'(a,f12.8)')'Fermi energy:',chempot
      write(dos_title(2),'(a,f12.8)')'LIZ Radius:',Neighbor%ShellRad(s)
      write(dos_title(3),'(a,i6)')'Num Energy Points:', nume_dos
!     ----------------------------------------------------------------
      call printDOSDATA(fu-ntypes-1,nvals,nume_dos,3,p_tdos,dos_title)
!     ----------------------------------------------------------------
      close(fu-ntypes-1)
   endif
!
   deallocate(aver_dos)
!
   end subroutine writeDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDOSDATA(fu,nvals,nume_dos,nt,p_dos,dos_title)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: fu, nvals, nume_dos, nt
!
   character (len=*), intent(in) :: dos_title(nt)
!
   integer (kind=IntKind) :: ie, id, i
!
   real (kind=RealKind), intent(in) :: p_dos(nvals,nume_dos)
!
   if (n_spin_cant == 2) then
      write(fu, '(94(''-''))')
   else if (n_spin_pola == 2) then
      write(fu, '(52(''-''))')
   else
      write(fu, '(32(''-''))')
   endif
   do i = 1, nt
      write(fu, '(a)')dos_title(i)
   enddo
   if (n_spin_cant == 2) then
      write(fu, '(94(''-''),/)')
      write(fu, '(a10,6a14)')'Energy', 'DOS_total', 'DOS_up', 'DOS_down', &
                             'MagDen_X', 'MagDen_Y', 'MagDen_Z'
      if (nvals < 7) then
         call ErrorHandler('writeDOS', 'nvals < 7', nvals)
      endif
   else if (n_spin_pola == 2) then
      write(fu, '(52(''-''))')
      write(fu, '(a10,3a14)')'Energy', 'DOS_total', 'DOS_up', 'DOS_down'
      if (nvals < 3) then
         call ErrorHandler('writeDOS', 'nvals < 3', nvals)
      endif
   else
      write(fu, '(32(''-''))')
      write(fu, '(a14,a14)')'Energy', 'DOS'
      if (nvals < 2) then
         call ErrorHandler('writeDOS', 'nvals < 2', nvals)
      endif
   endif
!
   do ie = 1, nume_dos
      if (n_spin_cant == 2) then
         write(fu, '(f10.6,6e14.6)')p_dos(1,ie),p_dos(2,ie), &
                                    p_dos(6,ie),p_dos(7,ie), &
                                    p_dos(3,ie),p_dos(4,ie), &
                                    p_dos(5,ie)
      else if (n_spin_pola == 2) then
         write(fu, '(f10.6,3e14.6)')p_dos(1,ie),             &
                                    p_dos(2,ie)+p_dos(3,ie), &
                                    p_dos(2,ie),p_dos(3,ie)
      else
         write(fu, '(f14.6,e18.6)')p_dos(1,ie), p_dos(2,ie)
      endif
   enddo
!
   end subroutine printDOSDATA
!  ===================================================================
!  The following are subroutines by xianglin for relativistic calculation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function returnRelSingleSiteDOS(info,e,aux,rfac,redundant) result(dos)
!  ===================================================================
!
!  This function returns single site DOS for a given energy on the real 
!  energy axis.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!
!  In the no-spin-polarized case, a factor of 2 is included in dos.
!  ===================================================================
   use PhysParamModule, only : Boltzmann
   use ScfDataModule, only : Temperature
   use GroupCommModule, only : GlobalSumInGroup
   use RadialGridModule, only : getGrid
   use StepFunctionModule, only : getVolumeIntegration
!   use SSSolverModule, only : computeDOS, getDOS
   use RelSSSolverModule, only : getRelSSDOS, SingleDiracScattering
!   use SSSolverModule, only : getOutsideDOS!, computePhaseShift, getPhaseShift
!
   implicit none
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind), intent(in), optional :: rfac
!
   logical, intent(in), optional :: redundant
   logical :: red
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: jmax_dos, kmax_phi, iend, n, kl, jl, ir, n0
   integer (kind=IntKind) :: is, id, print_dos, atom
!
   real (kind=RealKind) :: dos, sfac, dos_mt, dos_out, tps, rmul, dos1, t0
   real (kind=RealKind), pointer :: ps(:)
   real (kind=RealKind) :: msgbuf(5,NumPEsInEGroup)
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3); print_dos = info(4)
   sfac= ONE !TWO/real(n_spin_pola,kind=RealKind)
   sfac = sfac*getFermiDiracFunc(e,chempot,Boltzmann*Temperature)
!
   if (present(rfac)) then
      rmul = rfac
   else
      rmul = ONE
   endif
!
   if (present(redundant)) then
      red = redundant
   else
      red = .false.
   endif
!
!   energy = adjustEnergy(is,e)
   energy = e
!
   t0 = getTime()
!  -------------------------------------------------------------------
   call SingleDiracScattering(id,energy)
!  -------------------------------------------------------------------
   Timing_SS = Timing_SS + (getTime() - t0)
   NumCalls_SS = NumCalls_SS + 1
!
!   if (AtomicNumber(id) == 0) then
!!     call computeDOS(add_highl_fec=.true.)
!      call computeDOS()
!   else
!      call computeDOS()
!   endif
!  -------------------------------------------------------------------
!   dos_out = sfac*getOutsideDOS()
!  commented since no phase shift is calculated relativistically right now
!  -------------------------------------------------------------------
!   call computePhaseShift()
!   ps => getPhaseShift()
!   kmax_phi = (lmax_phi(id)+1)**2
!   tps = ZERO
!   do kl = 1, kmax_phi
!      tps = tps + ps(kl)
!   enddo
!!
!   if (print_dos > 0) then
!      msgbuf = ZERO
!      msgbuf(1,MyPEinEGroup+1) = real(energy)
!      msgbuf(2,MyPEinEGroup+1) = dos
!      msgbuf(3,MyPEinEGroup+1) = dos_mt
!      msgbuf(4,MyPEinEGroup+1) = dos_out
!      msgbuf(5,MyPEinEGroup+1) = tps
!!     ----------------------------------------------------------------
!      call GlobalSumInGroup(eGID,msgbuf,5,NumPEsInEGroup)
!!     ----------------------------------------------------------------
!      if ( node_print_level >= 0) then
!         do n = 1, NumPEsInEGroup
!!           write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy),dos,dos_mt,dos_out,tps
!            write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')msgbuf(1:5,n)
!         enddo
!      endif
!   endif
!
   dos_r_jl => getRelSSDOS(1,id) !is 1:4
   iend = size(dos_r_jl,1); jmax_dos = size(dos_r_jl,2)
   Grid => getGrid(id)
   do is = 1, n_spin_cant*n_spin_cant
      n0=(is-1)*(LastValue(id)%NumRs*jmax_dos+4)
      dos_r_jl => getRelSSDOS(is,id) !is 1:4
      do jl = 1, jmax_dos
         n = (jl-1)*LastValue(id)%NumRs
         do ir = 1, LastValue(id)%NumRs
            aux(n0+n+ir) = rmul*sfac*dos_r_jl(ir,jl)
         enddo
!        ----------------------------------------------------------------
!        call zcopy(LastValue(id)%NumRs,dos_r_jl(1,jl),1,aux(n+1),1)
!        ----------------------------------------------------------------
      enddo
      n = LastValue(id)%NumRs*jmax_dos
      dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                     jmax_dos, 2, dos_r_jl, dos_mt )
      dos_mt = sfac*dos_mt
!      print*,"single-site dos at real axis", "  is=", is
!      print*,"energy=",real(energy), "dos=", dos, "dos_mt=", dos_mt
      aux(n0+n+1) = rmul*dos
      aux(n0+n+2) = rmul*dos_mt
      aux(n0+n+3) = rmul*dos*energy
      aux(n0+n+4) = CZERO!rmul*dos_out
      if (is==1) then
         dos1=dos
      endif
   enddo
!
   dos=dos1
   end function returnRelSingleSiteDOS
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function returnRelSingleSiteDOS_pole(info,e,aux,rfac) result(dos)
!  ===================================================================
!  This function is not used because the DOS calculated is unstable
!
!  This function returns single site DOS for a given energy on the real 
!  energy axis.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!
!  In the no-spin-polarized case, a factor of 2 is included in dos.
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
   use RadialGridModule, only : getGrid
   use StepFunctionModule, only : getVolumeIntegration
!   use SSSolverModule, only : computeDOS, getDOS
   use RelSSSolverModule, only :getRelSSDOS,getRelSSGreen, getRelSSGreenOut,&
                SingleDiracScattering
!   use SSSolverModule, only : getOutsideDOS!, computePhaseShift, getPhaseShift
!
   implicit none
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind), intent(in), optional :: rfac
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: jmax_dos, kmax_phi, iend, n, kl, jl, ir, n0
   integer (kind=IntKind) :: is, id, print_dos, atom
!
   real (kind=RealKind) :: dos, sfac, dos_mt, dos_out, tps, rmul, dos1, dos_end, t0
   real (kind=RealKind), pointer :: ps(:)
   real (kind=RealKind) :: msgbuf(5,NumPEsInEGroup)
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3); print_dos = info(4)
   sfac= ONE !TWO/real(n_spin_pola,kind=RealKind)
!
!   if (present(rfac)) then
!      rmul = rfac
!   else
!      rmul = ONE
!   endif
!
!   energy = adjustEnergy(is,e)
!   energy = cmplx(e,0.d0)!cmplx(e,TEN2m6,kind=CmplxKind)
    energy = cmplx(e,TEN2m6,kind=CmplxKind)
!
   t0 = getTime()
!  -------------------------------------------------------------------
    call SingleDiracScattering(id,energy)
!  -------------------------------------------------------------------
   Timing_SS = Timing_SS + (getTime() - t0)
   NumCalls_SS = NumCalls_SS + 1
!
!   if (AtomicNumber(id) == 0) then
!!     call computeDOS(add_highl_fec=.true.)
!      call computeDOS()
!   else
!      call computeDOS()
!   endif
!  -------------------------------------------------------------------
!   dos_out = sfac*getOutsideDOS()
!  commented since no phase shift is calculated relativistically right now
!  -------------------------------------------------------------------
!   call computePhaseShift()
!   ps => getPhaseShift()
!   kmax_phi = (lmax_phi(id)+1)**2
!   tps = ZERO
!   do kl = 1, kmax_phi
!      tps = tps + ps(kl)
!   enddo
!!
!   if (print_dos > 0) then
!      msgbuf = ZERO
!      msgbuf(1,MyPEinEGroup+1) = real(energy)
!      msgbuf(2,MyPEinEGroup+1) = dos
!      msgbuf(3,MyPEinEGroup+1) = dos_mt
!      msgbuf(4,MyPEinEGroup+1) = dos_out
!      msgbuf(5,MyPEinEGroup+1) = tps
!!     ----------------------------------------------------------------
!      call GlobalSumInGroup(eGID,msgbuf,5,NumPEsInEGroup)
!!     ----------------------------------------------------------------
!      if ( node_print_level >= 0) then
!         do n = 1, NumPEsInEGroup
!!           write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy),dos,dos_mt,dos_out,tps
!            write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')msgbuf(1:5,n)
!         enddo
!      endif
!   endif
!
   dos_r_jl => getRelSSDOS(1,id) !is 1:4
   iend = size(dos_r_jl,1); jmax_dos = size(dos_r_jl,2)
   Grid => getGrid(id)
!-------------------------------------------------
!  find the correct rmul
!-------------------------------------------------
   dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
         jmax_dos, 2, dos_r_jl, dos_mt, truncated=.false. )
!
!  dos_out = int_0^infty A*r^2*exp(-2*k*x)*Y_00; dos(iend,1)=A*r^2 not used since it's wrong
   dos_end = real(dos_r_jl(iend,1),RealKind)
   if (e > 0) then
      call ErrorHandler("returnRelSingleSiteDOS_pole","unimplemented for poles on positive axis")
   endif
   dos_out = -AIMAG(getRelSSGreenOut(id))
!  SQRT(4.d0*PI)*dos_end/( 4.d0*(SQRT(-e))**3*(Grid%r_mesh(iend))**2 )
   rmul = 1.d0/(dos + dos_out)
   if (MyPE == 0) then
      print*, "dos=", dos, " dos_out=",dos_out, getRelSSGreenOut(id), " dos_total", dos+dos_out
      print*,"e=", e," number of electrons in this resonance: ", dos/(dos+dos_out)
   endif
!
   do is = 1, n_spin_cant*n_spin_cant
      n0=(is-1)*(LastValue(id)%NumRs*jmax_dos+4)
      dos_r_jl => getRelSSDOS(is,id) !is 1:4
      do jl = 1, jmax_dos
         n = (jl-1)*LastValue(id)%NumRs
         do ir = 1, LastValue(id)%NumRs
            aux(n0+n+ir) = rmul*sfac*dos_r_jl(ir,jl)
         enddo
!     ----------------------------------------------------------------
!     call zcopy(LastValue(id)%NumRs,dos_r_jl(1,jl),1,aux(n+1),1)
!     ----------------------------------------------------------------
      enddo
      n = LastValue(id)%NumRs*jmax_dos
      dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                     jmax_dos, 2, dos_r_jl, dos_mt )
      dos_mt = sfac*dos_mt
      aux(n0+n+1) = rmul*dos
      aux(n0+n+2) = rmul*dos_mt
      aux(n0+n+3) = rmul*dos*real(energy,RealKind)
      aux(n0+n+4) = CZERO!rmul*dos_out
      if (is==1) then
         dos1=dos
         if (MyPE == 0) then
            print*,"aux(n0+n+1)=",aux(n0+n+1)
         endif
      endif
   enddo
!
   dos=dos1
!   if (MyPE == 0 .and. e<-.4 ) then
!      open (104,file="dos_r_jl1",action="write")
!      do jl=1,jmax_dos
!         do ir=1,iend
!            write(104,*) Grid%r_mesh(ir), real(dos_r_jl(ir,jl)*rmul)
!         enddo
!      enddo
!      close (104)
!   else
!      open (105,file="dos_r_jl2",action="write")
!      do jl=1,jmax_dos
!         do ir=1,iend
!            write(105,*) Grid%r_mesh(ir), real(dos_r_jl(ir,jl)*rmul)
!         enddo
!      enddo
!      close (105)
!   endif
   end function returnRelSingleSiteDOS_pole
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function returnRelSingleSiteDOSinCP(info,e,aux,cfac) result(dos)
!  ===================================================================
!
!  This function returns single site DOS for a given energy on the real 
!  energy axis.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!
!  In the no-spin-polarized case, a factor of 2 is included in dos.
!  ===================================================================
   use PotentialTypeModule, only : isASAPotential
!
   use RadialGridModule, only : getGrid
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use RelSSSolverModule, only : getRelSSGreen, SingleDiracScattering
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in), optional :: cfac
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: is, id, ks, atom
   integer (kind=IntKind) :: kmax, jmax, iend, n, ir, jl, l, m, kl, klc, p0, n0
!
   real (kind=RealKind) :: dos, sfac, t0
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy, cmul, greenint, greenint_mt, ede
   complex (kind=CmplxKind), pointer :: green(:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:), p1(:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3)
   sfac= ONE !TWO/real(n_spin_pola,kind=RealKind)
!
   if (present(cfac)) then
      cmul = cfac
   else
      cmul = CONE
   endif
!
!  ===================================================================
!  For energy in the upper complex plane, use spherical solver
!  together with irregular solution
!  ===================================================================
   energy = adjustEnergy(is,e)
!
   t0 = getTime()
!  -------------------------------------------------------------------
   call SingleDiracScattering(id,energy)
!  -------------------------------------------------------------------
   Timing_SS = Timing_SS + (getTime() - t0)
   NumCalls_SS = NumCalls_SS + 1
!
!  -------------------------------------------------------------------
   green=>getRelSSGreen(1,id)
!  -------------------------------------------------------------------
!
   iend = size(green,1); kmax = size(green,2); jmax = jofk(kmax)
   Grid => getGrid(id)
!  -------------------------------------------------------------------
   greenint = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,            &
                                         kmax, 2, green, greenint_mt )
!  -------------------------------------------------------------------
   greenint_mt = sfac*greenint_mt
!
   dos = real(SQRTm1*greenint/PI,kind=RealKind)
!
   if (isASAPotential()) then
      greenint = greenint_mt
   endif
!
!
   if (iharris <= 1) then
      ede = energy
   else
      ede = (energy-chempot)
   endif
!
   n = iend*jmax
   do ks = 1, n_spin_cant*n_spin_cant
      n0=(ks-1)*(iend*jmax+4)
      green => getRelSSGreen(ks,id)
      p1 => aux(n0+1:n0+n)
      dos_r_jl => aliasArray2_c(p1,iend,jmax)
      do jl = 1, jmax
         l = lofj(jl); m = mofj(jl)
         kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
         pcv_x => green(1:iend,kl)
         pcv_y => green(1:iend,klc)
         pcv_z => dos_r_jl(1:iend,jl)
         pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
      enddo
   enddo
   aux(n0+n+1) = SQRTm1*cmul*greenint/PI
   aux(n0+n+2) = SQRTm1*cmul*greenint_mt/PI
   aux(n0+n+3) = SQRTm1*cmul*ede*greenint/PI
!  -------------------------------------------------------------------
   greenint = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,       &
                                         kmax, 2, green, greenint_mt, &
                                         truncated=.false. ) - greenint
!  -------------------------------------------------------------------
   aux(n+4) = SQRTm1*cmul*greenint/PI
!  print*,"aux(n0+n+1)=",aux(n0+n+1)
!  print*,"aux(n0+n+2)=",aux(n0+n+2)
!  print*,"aux(n0+n+3)=",aux(n0+n+3)
!  print*,"aux(n0+n+4)=",aux(n0+n+4)
!
   end function returnRelSingleSiteDOSinCP
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function returnRelMultipleSiteDOS(info,e,dos_array,cfac) result(dos)
!  ===================================================================
   use MathParamModule, only : CONE, ZERO, PI
!
   use RelMSSolverModule, only : getRelMSGreenFunction!, getMSGreenMatrix
!
   use PotentialTypeModule, only : isASAPotential
!
   use RadialGridModule, only : getGrid
!
   use AtomModule, only : getLocalEvecNew
!
   use StepFunctionModule, only : getVolumeIntegration
!
!   use SpinRotationModule, only : transformDensityMatrix
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in), optional :: cfac
!
   complex (kind=CmplxKind), intent(inout), target :: dos_array(*)
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: is, id, atom
   integer (kind=IntKind) :: ks, jid, j, l, m, jl, kl, klc, p0, n, i
   integer (kind=IntKind) :: jmax, kmax, iend, ns_sqr, gform
!
   real (kind=RealKind) :: r0j(3), evec_new(3), sfac, dos(2)
   real (kind=RealKind), pointer :: r_mesh(:), pra_y(:,:)
!
   complex (kind=CmplxKind) :: dosmt_pola(2)
   complex (kind=CmplxKind) :: dosws_pola(2), ede, energy, cmul
!
   complex (kind=CmplxKind) :: greenint(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: greenint_mt(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: gm(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: gm_mt(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind), pointer :: green_matrix(:,:,:)
   complex (kind=CmplxKind), pointer :: green(:,:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:), pca_x(:,:), pca_y(:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:), p1(:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3)
   energy = adjustEnergy(is,e)
!
   if (present(cfac)) then
      cmul = cfac
   else
      cmul = CONE
   endif
!
!  =================================================================
!  Note: green is the multiple scattering term of the Green function
!        multiplied by r^2 and is expanded in spherical harmonics.
!
!        For nonspin-polarized case, green, greenint, and
!        greenint_mt also contain a factor 2.0.
!
!        green and greenint are in local spin reference frame
!  =================================================================
!   if (RelativisticFlag == 2) then
!      stop 'Not implemented yet!'
!     --------------------------------------------------------------
  !   green  => getRelGreenFunction(id)
!     --------------------------------------------------------------
!   else
!     --------------------------------------------------------------
      green  => getRelMSGreenFunction(id)
!     --------------------------------------------------------------
!   endif
!
   sfac= ONE !TWO/real(n_spin_pola,kind=RealKind)
   ns_sqr = n_spin_cant*n_spin_cant
   iend = size(green,1); kmax = size(green,2); jmax = jofk(kmax)
!
   if (iend /= LastValue(id)%NumRs) then
      call ErrorHandler('returnRelMultipleSiteDOS','Inconsistent r-mesh size of green and dos_r_jl', &
                        iend,LastValue(id)%NumRs)
   else if (jmax /= LastValue(id)%jmax) then
      call ErrorHandler('returnRelMultipleSiteDOS','Inconsistent jmax size of green and dos_r_jl', &
                        jmax,LastValue(id)%jmax)
   endif
!
   Grid => getGrid(id)
!
   do ks = 1, ns_sqr
!     --------------------------------------------------------------
      greenint(ks) = cmul*sfac*getVolumeIntegration( id, iend, Grid%r_mesh, &
                                                     kmax, 2, green(:,:,ks), greenint_mt(ks) )
!     --------------------------------------------------------------
      greenint_mt(ks) = cmul*sfac*greenint_mt(ks)
   enddo
!
!   if (isASAPotential()) then
!      gm_mt = greenint_mt
!      if (n_spin_cant == 2) then
!        ------------------------------------------------------------
!         call transformDensityMatrix(id,gm_mt)
!        ------------------------------------------------------------
!      endif
!      greenint = greenint_mt
!      gm = gm_mt
!   else
      gm = greenint
      gm_mt = greenint_mt
!      if (n_spin_cant == 2) then
!        ------------------------------------------------------------
!         call transformDensityMatrix(id,gm)
!         call transformDensityMatrix(id,gm_mt)
!        ------------------------------------------------------------
!      endif
!   endif
!
   dosmt_pola = CZERO; dosws_pola = CZERO
   if(n_spin_cant == 2) then
!     ---------------------------------------------------------------
      evec_new(1:3) = getLocalEvecNew(id)
      dosmt_pola(1) = SQRTm1*HALF*(gm_mt(1)+gm_mt(2)*evec_new(1)+gm_mt(3)*evec_new(2)+gm_mt(4)*evec_new(3))/PI
      dosmt_pola(2) = SQRTm1*HALF*(gm_mt(1)-gm_mt(2)*evec_new(1)-gm_mt(3)*evec_new(2)-gm_mt(4)*evec_new(3))/PI
      dosws_pola(1) = SQRTm1*HALF*(gm(1)+gm(2)*evec_new(1)+gm(3)*evec_new(2)+gm(4)*evec_new(3))/PI
      dosws_pola(2) = SQRTm1*HALF*(gm(1)-gm(2)*evec_new(1)-gm(3)*evec_new(2)-gm(4)*evec_new(3))/PI
   else
      dosmt_pola(1) = SQRTm1*greenint_mt(1)/PI
      dosws_pola(1) = SQRTm1*greenint(1)/PI
   endif

   if (node_print_level >= 0) then
      write(6,'(/,''returnRelMultipleSiteDOS: energy ='',2f18.12,'', id ='',i4)')energy,id
      write(6,'(  ''                       Int [Z*(Tau-t)*Z] on MT         ='',2f18.12)') &
            dosmt_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
      write(6,'(  ''                       Int [Z*(Tau-t)*Z] on VP         ='',2f18.12)') &
            dosws_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
      write(6,'(  ''                       Im{-Int [Z*(Tau-t)*Z]/pi on MT} ='',f18.12)')  &
            real(dosmt_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
      write(6,'(  ''                       Im{-Int [Z*(Tau-t)*Z]/pi on VP} ='',f18.12)')  &
            real(dosws_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
   endif ! print_level

   dos = real(dosws_pola,kind=RealKind)

   if (iharris /= 0) then
      call ErrorHandler('returnRelMultipleSiteDOS','Relativistic Harris calculation not implemented')
   endif
!
   if (iharris <= 1) then !harris energy calculation,
      ede = energy
   else
      ede = (energy-chempot)
   endif
!
   n = iend*jmax
   p0 = 0
   do ks = 1, ns_sqr
      p1 => dos_array(p0+1:p0+n)
      dos_r_jl => aliasArray2_c(p1, iend, jmax)
      do jl = 1, jmax
         l = lofj(jl); m = mofj(jl)
         kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
         pcv_x => green(:,kl,ks)
         pcv_y => green(:,klc,ks)
         pcv_z => dos_r_jl(:,jl)
         pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
      enddo
      p0 = p0 + n
      dos_array(p0+1)  = SQRTm1*greenint(ks)/PI
      dos_array(p0+2) = SQRTm1*greenint_mt(ks)/PI
      dos_array(p0+3) = SQRTm1*ede*greenint(ks)/PI
      dos_array(p0+4) = CZERO
      p0 = p0 + 4
!      if (isDensityMatrixNeeded) then
!         kmax = size(green_matrix,1)
!         pca_x => green_matrix(:,:,ks)
!         p1 => dos_array(p0+1:p0+kmax*kmax)
!         pca_y => aliasArray2_c(p1,kmax,kmax)
!         pca_y = cmul*SQRTm1*pca_x
!         p0 = p0 + kmax*kmax
!      endif
   enddo
   end function returnRelMultipleSiteDOS

!  ===================================================================
   subroutine calRelIntegratedDOS(efermi)
!  ===================================================================
!   use PublicParamDefinitionsModule, only : ButterFly
!
   use PhysParamModule, only : Boltzmann
!
   use MPPModule, only : endMPP
!
   use GroupCommModule, only : GlobalMaxInGroup
!
   use ScfDataModule, only : ErBottom, ErTop, EiBottom, EiTop, Temperature
   use ScfDataModule, only : NumEs, ContourType, eGridType, isReadEmesh, getEmeshFileName
   use ScfDataModule, only : isKKR, isSSIrregularSolOn
   use ScfDataModule, only : NumSS_IntEs
!
   use AtomModule, only : getLocalSpeciesContent
   !
   use PotentialTypeModule, only : isFullPotential
   !
   use KreinModule, only : isLloydOn
   !
   use ContourModule, only : initContour, endContour, setupContour
   use ContourModule, only : printContour
!
!   use MSSolverModule, only : computeMSGreenFunction
   use RelMSSolverModule, only : computeRelMST
!
!  use ValenceDensityModule, only : printValenceDensity
!
   implicit none
!
   logical :: useIrregularSolution = .false.
   logical :: relativity
!
   integer (kind=IntKind) :: id, BadFermiEnergy, is, info(4), ns, ia
   integer (kind=IntKind), parameter :: MaxIterations = 20
   !
   real (kind=RealKind), intent(inout) :: efermi
   real (kind=RealKind) :: efermi_old, kBT, rfac
   real (kind=RealKind) :: Lloyd_factor(n_spin_pola), ssDOS, msDOS(2)
   real (kind=RealKind) :: int_dos(n_spin_cant*n_spin_pola,LocalNumAtoms)
   real (kind=RealKind) :: last_dos(n_spin_cant*n_spin_pola,LocalNumAtoms)
!
   complex (kind=CmplxKind) :: eLast
   if (isLloydOn()) then
      call ErrorHandler("calRelIntegratedDOS","isLloydOn=.ture., unimplemented case")
   endif
!
   kBT = Temperature*Boltzmann
!
!  ===================================================================
!  Initialize the enegy contour module
!  ===================================================================
   if (isReadEmesh()) then
!     ----------------------------------------------------------------
      call initContour(getEmeshFileName(), stop_routine, max_print_level)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call initContour( ContourType, eGridType, NumEs, Temperature,   &
                        stop_routine, maxval(print_level(1:LocalNumAtoms)), .true. )
!      ----------------------------------------------------------------
!
!     ================================================================
!     set up the enegy contour that ends at (efermi,0) for Gaussian grids
!     on a semi-circle contour or ends at (efermi,Eibottom) for uniform
!     grids on arectangular contour.
!     New feature: No extra last energy point is added!
!     ----------------------------------------------------------------
      call setupContour( ErBottom, efermi, EiBottom, EiTop )
!     ----------------------------------------------------------------
   endif

   if (node_print_level >= 0) then
!     ----------------------------------------------------------------
      call printContour()
!     ----------------------------------------------------------------
   endif

   do id =  1,LocalNumAtoms
!     ----------------------------------------------------------------
      call zeroElectroStruct(IntegrValue(id))
!     ----------------------------------------------------------------
   enddo
   !
   call calMultipleScatteringIDOS(useIrregularSolution, relativity=.true.)
   if (.not.useIrregularSolution) then
!     ================================================================
!     compute the single-site DOS from e<0 (shallow bound states) 
!     ================================================================
      if ( ErBottom < ZERO .and. isPole ) then
         do id =  1,LocalNumAtoms
            call zeroElectroStruct(ssIntegrValue(id))
         enddo
!        -------------------------------------------------------------
         call calSingleScatteringIDOS(LowerContour=.true., relativity=.true.) !xianglin
!        -------------------------------------------------------------
         do id =  1,LocalNumAtoms
!           ----------------------------------------------------------
            call addElectroStruct(CONE,ssIntegrValue(id),IntegrValue(id))
!           ----------------------------------------------------------
         enddo
         write(6,'(/,a)')'Rel: IDOS of the SS_Pole (e<0) term'
         do id =  1,LocalNumAtoms 
            do ia = 1, ssIntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', MS IDOS_mt_Pole = ',real(ssIntegrValue(id)%dos_mt(is,ia)),&
                          ', MS IDOS_ws_Pole = ',real(ssIntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
!     ================================================================
!     compute the DOS arising from the single site scattering term along
!     the real energy axis from E = 0 upto E = efermi
!     ================================================================
      do id =  1,LocalNumAtoms
!        -------------------------------------------------------------
         call zeroElectroStruct(ssIntegrValue(id))
!        -------------------------------------------------------------
      enddo
      call calSingleScatteringIDOS(Ebegin=0.00001d0,                  &
                                   Eend=efermi+8.0d0*log(10.0d0)*kBT, &
                                   relativity=.true.)
      if ( node_print_level >= 0) then
         write(6,'(/,a)')'Rel: IDOS of the MS term'
         do id =  1,LocalNumAtoms 
            do ia = 1, IntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', MS IDOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)),&
                          ', MS IDOS_ws = ',real(IntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
         write(6,'(/,a)')'Rel: IDOS of the SS (e>0) term'
         do id =  1,LocalNumAtoms 
            do ia = 1, ssIntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', SS IDOS_mt = ',real(ssIntegrValue(id)%dos_mt(is,ia)),&
                          ', SS IDOS_ws = ',real(ssIntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
      do id =  1,LocalNumAtoms
!        -------------------------------------------------------------
         call addElectroStruct(CONE,ssIntegrValue(id),IntegrValue(id))
!        -------------------------------------------------------------
      enddo
      if ( node_print_level >= 0) then
         write(6,'(/,a)')'Rel: After adding the IDOS of the SS term to the IDOS of the MS term'
         do id =  1,LocalNumAtoms 
            do ia = 1, IntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia ',ia,', is = ',is,   &
                          ', MS+SS IDOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)),&
                          ', MS+SS IDOS_ws = ',real(IntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
   else
      call ErrorHandler("calRelIntegratedDOS","useIrregularSolution not implemented")
   endif
!
!  note transformed after add ssIntegerValue and IntegerValue, different from Yang's old code
   if (n_spin_cant == 2 .and. RelativisticFlag .ne. 2) then !no transform for relativistic calculation
      do id =  1,LocalNumAtoms
!        -------------------------------------------------------------
         call transformElectroStruct(id,IntegrValue(id))
!        -------------------------------------------------------------
      enddo
   endif

!  ===================================================================
!  Iterate the ending point of the energy contour to find the Fermi energy
!  ===================================================================
   BadFermiEnergy = 1
   LOOP_LastE: do while (BadFermiEnergy > 0 .and. BadFermiEnergy <= MaxIterations)
   !     ===============================================================
   !     Solve the multiple scattering problem for e = eLast, which is
   !     set to be (efermi,0.001) in KKR case.
   !     ===============================================================
      if (isKKR()) then
         eLast = cmplx(efermi,0.005d0,kind=CmplxKind)
      else
         eLast = cmplx(efermi,0.000d0,kind=CmplxKind)
      endif
   !
   !     ===============================================================
   !     Calculate the DOS of the multiple scattering term for e = eLast,
   !     store the data in LastValue structure, and add the single site
   !     results (at efermi, not eLast) to it.
   !          LastValue(id)%dos = is the MST part of the DOS at (chempot,eib) for atom id
   !     ===============================================================
      if ( node_print_level >= 0) then
         write(6,'(/,a,2d15.8)')'M.S. Term DOS at the last energy: ',eLast
      endif
      do id =  1,LocalNumAtoms
   !        -------------------------------------------------------------
         call zeroElectroStruct(LastValue(id))
   !        -------------------------------------------------------------
      enddo
!
      do is = 1, n_spin_pola/n_spin_cant
!         if (isFullPotential()) then
!           =========================================================
!           In this case, the multiple scattering module returns DOS of
!           Z*(tau-t)*Z, instead of Z*tau*Z-Z*J. One needs to add the single
!           site DOS to the result.
!           ----------------------------------------------------------
            call computeRelMST(adjustEnergy(is,eLast))
!           ----------------------------------------------------------
!         else
!           =========================================================
!           In this case, the multiple scattering module returns DOS of
!           Z*tau*Z-Z*J.
!           ----------------------------------------------------------
!            call computeMSGreenFunction(is,adjustEnergy(is,eLast),    &
!                                        add_Gs=.true.,isSphSolver=.true.)
!           ----------------------------------------------------------
!         endif
         do id =  1,LocalNumAtoms
            info(1) = is; info(2) = id; info(3) = -1; info(4) = 1
!           ----------------------------------------------------------
            msDOS = returnRelMultipleSiteDOS(info,eLast,wk_dos)
            call calElectroStruct(info,4,wk_dos,LastValue(id))
!           ----------------------------------------------------------
            if ( node_print_level >= 0) then
               do ia = 1, LastValue(id)%NumSpecies
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', M.Site DOS_mt = ',real(LastValue(id)%dos_mt(is,ia)), &
                          ', M.Site DOS_ws = ',real(LastValue(id)%dos(is,ia))
               enddo
            endif
         enddo
      enddo
!      if (isFullPotential()) then
!        ============================================================
!        In this case, the multiple scattering module returns DOS of
!        Z*(tau-t)*Z, instead of Z*tau*Z-Z*J. One needs to add the single
!        site DOS to the result.
!        Solve the single scattering problem for e = efermi.
!        ssLastValue(id)%dos = is the single site part of the DOS at efermi on the
!                              real energy axis for atom id
!        ============================================================
         do id =  1,LocalNumAtoms
!           ----------------------------------------------------------
            call zeroElectroStruct(ssLastValue(id))
!           ----------------------------------------------------------
            do is = 1, 1!n_spin_pola*n_spin_pola !major change by xianglin
               info(1) = is; info(2) = id; info(3) = -1; info(4) = -1
!              -------------------------------------------------------
               ssDOS = returnRelSingleSiteDOS(info,efermi,wk_dos)
               call calElectroStruct(info,4,wk_dos,ssLastValue(id))
!              -------------------------------------------------------
               if ( node_print_level >= 0) then
                  write(6,'(/,a,d15.8)')'S.S. Term DOS at the last energy: ',efermi
                  do ia = 1, ssLastValue(id)%NumSpecies
                     write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', S.Site DOS_mt = ',real(ssLastValue(id)%dos_mt(is,ia)), &
                          ', S.Site DOS_ws = ',real(ssLastValue(id)%dos(is,ia))
                  enddo
               endif
            enddo
!
!           =========================================================
!           Now, add the single site DOS to the multiple scattering DOS term at last
!           energy point
!           ----------------------------------------------------------
            call addElectroStruct(ONE,ssLastValue(id),LastValue(id))
!           ----------------------------------------------------------
!
            if (n_spin_cant == 2 .and. RelativisticFlag .ne. 2) then
!              -------------------------------------------------------
               call transformElectroStruct(id,LastValue(id))
!              -------------------------------------------------------
            endif

            if ( node_print_level >= 0) then
               write(6,'(/,a,2d15.8)')'M.S.+S.S. DOS at the last energy: ',eLast
               do ia = 1, LastValue(id)%NumSpecies
                  do is = 1, n_spin_pola*n_spin_cant
                     write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                             ', MS+SS  DOS_mt = ',real(LastValue(id)%dos_mt(is,ia)), &
                             ', MS+SS  DOS_ws = ',real(LastValue(id)%dos(is,ia))
                  enddo
               enddo
            endif
         enddo
!      endif
!
      last_dos = ZERO; int_dos = ZERO
      do id = 1, LocalNumAtoms
         do ia = 1, LastValue(id)%NumSpecies
            rfac = getLocalSpeciesContent(id,ia)
            do is = 1, n_spin_cant*n_spin_pola
               last_dos(is,id) =  last_dos(is,id) + rfac*real(LastValue(id)%dos(is,ia),kind=RealKind)
               int_dos(is,id) =  int_dos(is,id) + rfac*real(IntegrValue(id)%dos(is,ia),kind=RealKind)
            enddo
         enddo
      enddo
!
      efermi_old = efermi
!     ================================================================
!     Compute the Fermi energy and efermi will be renewed.
!     ----------------------------------------------------------------
      call mufind(efermi,efermi_old,int_dos,last_dos,BadFermiEnergy,Lloyd_factor)
!     ----------------------------------------------------------------
      do id = 1, LocalNumAtoms
!        -------------------------------------------------------------
         call addElectroStruct(efermi-efermi_old,LastValue(id),IntegrValue(id))
!        -------------------------------------------------------------
      enddo
      if ( node_print_level >= 0) then
         write(6,'(/,a)')'After mufind, the integrated DOS are:'
         do id =  1,LocalNumAtoms
            do ia = 1, IntegrValue(id)%NumSpecies
               do is = 1, n_spin_pola*n_spin_cant
                  write(6,'(3(a,i2),2(a,d15.8))')'id = ',id,', ia = ',ia,', is = ',is, &
                          ', MS+SS  DOS_mt = ',real(IntegrValue(id)%dos_mt(is,ia)), &
                          ', MS+SS  DOS_ws = ',real(IntegrValue(id)%dos(is,ia))
               enddo
            enddo
         enddo
      endif
!
      if ( .not.isIterateEfOn .or. isLloydOn() ) then
         exit LOOP_LastE
      else
!        -------------------------------------------------------------
         call GlobalMaxInGroup(aGID,BadFermiEnergy)
!        -------------------------------------------------------------
      endif
   enddo Loop_LastE
!
!  -------------------------------------------------------------------
!  call updateValenceDOS(efermi,Lloyd_factor)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  end Energy contour modules
!  -------------------------------------------------------------------
   call endContour()
!  -------------------------------------------------------------------
!
   if (stop_routine == 'calRelIntegratedDOS') then
      call syncAllPEs()
      call endMPP()
      call StopHandler('calRelIntegratedDOS')
   endif
!
   end subroutine calRelIntegratedDOS

!  ===================================================================
   subroutine calQuadraticPoles(site,LocalNumSpecies,lkkr,lphi,       &
                                n_spin,eb,et,ldp,NumPoles,Poles,      &
                                PanelOnZero, isRel)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                               Ten2m8, FOURTH, CZERO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : solveSingleScattering, getJostMatrix, getSineMatrix
   use RelSSSolverModule, only : SingleDiracScattering, getRelJostMatrix
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix, &
                                     solveQuadraticEquation, getEigenValue,   &
                                     solveLinearEquation
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: site,LocalNumSpecies
   integer (kind=IntKind), intent(in) :: lkkr
   integer (kind=IntKind), intent(in) :: lphi
   integer (kind=IntKind), intent(in) :: n_spin,ldp
   integer (kind=IntKind), intent(out) :: NumPoles(LocalNumSpecies,n_spin)
!
   integer (kind=IntKind) :: ia, ie, is, iw, n, kmax_kkr, NumWindows, info
   integer (kind=IntKind) :: i, j, l, m, kl, lp, mp, klp, nv
!
   logical, optional, intent(in) :: PanelOnZero
   logical, optional, intent(in) :: isRel
   logical :: Relativity
   logical :: isZeroInterval = .false.
   logical, parameter :: isGeneral = .false.
!
!   real (kind=RealKind), parameter :: Delta = 0.002d0
   real (kind=RealKind) :: WindowWidth != 5.0*Delta!20.*Delta!4.05*Delta
   real (kind=RealKind) :: Delta, t0
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0
   real (kind=RealKind), intent(out) :: Poles(ldp,LocalNumSpecies,n_spin)
   real (kind=RealKind) :: Poles_loc(ldp,LocalNumSpecies,n_spin)
   integer (kind=IntKind) :: NumPoles_all(LocalNumSpecies,n_spin,NumPEsInEGroup)
!
   complex (kind=CmplxKind) :: e, kappa, a2l, a2lp
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), allocatable :: s0(:,:,:), s1(:,:,:), s2(:,:,:)
   complex (kind=CmplxKind), allocatable :: s1p(:,:,:)
   complex (kind=CmplxKind), pointer :: pv(:)
!
   if (eb > et) then
      call ErrorHandler('calQuadraticPoles','eb > et',eb,et)
   endif
!
   WindowWidth = (et-eb)/(nw_Pole_new)
   NumWindows = nwPG

   Delta=0.2d0*WindowWidth
!
   if (.not.present(PanelOnZero)) then
      isZeroInterval = .false.
   endif
   if (.not. present(isRel)) then
      Relativity = .false.
   else
      Relativity = isRel
   endif
!
   if(isZeroInterval) then
      NumWindows = 1
!   else
!      NumWindows = int((et-eb)/WindowWidth)
   endif
!
   if(isZeroInterval) then
!     de = FOURTH*(et-eb)
      de = HALF*(et-eb)
   else
      de = Delta
   endif
   de2 = de*TWO; dede2 = de*de*TWO
!
   NumPoles = 0
   if (Relativity) then ! Needs to be fixed for the CPA case
      kmax_kkr = 2*(lkkr+1)**2
      allocate( s0(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s1(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s2(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s1p(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      call initQuadraticMatrix(kmax_kkr,isGeneral)
!
      do is = 1,n_spin
         if (2*(lkkr+1)**2 /= kmax_kkr) then
            deallocate(s0, s1, s2, s1p)
            kmax_kkr =2*(lkkr+1)**2
            allocate( s0(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s1(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s2(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s1p(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            call endQuadraticMatrix()
            call initQuadraticMatrix(kmax_kkr)
         endif
!
         do iw = MyPEinEGroup*nwPG+1, nwPG*(MyPEinEGroup+1)
            !Epoints calculated by local E processor
            w0 = eb + (iw-1)*WindowWidth
            e0 = w0 + (HALF)*WindowWidth
            if (.not.isZeroInterval.and.(abs(e0) < Ten2m6 .or.        &
                abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6)) then
               write(6,*)'e0 shifted'
               e0 = e0 - HALF*de
            else if(isZeroInterval) then
               e0 = ZERO
            endif
!
            if(isZeroInterval) then
               s0 = CZERO
            else
               e = cmplx(e0,ZERO,kind=CmplxKind)
               kappa = sqrt(e)
!
               t0 = getTime()
!              ----------------------------------------------------------
               call SingleDiracScattering (site,e)
!              ----------------------------------------------------------
               Timing_SS = Timing_SS + (getTime() - t0)
               NumCalls_SS = NumCalls_SS + 1
!
               ia = 1 ! Temporary fix for the CPA case
!              ----------------------------------------------------------
               sin_mat => getRelJostMatrix(site)
!              ----------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0(1,1,ia),1)
!              ----------------------------------------------------------
            endif
!
            e = cmplx(e0+de,ZERO,kind=CmplxKind)
!
            t0 = getTime()
!           -------------------------------------------------------------
            call SingleDiracScattering (site,e)
!           -------------------------------------------------------------
            Timing_SS = Timing_SS + (getTime() - t0)
            NumCalls_SS = NumCalls_SS + 1
!
            ia = 1 ! Temporary fix for the CPA case
!           -------------------------------------------------------------
            sin_mat => getRelJostMatrix(site)
!           -------------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s2(1,1,ia),1)
!           -------------------------------------------------------------
!
            e = cmplx(e0-de,ZERO,kind=CmplxKind)
!
            t0 = getTime()
!           -------------------------------------------------------------
            call SingleDiracScattering (site,e)
!           -------------------------------------------------------------
            Timing_SS = Timing_SS + (getTime() - t0)
            NumCalls_SS = NumCalls_SS + 1
!
            ia = 1 ! Temporary fix for the CPA case
!           -------------------------------------------------------------
            sin_mat => getRelJostMatrix(site)
!           -------------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s1p(1,1,ia),1)
!           -------------------------------------------------------------
!
            s1 = (s2 - s1p)/de2
            s2 = (s2 + s1p - TWO*s0)/dede2
!
            do ia = 1, LocalNumSpecies
               if(isZeroInterval) then
!              -------------------------------------------------------
                  call solveLinearEquation(s1(:,:,ia),s2(:,:,ia),info)
!              -------------------------------------------------------
               else
!              -------------------------------------------------------
                  call solveQuadraticEquation(s0(:,:,ia),s1(:,:,ia),s2(:,:,ia),info)
!              -------------------------------------------------------
               endif
!
               if (info /= 0) then
                  stop 'Error in s0, s1, s2'
               endif
!
!              ----------------------------------------------------------
               pv => getEigenValue(nv)
!              ----------------------------------------------------------
               if ( node_print_level >= 0) then
                  write(6,'("calQuadraticPoles: iw",t40,"=", i5,i5)')iw
               endif
               do ie = 1, kmax_kkr*2
                  if (abs(aimag(pv(ie))) < Ten2m4) then
!                 if (aimag(pv(ie)) > ZERO .and. real(pv(ie),kind=RealKind) + e0 > ZERO) then
!                  if (aimag(sqrt(pv(ie)+e0)) < ZERO) then
                     pe = real(pv(ie),kind=RealKind) + e0
!                if (pe >= w0 .and. pe <= w0+WindowWidth .and. pe > 0.001d0) then
                     if (pe >= w0 .and. pe <= w0+WindowWidth ) then
                        NumPoles(ia,is) = NumPoles(ia,is) + 1
                        n = NumPoles(ia,is)
                        Poles_loc(n,ia,is) = pe
                     endif
                  endif
               enddo ! ie
            enddo ! ia
         enddo ! iw
      enddo !is loop
   else !relativistic or not
      kmax_kkr = (lkkr+1)**2
      allocate( s0(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s1(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s2(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s1p(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      call initQuadraticMatrix(kmax_kkr,isGeneral)
!
      do is = 1,n_spin
         if ((lkkr+1)**2 /= kmax_kkr) then
            deallocate(s0, s1, s2, s1p)
            kmax_kkr = (lkkr+1)**2
            allocate( s0(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s1(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s2(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s1p(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            call endQuadraticMatrix()
            call initQuadraticMatrix(kmax_kkr)
         endif
!
         do iw = MyPEinEGroup*nwPG+1, nwPG*(MyPEinEGroup+1)
            w0 = eb + (iw-1)*WindowWidth
            e0 = w0 + (HALF)*WindowWidth
            if (.not.isZeroInterval.and.(abs(e0) < Ten2m6 .or.        &
                abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6)) then
               write(6,*)'e0 shifted'
               e0 = e0 - HALF*de
            else if(isZeroInterval) then
               e0 = ZERO
            endif
!
            if(isZeroInterval) then
               s0 = CZERO
            else
               e = cmplx(e0,ZERO,kind=CmplxKind)
               kappa = sqrt(e)
!
               t0 = getTime()
!              ----------------------------------------------------------
               call solveSingleScattering(is, site, e, CZERO)
!              ----------------------------------------------------------
               Timing_SS = Timing_SS + (getTime() - t0)
               NumCalls_SS = NumCalls_SS + 1
!
               do ia = 1, LocalNumSpecies
                  sin_mat => getJostMatrix(spin=is,site=site,atom=ia)
!                 -------------------------------------------------------
                  call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0(1,1,ia),1)
!                 -------------------------------------------------------
               enddo
            endif
!
            e = cmplx(e0+de,ZERO,kind=CmplxKind)
!
            t0 = getTime()
!           -------------------------------------------------------------
            call solveSingleScattering(is, site, e, CZERO)
!           -------------------------------------------------------------
            Timing_SS = Timing_SS + (getTime() - t0)
            NumCalls_SS = NumCalls_SS + 1
!
            do ia = 1, LocalNumSpecies
               sin_mat => getJostMatrix(spin=is,site=site,atom=ia)
!              ----------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s2(1,1,ia),1)
!              ----------------------------------------------------------
            enddo
!
            e = cmplx(e0-de,ZERO,kind=CmplxKind)
!
            t0 = getTime()
!           -------------------------------------------------------------
            call solveSingleScattering(is, site, e, CZERO)
!           -------------------------------------------------------------
            Timing_SS = Timing_SS + (getTime() - t0)
            NumCalls_SS = NumCalls_SS + 1
!
            do ia = 1, LocalNumSpecies
               sin_mat => getJostMatrix(spin=is,site=site,atom=ia)
!              ----------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s1p(1,1,ia),1)
!              ----------------------------------------------------------
            enddo
!
            s1 = (s2 - s1p)/de2
            s2 = (s2 + s1p - TWO*s0)/dede2
!
            do ia = 1, LocalNumSpecies
               if(isZeroInterval) then
!                 -------------------------------------------------------
                  call solveLinearEquation(s1(:,:,ia),s2(:,:,ia),info)
!                 -------------------------------------------------------
               else
!                 -------------------------------------------------------
                  call solveQuadraticEquation(s0(:,:,ia),s1(:,:,ia),s2(:,:,ia),info)
!                 -------------------------------------------------------
               endif
!
               if (info /= 0) then
                  stop 'Error in s0, s1, s2'
               endif
!
!              ----------------------------------------------------------
               pv => getEigenValue(nv)
!              ----------------------------------------------------------
               if ( node_print_level >= 0) then
                  write(6,'("calQuadraticPoles: iw",t40,"=", i5,i5)')iw
               endif
               do ie = 1, kmax_kkr*2
!                 if (abs(aimag(pv(ie))) < Ten2m8) then
!                 if (aimag(pv(ie)) > ZERO .and. real(pv(ie),kind=RealKind) + e0 > ZERO) then
                  if (aimag(sqrt(pv(ie)+e0)) < ZERO) then
                     pe = real(pv(ie),kind=RealKind) + e0
                     if (pe >= w0 .and. pe <= w0+WindowWidth .and. pe > 0.01d0) then
                        NumPoles(ia,is) = NumPoles(ia,is) + 1
                        n = NumPoles(ia,is)
                        Poles_loc(n,ia,is) = pe
!      write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
                     endif
                  endif
               enddo ! ie
            enddo ! ia
         enddo ! iw
      enddo ! is
   endif !relativistic or not 
   NumPoles_all=0
   NumPoles_all(:,:,MyPEinEGroup+1)=NumPoles(:,:)
!  print*,"NumPoles(nd,ns), before global sum", NumPoles(:,:)
!  print*,"NumPoles_all(nd,ns,NumPEsInEGroup), before global sum", NumPoles_all(:,:,:)
   call GlobalSumInGroup(eGID, NumPoles_all,LocalNumSpecies,n_spin,NumPEsInEGroup)
   call GlobalSumInGroup(eGID, NumPoles,LocalNumSpecies,n_spin)
!  print*,"NumPoles(nd,ns), after global sum", NumPoles(ia,is)
!  print*,"NumPoles_all(nd,ns,NumPEsInEGroup), after global sum", NumPoles_all(:,:,:)
   do is = 1,n_spin
      do ia = 1,LocalNumSpecies
         if ( NumPoles(ia,is) >  2*(lmax_kkr_max+1)**2 ) then
            call ErrorHandler('calQuadraticPoles','NumPoles > 2*(lmax_kkr_max+1)**2',&
                              NumPoles(ia,is),2*(lmax_kkr_max+1)**2)
         endif
         lp = 0
         do i = 1, MyPEinEGroup
            lp = lp + NumPoles_all(ia,is,i)
         enddo ! find the total number of poles before local E points
         do kl =1, NumPoles_all(ia,is,MyPEinEGroup+1)
            Poles(lp+kl,ia,is)=Poles_loc(kl,ia,is)
         enddo
      enddo
   enddo
   call GlobalSumInGroup(eGID, Poles,ldp,LocalNumSpecies,n_spin)
!
   if ( node_print_level >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,25x,a)')'*********************************'
      write(6,'( 25x,a )')'* Output from calQuadraticPoles *'
      write(6,'(25x,a,/)')'*********************************'
!
      do is = 1,n_spin
         do ia = 1,LocalNumSpecies
            write(6,'(/,"NumPEsInEGroup",t40,"=", i5)')NumPEsInEGroup
            write(6,'("MyPEinEGroup",t40,"=", i5)')MyPEinEGroup
            write(6,'("calQuadraticPoles: Num. of Windows",t40,"=", i5)')nw_Pole_new
            write(6,'("pole-searching within Window number",t40,"=", i5,i5)') &
                 MyPEinEGroup*nwPG+1, nwPG*(MyPEinEGroup+1)
            write(6,'("ia",t40,"=", i5)')ia
            write(6,'("is",t40,"=", i5)')is
            write(6,'("NumPoles",t40,"=", i5)')NumPoles(ia,is)
            do i = 1, NumPoles(ia,is)
               write(6,'(/,"Pole",i5,t40,"=", f18.12)')i, Poles(i,ia,is)
            enddo
         enddo
      enddo
      write(6,'(80(''-''))')
   endif
!
!   if (MyPE == 0) then
!      do is = 1,n_spin
!         do ia = 1,LocalNumSpecies
!            print*,"calQuadraticPoles"
!            print*, "ia=", ia, "  is=",is
!            print*,"NumPoles(ia,is)=", NumPoles(ia,is)
!            do i = 1,NumPoles(ia,is)
!               print*, i," Poles=",Poles(i,ia,is)
!            enddo
!            print*,"NumPEsInEGroup=", NumPEsInEGroup," MyPEinEGroup=", MyPEinEGroup
!            print*, "NumPoles_all="
!            print*, NumPoles_all(ia,is,:)
!         enddo
!      enddo
!   endif
!
   nullify(pv)
   call endQuadraticMatrix()
!
   deallocate( s0, s1, s2, s1p )
!
   end subroutine calQuadraticPoles

!  ===================================================================
   subroutine calQuadraticPoles_cmplx(site,LocalNumSpecies,lkkr,lphi, &
                                n_spin,eb,et,ldp,NumPoles,Poles,      &
                                PanelOnZero, isRel)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                               Ten2m8, FOURTH, CZERO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : solveSingleScattering, getJostMatrix, getSineMatrix
   use RelSSSolverModule, only : SingleDiracScattering, getRelJostMatrix
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix, &
                                     solveQuadraticEquation, getEigenValue,   &
                                     solveLinearEquation
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, LocalNumSpecies
   integer (kind=IntKind), intent(in) :: lkkr
   integer (kind=IntKind), intent(in) :: lphi
   integer (kind=IntKind), intent(in) :: n_spin,ldp
   integer (kind=IntKind), intent(out) :: NumPoles(LocalNumSpecies,n_spin)
!  Note the "NumPoles" is just a local name, and actually is NumPoles_plus when this subroutine is called
!
   integer (kind=IntKind) :: ia, ie, is, iw, n, kmax_kkr, NumWindows, info
   integer (kind=IntKind) :: i, j, l, m, kl, lp, mp, klp, nv
!
   logical, optional, intent(in) :: PanelOnZero
   logical, optional, intent(in) :: isRel
   logical :: Relativity
   logical :: isZeroInterval = .false.
   logical, parameter :: isGeneral = .false.
!
!   real (kind=RealKind), parameter :: Delta = 0.002d0
   real (kind=RealKind) :: WindowWidth != 5.0*Delta!20.*Delta!4.05*Delta
   real (kind=RealKind) :: Delta, t0
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0
   complex (kind=CmplxKind), intent(out) :: Poles(ldp,LocalNumSpecies,n_spin)
!  Note the "Poles" is just a local name, and actually is Poles_plus when this subroutine is called
   complex (kind=CmplxKind) :: Poles_loc(ldp,LocalNumSpecies,n_spin)
   integer (kind=IntKind) :: NumPoles_all(LocalNumSpecies,n_spin,NumPEsInEGroup)
!
   complex (kind=CmplxKind) :: e, kappa, a2l, a2lp
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), allocatable :: s0(:,:,:), s1(:,:,:), s2(:,:,:)
   complex (kind=CmplxKind), allocatable :: s1p(:,:,:)
   complex (kind=CmplxKind), pointer :: pv(:)
!
   if (eb > et) then
      call ErrorHandler('calQuadraticPoles','eb > et',eb,et)
   endif
!
   WindowWidth = (et-eb)/(nw_Pole_new_plus)
   NumWindows = nwPG_plus
   Delta=0.2d0*WindowWidth
!
   if (.not.present(PanelOnZero)) then
      isZeroInterval = .false.
   endif
   if (.not. present(isRel)) then
      Relativity = .false.
   else
      Relativity = isRel
   endif
!
   if(isZeroInterval) then
      NumWindows = 1
!   else
!      NumWindows = int((et-eb)/WindowWidth)
   endif
!
   if(isZeroInterval) then
!     de = FOURTH*(et-eb)
      de = HALF*(et-eb)
   else
      de = Delta
   endif
   de2 = de*TWO; dede2 = de*de*TWO
!
   NumPoles = 0
   if (Relativity) then
      kmax_kkr = 2*(lkkr+1)**2
      allocate( s0(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s1(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s2(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s1p(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      call initQuadraticMatrix(kmax_kkr,isGeneral)
!
      do is = 1,n_spin
         if (2*(lkkr+1)**2 /= kmax_kkr) then
            deallocate(s0, s1, s2, s1p)
            kmax_kkr =2*(lkkr+1)**2
            allocate( s0(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s1(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s2(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s1p(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            call endQuadraticMatrix()
            call initQuadraticMatrix(kmax_kkr)
         endif
!
         do iw = MyPEinEGroup*nwPG_plus+1, nwPG_plus*(MyPEinEGroup+1)
            w0 = eb + (iw-1)*WindowWidth
            e0 = w0 + (HALF)*WindowWidth
            if (.not.isZeroInterval.and.(abs(e0) < Ten2m6 .or.        &
                abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6)) then
                write(6,*)'e0 shifted'
               e0 = e0 - HALF*de
            else if(isZeroInterval) then
               e0 = ZERO
            endif
!
            if (isZeroInterval) then
               s0 = CZERO
            else
               e = cmplx(e0,ZERO,kind=CmplxKind)
               kappa = sqrt(e)
!
               t0 = getTime()
!              -------------------------------------------------------
               call SingleDiracScattering (site,e)
!              -------------------------------------------------------
               Timing_SS = Timing_SS + (getTime() - t0)
               NumCalls_SS = NumCalls_SS + 1
!
               ia = 1 ! A temporary fix for the CPA case
!              -------------------------------------------------------
               sin_mat => getRelJostMatrix(site)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0(1,1,ia),1)
   !           -------------------------------------------------------
            endif
!
            e = cmplx(e0+de,ZERO,kind=CmplxKind)
!
            t0 = getTime()
!           ----------------------------------------------------------
            call SingleDiracScattering (site,e)
!           ----------------------------------------------------------
            Timing_SS = Timing_SS + (getTime() - t0)
            NumCalls_SS = NumCalls_SS + 1
!
            ia = 1 ! A temporary fix for the CPA case in relativistic MST
!           ----------------------------------------------------------
            sin_mat => getRelJostMatrix(site)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s2(1,1,ia),1)
!           ----------------------------------------------------------
!
            e = cmplx(e0-de,ZERO,kind=CmplxKind)
!
            t0 = getTime()
!           ----------------------------------------------------------
            call SingleDiracScattering (site,e)
!           ----------------------------------------------------------
            Timing_SS = Timing_SS + (getTime() - t0)
            NumCalls_SS = NumCalls_SS + 1
!
            ia = 1 ! A temporary fix for the CPA case
!           ----------------------------------------------------------
            sin_mat => getRelJostMatrix(site)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s1p(1,1,ia),1)
!           ----------------------------------------------------------
!
            s1 = (s2 - s1p)/de2
            s2 = (s2 + s1p - TWO*s0)/dede2
!
            do ia = 1, LocalNumSpecies
               if(isZeroInterval) then
!                 ----------------------------------------------------
                  call solveLinearEquation(s1(:,:,ia),s2(:,:,ia),info)
!                 ----------------------------------------------------
               else
!                 ----------------------------------------------------
                  call solveQuadraticEquation(s0(:,:,ia),s1(:,:,ia),s2(:,:,ia),info)
!                 ----------------------------------------------------
               endif
   !
               if (info /= 0) then
                  stop 'Error in s0, s1, s2'
               endif
   !
   !           -------------------------------------------------------
               pv => getEigenValue(nv)
   !           -------------------------------------------------------
               if ( node_print_level >= 0) then
                  write(6,'("calQuadraticPoles_cmplx: iw",t40,"=", i5,i5)')iw
               endif
               do ie = 1, kmax_kkr*2
                  if (abs(aimag(pv(ie))) < 0.1d0) then
   !              if (aimag(pv(ie)) > ZERO .and. real(pv(ie),kind=RealKind) + e0 > ZERO) then
   !               if (aimag(sqrt(pv(ie)+e0)) < ZERO) then
                     pe = real(pv(ie),kind=RealKind) + e0
                     if (pe >= w0 .and. pe <= w0+WindowWidth .and. aimag(sqrt(pv(ie)+e0)) < ZERO) then
                   ! if (pe >= w0 .and. pe <= w0+WindowWidth .and. pe > 0.003d0 .and. & 
                              ! aimag(sqrt(pv(ie)+e0)) < ZERO ) then
                        NumPoles(ia,is) = NumPoles(ia,is) + 1
                        n = NumPoles(ia,is)
                        Poles_loc(n,ia,is) = pv(ie)+e0
!   if (MyPE == 0) then
!      write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
!   endif
                     endif
                  endif
               enddo ! ie
            enddo ! ia
         enddo ! iw
      enddo ! is
   else !relativistic or not
      kmax_kkr = (lkkr+1)**2
      allocate( s0(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s1(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s2(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      allocate( s1p(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
      call initQuadraticMatrix(kmax_kkr,isGeneral)
!
      do is = 1,n_spin
         if ((lkkr+1)**2 /= kmax_kkr) then
            deallocate(s0, s1, s2, s1p)
            kmax_kkr = (lkkr+1)**2
            allocate( s0(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s1(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s2(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            allocate( s1p(1:kmax_kkr,1:kmax_kkr,1:LocalNumSpecies) )
            call endQuadraticMatrix()
            call initQuadraticMatrix(kmax_kkr)
         endif
!
         do iw = MyPEinEGroup*nwPG_plus+1, nwPG_plus*(MyPEinEGroup+1)
            w0 = eb + (iw-1)*WindowWidth
            e0 = w0 + (HALF)*WindowWidth
            if (.not.isZeroInterval.and.(abs(e0) < Ten2m6 .or.        &
                abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6)) then
               write(6,*)'e0 shifted'
               e0 = e0 - HALF*de
            else if(isZeroInterval) then
               e0 = ZERO
            endif
!
            if (isZeroInterval) then
               s0 = CZERO
            else
               e = cmplx(e0,ZERO,kind=CmplxKind)
               kappa = sqrt(e)
!
               t0 = getTime()
!              -------------------------------------------------------
               call solveSingleScattering(is, site, e, CZERO)
!              -------------------------------------------------------
               Timing_SS = Timing_SS + (getTime() - t0)
               NumCalls_SS = NumCalls_SS + 1
!
               do ia = 1, LocalNumSpecies
                  sin_mat => getJostMatrix(spin=is,site=site,atom=ia)
!                 ----------------------------------------------------
                  call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s0(1,1,ia),1)
!                 ----------------------------------------------------
               enddo
            endif
!
            e = cmplx(e0+de,ZERO,kind=CmplxKind)
!
            t0 = getTime()
!           ----------------------------------------------------------
            call solveSingleScattering(is, site, e, CZERO)
!           ----------------------------------------------------------
            Timing_SS = Timing_SS + (getTime() - t0)
            NumCalls_SS = NumCalls_SS + 1
!
            do ia = 1, LocalNumSpecies
               sin_mat => getJostMatrix(spin=is,site=site,atom=ia)
!              -------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s2(1,1,ia),1)
!              -------------------------------------------------------
            enddo
!
            e = cmplx(e0-de,ZERO,kind=CmplxKind)
!
            t0 = getTime()
!           -------------------------------------------------------------
            call solveSingleScattering(is, site, e, CZERO)
!           -------------------------------------------------------------
            Timing_SS = Timing_SS + (getTime() - t0)
            NumCalls_SS = NumCalls_SS + 1
!
            do ia = 1, LocalNumSpecies
               sin_mat => getJostMatrix(spin=is,site=site,atom=ia)
!              ----------------------------------------------------------
               call zcopy(kmax_kkr*kmax_kkr,sin_mat,1,s1p(1,1,ia),1)
!              ----------------------------------------------------------
            enddo
!
            s1 = (s2 - s1p)/de2
            s2 = (s2 + s1p - TWO*s0)/dede2
!
            do ia = 1, LocalNumSpecies
               if(isZeroInterval) then
!                 -------------------------------------------------------
                  call solveLinearEquation(s1(:,:,ia),s2(1,1,ia),info)
!                 -------------------------------------------------------
               else
!                 -------------------------------------------------------
                  call solveQuadraticEquation(s0(:,:,ia),s1(:,:,ia),s2(:,:,ia),info)
!                 -------------------------------------------------------
               endif
!
               if (info /= 0) then
                  stop 'Error in s0, s1, s2'
               endif
!
!              ----------------------------------------------------------
               pv => getEigenValue(nv)
!              ----------------------------------------------------------
               if ( node_print_level >= 0) then
                  write(6,'("calQuadraticPoles_cmplx: iw",t40,"=", i5,i5)')iw
               endif
               do ie = 1, kmax_kkr*2
!                 if (abs(aimag(pv(ie))) < Ten2m8) then
!                 if (aimag(pv(ie)) > ZERO .and. real(pv(ie),kind=RealKind) + e0 > ZERO) then
                  if (aimag(sqrt(pv(ie)+e0)) < ZERO) then
                     pe = real(pv(ie),kind=RealKind) + e0
                     if (pe >= w0 .and. pe <= w0+WindowWidth .and. pe > 0.01d0) then
                        NumPoles(ia,is) = NumPoles(ia,is) + 1
                        n = NumPoles(ia,is)
                        Poles_loc(n,ia,is) = pv(ie)+e0
!    write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
                     endif
                  endif
               enddo ! ie
            enddo ! ia
         enddo ! iw
      enddo ! is
   endif !relativistic or not 
!
   NumPoles_all=0
   NumPoles_all(:,:,MyPEinEGroup+1)=NumPoles(:,:)
!  print*,"NumPoles(nd,ns), before global sum complex", NumPoles(:,:)
!  print*,"NumPoles_all(nd,ns,NumPEsInEGroup), before global sum", NumPoles_all(:,:,:)
   call GlobalSumInGroup(eGID, NumPoles_all,LocalNumSpecies,n_spin,NumPEsInEGroup)
   call GlobalSumInGroup(eGID, NumPoles,LocalNumSpecies,n_spin)
!  print*,"NumPoles(nd,ns), after global sum", NumPoles(:,:)
!  print*,"NumPoles_all(nd,ns,NumPEsInEGroup), after global sum complex", NumPoles_all(:,:,:)
   do is = 1,n_spin
      do ia = 1,LocalNumSpecies
         if ( NumPoles(ia,is) >  2*(lmax_kkr_max+1)**2 ) then
            call ErrorHandler('calQuadraticPoles_cmplx','NumPoles > 2*(lmax_kkr_max+1)**2',&
                              NumPoles(ia,is),2*(lmax_kkr_max+1)**2)
         endif
         lp = 0
         do i = 1, MyPEinEGroup
            lp = lp + NumPoles_all(ia,is,i)
         enddo ! find the total number of poles before local E points
         do kl = 1, NumPoles_all(ia,is,MyPEinEGroup+1)
            Poles(lp+kl,ia,is)=Poles_loc(kl,ia,is)
         enddo
      enddo
   enddo
   call GlobalSumInGroup(eGID, Poles,ldp,LocalNumSpecies,n_spin)
!
   if ( node_print_level >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,25x,a)')'***************************************'
      write(6,'( 25x,a )')'* Output from calQuadraticPoles_cmplx *'
      write(6,'(25x,a,/)')'***************************************'
!
      do is = 1,n_spin
         do ia = 1,LocalNumSpecies
            write(6,'(/,"NumPEsInEGroup",t40,"=", i5)') NumPEsInEGroup
            write(6,'("MyPEinEGroup",t40,"=", i5)')MyPEinEGroup
            write(6,'("calQuadraticPoles_cmplx: Num. of Windows",t40,"=", i5)')nw_Pole_new_plus
            write(6,'("pole-searching within Window number",t40,"=", i5,i5)') &
                 MyPEinEGroup*nwPG_plus+1, nwPG_plus*(MyPEinEGroup+1)
            write(6,'("ia",t40,"=", i5)')ia
            write(6,'("is",t40,"=", i5)')is
            write(6,'("NumPoles",t40,"=", i5)')NumPoles(ia,is)
            do i = 1, NumPoles(ia,is)
               write(6,'("Pole",i5,t40,"=", 2f18.12)')i, Poles(i,ia,is)
            enddo
         enddo
      enddo
      write(6,'(80(''-''))')
   endif
!   if (MyPE == 0) then
!      do is = 1,n_spin
!         do ia = 1,LocalNumSpecies
!            print*,"calQuadraticPoles_cmplx"
!            print*, "ia=", ia, "  is=",is
!            print*,"NumPoles(ia,is)=", NumPoles(ia,is)
!            do i = 1,NumPoles(ia,is)
!               print*, i," Poles=",Poles(i,ia,is)
!            enddo
!            print*,"NumPEsInEGroup=", NumPEsInEGroup," MyPEinEGroup=", MyPEinEGroup
!            print*, "NumPoles_all="
!            print*, NumPoles_all(ia,is,:)
!         enddo
!      enddo
!   endif
   nullify(pv)
   call endQuadraticMatrix()
!
   deallocate( s0, s1, s2, s1p )
!
   end subroutine calQuadraticPoles_cmplx
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFermiDiracFunc(z,mu,kBT) result(fd)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: mu
   real (kind=RealKind), intent(in) :: kBT
   real (kind=RealKind), parameter :: temp_tol = TEN2m6
!
   real (kind=RealKind), intent(in) :: z
   real (kind=RealKind) :: fd
!
   if (kBT < temp_tol) then
      if (z <= mu) then
         fd = ONE
      else
         fd = ZERO
      endif
   else
      fd = ONE/(ONE+exp((z-mu)/kBT))
   endif
!
   end function getFermiDiracFunc
!  ===================================================================
end module GFMethodModule
