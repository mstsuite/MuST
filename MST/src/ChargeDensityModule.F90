module ChargeDensityModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use GroupCommModule, only : GlobalSumInGroup, GlobalMaxInGroup
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use MathParamModule, only : ZERO, ONE, TWO, FOUR, CZERO, CONE, HALF, &
                               PI2, PI4, THIRD, Ten2m2, Ten2m6, Ten2m7, Ten2m8, TEN2m4, &
                               Ten2m9, Ten2m10, TEN2m12, TEN2m5, SQRT_PI, Y0, SQRTm1
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
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
   use IntegerFactorsModule, only : lofj, mofj
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use SystemVolumeModule, only : getSystemVolume
!
   use RadialGridModule, only : getNumRmesh, getGrid
!
   use AtomModule, only : getLocalAtomicNumber
!
   use PolyhedraModule, only : getVolume
!
   use PotentialTypeModule, only : isASAPotential
!
   use CoreStatesModule, only : getDeepCoreDensity, getSemiCoreDensity
!
public:: initChargeDensity,         &
         endChargeDensity,          &
         getSphChargeDensity,       &
         getChargeDensity,          &
         getSphMomentDensity,       &
         getMomentDensity,          &
         getMultipoleMoment,        &
         constructChargeDensity,    &
         getRhoLmax,                &
         getPseudoRad,              &
         getPseudoNumRPts,          &
         printChargeDensity,        &
         printChargeDensity_L,      &
         printMomentDensity_L,      &
         isChargeComponentZero,     &
         getChargeComponentFlag,    &
         getNeutralChargeDensity,   &
         isSphericalChargeDensity,  &
         getVPCharge,               &
         getVPMomSize,              &
         setVPCharge,               &
         setVPMomSize,              &
         updateTotalDensity,        &
         getChargeDensityAtPoint,   &
         getMomentDensityAtPoint
!
   interface getChargeDensity
      module procedure getChrgDen_L, getChrgDen_Lj
   end interface getChargeDensity
!
   interface getMomentDensity
      module procedure getMomDen_L, getMomDen_Lj
   end interface getMomentDensity
!
   interface getMultipoleMoment
      module procedure getMM, getMM_j, getMM_idj
   end interface getMultipoleMoment
!
   interface printChargeDensity
      module procedure printDensity_r, printDensity_c
   end interface printChargeDensity
!
private
!  ==================================================================
!  Charge densities are defined as
!  rho = (rho_tilda + rho0_neutral) + (rho_pseudo - rho0_neutral)
!  ------------------------------------------------------------------
!
   type ChargeDensityStruct
      integer (kind=IntKind) :: L_ID
      integer (kind=IntKind) :: G_ID
      integer (kind=IntKind) :: lmax
      integer (kind=IntKind) :: jmax
      integer (kind=IntKind) :: jend
      integer (kind=IntKind) :: NumRPts
      integer (kind=IntKind) :: NumRPts_Pseudo
      integer (kind=IntKind) :: NumSpecies
!
      real (kind=RealKind) :: RcutPseudo
      real (kind=RealKind) :: ChargeVP
!
      integer (kind=IntKind), pointer :: ChargeCompFlag(:)
!
      real (kind=RealKind), pointer :: r_mesh(:)
!
      real (kind=RealKind), pointer :: rhoSph_Total(:,:)
      real (kind=RealKind), pointer :: momSph_Total(:,:)
      real (kind=RealKind), pointer :: der_rhoSph_Total(:,:)
      real (kind=RealKind), pointer :: der_momSph_Total(:,:)
!
      real (kind=RealKind), pointer :: rhoSph_TotalOld(:,:)
      real (kind=RealKind), pointer :: momSph_TotalOld(:,:)
!
      complex (kind=CmplxKind), pointer :: rhoL_Total(:,:,:)
      complex (kind=CmplxKind), pointer :: momL_Total(:,:,:)
      complex (kind=CmplxKind), pointer :: der_rhoL_Total(:,:,:)
      complex (kind=CmplxKind), pointer :: der_momL_Total(:,:,:)
!
      complex (kind=CmplxKind), pointer :: rhoL_TotalOld(:,:,:)
      complex (kind=CmplxKind), pointer :: momL_TotalOld(:,:,:)
!
      real (kind=RealKind), pointer :: rhoSph_Pseudo(:,:)
      real (kind=RealKind), pointer :: momSph_Pseudo(:,:)
!
      complex (kind=CmplxKind), pointer :: rhoL_Pseudo(:,:,:)
      complex (kind=CmplxKind), pointer :: momL_Pseudo(:,:,:)
!
      real (kind=RealKind), pointer :: rhoSph_Tilda(:,:)
      real (kind=RealKind), pointer :: momSph_Tilda(:,:)
!
      complex (kind=CmplxKind), pointer :: rhoL_Tilda(:,:,:)
      complex (kind=CmplxKind), pointer :: momL_Tilda(:,:,:)
!
      real (kind=RealKind), pointer :: rho_Core(:,:,:)
      real (kind=RealKind), pointer :: rho_SemiCore(:,:,:)
      real (kind=RealKind), pointer :: der_rho_Core(:,:,:)
      real (kind=RealKind), pointer :: der_rho_SemiCore(:,:,:)
!
      complex (kind=CmplxKind), pointer :: rhoL_Valence(:,:,:)
      complex (kind=CmplxKind), pointer :: momL_Valence(:,:,:,:)
!
      complex (kind=CmplxKind), pointer :: rhoL_ValPseudo(:,:,:)
!
      complex (kind=CmplxKind), pointer :: multipole_mom(:,:)
!
   end type ChargeDensityStruct
!
   logical :: Initialized = .false.
   logical :: isSphericalCharge = .true.
   logical :: isChargeSymmOn = .false.
   logical :: isFitChargeOn = .false.
   logical :: isMTFP = .false.
   logical :: rad_derivative = .false.
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: MaxLocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: jmax_max, lmax_max, NumRPts_max
   integer (kind=IntKind) :: jmax_glb
   integer (kind=IntKind), allocatable :: print_level(:)
!
   integer (kind=IntKind) :: Fit_Method = 1
!
!   integer (kind=IntKind), allocatable :: lofj(:)
!   integer (kind=IntKind), allocatable :: mofj(:)
!  integer (kind=IntKind), allocatable :: flag_jl(:)
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind) :: MM_size
   integer (kind=IntKind), allocatable, target :: MM_table_line(:)
!
   integer (kind=IntKind), parameter :: n_inter = 5 ! order of polynomial
                                                    ! interpolation
   real(kind=RealKind) :: rho0_neutral, rhoint
   real(kind=RealKind), parameter :: rho_tol = TEN2m6
!
   complex (kind=CmplxKind), pointer :: MultipoleMom(:,:)
!
   real (kind=RealKind), allocatable, target :: sqrt_r(:)
   real (kind=RealKind), allocatable, target :: den_r(:), den_r1(:)
   complex(kind=CmplxKind), allocatable, target :: den_in(:), den_out(:)
   complex(kind=CmplxKind), allocatable, target :: ws_ylm(:), ws_grad_ylm(:,:)
   complex(kind=CmplxKind), allocatable, target :: ws_den(:,:), ws_der_den(:,:)
!
   type (ChargeDensityStruct), allocatable, target :: ChargeDensityList(:)
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
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initChargeDensity( na, gindex, lmax_rho, pola, cant, iprint, isGGA)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
   use Atom2ProcModule, only : getMaxLocalNumAtoms
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid
   use SystemModule, only : getNumAtoms, getLmaxRho, getNumAlloyElements
   use AtomModule, only : getRcutPseudo, getLocalNumSpecies
   use PotentialTypeModule, only : isFullPotential, isMuffinTinFullPotential
   use ScfDataModule, only : isChargeSymm, isFittedChargeDen
   use ValenceDensityModule, only : getValenceElectronDensity, getValenceMomentDensity
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: gindex(na)
   integer (kind=IntKind), intent(in) :: lmax_rho(na)
   integer (kind=IntKind), intent(in) :: pola, cant
   integer (kind=IntKind), intent(in) :: iprint(na)
!
   logical, intent(in), optional :: isGGA
!
   integer (kind=IntKind) :: id, ia, nr, ns, l, jl, m, lmax, jmax, ir, jmt, ig
   integer (kind=IntKind), allocatable :: DataSize(:), DataSizeL(:)
!
   real (kind=RealKind) :: r_ps
   real (kind=RealKind), pointer :: r_mesh(:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   type (GridStruct), pointer :: Grid
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
   LocalNumAtoms = na
   n_spin_pola = pola
   n_spin_cant = cant
   ns = n_spin_pola
   GlobalNumAtoms = getNumAtoms()
   MaxLocalNumAtoms = getMaxLocalNumAtoms()
   isFitChargeOn = isFittedChargeDen()
   isMTFP = isMuffinTinFullPotential()
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   allocate( print_level(LocalNumAtoms) )
   print_level(:) = iprint(:)
!
   if (present(isGGA)) then
      rad_derivative = isGGA
   else
      rad_derivative = .false.
   endif
!
   lmax_max = maxval(lmax_rho)
   jmax_max = ((lmax_max+1)*(lmax_max+2))/2
!
   if ( lmax_max /= 0 .or. isFullPotential() ) then
      isSphericalCharge = .false.
   endif
!
   allocate( ws_ylm((lmax_max+1)*(lmax_max+1)) )
   allocate( ws_grad_ylm((lmax_max+1)*(lmax_max+1),3) )
   allocate( ws_den(n_inter,jmax_max), ws_der_den(n_inter,jmax_max) )
!
   allocate( DataSize(LocalNumAtoms) )
   do id = 1,LocalNumAtoms
      DataSize(id)  = (getNumRmesh(id)+1)*getLocalNumSpecies(id)
   enddo
!
   isChargeSymmOn = isChargeSymm()
!  ==================================================================
!  generate Storage space for spherical charge densities
!  ==================================================================
   if (.not.isDataStorageExisting('NewSphericalElectronDensity')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'NewSphericalElectronDensity',  &
                             DataSize(1:LocalNumAtoms),RealType)  ! en extra space is used
                                                                  ! storing ValenceVPCharge(1)
!     ----------------------------------------------------------------
      call setDataStorage2Value('NewSphericalElectronDensity',ZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         nr = getNumRmesh(id)+1
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'NewSphericalElectronDensity',nr)
!        -------------------------------------------------------------
      enddo
      if (rad_derivative) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'SphericalElectronDensityDer',  &
                                DataSize(1:LocalNumAtoms),RealType)
!        -------------------------------------------------------------
         call setDataStorage2Value('SphericalElectronDensityDer',ZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)+1
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'SphericalElectronDensityDer',nr)
!           ----------------------------------------------------------
         enddo
      endif
   endif
!
   if (n_spin_pola == 2) then
      if (.not.isDataStorageExisting('NewSphericalMomentDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'NewSphericalMomentDensity', &
                                DataSize(1:LocalNumAtoms),RealType) ! en extra space is used
                                                      ! storing ValenceVPCharge(2)
!        -------------------------------------------------------------
         call setDataStorage2Value('NewSphericalMomentDensity',ZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)+1
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'NewSphericalMomentDensity',nr)
!           ----------------------------------------------------------
         enddo
         if (rad_derivative) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'SphericalMomentDensityDer', &
                                   DataSize(1:LocalNumAtoms),RealType)
!           ----------------------------------------------------------
            call setDataStorage2Value('SphericalMomentDensityDer',ZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)+1
!              -------------------------------------------------------
               call setDataStorageLDA(id,'SphericalMomentDensityDer',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
   endif
!
   if (.not.isDataStorageExisting('OldSphericalElectronDensity')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'OldSphericalElectronDensity',  &
                             DataSize(1:LocalNumAtoms),RealType)  ! en extra space is used
                                                 ! storing ValenceVPCharge(1)
!     ----------------------------------------------------------------
      call setDataStorage2Value('OldSphericalElectronDensity',ZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         nr = getNumRmesh(id)+1
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'OldSphericalElectronDensity',nr)
!        -------------------------------------------------------------
      enddo
   endif
!
   if (n_spin_pola == 2) then
      if (.not.isDataStorageExisting('OldSphericalMomentDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'OldSphericalMomentDensity', &
                                DataSize(1:LocalNumAtoms),RealType) ! en extra space is used
                                                      ! storing ValenceVPCharge(2)
!        -------------------------------------------------------------
         call setDataStorage2Value('OldSphericalMomentDensity',ZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)+1
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'OldSphericalMomentDensity',nr)
!           ----------------------------------------------------------
         enddo
      endif
   endif
!
   if ( isSphericalCharge ) then
      do id = 1,LocalNumAtoms
         DataSize(id)  = getNumRmesh(id)*getLocalNumSpecies(id)
      enddo
!
      if (.not.isDataStorageExisting('NewL-ElectronDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'NewL-ElectronDensity',  &
                                DataSize,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('NewL-ElectronDensity',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'NewL-ElectronDensity',nr)
!           ----------------------------------------------------------
         enddo
         if (rad_derivative) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'L-ElectronDensityDer',  &
                                   DataSize,ComplexType)
!           ----------------------------------------------------------
            call setDataStorage2Value('L-ElectronDensityDer',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'L-ElectronDensityDer',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
!
      if (.not.isDataStorageExisting('OldL-ElectronDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'OldL-ElectronDensity',  &
                                DataSize,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('OldL-ElectronDensity',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'OldL-ElectronDensity',nr)
!           ----------------------------------------------------------
         enddo
      endif
!
      if (n_spin_pola == 2) then
         if (.not.isDataStorageExisting('NewL-MomentDensity')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'NewL-MomentDensity', &
                                   DataSize,ComplexType)
!           ----------------------------------------------------------
            call setDataStorage2Value('NewL-MomentDensity',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'NewL-MomentDensity',nr)
!              -------------------------------------------------------
            enddo
            if (rad_derivative) then
!              -------------------------------------------------------
               call createDataStorage(LocalNumAtoms,'L-MomentDensityDer', &
                                      DataSize,ComplexType)
!              -------------------------------------------------------
               call setDataStorage2Value('L-MomentDensityDer',CZERO)
!              -------------------------------------------------------
               do id = 1,LocalNumAtoms
                  nr = getNumRmesh(id)
!                 ----------------------------------------------------
                  call setDataStorageLDA(id,'L-MomentDensityDer',nr)
!                 ----------------------------------------------------
               enddo
            endif
         endif
      endif
!
      if (n_spin_pola == 2) then
         if (.not.isDataStorageExisting('OldL-MomentDensity')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'OldL-MomentDensity', &
                                   DataSize,ComplexType)
!           ----------------------------------------------------------
            call setDataStorage2Value('OldL-MomentDensity',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'OldL-MomentDensity',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
!
   endif
!
   allocate( ChargeDensityList(LocalNumAtoms) )
!
   NumRPts_max = 0
   do id = 1,LocalNumAtoms
!
      Grid => getGrid(id)
      p_CDL => ChargeDensityList(id)
      p_CDL%L_ID = id
      p_CDL%G_ID = gindex(id)
      p_CDL%lmax = lmax_rho(id)
      p_CDL%jmax = ((lmax_rho(id)+1)*(lmax_rho(id)+2))/2
      p_CDL%NumRPts = Grid%jend  !May28_2014 getNumRmesh(id)
      p_CDL%jend    = Grid%jend
      p_CDL%r_mesh => getRmesh(id)
      nr = p_CDL%NumRPts
!
      p_CDL%NumSpecies = getLocalNumSpecies(id)
!
      p_CDL%rhoSph_Total => getDataStorage( id,'NewSphericalElectronDensity', &
                                            nr+1, getLocalNumSpecies(id), RealMark )
      if (n_spin_pola == 2) then
         p_CDL%momSph_Total => getDataStorage( id,'NewSphericalMomentDensity', &
                                               nr+1, getLocalNumSpecies(id), RealMark )
      endif
!
      if (rad_derivative) then
         p_CDL%der_rhoSph_Total => getDataStorage( id,'SphericalElectronDensityDer', &
                                                   nr+1, getLocalNumSpecies(id), RealMark )
         if (n_spin_pola == 2) then
            p_CDL%der_momSph_Total => getDataStorage( id,'SphericalMomentDensityDer',&
                                                      nr+1, getLocalNumSpecies(id), RealMark )
         endif
      endif
!
      p_CDL%rhoSph_TotalOld => getDataStorage( id,'OldSphericalElectronDensity', &
                                               nr+1, getLocalNumSpecies(id), RealMark )
      if (n_spin_pola == 2) then
         p_CDL%momSph_TotalOld => getDataStorage( id,'OldSphericalMomentDensity', &
                                                  nr+1, getLocalNumSpecies(id), RealMark )
      endif
!
      if ( isSphericalCharge ) then
         p_CDL%rhoL_Total => getDataStorage( id,'NewL-ElectronDensity', &
                                             nr, p_CDL%jmax, getLocalNumSpecies(id), ComplexMark )
         if (n_spin_pola == 2) then
            p_CDL%momL_Total => getDataStorage( id,'NewL-MomentDensity', &
                                                nr, p_CDL%jmax, getLocalNumSpecies(id), ComplexMark )
         endif
!
         if (rad_derivative) then
            p_CDL%der_rhoL_Total => getDataStorage( id,'L-ElectronDensityDer', &
                                                    nr, p_CDL%jmax, getLocalNumSpecies(id), ComplexMark )
            if (n_spin_pola == 2) then
               p_CDL%der_momL_Total => getDataStorage( id,'L-MomentDensityDer', &
                                                       nr, p_CDL%jmax, getLocalNumSpecies(id), ComplexMark )
            endif
         endif
!
         p_CDL%rhoL_TotalOld => getDataStorage( id,'OldL-ElectronDensity', &
                                                nr, p_CDL%jmax, getLocalNumSpecies(id), ComplexMark )
         do ia = 1, p_CDL%NumSpecies
            p_CDL%rhoL_TotalOld(1:nr,1,ia) = p_CDL%rhoSph_TotalOld(1:nr,ia)/Y0
         enddo
         if (n_spin_pola == 2) then
            p_CDL%momL_TotalOld => getDataStorage( id,'OldL-MomentDensity', &
                                                   nr, p_CDL%jmax, getLocalNumSpecies(id), ComplexMark )
         endif
!
         p_CDL%rhoL_Valence => getValenceElectronDensity(id)
         if (n_spin_pola == 2) then
            p_CDL%momL_Valence => getValenceMomentDensity(id)
         endif
      endif
!
      allocate( p_CDL%ChargeCompFlag(1:p_CDL%jmax) )
      nullify( p_CDL%multipole_mom )
!
      NumRPts_max = max(NumRPts_max,nr)
!
   enddo
!
   if ( .not.isSphericalCharge ) then
      allocate( DataSizeL(LocalNumAtoms) )
      do id = 1,LocalNumAtoms
         p_CDL => ChargeDensityList(id)
         DataSizeL(id) = getNumRmesh(id)*p_CDL%jmax*getLocalNumSpecies(id)
      enddo
!
      if (.not.isDataStorageExisting('NewL-ElectronDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'NewL-ElectronDensity',  &
                                DataSizeL,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('NewL-ElectronDensity',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'NewL-ElectronDensity',nr)
!           ----------------------------------------------------------
         enddo
         if (rad_derivative) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'L-ElectronDensityDer',&
                                   DataSizeL,ComplexType)
!           ----------------------------------------------------------
            call setDataStorage2Value('L-ElectronDensityDer',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'L-ElectronDensityDer',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
!
      if (.not.isDataStorageExisting('OldL-ElectronDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'OldL-ElectronDensity',  &
                                DataSizeL,ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('OldL-ElectronDensity',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'OldL-ElectronDensity',nr)
!           ----------------------------------------------------------
         enddo
      endif
!
      if (n_spin_pola == 2) then
         if (.not.isDataStorageExisting('NewL-MomentDensity')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'NewL-MomentDensity', &
                                   DataSizeL,ComplexType)
!           ----------------------------------------------------------
            call setDataStorage2Value('NewL-MomentDensity',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'NewL-MomentDensity',nr)
!              -------------------------------------------------------
            enddo
            if (rad_derivative) then
!              -------------------------------------------------------
               call createDataStorage(LocalNumAtoms,'L-MomentDensityDer', &
                                      DataSizeL,ComplexType)
!              -------------------------------------------------------
               call setDataStorage2Value('L-MomentDensityDer',CZERO)
!              -------------------------------------------------------
               do id = 1,LocalNumAtoms
                  nr = getNumRmesh(id)
!                 ----------------------------------------------------
                  call setDataStorageLDA(id,'L-MomentDensityDer',nr)
!                 ----------------------------------------------------
               enddo
            endif
         endif
      endif
!
      if (n_spin_pola == 2) then
         if (.not.isDataStorageExisting('OldL-MomentDensity')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'OldL-MomentDensity', &
                                DataSizeL,ComplexType) ! en extra space is used
                                                        ! storing ValenceVPCharge(2)
!           ----------------------------------------------------------
            call setDataStorage2Value('OldL-MomentDensity',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'OldL-MomentDensity',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
!
      if (.not.isMTFP) then
         if (.not.isDataStorageExisting('PseudoElectronDensity')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'PseudoElectronDensity',  &
                                   DataSize,RealType)  ! en extra space is used
                                                       ! storing ValenceVPCharge(1)
!           ----------------------------------------------------------
            call setDataStorage2Value('PseudoElectronDensity',ZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'PseudoElectronDensity',nr)
!              -------------------------------------------------------
            enddo
         endif
!
         if (n_spin_pola == 2) then
            if (.not.isDataStorageExisting('PseudoMomentDensity')) then
!              -------------------------------------------------------
               call createDataStorage(LocalNumAtoms,'PseudoMomentDensity', &
                                      DataSize,RealType) ! en extra space is used
                                                         ! storing ValenceVPCharge(2)
!              -------------------------------------------------------
               call setDataStorage2Value('PseudoMomentDensity',ZERO)
!              -------------------------------------------------------
               do id = 1,LocalNumAtoms
                  nr = getNumRmesh(id)
!                 ----------------------------------------------------
                  call setDataStorageLDA(id,'PseudoMomentDensity',nr)
!                 ----------------------------------------------------
               enddo
            endif
         endif
!
         if (.not.isDataStorageExisting('L-PseudoElectronDensity')) then
!           ---------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'L-PseudoElectronDensity',  &
                                   DataSizeL,ComplexType)  ! en extra space is used
                                                          ! storing ValenceVPCharge(1)
!           ----------------------------------------------------------
            call setDataStorage2Value('L-PseudoElectronDensity',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'L-PseudoElectronDensity',nr)
!              -------------------------------------------------------
            enddo
         endif
!
         if (n_spin_pola == 2) then
            if (.not.isDataStorageExisting('L-PseudoMomentDensity')) then
!              -------------------------------------------------------
               call createDataStorage(LocalNumAtoms,'L-PseudoMomentDensity', &
                                      DataSizeL,ComplexType) ! en extra space is used
                                                           ! storing ValenceVPCharge(2)
!              -------------------------------------------------------
               call setDataStorage2Value('L-PseudoMomentDensity',CZERO)
!              -------------------------------------------------------
               do id = 1,LocalNumAtoms
                  nr = getNumRmesh(id)
!                 ----------------------------------------------------
                  call setDataStorageLDA(id,'L-PseudoMomentDensity',nr)
!                 ----------------------------------------------------
               enddo
            endif
         endif
!
         if (.not.isDataStorageExisting('L-PseudoValElectronDensity')) then
!           ----------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'L-PseudoValElectronDensity',  &
                                   DataSizeL,ComplexType)  ! en extra space is used
                                                           ! storing ValenceVPCharge(1)
!           ----------------------------------------------------------
            call setDataStorage2Value('L-PseudoValElectronDensity',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'L-PseudoValElectronDensity',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
!
      if (.not.isDataStorageExisting('TildaElectronDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'TildaElectronDensity',  &
                                DataSize,RealType)  ! en extra space is used
                                                 ! storing ValenceVPCharge(1)
!        -------------------------------------------------------------
         call setDataStorage2Value('TildaElectronDensity',ZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'TildaElectronDensity',nr)
!           ----------------------------------------------------------
         enddo
      endif
!
      if (n_spin_pola == 2) then
         if (.not.isDataStorageExisting('TildaMomentDensity')) then
!           ---------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'TildaMomentDensity', &
                                   DataSize,RealType) ! en extra space is used
                                                      ! storing ValenceVPCharge(2)
!           ----------------------------------------------------------
            call setDataStorage2Value('TildaMomentDensity',ZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'TildaMomentDensity',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
!
      if (.not.isDataStorageExisting('L-TildaElectronDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'L-TildaElectronDensity',  &
                                DataSizeL,ComplexType)  ! en extra space is used
                                                 ! storing ValenceVPCharge(1)
!        -------------------------------------------------------------
         call setDataStorage2Value('L-TildaElectronDensity',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            nr = getNumRmesh(id)
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'L-TildaElectronDensity',nr)
!           ----------------------------------------------------------
         enddo
      endif
!
      if (n_spin_pola == 2) then
         if (.not.isDataStorageExisting('L-TildaMomentDensity')) then
!           ---------------------------------------------------------
            call createDataStorage(LocalNumAtoms,'L-TildaMomentDensity', &
                                   DataSizeL,ComplexType) ! en extra space is used
                                                   ! storing ValenceVPCharge(2)
!           ----------------------------------------------------------
            call setDataStorage2Value('L-TildaMomentDensity',CZERO)
!           ----------------------------------------------------------
            do id = 1,LocalNumAtoms
               nr = getNumRmesh(id)
!              -------------------------------------------------------
               call setDataStorageLDA(id,'L-TildaMomentDensity',nr)
!              -------------------------------------------------------
            enddo
         endif
      endif
!
      allocate(MM_table_line(GlobalNumAtoms))
      MM_size = 0
      do ig = 1, GlobalNumAtoms
         MM_table_line(ig) = MM_size
         MM_size = MM_size + getNumAlloyElements(ig)
      enddo
      jmax_glb = jmax_max
      call GlobalMaxInGroup(GroupID,jmax_glb)
      if ( .not.isDataStorageExisting('MultipoleMoments') ) then
         call createDataStorage('MultipoleMoments', &
                                jmax_glb*MM_size, ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('MultipoleMoments',CZERO)
!        -------------------------------------------------------------
      endif
!
      MultipoleMom => getDataStorage( 'MultipoleMoments',    &
                                      MM_size, jmax_glb, ComplexMark )
!
      do id = 1,LocalNumAtoms
         p_CDL => ChargeDensityList(id)
         nr = p_CDL%NumRPts
         r_ps = getRcutPseudo(id)
         Grid => getGrid(id)
         jmt  = Grid%jmt
         r_mesh => getRmesh(id)
         if (isMTFP) then
            r_ps = p_CDL%r_mesh(jmt)
            ir = jmt
         else 
            if ( r_ps <= ONE ) then ! If r_ps < 1.0, it is regarded as r_ps/Rmt
               r_ps = r_ps*p_CDL%r_mesh(jmt)
            endif
            ir = 1
!           ----------------------------------------------------------
            call hunt(nr,p_CDL%r_mesh,r_ps,ir)
!           ----------------------------------------------------------
         endif
         if ( print_level(id)>=0 ) then
            write(6,'(/,a,2(1x,f12.8))') "ChargeDensity:: Input:  R_pseudo/R_mt, R_pseudo = ", &
                  r_ps/r_mesh(jmt), r_ps
         endif
         p_CDL%NumRPts_Pseudo = ir
         p_CDL%RcutPseudo = p_CDL%r_mesh(ir)
         if ( print_level(id)>=0 ) then
            write(6,'(a,i5,2(1x,f12.8))')"ChargeDensity:: Actual: ir-pseudo, R_pseudo/R_mt, R_pseudo = ", &
                     p_CDL%NumRPts_Pseudo, p_CDL%RcutPseudo/p_CDL%r_mesh(jmt), p_CDL%RcutPseudo
         endif
!
         jmax = p_CDL%jmax
         lmax = p_CDL%lmax
         p_CDL%rhoL_Total => getDataStorage( id,'NewL-ElectronDensity', &
                                             nr, jmax, getLocalNumSpecies(id), ComplexMark )
         if (n_spin_pola == 2) then
            p_CDL%momL_Total => getDataStorage( id,'NewL-MomentDensity', &
                                                nr, jmax, getLocalNumSpecies(id), ComplexMark )
         endif
!
         if (rad_derivative) then
            p_CDL%der_rhoL_Total => getDataStorage( id,'L-ElectronDensityDer', &
                                                    nr, jmax, getLocalNumSpecies(id), ComplexMark )
            if (n_spin_pola == 2) then
               p_CDL%der_momL_Total => getDataStorage( id,'L-MomentDensityDer', &
                                                       nr, jmax, getLocalNumSpecies(id), ComplexMark )
            endif
         endif
!
         p_CDL%rhoL_TotalOld => getDataStorage( id,'OldL-ElectronDensity', &
                                             nr, jmax, getLocalNumSpecies(id), ComplexMark )
         do ia = 1, p_CDL%NumSpecies
            p_CDL%rhoL_TotalOld(1:nr,1,ia) = p_CDL%rhoSph_TotalOld(1:nr,ia)/Y0
         enddo
         if (n_spin_pola == 2) then
            p_CDL%momL_TotalOld => getDataStorage( id,'OldL-MomentDensity', &
                                                nr, jmax, getLocalNumSpecies(id), ComplexMark )
         endif
!
         if (.not.isMTFP) then
            p_CDL%rhoSph_Pseudo => getDataStorage( id,'PseudoElectronDensity', &
                                                   nr, getLocalNumSpecies(id), RealMark )
            if (n_spin_pola == 2) then
               p_CDL%momSph_Pseudo => getDataStorage( id,'PseudoMomentDensity', &
                                                      nr, getLocalNumSpecies(id), RealMark )
            endif
!
            p_CDL%rhoL_Pseudo => getDataStorage( id,'L-PseudoElectronDensity', &
                                                 nr, jmax, getLocalNumSpecies(id), ComplexMark )
            if (n_spin_pola == 2) then
               p_CDL%momL_Pseudo => getDataStorage( id,'L-PseudoMomentDensity', &
                                                    nr, jmax, getLocalNumSpecies(id), ComplexMark )
            endif
!
            p_CDL%rhoL_ValPseudo => getDataStorage( id,'L-PseudoValElectronDensity', &
                                                    nr, jmax, getLocalNumSpecies(id), ComplexMark )
         endif
!
         p_CDL%rhoSph_Tilda => getDataStorage( id,'TildaElectronDensity', &
                                               nr, getLocalNumSpecies(id), RealMark )
         if (n_spin_pola == 2) then
            p_CDL%momSph_Tilda => getDataStorage( id,'TildaMomentDensity', &
                                                  nr, getLocalNumSpecies(id), RealMark )
         endif
!
         p_CDL%rhoL_Tilda => getDataStorage( id,'L-TildaElectronDensity', &
                                             nr, jmax, getLocalNumSpecies(id), ComplexMark )
         if (n_spin_pola == 2) then
            p_CDL%momL_Tilda => getDataStorage( id,'L-TildaMomentDensity', &
                                                 nr, jmax, getLocalNumSpecies(id), ComplexMark )
         endif
!
         p_CDL%rhoL_Valence => getValenceElectronDensity(id)
         if (n_spin_pola == 2) then
            p_CDL%momL_Valence => getValenceMomentDensity(id)
         endif
!
         allocate( p_CDL%multipole_mom(jmax,getLocalNumSpecies(id)) )
!
         p_CDL%ChargeCompFlag = 0
         p_CDL%ChargeCompFlag(1) = 1
!
      enddo
!
      nullify(Grid)
      deallocate( DataSizeL )
!
   endif   ! non-spherical charge
!
   deallocate( DataSize )
!
   nullify(p_CDL)
!
   allocate( sqrt_r(0:NumRPts_max), den_in(0:NumRPts_max), den_out(0:NumRPts_max) )
   allocate( den_r(0:NumRPts_max),den_r1(0:NumRPts_max) )
!  -------------------------------------------------------------------
   call initIntegerFactors(lmax_max)
!  -------------------------------------------------------------------
!
   Initialized = .true.
!
   end subroutine initChargeDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endChargeDensity()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
!
   implicit none
!
   integer (kind=IntKind) :: id
   type (ChargeDensityStruct), pointer :: p_CDL
!
   do id = 1, LocalNumAtoms
!
      p_CDL => ChargeDensityList(id)
      nullify( p_CDL%r_mesh )
      nullify( p_CDL%rhoSph_Total, p_CDL%momSph_Total )
      nullify( p_CDL%rhoSph_TotalOld, p_CDL%momSph_TotalOld )
!
      nullify( p_CDL%rhoL_Total, p_CDL%momL_Total )
      nullify( p_CDL%rhoL_TotalOld, p_CDL%momL_TotalOld )
!
      nullify( p_CDL%rhoSph_Pseudo, p_CDL%momSph_Pseudo )
      nullify( p_CDL%rhoL_Pseudo, p_CDL%momL_Pseudo )
!
      nullify( p_CDL%rhoSph_Tilda, p_CDL%momSph_Tilda )
      nullify( p_CDL%rhoL_Tilda, p_CDL%momL_Tilda )
!
      nullify( p_CDL%rhoL_Valence, p_CDL%momL_Valence )
!
      nullify( p_CDL%rhoL_ValPseudo )
!
      if (rad_derivative) then
         nullify( p_CDL%der_rhoSph_Total, p_CDL%der_momSph_Total )
         nullify( p_CDL%der_rhoL_Total, p_CDL%der_momL_Total )
      endif
!
      if ( associated(p_CDL%multipole_mom) ) then
         deallocate( p_CDL%multipole_mom )
      endif
!
      deallocate( p_CDL%ChargeCompFlag )
!
   enddo
   nullify( p_CDL, MultipoleMom )
!
   if ( .not.isSphericalCharge ) then
      deallocate( MM_table_line )
   endif
!
   deallocate( ChargeDensityList )
   deallocate( ws_den, ws_der_den, ws_ylm, ws_grad_ylm )
   deallocate( sqrt_r, den_in, den_out, den_r, den_r1 )
!  -------------------------------------------------------------------
   call endIntegerFactors()
!  -------------------------------------------------------------------
   Initialized = .false.
!
   end subroutine endChargeDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphChargeDensity(densityType, id, ia, p_grad_den) result(p_den)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: densityType
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   integer (kind=IntKind) :: nr
!
   real (kind=RealKind), pointer :: p_den(:)
   real (kind=RealKind), pointer, optional :: p_grad_den(:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   if ( .not.Initialized ) then
      call ErrorHandler("getChrgDen_Sph"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   else if ( id<1 .or. id> LocalNumAtoms ) then
      call ErrorHandler("getSphChargeDensity"," Invalid local atom index",id)
   else if (isMTFP .and. nocaseCompare(densityType,"Pseudo") ) then
      call ErrorHandler("getSphChargeDensity"," Invalid potential type",'Pseudo')
   endif
!
   p_CDL => ChargeDensityList(id)
   nr = p_CDL%NumRPts+1
!
   if ( nocaseCompare(densityType,"Pseudo") ) then
      p_den => p_CDL%rhoSph_Pseudo(1:nr,ia)
   else if ( nocaseCompare(densityType,"Tilda") ) then
      p_den => p_CDL%rhoSph_Tilda(1:nr,ia)
   else if ( nocaseCompare(densityType,"TotalNew") ) then
      p_den => p_CDL%rhoSph_Total(1:nr,ia)
      if ( present(p_grad_den) ) then
         p_grad_den => p_CDL%der_rhoSph_Total(1:nr,ia)
      endif
   else if ( nocaseCompare(densityType,"TotalOld") ) then
      p_den => p_CDL%rhoSph_TotalOld(1:nr,ia)
   else
      call ErrorHandler("getChrgDen_Sph","Undefined charge type")
   endif
!
   end function getSphChargeDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getChrgDen_L(densityType, id, ia, p_grad_den) result(p_den)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: densityType
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   integer (kind=IntKind) :: nr, jmax
!
   complex (kind=CmplxKind), pointer :: p_den(:,:)
   complex (kind=CmplxKind), pointer, optional :: p_grad_den(:,:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   if ( .not.Initialized ) then
      call ErrorHandler("getChrgDen_L"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   else if ( id<1 .or. id> LocalNumAtoms ) then
      call ErrorHandler("getChrgDen_L"," Invalid local atom index",id)
   else if (isMTFP .and. nocaseCompare(densityType,"Pseudo") ) then
      call ErrorHandler("getChrgDen_L"," Invalid potential type",'Pseudo')
   else if (isMTFP .and. densityType(1:5)=="ValPs" ) then
      call ErrorHandler("getChrgDen_L"," Invalid potential type",'ValPs')
   endif
!
   p_CDL => ChargeDensityList(id)
   jmax = p_CDL%jmax
   nr   = p_CDL%NumRPts
!
   if ( nocaseCompare(densityType,"Pseudo") ) then
         p_den => p_CDL%rhoL_Pseudo(1:nr,1:jmax,ia)
   else if ( nocaseCompare(densityType,"Tilda") ) then
         p_den => p_CDL%rhoL_Tilda(1:nr,1:jmax,ia)
   else if (nocaseCompare(densityType,"ValPs") ) then
         p_den => p_CDL%rhoL_ValPseudo(1:nr,1:jmax,ia)
   else if ( nocaseCompare(densityType,"Valence") ) then
         p_den => p_CDL%rhoL_Valence(1:nr,1:jmax,ia)
   else if ( nocaseCompare(densityType,"TotalNew") ) then
         p_den => p_CDL%rhoL_Total(1:nr,1:jmax,ia)
      if ( present(p_grad_den) ) then
         p_grad_den => p_CDL%der_rhoL_Total(1:nr,1:jmax,ia)
      endif
   else if ( nocaseCompare(densityType,"TotalOld") ) then
         p_den => p_CDL%rhoL_TotalOld(1:nr,1:jmax,ia)
   else
      call ErrorHandler("getChrgDen_L","Undefined charge type")
   endif
!
   end function getChrgDen_L
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getChrgDen_Lj(densityType, id, ia, jl, p_grad_den) result(p_den)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: densityType
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in) :: jl
!
   integer (kind=IntKind) :: nr, jmax
!
   complex (kind=CmplxKind), pointer :: p_den(:)
   complex (kind=CmplxKind), pointer, optional :: p_grad_den(:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   if ( .not.Initialized ) then
      call ErrorHandler("getChrgDen_Lj"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   else if ( id<1 .or. id> LocalNumAtoms ) then
      call ErrorHandler("getChrgDen_Lj"," Invalid local atom index",id)
   else if (isMTFP .and. nocaseCompare(densityType,"Pseudo") ) then
      call ErrorHandler("getChrgDen_Lj"," Invalid potential type",'Pseudo')
   else if (isMTFP .and. nocaseCompare(densityType,"ValPs") ) then
      call ErrorHandler("getChrgDen_Lj"," Invalid potential type",'ValPs')
   endif
!
   p_CDL => ChargeDensityList(id)
   jmax = p_CDL%jmax
   nr   = p_CDL%NumRPts
!
   if ( jl<1 .or. jl>jmax ) then
      call ErrorHandler("getChrgDen_Lj"," Invalid jl ", jl)
   endif
!
   if ( nocaseCompare(densityType,"Pseudo") ) then
         p_den => p_CDL%rhoL_Pseudo(1:nr,jl,ia)
   else if ( nocaseCompare(densityType,"Tilda") ) then
         p_den => p_CDL%rhoL_Tilda(1:nr,jl,ia)
   else if ( nocaseCompare(densityType,"ValPs") ) then
         p_den => p_CDL%rhoL_ValPseudo(1:nr,jl,ia)
   else if ( nocaseCompare(densityType,"Valence") ) then
         p_den => p_CDL%rhoL_Valence(1:nr,jl,ia)
   else if ( nocaseCompare(densityType,"TotalNew") ) then
         p_den => p_CDL%rhoL_Total(1:nr,jl,ia)
      if ( present(p_grad_den) ) then
         p_grad_den => p_CDL%der_rhoL_Total(1:nr,jl,ia)
      endif
   else if ( nocaseCompare(densityType,"TotalOld") ) then
         p_den => p_CDL%rhoL_TotalOld(1:nr,jl,ia)
   else
      call ErrorHandler("getChrgDen_Lj","Undefined charge type")
   endif
!
   end function getChrgDen_Lj
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getChargeDensityAtPoint(densityType, id, ia, posi, jmax_in, n_mult, grad) result(rho)
!  ===================================================================
   use MathParamModule, only : ZERO, TWO
!
   use SphericalHarmonicsModule, only : calYlm
!
   use InterpolationModule, only : PolyInterp
!
   use PolyhedraModule, only: getPointLocationFlag
!
   implicit none
!
   character(len=*), intent(in) :: densityType
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in), optional :: jmax_in
   integer (kind=IntKind), intent(in), optional :: n_mult
!
   real (kind=RealKind), intent(in) :: posi(3)
   real (kind=RealKind), intent(out), optional :: grad(3)
!
   integer (kind=IntKind) :: VP_point
   integer (kind=IntKind) :: jl, ir, is, irp, l, kl, i
   integer (kind=IntKind) :: iend, kmax, nr, jmax
!
   real (kind=RealKind) :: rho, err, r, fact, er(3)
!
   complex (kind=CmplxKind) :: rho_in, der_rho_in
   complex (kind=CmplxKind), pointer :: ylm(:)
   complex (kind=CmplxKind), pointer :: den_l(:,:)
   complex (kind=CmplxKind), pointer :: grad_ylm(:,:)
   complex (kind=CmplxKind), pointer :: der_den_l(:,:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   logical :: takeGradient = .false.
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
   if (.not.Initialized) then
      call ErrorHandler('getChargeDensityAtPoint',                    &
                        'Carge Density Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getChargeDensityAtPoint','Invalid atom index',id)
   else if (isMTFP .and. nocaseCompare(densityType,"Pseudo") ) then
      call ErrorHandler("getChrgDen_Lj"," Invalid potential type",'Pseudo')
   else if (isMTFP .and. nocaseCompare(densityType,"ValPs") ) then
      call ErrorHandler("getChrgDen_Lj"," Invalid potential type",'ValPs')
   else if (.not.nocaseCompare(densityType,"TotalNew") .and. present(grad)) then
      call ErrorHandler('getChrgDen_Lj','The gradient of the quantity is not implemented',densityType)
   endif
!
   p_CDL => ChargeDensityList(id)
!
   if ( present(jmax_in) ) then
      if (jmax_in < 1 .or. jmax_in > p_CDL%jmax) then
         call ErrorHandler('getChargeDensityAtPoint','Invalid jmax value',jmax_in)
      endif
      jmax = jmax_in
   else
      jmax = p_CDL%jmax
   endif
!
   if (nocaseCompare(densityType,"TotalNew") .and. present(grad) .and.  rad_derivative) then
      takeGradient = .true.
   else
      takeGradient = .false.
   endif
!
   nr = p_CDL%NumRPts
   r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
   if ( r > p_CDL%r_mesh(nr)+TEN2m8 ) then
      rho = ZERO
      return
   endif
   if (r > TEN2m8) then
      er = posi/r  ! unit r-vector
   else
      er(1) = ZERO; er(2) = ZERO; er(3) = ONE
   endif
!
   kmax = (p_CDL%lmax+1)*(p_CDL%lmax+1)
   ylm => ws_ylm(1:kmax)
   den_l => ws_den(1:n_inter,1:p_CDL%jmax)
!
   if (takeGradient) then
      grad_ylm => ws_grad_ylm
      der_den_l => ws_der_den
!     ----------------------------------------------------------------
      call calYlm(posi,p_CDL%lmax,ylm,grad_ylm)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call calYlm(posi,p_CDL%lmax,ylm)
!     ----------------------------------------------------------------
   endif
!
   iend = p_CDL%NumRPts
!  -------------------------------------------------------------------
   call hunt(iend,p_CDL%r_mesh(1:iend),r,ir)
!  -------------------------------------------------------------------
   if (ir > iend-(n_inter-1)/2) then
      irp=iend-n_inter+1
   else if (2*ir+1 > n_inter) then
      irp=ir-(n_inter-1)/2
   else
      irp=1
   endif
!
   den_l = CZERO
   if (takeGradient) then
      der_den_l = CZERO
   endif
!
   if ( present(n_mult) ) then
      if ( n_mult<1 ) then 
         call ErrorHandler("getChargeDensityAtPoint",'Invalid n_mult',n_mult)
      endif
      VP_point = getPointLocationFlag(id, posi(1), posi(2), posi(3))
      if ( nocaseCompare(densityType,"Pseudo") ) then
         if ( VP_point == 1 .or. VP_point == 0 ) then
            do jl = 1,p_CDL%jmax
               den_l(1:n_inter,jl) = p_CDL%rhoL_Pseudo(irp:irp+n_inter-1,jl,ia)
            enddo
         else  ! ouside Voronoi polyhedra, overlaping semicore charges
            rho = ZERO
            return
         endif
      else if ( nocaseCompare(densityType,"ValPs") ) then
         if ( VP_point >= 0 ) then
            do jl = 1,p_CDL%jmax
               den_l(1:n_inter,jl) = p_CDL%rhoL_ValPseudo(irp:irp+n_inter-1,jl,ia)
            enddo
         else
            rho = ZERO
            return
         endif
      else if ( nocaseCompare(densityType,"Tilda") ) then
         do jl = 1,p_CDL%jmax
            den_l(1:n_inter,jl) = p_CDL%rhoL_Tilda(irp:irp+n_inter-1,jl,ia)
         enddo
      else if ( nocaseCompare(densityType,"Valence") ) then
         if ( VP_point >= 0 ) then
            do jl = 1,p_CDL%jmax
               den_l(1:n_inter,jl) = p_CDL%rhoL_Valence(irp:irp+n_inter-1,jl,ia)
            enddo
         else
            rho = ZERO
            return
         endif
      else if ( nocaseCompare(densityType,"TotalNew") ) then
         if ( VP_point == 1 .or. VP_point == 0 ) then
            do jl = 1,p_CDL%jmax
               den_l(1:n_inter,jl) = p_CDL%rhoL_Total(irp:irp+n_inter-1,jl,ia)
            enddo
         else  ! ouside Voronoi polyhedra, overlaping semicore charges
!           rho = ZERO
!           return
            fact = real(n_mult,kind=RealKind)/Y0
            do is = 1,n_spin_pola
               do ir = 1, n_inter
                  den_l(ir,1) = den_l(ir,1) +                                 &
                                cmplx(fact*(p_CDL%rho_Core(irp+ir-1,is,ia) +  &
                                     p_CDL%rho_SemiCore(irp+ir-1,is,ia)),ZERO,&
                                     kind=CmplxKind)
               enddo
            enddo
         endif
         if (VP_point == 0) then ! For points on the cell bounday, since den_l contains
                                 ! the core density, it needs to be corrected
            fact = real(n_mult-1,kind=RealKind)/Y0
            do is = 1,n_spin_pola
               do ir = 1, n_inter
                  den_l(ir,1) = den_l(ir,1) +                                 &
                                cmplx(fact*(p_CDL%rho_Core(irp+ir-1,is,ia) +  &
                                     p_CDL%rho_SemiCore(irp+ir-1,is,ia)),ZERO,&
                                     kind=CmplxKind)
               enddo
            enddo
         endif
         if ( takeGradient ) then
            if ( VP_point == 1 .or. VP_point == 0 ) then
               do jl = 1,p_CDL%jmax
                  der_den_l(1:n_inter,jl) = p_CDL%der_rhoL_Total(irp:irp+n_inter-1,jl,ia)
               enddo
            else  ! ouside Voronoi polyhedra, overlaping semicore charges
               fact = real(n_mult,kind=RealKind)/Y0
               do is = 1,n_spin_pola
                  do ir = 1, n_inter
                     der_den_l(ir,1) = der_den_l(ir,1) +                                   &
                                       cmplx(fact*(p_CDL%der_rho_Core(irp+ir-1,is,ia) +    &
                                             p_CDL%der_rho_SemiCore(irp+ir-1,is,ia)),ZERO, &
                                             kind=CmplxKind)
                  enddo
               enddo
            endif
            if (VP_point == 0) then ! For points on the cell bounday, since den_l contains
                                    ! the core density, it needs to be corrected
               fact = real(n_mult-1,kind=RealKind)/Y0
               do is = 1,n_spin_pola
                  do ir = 1, n_inter
                     der_den_l(ir,1) = der_den_l(ir,1) +                                   &
                                       cmplx(fact*(p_CDL%der_rho_Core(irp+ir-1,is,ia) +    &
                                             p_CDL%der_rho_SemiCore(irp+ir-1,is,ia)),ZERO, &
                                             kind=CmplxKind)
                  enddo
               enddo
            endif
         endif
      else if ( nocaseCompare(densityType,"TotalOld") ) then
         if ( VP_point == 1 .or. VP_point == 0 ) then
            do jl = 1,p_CDL%jmax
               den_l(1:n_inter,jl) = p_CDL%rhoL_TotalOld(irp:irp+n_inter-1,jl,ia)
            enddo
         else ! include core density outside the atomic cell
            fact = real(n_mult,kind=RealKind)/Y0
            jmax=1
            do is = 1,n_spin_pola
               do ir = 1, n_inter
                  den_l(ir,1) = den_l(ir,1) +                                 &
                                cmplx(fact*(p_CDL%rho_Core(irp+ir-1,is,ia) +  &
                                     p_CDL%rho_SemiCore(irp+ir-1,is,ia)),ZERO,&
                                     kind=CmplxKind)
               enddo
            enddo
         endif
         if (VP_point == 0) then ! For points on the cell bounday, since den_l contains
                                 ! the core density, it needs to be corrected
            fact = real(n_mult-1,kind=RealKind)/Y0
            do is = 1,n_spin_pola
               do ir = 1, n_inter
                  den_l(ir,1) = den_l(ir,1) +                                 &
                                cmplx(fact*(p_CDL%rho_Core(irp+ir-1,is,ia) +  &
                                     p_CDL%rho_SemiCore(irp+ir-1,is,ia)),ZERO,&
                                     kind=CmplxKind)
               enddo
            enddo
         endif
      else
         call ErrorHandler("getChargeDensityAtPoint","Undefined charge type")
      endif
   else
      if ( nocaseCompare(densityType,"Pseudo") ) then
         do jl = 1,p_CDL%jmax
            den_l(1:n_inter,jl) = p_CDL%rhoL_Pseudo(irp:irp+n_inter-1,jl,ia)
         enddo
      else if ( nocaseCompare(densityType,"ValPs") ) then
         do jl = 1,p_CDL%jmax
            den_l(1:n_inter,jl) = p_CDL%rhoL_ValPseudo(irp:irp+n_inter-1,jl,ia)
         enddo
      else if ( nocaseCompare(densityType,"Tilda") ) then
         do jl = 1,p_CDL%jmax
            den_l(1:n_inter,jl) = p_CDL%rhoL_Tilda(irp:irp+n_inter-1,jl,ia)
         enddo
      else if ( nocaseCompare(densityType,"Valence") ) then
         do jl = 1,p_CDL%jmax
            den_l(1:n_inter,jl) = p_CDL%rhoL_Valence(irp:irp+n_inter-1,jl,ia)
         enddo
      else if ( nocaseCompare(densityType,"TotalNew") ) then
         do jl = 1,p_CDL%jmax
            den_l(1:n_inter,jl) = p_CDL%rhoL_Total(irp:irp+n_inter-1,jl,ia)
         enddo
         if ( takeGradient ) then
            do jl = 1,p_CDL%jmax
               der_den_l(1:n_inter,jl) = p_CDL%der_rhoL_Total(irp:irp+n_inter-1,jl,ia)
            enddo
         endif
      else if ( nocaseCompare(densityType,"TotalOld") ) then
         do jl = 1,p_CDL%jmax
            den_l(1:n_inter,jl) = p_CDL%rhoL_TotalOld(irp:irp+n_inter-1,jl,ia)
         enddo
      else
         call ErrorHandler("getChargeDensityAtPoint","Undefined charge type")
      endif
   endif
!
   rho = ZERO
   if (takeGradient) then
      grad = ZERO
   endif
   do jl = jmax,1,-1
      l = lofj(jl)
!     ----------------------------------------------------------------
      call PolyInterp(n_inter, p_CDL%r_mesh(irp:irp+n_inter-1), &
                      den_l(1:n_inter,jl), r, rho_in, err)
!     ----------------------------------------------------------------
      if ( takeGradient ) then
!        -------------------------------------------------------------
         call PolyInterp(n_inter, p_CDL%r_mesh(irp:irp+n_inter-1), &
                         der_den_l(1:n_inter,jl), r, der_rho_in, err)
!        -------------------------------------------------------------
      endif
      kl = (l+1)*(l+1)-l+mofj(jl)
      if (mofj(jl) == 0) then
         fact = ONE
      else
         fact = TWO
      endif
      rho = rho + fact*real(rho_in*ylm(kl),RealKind)
      if ( takeGradient ) then
         do i = 1, 3
            grad(i) = grad(i) + fact*real(der_rho_in*ylm(kl)*er(i)+  &
                                          rho_in*grad_ylm(kl,i),RealKind)
         enddo
      endif
   enddo
!
   end function getChargeDensityAtPoint
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphMomentDensity(densityType, id, ia, p_grad_mom) result(p_mom)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: densityType
!
   integer (kind=IntKind) :: id, ia
!
   integer (kind=IntKind) :: nr
!
   real (kind=RealKind), pointer :: p_mom(:)
   real (kind=RealKind), pointer, optional :: p_grad_mom(:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   if ( .not.Initialized ) then
      call ErrorHandler("getMomDen_Sph"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   else if ( id<1 .or. id> LocalNumAtoms ) then
      call ErrorHandler("getMomDen_Sph"," Invalid local atom index",id)
   else if (isMTFP .and. nocaseCompare(densityType,"Pseudo") ) then
      call ErrorHandler("getSphMomentDensity"," Invalid potential type",'Pseudo')
   endif
!
   p_CDL => ChargeDensityList(id)
   nr = p_CDL%NumRPts+1
!
   if ( nocaseCompare(densityType,"Pseudo") ) then
      p_mom => p_CDL%momSph_Pseudo(1:nr,ia)
   else if ( nocaseCompare(densityType,"Tilda") ) then
      p_mom => p_CDL%momSph_Tilda(1:nr,ia)
   else if ( nocaseCompare(densityType,"TotalNew") ) then
      p_mom => p_CDL%momSph_Total(1:nr,ia)
      if ( present(p_grad_mom) ) then
         p_grad_mom => p_CDL%der_momSph_Total(1:nr,ia)
      endif
   else if ( nocaseCompare(densityType,"TotalOld") ) then
      p_mom => p_CDL%momSph_TotalOld(1:nr,ia)
   else
      call ErrorHandler("getMomDen_Sph","Undefined moment type")
   endif
!
   end function getSphMomentDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMomDen_L( momType, id, ia, p_grad_mom) result(p_mom)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: momType
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   integer (kind=IntKind) :: nr, jmax
!
   complex (kind=CmplxKind), pointer :: p_mom(:,:)
   complex (kind=CmplxKind), pointer, optional :: p_grad_mom(:,:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   if ( .not.Initialized ) then
      call ErrorHandler("getMomDen_L"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   else if ( id<1 .or. id> LocalNumAtoms ) then
      call ErrorHandler("getMomDen_L"," Invalid local atom index",id)
   else if (isMTFP .and. nocaseCompare(momType,"Pseudo") ) then
      call ErrorHandler("getMomDen_L"," Invalid potential type",'Pseudo')
   endif
!
   p_CDL => ChargeDensityList(id)
   jmax = p_CDL%jmax
   nr   = p_CDL%NumRPts
!
   if ( nocaseCompare(momType,"Pseudo") ) then
      p_mom => p_CDL%momL_Pseudo(1:nr,1:jmax,ia)
   else if ( nocaseCompare(momType,"Tilda") ) then
      p_mom => p_CDL%momL_Tilda(1:nr,1:jmax,ia)
   else if ( nocaseCompare(momType,"Valence") ) then
      if (n_spin_cant == 2) then
         call WarningHandler('getMomDen_L','Only show z-component of the moment density')
         p_mom => p_CDL%momL_Valence(1:nr,1:jmax,3,ia)
      else
         p_mom => p_CDL%momL_Valence(1:nr,1:jmax,1,ia)
      endif
   else if ( nocaseCompare(momType,"TotalNew") ) then
      p_mom => p_CDL%momL_Total(1:nr,1:jmax,ia)
      if ( present(p_grad_mom) ) then
         p_grad_mom => p_CDL%der_momL_Total(1:nr,1:jmax,ia)
      endif
   else if ( nocaseCompare(momType,"TotalOld") ) then
      p_mom => p_CDL%momL_TotalOld(1:nr,1:jmax,ia)
   else
      call ErrorHandler("getMomDen_L","Undefined moment type")
   endif
!
   end function getMomDen_L
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMomDen_Lj(densityType, id, ia, jl, p_grad_mom) result(p_mom)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: densityType
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in) :: jl
!
   integer (kind=IntKind) :: nr, jmax
!
   complex (kind=CmplxKind), pointer :: p_mom(:)
   complex (kind=CmplxKind), pointer, optional :: p_grad_mom(:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   if ( .not.Initialized ) then
      call ErrorHandler("getMomDen_Lj"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   else if ( id<1 .or. id> LocalNumAtoms ) then
      call ErrorHandler("getMomDen_Lj"," Invalid local atom index",id)
   else if (isMTFP .and. nocaseCompare(densityType,"Pseudo") ) then
      call ErrorHandler("getMomDen_Lj"," Invalid potential type",'Pseudo')
   endif
!
   p_CDL => ChargeDensityList(id)
   jmax = p_CDL%jmax
   nr   = p_CDL%NumRPts
!
   if ( jl<1 .or. jl>jmax ) then
      call ErrorHandler("getMomDen_Lj"," Invalid jl ", jl)
   endif
!
   if ( nocaseCompare(densityType,"Pseudo") ) then
      p_mom => p_CDL%momL_Pseudo(1:nr,jl,ia)
   else if ( nocaseCompare(densityType,"Tilda") ) then
      p_mom => p_CDL%momL_Tilda(1:nr,jl,ia)
   else if ( nocaseCompare(densityType,"Valence") ) then
      if (n_spin_cant == 2) then
         call WarningHandler('getMomDen_L','Only show z-component of the moment density')
         p_mom => p_CDL%momL_Valence(1:nr,jl,3,ia)
      else
         p_mom => p_CDL%momL_Valence(1:nr,jl,1,ia)
      endif
   else if ( nocaseCompare(densityType,"TotalNew") ) then
      p_mom => p_CDL%momL_Total(1:nr,jl,ia)
      if ( present(p_grad_mom) ) then
         p_grad_mom => p_CDL%der_momL_Total(1:nr,jl,ia)
      endif
   else if ( nocaseCompare(densityType,"TotalOld") ) then
      p_mom => p_CDL%momL_TotalOld(1:nr,jl,ia)
   else
      call ErrorHandler("getMomDen_Lj","Undefined moment type")
   endif
!
   end function getMomDen_Lj
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMomentDensityAtPoint(momentType, id, ia, posi, jmax_in, n_mult, grad) result(mom)
!  ===================================================================
   use MathParamModule, only : ZERO, TWO
!
   use SphericalHarmonicsModule, only : calYlm
!
   use InterpolationModule, only : PolyInterp
!
   use PolyhedraModule, only: getPointLocationFlag
!
   implicit none
!
   character(len=*), intent(in) :: momentType
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in), optional :: jmax_in
   integer (kind=IntKind), intent(in), optional :: n_mult
!
   integer (kind=IntKind) :: VP_point
   integer (kind=IntKind) :: jl, ir, irp, l, kl, iend, i
   integer (kind=IntKind) :: kmax, nr, is, isig, jmax
!
   real (kind=RealKind), intent(in) :: posi(3)
   real (kind=RealKind), intent(out), optional :: grad(3)
   real (kind=RealKind) :: mom, err, r, fact, er(3)
!
   complex (kind=CmplxKind) :: mom_in, der_mom_in
   complex (kind=CmplxKind), pointer :: ylm(:), grad_ylm(:,:)
   complex (kind=CmplxKind), pointer :: mom_l(:,:), der_mom_l(:,:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   logical :: takeGradient = .false.
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
   if (.not.Initialized) then
      call ErrorHandler('getMomDen_r',                    &
                        'Carge Density Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getMomDen_r','Invalid atom index',id)
   else if (isMTFP .and. nocaseCompare(momentType,"Pseudo") ) then
      call ErrorHandler("getMomDen_r"," Invalid potential type",'Pseudo')
   else if (.not.nocaseCompare(momentType,"TotalNew") .and. present(grad)) then
      call ErrorHandler('getMomDen_Lj','The gradient of the quantity is not implemented',momentType)
   endif
!
   p_CDL => ChargeDensityList(id)
!
   if ( present(jmax_in) ) then
      if (jmax_in < 1 .or. jmax_in > p_CDL%jmax) then
         call ErrorHandler('getMomDen_r','Invalid jmax value',jmax_in)
      endif
      jmax = jmax_in
   else
      jmax = p_CDL%jmax
   endif
!
   if (nocaseCompare(momentType,"TotalNew") .and. present(grad) .and.  rad_derivative) then
      takeGradient = .true.
   else
      takeGradient = .false.
   endif
!
   nr = p_CDL%NumRPts
   r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
   if (r > p_CDL%r_mesh(nr)+TEN2m6 ) then
      mom = ZERO
      return
   endif
   if (r > TEN2m8) then
      er = posi/r  ! unit r-vector
   else
      er(1) = ZERO; er(2) = ZERO; er(3) = ONE
   endif
!
   kmax = (p_CDL%lmax+1)*(p_CDL%lmax+1)
   ylm => ws_ylm(1:kmax)
   mom_l => ws_den(1:n_inter,1:p_CDL%jmax)
   if (takeGradient) then
      grad_ylm => ws_grad_ylm
      der_mom_l => ws_der_den
!     ----------------------------------------------------------------
      call calYlm(posi,p_CDL%lmax,ylm,grad_ylm)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call calYlm(posi,p_CDL%lmax,ylm)
!     ----------------------------------------------------------------
   endif
!
   iend = p_CDL%NumRPts
!  -------------------------------------------------------------------
   call hunt(iend,p_CDL%r_mesh(1:iend),r,ir)
!  -------------------------------------------------------------------
   if (ir > iend-(n_inter-1)/2) then
      irp=iend-n_inter+1
   else if (2*ir+1 > n_inter) then
      irp=ir-(n_inter-1)/2
   else
      irp=1
   endif
!
   mom_l = CZERO
   if ( present(n_mult) ) then
      if ( n_mult<1 ) then
         call ErrorHandler("getMomentDensityAtPoint",'Invalid n_mult',n_mult)
      endif
      VP_point = getPointLocationFlag(id, posi(1), posi(2), posi(3))
      if ( nocaseCompare(momentType,"TotalNew") ) then
         if ( VP_point == 1 .or. VP_point == 0) then
            do jl = 1,jmax
               mom_l(1:n_inter,jl) = p_CDL%momL_Total(irp:irp+n_inter-1,jl,ia)
            enddo
            if (takeGradient) then
               do jl = 1,jmax
                  der_mom_l(1:n_inter,jl) = p_CDL%der_momL_Total(irp:irp+n_inter-1,jl,ia)
               enddo
            endif
            if (VP_point == 0) then ! For points on the cell bounday, since mom_l contains
                                    ! the core density, it needs to be corrected
               fact = real(n_mult-1,kind=RealKind)/Y0
               do is = 1,n_spin_pola
                  isig = 3-2*is
                  do ir = 1, n_inter
                     mom_l(ir,1) = mom_l(ir,1) +                                       &
                                   isig*cmplx(fact*(p_CDL%rho_Core(irp+ir-1,is,ia) +   &
                                              p_CDL%rho_SemiCore(irp+ir-1,is,ia)),ZERO,&
                                              kind=CmplxKind)
                  enddo
                  if (takeGradient) then
                     do ir = 1, n_inter
                        der_mom_l(ir,1) = der_mom_l(ir,1) +                            &
                                          isig*cmplx(fact*(p_CDL%der_rho_Core(irp+ir-1,is,ia) +    &
                                                           p_CDL%der_rho_SemiCore(irp+ir-1,is,ia)),&
                                                     ZERO, kind=CmplxKind)
                     enddo
                  endif
               enddo
            endif
         else
            fact = real(n_mult,kind=RealKind)/Y0
            do is = 1,n_spin_pola
               isig = 3-2*is
               do ir = 1, n_inter
                  mom_l(ir,1) = mom_l(ir,1) + isig*fact*(p_CDL%rho_Core(irp+ir-1,is,ia) +  &
                                                         p_CDL%rho_SemiCore(irp+ir-1,is,ia))
               enddo
               if (takeGradient) then
                  do ir = 1, n_inter
                     der_mom_l(ir,1) = der_mom_l(ir,1) + isig*fact*(p_CDL%der_rho_Core(irp+ir-1,is,ia) +  &
                                                                    p_CDL%der_rho_SemiCore(irp+ir-1,is,ia))
                  enddo
               endif
            enddo
         endif
      else
         call ErrorHandler("getMomentDensityAtPoint","Undefined charge type")
      endif
   else
      if ( nocaseCompare(momentType,"Pseudo") ) then
         do jl = 1,jmax
            mom_l(1:n_inter,jl) = p_CDL%momL_Pseudo(irp:irp+n_inter-1,jl,ia)
         enddo
      else if ( nocaseCompare(momentType,"Tilda") ) then
         do jl = 1,p_CDL%jmax
            mom_l(1:n_inter,jl) = p_CDL%momL_Tilda(irp:irp+n_inter-1,jl,ia)
         enddo
      else if ( nocaseCompare(momentType,"TotalNew") ) then
         do jl = 1,jmax
            mom_l(1:n_inter,jl) = p_CDL%momL_Total(irp:irp+n_inter-1,jl,ia)
         enddo
         if ( takeGradient ) then
            do jl = 1,jmax
               der_mom_l(1:n_inter,jl) = p_CDL%der_momL_Total(irp:irp+n_inter-1,jl,ia)
            enddo
         endif
      else if ( nocaseCompare(momentType,"TotalOld") ) then
         do jl = 1,jmax
            mom_l(1:n_inter,jl) = p_CDL%momL_TotalOld(irp:irp+n_inter-1,jl,ia)
         enddo
      else
         call ErrorHandler("getMomentDensityAtPoint","Undefined charge type")
      endif
   endif
!
   mom = ZERO
   if (takeGradient) then
      grad = ZERO
   endif
   do jl = jmax,1,-1
      l = lofj(jl)
!     ----------------------------------------------------------------
      call PolyInterp(n_inter, p_CDL%r_mesh(irp:irp+n_inter-1),       &
                      mom_l(1:n_inter,jl), r, mom_in, err)
!     ----------------------------------------------------------------
      if ( takeGradient ) then
!        -------------------------------------------------------------
         call PolyInterp(n_inter, p_CDL%r_mesh(irp:irp+n_inter-1),    &
                         der_mom_l(1:n_inter,jl), r, der_mom_in, err)
!        -------------------------------------------------------------
      endif
      kl = (l+1)*(l+1)-l+mofj(jl)
      if ( mofj(jl) == 0 ) then
         fact = ONE
      else
         fact = TWO
      endif
      mom = mom + fact*real(mom_in*ylm(kl),RealKind)
      if ( takeGradient ) then
         do i = 1, 3
            grad(i) = grad(i) + fact*real(der_mom_in*ylm(kl)*er(i)+   &
                                          mom_in*grad_ylm(kl,i),RealKind)
         enddo
      endif
   enddo
!
   end function getMomentDensityAtPoint
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMM(table_line)                             result(p_mm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer, intent(out) :: table_line(:)
   complex (kind=CmplxKind), pointer :: p_mm(:,:)
!
   if ( .not.Initialized ) then
      call ErrorHandler("getMultipoleMoment"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   endif
!
   table_line => MM_table_line
   p_mm => MultipoleMom(:,:)
!
   end function getMM
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMM_idj( ig, ia, jl )                       result(p_mm)
!  ===================================================================
   use SystemModule, only : getNumAlloyElements
   implicit none
!
   integer (kind=IntKind), intent(in) :: ig, ia, jl
   integer (kind=IntKind) :: lig
!
   complex (kind=CmplxKind) :: p_mm
!
   if ( .not.Initialized ) then
      call ErrorHandler("getMultipoleMoment"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   else if ( ig<1 .or. ig> GlobalNumAtoms ) then
      call ErrorHandler("getMultipoleMoment"," Invalid global site index",ig)
   else if ( ia<1 .or. ia> getNumAlloyElements(ig) ) then
      call ErrorHandler("getMultipoleMoment"," Invalid atom index on global site",ia,ig)
   else if ( jl<1 .or. jl> jmax_glb ) then
      call ErrorHandler("getMultipoleMoment"," Invalid j index",jl,jmax_glb)
   endif
!
   lig = MM_table_line(ig) + ia
   p_mm = MultipoleMom( lig, jl )
!
   end function getMM_idj
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMM_j( jl, table_line )                     result(p_mm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) ::  jl
   integer (kind=IntKind), pointer, intent(out) :: table_line(:)
!
   complex (kind=CmplxKind), pointer :: p_mm(:)
!
   if ( .not.Initialized ) then
      call ErrorHandler("getMultipoleMoment"," Module ChargeDensity"// &
                        " needs to be in itialized first")
   else if ( jl<1 .or. jl> jmax_glb ) then
      call ErrorHandler("getMultipoleMoment"," Invalid jl index",jl,jmax_glb)
   endif
!
   table_line => MM_table_line
   p_mm => MultipoleMom(:,jl)
!
   end function getMM_j
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isChargeComponentZero(id,jl) result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: jl
!
   logical :: flag
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('isChargeComponentZero','invalid id',id)
   else if (jl < 1 .or. jl > ChargeDensityList(id)%jmax) then
      flag = .true.
   else if (ChargeDensityList(id)%ChargeCompFlag(jl) == 0) then
      flag = .true.
   else
      flag = .false.
   endif
!
   end function isChargeComponentZero
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getChargeComponentFlag(id)                    result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer, pointer :: flag(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('isChargeComponentZero','invalid id',id)
   endif
!
   flag => ChargeDensityList(id)%ChargeCompFlag
!
   end function getChargeComponentFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine constructChargeDensity()
!  ===================================================================
   use KindParamModule, only : IntKind
!
   use GroupCommModule, only : getGroupID, GlobalSumInGroup
!
   use WriteFunctionModule, only : writeFunction
!
   use DerivativeModule, only : derv5
!
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid
!
!  use DataServiceCenterModule, only : getDataStorage, RealMark,      &
!                                      ComplexMark
   use CoreStatesModule, only : getDeepCoreDensity, getSemiCoreDensity
   use CoreStatesModule, only : getDeepCoreDensityDerivative,         &
                                getSemiCoreDensityDerivative
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use AtomModule, only : getLocalEvecNew, getLocalAtomicNumber
   use AtomModule, only : getLocalSpeciesContent
!
   use SystemSymmetryModule, only : getSymmetryFlags
!
   use PotentialTypeModule, only : isASAPotential
!
   use StepFunctionModule, only : getVolumeIntegration
   use IntegrationModule, only : calIntegration
!
   use SystemVolumeModule, only : getSystemVolume, getTotalInterstitialVolume
!
   use PolyhedraModule, only : getVolume
!
   use ValenceDensityModule, only : getValenceVPCharge, getValenceVPMoment, &
                                    getValenceSphericalElectronDensity,     &
                                    getValenceSphericalMomentDensity,       &
                                    getValenceElectronDensity, getValenceMomentDensity
!
   implicit none
!
   character(len=40) :: file_den
   character (len=6) :: denFlag, jlFlag
   character (len=2) :: specFlag
!
   integer (kind=IntKind) :: id, ia, nr, nr_ps, ir, is, ns, jl, jmax, jmt, l, jend, kmax, id_glb, lig
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind) :: evec(3), mvec(3), msgbuf(3)
   real (kind=RealKind) :: corr, qint(2), volume, rho_r,rho_i, q_tmp, q_tmp_mt, qlost, omega_vp, omega_mt
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), pointer :: rho0(:), mom0(:)
   real (kind=RealKind), pointer :: rho2p_r(:), rho2_r(:,:), rho3_r(:,:,:)
   real (kind=RealKind), pointer :: mom2_r(:,:), mom3_r(:,:,:)
   real (kind=RealKind), parameter :: PI8 = PI4*TWO
!
   complex (kind=CmplxKind) :: cfact
   complex (kind=CmplxKind), pointer :: moml(:,:), rhol(:)
   complex (kind=CmplxKind), pointer :: mom4_c(:,:,:,:)
   complex (kind=CmplxKind), pointer :: p_MultipoleMom(:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
   type (GridStruct), pointer :: Grid
!
!  ===========================================
!  Construct the Total Denity(Spherical) first
!  ===========================================
   if ( maxval(print_level) >= 0 ) then
      write(6,'(/,16x,a  )')'****************************************'
      write(6,'(  16x,a  )')'*  Output from constructChargeDensity  *'
      write(6,'(  16x,a  )')'*    in Module ChargeDensityModule     *'
      write(6,'(  16x,a,/)')'****************************************'
      call FlushFile(6)
   endif
!
   cfact = -SQRTm1
   volume = getSystemVolume()
   qint  = ZERO
   qlost = ZERO
   omega_vp = ZERO
   omega_mt = ZERO
   do id = 1,LocalNumAtoms
      Grid => getGrid(id)
      jmt = Grid%jmt
      jend = Grid%jend
      p_CDL => ChargeDensityList(id)
      nr = p_CDL%NumRPts
      nr_ps = p_CDL%NumRPts_Pseudo
      ns = 2*n_spin_cant-1
      jmax = p_CDL%jmax
      kmax = (p_CDL%lmax+1)**2
      r_mesh => ChargeDensityList(id)%r_mesh(1:nr)
      sqrt_r(0) = ZERO
      do ir = 1,nr
         sqrt_r(ir) = sqrt(r_mesh(ir))
      enddo
!
!     =============================
!     Total Spherical Densities
!     =============================
!
      p_CDL%rho_Core => getDeepCoreDensity(id)
      p_CDL%rho_SemiCore => getSemiCoreDensity(id)
      if (rad_derivative) then
         p_CDL%der_rho_Core => getDeepCoreDensityDerivative(id)
         p_CDL%der_rho_SemiCore => getSemiCoreDensityDerivative(id)
      endif
      p_CDL%ChargeCompFlag(1) = 1
      do ia = 1, p_CDL%NumSpecies
         rho2p_r => getValenceSphericalElectronDensity(id,ia)
         evec(1:3) = getLocalEvecNew(id)
         rho0 => p_CDL%rhoSph_Total(1:nr+1,ia)
         rho0(1:nr+1) = ZERO
         rho0(1:nr) = rho0(1:nr) + rho2p_r(1:nr)
!
         if ( print_level(id) >= 0 ) then
            q_tmp = getVolumeIntegration( id, jend, r_mesh, 0,           &
                                          rho2p_r(1:nr), q_tmp_mt )
            write(6,'(a,i4,a,i4,a)')'+++++ id = ',id,',   ia = ',ia,' +++++'
            write(6,'(a,f20.14)') "Charge Density :: Spherical Val Charge in MT   = ",  &
                       q_tmp_mt
            write(6,'(a,f20.14)') "Charge Density :: Spherical Val Charge in VP   = ",  &
                       q_tmp
         endif
!
         rho2_r => p_CDL%rho_Core(:,:,ia)
         rho0(1:nr) = rho0(1:nr) + rho2_r(1:nr,1)
         if ( n_spin_pola == 2 ) then
            rho0(1:nr) = rho0(1:nr) + rho2_r(1:nr,2)
         endif
!
         if ( print_level(id) >= 0 ) then
!
            den_r(0) = ZERO
            do ir = 1,nr
               den_r(ir) = rho2_r(ir,1)
            enddo
            if ( n_spin_pola == 2 ) then
               do ir = 1,nr
                  den_r(ir) = den_r(ir) + rho2_r(ir,2)
               enddo
            endif
!
            q_tmp = getVolumeIntegration( id, jend, r_mesh(1:jend), 0,    &
                                          den_r(1:jend), q_tmp_mt )
!
            den_r(0) = ZERO
            do ir = 1,nr
               den_r(ir) = rho2_r(ir,1)*r_mesh(ir)
            enddo
            if ( n_spin_pola == 2 ) then
               do ir = 1,nr
                  den_r(ir) = den_r(ir) + rho2_r(ir,2)*r_mesh(ir)
               enddo
            endif
!
            call calIntegration(nr+1,sqrt_r(0:nr),den_r(0:nr),den_r1(0:nr),3)
!
            write(6,'(a,f20.14)') "Charge Density :: DeepCore Charge in MT-Sphere = ", &
                                          den_r1(jmt)*PI8
            write(6,'(a,f20.14)') "Charge Density :: DeepCore Charge in WS-Sphere = ", &
                                          den_r1(jend)*PI8
            write(6,'(a,f20.14)') "Charge Density :: DeepCore Charge in EX-Sphere = ", &
                                          den_r1(nr)*PI8
            write(6,'(a,f20.14)') "Charge Density :: DeepCore Charge in MT        = ",  q_tmp_mt
            write(6,'(a,f20.14)') "Charge Density :: DeepCore Charge in VP        = ",  q_tmp
         endif
!
         rho2_r => p_CDL%rho_SemiCore(1:nr,1:n_spin_pola,ia)
         rho0(1:nr) = rho0(1:nr) + rho2_r(1:nr,1)
         if ( n_spin_pola == 2 ) then
            rho0(1:nr) = rho0(1:nr) + rho2_r(1:nr,2)
         endif
!
         if ( print_level(id) >= 0 ) then
!
            den_r(0)    = ZERO
            den_r(1:nr) = rho2_r(1:nr,1)
            if ( n_spin_pola == 2 ) then
               den_r(1:nr) = den_r(1:nr) + rho2_r(1:nr,2)
            endif
!
            q_tmp = getVolumeIntegration( id, jend, r_mesh(1:jend), 0,           &
                                          den_r(1:jend), q_tmp_mt )
            den_r(0)    = ZERO
            do ir = 1,nr
               den_r(ir) = rho2_r(ir,1)*r_mesh(ir)
            enddo
!
           if ( n_spin_pola == 2 ) then
               do ir = 1,nr
                  den_r(ir) = den_r(ir) + rho2_r(ir,2)*r_mesh(ir)
               enddo
            endif
!
            call calIntegration(nr+1,sqrt_r(0:nr),den_r(0:nr),den_r1(0:nr),3)
!
            write(6,'(a,f20.14)') "Charge Density :: SemiCore Charge in MT-Sphere = ", den_r1(jmt)*PI8
            write(6,'(a,f20.14)') "Charge Density :: SemiCore Charge in WS-Sphere = ", den_r1(jend)*PI8
            write(6,'(a,f20.14)') "Charge Density :: SemiCore Charge in EX-Sphere = ", den_r1(nr)*PI8
            write(6,'(a,f20.14)') "Charge Density :: SemiCore Charge in MT        = ", q_tmp_mt
            write(6,'(a,f20.14)') "Charge Density :: SemiCore Charge in VP        = ", q_tmp
!
         endif
!
!        =============================================================
!        The following code is added to take care small charge density
!        numerical noise in vaccume cell that could potentially
!        cause exchange correlation potential calculation to fail.
!        05/24/17 - Hopefully there will be better way to address
!        this problem in the future.
!        =============================================================
         do ir = 1, nr
            if (rho0(ir) < ZERO) then
               if (rho0(ir) > -TEN2m5 .or. getLocalAtomicNumber(id,ia) == 0) then
                  rho0(ir) = ZERO
               else
                  call ErrorHandler('constructChargeDensity','rho0(ir) < 0',rho0(ir),.true.)
               endif
            endif
         enddo
!        =============================================================
!
!        do ir = jend+1,nr
!           rho0(ir) = ZERO
!        enddo
         rho0(nr+1) = getValenceVPCharge(id,ia)
!
         do ir = 1,nr
            p_CDL%rhoL_Total(ir,1,ia) = cmplx(rho0(ir)/Y0,ZERO,kind=CmplxKind)
         enddo
!         p_CDL%rhoL_Total(jend+1:nr,1,ia) = CZERO
!
         if (isASAPotential()) then
            Grid => getGrid(id)
            jend = Grid%jend
            corr = rho0(1)*r_mesh(1)*r_mesh(1)*r_mesh(1)*PI2
            qint = qint + getLocalAtomicNumber(id,ia) -                  &
                   getVolumeIntegration(id,jend,r_mesh(1:jend),0,        &
                                        rho0(1:jend),truncated=.false.) + corr
         endif
!
         if ( n_spin_pola == 2 ) then
!
            mom0 => p_CDL%momSph_Total(1:nr+1,ia)
            mom0 = ZERO
            rho2_r => p_CDL%rho_Core(:,:,ia)
            mom0(1:nr) = mom0(1:nr) + rho2_r(1:nr,1) - rho2_r(1:nr,2)
            rho2_r => p_CDL%rho_SemiCore(:,:,ia)
            mom0(1:nr) = mom0(1:nr) + rho2_r(1:nr,1) - rho2_r(1:nr,2)
            mom2_r => getValenceSphericalMomentDensity(id,ia)
            mvec = getValenceVPMoment(id,ia)
            if ( n_spin_cant==1 ) then
               mom0(1:nr) = mom0(1:nr) + mom2_r(1:nr,1)
               mom0(nr+1) = mvec(3)
            else
               do is = 1,3
                  mom0(1:nr) = mom0(1:nr) + evec(is)*mom2_r(1:nr,is)
               enddo
               mom0(nr+1) = sqrt(mvec(1)**2+mvec(2)**2+mvec(3)**2)
            endif
!
!           ==========================================================
!           The following code is added to take care small moment density
!           numerical noise in vaccume cell that could potentially
!           cause exchange correlation potential calculate to fail.
!           05/24/17 - Hopefully there will be better way to address
!           this problem in the future.
!           ==========================================================
            do ir = 1,nr
               if (rho0(ir) < abs(mom0(ir))) then
                  if (rho0(ir) < TEN2m5 .or. getLocalAtomicNumber(id,ia) == 0) then
                     mom0(ir) = ZERO
                  else
!                    call ErrorHandler('constructChargeDensity','rho0 < mom0',rho0(ir),mom0(ir),.true.)
                     call WarningHandler('constructChargeDensity','rho0 < mom0',rho0(ir),mom0(ir),.true.)
                     mom0(ir) = sign(rho0(ir),mom0(ir)) ! Added on 01/12/2019 by YW
                  endif
               endif
            enddo
!           ==========================================================
!           do ir = jend+1,nr
!              mom0(ir) = ZERO 
!           enddo
            do ir = 1,nr
               p_CDL%momL_Total(ir,1,ia) = cmplx(mom0(ir)/Y0,ZERO,kind=CmplxKind)
            enddo
         endif
!
         if ( print_level(id) >= 0 ) then
!
            den_r(0) = CZERO
            do ir = 1,nr
               den_r(ir) = rho0(ir)*r_mesh(ir)
            enddo
!
            call calIntegration(nr+1,sqrt_r(0:nr),den_r(0:nr),den_r1(0:nr),3)
!
            write(6,'(a,f20.14)') "Charge Density :: Total Spherical Charge in MT = ", &
                                          den_r1(jmt)*PI8
            write(6,'(a,f20.14)') "Charge Density :: Total Spherical Charge in BS = ", &
                                          den_r1(jend)*PI8
            den_r(0)    = ZERO
            do ir = 1,nr
               den_r(ir) = rho0(ir)
            enddo
            q_tmp = getVolumeIntegration( id, jend, r_mesh(1:jend), 0, den_r(1:jend) )
            write(6,'(a,f20.14)') "Charge Density :: Total Spherical Charge in VP = ", q_tmp
         endif
!
         if ( .not.isSphericalCharge ) then
!
!           =========================================================
!           L Density components - should have been already computed in
!                                  ValenceStatesModule for l>1 and for l=0
!                                  component just copy the total spherical
!                                  charge into the real part.
!           =========================================================
            flag_jl => p_CDL%ChargeCompFlag(1:jmax)
!
            if ( .not.isChargeSymmOn ) then
               flag_jl = 0
               flag_jl(1) = 1
               do jl = 2,jmax
                  l = lofj(jl)
                  LOOP_NR: do ir = 1,nr
!                    if ( abs(p_CDL%rhoL_Valence(ir,jl,ia)) > TEN2m12 ) then
                     if ( abs(p_CDL%rhoL_Valence(ir,jl,ia)) > rho_tol ) then
                        flag_jl(jl) = 1
                        exit LOOP_NR
                     endif
                  enddo LOOP_NR
                  if ( flag_jl(jl) == 0 ) then
                     p_CDL%rhoL_Valence(1:nr,jl,ia) = CZERO
                  else
                     if ( mofj(jl) /= 0 ) then
                        do ir = 1,nr
                           rho_r = real(p_CDL%rhoL_Valence(ir,jl,ia), kind=RealKind)
                           rho_i = real(cfact*p_CDL%rhoL_Valence(ir,jl,ia), kind=RealKind)
                           if ( (abs(rho_r) /= ZERO) .and. (abs(rho_i)/abs(rho_r) < rho_tol) ) then
                              p_CDL%rhoL_Valence(ir,jl,ia) = cmplx(rho_r,ZERO,kind=CmplxKind)
                           else if ( (abs(rho_i) /= ZERO) .and. (abs(rho_r)/abs(rho_i) < rho_tol) ) then
                              p_CDL%rhoL_Valence(ir,jl,ia) = cmplx(ZERO,rho_i,kind=CmplxKind)
                              flag_jl(jl) = 2
                           else
                              flag_jl(jl) = 3
                           endif
                        enddo
                     else
                        do ir = 1,nr
                           rho_r = real(p_CDL%rhoL_Valence(ir,jl,ia),kind=RealKind)
                           p_CDL%rhoL_Valence(ir,jl,ia) = cmplx(rho_r,ZERO,kind=CmplxKind)
                        enddo
                     endif
                  endif
               enddo
!              ======================================================
!              Eliminates numerical noise of valence charge density for small r
!              ======================================================
               if ( isFitChargeOn ) then
!                 ---------------------------------------------------
                  call rtol_fit( id, ia )
!                 ---------------------------------------------------
               endif
            else
               flag_jl => getSymmetryFlags(id)
               do jl = 1,jmax
                  p_CDL%ChargeCompFlag(jl) = flag_jl(jl)
               enddo
               flag_jl => p_CDL%ChargeCompFlag(1:jmax)
               do jl = 1,jmax
                  if ( flag_jl(jl) ==0 ) then
                     p_CDL%rhoL_Valence(1:nr,jl,ia) = CZERO
                  else if ( flag_jl(jl) ==1 ) then
                     do ir = 1,nr
                        rho_r = real(p_CDL%rhoL_Valence(ir,jl,ia), kind=RealKind)
                        p_CDL%rhoL_Valence(ir,jl,ia) = cmplx(rho_r,ZERO,kind=CmplxKind)
                     enddo
                  else if ( flag_jl(jl) ==2 ) then
                     do ir = 1,nr
                        rho_i = real(cfact*p_CDL%rhoL_Valence(ir,jl,ia),kind=RealKind)
                        p_CDL%rhoL_Valence(ir,jl,ia) = cmplx(ZERO,rho_i,kind=CmplxKind)
                     enddo
                  endif
               enddo
            endif
!
            if ( print_level(id) >= 0 ) then
               q_tmp = getVolumeIntegration( id, nr, r_mesh(1:nr), kmax, jmax,  &
                                           0, p_CDL%rhoL_Valence(1:nr,1:jmax,ia),   &
                                           q_tmp_mt)
               write(6,'(a,f20.14)') "Charge Density :: Valence Charge in MT-Sphere  = ", &
                                                 q_tmp_mt
               write(6,'(a,f20.14)') "Charge Density :: Valence Charge in VP         = ", &
                                                 q_tmp
            endif
!
            do jl = 2,jmax
               p_CDL%rhoL_Total(1:nr,jl,ia) =  p_CDL%rhoL_Valence(1:nr,jl,ia)
               p_CDL%rhoL_Total(jend+1:nr,jl,ia) = CZERO
            enddo
            rho0 => p_CDL%rhoSph_Total(1:nr+1,ia)
            do ir = 1,nr
               p_CDL%rhoL_Total(ir,1,ia) = cmplx(rho0(ir)/Y0,ZERO,kind=CmplxKind)
            enddo
            p_CDL%rhoL_Total(jend+1:nr,1,ia) = CZERO
!
!           =========================================================
!           The total momemt components have to be summed up over the evec
!           directions because it is stored in ValenceStatesModule by the 
!           evec directions
!           =========================================================
!
            if ( n_spin_pola == 2 ) then
               momL => p_CDL%momL_Total(1:nr,1:jmax,ia)
               momL = CZERO
               mom4_c => getValenceMomentDensity(id)
               if ( n_spin_cant==1 ) then
                  do jl = 1,jmax
                     momL(1:nr,jl) = momL(1:nr,jl) + mom4_c(1:nr,jl,1,ia)
                  enddo
               else
                  do is = 1,3
                     do jl = 1,jmax
                        momL(1:nr,jl) = momL(1:nr,jl) + &
                                        mom4_c(1:nr,jl,is,ia)*evec(is)
                     enddo
                  enddo
               endif
               do ir = 1,nr
                  momL(ir,1) = cmplx(p_CDL%momSph_Total(ir,ia)/Y0,ZERO,kind=CmplxKind)
               enddo
!
!              =======================================================
!              The following code is added to take care small charge density
!              numerical noise in vaccume cell that could potentially
!              cause exchange correlation potential calculate to fail. 05/29/17 - Hopefully 
!              there will be better way to address this problem in the future.
!              =======================================================
!              if (getLocalAtomicNumber(id,ia) == 0) then
!                 do jl = 1,jmax
!                    do ir = 1, nr
!                       if (abs(p_CDL%rhoL_Total(ir,jl,ia)) < TEN2m5) then
!                          momL(ir,jl) = CZERO
!                       endif
!                    enddo
!                 enddo
!              endif
!              =======================================================
            endif
!
            q_tmp = getVolumeIntegration( id, nr, r_mesh, kmax, jmax,   &
                                            0, p_CDL%rhoL_Total(1:nr,1:jmax,ia), q_tmp_mt )
            qlost = qlost+(getLocalAtomicNumber(id,ia) - q_tmp)
            if (print_level(id) >= 0) then
               write(6,'(a,f20.14)')'Checking -- Missing charge in atomic cell      = ',qlost
            endif
!qlost = qlost+(getLocalAtomicNumber(ia) - q_tmp_mt)
!if (print_level(id) >= 0) then
!   write(6,'(a,f18.14)')'Checking --- Interstitial charge  = ',qlost
!endif
!
!           ==========================================================
!           Use getVolumeIntegration to calculate the volume will help
!           offsetting the numerical error in the determination of
!           interstitial charge density.
!           ==========================================================
            den_r = ONE
            omega_vp = omega_vp + getVolumeIntegration( id, nr, r_mesh, 0, den_r(1:nr), q_tmp_mt )
            omega_mt = omega_mt + q_tmp_mt
!           ==========================================================
         else
            omega_vp = omega_vp + getVolume(id)
            omega_mt = omega_mt + PI4*r_mesh(jmt)**3*THIRD
         endif
!
         if (rad_derivative) then
            rho2p_r => getValenceSphericalElectronDensity(id,ia,isDerivative=.true.)
            rho0 => p_CDL%der_rhoSph_Total(:,ia)
            rho0 = ZERO
            rho0(1:nr) = rho2p_r(1:nr)
!
            rho2_r => p_CDL%der_rho_Core(:,:,ia)
            rho0(1:nr) = rho0(1:nr) + rho2_r(1:nr,1)
            if ( n_spin_pola == 2 ) then
               rho0(1:nr) = rho0(1:nr) + rho2_r(1:nr,2)
            endif
!
            rho2_r => p_CDL%der_rho_SemiCore(:,:,ia)
            rho0(1:nr) = rho0(1:nr) + rho2_r(1:nr,1)
            if ( n_spin_pola == 2 ) then
               rho0(1:nr) = rho0(1:nr) + rho2_r(1:nr,2)
            endif
!
            do ir = 1,nr
               p_CDL%der_rhoL_Total(ir,1,ia) = cmplx(rho0(ir)/Y0,ZERO,kind=CmplxKind)
            enddo
!
            if ( n_spin_pola == 2 ) then
               mom0 => p_CDL%der_momSph_Total(:,ia)
               mom0 = ZERO
               rho2_r => p_CDL%der_rho_Core(:,:,ia)
               mom0(1:nr) = mom0(1:nr) + rho2_r(1:nr,1) - rho2_r(1:nr,2)
               rho2_r => p_CDL%der_rho_SemiCore(:,:,ia)
               mom0(1:nr) = mom0(1:nr) + rho2_r(1:nr,1) - rho2_r(1:nr,2)
!
               mom2_r => getValenceSphericalMomentDensity(id,ia,isDerivative=.true.)
               mvec = getValenceVPMoment(id,ia)
               if ( n_spin_cant==1 ) then
                  mom0(1:nr) = mom0(1:nr) + mom2_r(1:nr,1)
               else
                  do is = 1,3
                     mom0(1:nr) = mom0(1:nr) + evec(is)*mom2_r(1:nr,is)
                  enddo
               endif
!
               do ir = 1,nr
                  p_CDL%der_momL_Total(ir,1,ia) = cmplx(mom0(ir)/Y0,ZERO,kind=CmplxKind)
               enddo
            endif
!
            if ( .not.isSphericalCharge ) then
               do jl = 2,jmax
                  rhol => getValenceElectronDensity(id,ia,jl,isDerivative=.true.)
                  p_CDL%der_rhoL_Total(1:nr,jl,ia) =  rhol(1:nr)
               enddo
               rho0 => p_CDL%der_rhoSph_Total(:,ia)
               do ir = 1,nr
                  p_CDL%der_rhoL_Total(ir,1,ia) = cmplx(rho0(ir)/Y0,ZERO,kind=CmplxKind)
               enddo
               if ( n_spin_pola == 2 ) then
                  momL => p_CDL%der_momL_Total(:,:,ia)
                  momL = CZERO
                  mom4_c => getValenceMomentDensity(id,isDerivative=.true.)
                  if ( n_spin_cant==1 ) then
                     do jl = 1,jmax
                        momL(1:nr,jl) = momL(1:nr,jl) + mom4_c(1:nr,jl,1,ia)
                     enddo
                  else
                     do is = 1, 3
                        do jl = 1, jmax
                           momL(1:nr,jl) = momL(1:nr,jl) + mom4_c(1:nr,jl,is,ia)*evec(is)
                        enddo
                     enddo
                  endif
                  do ir = 1,nr
                     momL(ir,1) = cmplx(p_CDL%der_momSph_Total(ir,ia)/Y0,ZERO,kind=CmplxKind)
                  enddo
               endif
            endif
!
!           ==========================================================
!           Checking the radial derivatives...
!           ==========================================================
            if (maxval(print_level) >= 1) then
               write(denFlag,'(i6)')100000+getGlobalIndex(id)
               denFlag(1:1) = 'n'
               write(specFlag,'(i2)')10+ia
               specFlag(1:1) = 'c'
               file_den = 'SphDensity'//'_'//denFlag//specFlag
!              -------------------------------------------------------
               call derv5(p_CDL%rhoSph_Total(:,ia),den_r(1:),r_mesh,nr)
!              -------------------------------------------------------
               call writeFunction(file_den,nr,r_mesh,                 &
                                  p_CDL%rhoSph_Total(:,ia),           &
                                  p_CDL%der_rhoSph_Total(:,ia),       &
                                  den_r(1:))
!              -------------------------------------------------------
               if (n_spin_pola == 2) then
                  file_den = 'SphMomDensity'//'_'//denFlag//specFlag
!                 ----------------------------------------------------
                  call derv5(p_CDL%momSph_Total(:,ia),den_r(1:),r_mesh,nr)
!                 ----------------------------------------------------
                  call writeFunction(file_den,nr,r_mesh,                 &
                                     p_CDL%momSph_Total(:,ia),           &
                                     p_CDL%der_momSph_Total(:,ia),       &
                                     den_r(1:))
!                 ----------------------------------------------------
               endif
!
               if ( .not.isSphericalCharge ) then
                  do jl = 1,jmax
                     if (p_CDL%ChargeCompFlag(jl) > 0) then
                        write(jlFlag,'(i4)')1000+jl
                        jlFlag(1:2) = "jl"
                        file_den = 'FullDensity'//'_'//denFlag//specFlag//jlFlag
                        den_r1(1:) = real(p_CDL%rhoL_Total(:,jl,ia))
!                       ----------------------------------------------
                        call derv5(den_r1(1:),den_r(1:),r_mesh,nr)
!                       ----------------------------------------------
                        call writeFunction(file_den,nr,r_mesh,        &
                                           den_r1(1:),                &
                                           real(p_CDL%der_rhoL_Total(:,jl,ia)),&
                                           den_r(1:))
!                       ----------------------------------------------
                        if (n_spin_pola == 2) then
                           file_den = 'FullMomDensity'//'_'//denFlag//specFlag//jlFlag
                           den_r1(1:) = real(p_CDL%momL_Total(:,jl,ia))
!                          -------------------------------------------
                           call derv5(den_r1(1:),den_r(1:),r_mesh,nr)
!                          -------------------------------------------
                           call writeFunction(file_den,nr,r_mesh,        &
                                              den_r1(1:),                &
                                              real(p_CDL%der_momL_Total(:,jl,ia)),&
                                              den_r(1:))
!                          -------------------------------------------
                        endif
                     endif
                  enddo
               endif
            endif
!           ==========================================================
         endif
!
      enddo
   enddo
   msgbuf(1) = qlost; msgbuf(2) = omega_vp; msgbuf(3) = omega_mt
   call GlobalSumInGroup(GroupID,msgbuf,3)
   qlost = msgbuf(1); omega_vp = msgbuf(2); omega_mt = msgbuf(3)
!  ===================================================================
!  Determine the missing charge density, which is to be added to the
!  interstitial region. -Yang Wang @ 12/17/2018.
!
!  Using (omega_vp-omega_mt) instead of getTotalInterstitialVolume()
!  helps reducing the numerical error caused by truncation scheme in 
!  the volume integration.
!  ===================================================================
   qlost = qlost/(omega_vp-omega_mt)
   if (maxval(print_level) >= 0) then
      write(6,'(a)')'+++++++++++++++++++++++++++++++++'
      write(6,'(a,f20.14)')'Checking -- Missing density in interstial area = ',qlost
   endif
!! qlost = ZERO
!
   if (isASAPotential()) then
      rhoint = qint(1)/getSystemVolume()
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,rhoint)
!     ----------------------------------------------------------------
      if (rhoint > ZERO) then
         do id=1,LocalNumAtoms
            p_CDL=>ChargeDensityList(id)
            Grid => getGrid(id)
            jend = Grid%jend
            do ia = 1, p_CDL%NumSpecies
               rho0 => p_CDL%rhoSph_Total(1:jend,ia)
               do ir=1,jend
                  rho0(ir)=rho0(ir)+rhoint ! if rhoint<0, may cause error
               enddo
            enddo
         enddo
      endif
   endif
!
   if ( .not.isSphericalCharge ) then
      do id = 1,LocalNumAtoms
         Grid => getGrid(id)
         jmt = Grid%jmt
         do ia = 1, p_CDL%NumSpecies
            p_CDL => ChargeDensityList(id)
            nr = p_CDL%NumRPts
            rho0 => p_CDL%rhoSph_Total(1:nr+1,ia)
!===
!This is for testing on 12/13/2018
!Needs to be deleted.
!p_CDL%rhoL_Total(1:jmt,1,ia) = p_CDL%rhoL_Total(1:jmt,1,ia) - qlost/Y0
!p_CDL%rhoL_Total(jmt+1:,1,ia) = qlost/Y0
!rho0(1:jmt)=rho0(1:jmt)-qlost
!rho0(jmt+1:)=qlost
!===
            p_CDL%rhoL_Total(jmt+1:nr,1,ia) = p_CDL%rhoL_Total(jmt+1:nr,1,ia) + qlost/Y0
            rho0(jmt+1:nr) = rho0(jmt+1:nr) + qlost
            q_tmp = getVolumeIntegration( id, nr, r_mesh, kmax, jmax,   &
                                            0, p_CDL%rhoL_Total(1:nr,1:jmax,ia), &
                                            q_tmp_mt )
            if (print_level(id) >= 0) then
               write(6,'(a,2i5,f20.14)') "id, ia, Total Charge in MT = ", &
                                                 id, ia, q_tmp_mt
               write(6,'(a,2i5,f20.14)') "id, ia, Total Charge in VP = ", &
                                                 id, ia, q_tmp
            endif
            rho0(nr+1) = q_tmp_mt
         enddo
      enddo
!     ========================================
!     Construct the Pseudo Denity
!     ========================================
      if (.not.isMTFP) then
         call calPseudoDensity("Total")
      endif
!===
!This is for testing on 12/13/2018
!Needs to be deleted.
!p_CDL%rhoSph_Pseudo = ZERO
!p_CDL%rhoL_Pseudo = CZERO
!===
!
!     ========================================
!     Construct the Tilda Denity
!     ========================================
!
!Yang qint = qint/volume
!     ---------------------------------------------------------------
!Yang call GlobalSumInGroup(GroupID,qint,2)
!     ---------------------------------------------------------------
!Yang rhoint = qint(1)
      do id = 1,LocalNumAtoms
         p_CDL => ChargeDensityList(id)
         nr = p_CDL%NumRPts
         nr_ps = p_CDL%NumRPts_Pseudo
         jmax = p_CDL%jmax
         flag_jl => p_CDL%ChargeCompFlag(1:jmax)
!
         do ia = 1, p_CDL%NumSpecies
            p_CDL%rhoSph_Tilda(1:nr,ia) = ZERO
            p_CDL%rhoL_Tilda(1:nr,1:jmax,ia) = ZERO
!           =========================================================
!           Notes on 11/22/2018 -Yang
!           The original code in the following lines calculates 
!           rhoL_Tilda and rhoSph_Tilda up to nr.
!           This appears to be a bug. I replaced nr with nr_ps.
!           =========================================================
            if (isMTFP) then
               p_CDL%rhoSph_Tilda(1:nr_ps,ia) = p_CDL%rhoSph_Total(1:nr_ps,ia)
               p_CDL%rhoL_Tilda(1:nr_ps,1,ia) = p_CDL%rhoL_Total(1:nr_ps,1,ia)
               do jl = 2,jmax
                  if ( flag_jl(jl) /= 0 ) then
                     p_CDL%rhoL_Tilda(1:nr_ps,jl,ia) = p_CDL%rhoL_Total(1:nr_ps,jl,ia)
                  endif
               enddo
            else
               p_CDL%rhoSph_Tilda(1:nr_ps,ia) = p_CDL%rhoSph_Total(1:nr_ps,ia)  -       &
                                                p_CDL%rhoSph_Pseudo(1:nr_ps,ia)
!
               p_CDL%rhoL_Tilda(1:nr_ps,1,ia) = p_CDL%rhoL_Total(1:nr_ps,1,ia) -        &
                                                p_CDL%rhoL_Pseudo(1:nr_ps,1,ia)
!===
!This is for testing on 12/13/2018
!Needs to be deleted.
!p_CDL%rhoL_Tilda(1:nr_ps,1,ia) = p_CDL%rhoL_Tilda(1:nr_ps,1,ia) - qlost/Y0
!p_CDL%rhoSph_Tilda(1:nr_ps,ia)=p_CDL%rhoSph_Tilda(1:nr_ps,ia)-qlost
!===
               do jl = 2,jmax
                  if ( flag_jl(jl) /= 0 ) then
                     p_CDL%rhoL_Tilda(1:nr_ps,jl,ia) = p_CDL%rhoL_Total(1:nr_ps,jl,ia)- &
                                                       p_CDL%rhoL_Pseudo(1:nr_ps,jl,ia)
                  endif
               enddo
            endif
!           =========================================================
         enddo
      enddo
!
!     ===============================================================
!     Calculate multipole moments and rho0_neutral, which is the density
!     to neutralize rhoL_Tilda density
!     ---------------------------------------------------------------
      call constructMultipoleMoments()
!     ---------------------------------------------------------------
!
      if (.not.isASAPotential()) then
         rhoint = rho0_neutral
      endif
!
      if ( print_level(1) >= 0 ) then
         write(6,'(a,f20.14)') "Charge Density :: Neutralizing Charge density  = ",  &
                    rho0_neutral
         write(6,'(a,f20.14)') "Charge Density :: Neutralizing Charge Per Atom = ",  &
                    rho0_neutral*volume/GlobalNumAtoms
         write(6,'(80(''-''),/)')
      endif
!
      if ( print_level(1) >= 0 ) then
         call printMultipoleMom()
      endif
!
!     ================================================================
!     Correct the pseudo charge density for charge neutrality
!     ================================================================
!
      if (isMTFP) then
         rho0_neutral = rho0_neutral*volume/getTotalInterstitialVolume()
         p_MultipoleMom => MultipoleMom(:,1)
         p_MultipoleMom = CZERO
         do id = 1, LocalNumAtoms
            p_CDL => ChargeDensityList(id)
            id_glb = getGlobalIndex(id)
            do ia = 1, p_CDL%NumSpecies
               lig = MM_table_line(id_glb) + ia
               do ir = 1, p_CDL%NumRPts_Pseudo
                  p_CDL%rhoL_Tilda(ir,1,ia) = p_CDL%rhoL_Tilda(ir,1,ia) + rho0_neutral/Y0
               enddo
               p_CDL%multipole_mom(1,ia) = p_CDL%multipole_mom(1,ia) +   &
                                           rho0_neutral/Y0*PI4*THIRD*p_CDL%RcutPseudo**3
               p_MultipoleMom(lig) = p_CDL%multipole_mom(1,ia)
            enddo
         enddo
!        -------------------------------------------------------------
         call GlobalSumInGroup(GroupID,p_MultipoleMom,MM_size)
!        -------------------------------------------------------------
      else
         do id = 1,LocalNumAtoms
            p_CDL => ChargeDensityList(id)
            jmax = p_CDL%jmax
            nr = p_CDL%NumRPts
            do ia = 1, p_CDL%NumSpecies
               p_CDL%rhoSph_Pseudo(1:nr,ia) = p_CDL%rhoSph_Pseudo(1:nr,ia) - rho0_neutral
               p_CDL%rhoL_Pseudo(1:nr,1,ia) = p_CDL%rhoL_Pseudo(1:nr,1,ia) - rho0_neutral/Y0
            enddo
         enddo
      endif
!===
!This is for testing on 12/13/2018
!Needs to be deleted.
!p_CDL%rhoSph_Pseudo = ZERO
!p_CDL%rhoL_Pseudo = CZERO
!===
!
   endif
!
   call FlushFile(6)
!
   end subroutine constructChargeDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calPseudoDensity(rho_type)
!  ===================================================================
!
   use AtomModule, only : getRcutPseudo
   use RadialGridModule, only : getGrid
!
   implicit none
!
   character (len=*), intent(in) :: rho_type
!
   integer (kind=IntKind) :: id, ia, ir, jl, m, l
   integer (kind=IntKind) :: jmax, NumRs
   integer (kind=IntKind), pointer :: flag_jl(:)
   integer (kind=IntKind), allocatable :: ir0(:)
!
   real (kind=RealKind) :: r0, r, r2, r3, rl, tol, efact
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), allocatable :: rhol_r(:),  rhol_c(:)
   real (kind=RealKind), allocatable :: drhol_r(:), ddrhol_r(:)
   real (kind=RealKind), allocatable :: drhol_c(:), ddrhol_c(:)
!
   complex (kind=CmplxKind) :: rhol_r0, drhol_r0, ddrhol_r0, cfact
   complex (kind=CmplxKind), pointer :: rhol(:,:,:), rhol_tmp(:,:)
   complex (kind=CmplxKind), allocatable :: ac(:), bc(:), cc(:), dc(:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   tol = TEN2m6
   cfact = -SQRTm1
!
   allocate( ir0(jmax_max) )
   allocate( rhol_r(NumRPts_max), rhol_c(NumRPts_max) )
   allocate( ac(jmax_max), bc(jmax_max), cc(jmax_max), dc(jmax_max) )
   allocate( drhol_r(NumRPts_max), ddrhol_r(NumRPts_max) )
   allocate( drhol_c(NumRPts_max), ddrhol_c(NumRPts_max) )
!
   do id = 1, LocalNumAtoms
      p_CDL => ChargeDensityList(id)
      NumRs = p_CDL%NumRPts
      jmax = p_CDL%jmax
      if (jmax>jmax_max .or. jmax<0 .or. NumRs>NumRPts_max .or. NumRs<0 ) then
         write(6,*) jmax_max,NumRPts_max
         call ErrorHandler("calPseudoDensity","Wrong jmax, NumRs", jmax,NumRs)
      endif
      p_CDL%rhoL_Pseudo = CZERO
      r_mesh => p_CDL%r_mesh(1:NumRs)
!
      if (rho_type(1:5)=="ValPs") then
!        do ia = 1, p_CDL%NumSpecies
!           do jl = 1,jmax
!              p_CDL%rhoL_ValPseudo(1:NumRs,jl,ia) = p_CDL%rhoL_Valence(1:NumRs,jl,ia)
!           enddo
!        enddo
         p_CDL%rhoL_ValPseudo = p_CDL%rhoL_Valence
!        rhol => p_CDL%rhoL_Valence(1:NumRs,1:jmax,1:p_CDL%NumSpecies)
         rhol => p_CDL%rhoL_Valence
      else
!        do ia = 1, p_CDL%NumSpecies
!           do jl = 1,jmax
!              p_CDL%rhoL_Pseudo(1:NumRs,jl,ia) = p_CDL%rhoL_Total(1:NumRs,jl,ia)
!           enddo
!        enddo
         p_CDL%rhoL_Pseudo = p_CDL%rhoL_Total
!        rhol => p_CDL%rhoL_Total(1:NumRs,1:jmax,1:p_CDL%NumSpecies)
         rhol => p_CDL%rhoL_Total
      endif
!
      do ir = 1,NumRs
         sqrt_r(ir) = sqrt(r_mesh(ir))
      enddo
!
      flag_jl => p_CDL%ChargeCompFlag(1:jmax)
!
!     get the radius for fitting
!
      ir = p_CDL%NumRPts_Pseudo
      ir0(1:jmax) = ir
      r0 = r_mesh(ir)
      p_CDL%RcutPseudo = r0
!
      do ia = 1, p_CDL%NumSpecies
         rhol_r = ZERO
         rhol_c = ZERO
         do jl = 1,jmax
            if ( flag_jl(jl) /= 0 ) then
               l =lofj(jl)
               do ir = 1,NumRs
                  rhol_r(ir) = real(rhol(ir,jl,ia),kind=RealKind)
                  rhol_c(ir) = real(cfact*rhol(ir,jl,ia),kind=RealKind)
               enddo
!              -------------------------------------------------------
               call newder(rhol_r,drhol_r,sqrt_r(1:),NumRs)
!              -------------------------------------------------------
               do ir = 1, NumRs
                  drhol_r(ir) = HALF*drhol_r(ir)/sqrt_r(ir)
               enddo
!              -------------------------------------------------------
               call newder(drhol_r,ddrhol_r,sqrt_r(1:),NumRs)
!              -------------------------------------------------------
               do ir = 1, NumRs
                  ddrhol_r(ir) = HALF*ddrhol_r(ir)/sqrt_r(ir)
               enddo
!              -------------------------------------------------------
               call newder(rhol_c,drhol_c,sqrt_r(1:),NumRs)
!              -------------------------------------------------------
               do ir = 1, NumRs
                  drhol_c(ir) = HALF*drhol_c(ir)/sqrt_r(ir)
               enddo
!              -------------------------------------------------------
               call newder(drhol_c,ddrhol_c,sqrt_r(1:),NumRs)
!              -------------------------------------------------------
               do ir = 1, NumRs
                  ddrhol_c(ir) = HALF*ddrhol_c(ir)/sqrt_r(ir)
               enddo
               ir = ir0(jl)
               rhol_r0   = cmplx(rhol_r(ir),rhol_c(ir),kind=CmplxKind)
               drhol_r0  = cmplx(drhol_r(ir),drhol_c(ir),kind=CmplxKind)
               ddrhol_r0 = cmplx(ddrhol_r(ir),ddrhol_c(ir),kind=CmplxKind)
!
               if ( Fit_Method == 1 ) then
!                 ----------------------------------------------------
                  call findQuadraticFit( l, r0, rhol_r0, drhol_r0,        &
                                         ac(jl), bc(jl), cc(jl) )
!                 ----------------------------------------------------
               else if ( Fit_Method == 2 ) then
!                 ----------------------------------------------------
                  call findCubicFit( l, r0, rhol_r0, drhol_r0, ddrhol_r0, &
                                     ac(jl), bc(jl), cc(jl), dc(jl) )
!                 ----------------------------------------------------
               else
                 call ErrorHandler('calPseudoDensity', &
                                   'Unknown Fitting Method', Fit_Method)
               endif
            endif
         enddo
         if (rho_type(1:5)=="ValPs") then
            rhol_tmp => p_CDL%rhoL_ValPseudo(1:NumRs,1:jmax,ia)
         else
            rhol_tmp => p_CDL%rhoL_Pseudo(1:NumRs,1:jmax,ia)
         endif
         if ( Fit_Method == 1 ) then
!           using the quadratic fitting method
            do ir = 1,ir0(1)-1
               r  = r_mesh(ir)
               r2 = r*r
               rl = ONE
               jl = 0
               do l = 0, p_CDL%lmax
                  do m = 0, l
                     jl = jl + 1
                     efact = ONE
!                     if (jl==1) then
!                        efact = 3.0d0/(3.0d0+exp(-6.0d0*r**7/(r0**7)))
!                     endif
                     if ( flag_jl(jl) /= 0 ) then
                        if ( ir <= ir0(1) ) then
                           rhol_tmp(ir,jl) = &
                                          efact*rl*(cc(jl)*r2+bc(jl)*r+ac(jl))
                        endif
                     endif
                  enddo
                  rl = rl*r
               enddo
            enddo
         else if ( Fit_Method == 2 ) then
!           using the cubic fitting method
            do ir = 1,ir0(1)-1
               r  = r_mesh(ir)
               r2 = r*r
               r3 = r*r*r
               rl = ONE
               jl = 0
               do l = 0, p_CDL%lmax
                  do m = 0, l
                     jl = jl + 1
                     efact = ONE
!                     if (jl==1) then
!                        efact = 1/(1+.5*exp((-4*r*r/(r0*r0))))
!                     endif
                     if ( flag_jl(jl) /= 0 ) then
                        if ( ir <= ir0(1) ) then
                           rhol_tmp(ir,jl) = &
                                   efact*rl*(dc(jl)*r3+cc(jl)*r2+bc(jl)*r+ac(jl))
                        endif
                     endif
                  enddo
                  rl = rl*r
               enddo
            enddo
         else
            call ErrorHandler('calPseudoDensity', &
                              'Unknown Fitting Method', Fit_Method)
         endif
         if (rho_type(1:5)/="ValPs") then
            do ir = 1,NumRs
               p_CDL%rhoSph_Pseudo(ir,ia) = real(p_CDL%rhoL_Pseudo(ir,1,ia),kind=RealKind)*Y0
            enddo
         endif
      enddo
   enddo
!
   deallocate( ir0 )
   deallocate( rhol_r, rhol_c )
   deallocate( ac, bc, cc, dc )
   deallocate( drhol_r, ddrhol_r, drhol_c, ddrhol_c )
!
   end subroutine calPseudoDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine constructMultipoleMoments()
!  ===================================================================
   use DataServiceCenterModule, only : getDataStorage, ComplexMark
   use IntegrationModule, only : calIntegration
   use AtomModule, only : getLocalAtomicNumber, getLocalSpeciesContent
   use ChemElementModule, only : getZtot
   use SystemModule, only : getAdditionalElectrons
   use SystemVolumeModule, only : getSystemVolume
   use Atom2ProcModule, only : getGlobalIndex
   use InterpolationModule, only : FitInterp
!
   implicit none
!
   integer (kind=IntKind) :: nr, nr_ps, id, ia, ir, id_glb, pe, jmt
   integer (kind=IntKind) :: jl, l, jmax, lig
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind) :: Z
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: dps
   complex (kind=CmplxKind), allocatable :: MM_tmp(:,:,:)
   complex (kind=CmplxKind), pointer :: p_mm(:)
   complex (kind=CmplxKind), pointer :: densityTilda_l(:,:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
   type (GridStruct), pointer :: Grid
!
   rho0_neutral = ZERO
!  ======================
!  Loop over local atoms
!  ======================
   do id = 1, LocalNumAtoms
      Grid => getGrid(id)
      jmt = Grid%jmt
      p_CDL => ChargeDensityList(id)
      jmax   =  p_CDL%jmax
      nr     =  p_CDL%NumRPts
      nr_ps  =  p_CDL%NumRPts_Pseudo
      r_mesh => p_CDL%r_mesh(1:nr)
!
      sqrt_r(0) = ZERO
      do ir = 1, nr
         sqrt_r(ir) = sqrt(r_mesh(ir))
      enddo
!
      do ia = 1, p_CDL%NumSpecies
         densityTilda_l => p_CDL%rhoL_Tilda(1:nr,1:jmax,ia)
!
         ir      = 0
         flag_jl => p_CDL%ChargeCompFlag(1:jmax)
!
         p_mm => p_CDL%multipole_mom(1:jmax,ia)
!
         p_mm = CZERO
         do jl = 1, jmax
            if ( flag_jl(jl) /= 0 ) then
               do ir = 1, nr_ps
                  den_in(ir) = densityTilda_l(ir,jl)
               enddo
!              -------------------------------------------------------
               call FitInterp( 4, sqrt_r(1:4), den_in(1:4), ZERO,         &
                               den_in(0), dps )
!              -------------------------------------------------------
               l = lofj(jl)
!              -------------------------------------------------------
               call calIntegration(nr_ps+1, sqrt_r(0:nr_ps), den_in(0:nr_ps), &
                                   den_out(0:nr_ps), 5+2*l )
!              -------------------------------------------------------
               p_mm(jl) = TWO*den_out(nr_ps)*PI4 
            endif
         enddo
         Z = getLocalAtomicNumber(id,ia)
!
         rho0_neutral = rho0_neutral + (Z -real(p_mm(1),kind=RealKind)*Y0)*getLocalSpeciesContent(id,ia)
      enddo
      if (print_level(id) >= 1) then
         write(6,'(/,a)')'In constructMultipoleMoments......'
         write(6,'(a,3i5,2x,f20.14)')'id,jmt,nr_ps, net charge: ',id,jmt,nr_ps,rho0_neutral
      endif
   enddo
!
   call GlobalSumInGroup(GroupID,rho0_neutral)
!
   rho0_neutral = rho0_neutral + getAdditionalElectrons()
   rho0_neutral = rho0_neutral/getSystemVolume()
   if (maxval(print_level) >= 1) then
      write(6,'(a,f20.14,/)')'Neutralizing density: ',rho0_neutral
   endif
!
   MultipoleMom = CZERO
   do id = 1,LocalNumAtoms
      id_glb = getGlobalIndex(id)
      p_CDL => ChargeDensityList(id)
      jmax   =  p_CDL%jmax
      do ia = 1, p_CDL%NumSpecies
         lig = MM_table_line(id_glb) + ia
         p_mm => p_CDL%multipole_mom(1:jmax,ia)
         do jl = 1,jmax
            MultipoleMom(lig,jl) = p_mm(jl)
         enddo
      enddo
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,MultipoleMom,MM_size,jmax_glb)
!  -------------------------------------------------------------------
!
   end subroutine constructMultipoleMoments
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVPCharge(dtype,id,ia)                       result(chg)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: dtype
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   integer (kind=IntKind) :: nr
!
   real (kind=IntKind) :: chg
!
   if ( id > LocalNumAtoms ) then
      call ErrorHandler("getVPCharge",'Invalid atom index',id)
   endif
!
   nr = ChargeDensityList(id)%NumRPts
   if ( dtype=='Old' ) then
      chg = ChargeDensityList(id)%rhoSph_TotalOld(nr+1,ia)
   else if ( dtype=="New" ) then
      chg = ChargeDensityList(id)%rhoSph_Total(nr+1,ia)
   endif
!
   end function getVPCharge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVPMomSize(dtype,id,ia)                       result(mom)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: dtype
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   integer (kind=IntKind) :: nr
!
   real (kind=IntKind) :: mom
!
   if ( id > LocalNumAtoms ) then
      call ErrorHandler("getVPMomSize",'Invalid atom index',id)
   else if ( n_spin_pola/=2 ) then
      call ErrorHandler("getVPMomSize",'Invalid spin parameter',n_spin_pola)
   endif
!
   nr = ChargeDensityList(id)%NumRPts
   if ( dtype=='Old' ) then
      mom = ChargeDensityList(id)%momSph_TotalOld(nr+1,ia)
   else if ( dtype=="New" ) then
      mom = ChargeDensityList(id)%momSph_Total(nr+1,ia)
   endif
!
   end function getVPMomSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setVPCharge(chg,id,ia)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), intent(in) :: chg
!
   integer (kind=IntKind) :: nr
!
   if ( id > LocalNumAtoms ) then
      call ErrorHandler("setVPCharge",'Invalid atom index',id)
   endif
!
   nr = ChargeDensityList(id)%NumRPts
   ChargeDensityList(id)%rhoSph_Total(nr+1,ia) = chg
!
   end subroutine setVPCharge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setVPMomSize(mom,id,ia)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), intent(in) :: mom
!
   integer (kind=IntKind) :: nr
!
   if ( id > LocalNumAtoms ) then
      call ErrorHandler("setVPMomSize",'Invalid atom index',id)
   else if ( n_spin_pola/=2 ) then
      call ErrorHandler("setVPMomSize",'Invalid spin parameter',n_spin_pola)
   endif
!
   nr = ChargeDensityList(id)%NumRPts
   ChargeDensityList(id)%momSph_Total(nr+1,ia) = mom
!
   end subroutine setVPMomSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateTotalDensity(setValenceVPCharge, setValenceVPMomentSize)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id, ia, ir, nr, jl
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind), pointer :: old_SphRho(:)
   real (kind=RealKind), pointer :: new_SphRho(:)
   real (kind=RealKind), pointer :: old_SphMom(:)
   real (kind=RealKind), pointer :: new_SphMom(:)
!
   complex (kind=CmplxKind), pointer :: old_rhoL(:)
   complex (kind=CmplxKind), pointer :: new_rhoL(:)
   complex (kind=CmplxKind), pointer :: old_momL(:)
   complex (kind=CmplxKind), pointer :: new_momL(:)
!
   type (ChargeDensityStruct), pointer :: p_CDL
!
   interface
      subroutine setValenceVPCharge(id,ia,vc)
         use KindParamModule, only : IntKind, RealKind  
         integer (kind=IntKind), intent(in) :: id, ia
         real (kind=RealKind), intent(in) :: vc
      end subroutine setValenceVPCharge
!
      subroutine setValenceVPMomentSize(id,ia,vm)
         use KindParamModule, only : IntKind, RealKind
         integer (kind=IntKind), intent(in) :: id, ia
         real (kind=RealKind), intent(in) :: vm
      end subroutine setValenceVPMomentSize
   end interface
!
   if ( isSphericalCharge ) then
      do id = 1, LocalNumAtoms
         p_CDL=>ChargeDensityList(id)
         nr = p_CDL%NumRPts+1
!
         do ia = 1, p_CDL%NumSpecies
            old_SphRho => p_CDL%rhoSph_TotalOld(1:nr,ia)
            new_SphRho => p_CDL%rhoSph_Total(1:nr,ia)
            old_SphRho(1:nr) = new_SphRho(1:nr)
!
            call setValenceVPCharge(id,ia,new_SphRho(nr))
            if (n_spin_pola == 2) then
               old_SphMom => p_CDL%momSph_TotalOld(1:nr,ia)
               new_SphMom => p_CDL%momSph_Total(1:nr,ia)
               old_SphMom(1:nr) = new_SphMom(1:nr)
               call setValenceVPMomentSize(id,ia,new_SphMom(nr))
            endif
         enddo
!
      enddo
   else
      do id = 1, LocalNumAtoms
         p_CDL=>ChargeDensityList(id)
         nr = p_CDL%NumRPts
!
!        Set the l components
!
         flag_jl => p_CDL%ChargeCompFlag(1:p_CDL%jmax)
!
         do ia = 1, p_CDL%NumSpecies
            do jl = 1,p_CDL%jmax
               if ( flag_jl(jl)/=0 ) then
                  old_rhoL => p_CDL%rhoL_TotalOld(1:nr,jl,ia)
                  new_rhoL => p_CDL%rhoL_Total(1:nr,jl,ia)
                  old_rhoL(1:nr) = new_rhoL(1:nr)
                  if (n_spin_pola == 2) then
                     old_momL => p_CDL%momL_TotalOld(1:nr,jl,ia)
                     new_momL => p_CDL%momL_Total(1:nr,jl,ia)
                     old_momL(1:nr) = new_momL(1:nr)
                  endif
               endif
            enddo
         enddo
!
!        Set the spherical component
!
         do ia = 1, p_CDL%NumSpecies
            old_SphRho => p_CDL%rhoSph_TotalOld(1:nr+1,ia)
            new_SphRho => p_CDL%rhoSph_Total(1:nr+1,ia)
            do ir = 1,nr
               new_SphRho(ir) = real(p_CDL%rhoL_Total(ir,1,ia),kind=RealKind)*Y0
               old_SphRho(ir) = new_SphRho(ir)
            enddo
            call setValenceVPCharge(id,ia,new_SphRho(nr+1))
!
            if (n_spin_pola == 2) then
               old_SphMom => p_CDL%momSph_TotalOld(1:nr+1,ia)
               new_SphMom => p_CDL%momSph_Total(1:nr+1,ia)
               do ir = 1,nr
                  new_SphMom(ir) = real(p_CDL%momL_Total(ir,1,ia),kind=RealKind)*Y0
                  old_SphMom(ir) = new_SphMom(ir)
               enddo
               call setValenceVPMomentSize(id,ia,new_SphMom(nr+1))
            endif
         enddo
!
      enddo
   endif
!
   end subroutine updateTotalDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSphericalChargeDensity()               result(isSphChrg)
!  ===================================================================
   implicit none
!
   logical :: isSphChrg
!
   isSphChrg = isSphericalCharge
!
   end function isSphericalChargeDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRhoLmax(id)                                   result(l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: l
!
   l = ChargeDensityList(id)%lmax
!
   end function getRhoLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPseudoRad(id)                               result(rad)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: rad
!
   rad = ChargeDensityList(id)%RcutPseudo
!
   end function getPseudoRad
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPseudoNumRPts(id)                            result(nr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: nr
!
   nr = ChargeDensityList(id)%NumRPts_Pseudo
!
   end function getPseudoNumRPts
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNeutralChargeDensity()                    result(rho_0)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: rho_0
!
   rho_0 = rho0_neutral
!
   end function getNeutralChargeDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine rtol_fit(id,ia)
!  ===================================================================
   use InterpolationModule, only : LeastSqFitInterp
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   integer (kind=IntKind) :: ir, ml, ll, jl, lmax, jmax, nRPts, im1_r0l
   integer (kind=IntKind), allocatable :: i_r0l(:)
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind) :: a,b, r0l, r0_max
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: rho_l(:,:)
   complex (kind=CmplxKind), allocatable :: rho_tmp(:)
!
   a = Ten2m2
   b = TWO
!
   lmax = ChargeDensityList(id)%lmax
   jmax = ChargeDensityList(id)%jmax
   if ( lmax==0 ) then
      return
   endif
!
   nRPts  =  ChargeDensityList(id)%NumRPts
   r0_max =  ChargeDensityList(id)%RcutPseudo
   r_mesh => ChargeDensityList(id)%r_mesh(1:nRPts)
   rho_l  => ChargeDensityList(id)%rhoL_Valence(1:nRPts,1:jmax,ia)
   flag_jl => ChargeDensityList(id)%ChargeCompFlag(1:jmax)
!
   if ( lmax==0 ) then
      return
   endif
!
   allocate( i_r0l(1:lmax) )
   allocate( rho_tmp(1:nRPts) )
!
   do ll = 1,lmax
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     For each l component of charge density set a radia from where 
!     the extrapolation to small r begins
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      r0l = a*(b**ll)
!     ----------------------------------------------------------------
      call hunt(nRPts,r_mesh(1:nRPts),r0l,ir)
!     ----------------------------------------------------------------
      if ( r_mesh(ir) <= r0_max ) then
         i_r0l(ll) = ir
         im1_r0l = ir-1
      else
         i_r0l(ll) = ChargeDensityList(id)%NumRPts_Pseudo
         im1_r0l = i_r0l(ll)-1
      endif
!      print *,"l =",ll,"R0l =", r_mesh(i_r0l(ll))
!
      LOOP_ML : do ml = 0,ll
         jl = (ll+1)*(ll+2)/2 - ll + ml
         if ( flag_jl(jl) == 0 ) then
            cycle LOOP_ML
         endif
         do ir = 1,nRPts
            rho_tmp(ir) = rho_l(ir,jl)/(r_mesh(ir)**ll)
         enddo
         call LeastSqFitInterp( 10, r_mesh(i_r0l(ll):nRPts),          &
                                rho_tmp(i_r0l(ll):nRPts), im1_r0l,    &
                                r_mesh(1:im1_r0l), rho_tmp(1:im1_r0l) )
         do ir = 1,im1_r0l
            rho_l(ir,jl) = rho_tmp(ir)*(r_mesh(ir)**ll)
         enddo
!
      enddo LOOP_ML
   enddo
!
   deallocate( i_r0l, rho_tmp )
!
   end subroutine rtol_fit
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDensity_r(denType, id, nrs, r_mesh, den_r_in)
!  ===================================================================
   use MPPModule, only : MyPE
   use MathParamModule, only: Ten2m8
!
   implicit none
!
   character (len=*), intent(in) :: denType
!
   integer (kind=IntKind), intent(in) :: id, nrs
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), intent(in), target  :: den_r_in(:)
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
   file_den(1:lenDenName) = trim(denType)//"_"
!
   if ( id == 1 ) then 
      write(file_den(lenDenName+1:lenDenName+6),'(i6)') offset+MyPE+id
      file_den(lenDenName+1:lenDenName+1) = 'n'
      funit = 55+MyPE+id
      open(unit=funit,file=trim(file_den),status='unknown')
      write(funit,'(a)') "#Ind     r_mesh   (lm):"
      write(funit,'(a)') "  "
      do ir = 1, nrs
         write(funit,'(i4,2(1x,d16.8))') id, r_mesh(ir), den_r_in(ir)
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
   subroutine printDensity_c(denType, id, nrs, r_mesh, den_c)
!  ===================================================================
   use MPPModule, only : MyPE
   use MathParamModule, only: Ten2m8
!
   implicit none
!
   character (len=*), intent(in) :: denType
!
   integer (kind=IntKind), intent(in) :: id, nrs
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), intent(in), target :: den_c(:)
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
   file_den(1:lenDenName) = trim(denType)//"_"
!
   if ( id == 1 ) then 
      write(file_den(lenDenName+1:lenDenName+6),'(i6)') offset+MyPE+id
      file_den(lenDenName+1:lenDenName+1) = 'n'
      funit = 55+MyPE+id
      open(unit=funit,file=trim(file_den),status='unknown')
      write(funit,'(a)') "#Ind     r_mesh   (lm):"
      write(funit,'(a)') "  "
      do ir = 1, nrs
         write(funit,'(i4,3(1x,d16.8))') id, r_mesh(ir), den_c(ir)
      enddo
      write(funit,'(a)') " "
      close(funit)
   endif
!
   end subroutine printDensity_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printMultipoleMom(id,ia)
!  ===================================================================
   use MPPModule, only : MyPE
   use MathParamModule, only: Ten2m8
   use SystemModule, only : getNumAlloyElements
   use Atom2ProcModule, only : getGlobalIndex
!
   implicit none
!
   integer (kind=IntKind), optional, intent(in) :: id, ia
!
   complex (kind=CmplxKind), pointer :: mmom(:,:)
!
   character(len=20) :: file_mm
   integer (kind=IntKind) :: jl, funit, ig, ja, l, m, lig
   integer (kind=IntKind) :: offset = 100000
   integer (kind=IntKind), allocatable :: flag_jl(:)
!
   if (MyPE /= 0) then
     return
   endif
   file_mm(1:20) = "                    "
   file_mm(1:13) = "MultipoleMom_"
   allocate( flag_jl(jmax_glb) )
   flag_jl(1:jmax_glb) = 0
!
   if ( present(id) ) then
      if (.not.present(ia)) then
         call ErrorHandler('printMultipoleMom','atom index on site is not specified')
      endif
      if (id < 1 .or. id > LocalNumAtoms) then
         call ErrorHandler('printMultipoleMom','local site index out of range',id)
      else if (ia < 1 .or. ia > ChargeDensityList(id)%NumSpecies) then
         call ErrorHandler('printMultipoleMom','atom index out of range',ia)
      endif
      ig = getGlobalIndex(id)
      write(file_mm(14:19),'(i6)') offset+ig*10+ia
      file_mm(14:14) = 'n'
      funit = 55+id
      open(unit=funit,file=trim(file_mm),status='unknown')
      write(funit,'(a,i4,a,i4,/,a)') "# Global Site ID::",ig,", Atom ID::",ia,"    jl   (lm):"
      write(funit,'(a)') "  "
      lig = MM_table_line(ig) + ia
      mmom => MultipoleMom
      do jl = 1,ChargeDensityList(id)%jmax
         if ( abs(mmom(lig,jl)) > Ten2m8) then
               flag_jl(jl) = 1
         endif
      enddo
      do jl = 1, ChargeDensityList(id)%jmax
         if ( flag_jl(jl) == 1 ) then
            l = lofj(jl)
            m = mofj(jl)
            write(funit,'(2i4,4x,2(1x,d16.8))') l, m, mmom(lig,jl)
         endif
      enddo
      write(funit,'(a)') " "
      close(funit)
   else if (present(ia)) then
      call ErrorHandler('printMultipoleMom','site index is not specified')
   else
      file_mm(14:19)="system"
      funit = 29
      open(unit=funit,file=trim(file_mm),status='unknown')
      mmom => MultipoleMom
      do jl = 1,jmax_glb
         FlagLOOP: do ig = 1,GlobalNumAtoms
            do ja = 1, getNumAlloyElements(ig)
               lig = MM_table_line(ig) + ja
               if ( abs(mmom(lig,jl)) > Ten2m8) then
                  flag_jl(jl) = 1
                  exit FlagLOOP
               endif
            enddo
         enddo FlagLOOP
      enddo
      write(funit,'(a,$)') "#AtomID(jl)::     "
      do jl = 1, jmax_glb
         if ( flag_jl(jl) == 1 ) then
            l = lofj(jl)
             m = mofj(jl)
             write(funit,'(a1,i2,1x,i2,a28,$)') "(",l, m, &
                           " )                          "
          endif
      enddo
      write(funit,'(a)') "  "
      do ig = 1, GlobalNumAtoms
         do ja = 1, getNumAlloyElements(ig)
            lig = MM_table_line(ig) + ja
            write(funit,'(i4,3x,i4,3x,$)') ig, ja
            do jl = 1,jmax_glb
               if (flag_jl(jl)==1) then
                  write(funit,'(2(1x,d16.8),$)') mmom(lig,jl)
               endif
            enddo
         enddo
         write(funit,'(a)') " "
      enddo
      close(funit)
   endif
   deallocate( flag_jl )
!
   end subroutine printMultipoleMom
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printChargeDensity_L( id, den_type, flag_rl )
!  ===================================================================
   use MPPModule, only : MyPE
   use MathParamModule, only: Ten2m8
!
   implicit none
!
   character (len=*), intent(in), optional :: den_type
!
   integer (kind=IntKind), intent(in), optional :: flag_rl
   integer (kind=IntKind), intent(in) :: id
!
   character(len=30) :: file_denl
   character(len=9)  :: char_lm
   integer (kind=IntKind) :: l, m, jl, jmax, ia, ir, NumRs
   integer (kind=IntKind) ::  funit, nc, len_name
   integer (kind=IntKind) :: offset = 100000
   integer (kind=IntKind) :: offset_at = 1000
   integer (kind=IntKind) :: offset_lm = 100
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: denl(:,:,:)
!
   if (isMTFP .and. nocaseCompare(den_type,"Pseudo") ) then
      call ErrorHandler('printChargeDensity_L','Invalid density type','Pseudo')
   endif
!
   file_denl(1:30) = " "
   if ( present(den_type) ) then
      len_name = len(den_type)
      if ( present(flag_rl) ) then
         if (flag_rl==0) then
            nc = 5+len_name
            file_denl(4+len_name:nc)  = "L_"
         else
            nc = 5+len_name+1
            file_denl(4+len_name:nc)  = "RL_"
         endif
      else
         nc = 5+len_name
         file_denl(4+len_name:nc)  = "L_"
      endif
      file_denl(1:3)  = "Rho"
      file_denl(4:3+len_name)  = den_type(1:len_name)
   else
      if ( present(flag_rl) ) then
         if (flag_rl==0) then
            nc = 5
            file_denl(1:nc)  = "RhoNL"
         else
            nc = 6
            file_denl(1:nc)  = "RhoNRL"
         endif
      else
         nc = 5
         file_denl(1:nc)  = "RhoNL"
      endif
   endif
   write(file_denl(nc+1:nc+6),'(i6)') offset+MyPE
   file_denl(nc+1:nc+1) = 'n'
   file_denl(nc+7:nc+7) = '_'
   write(file_denl(nc+8:nc+11),'(i4)') offset_at+id
   file_denl(nc+8:nc+8) = 'a'
   funit = 55+MyPE+id
   open(unit=funit,file=trim(file_denl),status='unknown')
!
   NumRs = ChargeDensityList(id)%NumRPts
   r_mesh => ChargeDensityList(id)%r_mesh(1:NumRs)
   jmax =  ChargeDensityList(id)%jmax
   if ( present(den_type) ) then
      if ( nocaseCompare(den_type,"Tilda") ) then
         denl => ChargeDensityList(id)%rhoL_Tilda
      else if ( nocaseCompare(den_type,"Pseudo") ) then
         denl => ChargeDensityList(id)%rhoL_Pseudo
      else if ( nocaseCompare(den_type,"Valence") ) then
         denl => ChargeDensityList(id)%rhoL_Valence
      else if ( nocaseCompare(den_type,"TotalNew") ) then
         denl => ChargeDensityList(id)%rhoL_Total
      else if ( nocaseCompare(den_type,"TotalOld") ) then
         denl => ChargeDensityList(id)%rhoL_TotalOld
      else
         denl => ChargeDensityList(id)%rhoL_Total
      endif
   else
      denl => ChargeDensityList(id)%rhoL_Total
   endif
!
   flag_jl => ChargeDensityList(id)%ChargeCompFlag
!
   write(funit,'(a,$)') "#Ind  Spec      r_mesh"
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
   do ia = 1, ChargeDensityList(id)%NumSpecies
      do ir = 1, NumRs
         write(funit,'(i4,1x,i4,1x,d16.8,$)') id, ia, r_mesh(ir)
         JlLoop: do jl = 1,jmax
            if ( flag_jl(jl) /= 0 ) then
               l = lofj(jl)
               if ( present(flag_rl) ) then
                  if (flag_rl==0) then
                      write(funit,'(1x,(2(1x,d16.8)),$)') denl(ir,jl,ia) ! /(r_mesh(ir)**l)
                  else
                      write(funit,'(1x,(2(1x,d16.8)),$)') denl(ir,jl,ia)/(r_mesh(ir)**l)
                  endif
               else
                  write(funit,'(1x,(2(1x,d16.8)),$)') denl(ir,jl,ia) ! /(r_mesh(ir)**l)
               endif
            endif
         enddo JlLoop
         write(funit,'(a)') " "
      enddo
   enddo
   nullify(flag_jl)
!
   close(funit)
!
   end subroutine printChargeDensity_L
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printMomentDensity_L( id, mom_type, flag_rl )
!  ===================================================================
   use MPPModule, only : MyPE
   use MathParamModule, only: Ten2m8
!
   implicit none
!
   character (len=*), intent(in), optional :: mom_type
!
   integer (kind=IntKind), intent(in), optional :: flag_rl
   integer (kind=IntKind), intent(in) :: id
!
   character(len=30) :: file_moml
   character(len=9)  :: char_lm
   integer (kind=IntKind) :: l, m, jl, jmax, ia, ir, NumRs
   integer (kind=IntKind) ::  funit, nc, len_name
   integer (kind=IntKind) :: offset = 100000
   integer (kind=IntKind) :: offset_at = 1000
   integer (kind=IntKind) :: offset_lm = 100
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: moml(:,:)
!
   if (isMTFP .and. nocaseCompare(mom_type,"Pseudo")) then
      call ErrorHandler('printMomentDensity_L','Invalid density type','Pseudo')
   endif
!
   file_moml(1:30) = " "
   if ( present(mom_type) ) then
      len_name = len(mom_type)
      if ( present(flag_rl) ) then
         if (flag_rl==0) then
            nc = 5+len_name
            file_moml(4+len_name:nc)  = "L_"
         else
            nc = 5+len_name+1
            file_moml(4+len_name:nc)  = "RL_"
         endif
      else
         nc = 5+len_name
         file_moml(4+len_name:nc)  = "L_"
      endif
      file_moml(1:3)  = "Mom"
      file_moml(4:3+len_name)  = mom_type(1:len_name)
   else
      if ( present(flag_rl) ) then
         if (flag_rl==0) then
            nc = 5
            file_moml(1:nc)  = "MomNL"
         else
            nc = 6
            file_moml(1:nc)  = "MomNRL"
         endif
      else
         nc = 5
         file_moml(1:nc)  = "MomNL"
      endif
   endif
   write(file_moml(nc+1:nc+6),'(i6)') offset+MyPE
   file_moml(nc+1:nc+1) = 'n'
   file_moml(nc+7:nc+7) = '_'
   write(file_moml(nc+8:nc+11),'(i4)') offset_at+id
   file_moml(nc+8:nc+8) = 'a'
   funit = 55+MyPE+id
   open(unit=funit,file=trim(file_moml),status='unknown')
!
   NumRs = ChargeDensityList(id)%NumRPts
   r_mesh => ChargeDensityList(id)%r_mesh(1:NumRs)
   jmax =  ChargeDensityList(id)%jmax
!
   flag_jl => ChargeDensityList(id)%ChargeCompFlag
!
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
!
   do ia = 1, ChargeDensityList(id)%NumSpecies
      if ( present(mom_type) ) then
         if ( nocaseCompare(mom_type,"Tilda") ) then
            moml => ChargeDensityList(id)%momL_Tilda(:,:,ia)
         else if ( nocaseCompare(mom_type,"Pseudo") ) then
            moml => ChargeDensityList(id)%momL_Pseudo(:,:,ia)
         else if ( nocaseCompare(mom_type,"Valence") ) then
            moml => ChargeDensityList(id)%momL_Valence(:,:,1,ia) ! spin-polarized case only
         else if ( nocaseCompare(mom_type,"TotalNew") ) then
            moml => ChargeDensityList(id)%momL_Total(:,:,ia)
         else if ( nocaseCompare(mom_type,"TotalOld") ) then
            moml => ChargeDensityList(id)%momL_TotalOld(:,:,ia)
         else
            moml => ChargeDensityList(id)%momL_Total(:,:,ia)
         endif
      else
         moml => ChargeDensityList(id)%momL_Total(:,:,ia)
      endif
      do ir = 1, NumRs
         write(funit,'(i4,1x,i4,1x,d16.8,$)') id, ia, r_mesh(ir)
         JlLoop: do jl = 1,jmax
            if ( flag_jl(jl) /= 0 ) then
               l = lofj(jl)
               if ( present(flag_rl) ) then
                  if (flag_rl==0) then
                      write(funit,'(1x,(2(1x,d16.8)),$)') moml(ir,jl) ! /(r_mesh(ir)**l)
                  else
                      write(funit,'(1x,(2(1x,d16.8)),$)') moml(ir,jl)/(r_mesh(ir)**l)
                  endif
               else
                  write(funit,'(1x,(2(1x,d16.8)),$)') moml(ir,jl) ! /(r_mesh(ir)**l)
               endif
            endif
         enddo JlLoop
         write(funit,'(a)') " "
      enddo
   enddo
   nullify(flag_jl)
!
   close(funit)
!
   end subroutine printMomentDensity_L
!  ==================================================================
!
end module ChargeDensityModule
