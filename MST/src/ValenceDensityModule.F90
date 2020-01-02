module ValenceDensityModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, PI, PI2, PI4, CZERO, CONE, TWO, HALF, SQRTm1, Y0
   use MathParamModule, only : TEN2m6, TEN2m5, TEN2m8, TEN2m4
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
   use IntegerFactorsModule, only : lofk, lofj, mofj, m1m, mofk
   use PublicTypeDefinitionsModule, only : GridStruct
   use TimerModule, only : getTime
   use MPPModule, only : MyPE, syncAllPEs
!
public :: initValenceDensity,        &
          endValenceDensity,         &
          printRho_L,                &
          printMom_L,                &
          printValenceDensity,       &
          getValenceVPCharge,        &
          setValenceVPCharge,        &
          getValenceMTCharge,        &
          setValenceMTCharge,        &
          getValenceVPMoment,        &
          setValenceVPMoment,        &
          getValenceMTMoment,        &
          setValenceMTMoment,        &
          getValenceVPMomentSize,    &
          setValenceVPMomentSize,    &
          getValenceMTMomentSize,    &
          getValenceElectronDensity, &
          setValenceElectronDensity, &
          getValenceElectronDensityAtPosi, &
          getValenceMomentDensity,   &
          setValenceMomentDensity,   &
          getValenceMomentDensityAtPosi, &
          getValenceMomentTorque,    &
          setValenceMomentTorque,    &
          getValenceSphericalElectronDensity,  &
          getValenceSphericalMomentDensity,  &
          getSphRho,                 &
          setSphRho,                 &
          getSphMom,                 &
          setSphMom,                 &
          getBandEnergy,             &
          getValenceKineticEnergy,   &
          getFermiEnergy,            &
          setFermiEnergy,            &
          getExchangeEnergy,         &  ! This is magnetic exchange energy
          getInterstitialValenceDensity,  &
          setInterstitialValenceDensity,  &
          getInterstitialValenceMoment,   &
          setInterstitialValenceMoment,   &
          updateValenceDensity,      &
          updateValenceEvec,         &
          updateValenceCharge,       &
          updateValenceEnergy,       &
          updateFermiEnergy,         &
          isDenComponentZero,        &
          getDenComponentFlag
!
   interface getValenceElectronDensity
      module procedure getved0, getved1, getved2
   end interface
!
   interface getValenceMomentDensity
      module procedure getvmd0, getvmd1, getvmd2, getvmd3
   end interface
!
   interface getSphRho
      module procedure getSphRho_is, getSphRho_tot, getSphRho_evec
   end interface
!
   interface getSphMom
      module procedure getSphMom_is, getSphMom_tot, getSphMom_evec
   end interface
!
   interface getBandEnergy
      module procedure getBandEnergy_0, getBandEnergy_1, getBandEnergy_2
   end interface
!
private
!
!  ===================================================================
!  lmax:     angular momentum l cutoff
!  jmax:     = (lmax+1)*(lmax+2)/2
!  rho_0:    spherical component of charge density data
!  mom_0:    spherical component of moment density data
!  Grid:     pointer to radial grid data
!  rho_l:    charge density data for 0 <= l <= lmax
!  mom_l:    moment density data for 0 <= l <= lmax
!  ===================================================================
   type ValenceDensityStruct
      integer (kind=IntKind) :: NumSpecies
      real (kind=RealKind)   :: NumValenceElectrons
      integer (kind=IntKind) :: NumRs ! the radial mesh size of charghe density
      integer (kind=IntKind) :: lmax
      integer (kind=IntKind) :: jmax
      integer (kind=IntKind), pointer :: DenCompFlag(:)
!
      real (kind=RealKind), pointer :: rho_0(:,:)
      real (kind=RealKind), pointer :: mom_0(:,:,:)
      complex (kind=CmplxKind), pointer :: rho_l(:,:,:)
      complex (kind=CmplxKind), pointer :: mom_l(:,:,:,:)
      real (kind=RealKind), pointer :: der_rho_0(:,:)
      real (kind=RealKind), pointer :: der_mom_0(:,:,:)
      complex (kind=CmplxKind), pointer :: der_rho_l(:,:,:)
      complex (kind=CmplxKind), pointer :: der_mom_l(:,:,:,:)
      real (kind=RealKind), allocatable :: eigensum(:,:)
      real (kind=RealKind), allocatable :: ChargeVP(:)
      real (kind=RealKind), allocatable :: ChargeMT(:)
      real (kind=RealKind), allocatable :: MomentVP(:,:)
      real (kind=RealKind), allocatable :: MomentMT(:,:)
      real (kind=RealKind), allocatable :: ExchangeEnergy(:)
!
      type (GridStruct), pointer :: Grid
   end type ValenceDensityStruct
!
   type (ValenceDensityStruct), allocatable, target :: Density(:)
!
   logical :: Initialized = .false.
   logical :: isValSymmOn = .false.
   logical :: rad_derivative = .false.
!
!  Flags to stop the code after the specific task is completed( input(external) flags )
!
   character (len = 50) :: stop_routine
!
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: NumVacancies
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind) :: node_print_level
   integer (kind=IntKind), allocatable :: print_level(:)
!
   integer (kind=IntKind), parameter :: n_inter = 5 ! order of polynomial
                                                    ! interpolation
!
   integer (kind=IntKind) :: SpinIndex
   integer (kind=IntKind) :: LocalAtomID
!
!  real (kind=RealKind), parameter :: rho_tol = TEN2m6
   real (kind=RealKind), parameter :: rho_tol = TEN2m4  ! This is changed back from 10^-6 to 10^-4 on 09-23-18
   real (kind=RealKind) :: FermiEnergy
   real (kind=RealKind) :: BandEnergy
   real (kind=RealKind) :: qvaltws
   real (kind=RealKind) :: qvaltmt
   real (kind=RealKind) :: rhoint
   real (kind=RealKind) :: mint
   real (kind=RealKind) :: mint_vec(3)
   real (kind=RealKind), allocatable, target :: rho0(:)
   real (kind=RealKind), allocatable, target :: mom0(:)
   real (kind=RealKind), allocatable, target :: NetValenceCharge(:)
   real (kind=RealKind), allocatable :: factorVacancy(:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initValenceDensity(na,vnum,lmax_rho_in,pola,cant,efermi, &
                                 istop,iprint,isGGA)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
!
   use OutputModule, only : getStandardOutputLevel
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   use RadialGridModule, only : getGrid, getNumRmesh
!
   use SystemModule, only : getNumAtoms, getNumVacancies, getLmaxMax
!
   use AtomModule, only : getLocalNumSpecies
!
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       setDataStorage2Value,  &
                                       setDataStorageLDA,     &
                                       isDataStorageExisting, &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark
!
   use ScfDataModule, only : isChargeSymm
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: lmax_rho_in(na)
   integer (kind=IntKind), intent(in) :: pola, cant
   integer (kind=IntKind), intent(in) :: iprint(na)
   integer (kind=IntKind) :: id, is, jend, iend, ns
   integer (kind=IntKind) :: jmax_rho, lmax_max
   integer (kind=IntKind), allocatable :: DataSize(:), DataSizeL(:)
!
   real (kind=RealKind), intent(in) :: vnum(na)
   real (kind=RealKind), intent(in) :: efermi
!
   logical, intent(in), optional :: isGGA
!
   type (GridStruct), pointer :: Grid
!
   if (pola == 1 .or. pola == 2) then
      n_spin_pola = pola
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initValenceDensity',                         &
                        'Invalid spin polarizing parameter',pola)
!     ----------------------------------------------------------------
   endif
!
   if (cant == 1 .or. cant == 2) then
      n_spin_cant = cant
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initValenceDensity',                         &
                        'Invalid spin canting parameter',cant)
!     ----------------------------------------------------------------
   endif
!
   if (n_spin_pola < n_spin_cant) then
!     ----------------------------------------------------------------
      call ErrorHandler('initValenceDensity','n_spin_cant > n_spin_pola')
!     ----------------------------------------------------------------
   endif
!
   stop_routine = istop
!
   if (present(isGGA)) then
      rad_derivative = isGGA
   else
      rad_derivative = .false.
   endif
!
   isValSymmOn = isChargeSymm();  ! isValSymmOn = .true.
   LocalNumAtoms = na
   GlobalNumAtoms = getNumAtoms()
   NumVacancies = getNumVacancies()
   if (GlobalNumAtoms == NumVacancies) then
      NumVacancies = 0
   endif
!
   allocate( Density(LocalNumAtoms) )
   allocate( print_level(LocalNumAtoms) )
   allocate( factorVacancy(LocalNumAtoms) )
!
   allocate( DataSize(LocalNumAtoms), DataSizeL(LocalNumAtoms) )
!
   lmax_max = 0
   do id = 1,LocalNumAtoms
      lmax_max=max(lmax_max, lmax_rho_in(id))
   enddo
   lmax_max = max(lmax_max,getLmaxMax())
!  ------------------------------------------------------------------
   call initIntegerFactors(lmax_max)
!  ------------------------------------------------------------------
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   do id = 1,LocalNumAtoms
      Grid => getGrid(id)
      jend = Grid%jend
      jmax_rho=(lmax_rho_in(id)+1)*(lmax_rho_in(id)+2)/2
      DataSize(id)  = jend*getLocalNumSpecies(id)
      DataSizeL(id) = jend*jmax_rho*getLocalNumSpecies(id)
   enddo
!
   node_print_level = getStandardOutputLevel()
   do id = 1, LocalNumAtoms
      print_level(id) = iprint(id)
   enddo
!
!  ==================================================================
!  generate Storage space for charge densities
!  ==================================================================
   if (.not.isDataStorageExisting('ValenceElectronDensity')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'ValenceElectronDensity',  &
                             DataSize,RealType)
!     ----------------------------------------------------------------
      call setDataStorage2Value('ValenceElectronDensity',ZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         Grid => getGrid(id)
         jend = Grid%jend
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'ValenceElectronDensity',jend)
!        -------------------------------------------------------------
      enddo
   endif
!
   if (n_spin_pola == 2) then
      if (.not.isDataStorageExisting('ValenceMomentDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'ValenceMomentDensity', &
                                DataSize*(2*n_spin_cant-1),RealType)
!        -------------------------------------------------------------
         call setDataStorage2Value('ValenceElectronDensity',ZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            Grid => getGrid(id)
            jend = Grid%jend
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'ValenceElectronDensity',jend)
!           ----------------------------------------------------------
         enddo
      endif
   endif
!
   if (.not.isDataStorageExisting('L-ValenceElectronDensity')) then
!     ----------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'L-ValenceElectronDensity',  &
                             DataSizeL,ComplexType)
!     ----------------------------------------------------------------
      call setDataStorage2Value('L-ValenceElectronDensity',CZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
         Grid => getGrid(id)
         jend = Grid%jend
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'L-ValenceElectronDensity',jend)
!        -------------------------------------------------------------
      enddo
   endif
!
   if (n_spin_pola == 2) then
      if (.not.isDataStorageExisting('L-ValenceMomentDensity')) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'L-ValenceMomentDensity', &
                               DataSizeL*(2*n_spin_cant-1),ComplexType)
!        -------------------------------------------------------------
         call setDataStorage2Value('L-ValenceElectronDensity',CZERO)
!        -------------------------------------------------------------
         do id = 1,LocalNumAtoms
            Grid => getGrid(id)
            jend = Grid%jend
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'L-ValenceElectronDensity',jend)
!           ----------------------------------------------------------
         enddo
      endif
   endif
!
   iend = 0
   do id = 1, LocalNumAtoms
      Grid => getGrid(id)
      jmax_rho=(lmax_rho_in(id)+1)*(lmax_rho_in(id)+2)/2
      jend = Grid%jend
      Density(id)%NumSpecies = getLocalNumSpecies(id)
      Density(id)%NumRs = jend
      ns = Density(id)%NumSpecies
!     ================================================================
      Density(id)%rho_l => getDataStorage( id,'L-ValenceElectronDensity',  &
                                           jend, jmax_rho, ns, ComplexMark )
      Density(id)%rho_l = CZERO
      allocate(Density(id)%DenCompFlag(jmax_rho))
      if (n_spin_pola == 2) then
         Density(id)%mom_l => getDataStorage( id,'L-ValenceMomentDensity', &
                                              jend, jmax_rho, 2*n_spin_cant-1, ns, ComplexMark )
         Density(id)%mom_l = CZERO
      endif
      allocate(Density(id)%eigensum(n_spin_pola,ns))
      allocate(Density(id)%ChargeVP(ns))
      allocate(Density(id)%ChargeMT(ns))
      allocate(Density(id)%MomentVP(3,ns))
      allocate(Density(id)%MomentMT(3,ns))
      allocate(Density(id)%ExchangeEnergy(ns))
      Density(id)%eigensum = ZERO
      Density(id)%ChargeVP = ZERO
      Density(id)%ChargeMT = ZERO
      Density(id)%MomentVP = ZERO
      Density(id)%MomentMT = ZERO
      Density(id)%ExchangeEnergy = ZERO
!     ================================================================
      Density(id)%rho_0 => getDataStorage( id,'ValenceElectronDensity', &
                                     jend, ns, RealMark )
      Density(id)%rho_0 = ZERO
      if (n_spin_pola == 2) then
         Density(id)%mom_0 => getDataStorage( id,'ValenceMomentDensity', &
                                              jend,2*n_spin_cant-1,ns,RealMark )
         Density(id)%mom_0 = ZERO
      endif
!     ================================================================
      if (rad_derivative) then
         allocate(Density(id)%der_rho_0(jend,ns),Density(id)%der_rho_l(jend,jmax_rho,ns))
         Density(id)%der_rho_0 = ZERO
         Density(id)%der_rho_l = CZERO
         if (n_spin_pola == 2) then
            allocate(Density(id)%der_mom_0(jend,2*n_spin_cant-1,ns))
            allocate(Density(id)%der_mom_l(jend,jmax_rho,2*n_spin_cant-1,ns))
            Density(id)%der_mom_0 = ZERO
            Density(id)%der_mom_l = CZERO
         endif
      endif
!     ================================================================
      Density(id)%Grid=>Grid
      Density(id)%lmax=lmax_rho_in(id)
      Density(id)%jmax=jmax_rho
      Density(id)%NumValenceElectrons = vnum(id)
!
      if ( vnum(id) /= 0 ) then
         factorVacancy(id)=ONE
      else
         factorVacancy(id)=ZERO
      endif
!
      iend=max(iend,jend)
   enddo
!
   deallocate( DataSize, DataSizeL )
!
   rhoint = ZERO
   mint = ZERO
   mint_vec(1:3) = ZERO
!
   FermiEnergy = efermi
!
   allocate(NetValenceCharge(GlobalNumAtoms))    ! Need to be moved to SystemModule
   NetValenceCharge(1:GlobalNumAtoms) = ZERO     ! Need to be moved to SystemModule
!
   allocate( rho0(iend), mom0(iend) )
!
   Initialized = .true.
!
   end subroutine initValenceDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endValenceDensity()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
!
   implicit none
   integer (kind=IntKind) :: id
!
   do id = 1, LocalNumAtoms
      nullify(Density(id)%rho_l)
      if (n_spin_pola == 2) then
         nullify(Density(id)%mom_l)
      endif
      nullify(Density(id)%rho_0)
      if (n_spin_pola == 2) then
         nullify(Density(id)%mom_0)
      endif
      nullify(Density(id)%Grid)
      deallocate(Density(id)%DenCompFlag)
      if (rad_derivative) then
         deallocate(Density(id)%der_rho_0, Density(id)%der_rho_l)
         if (n_spin_pola == 2) then
            deallocate(Density(id)%der_mom_0, Density(id)%der_mom_l)
         endif
      endif
   enddo
!
   deallocate( Density )
   deallocate(rho0, mom0, factorVacancy)
!  -------------------------------------------------------------------
   call endIntegerFactors()
!  -------------------------------------------------------------------
!
   deallocate(NetValenceCharge, print_level)
!
   Initialized = .false.
!
   end subroutine endValenceDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printValenceDensity()
!  ===================================================================
   use IntegrationModule, only : calIntegration
!
   use AtomModule, only : getLocalEvecNew, getLocalSpeciesContent
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use PolyhedraModule, only : getNumPlanes, getPlane
   use PolyhedraModule, only : getNumCorners, getCorner
!
   implicit none
!
   integer (kind=Intkind) :: id, is, ip, ic, nump, numc, j, ia
   integer (kind=Intkind) :: nr, klmax, jlmax
!
   real (kind=RealKind), pointer :: pp(:,:), pc(:,:)
   real (kind=RealKind) :: rho, mom, r
   real (kind=RealKind) :: intr(10000)
!
   real (kind=RealKind) :: q_VP, q_MT
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
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('printValenceDensity','valence module not initialized')
!     ----------------------------------------------------------------
   endif
!
   if (maxval(print_level) >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,23x,a)')'***********************************'
      write(6,'( 23x,a )')'* Output from printValenceDensity *'
      write(6,'(23x,a,/)')'***********************************'
!
      write(6,'(2x,a,t40,a,f18.11)')'Fermi Energy (Ryd)','=',FermiEnergy
      write(6,'(2x,a,t40,a,f18.11)')'GLOBAL VP Int[n(e)]/GlobalNumAtoms','=', &
                                    qvaltws/real(GlobalNumAtoms-NumVacancies,RealKind)
!     write(6,'(2x,''GLOBAL Int[e*n(e)]/GlobalNumAtoms'', t40,''='',f18.11)') &
!                                   BandEnergy/real(GlobalNumAtoms-NumVacancies,RealKind)
      do id =1, LocalNumAtoms
         nr = Density(id)%NumRs
         jlmax = Density(id)%jmax; klmax = (Density(id)%lmax+1)**2
         if (print_level(id) >= 0) then
            do ia = 1, Density(id)%NumSpecies
!
               write(6,'(/,2x,60(''-''))')
               write(6,'(2x,a,i4,a,i4)')'Local Atom Index =',id,', Species Index',ia
               do is=1,n_spin_pola
                  write(6,'(4x,a,i2,a,t40,a,f18.11)')                          &
                        'Spin =',is,', LOCAL Int[e*n(e)]','=',Density(id)%eigensum(is,ia)
               enddo
!              -------------------------------------------------------
               call calIntegration(Density(id)%Grid%jmt,Density(id)%Grid%r_mesh,    &
                                   Density(id)%rho_0(:,ia),intr,2)
!              -------------------------------------------------------
               write(6,'(4x,a,t40,a,f18.11)')'Integrated Electron Density in MT','=', &
                                             intr(Density(id)%Grid%jmt)*PI4
!              -------------------------------------------------------
               q_VP=getVolumeIntegration( id, nr, Density(id)%Grid%r_mesh(1:nr),&
                                          klmax, jlmax, 0, Density(id)%rho_l(1:nr,1:jlmax,ia), q_MT )
!              -------------------------------------------------------
               write(6,'(4x,a,t40,a,2f18.11)')                                  &
                     'Integrated Electron Density in VP','=',q_VP, q_MT
               write(6,'(4x,a,t40,a,3f18.11)')                                 &
                     'Integrated Moment Density in MT','=',Density(id)%MomentMT(:,ia)
               write(6,'(4x,a,t40,a,3f18.11)')                                 &
                     'Integrated Moment Density in VP','=',Density(id)%MomentVP(:,ia)
!              -------------------------------------------------------
               if (n_spin_pola == 1) then
                  write(6,'(/,15x,a)')                                         &
                          ' Electron Density on Selected Cell Boundary Points'
                  write(6,'(12x,55(''=''))')
                  write(6,'(12x,a)')                                           &
                      '      X         Y         Z         Electron Density'
                  write(6,'(12x,55(''-''))')
               else
                  write(6,'(/,10x,a)')                                         &
                          'Electron and Moment Density on Selected Cell Boundary Points'
                  write(6,'(1x,78(''=''))')
                  write(6,'(1x,a)')                                            &
                      '      X         Y         Z         Electron Density         Moment Density'
                  write(6,'(1x,78(''-''))')
               endif
!
               nump = getNumPlanes(id)
               pp => getPlane(id)
               do ip = 1, nump
                  r = sqrt(pp(1,ip)*pp(1,ip)+pp(2,ip)*pp(2,ip)+pp(3,ip)*pp(3,ip))
                  if (r <= Density(id)%Grid%rend) then
!                    rho = getValenceElectronDensityAtPosi(id,1,pp(1:3,ip))
                     rho = getValueAtPosi(Density(id)%Grid%r_mesh,pp(1:3,ip),Density(id)%rho_l(:,:,ia))
                     if (n_spin_pola == 1) then
                        write(6,'(12x,3f10.5,4x,d20.13)')pp(1:3,ip),rho
                     else
                        mom = getValenceMomentDensityAtPosi(id,1,pp(1:3,ip))
                        write(6,'(1x,3f10.5,2(4x,d20.13))')pp(1:3,ip),rho,mom
                     endif
                  endif
               enddo
!
               numc = getNumCorners(id)
               pc => getCorner(id)
               do ic = 1, numc
                  r = sqrt(pc(1,ic)*pc(1,ic)+pc(2,ic)*pc(2,ic)+pc(3,ic)*pc(3,ic))
                  if (r <= Density(id)%Grid%rend) then
!                    rho = getValenceElectronDensityAtPosi(id,1,pc(1:3,ic))
                     rho = getValueAtPosi(Density(id)%Grid%r_mesh,pc(1:3,ic),Density(id)%rho_l(:,:,ia))
                     if (n_spin_pola == 1) then
                        write(6,'(12x,3f10.5,4x,d20.13)')pc(1:3,ic),rho
                     else
                        mom = getValenceMomentDensityAtPosi(id,1,pc(1:3,ic))
!                       mom = getValueAtPosi(Density(id)%Grid%r_mesh,pc(1:3,ic),Density(id)%mom_l(:,:,:,ia))
                        write(6,'(1x,3f10.5,2(4x,d20.13))')pc(1:3,ic),rho,mom
                     endif
                  endif
               enddo
!
               if (n_spin_pola == 1) then
                  write(6,'(12x,55(''=''),/)')
               else
                  write(6,'(1x,78(''=''),/)')
               endif
            enddo
!
         endif
      enddo
      write(6,'(/,80(''-''))')
   endif
!
   end subroutine printValenceDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceVPCharge(id,ia) result(vc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind) :: vc
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceVPCharge','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceVPCharge','Invalid species index',ia)
   endif
!
   vc = Density(id)%ChargeVP(ia)
!
   end function getValenceVPCharge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setValenceVPCharge(id,ia,vc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind), intent(in) :: vc
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setValenceVPCharge','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('setValenceVPCharge','Invalid species index',ia)
   endif
!
   Density(id)%ChargeVP(ia) = vc
!
   end subroutine setValenceVPCharge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceMTCharge(id,ia) result(vc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind) :: vc
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceMTCharge','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceMTCharge','Invalid species index',ia)
   endif
!
   vc = Density(id)%ChargeMT(ia)
!
   end function getValenceMTCharge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setValenceMTCharge(id,ia,vc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind), intent(in) :: vc
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setValenceMTCharge','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('setValenceMTCharge','Invalid species index',ia)
   endif
!
   Density(id)%ChargeMT(ia) = vc
!
   end subroutine setValenceMTCharge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceVPMoment(id,ia) result(vm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind) :: vm(3)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceVPMoment','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceVPCharge','Invalid species index',ia)
   endif
!
   vm(:) = Density(id)%MomentVP(:,ia)
!
   end function getValenceVPMoment
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setValenceVPMoment(id,ia,vm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind), intent(in) :: vm(3)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setValenceVPMoment','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('setValenceVPMoment','Invalid species index',ia)
   endif
!
   Density(id)%MomentVP(:,ia) = vm(:)
!
   end subroutine setValenceVPMoment
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceVPMomentSize(id,ia) result(vm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind) :: vm
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceVPMomentSize','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceVPmomentSize','Invalid species index',ia)
   endif
!
   vm = sqrt(Density(id)%MomentVP(1,ia)**2+Density(id)%MomentVP(2,ia)**2    &
                                          +Density(id)%MomentVP(3,ia)**2)
!
   end function getValenceVPMomentSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setValenceVPMomentSize(id,ia,vm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind), intent(in) :: vm
   real (kind=RealKind) :: vm0
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setValenceVPMomentSize','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('setValenceVPMomentSize','Invalid species index',ia)
   endif
!
   vm0 = sqrt(Density(id)%MomentVP(1,ia)**2+Density(id)%MomentVP(2,ia)**2    &
                                           +Density(id)%MomentVP(3,ia)**2)
   if (vm0 < TEN2m8) then
      Density(id)%MomentVP(:,ia) = ZERO
      Density(id)%MomentVP(3,ia) = vm   ! Assuming z-polarization
   else
      vm0 = vm/vm0
      Density(id)%MomentVP(:,ia) = Density(id)%MomentVP(:,ia)*vm0   ! Assuming z-polarization
   endif
!
   end subroutine setValenceVPMomentSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceMTMoment(id,ia) result(vm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind) :: vm(3)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceMTMoment','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceMTMoment','Invalid species index',ia)
   endif
!
   vm(:) = Density(id)%MomentMT(:,ia)
!
   end function getValenceMTMoment
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setValenceMTMoment(id,ia,vm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind), intent(in) :: vm(3)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setValenceMTMoment','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('setValenceMTMoment','Invalid species index',ia)
   endif
!
   Density(id)%MomentMT(:,ia) = vm(:)
!
   end subroutine setValenceMTMoment
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceMTMomentSize(id,ia) result(vm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   real (kind=RealKind) :: vm
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceMTMomentSize','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceMTMomentSize','Invalid species index',ia)
   endif
!
   vm = sqrt(Density(id)%MomentMT(1,ia)**2+Density(id)%MomentMT(2,ia)**2    &
                                          +Density(id)%MomentMT(3,ia)**2)
!
   end function getValenceMTMomentSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFermiEnergy() result(ef)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: ef
!
   ef=FermiEnergy
   end function getFermiEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setFermiEnergy(ef)
!  ===================================================================
   implicit none
   real (kind=RealKind), intent(in) :: ef
!
   FermiEnergy = ef
!
   end subroutine setFermiEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphRho_evec(id,ia,is,evec,isDerivative) result(p_rho0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
   integer (kind=IntKind) :: iend, isig, js
!
   real (kind=RealKind), intent(in) :: evec(3)
   real (kind=RealKind), pointer :: p_rho0(:)
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphRho','Invalid atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSphRho','Invalid spin index',is)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   else if (n_spin_cant == 1) then
      call ErrorHandler('getSphRho','This is not spin-canted case')
   else if (abs(sqrt(evec(1)**2+evec(2)**2+evec(3)**2)-ONE) > TEN2m6) then
      call ErrorHandler('getSphRho','evec is not unit vector')
   endif
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphRho_evec','Invalid atom index',id)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getSphRho_evec',                             &
                        'ValenceDensityModule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
   isig = 3-2*is
   rho0(1:iend) = ZERO
   if (deriv) then
      do js = 1, 3
         rho0(1:iend) = rho0(1:iend) + evec(js)*Density(id)%der_mom_0(1:iend,js,ia)
      enddo
      rho0(1:iend) = (Density(id)%der_rho_0(1:iend,ia)+isig*rho0(1:iend))*HALF
   else
      do js = 1, 3
         rho0(1:iend) = rho0(1:iend) + evec(js)*Density(id)%mom_0(1:iend,js,ia)
      enddo
      rho0(1:iend) = (Density(id)%rho_0(1:iend,ia)+isig*rho0(1:iend))*HALF
   endif
   p_rho0=>rho0(:)
!
   end function getSphRho_evec
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphRho_is(id,ia,is,isDerivative) result(p_rho0)
!  ===================================================================
   use AtomModule, only : getLocalEvecOld
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
   integer (kind=IntKind) :: iend, isig, js
!
   real (kind=RealKind), pointer :: p_rho0(:)
   real (kind=RealKind) :: evec(3)
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphRho_is','Invalid atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSphRho_is','Invalid spin index',is)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphRho_is','Invalid atom index',id)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getSphRho_is',                               &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
   if (n_spin_pola == 1) then
      if (deriv) then
         p_rho0=>Density(id)%der_rho_0(:,ia)
      else
         p_rho0=>Density(id)%rho_0(:,ia)
      endif
   else
      if (n_spin_cant == 1) then
         isig = 3-2*is
         if (deriv) then
            rho0(1:iend)=( Density(id)%der_rho_0(1:iend,ia)           &
                           +isig*Density(id)%der_mom_0(1:iend,1,ia) )*HALF
         else
            rho0(1:iend)=( Density(id)%rho_0(1:iend,ia)               &
                           +isig*Density(id)%mom_0(1:iend,1,ia) )*HALF
         endif
      else
         evec(1:3) = getLocalEvecOld(id)
         isig = 3-2*is
         rho0(1:iend) = ZERO
         if (deriv) then
            do js = 1, 3
               rho0(1:iend) = rho0(1:iend) +                             &
                              evec(js)*Density(id)%der_mom_0(1:iend,js,ia)
            enddo
            rho0(1:iend)=(Density(id)%der_rho_0(1:iend,ia)+isig*rho0(1:iend))*HALF
         else
            do js = 1, 3
               rho0(1:iend) = rho0(1:iend) +                             &
                              evec(js)*Density(id)%mom_0(1:iend,js,ia)
            enddo
            rho0(1:iend)=(Density(id)%rho_0(1:iend,ia)+isig*rho0(1:iend))*HALF
         endif
      endif
      p_rho0=>rho0(:)
   endif
!
   end function getSphRho_is
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphRho_tot(id,ia,isDerivative) result(p_rho0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: iend
!
   real (kind=RealKind), pointer :: p_rho0(:)
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphRho','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getSphRho_tot',                              &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
!
   if (deriv) then
      p_rho0=>Density(id)%der_rho_0(:,ia)
   else
      p_rho0=>Density(id)%rho_0(:,ia)
   endif
!
   end function getSphRho_tot
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceSphericalElectronDensity(id,ia,isDerivative) result(p_rho0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: iend
!
   real (kind=RealKind), pointer :: p_rho0(:)
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValSphDen','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceSphericalElectronDensity',                 &
                        'ValenceDensityModule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
!
   if (deriv) then
      p_rho0=>Density(id)%der_rho_0(:,ia)
   else
      p_rho0=>Density(id)%rho_0(:,ia)
   endif
!
   end function getValenceSphericalElectronDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceSphericalMomentDensity(id,ia,isDerivative) result(p_mom0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: iend
!
   real (kind=RealKind), pointer :: p_mom0(:,:)
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValSphDen','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceSphericalMomentDensity',           &
                        'ValenceDensityModule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
!
   if (deriv) then
      p_mom0=>Density(id)%der_mom_0(:,:,ia)
   else
      p_mom0=>Density(id)%mom_0(:,:,ia)
   endif
!
   end function getValenceSphericalMomentDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphMom_tot(id,ia,isDerivative) result(p_mom0)
!  ===================================================================
   use AtomModule, only : getLocalEvecOld
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: iend, is
!
   real (kind=RealKind), pointer :: p_mom0(:)
   real (kind=RealKind) :: evec(3)
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphMom','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getSphMom_tot',                              &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
   if (n_spin_pola == 1) then
      mom0(1:iend)=ZERO
      p_mom0 => mom0(:)
   else if (n_spin_cant == 1) then
      if (deriv) then
         p_mom0 => Density(id)%der_mom_0(:,1,ia)
      else
         p_mom0 => Density(id)%mom_0(:,1,ia)
      endif
   else
      evec(1:3) = getLocalEvecOld(id)
      mom0(1:iend) = ZERO
      if (deriv) then
         do is = 1,3
            mom0(1:iend)=mom0(1:iend)+Density(id)%der_mom_0(1:iend,is,ia)*evec(is)
         enddo
      else
         do is = 1,3
            mom0(1:iend)=mom0(1:iend)+Density(id)%mom_0(1:iend,is,ia)*evec(is)
         enddo
      endif
      p_mom0 => mom0(:)
   endif
!
   end function getSphMom_tot
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphMom_is(id,ia,is,isDerivative) result(p_mom0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
   integer (kind=IntKind) :: iend
!
   real (kind=RealKind), pointer :: p_mom0(:)
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphMom','Invalid atom index',id)
   else if(is < 1 .or. is > 2*n_spin_cant-1) then
      call ErrorHandler('getSphMom','Invalid spin index',is)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   else if (n_spin_pola == 1) then
      return
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getSphMom_tot',                              &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
   if (deriv) then
      p_mom0=>Density(id)%der_mom_0(:,is,ia)
   else
      p_mom0=>Density(id)%mom_0(:,is,ia)
   endif
!
   end function getSphMom_is
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphMom_evec(id,ia,evec,isDerivative) result(p_mom0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: iend, is
!
   real (kind=RealKind), intent(in) :: evec(3)
   real (kind=RealKind), pointer :: p_mom0(:)
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphMom','Invalid atom index',id)
   else if (n_spin_pola == 1) then
      call ErrorHandler('getSphMom','This is not spin-canted case')
   else if (n_spin_cant == 1) then
      call ErrorHandler('getSphMom','This is not spin-canted case')
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   else if (abs(sqrt(evec(1)**2+evec(2)**2+evec(3)**2)-ONE) > TEN2m6) then
      call ErrorHandler('getSphMom','evec is not unit vector')
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getSphMom_tot',                              &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
   mom0(1:iend) = ZERO
   if (deriv) then
      do is = 1,3
         mom0(1:iend)=mom0(1:iend)+Density(id)%der_mom_0(1:iend,is,ia)*evec(is)
      enddo
   else
      do is = 1,3
         mom0(1:iend)=mom0(1:iend)+Density(id)%mom_0(1:iend,is,ia)*evec(is)
      enddo
   endif
   p_mom0 => mom0(:)
!
   end function getSphMom_evec
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSphRho(id,ia,in_rho0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: iend
!
   real (kind=RealKind), intent(in) :: in_rho0(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setSphRho','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   iend = Density(id)%Grid%jend
   Density(id)%rho_0(1:iend,ia)=in_rho0(1:iend)
   Density(id)%rho_l(1:iend,1,ia)=in_rho0(1:iend)/Y0
!
   end subroutine setSphRho
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSphMom(id,ia,is,in_mom0)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
   integer (kind=IntKind) :: iend
!
   real (kind=RealKind), intent(in) :: in_mom0(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setSphMom','Invalid atom index',id)
   else if(is < 1 .or. is > 2*n_spin_cant-1) then
      call ErrorHandler('setSphMom','Invalid spin index',is)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   iend = Density(id)%Grid%jend
   Density(id)%mom_0(1:iend,is,ia)=in_mom0(1:iend)
   Density(id)%mom_l(1:iend,1,is,ia)=in_mom0(1:iend)/Y0
!
   end subroutine setSphMom
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getved0(id,ia,jl,isDerivative) result(rho_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind) :: iend
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   complex (kind=CmplxKind), pointer :: rho_l(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceElectronDensity','Invalid atom index',id)
   else if(jl < 1 .or. jl > Density(id)%jmax) then
      call ErrorHandler('getValenceElectronDensity','Invalid jl index',jl)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceElectronDensity',                  &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
!
   if (deriv) then
      rho_l=>Density(id)%der_rho_l(:,jl,ia)
   else
      rho_l=>Density(id)%rho_l(:,jl,ia)
   endif
!
   end function getved0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getved1(id,ia,isDerivative) result(rho_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: iend, jmax
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   complex (kind=CmplxKind), pointer :: rho_l(:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceElectronDensity','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceElectronDensity',                  &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
   jmax = Density(id)%jmax
!
   if (deriv) then
      rho_l=>Density(id)%der_rho_l(:,:,ia)
   else
      rho_l=>Density(id)%rho_l(:,:,ia)
   endif
!
   end function getved1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getved2(id,isDerivative) result(rho_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: iend, jmax
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   complex (kind=CmplxKind), pointer :: rho_l(:,:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceElectronDensity','Invalid atom index',id)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceElectronDensity',                  &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   iend = Density(id)%Grid%jend
   jmax = Density(id)%jmax
!
   if (deriv) then
      rho_l=>Density(id)%der_rho_l(:,:,:)
   else
      rho_l=>Density(id)%rho_l(:,:,:)
   endif
!
   end function getved2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceElectronDensityAtPosi(title,id,ia,is,posi) result(rho)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   use InterpolationModule, only : PolyInterp
!
   implicit none
!
   character (len=*), intent(in) :: title
   integer (kind=IntKind), intent(in) :: id, ia, is
   integer (kind=IntKind) :: jl, ir, irp, l, kl, iend, kmax
!
   real (kind=RealKind), intent(in) :: posi(3)
   real (kind=RealKind) :: rho, mom, err, r
!
   complex (kind=CmplxKind) :: rho_in
   complex (kind=CmplxKind), pointer :: rho_l(:,:), mom_l(:,:)
   complex (kind=CmplxKind), allocatable :: ylm(:)
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
!
      function getValueAtPosi(r_mesh,posi,val_l) result(val)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         real (kind=RealKind), intent(in) :: r_mesh(:), posi(3)
         real (kind=RealKind) :: val
         complex (kind=CmplxKind), intent(in) :: val_l(:,:)
      end function getValueAtPosi
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('getValenceElectronDensityAtPosi',                 &
                        'Valence Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceElectronDensityAtPosi','Invalid atom index',id)
   else if (is < 1 .or. is > 2) then
      call ErrorHandler('getSphRho','Invalid spin index',is)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getSphRho','Invalid species index',ia)
   endif
!
   rho = getValueAtPosi(Density(id)%Grid%r_mesh,posi,Density(id)%rho_l(:,:,ia))
   if (n_spin_pola == 2) then  ! Needs some work for spin-canting case
      mom = getValueAtPosi(Density(id)%Grid%r_mesh,posi,Density(id)%mom_l(:,:,1,ia))
      rho = HALF*(rho + (3-2*is)*mom)
   endif
   return
!
   r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
   if (r > Density(id)%Grid%rend) then
!      call WarningHandler('getValenceElectronDensityAtPosi','r > rend',   &
!                          r,Density(id)%Grid%rend)
      rho = ZERO
      return
   endif
!
   kmax = (Density(id)%lmax+1)**2
   allocate( ylm(kmax) )
!  -------------------------------------------------------------------
   call calYlm(posi,Density(id)%lmax,ylm)
!  -------------------------------------------------------------------
!
   iend = Density(id)%Grid%jend
!  -------------------------------------------------------------------
   call hunt(iend,Density(id)%Grid%r_mesh(1:iend),r,ir)
!  -------------------------------------------------------------------
   if (ir > iend-(n_inter-1)/2) then
      irp=iend-n_inter+1
   else if (2*ir+1 > n_inter) then
      irp=ir-(n_inter-1)/2
   else
      irp=1
   endif
!
   rho_l=>Density(id)%rho_l(:,:,ia)
!
   rho = ZERO
!
   do jl = 1,Density(id)%jmax
      l = lofj(jl)
!     ----------------------------------------------------------------
      call PolyInterp(n_inter, Density(id)%Grid%r_mesh(irp:irp+n_inter-1), &
                      rho_l(irp:irp+n_inter-1,jl), r, rho_in, err)
!     ----------------------------------------------------------------
      kl = (l+1)*(l+1)-l+mofj(jl)
      if (mofj(jl) == 0) then
         rho = rho + real(rho_in*ylm(kl),RealKind)
      else
         rho = rho + TWO*real(rho_in*ylm(kl),RealKind)
      endif
   enddo
!
   deallocate( ylm )
!
   end function getValenceElectronDensityAtPosi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getvmd0(id,ia,is,jl,isDerivative) result(mom_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: jl
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   complex (kind=CmplxKind), pointer :: mom_l(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceMomentDensity','Invalid atom index',id)
   else if(is < 1 .or. is > 2*n_spin_pola-1) then
      call ErrorHandler('getValenceMomentDensity','Invalid spin index',is)
   else if(jl < 1 .or. jl > Density(id)%jmax) then
      call ErrorHandler('getValenceMomentDensity','Invalid jl index',jl)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceMomentDensity','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceMomentDensity',                  &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   if (deriv) then
      mom_l=>Density(id)%der_mom_l(:,jl,is,ia)
   else
      mom_l=>Density(id)%mom_l(:,jl,is,ia)
   endif
!
   end function getvmd0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getvmd1(id,ia,is,isDerivative) result(mom_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in) :: is
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   complex (kind=CmplxKind), pointer :: mom_l(:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceMomentDensity','Invalid atom index',id)
   else if(is < 1 .or. is > 2*n_spin_pola-1) then
      call ErrorHandler('getValenceMomentDensity','Invalid spin index',is)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceMomentDensity','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceMomentDensity',                  &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   if (deriv) then
      mom_l=>Density(id)%der_mom_l(:,:,is,ia)
   else
      mom_l=>Density(id)%mom_l(:,:,is,ia)
   endif
!
   end function getvmd1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getvmd2(id,ia,isDerivative) result(mom_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   complex (kind=CmplxKind), pointer :: mom_l(:,:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceMomentDensity','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceMomentDensity','Invalid species index',ia)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceMomentDensity',                  &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   if (deriv) then
      mom_l=>Density(id)%der_mom_l(:,:,:,ia)
   else
      mom_l=>Density(id)%mom_l(:,:,:,ia)
   endif
!
   end function getvmd2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getvmd3(id,isDerivative) result(mom_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   logical, intent(in), optional :: isDerivative
   logical :: deriv
!
   complex (kind=CmplxKind), pointer :: mom_l(:,:,:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceMomentDensity','Invalid atom index',id)
   endif
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv .and. .not.rad_derivative) then
      call ErrorHandler('getValenceMomentDensity',                  &
                        'ValenceDensityMOdule is not initialized with GGA enabled')
   endif
!
   if (deriv) then
      mom_l=>Density(id)%der_mom_l(:,:,:,:)
   else
      mom_l=>Density(id)%mom_l(:,:,:,:)
   endif
!
   end function getvmd3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceMomentDensityAtPosi(id,ia,posi) result(mom)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   use InterpolationModule, only : PolyInterp
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: jl, ir, irp, l, kl, is, iend, kmax
!
   real (kind=RealKind), intent(in) :: posi(3)
   real (kind=RealKind) :: momv(3),mom, err, r, x
!
   complex (kind=CmplxKind) :: momv_in(3)
   complex (kind=CmplxKind), pointer :: mom_l(:,:,:)
   complex (kind=CmplxKind), allocatable :: ylm(:)
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
      call ErrorHandler('getValenceMomentDensityAtPosi',              &
                        'Valence Module is not initialized')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceMomentDensityAtPosi','Invalid atom index',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceMomentDensity','Invalid species index',ia)
   else if (n_spin_pola == 1) then
      mom = ZERO
      return
   endif
!
   r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
   if (r > Density(id)%Grid%rend) then
!     call WarningHandler('getValenceMomentDensityAtPosi','r > rend', &
!                         r,Density(id)%Grid%rend)
      mom = 0
      return
   endif
   x = log(r)
!
   kmax = (Density(id)%lmax+1)**2
   allocate( ylm(kmax) )
!  -------------------------------------------------------------------
   call calYlm(posi,Density(id)%lmax,ylm)
!  -------------------------------------------------------------------
!
   iend = Density(id)%Grid%jend
!  -------------------------------------------------------------------
   call hunt(iend,Density(id)%Grid%r_mesh(1:iend),r,ir)
!  -------------------------------------------------------------------
   if (ir > iend-(n_inter-1)/2) then
      irp=iend-n_inter+1
   else if (2*ir+1 > n_inter) then
      irp=ir-(n_inter-1)/2
   else
      irp=1
   endif
!
   mom_l=>Density(id)%mom_l(1:iend,1:Density(id)%jmax,1:2*n_spin_cant-1,ia)
!
   momv = ZERO
!
   do is = 1, 2*n_spin_cant-1
!     ----------------------------------------------------------------
!      call PolyInterp(n_inter, Density(id)%Grid%r_mesh(irp:irp+n_inter-1), &
!                   Density(id)%mom_0(irp:irp+n_inter-1,is,ia), r, momv(is), err)
!     ----------------------------------------------------------------
      do jl = 1,Density(id)%jmax
         l = lofj(jl)
!        -------------------------------------------------------------
         call PolyInterp(n_inter,                                    &
                         Density(id)%Grid%r_mesh(irp:irp+n_inter-1), &
                         mom_l(irp:irp+n_inter-1,jl,is), r,          &
                         momv_in(is), err)
!        -------------------------------------------------------------
         kl = (l+1)*(l+1)-l+mofj(jl)
         if (mofj(jl) == 0) then
            momv(is) = momv(is) + real(momv_in(is)*ylm(kl),RealKind)
         else
            momv(is) = momv(is) + TWO*real(momv_in(is)*ylm(kl),RealKind)
         endif
      enddo
   enddo
!
   deallocate( ylm )
!
   if (n_spin_cant == 1) then
      mom = abs(momv(1))
   else
      mom = sqrt(momv(1)*momv(1)+momv(2)*momv(2)+momv(3)*momv(3))
   endif
!
   end function getValenceMomentDensityAtPosi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setValenceElectronDensity(id,ia,jl,rho_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind) :: iend
!
   complex (kind=CmplxKind), intent(in) :: rho_l(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setValenElectronDensity','Invalid atom index',id)
   else if(jl < 1 .or. jl > Density(id)%jmax) then
      call ErrorHandler('setValenElectronDensity','Invalid jl index',jl)
   else if(ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('setValenMomentDensity','Invalid species index',ia)
   endif
!
   iend = Density(id)%Grid%jend
   if (jl == 1) then
      Density(id)%rho_0(1:iend,ia)=rho_l(1:iend)*Y0
   endif
   Density(id)%rho_l(1:iend,jl,ia)=rho_l(1:iend)
!
   end subroutine setValenceElectronDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setValenceMomentDensity(id,ia,is,jl,mom_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind) :: iend
!
   complex (kind=CmplxKind), intent(in) :: mom_l(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setValenMomentDensity','Invalid atom index',id)
   else if(is < 1 .or. is > 2*n_spin_cant-1) then
      call ErrorHandler('setValenMomentDensity','Invalid spin index',is)
   else if(jl < 1 .or. jl > Density(id)%jmax) then
      call ErrorHandler('setValenMomentDensity','Invalid jl index',jl)
   else if(ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('setValenMomentDensity','Invalid species index',ia)
   endif
!
   iend = Density(id)%Grid%jend
   if (jl == 1) then
      Density(id)%mom_0(1:iend,is,ia)=mom_l(1:iend)
   endif
   Density(id)%mom_l(1:iend,jl,is,ia)=mom_l(1:iend)
!
   end subroutine setValenceMomentDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRhoLmax(id) result(lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: lmax
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getRhoLmax','Invalid atom index',id)
   endif
!
   lmax = Density(id)%lmax
   end function getRhoLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBandEnergy_0() result(ev)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: ev
!
   ev = BandEnergy
!
   end function getBandEnergy_0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBandEnergy_1(id,is) result(ev)
!  ===================================================================
   use AtomModule, only : getLocalSpeciesContent
   implicit none
   integer (kind=IntKind), intent(in) :: id, is
   integer (kind=IntKind) :: ia
   real (kind=RealKind) :: ev
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getBandEnergy','Invalid atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getBandEnergy','Invalid spin index',is)
   endif
!
   ev = ZERO
   do ia = 1, Density(id)%NumSpecies
      ev = ev + getLocalSpeciesContent(id,ia)*Density(id)%eigensum(is,ia)
   enddo
!
   end function getBandEnergy_1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBandEnergy_2(id,ia,is) result(ev)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, ia, is
   real (kind=RealKind) :: ev
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getBandEnergy','Invalid atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getBandEnergy','Invalid spin index',is)
   else if(ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getBandEnergy','Invalid species index',ia)
   endif
!
   ev = Density(id)%eigensum(is,ia)
!
   end function getBandEnergy_2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getValenceKineticEnergy(id,ia,is) result(ke)
!  ===================================================================
   use PotentialModule, only : getPotential, getTruncatedPotential
   use PotentialModule, only : getPotLmax, getTruncPotLmax
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use AtomModule, only : getLocalEvecOld
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is
   integer (kind=IntKind) :: jend, ir, jl, nr, ns
   integer (kind=IntKind) :: lmax_rho, jmax_rho, lmax_pot, jmax_pot
   integer (kind=IntKind) :: lmax_prod, jmax_prod, kmax_prod
!
   real (kind=RealKind) :: ke, rhovint, vint_mt, evec(3)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: rho_val(:,:)
   complex (kind=CmplxKind), pointer :: mom_val(:,:,:)
   complex (kind=CmplxKind), pointer :: v_old(:,:)
   complex (kind=CmplxKind), allocatable :: rho_tmp(:,:), v_tmp(:,:), prod(:,:)
!
   type (GridStruct), pointer :: Grid
!
   interface
      subroutine computeProdExpan(n,lf,f,lg,g,lh,h)
        use KindParamModule, only : IntKind, CmplxKind
        integer (kind=IntKind), intent(in) :: n, lf, lg, lh
        complex (kind=CmplxKind), intent(in) :: f(:,:), g(:,:)
        complex (kind=CmplxKind), intent(out) :: h(:,:)
      end subroutine computeProdExpan
   end interface
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getValenceKineticEnergy','Invalid atom index',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getValenceKineticEnergy','Invalid spin index',is)
   else if(ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getValenceKineticEnergy','Invalid species index',ia)
   endif
!
   lmax_rho = Density(id)%lmax
   jmax_rho = Density(id)%jmax
   rho_val=>Density(id)%rho_l(:,:,ia)
   Grid => Density(id)%Grid
   jend = Grid%jend
   r_mesh => Grid%r_mesh(1:jend)
!
   if (n_spin_pola == 2) then
      mom_val => Density(id)%mom_l(:,:,:,ia)
   endif
!
   lmax_pot = getPotLmax(id)
   jmax_pot = (lmax_pot+1)*(lmax_pot+2)/2
   lmax_prod = lmax_rho+lmax_pot
   jmax_prod = (lmax_prod+1)*(lmax_prod+2)/2
   kmax_prod = (lmax_prod+1)**2
   allocate( rho_tmp(jend,jmax_rho), prod(jend,jmax_prod) )
!
   if ( n_spin_pola == 1 ) then
      do jl = 1,jmax_rho
         do ir = 1,jend
            rho_tmp(ir,jl) = rho_val(ir,jl)
         enddo
      enddo
   else if (n_spin_cant == 1) then
      do jl = 1,jmax_rho
         do ir = 1,jend
            rho_tmp(ir,jl) = (rho_val(ir,jl) + (3-2*is)*mom_val(ir,jl,1))*HALF
         enddo
      enddo
   else  
      evec = getLocalEvecOld(id)
      ns = 3-2*is
      do jl = 1,jmax_rho
         do ir = 1,jend
            rho_tmp(ir,jl) = ( rho_val(ir,jl) +                                          &
                               ns*(  evec(1)*mom_val(ir,jl,1) + evec(2)*mom_val(ir,jl,2) &
                                   + evec(3)*mom_val(ir,jl,3) ) )*HALF
         enddo
      enddo
   endif
!
   v_old => getPotential(id,ia,is)
!  -------------------------------------------------------------------
   call computeProdExpan(jend, lmax_rho, rho_tmp, lmax_pot, v_old, lmax_prod, prod)
   rhovint = getVolumeIntegration( id, jend, r_mesh, kmax_prod, jmax_prod, 0, prod, vint_mt)
!  -------------------------------------------------------------------
   ke = Density(id)%eigensum(is,ia) - rhovint
!
   deallocate(rho_tmp, prod)
!
   end function getValenceKineticEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterstitialValenceDensity() result(r0)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: r0
!
   r0 = rhoint
!
   end function getInterstitialValenceDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterstitialValenceMoment() result(m0)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: m0
!
   m0 = mint
!
   end function getInterstitialValenceMoment
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalNetValenceChargeTable() result(p_tble)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: p_tble(:)
!
   p_tble => NetValenceCharge(1:GlobalNumAtoms)
!
   end function getGlobalNetValenceChargeTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchangeEnergy(id,ia) result(exc_en)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind) :: exc_en
!
   if ( id < 1 .or. id > LocalNumAtoms ) then
      call ErrorHandler('getExchangeEnergy','Local atom index is out of range',id)
   else if (ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('getExchangeEnergy','Invalid species index',ia)
   endif
!
   exc_en = Density(id)%ExchangeEnergy(ia)
!
   end function getExchangeEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateValenceCharge(id,ia,dos,dos_mt)
!  ===================================================================
   use AtomModule, only : getLocalEvecNew
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), intent(in) :: dos(*), dos_mt(*)
   real (kind=RealKind) :: evec(3)
!
   if ( id < 1 .or. id > LocalNumAtoms ) then
      call ErrorHandler('updateValenceCharge','Local atom index is out of range',id)
   else if(ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('updateValenceCharge','Invalid species index',ia)
   endif
!
   evec = getLocalEvecNew(id)
   if (n_spin_cant == 2) then
      Density(id)%ChargeVP(ia) = dos(1)
      Density(id)%ChargeMT(ia) = dos_mt(1)
      Density(id)%MomentVP(1:3,ia) = dos(2:4)
      Density(id)%MomentMT(1:3,ia) = dos_mt(2:4)
   else if(n_spin_pola == 2) then
      Density(id)%ChargeVP(ia) = dos(1) + dos(2)
      Density(id)%ChargeMT(ia) = dos_mt(1) + dos_mt(2)
      Density(id)%MomentVP(:,ia) = ZERO
      Density(id)%MomentMT(:,ia) = ZERO
      Density(id)%MomentVP(3,ia) = dos(1) - dos(2)
      Density(id)%MomentMT(3,ia) = dos_mt(1) - dos_mt(2)
   else
      Density(id)%ChargeVP(ia) = dos(1)
      Density(id)%ChargeMT(ia) = dos_mt(1)
      Density(id)%MomentVP(:,ia) = ZERO
      Density(id)%MomentMT(:,ia) = ZERO
   endif
   if (print_level(id) >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,23x,a)')'***********************************'
      write(6,'( 23x,a )')'* Output from updateValenceCharge *'
      write(6,'(23x,a,/)')'***********************************'
!
      write(6,'(''id = '',i3,'', ia = '',i3,'', Charge_MT = '',f22.16)') id,ia,Density(id)%ChargeMT(ia)
      write(6,'(20x,''Charge_VP = '',f22.16)') Density(id)%ChargeVP(ia)
      if ( n_spin_pola/=1 ) then
         write(6,'(10x,''Moment_MT ='',3f22.16)') Density(id)%MomentMT(:,ia)
         write(6,'(10x,''Moment_VP ='',3f22.16)') Density(id)%MomentVP(:,ia)
      endif
   endif
!
   end subroutine updateValenceCharge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateValenceEnergy(id,ia,evalsum,exc)
!  ===================================================================
   use AtomModule, only : getLocalEvecOld
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind), intent(in) :: evalsum(:),exc
   real (kind=RealKind) :: ev0, evm
   real (kind=RealKind) :: evec(3)
!
   if ( id < 1 .or. id > LocalNumAtoms ) then
      call ErrorHandler('updateValenceEnergy','Local atom index is out of range',id)
   else if(ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('updateValenceEnergy','Invalid species index',ia)
   endif
!
   evec = getLocalEvecOld(id) ! Should use EvecOld, not EvecNew
   if (n_spin_cant == 2) then
      ev0 = evalsum(1)
      evm = evalsum(2)*evec(1)+evalsum(3)*evec(2)+evalsum(4)*evec(3)
      Density(id)%eigensum(1,ia) = HALF*(ev0+evm)
      Density(id)%eigensum(2,ia) = HALF*(ev0-evm)
   else if(n_spin_pola == 2) then
      Density(id)%eigensum(1,ia) = evalsum(1)
      Density(id)%eigensum(2,ia) = evalsum(2)
   else
      Density(id)%eigensum(1,ia) = evalsum(1)
   endif
!
   Density(id)%ExchangeEnergy(ia) = exc
!
   end subroutine updateValenceEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateFermiEnergy(ef,zvaltss)
!  ===================================================================
!
!  Note: This routine should be called after updateValenceEnergy!!!
!
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use SystemVolumeModule, only : getSystemVolume, getTotalInterstitialVolume
!
   use ScfDataModule, only : Harris
!
   use AtomModule, only : getLocalEvecNew, getLocalSpeciesContent
!
   use PotentialTypeModule, only : isASAPotential, isFullPotential,   &
                                   isMuffinTinFullPotential, isMuffinTinPotential
!
   use PolyhedraModule, only : getVolume, getInscrSphVolume
!
   implicit none
!
   integer (kind=IntKind) :: id, ig, i, iend, ia
!
   real (kind=RealKind), intent(in) :: ef, zvaltss
   real (kind=RealKind) :: qvalws
   real (kind=RealKind) :: qvalmt
   real (kind=RealKind) :: qint, zvaltss_avr, cfac, charge
   real (kind=RealKind) :: evec(3)
   real (kind=RealKind) :: wspace(8)
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), parameter :: qtol = TEN2m6
!
   character (len=17) :: sname='updateFermiEnergy'
!
   FermiEnergy = ef
!
!  ===================================================================
!  Also update other global quantities......
!  ===================================================================
!  zvaltss_avr = zvaltss/real(GlobalNumAtoms-NumVacancies,RealKind)
   zvaltss_avr = zvaltss/real(GlobalNumAtoms,RealKind)
   qvalws = ZERO
   qvalmt = ZERO
   mint_vec(1:3) = ZERO 
   BandEnergy = ZERO
   do id = 1, LocalNumAtoms
      evec(1:3) = getLocalEvecNew(id)
      do ia = 1, Density(id)%NumSpecies
         cfac = getLocalSpeciesContent(id,ia)
         if(n_spin_pola == 2) then
            BandEnergy = BandEnergy + cfac*(Density(id)%eigensum(1,ia) + Density(id)%eigensum(2,ia))
         else
            BandEnergy = BandEnergy + cfac*(Density(id)%eigensum(1,ia))
         endif
!        qvalws = qvalws + (Density(id)%ChargeVP-factorVacancy(id)*zvaltss_avr)
!        qvalmt = qvalmt + (Density(id)%ChargeMT-factorVacancy(id)*zvaltss_avr)
         qvalws = qvalws + cfac*(Density(id)%ChargeVP(ia)-zvaltss_avr)
         qvalmt = qvalmt + cfac*(Density(id)%ChargeMT(ia)-zvaltss_avr)
         mint_vec(1:3) = mint_vec(1:3)                                &
                       + cfac*(Density(id)%MomentVP(1:3,ia) - Density(id)%MomentMT(1:3,ia))
      enddo
   enddo
!           
   wspace(1)=BandEnergy
   wspace(2)=qvalws
   wspace(3)=qvalmt
   wspace(4:6)=mint_vec(1:3)
   if(isASAPotential()) then
      wspace(7) = (wspace(2)-wspace(3))/(getSystemVolume()*n_spin_pola)
   else  
      wspace(7) = (wspace(2)-wspace(3))/getTotalInterstitialVolume()
   endif 
   wspace(8) = (wspace(2)-wspace(3))
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,wspace,8)
!  -------------------------------------------------------------------
!        
   BandEnergy=wspace(1)
   qvaltws=wspace(2)+zvaltss
   qvaltmt=wspace(3)+zvaltss
   mint_vec(1:3)=wspace(4:6)
   rhoint = wspace(7)
!     
!  ===================================================================
!  check that average integrated valence DOS = average valence.....
!  ===================================================================
   if (abs(qvaltws-zvaltss) > qtol) then
      if (node_print_level >= 0) then
         write(6,'(/,''updateFermiEnergy:: Trouble: VP[q,z]:'',2d16.8)')qvaltws,zvaltss
      endif
   endif
!
   qint = wspace(8)
   if (abs(qint) < TEN2m8) then
      qint = ZERO
   else if (qint < ZERO) then
      call WarningHandler('updateFermiEnergy','negative qint',qint)
   endif
!
   if(isASAPotential()) then
       if (rhoint > TEN2m8) then
          do id = 1, LocalNumAtoms
             iend = Density(id)%Grid%jend
             r_mesh => Density(id)%Grid%r_mesh
             do ia = 1, Density(id)%NumSpecies
                do i=1, iend
                   Density(id)%rho_0(i,ia) = Density(id)%rho_0(i,ia) + rhoint
                   Density(id)%rho_l(i,1,ia) = Density(id)%rho_l(i,1,ia) + sqrt(PI4)*rhoint
                enddo
             enddo
          enddo
       endif
       qint = ZERO
       rhoint = ZERO
       mint = ZERO
       mint_vec = ZERO
   else
       if(n_spin_pola.eq.2) then
           mint=sqrt(mint_vec(1)*mint_vec(1)+mint_vec(2)*mint_vec(2)+mint_vec(3)*mint_vec(3))
       else
           mint=ZERO
       endif
       if(mint < TEN2m5) then
           mint=ZERO
           mint_vec=ZERO
       endif
   endif
!
   do id = 1, LocalNumAtoms
      if(print_level(id) >= 0) then
          write(6,'(/,a,t40,a,f22.16)')'Number of Electrons outside MT-sphere','=',qint
          write(6,'(a,t40,a,f22.16)')'Magnetic moment outside MT-sphere','=',mint
          write(6,'(a,t40,a,f22.16)')'Electron density outside MT-sphere','=',rhoint
      endif
   enddo
!
   NetValenceCharge(1:GlobalNumAtoms) = ZERO
   do id = 1, LocalNumAtoms
       ig = getGlobalIndex(id)
       charge = ZERO
       do ia = 1, Density(id)%NumSpecies
          cfac = getLocalSpeciesContent(id,ia)
          if (isFullPotential()) then
             charge = charge + cfac*Density(id)%ChargeVP(ia)
          else
             charge = charge + cfac*Density(id)%ChargeMT(ia)
          endif
      enddo
      if (isMuffinTinPotential() .or. isMuffinTinFullPotential())  then
          NetValenceCharge(ig) = Density(id)%NumValenceElectrons-charge  &
                                -rhoint*(getVolume(id)-getInscrSphVolume(id))
      else
          NetValenceCharge(ig) = Density(id)%NumValenceElectrons-charge
      endif
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,NetValenceCharge,GlobalNumAtoms)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Add [N*FermiEnergy] to band energy: fixed chem. potl calc...........
!  ===================================================================
   if(Harris > 1) then
      BandEnergy=BandEnergy + HALF*(qvaltws-zvaltss)*FermiEnergy
   endif
!
   if(stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine updateFermiEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateValenceDensity(id,ia,dos_r_jl,isDerivative)
!  ===================================================================
!  This code needs modified to allow for dos_r_jl has an extra
!  dimention for multi chemical species
   use MPPModule, only : MyPE, syncAllPEs
!
   use SystemSymmetryModule, only : getSymmetryFlags
!
   implicit none
!
   character (len=20) :: sname='updateValenceDensity'
!
   logical :: isZERO, deriv
   logical, intent(in), optional :: isDerivative
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind) :: is, ir, jl, NumRs, izamax
   integer (kind=IntKind) :: jmax_rho
   integer (kind=IntKind), pointer :: flags_jl(:)
!
   real (kind=RealKind) :: rfac, rho_r, rho_i, mom_r, mom_i
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), pointer :: rho_0(:), mom_0(:,:)
!
   complex (kind=CmplxKind), intent(in) :: dos_r_jl(:,:,:)
   complex (kind=CmplxKind), pointer :: rho_l(:,:), mom_l(:,:,:)
!
   if ( id < 1 .or. id > LocalNumAtoms ) then
      call ErrorHandler('updateValenceDensity','Local atom index is out of range',id)
   else if(ia < 1 .or. ia > Density(id)%NumSpecies) then
      call ErrorHandler('updateValenceDensity','Invalid species index',ia)
   endif
!
   NumRs = Density(id)%NumRs
   r_mesh => Density(id)%Grid%r_mesh(1:NumRs)
   jmax_rho = Density(id)%jmax
!
   if (present(isDerivative)) then
      deriv = isDerivative
   else
      deriv = .false.
   endif
!
   if (deriv) then
      rho_0 => Density(id)%der_rho_0(:,ia)
      rho_l => Density(id)%der_rho_l(:,:,ia)
      if (n_spin_pola == 2) then
         mom_0 => Density(id)%der_mom_0(:,:,ia)
         mom_l => Density(id)%der_mom_l(:,:,:,ia)
      endif
   else
      rho_0 => Density(id)%rho_0(:,ia)
      rho_l => Density(id)%rho_l(:,:,ia)
      if (n_spin_pola == 2) then
         mom_0 => Density(id)%mom_0(:,:,ia)
         mom_l => Density(id)%mom_l(:,:,:,ia)
      endif
   endif
!
   rfac = ONE/(sqrt(PI4))
   if (n_spin_pola == 1 .or. n_spin_cant == 2) then
      do ir = 1, NumRs
         rho_0(ir) = real(dos_r_jl(ir,1,1),kind=RealKind)*rfac/(r_mesh(ir)*r_mesh(ir))
      enddo
   else
      do ir = 1, NumRs
         rho_0(ir) = real(dos_r_jl(ir,1,1)+dos_r_jl(ir,1,2),kind=RealKind)*rfac/(r_mesh(ir)*r_mesh(ir))
      enddo
   endif
!
   do jl = 1, jmax_rho
      if (n_spin_pola == 1 .or. n_spin_cant == 2) then
         do ir = 1, NumRs
            rho_l(ir,jl) = dos_r_jl(ir,jl,1)/(r_mesh(ir)*r_mesh(ir))
         enddo
      else
         do ir = 1, NumRs
            rho_l(ir,jl) = (dos_r_jl(ir,jl,1)+dos_r_jl(ir,jl,2))/(r_mesh(ir)*r_mesh(ir))
         enddo
      endif
!
      isZERO = .true.
      LOOP_ir: do ir = 1, NumRs
         if (abs(rho_l(ir,jl)) > rho_tol) then
            isZERO = .false.
            exit LOOP_ir
         endif
      enddo LOOP_ir
      if (isZERO) then
!        if (getSymmetryFlags(id,jl) > 0) then
!           call WarningHandler('updateValenceDensity','Non-zero valence density in this (l,m) channel is expected',lofj(jl),mofj(jl))
!        endif
         Density(id)%DenCompFlag(jl) = 0
      else
         Density(id)%DenCompFlag(jl) = 1
         if (getSymmetryFlags(id,jl) == 0) then
            call WarningHandler('updateValenceDensity','Zero valence density in this (l,m) channel is expected', &
                                lofj(jl),mofj(jl))
            ir = izamax(NumRs,rho_l(:,jl),1)
            write(6,'(a,3(i5,2x),2d16.8)')'Maximum density value in this channel = ',MyPE,ir,jl,rho_l(ir,jl)
            if ( isValSymmOn ) then
               Density(id)%DenCompFlag(jl) = 0
            endif
         endif
      endif
      if (Density(id)%DenCompFlag(jl) == 0) then
         do ir = 1, NumRs
            rho_l(ir,jl) = CZERO
         enddo
      endif
   enddo
!
   if (n_spin_cant == 2) then
      mom_l = CZERO
      do is=1,3
         do ir = 1, NumRs
            mom_0(ir,is) = real(dos_r_jl(ir,1,is+1),kind=RealKind)*rfac/(r_mesh(ir)*r_mesh(ir))
         enddo
         do jl = 1, jmax_rho
            if (Density(id)%DenCompFlag(jl) > 0) then
               do ir = 1, NumRs
                  mom_l(ir,jl,is) = dos_r_jl(ir,jl,is+1)/(r_mesh(ir)*r_mesh(ir))
               enddo
            endif
         enddo
      enddo
   else if (n_spin_pola == 2) then
      mom_l = CZERO
      do ir = 1, NumRs
         mom_0(ir,1) = real(dos_r_jl(ir,1,1)-dos_r_jl(ir,1,2),kind=RealKind)*rfac/(r_mesh(ir)*r_mesh(ir))
      enddo
      do jl = 1, jmax_rho
         if (Density(id)%DenCompFlag(jl) > 0) then
            do ir = 1, NumRs
               mom_l(ir,jl,1) = (dos_r_jl(ir,jl,1)-dos_r_jl(ir,jl,2))/(r_mesh(ir)*r_mesh(ir))
            enddo
         endif
      enddo
   endif
!
   if ( isValSymmOn ) then
      flags_jl => getSymmetryFlags(id)
      do jl = 1, jmax_rho
         if ( flags_jl(jl) ==0 ) then
            rho_l(1:NumRs,jl) = CZERO
            Density(id)%DenCompFlag(jl) = 0
         else if ( flags_jl(jl) == 1 ) then
            do ir = 1,NumRs
               rho_r = real(rho_l(ir,jl), kind=RealKind )
               rho_l(ir,jl) = cmplx(rho_r,ZERO,kind=CmplxKind)
            enddo
         else if ( flags_jl(jl) == 2 ) then
            do ir = 1,NumRs
               rho_i = real(-sqrtm1*rho_l(ir,jl), kind=RealKind )
               rho_l(ir,jl) = cmplx(ZERO, rho_i, kind=CmplxKind)
            enddo
         endif
      enddo
      if (n_spin_pola == 2) then
         do is = 1, 2*n_spin_cant-1
            do jl = 1, jmax_rho
               if ( flags_jl(jl) ==0 ) then
                  mom_l(1:NumRs,jl,is) = CZERO
               else if ( flags_jl(jl) == 1 ) then
                  do ir = 1,NumRs
                     mom_r = real( mom_l(ir,jl,is), kind=RealKind )
                     mom_l(ir,jl,is) = cmplx(mom_r,ZERO,kind=CmplxKind)
                  enddo
               else if ( flags_jl(jl) == 2 ) then
                  do ir = 1,NumRs
                     mom_i = real(-sqrtm1*mom_l(ir,jl,is), kind=RealKind)
                     mom_l(ir,jl,is) = cmplx(ZERO, mom_i, kind=CmplxKind)
                  enddo
               endif
            enddo
         enddo
      endif
   endif
!
   if(stop_routine == sname) then
      call StopHandler(sname)
   endif
!
   end subroutine updateValenceDensity
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printRho_L(id,ia,aux_name)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: id, ia
!
   character(len=30) :: file_rhol
   character(len=*), intent(in), optional :: aux_name
   integer (kind=IntKind) :: i, l, m, jl, jmax, ir, NumRs, funit, count_flag
   integer (kind=IntKind) :: offset = 100000, alen, nc, na, nb, ja
   integer (kind=IntKind), allocatable :: flag_jl(:)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: rhol(:,:)
!
   if (MyPE /= 0) then
     return
   endif
   alen=len(file_rhol)
   file_rhol(1:alen) = " "
   file_rhol(1:5) = "RhoL_"
   if (present(aux_name)) then
      alen = len(trim(aux_name))
      file_rhol(6:6+alen-1) = aux_name(1:alen)
      nc = 6 + alen
   else
      nc = 6
   endif
!
   if (present(id)) then
      na = id; nb = id
   else
      na = 1; nb = LocalNumAtoms
   endif
!
   if (present(ia)) then
      ja = ia
   else
      ja = 1
   endif
!
   if ( na < 1 .or. na > LocalNumAtoms ) then
      call ErrorHandler('printRho_L','Local atom index is out of range',na)
   else if(ja < 1 .or. ja > Density(id)%NumSpecies) then
      call ErrorHandler('printRho_L','Invalid species index',ja)
   endif
!
   do i = na, nb
      write(file_rhol(nc:nc+5),'(i6)') offset+MyPE+i
      file_rhol(nc:nc) = 'n'
      funit = 55+MyPE+i
      open(unit=funit,file=trim(file_rhol),status='unknown')
      NumRs = Density(i)%NumRs
      r_mesh => Density(i)%Grid%r_mesh(1:NumRs)
      rhol => Density(i)%rho_l(:,:,ja)
      jmax =  Density(i)%jmax
      allocate (flag_jl(jmax))
      flag_jl = 0
      count_flag = 0
      do jl = 1,jmax
         FlagLOOP: do ir = 1,NumRs
            if ( abs(rhol(ir,jl)) > TEN2m8) then
               flag_jl(jl) = 1
               count_flag = count_flag+1
               exit FlagLOOP
            endif
         enddo FlagLOOP
      enddo
      write(funit,'(a,$)') " Ind     r_mesh   (lm):"
      do jl = 1,jmax
         if ( flag_jl(jl) > 0 ) then
            l = lofj(jl)
            m = mofj(jl)
            write(funit,'(11x,a2,i2,1x,i2,a2,11x,$)') "( ", l, m," )"
         endif
!        write(funit,'(i4,$)') flag_jl(j)
      enddo
      write(funit,'(a)') "  "
      do ir = 1, NumRs
         write(funit,'(i4,1x,d16.8,$)') i, r_mesh(ir)
         JlLoop: do jl = 1,jmax
            if ( flag_jl(jl) > 0 ) then
               write(funit,'(1x,(2(1x,d16.8)),$)') rhol(ir,jl)
            endif
         enddo JlLoop
         write(funit,'(a)') " "
      enddo
      deallocate(flag_jl)
      close(funit)
   enddo
!
   end subroutine printRho_L
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printMom_L(id,ia,aux_name)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) , optional :: id, ia
!
   character(len=*), intent(in), optional :: aux_name
   character(len=11) :: file_moml
   integer (kind=IntKind) :: i, j, l, m, is, jl, jmax, ir, NumRs, funit, count_flag
   integer (kind=IntKind) :: offset = 100000, alen, nc, na, nb, ja
   integer (kind=IntKind), allocatable :: flag_jl(:)
!
   real (kind=RealKind), pointer :: r_mesh(:), moml0(:,:)
!
   complex (kind=CmplxKind), pointer :: moml(:,:,:)
!
   if (MyPE /= 0) then
     return
   endif
   alen=len(file_moml)
   file_moml(1:alen) = " "
   file_moml(1:5) = "MomL_"
   if (present(aux_name)) then
      alen = len(trim(aux_name))
      file_moml(6:6+alen-1) = aux_name(1:alen)
      nc = 6 + alen
   else
      nc = 6
   endif
!
   if (present(id)) then
      na = id; nb = id
   else
      na = 1; nb = LocalNumAtoms
   endif
!
   if (present(ia)) then
      ja = ia
   else
      ja = 1
   endif
!
   if ( na < 1 .or. na > LocalNumAtoms ) then
      call ErrorHandler('printMom_L','Local atom index is out of range',na)
   else if(ja < 1 .or. ja > Density(id)%NumSpecies) then
      call ErrorHandler('printMom_L','Invalid species index',ja)
   endif
!
   do i = na, nb
      write(file_moml(nc:nc+5),'(i6)') offset+MyPE+i
      file_moml(nc:nc) = 'n'
      funit = 55+MyPE+i
      open(unit=funit,file=trim(file_moml),status='unknown')
      NumRs = Density(i)%NumRs
      r_mesh => Density(i)%Grid%r_mesh(1:NumRs)
      moml => Density(i)%mom_l(:,:,:,ja)
      moml0 => Density(i)%mom_0(:,:,ja)
      jmax =  Density(i)%jmax
      allocate (flag_jl(jmax))
      flag_jl = 0
      count_flag = 0
      do is = 1,2*n_spin_cant-1
         do jl = 1,jmax
            FlagLOOP: do ir = 1,NumRs
               if ( abs(moml(ir,jl,is)) > TEN2m8) then
                  flag_jl(jl) = 1
                  count_flag = count_flag+1
                  exit FlagLOOP
               endif
            enddo FlagLOOP
         enddo
      enddo
      write(funit,'(a,$)') " Ind     r_mesh   (lm):"
      do j = 1,jmax
         if ( flag_jl(jl) > 0 ) then
            l = lofj(jl)
            m = mofj(jl)
            write(funit,'(11x,a2,i2,1x,i2,a2,11x,$)') "( ", l, m," )"
         endif
!        write(funit,'(i4,$)') flag_jl(j)
      enddo
      write(funit,'(a)') "  "
      do ir = 1, NumRs
         write(funit,'(i4,1x,d16.8,$)') i, r_mesh(ir)
         do is = 1,2*n_spin_cant-1
            JlLoop: do jl = 1,jmax
               if ( flag_jl(jl) > 0 ) then
                  write(funit,'(1x,(2(1x,d16.8)),$)') moml(ir,jl,is)
               endif
            enddo JlLoop
         enddo
         write(funit,'(a)') " "
      enddo
      deallocate(flag_jl)
      close(funit)
   enddo
!
   end subroutine printMom_L
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateValenceEvec(id,torque)
!  ===================================================================
   use Atom2ProcModule, only : getGlobalIndex
!
   use SystemModule, only : setMomentDirection, setMomentDirectionOld
!
   use AtomModule, only : getMixingParam4Evec, getLocalEvecOld,       &
                          setLocalEvecNew, setLocalMoment,            &
                          getLocalEvecOut, getLocalSpeciesContent
!
   use ScfDataModule, only : getCantedMomentTorqueFactor
!
   use PotentialTypeModule, only : isFullPotential
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), intent(in) :: torque(3)
!
   integer (kind=IntKind) :: i, ig, ia
!
   real (kind=RealKind) :: evec_new(3)
   real (kind=RealKind) :: moment(3)
   real (kind=RealKind) :: alpev
   real (kind=RealKind) :: ctq
   real (kind=RealKind) :: evec_old(3)
   real (kind=RealKind) :: emag
   real (kind=RealKind), parameter :: mtol = TEN2m8
!
   if ( id < 1 .or. id > LocalNumAtoms ) then
      call ErrorHandler('updateValenceEvec','Local atom index is out of range',id)
   else if (n_spin_cant == 1) then
      return
   endif
!
   ig = getGlobalIndex(id)
   alpev = getMixingParam4Evec(id)
   ctq = getCantedMomentTorqueFactor()
   evec_old(1:3) = getLocalEvecOld(id)
!
   moment = ZERO
   do ia = 1, Density(id)%NumSpecies
      if ( .not.isFullPotential() ) then
         moment = moment + getLocalSpeciesContent(id,ia)*getValenceMTMoment(id,ia)
      else
         moment = moment + getLocalSpeciesContent(id,ia)*getValenceVPMoment(id,ia)
      endif
   enddo
!
   call setLocalMoment(id,moment)
!
!  ===================================================================
!  determine evec_new according to moment..........................
!  ===================================================================
   emag=sqrt(moment(1)*moment(1)+moment(2)*moment(2)+moment(3)*moment(3))
   if(node_print_level.ge.0) then
       write(6,'(/,a)')'updateValenceEvec:'
       write(6,'(4x,a,t40,''='',3f11.6)')'magnetic moment vector',    &
            (moment(i),i=1,3)
       write(6,'(4x,a,t40,''='',1f11.6)')'magnetic moment magnitude', &
            emag
   endif
   if(emag > mtol) then
!     ================================================================
!        evec is the new moment orientation:
!             evec(1) = the x-component of e-vector
!             evec(2) = the y-component of e-vector
!             evec(3) = the z-component of e-vector
!        it is determined by the total moment inside the muffin-tin
!        sphare
!     ================================================================
      evec_new(1:3)=moment(1:3)/emag
   else
      evec_new(1:3)=evec_old(1:3)
   endif
!
!  ===================================================================
!     add precession motion to evec_new...............................
!     Note, an arbitary constant ctq is multiplied to torque..........
!  ===================================================================
   evec_new(1:3)=evec_new(1:3)+ctq*torque(1:3)
   emag=sqrt(evec_new(1)*evec_new(1)+evec_new(2)*evec_new(2)          &
             +evec_new(3)*evec_new(3))
!
   if(emag < mtol) then
      call ErrorHandler('getevec','emag = 0',emag)
   endif
!
   evec_new(1:3)=evec_new(1:3)/emag
!
   if(node_print_level.ge.0) then
      write(6,'(4x,a,t40,''='',3f11.6)')                              &
            'Moment direc. old          ',evec_old(1),evec_old(2),evec_old(3)
      write(6,'(4x,a,t40,''='',3f11.6)')                              &
            'Moment direc. before mixing',evec_new(1),evec_new(2),evec_new(3)
   endif
!
!  ===================================================================
!     perform simple mixing of evec_new and evec_old, and redefine
!     evec_new........................................................
!  ===================================================================
   evec_new(1:3)=alpev*evec_new(1:3)+(one-alpev)*evec_old(1:3)
   emag=sqrt(evec_new(1)*evec_new(1)+evec_new(2)*evec_new(2)          &
             +evec_new(3)*evec_new(3))
!
   if(emag < mtol) then
      call ErrorHandler('getevec','emag = 0',emag)
   endif
!
   evec_new(1:3)=evec_new(1:3)/emag
!
!  -------------------------------------------------------------------
   call setLocalEvecNew(id,evec_new)
!  -------------------------------------------------------------------
   if(node_print_level.ge.0) then
      write(6,'(4x,a,t40,''='',3f11.6)')                              &
            'Moment direc. new        ',evec_new(1),evec_new(2),evec_new(3)
   endif
!  -------------------------------------------------------------------
   evec_new(1:3) = getLocalEvecOut(id)
!  -------------------------------------------------------------------
   call setMomentDirection(ig,evec_new)
!  -------------------------------------------------------------------
   call setMomentDirectionOld(ig,evec_old)
!  -------------------------------------------------------------------
!
   if(node_print_level.ge.0) then
      write(6,'(4x,a,t40,''='',3f11.6)')                              &
            'New moment orientation',evec_new(1),evec_new(2),evec_new(3)
   endif
!
   end subroutine updateValenceEvec
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isDenComponentZero(id,is,jl) result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is
   integer (kind=IntKind), intent(in) :: jl
!
   logical :: flag
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('isDenComponentZero','invalid id',id)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('isDenComponentZero','invalid spin index',is)
   else if (jl < 1 .or. jl > Density(id)%jmax) then
      flag = .true.
   else if (Density(id)%DenCompFlag(jl) == 0) then
      flag = .true.
   else
      flag = .false.
   endif
!
   end function isDenComponentZero
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDenComponentFlag(id)                    result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer, pointer :: flag(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDenComponentZero','invalid id',id)
   endif
!
   flag => Density(id)%DenCompFlag(:)
!
   end function getDenComponentFlag
!  ===================================================================
end module ValenceDensityModule
