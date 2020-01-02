!  *******************************************************************
!  * Note: VPCellElectronTables should be the integration of total   *
!  *       electron densidy over VP-cell plus the amount of the core *
!  *       electrons of the neighboring atoms which are extended into*
!  *       the local cell.                                           *
!  *       In the present version of this code, VPCellElectronTables *
!  *       takes the value fed in from the update function. The value*
!  *       in fact is the integration of the valence electron density*
!  *       over VP-cell plus the untruncated total core electrons.   *
!  *******************************************************************
module ChargeDistributionModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
   use MathParamModule, only : ZERO, HALF, ONE, TWO, PI2, PI4, THIRD, TEN2m8
!
public :: initChargeDistribution,           &
          endChargeDistribution,            &
          updateChargeDistribution,         &
          getGlobalTableLine,               &
          getGlobalOnSiteElectronTable,     &   ! Qmt_i + rho_0*Omega0_i
          getGlobalMTSphereElectronTable,   &   ! Qmt_i
          getGlobalVPCellElectronTable,     &   ! Qvp_i
          getGlobalNetMomentTable,          &
          getGlobalMTSphereMomentTable,     &
          getGlobalVPCellMomentTable,       &
          getGlobalExchangeEnergyTable,     &
          getInterstitialElectronDensity,   &   ! rho_0
          getInterstitialMomentDensity,     &
          getGlobalOnSiteElectronTableOld,  &   ! Qmt_i + rho_0*Omega0_i
          getGlobalMTSphereElectronTableOld,&   ! Qmt_i
          getGlobalVPCellElectronTableOld,  &   ! Qvp_i
          getInterstitialElectronDensityOld,&   ! rho_0
          getAverageMoment,                 &
          printChargeDistribution
!
private
   logical :: Initialized = .false.
   logical :: Updated = .false.
!
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: MaxLocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind), allocatable, target :: global_table_line(:) ! Atom address in the global table
   integer (kind=IntKind), allocatable, target :: local_array_line(:)  ! Atom address in the local array
   integer (kind=IntKind) :: global_table_size, local_array_size
!
   real (kind=RealKind) :: mom_aver = ZERO
   real (kind=RealKind) :: mdenint
   real (kind=RealKind), allocatable, target :: Table_Wkspace(:,:)
   real (kind=RealKind), pointer :: NetMomentTable(:)
   real (kind=RealKind), pointer :: MTSphereMomentTable(:)
   real (kind=RealKind), pointer :: VPCellMomentTable(:)
   real (kind=RealKind), pointer :: ExchangeEnergyTable(:)
!
   real (kind=RealKind) :: rhoint
   real (kind=RealKind), pointer :: OnSiteElectronTable(:)
   real (kind=RealKind), pointer :: MTSphereElectronTable(:)
   real (kind=RealKind), pointer :: VPCellElectronTable(:)
!
   real (kind=RealKind) :: rhoint_old
   real (kind=RealKind), pointer :: OnSiteElectronTableOld(:)
   real (kind=RealKind), pointer :: MTSphereElectronTableOld(:)
   real (kind=RealKind), pointer :: VPCellElectronTableOld(:)
!
!  ExEn here is the exchange splitting energy between spin-ip and spin down states
   real (kind=RealKind), allocatable :: Vint(:), Qmt(:), Qvp(:), Mmt(:), Mvp(:), ExEn(:)
   real (kind=RealKind), allocatable :: Qmt_old(:), Qvp_old(:)
   real (kind=RealKind), allocatable :: memtemp(:,:), q_mix(:)
!
   type EMDStruct
      real (kind=RealKind) :: VPCharge
      real (kind=RealKind) :: VPMoment
      real (kind=RealKind), pointer :: rho0(:)
      real (kind=RealKind), pointer :: mom0(:)
   end type EMDStruct
!
   type ChargeListStruct
      integer (kind=IntKind) :: jmt
      integer (kind=IntKind) :: NumRs
      integer (kind=IntKind) :: NumSpecies
      real (kind=RealKind) :: rmt
      real (kind=RealKind), pointer :: r_mesh(:)
      type (EMDStruct), pointer :: EMD(:)
      type (ChargeListStruct), pointer :: next
   end type ChargeListStruct
!
   type (ChargeListStruct), target :: ChargeList
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initChargeDistribution(na,nt,ns)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalSumInGroup
   use Atom2ProcModule, only : getMaxLocalNumAtoms, getGlobalIndex
   use Atom2ProcModule, only : getLocalIndex, getAtom2ProcInGroup
   use AtomModule, only : getMixingParam4Charge, getLocalNumSpecies
   use AtomModule, only : getLocalSpeciesContent
   use StepFunctionModule, only : getVolumeIntegration
   use PolyhedraModule, only : getVolume
   use SystemModule, only : getAtomicNumber, getNumAlloyElements
   use PotentialTypeModule, only : isASAPotential
   implicit none
!
   integer (kind=IntKind), intent(in) :: nt
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: ns
!
   type (ChargeListStruct), pointer :: p_DL
!
   integer (kind=IntKind) :: ig, ip, iq, id, jmt, NumRs
   integer (kind=IntKind) :: lig, lid, ia
!
   real (kind=RealKind) :: rmt, corr, qint_old, GVint
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), allocatable :: rho0(:)
   real (kind=RealKind) :: sums(2)
!
   GlobalNumAtoms = nt
   LocalNumAtoms = na
   MaxLocalNumAtoms = getMaxLocalNumAtoms()
   n_spin_pola = ns
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   allocate(global_table_line(GlobalNumAtoms))
   global_table_line = 0
   global_table_size = 0
   do ig = 1, GlobalNumAtoms
      global_table_line(ig) = global_table_size
      global_table_size = global_table_size + getNumAlloyElements(ig)
   enddo
!
   if ( n_spin_pola==2 ) then
      allocate( Table_Wkspace(1:global_table_size,10) )
      OnSiteElectronTableOld   => Table_Wkspace(1:global_table_size,1)
      MTSphereElectronTableOld => Table_Wkspace(1:global_table_size,2)
      VPCellElectronTableOld   => Table_Wkspace(1:global_table_size,3)
      OnSiteElectronTable      => Table_Wkspace(1:global_table_size,4)
      MTSphereElectronTable    => Table_Wkspace(1:global_table_size,5)
      VPCellElectronTable      => Table_Wkspace(1:global_table_size,6)
      NetMomentTable      => Table_Wkspace(1:global_table_size,7)
      MTSphereMomentTable => Table_Wkspace(1:global_table_size,8)
      VPCellMomentTable   => Table_Wkspace(1:global_table_size,9)
      ExchangeEnergyTable => Table_Wkspace(1:global_table_size,10)
   else
      allocate( Table_Wkspace(1:global_table_size,1:6) )
      OnSiteElectronTableOld   => Table_Wkspace(1:global_table_size,1)
      MTSphereElectronTableOld => Table_Wkspace(1:global_table_size,2)
      VPCellElectronTableOld   => Table_Wkspace(1:global_table_size,3)
      OnSiteElectronTable      => Table_Wkspace(1:global_table_size,4)
      MTSphereElectronTable    => Table_Wkspace(1:global_table_size,5)
      VPCellElectronTable      => Table_Wkspace(1:global_table_size,6)
      nullify( NetMomentTable )
      nullify( MTSphereMomentTable, VPCellMomentTable, ExchangeEnergyTable )
   endif
!
   allocate(local_array_line(LocalNumAtoms))
   local_array_line = 0
   local_array_size = 0
   do id = 1, LocalNumAtoms
      local_array_line(id) = local_array_size
      local_array_size = local_array_size + getLocalNumSpecies(id)
   enddo
!
   allocate( Vint(LocalNumAtoms), q_mix(LocalNumAtoms) )
   allocate( Qmt(local_array_size), Qmt_old(local_array_size) )
   allocate( Qvp(local_array_size), Qvp_old(local_array_size) )
   allocate( Mmt(local_array_size) )
   allocate( Mvp(local_array_size) )
   allocate( ExEn(local_array_size) )
!
   OnSiteElectronTable = ZERO
   MTSphereElectronTable = ZERO
   VPCellElectronTable = ZERO
   rhoint = ZERO
!
   OnSiteElectronTableOld = ZERO
   MTSphereElectronTableOld = ZERO
   VPCellElectronTableOld = ZERO
   rhoint_old = ZERO
!
   if ( n_spin_pola==2 ) then
     NetMomentTable = ZERO
     MTSphereMomentTable = ZERO
     VPCellMomentTable = ZERO
     ExchangeEnergyTable = ZERO
   endif
   mdenint = ZERO
!
   p_DL => ChargeList
   do id=2,LocalNumAtoms
      allocate(p_DL%next)
      p_DL => p_DL%next
   enddo
   nullify( p_DL%next )
!
   Vint = ZERO
   Qmt = ZERO
   Qvp = ZERO
   Mmt = ZERO
   Mvp = ZERO
   ExEn = ZERO
!
!  ===================================================================
!  setup the initial Qmt_old, Qvp_old
!  -------------------------------------------------------------------
   call setupChargeList('Old')
!  -------------------------------------------------------------------
   p_DL => ChargeList
   jmt = 0
   do id = 1, LocalNumAtoms
      jmt = max(jmt, p_DL%jmt)
      p_DL => p_DL%next
   enddo
!
   allocate( rho0(1:jmt) )
   qint_old = ZERO
   GVint = ZERO
   p_DL => ChargeList
   do id = 1, LocalNumAtoms
!     ----------------------------------------------------------------
      q_mix(id) = getMixingParam4Charge(id)
!     ----------------------------------------------------------------
      jmt = p_DL%jmt
      NumRs = p_DL%NumRs
      rmt = p_DL%rmt
      Vint(id) = getVolume(id) - PI4*rmt*rmt*rmt*THIRD
      GVint = GVint + Vint(id)
      r_mesh => p_DL%r_mesh(1:NumRs)
      ig = getGlobalIndex(id)
      do ia = 1, p_DL%NumSpecies
         lid = local_array_line(id) + ia
!        -------------------------------------------------------------
         call dcopy(jmt,p_DL%EMD(ia)%rho0,1,rho0,1)
!        -------------------------------------------------------------
         corr = rho0(1)*r_mesh(1)*r_mesh(1)*r_mesh(1)*PI2
         Qmt_old(lid) = getVolumeIntegration(id,jmt,r_mesh,0,rho0,truncated=.false.) + corr
         Qvp_old(lid) = p_DL%EMD(ia)%VPCharge
         qint_old = qint_old + getLocalSpeciesContent(id,ia)*(getAtomicNumber(ig,ia)-Qmt_old(lid))
      enddo
      p_DL => p_DL%next
   enddo
   deallocate( rho0 )
!
   sums(1) = GVint
   sums(2) = qint_old
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID, sums, 2)
!  -------------------------------------------------------------------
   GVint = sums(1)
   qint_old = sums(2)
!
   if (isASAPotential()) then
      rhoint_old = ZERO
   else
      rhoint_old = qint_old/GVint
   endif
!
   do id = 1, LocalNumAtoms
      ig = getGlobalIndex(id)
      do ia = 1, getLocalNumSpecies(id)
         lig = global_table_line(ig) + ia
         lid = local_array_line(id) + ia
         OnSiteElectronTableOld(lig) = Qmt_old(lid)+rhoint_old*Vint(id)
         MTSphereElectronTableOld(lig) = Qmt_old(lid)
         VPCellElectronTableOld(lig) = Qvp_old(lid)
      enddo
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,Table_Wkspace,global_table_size,3)
!  -------------------------------------------------------------------
!  ===================================================================
!
   Initialized = .true.
   Updated = .false.
!
   end subroutine initChargeDistribution
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endChargeDistribution()
!  ===================================================================
   implicit none
!
   nullify( OnSiteElectronTableOld )
   nullify( MTSphereElectronTableOld )
   nullify( VPCellElectronTableOld )
   nullify( OnSiteElectronTable )
   nullify( MTSphereElectronTable )
   nullify( VPCellElectronTable )
   if ( n_spin_pola==2 ) then
      nullify( NetMomentTable )
      nullify( MTSphereMomentTable )
      nullify( VPCellMomentTable )
      nullify( ExchangeEnergyTable )
   endif
   deallocate( Table_Wkspace, global_table_line, local_array_line)
   deallocate( Vint, Qmt, Qvp, Mmt, Mvp, ExEn, Qmt_old, Qvp_old, q_mix )
!
   Initialized = .false.
   Updated = .false.
!
   end subroutine endChargeDistribution
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateChargeDistribution(getExchangeEnergy)
!  ===================================================================
   use MPPModule, only : MyPE
   use GroupCommModule, only : GlobalSumInGroup
   use InterpolationModule, only : FitInterp
   use IntegrationModule, only : calIntegration
   use Atom2ProcModule, only : getLocalIndex, getAtom2ProcInGroup, getGlobalIndex
   use PolyhedraModule, only : getVolume
   use StepFunctionModule, only : getVolumeIntegration
   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies
   use AtomModule, only : getLocalSpeciesContent
   use SystemModule, only : getNumAlloyElements, getAlloyElementContent
   use PotentialTypeModule, only : isASAPotential
   implicit none
!
   type (ChargeListStruct), pointer :: p_DL
!
   integer (kind=IntKind) :: id, ig, ip, ia, lid, lig
   integer (kind=IntKind) :: jmt, NumRs
#ifdef No_BLAS
   integer (kind=IntKind) :: j
#endif
!
   real (kind=RealKind) :: rmt, GVint, corr
   real (kind=RealKind) :: qint, mint, qint_old
   real (kind=RealKind) :: sums(4)
   real (kind=RealKind), parameter :: PI8 = PI4*TWO
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), allocatable :: rho0(:)
!
   interface
      function getExchangeEnergy(id,ia) result(exc_en)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: id, ia
         real (kind=RealKind) :: exc_en
      end function getExchangeEnergy
   end interface
!
!  -------------------------------------------------------------------
   call setupChargeList('New')
!  -------------------------------------------------------------------
!
   GVint = ZERO
   jmt = 0
   p_DL => ChargeList
   do id = 1, LocalNumAtoms
      jmt = max(jmt, p_DL%jmt)
      rmt = p_DL%rmt
      Vint(id) = getVolume(id) - PI4*rmt*rmt*rmt*THIRD
      GVint = GVint + Vint(id)
      p_DL => p_DL%next
   enddo
   allocate( rho0(1:jmt) )
   sums(1) = GVint
!
   qint = ZERO
   qint_old = ZERO
   mint = ZERO
   p_DL => ChargeList
   do id = 1, LocalNumAtoms
      jmt = p_DL%jmt
      NumRs = p_DL%NumRs
      rmt = p_DL%rmt
      r_mesh => p_DL%r_mesh(1:NumRs)
      do ia = 1, p_DL%NumSpecies
         lid = local_array_line(id)+ia
!        -------------------------------------------------------------
         call dcopy(jmt,p_DL%EMD(ia)%rho0,1,rho0,1)
!        -------------------------------------------------------------
         corr = rho0(1)*r_mesh(1)*r_mesh(1)*r_mesh(1)*PI2
         Qmt(lid) = getVolumeIntegration(id,jmt,r_mesh,0,rho0,truncated=.false.) + corr
         Qvp(lid) = p_DL%EMD(ia)%VPCharge
         qint = qint + getLocalSpeciesContent(id,ia)*(getLocalAtomicNumber(id,ia)-Qmt(lid))
!        =============================================================
!        Mixing the new Qmt, Qvp with the old ones.
!        =============================================================
         Qmt_old(lid) = q_mix(id)*Qmt(lid) + (ONE-q_mix(id))*Qmt_old(lid)
         Qvp_old(lid) = q_mix(id)*Qvp(lid) + (ONE-q_mix(id))*Qmt_old(lid)
         qint_old = qint_old +                                        &
                  getLocalSpeciesContent(id,ia)*(getLocalAtomicNumber(id,ia)-Qmt_old(lid))
!        =============================================================
         if (n_spin_pola == 2) then
!           ----------------------------------------------------------
            call dcopy(jmt,p_DL%EMD(ia)%mom0,1,rho0,1)
!           ----------------------------------------------------------
            corr = rho0(1)*r_mesh(1)*r_mesh(1)*r_mesh(1)*PI2
            Mmt(lid) = getVolumeIntegration(id,jmt,r_mesh,0,rho0) + corr
            Mvp(lid) = p_DL%EMD(ia)%VPMoment
            mint = mint + getLocalSpeciesContent(id,ia)*(Mvp(lid) - Mmt(lid))
            ExEn(lid) = getExchangeEnergy(id,ia) ! This is the exchange field that
                                                 ! causes spin-up and spin-down DOS peak to split
         else
            Mmt(lid) = ZERO
            Mvp(lid) = ZERO
            ExEn(lid) = ZERO
         endif
      enddo
      p_DL => p_DL%next
   enddo
!
   deallocate( rho0 )
!
   sums(2) = qint
   sums(3) = qint_old
   if (n_spin_pola == 2) then
      sums(4) = mint
   else
      sums(4) = ZERO
   endif
!
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID, sums, n_spin_pola+2)
!  -------------------------------------------------------------------
   GVint = sums(1)
   qint = sums(2)
   qint_old = sums(3)
   mint = sums(4)
!
   if (isASAPotential()) then
      rhoint = ZERO
      rhoint_old = ZERO
      mdenint = ZERO
   else
      rhoint = qint/GVint
      rhoint_old = qint_old/GVint
      mdenint = mint/GVint
   endif
!
   if (MyPE == 0) then
      write(6,'(/,a,f20.14,/)')'updateChargeDistribution :: Interstitial rho0  = ', &
            rhoint
   endif
!
   Table_Wkspace = ZERO
   do id = 1, LocalNumAtoms
      ig = getGlobalIndex(id)
      do ia = 1, getLocalNumSpecies(id)
         lig = global_table_line(ig) + ia
         lid = local_array_line(id) + ia
         Table_Wkspace(lig,1) = Qmt(lid)+rhoint*Vint(id)
         Table_Wkspace(lig,2) = Qmt(lid)
         Table_Wkspace(lig,3) = Qvp(lid)
         Table_Wkspace(lig,4) = Qmt_old(lid)+rhoint_old*Vint(id)
         Table_Wkspace(lig,5) = Qmt_old(lid)
         Table_Wkspace(lig,6) = Qvp_old(lid)
         if ( n_spin_pola==2 ) then
            Table_Wkspace(lig,7) = Mmt(lid) + mdenint*Vint(id)
            Table_Wkspace(lig,8) = Mmt(lid)
            Table_Wkspace(lig,9) = Mvp(lid)
            Table_Wkspace(lig,10) = ExEn(lid)
         endif
      enddo
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,Table_Wkspace,global_table_size,6+4*(n_spin_pola-1))
!  -------------------------------------------------------------------
!
   mom_aver = ZERO
   if ( n_spin_pola==2 ) then
      do ig = 1, GlobalNumAtoms
         do ia = 1, getNumAlloyElements(ig)
            lig = global_table_line(ig) + ia
            mom_aver = mom_aver + getAlloyElementContent(ig,ia)*VPCellMomentTable(lig)
         enddo
      enddo
   endif
   mom_aver = mom_aver/GlobalNumAtoms
!
   Updated = .true.
!
   end subroutine updateChargeDistribution
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalTableLine() result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: p(:)
!
   p => global_table_line(1:GlobalNumAtoms)
!
   end function getGlobalTableLine
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalOnSiteElectronTable() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   r => OnSiteElectronTable(1:global_table_size)
!
   end function getGlobalOnSiteElectronTable
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalMTSphereElectronTable() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   r => MTSphereElectronTable(1:global_table_size)
!
   end function getGlobalMTSphereElectronTable
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalVPCellElectronTable() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   r => VPCellElectronTable(1:global_table_size)
!
   end function getGlobalVPCellElectronTable
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalOnSiteElectronTableOld() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   r => OnSiteElectronTableOld(1:global_table_size)
!
   end function getGlobalOnSiteElectronTableOld
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalMTSphereElectronTableOld() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   r => MTSphereElectronTableOld(1:global_table_size)
!
   end function getGlobalMTSphereElectronTableOld
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalVPCellElectronTableOld() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   r => VPCellElectronTableOld(1:global_table_size)
!
   end function getGlobalVPCellElectronTableOld
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalNetMomentTable() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   if (n_spin_pola == 2) then
      r => NetMomentTable(1:global_table_size)
   else
      nullify(r)
   endif
!
   end function getGlobalNetMomentTable
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalMTSphereMomentTable() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   if (n_spin_pola == 2) then
      r => MTSphereMomentTable(1:global_table_size)
   else
      nullify(r)
   endif
!
   end function getGlobalMTSphereMomentTable
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalVPCellMomentTable() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   if (n_spin_pola == 2) then
      r => VPCellMomentTable(1:global_table_size)
   else
      nullify(r)
   endif
!
   end function getGlobalVPCellMomentTable
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function  getGlobalExchangeEnergyTable() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: r(:)
!
   r => ExchangeEnergyTable(1:global_table_size)
!
   end function  getGlobalExchangeEnergyTable
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterstitialElectronDensity() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: r
!
   r = rhoint
!
   end function getInterstitialElectronDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterstitialElectronDensityOld() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: r
!
   r = rhoint_old
!
   end function getInterstitialElectronDensityOld
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterstitialMomentDensity() result(r)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: r
!
   r = mdenint
!
   end function getInterstitialMomentDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAverageMoment() result(m)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: m
!
   m = mom_aver
!
   end function getAverageMoment
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printChargeDistribution(iter,fu)
!  ===================================================================
   use Atom2ProcModule, only : getAtom2ProcInGroup
   use SystemModule, only : getAtomicNumber, getAtomName, getNumAlloyElements
   implicit none
!
   logical :: FileExist
!
   integer (kind=IntKind), intent(in) :: iter
   integer (kind=IntKind), intent(in) :: fu
   integer (kind=IntKind) :: ig, ia, lig
!
   if (.not.Initialized) then
      call ErrorHandler('printChargeDistributionTable',               &
                        'Module is not initialized')
   else if (.not.Updated) then
      call ErrorHandler('printChargeDistributionTable','Updated Table is empty')
   endif
!
   if (fu /= 6) then
      inquire(file='ChargeTable',exist=FileExist)
      if (FileExist) then
         open(unit=fu,file='ChargeTable',form='formatted',status='old', &
              position='append')
      else
         open(unit=fu,file='ChargeTable',form='formatted',status='unknown')
         write(fu,'(/,80(''-''))')
         write(fu,'(/,20x,a)')'***************************************'
         write(fu,'( 20x,a )')'* Output from printChargeDistribution *'
         write(fu,'(20x,a,/)')'***************************************'
      endif
   else
      write(fu,'(/,80(''-''))')
      write(fu,'(/,20x,a)')'***************************************'
      write(fu,'( 20x,a )')'* Output from printChargeDistribution *'
      write(fu,'(20x,a,/)')'***************************************'
      write(fu,'(2x,''Interstitial Electron Density ='',f17.8)') rhoint
      write(fu,'(2x,''Interstitial Moment   Density ='',f17.8)') mdenint
   endif
!
   write(fu,'(/,a,i5)')'# ITERATION :',iter
   write(fu,'(80(''=''))')
   if (n_spin_pola == 1) then
      write(fu,'(a)')' Atom   Index       Q         Qmt        Qvp          dQ'
      write(fu,'(80(''=''))')
      do ig = 1, GlobalNumAtoms
         do ia = 1, getNumAlloyElements(ig)
            lig = global_table_line(ig) + ia
            write(fu,'(2x,a3,2x,i7,4(2x,f9.5))')getAtomName(ig,ia),ig, &
                  OnSiteElectronTable(lig),MTSphereElectronTable(lig), &
                  VPCellElectronTable(lig),                            &
                  OnSiteElectronTable(lig)-getAtomicNumber(ig,ia)
         enddo
      enddo
   else
      write(fu,'(a)')      &
' Atom   Index       Q         Qmt        Qvp          dQ        Mmt        Mvp'
      write(fu,'(80(''=''))')
      do ig = 1, GlobalNumAtoms
         do ia = 1, getNumAlloyElements(ig)
            lig = global_table_line(ig) + ia
            write(fu,'(2x,a3,2x,i7,6(2x,f9.5))')getAtomName(ig,ia),ig, &
                  OnSiteElectronTable(lig),MTSphereElectronTable(lig), &
                  VPCellElectronTable(lig),                            &
                  OnSiteElectronTable(lig)-getAtomicNumber(ig,ia),     &
                  MTSphereMomentTable(lig),VPCellMomentTable(lig)
!                 NetMomentTable(lig)
         enddo
      enddo
   endif
   write(fu,'(80(''=''))')
!
   if (fu /= 6) then
      close(unit=fu)
   endif
!
   end subroutine printChargeDistribution
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupChargeList(stage)
!  ===================================================================
   use KindParamModule, only : IntKind
!
   use MathParamModule, only : ZERO
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid
!
   use DataServiceCenterModule, only : getDataStorage, RealType, RealMark
!
   use CoreStatesModule, only : getCoreVPCharge, getCoreVPMoment
!
   use AtomModule, only : getLocalNumSpecies
!
   implicit none
!
   character (len=3), intent(in) :: stage
   character (len=27) :: rho_storage
   character (len=25) :: mom_storage
!
   integer (kind=IntKind) :: id, nr, ia
!
   real (kind=RealKind), pointer :: pr(:), rho0(:,:), mom0(:,:)
!
   type (GridStruct), pointer :: Grid
   type (ChargeListStruct), pointer :: p_DL
!
   rho_storage = stage//'SphericalElectronDensity'
   mom_storage = stage//'SphericalMomentDensity'
!
   p_DL => ChargeList
   do id = 1, LocalNumAtoms
      nr = getNumRmesh(id)
      p_DL%NumSpecies = getLocalNumSpecies(id)
      p_DL%r_mesh => getRmesh(id)
      Grid=>getGrid(id)
      p_DL%jmt = Grid%jmt
      p_DL%NumRs = nr
      p_DL%rmt = Grid%rmt
      allocate(p_DL%EMD(p_DL%NumSpecies))
      pr => getDataStorage(id,rho_storage,(nr+1)*getLocalNumSpecies(id),RealMark)
      rho0 => aliasArray2_r(pr,nr+1,getLocalNumSpecies(id))
      if (n_spin_pola == 2) then
         pr => getDataStorage(id,mom_storage,(nr+1)*getLocalNumSpecies(id),RealMark)
         mom0 => aliasArray2_r(pr,nr+1,getLocalNumSpecies(id))
      endif
      do ia = 1, p_DL%NumSpecies
         p_DL%EMD(ia)%rho0 => rho0(:,ia)
         p_DL%EMD(ia)%VPCharge = p_DL%EMD(ia)%rho0(nr+1) !+ getCoreVPCharge(id)
         if (n_spin_pola == 2) then
            p_DL%EMD(ia)%mom0 => mom0(:,ia)
            p_DL%EMD(ia)%VPMoment = p_DL%EMD(ia)%mom0(nr+1) !+ getCoreVPMoment(id)
         else
            p_DL%EMD(ia)%VPMoment = ZERO
            nullify( p_DL%EMD(ia)%mom0 )
         endif
      enddo
      p_DL => p_DL%next
   enddo
!
   end subroutine setupChargeList
!  ===================================================================
end module ChargeDistributionModule
