module CPAMediumModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initCPAMedium,            &
          endCPAMedium,             &
          computeCPAMedium,         &
          getNumCPAMediums,         &
          getNumSpeciesInCPAMedium, &
          getCPAMatrix,             &
          getImpurityMatrix
!
private
   integer (kind=IntKind) :: MaxIterations = 50
   integer (kind=IntKind) :: iteration = 0
!
   integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
   integer (kind=IntKind) :: nSpinCant, NumCPAMediums = 0
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind), allocatable :: MediumIndex(:)
   integer (kind=IntKind) :: kmax_kkr_max
   integer (kind=IntKind) :: ndim_Tmat
   integer (kind=IntKind) :: print_instruction
!
   type CPAMatrixStruct
      real (kind=RealKind) :: content
      complex (kind=CmplxKind), pointer :: tmat_a(:,:) ! In global spin framework
      complex (kind=CmplxKind), pointer :: tau_a(:,:,:)  ! In local spin framework
      complex (kind=CmplxKind), pointer :: kau_a(:,:,:) ! In local spin framework
   end type CPAMatrixStruct
!
   type CPAMediumStruct ! In global spin framework
      integer (kind=IntKind) :: dsize
      integer (kind=IntKind) :: local_index
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: num_species
      complex (kind=CmplxKind), pointer :: tau_c(:,:,:) ! In local spin framework
      complex (kind=CmplxKind), pointer :: Tcpa(:,:)
      complex (kind=CmplxKind), pointer :: TcpaInv(:,:)
      complex (kind=CmplxKind), pointer :: Tcpa_old(:,:)
      complex (kind=CmplxKind), pointer :: TcpaInv_old(:,:)
      type (CPAMatrixStruct), pointer :: CPAMatrix(:)
   end type CPAMediumStruct
!
   type (CPAMediumStruct), allocatable :: CPAMedium(:)
!
   logical :: isRelativistic = .false.
!
   character (len=50) :: stop_routine
!
   complex (kind=CmplxKind), allocatable, target :: Tcpa(:)
   complex (kind=CmplxKind), allocatable, target :: TcpaInv(:)
   complex (kind=CmplxKind), allocatable, target :: Tcpa_old(:)
   complex (kind=CmplxKind), allocatable, target :: TcpaInv_old(:)
!  complex (kind=CmplxKind), allocatable, target :: Tau(:)
   complex (kind=CmplxKind), allocatable, target :: TauA(:)
   complex (kind=CmplxKind), allocatable, target :: KauA(:)
   complex (kind=CmplxKind), allocatable, target :: Tmat_global(:,:)
   complex (kind=CmplxKind), allocatable, target :: WORK0(:), WORK1(:), WORK2(:)
!
   real (kind=RealKind) :: CPA_tolerance = TEN2m8
   real (kind=RealKind) :: CPA_alpha = 0.10d0
   real (kind=RealKind) :: CPA_slow_alpha = 0.02d0
   real (kind=RealKind) :: CPA_switch_param = 0.003d0
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initCPAMedium(cant, lmax_kkr, rel, cpa_mix_type, cpa_max_iter,  &
                            cpa_mix_0, cpa_mix_1, cpa_eswitch, cpa_tol, istop, iprint)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : getNumGroups, getGroupLabel, isGroupExisting
!
   use CrystalMatrixModule, only : initCrystalMatrix
   use MediumHostModule, only  : getNumSites, getLocalNumSites, getNumSpecies, &
                                 getLmaxKKR, getGlobalSiteIndex, getSpeciesContent
!
   use AccelerateCPAModule, only : initAccelerateCPA
!
   use EmbeddedClusterModule, only : initEmbeddedCluster
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: cant, rel, cpa_mix_type, cpa_max_iter
   integer (kind=IntKind), intent(in) :: iprint(:)
   integer (kind=IntKind), intent(in) :: lmax_kkr(:)
   integer (kind=IntKind) :: i, ig, kmaxi, ic, n, NumImpurities
   integer (kind=IntKind) :: aid, num, dsize, lmaxi, nsize
!
   real (kind=RealKind), intent(in) :: cpa_mix_0, cpa_mix_1, cpa_eswitch, cpa_tol
!
   character (len=14) :: sname = "initCPAMedium"
!
   complex (kind=CmplxKind), pointer :: p1(:)
!
   GlobalNumSites = getNumSites()
   LocalNumSites = getLocalNumSites()
   nSpinCant = cant
!
!  ===================================================================
!  This needs to be further thought through for whether it needs to 
!  create a Medium Cell communicator.
!  ===================================================================
   if (isGroupExisting('Medium Cell')) then
      GroupID = getGroupID('Medium Cell')
      NumPEsInGroup = getNumPEsInGroup(GroupID)
      MyPEinGroup = getMyPEinGroup(GroupID)
   else if (isGroupExisting('Unit Cell')) then
      GroupID = getGroupID('Unit Cell')
      NumPEsInGroup = getNumPEsInGroup(GroupID)
      MyPEinGroup = getMyPEinGroup(GroupID)
   else
      GroupID = -1
      NumPEsInGroup = 1
      MyPEinGroup = 0
   endif
!
   if (rel>1) then
      isRelativistic = .true.
   else
      isRelativistic = .false.
   endif
!
   MaxIterations = cpa_max_iter
   CPA_tolerance = cpa_tol
   CPA_alpha = cpa_mix_0
   CPA_slow_alpha = cpa_mix_1
   CPA_switch_param = cpa_eswitch
!
   kmax_kkr_max = 0
   do ig = 1, GlobalNumSites
      lmaxi = getLmaxKKR(ig)
      kmaxi = (lmaxi+1)**2
      kmax_kkr_max = max( kmax_kkr_max,kmaxi )
   enddo
!
   allocate(WORK0(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant))
   allocate(WORK1(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant))
   allocate(WORK2(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant))
   allocate(MediumIndex(LocalNumSites))
!
   NumCPAMediums = 0
   NumImpurities = 0
   do i = 1, LocalNumSites
      ig = getGlobalSiteIndex(i)
      if (getNumSpecies(ig) > 1) then
         NumCPAMediums = NumCPAMediums + 1
         NumImpurities = NumImpurities + getNumSpecies(ig)
      endif
   enddo
!
   if (NumCPAMediums > 0) then
      allocate(CPAMedium(NumCPAMediums))
      allocate(Tcpa(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(TcpaInv(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(Tcpa_old(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(TcpaInv_old(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
!     allocate(Tau(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(TauA(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumImpurities))
      allocate(KauA(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumImpurities))
      if (nSpinCant == 2) then
         allocate(Tmat_global(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant,NumImpurities))
      endif
      Tcpa = CZERO; TcpaInv = CZERO
!     Tau = CZERO
      TauA=CZERO
      KauA = CZERO
   endif
!
   n = 0
   NumImpurities = 0
   MediumIndex = 0
   ndim_Tmat = 0; aid = 0
   do i = 1, LocalNumSites
      ig = getGlobalSiteIndex(i)
      num = getNumSpecies(ig)
      if (num > 1) then
         n = n + 1
         MediumIndex(i) = n
         CPAMedium(n)%num_species = num
         CPAMedium(n)%local_index = i
         CPAMedium(n)%global_index = ig
         lmaxi = getLmaxKKR(ig)
         dsize = (lmaxi+1)**2
         CPAMedium(n)%dsize = dsize
         nsize = dsize*nSpinCant
!        CPAMedium(n)%Tcpa => aliasArray2_c(Tcpa(ndim_Tmat+1:),nsize,nsize)
         p1 => Tcpa(ndim_Tmat+1:)
         CPAMedium(n)%Tcpa => aliasArray2_c(p1,nsize,nsize)
!        CPAMedium(n)%TcpaInv => aliasArray2_c(TcpaInv(ndim_Tmat+1:),nsize,nsize)
         p1 => TcpaInv(ndim_Tmat+1:)
         CPAMedium(n)%TcpaInv => aliasArray2_c(p1,nsize,nsize)
!        CPAMedium(n)%Tcpa_old => aliasArray2_c(Tcpa_old(ndim_Tmat+1:),nsize,nsize)
         p1 => Tcpa_old(ndim_Tmat+1:)
         CPAMedium(n)%Tcpa_old => aliasArray2_c(p1,nsize,nsize)
!        CPAMedium(n)%TcpaInv_old => aliasArray2_c(TcpaInv_old(ndim_Tmat+1:),nsize,nsize)
         p1 => TcpaInv_old(ndim_Tmat+1:)
         CPAMedium(n)%TcpaInv_old => aliasArray2_c(p1,nsize,nsize)
!        CPAMedium(n)%tau_c => aliasArray3_c(Tau(ndim_Tmat+1:),dsize,dsize,nSpinCant*nSpinCant)
         nullify(CPAMedium(n)%tau_c)
         ndim_Tmat = ndim_Tmat + nsize*nsize
         allocate(CPAMedium(n)%CPAMatrix(num))
         do ic = 1, num
            NumImpurities = NumImpurities + 1
            CPAMedium(n)%CPAMatrix(ic)%content = getSpeciesContent(ic,ig)
            if (nSpinCant == 2) then ! Spin-canted case, use the locally allocated space
!              CPAMedium(n)%CPAMatrix(ic)%tmat_a => aliasArray2_c(Tmat_global(:,NumImpurities),nsize,nsize)
               p1 => Tmat_global(:,NumImpurities)
               CPAMedium(n)%CPAMatrix(ic)%tmat_a => aliasArray2_c(p1,nsize,nsize)
            else ! In (non-)spin-polarized case, use the Tmat space in SSSolverModule
               nullify(CPAMedium(n)%CPAMatrix(ic)%tmat_a)
            endif
            p1 => TauA(aid+1:)
            CPAMedium(n)%CPAMatrix(ic)%tau_a => aliasArray3_c(p1,dsize,dsize,nSpinCant*nSpinCant)
!           nullify(CPAMedium(n)%CPAMatrix(ic)%tau_a)
            p1 => KauA(aid+1:)
            CPAMedium(n)%CPAMatrix(ic)%kau_a => aliasArray3_c(p1,dsize,dsize,nSpinCant*nSpinCant)
            aid = aid + nsize*nsize
         enddo
      endif
   enddo
!
   iteration = 0
   print_instruction = maxval(iprint)
!
!  -------------------------------------------------------------------
   call initCrystalMatrix(LocalNumSites, cant, lmax_kkr, rel, istop, iprint)
!  -------------------------------------------------------------------
   call initEmbeddedCluster(cant, istop, print_instruction)
!  -------------------------------------------------------------------
   call initAccelerateCPA(cpa_mix_type, cpa_max_iter, cpa_mix_0, cpa_tol,  &
                          kmax_kkr_max)
!  -------------------------------------------------------------------
!
   end subroutine initCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCPAMedium()
!  ===================================================================
   use CrystalMatrixModule, only : endCrystalMatrix
   use AccelerateCPAModule, only : endAccelerateCPA
   use EmbeddedClusterModule, only : endEmbeddedCluster
!
   implicit none
!
   integer (kind=IntKind) :: n, ic
!
   deallocate(WORK0, WORK1, WORK2, MediumIndex)
!
   do n = 1, NumCPAMediums
      nullify(CPAMedium(n)%tau_c, CPAMedium(n)%Tcpa, CPAMedium(n)%TcpaInv)
      nullify(CPAMedium(n)%Tcpa_old, CPAMedium(n)%TcpaInv_old)
      do ic = 1, CPAMedium(n)%num_species
         nullify(CPAMedium(n)%CPAMatrix(ic)%tau_a, CPAMedium(n)%CPAMatrix(ic)%tmat_a)
         nullify(CPAMedium(n)%CPAMatrix(ic)%kau_a)
      enddo
      deallocate(CPAMedium(n)%CPAMatrix)
   enddo
!
   if (NumCPAMediums > 0) then
      deallocate(Tcpa, TcpaInv, Tcpa_old, TcpaInv_old, CPAMedium)
!     deallocate(Tau)
      deallocate(KauA, TauA)
      if (nSpinCant ==2) then
         deallocate(Tmat_global)
      endif
   endif
!
!  -------------------------------------------------------------------
   call endAccelerateCPA()
!  -------------------------------------------------------------------
   call endEmbeddedCluster()
!  -------------------------------------------------------------------
   call endCrystalMatrix()
!  -------------------------------------------------------------------
!
   NumCPAMediums = 0
   iteration = 0
!
   end subroutine endCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCPAMediums() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumCPAMediums
!
   end function getNumCPAMediums
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeCPAMedium(e)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
!
   use SSSolverModule, only : getScatteringMatrix
!
!  use CrystalMatrixModule, only : calCrystalMatrix
!  use CrystalMatrixModule, only : getCPAMediumTau => getTau
!
   use AccelerateCPAModule, only : initializeAcceleration, accelerateCPA
   use AccelerateCPAModule, only : setAccelerationParam
!
   use EmbeddedClusterModule, only : setupHostMedium, getTau
!
   use WriteMatrixModule,  only : writeMatrix
!
   use MatrixModule,  only : computeUAU
!
   implicit none
!
   character (len=12) :: description
!
   integer (kind=IntKind) :: ia, id, nsize, n, dsize
   integer (kind=IntKind) :: ns, is, js
   integer (kind=IntKind) :: site_config(LocalNumSites)
!
!  ===================================================================
!  CPA iteration acceleration parameters...
!  These parameters are taken from the mkkrcpa code
!  ===================================================================
   integer (kind=IntKind), parameter :: ipits = 4
   real (kind=RealKind), parameter ::  ctol=1.0d-08, cmix=0.15d0, cw0=5.0d-03
!  ===================================================================
!
   real (kind=RealKind) :: err, max_err
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind), pointer :: tau_a(:,:,:)
   complex (kind=CmplxKind), pointer :: mat_a(:,:)
   complex (kind=CmplxKind), pointer :: kau_a(:,:)
   complex (kind=CmplxKind), pointer :: Jost(:,:), OH(:,:), Tinv(:,:), Jinv(:,:)
   complex (kind=CmplxKind), pointer :: SinvL(:,:), SinvR(:,:)
!  
   if (nSpinCant == 1) then
      do n = 1, NumCPAMediums
         id = CPAMedium(n)%local_index
         do ia = 1, CPAMedium(n)%num_species
            CPAMedium(n)%CPAMatrix(ia)%tmat_a =>                      &
               getScatteringMatrix('T-Matrix',spin=1,site=id,atom=ia)
         enddo
      enddo
   else
      call ErrorHandler('computeCPAMedium',                           &
               'Spin-canting calculation needs to be checked in CPA case')
   endif
!
   kappa = sqrt(e)
!
   do n = 1, NumCPAMediums
!     ----------------------------------------------------------------
      call averageTMatrix(n)
!     ----------------------------------------------------------------
   enddo
!
!  ===================================================================
!  Note that:
!       CPAMedium(n)%Tcpa(:,:,:) is an alias of Tcpa(:)
!       CPAMedium(n)%Tcpa_old(:,:,:) is an alias of Tcpa_old(:)
!       CPAMedium(n)%TcpaInv_old(:,:,:) is an alias of TcpaInv_old(:)
!  ===================================================================
   if (aimag(e) >= CPA_switch_param .or. real(e) < ZERO) then
!     ----------------------------------------------------------------
      call setAccelerationParam(CPA_alpha)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call setAccelerationParam(CPA_slow_alpha)
!     ----------------------------------------------------------------
   endif
!
   do id = 1, LocalNumSites
      if (MediumIndex(id) > 0) then ! This is a random alloy sublattice site
         site_config(id) = 0 ! Set the site to be the CPA medium site
      endif
   enddo
!
   iteration = 0
   LOOP_iter: do while (iteration < MaxIterations)
      iteration = iteration + 1
!     write(6,'(a,i5)')'At iteration: ',iteration
      max_err = ZERO
!     ----------------------------------------------------------------
      call setupHostMedium(e,getSingleSiteTmat,configuration=site_config)   
!     ----------------------------------------------------------------
      do n = 1, NumCPAMediums
         dsize = CPAMedium(n)%dsize
         nsize = dsize*nSpinCant
         CPAMedium(n)%Tcpa_old = CPAMedium(n)%Tcpa
         CPAMedium(n)%TcpaInv_old = CPAMedium(n)%TcpaInv
!        -------------------------------------------------------------
!        call writeMatrix('tin',CPAMedium(n)%Tcpa,nsize,nsize,TEN2m8)
!        -------------------------------------------------------------
         call initializeAcceleration(CPAMedium(n)%Tcpa,nsize*nsize,iteration)
!        -------------------------------------------------------------
!
!        =============================================================
!        calculate Tau00, Tauij, of the medium made of tmat_c's
!        -------------------------------------------------------------
!        call calCrystalMatrix(e,getSingleSiteTmat,use_tmat=.true.,   &
!                              tau_needed=.true.,configuration=site_config)
!        call calCrystalMatrix(e,getSingleSiteMatrix,                 &
!                              tau_needed=.true.,configuration=site_config)
!        -------------------------------------------------------------
!
!        call writeMatrix('t-cpa',CPAMedium(n)%Tcpa,nsize,nsize,TEN2m8)
!        =============================================================
!        Assume that the CPA mediums on each sub-lattice are not correlated,
!        i.e., we are taking the single site approximation.
!        -------------------------------------------------------------
         call iterateCPAMedium(n)
!        -------------------------------------------------------------
         call checkCPAMedium(n,err)
!        -------------------------------------------------------------
         if (print_instruction >= 0) then
            write(6,'(a,2i4,2x,d15.8)')'In computeCPAMedium: iter, medium, err = ', &
                  iteration, n, err   
         endif
!
!        -------------------------------------------------------------
!        call writeMatrix('tcout',CPAMedium(n)%Tcpa,nsize,nsize,TEN2m8)
!        -------------------------------------------------------------
         call accelerateCPA(CPAMedium(n)%Tcpa,nsize*nsize,iteration)
!        -------------------------------------------------------------
         CPAMedium(n)%TcpaInv = CPAMedium(n)%Tcpa
!        -------------------------------------------------------------
         call MtxInv_LU(CPAMedium(n)%TcpaInv,nsize)
!        -------------------------------------------------------------
         max_err = max(max_err,err)
      enddo
!
      if (max_err < CPA_tolerance) then
         exit LOOP_iter
      endif
   enddo LOOP_iter
!
   if (print_instruction >= 1) then
      do n = 1, NumCPAMediums
         dsize = CPAMedium(n)%dsize
         nsize = dsize*nSpinCant
         do ia = 1, CPAMedium(n)%num_species
            write(description,'(a,i2)')'t-matrix:',ia
            call writeMatrix(description,                                &
                             CPAMedium(n)%CPAMatrix(ia)%tmat_a,nsize,nsize,TEN2m8)
         enddo
         call writeMatrix('Final t-cpa',CPAMedium(n)%Tcpa,nsize,nsize,TEN2m8)
         call writeMatrix('Final t-cpa-inv',CPAMedium(n)%TcpaInv,nsize,nsize,TEN2m8)
      enddo
   endif
!  
!  ===================================================================
!  Calculate Tau_a and Kau_a for each species in local spin framework
!  ===================================================================
   do n = 1, NumCPAMediums
      id = CPAMedium(n)%local_index
      dsize = CPAMedium(n)%dsize
      nsize = dsize*nSpinCant
      do ia = 1, CPAMedium(n)%num_species
         tau_a => getTau(local_id=1) ! Associated tau_a space with the space
                                     ! allocated in EmbeddedCluster module
!        =============================================================
!        Substitute one CPA medium site by a real atom. The returning
!        xmat_a is the tau_a matrix.
!        Note: This needs to be carefully checked in the spin-canted case
!              for which tau_a needs to be transformed from the global
!              spin framework to the local spin framework in subroutine
!              substituteTcByTa
!        -------------------------------------------------------------
         call substituteTcByTa(id,ia,spin=1,mat_a=tau_a(:,:,1))
!        -------------------------------------------------------------
         CPAMedium(n)%CPAMatrix(ia)%tau_a = tau_a
         ns = 0
         do js = 1, nSpinCant
            Jinv => getScatteringMatrix('JostInv-Matrix',spin=js,site=id,atom=ia)
            Tinv => getScatteringMatrix('TInv-Matrix',spin=js,site=id,atom=ia)
            SinvL => aliasArray2_c(WORK0,dsize,dsize)
!           ==========================================================
!           S^{-1} = Jost^{-1}*tmat_a^{-1}/kappa
!           ----------------------------------------------------------
            call zgemm( 'n', 'n', dsize, dsize, dsize, CONE/kappa,    &
                        Jinv, dsize, Tinv, dsize, CZERO, SinvL, dsize)
!           ----------------------------------------------------------
            do is = 1, nSpinCant
               Jost => getScatteringMatrix('Jost-Matrix',spin=is,site=id,atom=ia)
               OH => getScatteringMatrix('OmegaHat-Matrix',spin=is,site=id,atom=ia)
               SinvR => aliasArray2_c(WORK1,dsize,dsize)
!              =======================================================
!              OmegaHat = S^{-1} * tmat_a * S^{-T*}/kappa
!              S^{-T*} = Jost*OmegaHat
!              -------------------------------------------------------
               call zgemm( 'n', 'n', dsize, dsize, dsize, CONE,       &
                           Jost, dsize, OH, dsize, CZERO, SinvR, dsize)
!              -------------------------------------------------------
!
               ns = ns + 1
               if (print_instruction >= 1) then
                  call writeMatrix('Tau_a',CPAMedium(n)%CPAMatrix(ia)%tau_a(:,:,ns), &
                                   dsize,dsize,TEN2m8)
               endif
               kau_a => CPAMedium(n)%CPAMatrix(ia)%kau_a(:,:,ns)
               mat_a => CPAMedium(n)%CPAMatrix(ia)%tau_a(:,:,ns)
!              =======================================================
!              kau_a = energy * S^{-1} * tau_a * S^{-T*}
!              -------------------------------------------------------
               call computeUAU(SinvL,dsize,dsize,SinvR,dsize,e,       &
                               mat_a,dsize,CZERO,kau_a,dsize,WORK2)
!              -------------------------------------------------------
               if (is == js) then
                  kau_a = kau_a - kappa*OH
               endif
            enddo
         enddo
      enddo
   enddo
!
   end subroutine computeCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageTMatrix(n)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
!
   use SpinRotationModule, only : rotateLtoG
!
   use SSSolverModule, only : getScatteringMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n ! n = CPA medium index
   integer (kind=IntKind) :: i, ic, ia, nsize
!
   complex (kind=CmplxKind), pointer :: tm1(:,:), tm2(:,:), tm0(:,:)
   complex (kind=CmplxKind), pointer :: tavr(:,:)
   complex (kind=CmplxKind) :: cfac
!
!  *******************************************************************
!  The averaging is performed in the global spin reference framework
!  *******************************************************************
!
   tavr => CPAMedium(n)%Tcpa
   tavr = CZERO
   ia = CPAMedium(n)%local_index
   nsize = CPAMedium(n)%dsize*nSpinCant
   do ic = 1, CPAMedium(n)%num_species
      tm0 => CPAMedium(n)%CPAMatrix(ic)%tmat_a
      if (nSpinCant == 2) then
         tm1 => getScatteringMatrix('T-Matrix',spin=1,site=ia,atom=ic)
         tm2 => getScatteringMatrix('T-Matrix',spin=2,site=ia,atom=ic)
         nsize = size(tm1,1)
!        -------------------------------------------------------------
         call rotateLtoG(i, nsize, nsize, tm1, tm2, tm0)
!        -------------------------------------------------------------
      endif
      cfac = CPAMedium(n)%CPAMatrix(ic)%content
!     ----------------------------------------------------------------
      call zaxpy(nsize*nsize,cfac,tm0,1,tavr,1)
!     ----------------------------------------------------------------
   enddo
   CPAMedium(n)%TcpaInv = CPAMedium(n)%Tcpa
!  -------------------------------------------------------------------
   call MtxInv_LU(CPAMedium(n)%TcpaInv,nsize)
!  -------------------------------------------------------------------
!
   end subroutine averageTMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumSpeciesInCPAMedium(site) result(num)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
   integer (kind=IntKind) :: num, n
!
   if (site < 1 .or. site > LocalNumSites) then
      call ErrorHandler('getNumSpeciesInCPAMedium',                   &
                        'The local atom index is out of range',site)
   endif
!
   n = MediumIndex(site)
   if (n == 0) then
      num = 1
   else
      num = CPAMedium(n)%num_species
   endif
!
   end function getNumSpeciesInCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCPAMatrix(matrix_type,site,atom,dsize,err) result(mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
   integer (kind=IntKind), intent(in), optional :: atom
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind), intent(out), optional :: err
   integer (kind=IntKind) :: n, ia
!
   character (len=*), intent(in) :: matrix_type
!
   complex (kind=CmplxKind), pointer :: mat(:,:)
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (site < 1 .or. site > LocalNumSites) then
      call ErrorHandler('getCPAMatrix','The local atom index is out of range',site)
   endif
!
   n = MediumIndex(site)
   if (n == 0) then
      nullify(mat)
      if (present(dsize)) then
         dsize = 0
      endif
      if (present(err)) then
         err = 1
      endif
      return
   endif
!
   if (present(atom)) then
      if (atom < 0 .or. atom > CPAMedium(n)%num_species) then
         call ErrorHandler('getCPAMatrix','The species index is out of range',atom)
      else
         ia = atom
      endif
   else
      ia = 0
   endif
!
   if (present(dsize)) then
      dsize = CPAMedium(n)%dsize
   endif
   if (present(err)) then
      err = 0
   endif
   if (nocaseCompare(matrix_type,'Tcpa')) then
      mat => CPAMedium(n)%Tcpa
   else if (nocaseCompare(matrix_type,'Tau')) then
      if (ia == 0) then
         mat => CPAMedium(n)%tau_c(:,:,1)
      else
         mat => CPAMedium(n)%CPAMatrix(ia)%tau_a(:,:,1)
      endif
   else if (nocaseCompare(matrix_type,'Tau_cpa')) then
      mat => CPAMedium(n)%tau_c(:,:,1)
   else 
      nullify(mat)
      call ErrorHandler('getCPAMatrix','The matrix type is invalid',matrix_type)
   endif
!
   end function getCPAMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getImpurityMatrix(matrix_type,site,atom,dsize,err) result(mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
   integer (kind=IntKind), intent(in) :: atom
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind), intent(out), optional :: err
   integer (kind=IntKind) :: n
!
   character (len=*), intent(in) :: matrix_type
!
   complex (kind=CmplxKind), pointer :: mat(:,:,:)
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (site < 1 .or. site > LocalNumSites) then
      call ErrorHandler('getImpurityMatrix','The local atom index is out of range',site)
   endif
!
   n = MediumIndex(site)
!
   if (atom < 0 .or. atom > CPAMedium(n)%num_species) then
      call ErrorHandler('getImpurityMatrix','The species index is out of range',atom)
   endif
!
   if (n == 0) then
      nullify(mat)
      if (present(dsize)) then
         dsize = 0
      endif
      if (present(err)) then
         err = 1
      endif
      return
   endif
!
   if (present(dsize)) then
      dsize = CPAMedium(n)%dsize
   endif
   if (present(err)) then
      err = 0
   endif
   if (nocaseCompare(matrix_type,'Tau_a')) then
      if (atom > 0) then
         mat => CPAMedium(n)%CPAMatrix(atom)%tau_a
      else
         mat => CPAMedium(n)%tau_c
      endif
   else if (nocaseCompare(matrix_type,'Kau_a') .and. atom > 0) then
      mat => CPAMedium(n)%CPAMatrix(atom)%kau_a
   else 
      nullify(mat)
      call ErrorHandler('getImpurityMatrix','The requested matrix does not exist', &
                         matrix_type)
   endif
!
   end function getImpurityMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine iterateCPAMedium(n)
!  ===================================================================
!
!  calculate Tcpa for the next iteration.
!
!  ===================================================================
   use SSSolverModule, only : getScatteringMatrix
!
   use CrystalMatrixModule, only : getCPAMediumTau => getTau
!
   use MatrixInverseModule, only : MtxInv_LU
!
   use MatrixModule, only : computeAprojB
!
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: id, ia, i, nsize
!
   complex (kind=CmplxKind) :: cfac
   complex (kind=CmplxKind), pointer :: xmat_a(:,:), xmat_c(:,:)
   complex (kind=CmplxKind), pointer :: xmat_proj(:,:)
!
!*********************************************************************
!  complex (kind=CmplxKind), allocatable, target :: tab(:,:), xab(:,:), wspace(:)
!  real (kind=RealKind) :: atcon(10), errcpa
!*********************************************************************
!
   id = CPAMedium(n)%local_index
   nsize = CPAMedium(n)%dsize*nSpinCant
!*********************************************************************
!! The following lines of the code are used for testing purpose.
!! Note mnewtc is the way of iterating tc in the MCPA code.
!! ===================================================================
!  allocate(tab(nsize*nsize,CPAMedium(n)%num_species))
!  allocate(xab(nsize*nsize,CPAMedium(n)%num_species))
!  allocate(wspace(3*nsize+4*nsize*nsize))
!  do ia = 1, CPAMedium(n)%num_species
!  call zcopy(nsize*nsize,CPAMedium(n)%CPAMatrix(ia)%tmat_a,1,tab(:,ia),1)
!  atcon(ia)=CPAMedium(n)%CPAMatrix(ia)%content
!  enddo
!  CPAMedium(n)%tau_c => getCPAMediumTau(id)
!  if (iteration == 1) then
!  call writeMatrix('Initial Tau_c',CPAMedium(n)%tau_c(:,:,1),        &
!                   CPAMedium(n)%dsize,CPAMedium(n)%dsize,TEN2m8)
!  endif
!  call mnewtc(CPAMedium(n)%tau_c,CPAMedium(n)%Tcpa,tab,xab,wspace,nsize,&
!              atcon,CPAMedium(n)%num_species,CPAMedium(n)%num_species,errcpa, &
!              0,'main      ')
!  CPAMedium(n)%TcpaInv=CPAMedium(n)%Tcpa
!  call MtxInv_LU(CPAMedium(n)%TcpaInv,nsize)
!  write(6,'(a,d15.8)')'After calling mnewtc, errcpa = ',errcpa
!  deallocate(tab,xab,wspace)
!  return
!*********************************************************************
   xmat_a => aliasArray2_c(WORK1,nsize,nsize)
   xmat_c => aliasArray2_c(WORK2,nsize,nsize)
   xmat_c = CZERO
   do ia = 1, CPAMedium(n)%num_species
!     ================================================================
!     Substitute one CPA medium site by a real atom. The returning
!     xmat_a is the X_a matrix if compute_X=.true., or the tau_a matrix
!     if compute_X=.false.
!     ----------------------------------------------------------------
      call substituteTcByTa(id,ia,spin=1,mat_a=xmat_a,compute_X=.true.)
!     ----------------------------------------------------------------
!
!     ================================================================
!     add content*Xmat_a to Xmat_c
!     ================================================================
      cfac = -CPAMedium(n)%CPAMatrix(ia)%content
!     ----------------------------------------------------------------
      call zaxpy(nsize*nsize,cfac,xmat_a,1,xmat_c,1)
!     ----------------------------------------------------------------
   enddo
   CPAMedium(n)%tau_c => getCPAMediumTau(id)
!  if (iteration == 1) then
!     call writeMatrix('Initial Tau_c',CPAMedium(n)%tau_c(:,:,1),     &
!                      CPAMedium(n)%dsize,CPAMedium(n)%dsize,TEN2m8)
!     call writeMatrix('Initial t_cpa',CPAMedium(n)%Tcpa,nsize,nsize,TEN2m8)
!  endif
   xmat_proj => aliasArray2_c(WORK0,nsize,nsize)
!  ===================================================================
!  This only works for non-spin-canted case.
!  -------------------------------------------------------------------
   call computeAprojB('L',nsize,xmat_c,CPAMedium(n)%tau_c(:,:,1),xmat_proj)
!  -------------------------------------------------------------------
   CPAMedium(n)%TcpaInv = CPAMedium(n)%TcpaInv - xmat_proj
   CPAMedium(n)%Tcpa = CPAMedium(n)%TcpaInv
!  -------------------------------------------------------------------
   call MtxInv_LU(CPAMedium(n)%Tcpa,nsize)
!  -------------------------------------------------------------------
!  call writeMatrix('Iterated t_cpa',CPAMedium(n)%Tcpa,            &
!                   nsize,nsize,TEN2m8)
!
   end subroutine iterateCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine substituteTcByTa(id,ia,spin,mat_a,compute_X)
!  ===================================================================
   use SSSolverModule, only : getScatteringMatrix
!
   use EmbeddedClusterModule, only : embedScatterInHostMedium,        &
                                     beginEmbedding, endEmbedding,    &
                                     calEmbeddedSiteMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, spin
!
   logical, intent(in), optional :: compute_X
   logical :: X
!
   complex (kind=CmplxKind), intent(out) :: mat_a(:,:)
   complex (kind=CmplxKind), pointer :: tmat_a_inv(:,:)
!
   if (present(compute_X)) then
      X = compute_X
   else
      X = .false.
   endif
!
   tmat_a_inv => getScatteringMatrix('TInv-Matrix',spin=spin,site=id,atom=ia)
!
!  ===================================================================
!  substitute tmat_c by tmat_a in the medium
!  Note: This should be performed in the global spin frame
!  -------------------------------------------------------------------
   call beginEmbedding()
   call embedScatterInHostMedium(id,tmat_a_inv,isInverse=.true.)
   call endEmbedding()
   call calEmbeddedSiteMatrix(id,mat_a,compute_X=X)
!  -------------------------------------------------------------------
!
   end subroutine substituteTcByTa
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSingleSiteTmat(sm_type,spin,site,atom,nsize) result(tmat)
!  ===================================================================
   use SpinRotationModule, only : rotateLtoG
!
   use SSSolverModule, only : getScatteringMatrix
!
   implicit none
!
   character (len=*), intent(in) :: sm_type
   integer (kind=IntKind), intent(in), optional :: spin, site, atom  ! Local atom index
   integer (kind=IntKind), intent(out), optional :: nsize
   integer (kind=IntKind) :: msize, ic, n
!
   complex (kind=CmplxKind), pointer :: tmat(:,:)  ! t-matrix
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (.not.nocaseCompare(sm_type,'T-Matrix') .and.                   &
       .not.nocaseCompare(sm_type,'TInv-Matrix')) then
      call ErrorHandler('getSingleSiteTmat',                          &
           'The scattering matrix type is not recognized in this case',sm_type)
   else if (.not.present(site) .or. .not.present(atom)) then
      call ErrorHandler('getSingleSiteTmat',                          &
           'site and atom arguments need to be passed in from the calling routine')
   endif
!
   n = MediumIndex(site)
   if (n > 0) then
      if (atom < 0 .or. atom > CPAMedium(n)%num_species) then
         call ErrorHandler('getSingleSiteTmat','The species index is out of range',atom)
      else if (atom == 0) then ! For a CPA site, it returns Tcpa or TcpaInv if atom = 0
         if (nocaseCompare(sm_type,'T-Matrix')) then
            tmat => CPAMedium(n)%Tcpa
         else
            tmat => CPAMedium(n)%TcpaInv
         endif
      else
         if (nocaseCompare(sm_type,'T-Matrix')) then
            tmat => CPAMedium(n)%CPAMatrix(atom)%tmat_a
         else
            tmat => getScatteringMatrix('TInv-Matrix',spin=1,site=site,atom=atom)
         endif
      endif
      if (present(nsize)) then
         nsize = CPAMedium(n)%dsize*nSpinCant
      endif
   else 
      if (atom == 0) then
         ic = 1   ! For a non-CPA site, index atom is irrelevant
      else
         ic = atom  ! For a CPA site, set ic to a real species index
      endif
      if (present(spin)) then
         tmat => getScatteringMatrix(sm_type,spin=spin,site=site,atom=ic,dsize=msize)
      else
         tmat => getScatteringMatrix(sm_type,spin=1,site=site,atom=ic,dsize=msize)
      endif
      if (present(nsize)) then
         nsize = msize
      endif
   endif
!
   end function getSingleSiteTmat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSingleSiteMatrix(sm_type,spin,site,atom,nsize) result(ss_mat)
!  ===================================================================
   use SpinRotationModule, only : rotateLtoG
!
   use SSSolverModule, only : getScatteringMatrix
!
   implicit none
!
   character (len=*), intent(in) :: sm_type
   integer (kind=IntKind), intent(in), optional :: spin, site, atom  ! Local atom index
   integer (kind=IntKind), intent(out), optional :: nsize
   integer (kind=IntKind) :: msize, ic, n
!
   complex (kind=CmplxKind), pointer :: ss_mat(:,:)
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (.not.present(site) .or. .not.present(atom)) then
      call ErrorHandler('getSingleSiteMatrix',                          &
           'site and atom arguments need to be passed in from the calling routine')
   endif
!
   n = MediumIndex(site)
   if (n > 0 .and. nocaseCompare(sm_type,'T-Matrix')) then
      if (atom < 0 .or. atom > CPAMedium(n)%num_species) then
         call ErrorHandler('getSingleSiteMatrix','The species index is out of range',atom)
      else if (atom == 0) then ! For a CPA site, it returns Tcpa if atom = 0
         ss_mat => CPAMedium(n)%Tcpa
      else
         ss_mat => CPAMedium(n)%CPAMatrix(atom)%tmat_a
      endif
      if (present(nsize)) then
         nsize = CPAMedium(n)%dsize*nSpinCant
      endif
   else 
      if (atom == 0) then
         ic = 1   ! For a non-CPA site, index atom is irrelevant
      else
         ic = atom  ! For a CPA site, set ic to a real species index
      endif
      if (present(spin)) then
         ss_mat => getScatteringMatrix(sm_type,spin=spin,site=site,atom=ic,dsize=msize)
      else
         ss_mat => getScatteringMatrix(sm_type,spin=1,site=site,atom=ic,dsize=msize)
      endif
      if (present(nsize)) then
         nsize = msize
      endif
   endif
!
   end function getSingleSiteMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkCPAMedium(n,err)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i, j, dsize, nsize
!
   real (kind=RealKind), intent(out) :: err
   real (kind=RealKind) :: trace
!
   complex (kind=CmplxKind), pointer :: tc_new(:,:), tc_old(:,:)
!
   tc_new => CPAMedium(n)%Tcpa
   tc_old => CPAMedium(n)%Tcpa_old
   dsize = CPAMedium(n)%dsize
   nsize = dsize*nSpinCant
!
   err = ZERO; trace = ZERO
   do j = 1, nsize
      do i = 1, nsize
!        err = err + abs(tc_new(i,j) - tc_old(i,j))
         err = max(err,abs(tc_new(i,j) - tc_old(i,j)))
      enddo
      trace = trace + abs(tc_old(j,j))
   enddo
   if (trace < TEN2m6) then
!     ----------------------------------------------------------------
      call ErrorHandler('checkCPAMedium','trace of Tcpa is too small',trace)
!     ----------------------------------------------------------------
   endif
!
!  err = err/trace
!
   end subroutine checkCPAMedium
!  ===================================================================
end module CPAMediumModule
