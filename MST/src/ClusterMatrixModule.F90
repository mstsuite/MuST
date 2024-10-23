!  *******************************************************************
!  *                                                                 *
!  *  version 1.0, Oct 31st, 2013                                    *
!  *     Implemented based on TauLSMSModule and TauKKRModule         *
!  *                                                                 *
!  *******************************************************************
module ClusterMatrixModule
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only  : sortNeighbors, getNeighbor
   use NeighborModule, only  : getNumReceives, getMaxReceives
   use MathParamModule, only : czero, cone, zero
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use WriteMatrixModule,  only : writeMatrix
!
public :: initClusterMatrix, &
          endClusterMatrix,  &
          calClusterMatrix,  &
          calClusterMatrixNonPeriodic, &
          getTau,            &
          checkIfNeighbor,   &
          getClusterTau,     &
          getNeighborTau,    &
          getKau
!
private
!
   integer (kind=IntKind) :: Relativity
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: MtxInvAlg = 2
   integer (kind=IntKind) :: dsize_max = 0
   integer (kind=IntKind) :: isize_max = 0
   integer (kind=IntKind) :: numnb_max = 0
   integer (kind=IntKind) :: MaxPrintLevel = -1
!
   integer (kind=IntKind) ::  lmax_max
   integer (kind=IntKind) ::  kmax_max
!
   integer (kind=IntKind), allocatable :: print_level(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: kmax_kkr(:)
   integer (kind=IntKind), allocatable :: kmax_phi(:)
   integer (kind=IntKind), allocatable :: ipvt(:)
   integer (kind=IntKind), allocatable :: dsize_l(:)
   integer (kind=IntKind), allocatable :: kmax_nrs(:)
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
!
   real (kind=RealKind), allocatable :: Position(:,:)
!
   complex (kind=CmplxKind), allocatable, target :: jinv_g(:)
   complex (kind=CmplxKind), allocatable, target :: sine_g(:)
   complex (kind=CmplxKind), allocatable, target :: oh_g(:)
   complex (kind=CmplxKind), allocatable, target :: BlockMatrix(:)
   complex (kind=CmplxKind), allocatable, target :: BigMatrix(:)
   complex (kind=CmplxKind), allocatable, target :: BigMatrixInv(:,:)
   complex (kind=CmplxKind), allocatable, target :: ClusterMatrix(:,:)
   complex (kind=CmplxKind), allocatable, target :: ClusterKau(:,:)
   complex (kind=CmplxKind), allocatable, target :: ClusterTau(:,:)
!
   type TauNeighborStruct
      complex (kind=CmplxKind), allocatable :: wau_l(:,:)
      complex (kind=CmplxKind), allocatable :: kau_l(:,:)
      complex (kind=CmplxKind), allocatable :: tau_l(:,:)
   end type TauNeighborStruct

   type TauijStruct
      integer (kind=IntKind) :: id_cluster  ! the index within LIZ cluster
      integer (kind=IntKind) :: id_global   ! the atom global index
      integer (kind=IntKind) :: isize
      integer (kind=IntKind) :: jsize
      complex (kind=CmplxKind), pointer :: kau_l(:,:,:)
      complex (kind=CmplxKind), pointer :: tau_l(:,:,:)
      type (TauNeighborStruct), allocatable :: neigh1j(:)
      type (TauNeighborStruct), allocatable :: neighj1(:)
      type (TauijStruct), pointer :: next
   end type TauijStruct
!
   type (TauijStruct), allocatable, target :: Tau00(:)
!
   type (NeighborStruct), pointer :: Neighbor
!
   complex (kind=CmplxKind), allocatable, target :: wsTau00L(:)
   complex (kind=CmplxKind), allocatable, target :: wsKau00L(:)
   complex (kind=CmplxKind), allocatable, target :: wsStoreA(:)
   complex (kind=CmplxKind), allocatable, target :: gijrel(:,:)
!
   logical :: Initialized = .false.
   logical :: InitializedFactors = .false.
!
   logical :: ExchSSSmatAllocated = .false.
   logical :: istauij_needed = .false.
!
   complex (kind=CmplxKind), allocatable, target :: remote_jsmtx(:,:)
   complex (kind=CmplxKind), allocatable, target :: local_jsmtx(:,:)
   complex (kind=CmplxKind), allocatable, target :: store_sine(:)
   complex (kind=CmplxKind), allocatable :: gij(:,:)
#ifdef OpenMPI
   complex (kind=CmplxKind), allocatable, target :: jsmtx_buf(:,:,:)
#else
   complex (kind=CmplxKind), allocatable :: jsmtx_buf(:,:)
#endif
!
   integer (kind=IntKind) :: NumNonLocalNeighbors
   integer (kind=IntKind) :: NumReceives
   integer (kind=IntKind) :: NumSends
   integer (kind=IntKind), pointer :: NonLocalNeighbors(:)
   integer (kind=IntKind), pointer :: TargetProc(:)
   integer (kind=IntKind), pointer :: SourceProc(:)
   integer (kind=IntKind), allocatable :: remote_kmax(:)
   integer (kind=IntKind), allocatable :: mapping_jsmtx(:,:)
   integer (kind=IntKind), allocatable :: recv_msgid(:), send_msgid(:)
!
   character (len=15) :: stop_routine
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initClusterMatrix(num_latoms,index,lmaxkkr,lmaxphi,     &
                                posi,cant,rel,istop,iprint)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   use ScfDataModule, only : Minv_alg, NumEs, tauij_needed
!
   use SystemModule, only  : getLmaxPhi
!
   use RSpaceStrConstModule, only : initRSpaceStrConst
!
   use MatrixBlockInversionModule, only : initMatrixBlockInv
   use WriteMatrixModule, only : writeMatrix
!
   integer (kind=IntKind), intent(in) :: num_latoms
   integer (kind=IntKind), intent(in) :: index(num_latoms)
   integer (kind=IntKind), intent(in) :: lmaxkkr(num_latoms)
   integer (kind=IntKind), intent(in) :: lmaxphi(num_latoms)
   integer (kind=IntKind), intent(in) :: cant
   integer (kind=IntKind), intent(in) :: rel
   integer (kind=IntKind), intent(in) :: iprint(num_latoms)
!
   character (len=*), intent(in) :: istop
!
   real (kind=RealKind), intent(in) :: posi(3,num_latoms)
!
   integer (kind=IntKind) :: lmax, i, j, max_nblck, kmaxns, bsz
   integer (kind=IntKind) :: wsTau00_size, pos_tau
   integer (kind=IntKind), allocatable :: nblck(:),kblocks(:,:)
!
   complex (kind=CmplxKind), pointer :: p1(:)
!
   n_spin_cant = cant
   Relativity = rel
   LocalNumAtoms = num_latoms
   stop_routine = istop
!
   if (tauij_needed == 1) then
     istauij_needed = .true.
   endif

   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   allocate( GlobalIndex(LocalNumAtoms), Position(3,LocalNumAtoms) )
   allocate( lmax_kkr(LocalNumAtoms), lmax_phi(LocalNumAtoms) )
   allocate( kmax_kkr(LocalNumAtoms), kmax_phi(LocalNumAtoms) )
   allocate( Tau00(LocalNumAtoms) )
   allocate( dsize_l(LocalNumAtoms), kmax_nrs(LocalNumAtoms) )
   allocate( print_level(LocalNumAtoms) )
!
   print_level(1:LocalNumAtoms) = iprint(1:LocalNumAtoms)
   MaxPrintLevel = maxval(print_level)
   if (MaxPrintLevel >= 0) then
      write(6,'(/,80("="))')
   endif
!
!  ===================================================================
!  Initialize the structure constant module.
!  -------------------------------------------------------------------
   call initRSpaceStrConst(getLmaxPhi(),istop,MaxPrintLevel)  
!  -------------------------------------------------------------------
!
   wsTau00_size = 0
   allocate( nblck(LocalNumAtoms) )
   lmax_max = 0
   kmax_max = 0
   do i=1,LocalNumAtoms
      GlobalIndex(i) = index(i)
      Position(1:3,i) = posi(1:3,i)
      lmax_kkr(i) = lmaxkkr(i)
      lmax_phi(i) = lmaxphi(i)
      kmax_kkr(i) = (lmaxkkr(i)+1)**2
      kmax_phi(i) = (lmaxphi(i)+1)**2
      wsTau00_size = wsTau00_size + kmax_kkr(i)*kmax_kkr(i)
      Tau00(i)%id_global = index(i)
      Tau00(i)%id_cluster = 0
      Tau00(i)%isize = kmax_kkr(i)*n_spin_cant
      Tau00(i)%jsize = kmax_kkr(i)*n_spin_cant
      nullify( Tau00(i)%tau_l, Tau00(i)%kau_l, Tau00(i)%next )
      Neighbor => getNeighbor(i)
      dsize_l(i) = kmax_kkr(i)
      kmax_nrs(i) = kmax_kkr(i)
      lmax_max = max(lmax_max,lmax_kkr(i))
      do j = 1,Neighbor%NumAtoms
         lmax = Neighbor%Lmax(j)
         lmax_max = max(lmax_max,lmax)
         dsize_l(i) = dsize_l(i) + (lmax+1)*(lmax+1)
         kmax_nrs(i) = max(kmax_nrs(i),(lmax+1)*(lmax+1))
      enddo
      kmax_max = max(kmax_max,kmax_nrs(i))
      dsize_l(i) = dsize_l(i)*n_spin_cant
      nblck(i) = Neighbor%NumAtoms+1
      if (print_level(i) >= 0) then
         write(6,'(a,i5)')           'My global index  = ',index(i)
         write(6,'(a,i5,'' x '',i5)')'My T-Matrix Size = ',            &
                       kmax_kkr(i)*n_spin_cant,kmax_kkr(i)*n_spin_cant
         write(6,'(a,i5)')           'LIZ cluster size = ',nblck(i)
         write(6,'(a,i5,'' x '',i5)')'Tau-Matrix Size  = ',            &
                                                 dsize_l(i),dsize_l(i)
      endif
   enddo
!
   if (MaxPrintLevel >= 0) then
      write(6,'(80("="),/)')
   endif
!
   max_nblck = maxval( nblck )
   max_nblck = max(max_nblck,2) ! This takes care the case when LIZ=1
   allocate( kblocks(max_nblck,LocalNumAtoms) )
!
   wsTau00_size = wsTau00_size*n_spin_cant*n_spin_cant
   allocate( wsKau00L(wsTau00_size), wsTau00L(wsTau00_size) )
   allocate( wsStoreA(2*kmax_max*kmax_max*n_spin_cant*n_spin_cant) )
!
   pos_tau = 0
   numnb_max = 0
   do i = 1,LocalNumAtoms
      kmaxns = kmax_kkr(i)*n_spin_cant
      kblocks(1,i) = kmaxns
      Neighbor => getNeighbor(i)
      if (istauij_needed) then
        allocate(Tau00(i)%neigh1j(Neighbor%NumAtoms+1))
        allocate(Tau00(i)%neighj1(Neighbor%NumAtoms+1))
      endif
      do j = 1,Neighbor%NumAtoms
         lmax = Neighbor%Lmax(j)
         lmax_max = max(lmax_max,lmax)
         kblocks(j+1,i) = (lmax+1)*(lmax+1)*n_spin_cant
         if (istauij_needed) then
           allocate(Tau00(i)%neigh1j(j)%tau_l(kmax_max, kmax_max))
           allocate(Tau00(i)%neigh1j(j)%kau_l(kmax_max, kmax_max))
           allocate(Tau00(i)%neigh1j(j)%wau_l(kmax_max, kmax_max))
           allocate(Tau00(i)%neighj1(j)%tau_l(kmax_max, kmax_max))
           allocate(Tau00(i)%neighj1(j)%kau_l(kmax_max, kmax_max))
           allocate(Tau00(i)%neighj1(j)%wau_l(kmax_max, kmax_max))
           Tau00(i)%neigh1j(j)%tau_l = CZERO
           Tau00(i)%neigh1j(j)%kau_l = CZERO
           Tau00(i)%neigh1j(j)%wau_l = CZERO
           Tau00(i)%neighj1(j)%tau_l = CZERO
           Tau00(i)%neighj1(j)%kau_l = CZERO
           Tau00(i)%neighj1(j)%wau_l = CZERO
         endif
      enddo
      if (istauij_needed) then
        allocate(Tau00(i)%neigh1j(Neighbor%NumAtoms+1)%tau_l(kmax_max, kmax_max))
        allocate(Tau00(i)%neigh1j(Neighbor%NumAtoms+1)%kau_l(kmax_max, kmax_max))
        allocate(Tau00(i)%neigh1j(Neighbor%NumAtoms+1)%wau_l(kmax_max, kmax_max))
        allocate(Tau00(i)%neighj1(Neighbor%NumAtoms+1)%tau_l(kmax_max, kmax_max))
        allocate(Tau00(i)%neighj1(Neighbor%NumAtoms+1)%kau_l(kmax_max, kmax_max))
        allocate(Tau00(i)%neighj1(Neighbor%NumAtoms+1)%wau_l(kmax_max, kmax_max))
      endif
!     Tau00(i)%tau_l => aliasArray3_c( wsTau00L(pos_tau+1:pos_tau+kmaxns*kmaxns),&
!                                      kmax_kkr(i), kmax_kkr(i), n_spin_cant*n_spin_cant )
      p1 => wsTau00L(pos_tau+1:pos_tau+kmaxns*kmaxns)
      Tau00(i)%tau_l => aliasArray3_c( p1, kmax_kkr(i), kmax_kkr(i), n_spin_cant*n_spin_cant )
!     Tau00(i)%kau_l => aliasArray3_c( wsKau00L(pos_tau+1:pos_tau+kmaxns*kmaxns),&
!                                      kmax_kkr(i), kmax_kkr(i), n_spin_cant*n_spin_cant )
      p1 => wsKau00L(pos_tau+1:pos_tau+kmaxns*kmaxns)
      Tau00(i)%kau_l => aliasArray3_c( p1, kmax_kkr(i), kmax_kkr(i), n_spin_cant*n_spin_cant )
      pos_tau = pos_tau + kmaxns*kmaxns
      numnb_max = max(numnb_max, Neighbor%NumAtoms)
   enddo
!
   dsize_max = maxval( dsize_l )
   isize_max = maxval( kmax_nrs )*n_spin_cant
!
   allocate( jinv_g(isize_max*isize_max) )
   allocate( sine_g(isize_max*isize_max) )
   allocate( store_sine(isize_max*isize_max) )
   allocate( BigMatrix(dsize_max*dsize_max) )
!
   jinv_g = CZERO; sine_g = CZERO
!
!  Allocate matrices for tauij calculation (i,j != 0)
   if (istauij_needed) then
      allocate( oh_g(isize_max*isize_max) )
      oh_g = CZERO
!    allocate (BigMatrixInv(dsize_max, dsize_max))
!    allocate (ClusterMatrix(LocalNumAtoms*kmax_max,LocalNumAtoms*kmax_max))
!    allocate (ClusterKau(LocalNumAtoms*kmax_max,LocalNumAtoms*kmax_max))
!    allocate (ClusterTau(LocalNumAtoms*kmax_max,LocalNumAtoms*kmax_max))
   endif
!
!  -------------------------------------------------------------------
   call initMatrixBlockInv( Minv_alg, LocalNumAtoms, NumEs, nblck,    &
                           max_nblck, kblocks, print_level ) 
!  -------------------------------------------------------------------
!
   allocate( ipvt(isize_max) )
   allocate( BlockMatrix(isize_max*isize_max) )
!
   deallocate( nblck, kblocks )
   nullify( Neighbor )
!
   if (Relativity > 1) then
      allocate( gijrel(kmax_max*n_spin_cant, kmax_max*n_spin_cant) )
   else
      allocate( gij(kmax_max, kmax_max) )
   endif
!
   ExchSSSmatAllocated = .false.
   Initialized = .true.
!
   end subroutine initClusterMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endClusterMatrix()
!  ===================================================================
   use MatrixBlockInversionModule, only : endMatrixBlockInv
!
   use RSpaceStrConstModule, only : endRSpaceStrConst
!
   implicit none
!
   integer (kind=IntKind) :: i, j
!
   if (.not.Initialized) then
      call ErrorHandler('endTauLSMS','Module not initialized')
   endif
!
   deallocate( GlobalIndex, Position, lmax_kkr, lmax_phi, kmax_kkr )
   deallocate( kmax_phi, kmax_nrs, dsize_l, print_level )
!
   do i=1,LocalNumAtoms
      if ( associated(Tau00(i)%tau_l) ) then
         nullify( Tau00(i)%tau_l )
      endif
   enddo
   deallocate( wsStoreA )
   deallocate( wsKau00L, wsTau00L )
!
   deallocate( ipvt, jinv_g, sine_g, store_sine )
   deallocate( BlockMatrix, BigMatrix )
   deallocate( Tau00 )
!
   if (istauij_needed) then
      deallocate(oh_g)
      istauij_needed = .false.
   endif
!
   call endMatrixBlockInv()
   nullify( Neighbor )
!
   if ( InitializedFactors ) then
      InitializedFactors = .false.
   endif
!
   if (ExchSSSmatAllocated) then
      call cleanExchSSSmat()
   endif
!
   if (Relativity > 1) then
      deallocate( gijrel )
   else
      deallocate( gij )
   endif
!
!  -------------------------------------------------------------------
   call endRSpaceStrConst()
!  -------------------------------------------------------------------
!
   Initialized = .false.
!
   end subroutine endClusterMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocateExchSSSmat(bsize)
!  ===================================================================
   use SendRecvTmatModule, only : getMaxNeighbors
   use SendRecvTmatModule, only : getNumNonLocalNeighbors, getNonLocalNeighbors
   use SendRecvTmatModule, only : getNumProcWiseSends, getNumProcWiseReceives
   use SendRecvTmatModule, only : getProcWiseTargetProcs, getProcWiseSourceProcs
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: bsize
   integer (kind=IntKind) :: m, n, i, j, iri, pe, gid
!
!  ===================================================================
!  We assume that bsize is the same on all processors. This needs to be
!  fixed in the future.
!  ===================================================================
   m = getMaxNeighbors()
   NumNonLocalNeighbors = getNumNonLocalNeighbors()
!
   NonLocalNeighbors => getNonLocalNeighbors()
   NumSends = getNumProcWiseSends()
   NumReceives = getNumProcWiseReceives()
   TargetProc => getProcWiseTargetProcs()
   SourceProc => getProcWiseSourceProcs()
!
   allocate( local_jsmtx(bsize+1,LocalNumAtoms) )
   allocate(recv_msgid(NumReceives),send_msgid(NumSends))
#ifdef OpenMPI
   allocate( jsmtx_buf(bsize+1,LocalNumAtoms,NumReceives) )
#else
   allocate( jsmtx_buf(bsize+1,LocalNumAtoms) )
#endif
   allocate( mapping_jsmtx(m,LocalNumAtoms) )
!
   mapping_jsmtx = -1
   local_jsmtx = CZERO
   jsmtx_buf = CZERO
!
   do i = 1, LocalNumAtoms
      Neighbor => getNeighbor(i)
      do iri=1, Neighbor%NumAtoms
         pe = Neighbor%ProcIndex(iri)
         if (pe /= MyPEinGroup) then
            gid = Neighbor%GlobalIndex(iri)
            LOOP_j: do j = 1, NumNonLocalNeighbors
               if (gid == NonLocalNeighbors(j)) then
                  mapping_jsmtx(iri,i) = j
                  exit LOOP_j
               endif
            enddo LOOP_j
         endif
      enddo
   enddo
   allocate( remote_kmax(NumNonLocalNeighbors), remote_jsmtx(bsize,NumNonLocalNeighbors) )
!
   ExchSSSmatAllocated = .true.
!
   end subroutine allocateExchSSSmat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cleanExchSSSmat()
!  ===================================================================
   implicit none
!
   NumNonLocalNeighbors = 0
   NumSends = 0
   NumReceives = 0
!
   deallocate(remote_kmax, remote_jsmtx, local_jsmtx, jsmtx_buf)
   deallocate(mapping_jsmtx)
   deallocate(recv_msgid,send_msgid)
!
   nullify(NonLocalNeighbors, SourceProc, TargetProc)
!
   ExchSSSmatAllocated = .false.
!
   end subroutine cleanExchSSSmat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine exchangeSSSMatrix(getSingleScatteringMatrix)
!  ===================================================================
   use MPPModule, only : AnyPE, MyPE
   use MPPModule, only : nbsendMessage, nbrecvMessage
   use MPPModule, only : recvMessage, isThereMessage, getSource, setMaxWaits
   use MPPModule, only : waitMessage, setCommunicator, resetCommunicator
!
   use GroupCommModule, only : getGroupCommunicator
!
   use DataServiceCenterModule, only : getDataStorage, ComplexMark
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use SpinRotationModule, only : rotateLtoG
!
   implicit none
!
   integer (kind=IntKind) :: comm, nfac
   integer (kind=IntKind) :: max_tsize, kmax_max_ns, kkr_ns, kli, klj
   integer (kind=IntKind) :: i, j, k, n, nl, nr, np, ns, nm, kmax, NumFills
   integer (kind=IntKind) :: proc, gid, t0size, t0size_ns, kkri_ns
!
   complex (kind=CmplxKind), pointer :: sm1(:,:), sm2(:,:), gmat(:,:)
   complex (kind=CmplxKind), pointer :: jm1(:,:), jm2(:,:), p1(:)
   complex (kind=CmplxKind), pointer :: om1(:,:), om2(:,:)
#ifdef OpenMPI
   complex (kind=CmplxKind), pointer :: p_trecv(:,:)
#endif
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('exchangeSSSMatrix','Module not initialized')
!     ----------------------------------------------------------------
   endif
!
   max_tsize = kmax_max*kmax_max*n_spin_cant*n_spin_cant
   kmax_max_ns = kmax_max*n_spin_cant
!
   if (istauij_needed) then
      nfac = 3
   else
      nfac = 2
   endif
!
   if (.not. ExchSSSmatAllocated) then
!     ----------------------------------------------------------------
      call allocateExchSSSmat(nfac*max_tsize) ! Allocate space for exchanging Jinv and Sine Matrices
                                              ! as well as the space for OmegaHat matrix if tauij is needed
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   comm = getGroupCommunicator(GroupID)
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!  -------------------------------------------------------------------
!
   do i = 1, LocalNumAtoms
!     ================================================================
!     Obtain the Jinv-matrix and the Sine-Matrix
!     ================================================================
      t0size = kmax_kkr(i)*kmax_kkr(i)
      if ( n_spin_cant == 2 ) then
!        =============================================================
!        calculate jinv_g and sine_g in global frame of reference.
!        =============================================================
         t0size_ns = t0size*n_spin_cant*n_spin_cant
         kkri_ns =  kmax_kkr(i)*n_spin_cant
         jm1 => getSingleScatteringMatrix('JostInv-Matrix',spin=1,site=i)
         jm2 => getSingleScatteringMatrix('JostInv-Matrix',spin=2,site=i)
         gmat => aliasArray2_c(jinv_g,kkri_ns,kkri_ns)
!        -------------------------------------------------------------
         call rotateLtoG(i, kmax_kkr(i), kmax_kkr(i), jm1, jm2, gmat)
!        -------------------------------------------------------------
         sm1 => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=i)
         sm2 => getSingleScatteringMatrix('Sine-Matrix',spin=2,site=i)
!        gmat => aliasArray2_c(sine_g(1:t0size_ns),kkri_ns,kkri_ns)
         p1 => sine_g(1:t0size_ns)
         gmat => aliasArray2_c(p1,kkri_ns,kkri_ns)
!        -------------------------------------------------------------
         call rotateLtoG(i, kmax_kkr(i), kmax_kkr(i), sm1, sm2, gmat)
!        -------------------------------------------------------------
         if (istauij_needed) then
            om1 => getSingleScatteringMatrix('OmegaHat-Matrix',spin=1,site=i)
            om2 => getSingleScatteringMatrix('OmegaHat-Matrix',spin=2,site=i)
            gmat => aliasArray2_c(oh_g,kkri_ns,kkri_ns)
!           ----------------------------------------------------------
            call rotateLtoG(i, kmax_kkr(i), kmax_kkr(i), om1, om2, gmat)
!           ----------------------------------------------------------
         endif
      else
         jm1 => getSingleScatteringMatrix('JostInv-Matrix',spin=1,site=i)
         sm1 => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=i)
!        -------------------------------------------------------------
         call zcopy( t0size, jm1, 1, jinv_g, 1 )
         call zcopy( t0size, sm1, 1, sine_g, 1 )
!        -------------------------------------------------------------
         if (istauij_needed) then
            om1 => getSingleScatteringMatrix('OmegaHat-Matrix',spin=1,site=i)
            call zcopy( t0size, om1, 1, oh_g, 1 )
         endif
      endif
!
!     ================================================================
!     Store the jinv-matrix and sine-matrix into the local communication buffer
!     ================================================================
      local_jsmtx(1,i) = cmplx( kmax_kkr(i), getGlobalIndex(i), kind=CmplxKind)
!     ----------------------------------------------------------------
      call zcopy(max_tsize,jinv_g,1,local_jsmtx(2,i),1)
      call zcopy(max_tsize,sine_g,1,local_jsmtx(2+max_tsize,i),1)
!     ----------------------------------------------------------------
      if (istauij_needed) then
!        -------------------------------------------------------------
         call zcopy(max_tsize,oh_g,1,local_jsmtx(2+2*max_tsize,i),1)
!        -------------------------------------------------------------
      endif
   enddo
!
   nullify(sm1,sm2,jm1,jm2,gmat)
!
   nr = 0; ns = 0; NumFills = 0
   remote_kmax = 0
!
#ifdef OpenMPI
   nm = max(NumReceives, NumSends)
   call setMaxWaits(nm)
   do k = 1, nm
      if (k <= NumReceives) then
         p_trecv => jsmtx_buf(1:nfac*max_tsize+1,1:LocalNumAtoms,k)
!        -------------------------------------------------------------
         recv_msgid(k)=nbrecvMessage(p_trecv,nfac*max_tsize+1,LocalNumAtoms,1010,AnyPE)
!        ------------------------------------------------------------
      endif
   enddo
!
   do k = 1, nm
      if (k < NumSends+1) then
         proc = TargetProc(k)
!        -------------------------------------------------------------
         send_msgid(k) = nbsendMessage(local_jsmtx,nfac*max_tsize+1,LocalNumAtoms,1010,proc)
!        -------------------------------------------------------------
      endif
   enddo
!
   do k = 1, nm
      if (nr < NumReceives) then
         p_trecv => jsmtx_buf(1:nfac*max_tsize+1,1:LocalNumAtoms,nr+1)
!        -------------------------------------------------------------
         call waitMessage(recv_msgid(k))
!        -------------------------------------------------------------
         nr = nr + 1
!
         do i = 1, LocalNumAtoms
            kmax = int(real(jsmtx_buf(1,i,nr),kind=RealKind))
            gid = int(aimag(jsmtx_buf(1,i,nr)))
            LOOP_j: do j = 1, NumNonLocalNeighbors
               if (gid == NonLocalNeighbors(j)) then
!                 ----------------------------------------------------
                  call zcopy(nfac*max_tsize,jsmtx_buf(2,i,nr),1,remote_jsmtx(1,j),1)
!                 ----------------------------------------------------
!                 remote_jsmtx(1:nfac*max_tsize,j) = jsmtx_buf(2:nfac*max_tsize+1,i,nr)
                  remote_kmax(j) = kmax
                  NumFills = NumFills + 1
                  exit LOOP_j
               endif
            enddo LOOP_j
         enddo
      else if (NumFills /= NumNonLocalNeighbors) then
!        -------------------------------------------------------------
         call ErrorHandler('exchangeSSSMatrix','NumFills /= NumNonLocalNeighbors', &
                           NumFills,NumNonLocalNeighbors)
!        -------------------------------------------------------------
      endif
   enddo
!
   do k = 1, nm
      if (k < NumSends+1) then
!        -------------------------------------------------------------
         call waitMessage(send_msgid(k))
!        -------------------------------------------------------------
         ns = ns + 1
      endif
   enddo
!
#else
   do while (nr < NumReceives .or. ns < NumSends)
      if (ns < NumSends) then
         ns = ns + 1
         proc = TargetProc(ns)
         send_msgid(ns) = nbsendMessage(local_jsmtx,nfac*max_tsize+1,LocalNumAtoms,1010,proc)
      endif
!
      if (nr < NumReceives) then
         if (isThereMessage(1010,AnyPE)) then
!           ----------------------------------------------------------
            call recvMessage(jsmtx_buf,nfac*max_tsize+1,LocalNumAtoms,1010,getSource())
!           ----------------------------------------------------------
            nr = nr + 1
!
            do i = 1, LocalNumAtoms
               kmax = int(real(jsmtx_buf(1,i),kind=RealKind))
               gid = int(aimag(jsmtx_buf(1,i)))
               LOOP_j: do j = 1, NumNonLocalNeighbors
                  if (gid == NonLocalNeighbors(j)) then
!                    -------------------------------------------------
                     call zcopy(nfac*max_tsize,jsmtx_buf(2,i),1,remote_jsmtx(1,j),1)
!                    -------------------------------------------------
!                    remote_jsmtx(1:nfac*max_tsize,j) = jsmtx_buf(2:nfac*max_tsize+1,i)
                     remote_kmax(j) = kmax
                     NumFills = NumFills + 1
                     exit LOOP_j
                  endif
               enddo LOOP_j
            enddo
         endif
      endif
   enddo
!
   if (NumFills /= NumNonLocalNeighbors) then
!     ----------------------------------------------------------------
      call ErrorHandler('exchangeSSSMatrix','NumFills /= NumNonLocalNeighbors', &
                        NumFills,NumNonLocalNeighbors)
!     ----------------------------------------------------------------
   endif
!
   do k = 1, NumSends
      call waitMessage(send_msgid(k))
   enddo
#endif
!
!  -------------------------------------------------------------------
   call resetCommunicator(sync=.true.)
!  -------------------------------------------------------------------
!
   end subroutine exchangeSSSMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calClusterMatrixNonPeriodic(energy,getSingleScatteringMatrix,tau_needed,kp)
!  ===================================================================
   use TimerModule, only : getTime
   use MPPModule, only : MyPE, syncAllPEs
!
   use SpinRotationModule, only : rotateGtoL
!
   use RSpaceStrConstModule, only : getStrConstMatrix
!
   use MatrixModule, only : setupUnitMatrix, computeUAUts
!
   use MatrixBlockInversionModule, only : invertMatrixBlock
   use MatrixInverseModule, only : MtxInv_LU
   use WriteMatrixModule, only : writeMatrix

   implicit none
!
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind), intent(in), optional :: kp
!
   logical, intent(in), optional :: tau_needed
!
   integer (kind=IntKind) :: my_atom, iri, irj, kli, klj, is, js, ns
   integer (kind=IntKind) :: lmaxi, lmaxj, kkri, kkrj, kkri_ns, kkrj_ns
   integer (kind=IntKind) :: nbi, nbj, sbi, sbj, i, j, kl, jl, my_atom2, nindex
   integer (kind=IntKind) :: nbr_kkri, nbr_kkri_ns, nbr_kkrj, nbr_kkrj_ns
   integer (kind=IntKind) :: tsize, kmax_max_ns, dsize, kkrsz_ns, kkrsz
!
   integer (kind=IntKind) :: lid_i, lid_j, pei, pej, jid, info, iid
!
   real (kind=RealKind) :: rij(3), posi(3), posj(3)
!
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind), pointer :: p_jinvi(:)
   complex (kind=CmplxKind), pointer :: p_sinej(:)
   complex (kind=CmplxKind), pointer :: OmegaHat(:,:)
   complex (kind=CmplxKind), pointer :: wau_g(:,:)
   complex (kind=CmplxKind), pointer :: wau_l(:,:,:)
   complex (kind=CmplxKind), pointer :: jig(:), ubmat(:,:)
   complex (kind=CmplxKind), pointer :: smi(:,:), smj(:,:), tm(:,:)
   complex (kind=CmplxKind), pointer :: kau_l(:,:), tau_l(:,:)
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface

   interface
      subroutine convertGijToRel(gij, bgij, kkr1, kkr2, ce)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: kkr1, kkr2
         complex (kind=CmplxKind), intent(in) :: gij(:,:)
         complex (kind=CmplxKind), intent(out) :: bgij(:,:)
         complex (kind=CmplxKind), intent(in) :: ce
      end subroutine convertGijToRel
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('calTauMatrix','Module not initialized')
   endif
!
!  -------------------------------------------------------------------
   call exchangeSSSMatrix(getSingleScatteringMatrix)
!  -------------------------------------------------------------------
!
   if (present(kp)) then
     kappa = kp
   else
     kappa = sqrt(energy)
   endif
   kkrsz = kmax_max*n_spin_cant
   tsize = kmax_max*kmax_max*n_spin_cant*n_spin_cant

   call setupUnitMatrix(LocalNumAtoms*kkrsz, ClusterMatrix)

   do my_atom = 1, LocalNumAtoms
     jig => wsStoreA
     posj(1:3) = Position(1:3,my_atom)
     lmaxj = lmax_kkr(my_atom)
     nbr_kkrj = kmax_kkr(my_atom)
     p_sinej => local_jsmtx(tsize+2:2*tsize+1,my_atom)
     kkrj  = (lmaxj+1)*(lmaxj+1)
     kkrj_ns  = kkrj*n_spin_cant
     nbr_kkrj_ns  = nbr_kkrj*n_spin_cant
     do my_atom2 = 1, LocalNumAtoms
       posi(1:3) = Position(1:3, my_atom2)
       lmaxi = lmax_kkr(my_atom2)
       nbr_kkri = kmax_kkr(my_atom2)
       p_jinvi => local_jsmtx(2:tsize+1,my_atom2)
       kkri = (lmaxi+1)*(lmaxi+1)
       kkri_ns = kkri*n_spin_cant
       nbr_kkri_ns = nbr_kkri*n_spin_cant
       if (my_atom .ne. my_atom2) then
         rij = posj - posi
         call getStrConstMatrix(lmaxi,lmaxj,kappa,rij,gij)
         do js = 1, n_spin_cant
           do is = 1, n_spin_cant
!            ----------------------------------------------
             call zgemm('n', 'n', kkri, kkrj, kkri, CONE,                      &
                        p_jinvi((js-1)*nbr_kkri_ns*nbr_kkri+(is-1)*nbr_kkri+1),&
                        nbr_kkri_ns, gij, kkri, CZERO,                         &
                        jig((js-1)*kkri_ns*kkrj+(is-1)*kkri+1), kkri_ns)
!            ----------------------------------------------
           enddo
         enddo
!        ------------------------------------------------------
         call zgemm('n', 'n', kkri_ns, kkrj_ns, kkrj_ns, -CONE/kappa, &
                    jig, kkri_ns, p_sinej, nbr_kkrj_ns, CZERO,        &
                    ClusterMatrix((my_atom2-1)*kkrsz+1:my_atom2*kkrsz, &
                     (my_atom-1)*kkrsz+1:my_atom*kkrsz), kkrsz)
!        ------------------------------------------------------
       endif 
     enddo
   enddo
   
   call MtxInv_LU(ClusterMatrix, LocalNumAtoms*kkrsz)
  
   do i = 1, LocalNumAtoms*kkrsz
     ClusterMatrix(i,i) = ClusterMatrix(i,i) - CONE
   enddo
   ClusterKau = CZERO
   ClusterTau = CZERO

   do my_atom = 1, LocalNumAtoms
!    ----------------------------------------------------------
     OmegaHat => getSingleScatteringMatrix('OmegaHat-Matrix',spin=js,site=my_atom)
!    ----------------------------------------------------------
     do my_atom2 = 1, LocalNumAtoms
       call zgemm('n', 'n', kmax_kkr(my_atom), kmax_kkr(my_atom), &
       kmax_kkr(my_atom), kappa, ClusterMatrix((my_atom2-1)*kkrsz+1:my_atom2*kkrsz, &
       (my_atom-1)*kkrsz+1:my_atom*kkrsz), kmax_kkr(my_atom), OmegaHat, &
       kmax_kkr(my_atom), CZERO, ClusterKau((my_atom2-1)*kkrsz+1:my_atom2*kkrsz, &
       (my_atom-1)*kkrsz+1:my_atom*kkrsz), kmax_kkr(my_atom))
     enddo
   enddo

   do my_atom = 1, LocalNumAtoms
     smj => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=my_atom)
     do my_atom2 = 1, LocalNumAtoms
       smi => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=my_atom2)
       call computeUAUts(smi,kmax_kkr(my_atom),kmax_kkr(my_atom), &
                         smj,kmax_kkr(my_atom), CONE/energy, &
         ClusterKau((my_atom2-1)*kkrsz+1:my_atom2*kkrsz, &
          (my_atom-1)*kkrsz+1:my_atom*kkrsz), kmax_kkr(my_atom),CZERO, &
         ClusterTau((my_atom2-1)*kkrsz+1:my_atom2*kkrsz, &
          (my_atom-1)*kkrsz+1:my_atom*kkrsz), kmax_kkr(my_atom),wsStoreA)
     enddo
   enddo

   do my_atom = 1, LocalNumAtoms
     tm => getSingleScatteringMatrix('T-Matrix',spin=1,site=my_atom)
     ClusterTau((my_atom-1)*kkrsz+1:my_atom*kkrsz, &
          (my_atom-1)*kkrsz+1:my_atom*kkrsz) = ClusterTau((my_atom-1)*kkrsz+1:my_atom*kkrsz, &
          (my_atom-1)*kkrsz+1:my_atom*kkrsz) + tm
   enddo

!  deallocate(ClusterMatrix, ClusterKau)

   end subroutine calClusterMatrixNonPeriodic   
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calClusterMatrix(energy,getSingleScatteringMatrix,tau_needed,kp)
!  ===================================================================
   use TimerModule, only : getTime
use MPPModule, only : MyPE, syncAllPEs
!
   use SpinRotationModule, only : rotateGtoL
!
   use RSpaceStrConstModule, only : getStrConstMatrix
!
   use MatrixModule, only : setupUnitMatrix, computeUAUts
!
   use MatrixBlockInversionModule, only : invertMatrixBlock
   use MatrixInverseModule, only : MtxInv_LU
   use WriteMatrixModule, only : writeMatrix
   !
   use NVTX
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind), intent(in), optional :: kp
!
   logical, intent(in), optional :: tau_needed
!
   integer (kind=IntKind) :: my_atom, iri, irj, kli, klj, is, js, ns
   integer (kind=IntKind) :: lmaxi, lmaxj, kkri, kkrj, kkri_ns, kkrj_ns
   integer (kind=IntKind) :: nbi, nbj, sbi, sbj, i, j, kl, jl, my_atom2, nindex
   integer (kind=IntKind) :: nbr_kkri, nbr_kkri_ns, nbr_kkrj, nbr_kkrj_ns
   integer (kind=IntKind) :: tsize, kmax_max_ns, dsize, kkrsz_ns
!
   integer (kind=IntKind) :: lid_i, lid_j, pei, pej, jid, info, iid
!
   real (kind=RealKind) :: rij(3), posi(3), posj(3), t0
!
   complex (kind=CmplxKind) :: kappa_1, kappa
   complex (kind=CmplxKind), pointer :: p_jinvi(:)
   complex (kind=CmplxKind), pointer :: p_sinej(:)
   complex (kind=CmplxKind), pointer :: pBigMatrix(:,:)
   complex (kind=CmplxKind), pointer :: pBlockMatrix(:,:)
   complex (kind=CmplxKind), pointer :: OmegaHat(:,:), OmegaHatj(:,:), p_oh(:)
   complex (kind=CmplxKind), pointer :: wau_g(:,:)
   complex (kind=CmplxKind), pointer :: wau_l(:,:,:)
   complex (kind=CmplxKind), pointer :: jig(:), ubmat(:,:)
   complex (kind=CmplxKind), pointer :: sm1(:,:), sm2(:,:), tm(:,:)
   complex (kind=CmplxKind), pointer :: smi(:,:), smj(:,:)
   complex (kind=CmplxKind), pointer :: kau_l(:,:), tau_l(:,:)
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
   interface
      subroutine convertGijToRel(gij, bgij, kkr1, kkr2, ce)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: kkr1, kkr2
         complex (kind=CmplxKind), intent(in) :: gij(:,:)
         complex (kind=CmplxKind), intent(out) :: bgij(:,:)
         complex (kind=CmplxKind), intent(in) :: ce
      end subroutine convertGijToRel
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('calTauMatrix','Module not initialized')
   endif
!
!  -------------------------------------------------------------------
   call exchangeSSSMatrix(getSingleScatteringMatrix)
!  -------------------------------------------------------------------
!
   if (present(kp)) then
     kappa_1 = kp
   else
     kappa_1 = sqrt(energy)
   endif
   kmax_max_ns = kmax_max*n_spin_cant
   tsize = kmax_max*kmax_max*n_spin_cant*n_spin_cant
!
   call nvtxStartRange("Step 1")
!
   do my_atom = 1, LocalNumAtoms
      Neighbor => getNeighbor(my_atom)
      kkrsz_ns = kmax_kkr(my_atom)*n_spin_cant
      jig => wsStoreA
!     ================================================================
!     Construct BigMatrix for local atom: my_atom
!     ================================================================
      dsize = dsize_l(my_atom)  ! This is the number of rows of BigMatrix
!     ================================================================
!     initialize BigMatrix so that it is a unit matrix to begin with..
!     ================================================================
      call setupUnitMatrix(dsize,BigMatrix)
!     ================================================================
!     loop over cluster neighbors to build BigMatrix [1-Jinv*G*S]
!
!     NOTE: The following code may fail if the 
!           number of neighbors associated with each atom is different.
!     ================================================================
      nbj=0
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE(rij,kappa_1,kappa,lmaxi,lmaxj,posj,p_sinej,p_jinvi,store_sine,posi,gij,gijrel) &
      !$ACC DEFAULT (PRESENT) &
      !$ACC REDUCTION(+:nbj)
      do irj=0,Neighbor%NumAtoms
!        =============================================================
!        get the Jinv-matrix and Sine-matrix and store the data in
!        the following way:
!             jinvi(nbr_kkrj_ns,nbr_kkrj_ns)
!             sinej(nbr_kkrj_ns,nbr_kkrj_ns)
!
!        Note: nbr_kkri and nbr_kkri are the actual kmax of the neighboring
!              atom, while kkri and kkrj are the reduced kmax of the neighboring 
!              atom which decreases as the distance from the center atom further
!              away.
!        =============================================================
         if (irj == 0) then
            lid_j = my_atom
            pej = MyPEinGroup
            posj(1:3) = Position(1:3,my_atom)
            lmaxj = lmax_kkr(my_atom)
         else
            lid_j = Neighbor%LocalIndex(irj)
            pej = Neighbor%ProcIndex(irj)
            posj(1:3) = Neighbor%Position(1:3,irj)
            lmaxj = Neighbor%Lmax(irj)
         endif
!
         if (pej == MyPEinGroup) then
            nbr_kkrj = kmax_kkr(lid_j)
!!!         p_sinej => local_jsmtx(tsize+2:2*tsize+1,lid_j)
         else
            jid = mapping_jsmtx(irj,my_atom)
            nbr_kkrj = remote_kmax(jid)
!!!         p_sinej => remote_jsmtx(tsize+1:2*tsize,jid)
         endif
         kkrj  = (lmaxj+1)*(lmaxj+1)
!        ==============================================================
!        nbr_kkrj = kmax of the remote atom
!        kkrj = kmax to be used in the construction of the "big matrix"
!             <= nbr_kkrj
!        ==============================================================
#ifndef _OPENACC
         if (nbr_kkrj < kkrj) then
!           -----------------------------------------------------------
            call ErrorHandler('calTauMatrix','remote atom kmax < kmaxj',nbr_kkrj,kkrj)
!           -----------------------------------------------------------
         endif
#endif
         kkrj_ns  = kkrj*n_spin_cant
         nbr_kkrj_ns  = nbr_kkrj*n_spin_cant
!
!        =============================================================
!        Rearrange the sine matrix storage if needed..................
         if (n_spin_cant == 2 .and. kkrj_ns /= nbr_kkrj_ns) then
            do js = 1, n_spin_cant
               sbi = (js-1)*nbr_kkrj_ns*nbr_kkrj
               sbj = (js-1)*kkrj_ns*kkrj
               do klj = 1, kkrj
                  do is = 1, n_spin_cant
                     do kli = 1, kkrj
                        store_sine(sbj+kli)=p_sinej(sbi+kli)
                     enddo
                     sbi = sbi + nbr_kkrj
                     sbj = sbj + kkrj
                  enddo
               enddo
            enddo
         else
            store_sine = CZERO
         endif
!        =============================================================
         nbi = 0
         do iri=0,Neighbor%NumAtoms
            if (iri == 0) then
               lid_i = my_atom
               pei = MyPEinGroup
               posi(1:3) = Position(1:3,my_atom)
               lmaxi = lmax_kkr(my_atom)
            else
               lid_i = Neighbor%LocalIndex(iri)
               pei = Neighbor%ProcIndex(iri)
               posi(1:3) = Neighbor%Position(1:3,iri)
               lmaxi = Neighbor%Lmax(iri)
            endif
!
            kkri = (lmaxi+1)*(lmaxi+1)
            if (pei == MyPEinGroup) then
               nbr_kkri = kmax_kkr(lid_i)
!!!            p_jinvi => local_jsmtx(2:tsize+1,lid_i)
            else
               iid = mapping_jsmtx(iri,my_atom)
               nbr_kkri = remote_kmax(iid)
!!!            p_jinvi => remote_jsmtx(1:tsize,iid)
            endif
!           ===========================================================
!           nbr_kkri = kmax of the remote atom
!           kkri = kmax to be used in the construction of the "big matrix"
!                <= nbr_kkri
!           ===========================================================
#ifndef _OPENACC
            if (nbr_kkri < kkri) then
!              --------------------------------------------------------
               call ErrorHandler('calTauMatrix','remote atom kmax < kmaxi',nbr_kkri,kkri)
!              --------------------------------------------------------
            endif
#endif
            kkri_ns = kkri*n_spin_cant
            nbr_kkri_ns = nbr_kkri*n_spin_cant
!
            if ( irj /= iri ) then
!              =======================================================
!              g(Rij) calculation
!              The data is stored as gij(kkri,kmax_phi_j)
!
!                       s,s"     ij'      s",s'
!              form Jinv  (e) * g  (e) * S   (e) for current Rij and
!                   ----i       -s"      -j'
!
!              store the result into BigMatrix => 1 - Jinv*G*S/kappa
!              =======================================================
!08/28/14      rij(1:3) = posi(1:3)-posj(1:3)
               rij = posj - posi
               kappa = kappa_1
               call getStrConstMatrix(lmaxi,lmaxj,kappa,rij,gij)
!
#ifndef _OPENACC
               if (Relativity > 1) then
!                 ----------------------------------------------------
                  call convertGijToRel(gij, gijrel, kkri, kkrj, energy)
!                 ----------------------------------------------------
                  call zgemm('n', 'n', kkri_ns, kkrj_ns, kkri_ns, CONE,&
                             p_jinvi, nbr_kkri_ns, gijrel, kkri_ns, CZERO, jig, kkri_ns)
!                 ----------------------------------------------------
               else
                  do js = 1, n_spin_cant
                     do is = 1, n_spin_cant
!                       ----------------------------------------------
                        call zgemm('n', 'n', kkri, kkrj, kkri, CONE,                      &
                                   p_jinvi((js-1)*nbr_kkri_ns*nbr_kkri+(is-1)*nbr_kkri+1),&
                                   nbr_kkri_ns, gij, kkri, CZERO,                         &
                                   jig((js-1)*kkri_ns*kkrj+(is-1)*kkri+1), kkri_ns)
!                       ----------------------------------------------
                     enddo
                  enddo
               endif
               if (n_spin_cant == 2 .and. kkrj_ns /= nbr_kkrj_ns) then
                  ! In this case, p_sinej data have been rearranged and copied
                  ! to store_sine
!                 ----------------------------------------------------
                  call zgemm('n', 'n', kkri_ns, kkrj_ns, kkrj_ns, -CONE/kappa_1, &
                             jig, kkri_ns, store_sine, kkrj_ns, CZERO,         &
                             BigMatrix(dsize*nbj+nbi+1), dsize)
!                 ----------------------------------------------------
               else
!                 ----------------------------------------------------
                  call zgemm('n', 'n', kkri_ns, kkrj_ns, kkrj_ns, -CONE/kappa_1, &
                             jig, kkri_ns, p_sinej, nbr_kkrj_ns, CZERO,        &
                             BigMatrix(dsize*nbj+nbi+1), dsize)
!                 ----------------------------------------------------
                 !call writeMatrix('BigMatrixIteration', BigMatrix, dsize*dsize)
               endif
#endif
            endif
            nbi = nbi + kkri_ns
         enddo
         nbj = nbj + kkrj_ns
      enddo
!     call writeMatrix('BigMatrix',BigMatrix,dsize_max*dsize_max)
      call nvtxEndRange()
      call nvtxStartRange("Step 2")
!
      if ( nbi /= nbj ) then
!        -------------------------------------------------------------
         call ErrorHandler(' calTauMatrix ::', 'nbi /= nbj', nbi, nbj)
!        -------------------------------------------------------------
      else if ( nbi /= dsize ) then
!        -------------------------------------------------------------
         call ErrorHandler(' calTauMatrix ::', 'nbi /= dsize', dsize)
!        -------------------------------------------------------------
      endif
!
      Tau00(my_atom)%kau_l = CZERO
      if (istauij_needed) then
        do j = 1, Neighbor%NumAtoms+1
          Tau00(my_atom)%neigh1j(j)%kau_l = CZERO
          Tau00(my_atom)%neigh1j(j)%wau_l = CZERO
          Tau00(my_atom)%neighj1(j)%kau_l = CZERO
          Tau00(my_atom)%neighj1(j)%wau_l = CZERO
        enddo
      endif
!
      if (Neighbor%NumAtoms > 0) then
!        =============================================================
!        invert [1 - jinv*gij*S/kappa] and obtain the kmax_kkr_ns*kmax_kkr_ns block only.
!        returned in pBlockMatrix is used to calculate wau_g as follows:
!        wau_g = [(1-pBlockMatrix)^{-1}-1]
!              = (1-pBlockMatrix)^{-1}*pBlockMatrix
!        =============================================================
         if (istauij_needed) then
           allocate(BigMatrixInv(dsize, dsize))
           BigMatrixInv = CZERO
           call zcopy(dsize*dsize, BigMatrix,1,BigMatrixInv,1) 
          !pBigMatrix => aliasArray2_c(BigMatrixInv,dsize,dsize)
           call MtxInv_LU(BigMatrixInv, dsize)
           do i = 1,dsize
             BigMatrixInv(i,i) = BigMatrixInv(i,i) - CONE
           enddo
           Neighbor => getNeighbor(my_atom)
           do j = 1, Neighbor%NumAtoms+1
              do kl = 1, kmax_kkr(my_atom)
                 do jl = 1, kmax_kkr(my_atom)
                    Tau00(my_atom)%neigh1j(j)%wau_l(jl,kl) = BigMatrixInv(jl, (j-1)*kkrsz_ns + kl)
                    Tau00(my_atom)%neighj1(j)%wau_l(jl,kl) = BigMatrixInv((j-1)*kkrsz_ns + jl, kl)
                 enddo
             enddo
          !  call zcopy(kkrsz_ns*kkrsz_ns, &
          !      pBigMatrix((j-1)*kkrsz_ns+1:j*kkrsz_ns,(j-1)*kkrsz_ns+1:j*kkrsz_ns), &
          !      1, Tau00(my_atom)%neighMat(j)%wau_l, 1)
           enddo
  !        call writeMatrix("BigMatrixInv00",Tau00(my_atom)%neighMat(1)%wau_l, kkrsz_ns, kkrsz_ns)
  !        call writeMatrix("BigMatrixInv10",pBigMatrix(kkrsz_ns+1:2*kkrsz_ns, 1:kkrsz_ns), kkrsz_ns, kkrsz_ns)
           deallocate(BigMatrixInv)
         endif
         BlockMatrix = CZERO
         pBlockMatrix => aliasArray2_c(BlockMatrix,kkrsz_ns,kkrsz_ns)
!        =============================================================
!        Note: For inversion algorithm <= 6, the following matrix inverse
!              routine appears problematic for lmax > 3 in the full-
!              potential case. The sympotom is: 
!              Different MPI processes give slightly different results.
!              It needs to be further checked
!        =============================================================
!write(6,'(a,i5,2f9.5,2x,4d16.8)')'BigMatrix =',MyPE,energy,          &
!      BigMatrix((dsize-1)*dsize-kkrsz_ns-1),BigMatrix(kkrsz_ns*dsize+2)
!call syncAllPEs()
!        -------------------------------------------------------------
!         print *,"In calClusterMatrix"
!         print *,"my_atom = ",my_atom
!         print *,"pBlockMatrix shape", shape(pBlockMatrix)
!         print *,"kkrsz_ns",kkrsz_ns
!         print *,"dsize = ",dsize
!         print *,"BigMatrix shape",shape(BigMatrix)
#ifdef ACCEL
         pBigMatrix => aliasArray2_c(BigMatrix, dsize, dsize)
         call invertMatrixLSMS_CUDA(my_atom, pBlockMatrix, kkrsz_ns,  &
                                 pBigMatrix, dsize )
!         print *,"CUDA pBlockMatrix(:, 1) = ",pBlockMatrix(:,1)
#else
         call invertMatrixBlock( my_atom, pBlockMatrix, kkrsz_ns, kkrsz_ns,  &
                                 BigMatrix, dsize, dsize )
#endif
!        -------------------------------------------------------------
!write(6,'(a,i5,2f9.5,2x,4d16.8)')'Block =',MyPE,energy,pBlockMatrix(1,1),pBlockMatrix(3,3)
!call syncAllPEs()
!
!        =============================================================
!        store [1 - BlockMatrix] in ubmat, which uses wsTau00L as temporary space
!        =============================================================
         ubmat => aliasArray2_c(wsTau00L, kkrsz_ns, kkrsz_ns)
!        -------------------------------------------------------------
         call setupUnitMatrix(kkrsz_ns,ubmat)
!        -------------------------------------------------------------
         ubmat = ubmat - pBlockMatrix
!
         wau_g => aliasArray2_c(wsStoreA, kkrsz_ns, kkrsz_ns)
!        =============================================================
!        Calculate wau_g = {[1-Jinv*G*S]**(-1)-1} : for central site only
!        solve the eigenvectors of the following linear equations:
!              [1 - pBlockMatrix]*wau_g = pBlockMatrix
!        to obtain wau_g
!        -------------------------------------------------------------
         call zcopy(kkrsz_ns*kkrsz_ns,pBlockMatrix,1,wau_g,1)
!        -------------------------------------------------------------
         call zgetrf(kkrsz_ns, kkrsz_ns, ubmat, kkrsz_ns, ipvt, info)
!        -------------------------------------------------------------
         call zgetrs( 'n', kkrsz_ns, kkrsz_ns, ubmat, kkrsz_ns, ipvt, &
                      wau_g, kkrsz_ns, info )
!        -------------------------------------------------------------
         nullify( pBlockMatrix )
         nullify( pBigMatrix )
!        -------------------------------------------------------------
!        call writeMatrix('wau_g', wau_g, kkrsz_ns, kkrsz_ns)
!        =============================================================
!        Rotate the wau_g-matrix from global frame to local frame of references
!        Note: if n_spin_cant = 1, wau_g and wau_l share the same space.
!        =============================================================
         if ( n_spin_cant==2 ) then
            wau_l => aliasArray3_c(wsTau00L, kmax_kkr(my_atom), kmax_kkr(my_atom), 4)
!           ----------------------------------------------------------
            call rotateGtoL(my_atom, kmax_kkr(my_atom), kmax_kkr(my_atom),&
                            wau_g, wau_l )
!           ----------------------------------------------------------
         else ! set wau_l = wau_g
            wau_l => aliasArray3_c(wsStoreA, kmax_kkr(my_atom), kmax_kkr(my_atom), 1)
         endif
!
!        =============================================================
!        Calculate kau_l = wau_l*OmegaHat*kappa in local space
!        =============================================================
         ns = 0
         do js = 1, n_spin_cant
!           ----------------------------------------------------------
            OmegaHat => getSingleScatteringMatrix('OmegaHat-Matrix',spin=js,site=my_atom)
!           ----------------------------------------------------------
            do is = 1, n_spin_cant
               ns = ns + 1
!              -------------------------------------------------------
               call zgemm('n', 'n', kmax_kkr(my_atom), kmax_kkr(my_atom), &
                          kmax_kkr(my_atom), kappa_1, wau_l(1,1,ns),        &
                          kmax_kkr(my_atom), OmegaHat, kmax_kkr(my_atom), &
                          CZERO, Tau00(my_atom)%kau_l(1,1,ns), kmax_kkr(my_atom))
!              -------------------------------------------------------
               if (istauij_needed) then ! Needs to be re-examined if the remote atom has a different kkr size
                 do j = 1, Neighbor%NumAtoms+1
                   if (j == 1) then
!                    --------------------------------------------------------
                     OmegaHatj => getSingleScatteringMatrix('OmegaHat-Matrix',spin=js,site=my_atom)
!                    --------------------------------------------------------
                   else if (Neighbor%ProcIndex(j-1) == MyPEinGroup) then
                     lid_j = Neighbor%LocalIndex(j-1)
!                    --------------------------------------------------------
                     OmegaHatj => getSingleScatteringMatrix('OmegaHat-Matrix',spin=js,site=lid_j)
!                    --------------------------------------------------------
                   else
                     jid = mapping_jsmtx(j-1,my_atom)
                     p_oh => remote_jsmtx(2*tsize+1:3*tsize,jid)
                     OmegaHatj => aliasArray2_c(p_oh,kmax_kkr(my_atom),kmax_kkr(my_atom))
                   endif
                   call zgemm('n', 'n', kmax_kkr(my_atom), kmax_kkr(my_atom), &
                    kmax_kkr(my_atom), kappa_1, Tau00(my_atom)%neigh1j(j)%wau_l, &
                    kmax_kkr(my_atom), OmegaHatj, kmax_kkr(my_atom), &
                    CZERO, Tau00(my_atom)%neigh1j(j)%kau_l, kmax_kkr(my_atom))
!                  ----------------------------------------------------------
                   call zgemm('n', 'n', kmax_kkr(my_atom), kmax_kkr(my_atom), &
                    kmax_kkr(my_atom), kappa_1, Tau00(my_atom)%neighj1(j)%wau_l, &
                    kmax_kkr(my_atom), OmegaHat, kmax_kkr(my_atom), &
                    CZERO, Tau00(my_atom)%neighj1(j)%kau_l, kmax_kkr(my_atom))
!                  ----------------------------------------------------------
                 enddo
               endif
!              =======================================================
!              As a check, set Tau00(my_atom)%kau_l = kappa*OmegaHat
!              the resulting Green function should be the single site 
!              Green function
!              -------------------------------------------------------
!              call zcopy(kmax_kkr(my_atom)*kmax_kkr(my_atom),OmegaHat,1, &
!                         Tau00(my_atom)%kau_l(1,1,ns),1)
!              Tau00(my_atom)%kau_l = kappa*Tau00(my_atom)%kau_l
!              =======================================================
            enddo
         enddo
      endif
   enddo
   call nvtxEndRange()
   call nvtxStartRange("Step 3")
!  ===================================================================
!  Determine Tau00 = S*Kau00*S^{*T}/energy + t_matrix
!  calculate tau00 if needed
!  ===================================================================
   if (present(tau_needed)) then
      if (tau_needed) then
         do my_atom = 1, LocalNumAtoms
            ns = 0
            do js = 1, n_spin_cant
               sm2 => getSingleScatteringMatrix('Sine-Matrix',spin=js,site=my_atom)
               tm => getSingleScatteringMatrix('T-Matrix',spin=js,site=my_atom)
               do is = 1, n_spin_cant
                  sm1 => getSingleScatteringMatrix('Sine-Matrix',spin=is,site=my_atom)
                  ns = ns + 1
                  kau_l => Tau00(my_atom)%kau_l(:,:,ns)
                  tau_l => Tau00(my_atom)%tau_l(:,:,ns)
!                 ----------------------------------------------------
                  call computeUAUts(sm1,kmax_kkr(my_atom),kmax_kkr(my_atom), &
                                    sm2,kmax_kkr(my_atom), CONE/energy,      &
                                    kau_l,kmax_kkr(my_atom),CZERO,           &
                                    tau_l,kmax_kkr(my_atom),wsStoreA)
!                 ----------------------------------------------------
                  if (is == js) then
                     tau_l = tau_l + tm
                  endif
               enddo
            enddo
         enddo
      endif
   endif

   if (istauij_needed) then ! Needs to be re-examined in the spin-canted case
     do my_atom = 1, LocalNumAtoms
       Neighbor => getNeighbor(my_atom)
       smi => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=my_atom)
       tm => getSingleScatteringMatrix('T-Matrix',spin=1,site=my_atom)
       do j = 1, Neighbor%NumAtoms+1
         if (j == 1) then
           smj => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=my_atom)
         else if (Neighbor%ProcIndex(j-1) == MyPEinGroup) then
           lid_j = Neighbor%LocalIndex(j-1)
           smj => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=lid_j)
         else
            jid = mapping_jsmtx(j-1,my_atom)
            p_sinej => remote_jsmtx(tsize+1:2*tsize,jid)
            smj => aliasArray2_c(p_sinej,kmax_kkr(my_atom),kmax_kkr(my_atom))
         endif
         wsStoreA = CZERO
         call computeUAUts(smi,kmax_kkr(my_atom),kmax_kkr(my_atom), &
                          smj,kmax_kkr(my_atom), CONE/energy, &
                Tau00(my_atom)%neigh1j(j)%kau_l, kmax_kkr(my_atom),CZERO, &
                Tau00(my_atom)%neigh1j(j)%tau_l,kmax_kkr(my_atom),wsStoreA)
         wsStoreA = CZERO
         call computeUAUts(smj,kmax_kkr(my_atom),kmax_kkr(my_atom), &
                          smi,kmax_kkr(my_atom), CONE/energy, &
                Tau00(my_atom)%neighj1(j)%kau_l, kmax_kkr(my_atom),CZERO, &
                Tau00(my_atom)%neighj1(j)%tau_l,kmax_kkr(my_atom),wsStoreA)
       enddo
       Tau00(my_atom)%neigh1j(1)%tau_l = Tau00(my_atom)%neigh1j(1)%tau_l + tm
       Tau00(my_atom)%neighj1(1)%tau_l = Tau00(my_atom)%neighj1(1)%tau_l + tm
     enddo
   endif
!
   call nvtxEndRange()

   nullify(p_jinvi, p_sinej, jig, wau_g, wau_l, ubmat, pBlockMatrix, pBigMatrix)
!
   end subroutine calClusterMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTau(local_id,global_i,global_j,dsize) result(ptau)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: local_id
   integer (kind=IntKind), intent(in), optional :: global_i, global_j
   integer (kind=IntKind), intent(out), optional :: dsize
!
   integer (kind=IntKind) :: n, k
   complex (kind=CmplxKind), pointer :: ptau(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getTau','Module not initialized')
   else if ( local_id < 1 .or. local_id > LocalNumAtoms) then
      call ErrorHandler('getTau','invalid local atom index',local_id)
   endif
!
   if (present(dsize)) then
      dsize = Tau00(local_id)%isize
   endif
!
   if (.not.present(global_i) .and. .not.present(global_j)) then
      ptau => Tau00(local_id)%tau_l
      return
   else if (.not.present(global_j)) then
      call ErrorHandler('getTau','both i and j indeces need to be specified')
   else if (global_i /= 1 .and. global_j /= 1) then
      call ErrorHandler('getTau','One of i or j needs to be 1',global_i,global_j)
   endif
!
   Neighbor => getNeighbor(local_id)
!
   if (global_j == 1 .and. global_i == 1) then
      call ErrorHandler('getTau','It is needs to be implemented for calculating tau00')
      ptau => Tau00(local_id)%tau_l
      return
   else if (global_i < 1 .or. global_i > Neighbor%NumAtoms) then
      call ErrorHandler('getTau','i index is out of range',global_i)
   else if (global_j < 1 .or. global_j > Neighbor%NumAtoms) then
      call ErrorHandler('getTau','j index is out of range',global_j)
   else
      call ErrorHandler('getTau','Only Tau00 is implemented.',global_i,global_j)
   endif
!
   end function getTau
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function checkIfNeighbor(id1, id2) result(is_neighbor)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: id1, id2
   real (kind=RealKind) :: R12
   real (kind=RealKind) :: pos1(3), pos2(3) 
   logical :: is_neighbor

   Neighbor => getNeighbor(id1)
   pos1 = Position(1:3,id1)
   pos2 = Position(1:3,id2)
   R12 = (pos1(1) - pos2(1))**2 + (pos1(2) - pos2(2))**2 + (pos1(3) - pos2(3))**2
   if (sqrt(R12) <= maxval(Neighbor%ShellRad)) then
     is_neighbor = .true.
   else
     is_neighbor = .false.
   endif 

   end function checkIfNeighbor
!  ===================================================================
!
!  *******************************************************************
!  
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function determineNeighborIndex(id1, id2) result(nindex)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: id1, id2
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: nindex
   real (kind=RealKind) :: pos2(3), posn(3)

   pos2 = Position(1:3, id2)
   Neighbor => getNeighbor(id1)
   do i = 1, Neighbor%NumAtoms
     posn = Neighbor%Position(1:3, i)
     if (posn(1) == pos2(1) .and. posn(2) == pos2(2) &
                          .and. posn(3) == pos2(3)) then
       nindex = i+1
     endif
   enddo

   if (id1 == id2) then
     nindex = 1
   endif

   end function determineNeighborIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNeighborTau(id1, nindex, order) result(ptau)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: id1, nindex, order
   complex (kind=CmplxKind) :: ptau(kmax_max,kmax_max)

!  print *, "Tau for (m,n) =", id2, id1

!  nindex = determineNeighborIndex(id2, id1)
   if (order == 0) then
     ptau = Tau00(id1)%neigh1j(nindex)%tau_l
   else
     ptau = Tau00(id1)%neighj1(nindex)%tau_l
   endif
!  call writeMatrix("Taustored", Tau00(id1)%neighMat(nindex)%tau_l, kmax_max, kmax_max)
!  call writeMatrix("ptau", ptau, kmax_max, kmax_max)

   end function getNeighborTau
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getClusterTau(id1, id2) result(ptau)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: id1, id2
   complex (kind=CmplxKind) :: ptau(kmax_max,kmax_max)

   ptau = ClusterTau((id1-1)*kmax_max+1:id1*kmax_max, &
                     (id2-1)*kmax_max+1:id2*kmax_max)

   end function getClusterTau
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKau(local_id,global_i,global_j,dsize) result(pkau)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: local_id
   integer (kind=IntKind), intent(in), optional :: global_i, global_j
   integer (kind=IntKind), intent(out), optional :: dsize
!
   integer (kind=IntKind) :: n, k
   complex (kind=CmplxKind), pointer :: pkau(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getKau','Module not initialized')
   else if ( local_id < 1 .or. local_id > LocalNumAtoms) then
      call ErrorHandler('getKau','invalid local atom index',local_id)
   endif
!
   if (present(dsize)) then
      dsize = Tau00(local_id)%isize
   endif
!
   if (.not.present(global_i) .and. .not.present(global_j)) then
      pkau => Tau00(local_id)%kau_l
      return
   else if (.not.present(global_j)) then
      call ErrorHandler('getKau','both i and j indeces need to be specified')
   else if (global_i /= 1 .and. global_j /= 1) then
      call ErrorHandler('getKau','Not implemented for i<>1 and/or j<>1.', &
                        global_i,global_j)
   endif
!
   pkau => Tau00(local_id)%kau_l
!
   end function getKau
!  ====================================================================
end module ClusterMatrixModule
