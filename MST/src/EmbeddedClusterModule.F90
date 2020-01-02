module EmbeddedClusterModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1
   use MathParamModule, only : two, half, third
   use MathParamModule, only : ten2m8, ten2m6, ten2m14, ten2m9, ten2m10
   use MathParamModule, only : sqrtm1, pi2, pi
   use TimerModule, only: getTime
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : MatrixBandStruct
!
public :: initEmbeddedCluster,     &
          endEmbeddedCluster,      &
          setupHostMedium,         &
          beginEmbedding,          &
          endEmbedding,            &
          embedScatterInHostMedium,&
          calEmbeddedSiteMatrix,   &
          calEmbeddedClusterMatrix,&
          getTau
!
private
   type (MatrixBandStruct), allocatable :: EmbeddedClusterMatrixBand(:)  ! Matrix is stored a band of columns
!
   type ScatterInsertStruct
      integer (kind=IntKind) :: embed_index
      integer (kind=IntKind) :: local_index
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: block_index
      integer (kind=IntKind) :: kmax_kkr
      complex (kind=CmplxKind), pointer :: diff_TmatInv(:,:)
      type (ScatterInsertStruct), pointer :: next
      type (ScatterInsertStruct), pointer :: prev
   end type ScatterInsertStruct
!
   type (ScatterInsertStruct), target :: EmbeddedCluster
   type (ScatterInsertStruct), pointer :: present_embed
!
   integer (kind=IntKind) :: NumEmbeddedAtoms
!
   complex (kind=CmplxKind), allocatable, target :: Tau_MatrixBand(:,:,:)
   complex (kind=CmplxKind), allocatable, target :: TMP_MatrixBand(:)
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: GlobalNumSitesInCluster, LocalNumSitesInCluster
   integer (kind=IntKind) :: nSpinCant
!
   integer (kind=IntKind), allocatable :: gid_array(:)! Relates a global index to the corresponding matrix block 
                                                      ! row index
!
   complex (kind=CmplxKInd), allocatable, target :: WORK(:)
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
!
   integer (kind=IntKind) :: host_kmax_kkr
!
   complex (kind=CmplxKind), allocatable, target :: t_host_inv(:,:) ! t-matrix inverse, in local spin frame, of the 
                                                                    ! host sites on local CPU
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initEmbeddedCluster(cant, istop, iprint)
!  ===================================================================
   use MPPModule, only : MyPE
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : syncAllPEsInGroup, isGroupExisting
   use MediumHostModule, only : getLmaxKKR, getNumSites, getLocalNumSites
   use MediumHostModule, only : getGlobalSiteIndex
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: cant    ! Spin canting index
   integer (kind=IntKind), intent(in) :: iprint
!
   character (len=14) :: sname = "initEmbeddedCluster"
!
   integer (kind=IntKind) :: i, j, ig, na, n, k, nk, lmax_kkr
   integer (kind=IntKind) :: kmax_kkr, kmax_kkr_ns, tsize
   integer (kind=IntKind) :: KKRMatrixSize
   integer (kind=IntKind) :: KKRMatrixSizeCant
   integer (kind=IntKind) :: BandSize
   integer (kind=IntKind) :: BandSizeCant
!
   complex (kind=CmplxKind), pointer :: p1(:)
!
   type (ScatterInsertStruct), pointer :: head_embed
!
   stop_routine = istop
!
   GlobalNumSitesInCluster = getNumSites()      ! This is the total number of sublattices
   LocalNumSitesInCluster = getLocalNumSites()  ! This is the number of sublattices on the local process
   nSpinCant = cant
   print_level= iprint
!
   if (nSpinCant > 1) then
!     ================================================================
!     In spin-canted case, the 3d dimension of the tau matrix is greater
!     than 1. Some code charges are needed.
!     ----------------------------------------------------------------
      call ErrorHandler('initEmbeddedCluster',                        &
                        'For spin-canted case, the code needs some work')
!     ----------------------------------------------------------------
   endif
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
!  ===================================================================
!
   allocate( EmbeddedClusterMatrixBand(LocalNumSitesInCluster) )
   allocate( gid_array(GlobalNumSitesInCluster) )
!
   BandSize = 0
   do j = 1, LocalNumSitesInCluster
      ig = getGlobalSiteIndex(j)
      EmbeddedClusterMatrixBand(j)%global_index = ig
      lmax_kkr = getLmaxKKR(ig)
      kmax_kkr = (lmax_kkr+1)**2
      kmax_kkr_ns = kmax_kkr*nSpinCant
      BandSize = BandSize + kmax_kkr
      EmbeddedClusterMatrixBand(j)%lmax_kkr = lmax_kkr
      EmbeddedClusterMatrixBand(j)%kmax_kkr = kmax_kkr
      allocate( EmbeddedClusterMatrixBand(j)%MatrixBlock(GlobalNumSitesInCluster) )
      EmbeddedClusterMatrixBand(j)%column_index = (j-1)*kmax_kkr_ns + 1
   enddo
!
   BandSizeCant = BandSize*nSpinCant
   KKRMatrixSize = 0
   tsize = 0
   do ig = 1, GlobalNumSitesInCluster
      lmax_kkr = getLmaxKKR(ig)
      kmax_kkr = (lmax_kkr+1)**2
      kmax_kkr_ns = kmax_kkr*nSpinCant
      KKRMatrixSize = KKRMatrixSize + kmax_kkr
      tsize = max(tsize,kmax_kkr_ns*kmax_kkr_ns)
   enddo
   KKRMatrixSizeCant = KKRMatrixSize*nSpinCant
!
   allocate ( Tau_MatrixBand(tsize,GlobalNumSitesInCluster,LocalNumSitesInCluster) )
   allocate ( TMP_MatrixBand(KKRMatrixSizeCant*BandSizeCant) )
   Tau_MatrixBand = CZERO; TMP_MatrixBand = CZERO
!
   allocate( WORK(tsize), t_host_inv(tsize,LocalNumSitesInCluster) )
   WORK = CZERO; t_host_inv = CZERO
!
!  ===================================================================
!  Set up the matrix block row/column index in such way that each
!  a band of columns of the matrix is stored on a processor
!      n = the matrix block row/column index in the band matrix
!     ig = the global index of the site that corresponds to the matrix block
!  ===================================================================
   n = 0
   nk = 0
   do k = 1, NumPEsInGroup
      na = getLocalNumSites(k-1)
      do i = 1, na
         n = n + 1
         ig = getGlobalSiteIndex(i,k-1)
         lmax_kkr = getLmaxKKR(ig)
         kmax_kkr = (lmax_kkr+1)**2
         kmax_kkr_ns = kmax_kkr*nSpinCant
         do j = 1, LocalNumSitesInCluster
            EmbeddedClusterMatrixBand(j)%MatrixBlock(n)%lmax_kkr = lmax_kkr
            EmbeddedClusterMatrixBand(j)%MatrixBlock(n)%kmax_kkr = kmax_kkr
            EmbeddedClusterMatrixBand(j)%MatrixBlock(n)%global_index = ig
            EmbeddedClusterMatrixBand(j)%MatrixBlock(n)%row_index = nk + 1
!           ==========================================================
!           Note: kau_l stores tau_a, the tau-matrix for the atom in a
!                 cluster embedded in effective medium
!                 tau_l stores tau_c, the tau-matrix for the effective
!                 medium.
!           ==========================================================
            p1 => Tau_MatrixBand(:,ig,j)
            EmbeddedClusterMatrixBand(j)%MatrixBlock(n)%kau_l =>      &
                   aliasArray3_c(p1,kmax_kkr,kmax_kkr,nSpinCant*nSpinCant)
            nullify(EmbeddedClusterMatrixBand(j)%MatrixBlock(n)%tau_l)
         enddo
         nk = nk + kmax_kkr_ns
         gid_array(ig) = n
      enddo
   enddo
!
   host_kmax_kkr = 0
!
!  ===================================================================
!  Set up the structure of embedded cluster linked list
!  ===================================================================
   head_embed => EmbeddedCluster
   allocate(head_embed%diff_TmatInv(kmax_kkr_ns,kmax_kkr_ns))
   nullify(head_embed%prev)
   nullify(head_embed%next)
   do j = 2, LocalNumSitesInCluster
      allocate(head_embed%next)
      head_embed%next%prev => head_embed
      head_embed => head_embed%next
      kmax_kkr_ns = EmbeddedClusterMatrixBand(j)%kmax_kkr*nSpinCant
      allocate(head_embed%diff_TmatInv(kmax_kkr_ns,kmax_kkr_ns))
      nullify(head_embed%next)
   enddo
   nullify(present_embed)
!
!  ===================================================================
   if ( print_level >0 .or. MyPE==0 ) then
      write(6,*) "===================================================="
      write(6,*) "init EmbeddedCluster Module ::  "
      write(6,'(a,i5)') "   LocalNumSites   : ", LocalNumSitesInCluster
      write(6,'(a,i5)') "   n_spin_cant     : ", nSpinCant
      write(6,'(a,i5)') "   KKR_row_dim     : ", KKRMatrixSize
      write(6,'(a,i5)') "   KKR_col_dim     : ", BandSize
      do j = 1,LocalNumSitesInCluster
         write(6,'(a,i5)') "   Local Site Index : ", j
         write(6,'(a,i5)') "       LocalCol_Ind : ", EmbeddedClusterMatrixBand(j)%column_index
         write(6,'(a,i5)') "       kmax         : ", EmbeddedClusterMatrixBand(j)%kmax_kkr
         write(6,'(a,i5)') "       lmax         : ", EmbeddedClusterMatrixBand(j)%lmax_kkr
         write(6,'(a,i5)') "       global index : ", EmbeddedClusterMatrixBand(j)%global_index
         if (GlobalNumSitesInCluster < 101) then
            do i = 1, GlobalNumSitesInCluster
               write(6,'(a,4i5)')"   block index, starting row, kmax, global Index : ",    &
                                 i, EmbeddedClusterMatrixBand(j)%MatrixBlock(i)%row_index, &
                                 EmbeddedClusterMatrixBand(j)%MatrixBlock(i)%kmax_kkr,     &
                                 EmbeddedClusterMatrixBand(j)%MatrixBlock(i)%global_index
            enddo
         endif
      enddo
      write(6,*) "===================================================="
   endif
!
   end subroutine initEmbeddedCluster
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endEmbeddedCluster()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, j
!
   type (ScatterInsertStruct), pointer :: head_embed
!
   do j = 1, LocalNumSitesInCluster
      do i = 1, GlobalNumSitesInCluster
         nullify( EmbeddedClusterMatrixBand(j)%MatrixBlock(i)%tau_l )
         nullify( EmbeddedClusterMatrixBand(j)%MatrixBlock(i)%kau_l )
      enddo
      deallocate( EmbeddedClusterMatrixBand(j)%MatrixBlock )
   enddo
   deallocate( EmbeddedClusterMatrixBand )
!
   head_embed => EmbeddedCluster
   do i = 1, NumEmbeddedAtoms-1
      deallocate(head_embed%diff_TmatInv)
      head_embed => head_embed%next
   enddo
   deallocate(head_embed%diff_TmatInv)
   do i = NumEmbeddedAtoms, 2, -1
      head_embed => head_embed%prev
      deallocate(head_embed%next)
      nullify(head_embed%next)
   enddo
   nullify(present_embed)
!
   deallocate( Tau_MatrixBand )
   deallocate( TMP_MatrixBand )
   deallocate( gid_array )
!
   deallocate( WORK, t_host_inv )
!
   host_kmax_kkr = 0
!
   end subroutine endEmbeddedCluster
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTau(local_id,global_i,global_j,dsize) result(tau)
!  ===================================================================
   implicit none
!
!  *******************************************************************
!  This function returns tau-matrix for the atoms, with local_id, in 
!  the cluster embedded in the effective.
!  *******************************************************************
!
   integer (kind=IntKind), intent(in), optional :: local_id ! Local index
   integer (kind=IntKind), intent(in), optional :: global_i, global_j
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind) :: ni, jd, ig, id
!
   complex (kind=CmplxKind), pointer :: tau(:,:,:)
!
   if (present(local_id)) then
      if (local_id < 1 .or. local_id > LocalNumSitesInCluster) then
         call ErrorHandler('getTau','invalid local index i',local_id)
      else if (present(global_j)) then
         call ErrorHandler('getTau','The function call is confusing')
      endif
!
      jd = local_id
      if (present(global_i)) then
         ig = global_i
      else
!        take the diagonal block associated with the atom with local_id
         ig = EmbeddedClusterMatrixBand(jd)%global_index
      endif
   else if (present(global_i) .and. present(global_j)) then
      if (global_i < 1 .or. global_i > GlobalNumSitesInCluster) then
         call ErrorHandler('getTau','invalid global index i',global_i)
      else if (global_j < 1 .or. global_j > GlobalNumSitesInCluster) then
         call ErrorHandler('getTau','invalid global index j',global_j)
      endif
!     take the appropriate matrix block associated with the atom specified
!     by the global index pair: (global_i, global_j)
      jd = 0
      LOOP_id: do id = 1, LocalNumSitesInCluster
         if (EmbeddedClusterMatrixBand(id)%global_index == global_j) then
            jd = id
            exit LOOP_id
         endif
      enddo LOOP_id
      if (jd == 0) then
         call ErrorHandler('getTau',                                  &
              'The global j index is not associated with a local atom',global_j)
      endif
      ig = global_i
   else if (present(global_i) .or. present(global_j)) then
      call ErrorHandler('getTau','The function call is confusing')
   else
!     take the upper-left corner block
      jd = 1
      ig = EmbeddedClusterMatrixBand(jd)%global_index
   endif
!
   ni = gid_array(ig)
!  ===================================================================
!  Note: kau_l stores tau_a, the tau-matrix for the atom in a cluster
!        embedded in the effective medium
!  tau is in the local spin framework.
!  ===================================================================
   tau => EmbeddedClusterMatrixBand(jd)%MatrixBlock(ni)%kau_l
!
   if (present(dsize)) then
      dsize = EmbeddedClusterMatrixBand(jd)%MatrixBlock(ni)%kmax_kkr*nSpinCant
   endif
!
   end function getTau
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printEmbeddedClusterMatrix(local_i,global_j,e)
!  ===================================================================
   use WriteMatrixModule, only  : writeMatrix
   implicit none
!
   integer (kind=IntKind), intent(in) :: local_i ! Local index
   integer (kind=IntKind), intent(in), optional :: global_j ! Global index
   integer (kind=IntKind) :: ni, ig, id, jd, js
   integer (kind=IntKind) :: kmaxi, kmaxj, kmaxi_ns, kmaxj_ns
!
   complex (kind=CmplxKind), intent(in), optional :: e
   complex (kind=CmplxKind), pointer :: tau(:,:,:)
!
   if (local_i < 1 .or. local_i > LocalNumSitesInCluster) then
      call ErrorHandler('printEmbeddedClusterMatrix',                 &
                        'invalid local index i',local_i)
   else if (present(global_j)) then
      if (global_j < 1 .or. global_j > GlobalNumSitesInCluster) then
         call ErrorHandler('printEmbeddedClusterMatrix',              &
                           'invalid global index j',global_j)
      else
         jd = global_j
      endif
   else
      jd = EmbeddedClusterMatrixBand(local_i)%global_index ! For printing the diagonal block
   endif
!
   ni = gid_array(jd)
   kmaxj = EmbeddedClusterMatrixBand(local_i)%kmax_kkr
   kmaxj_ns = kmaxj*nSpinCant
   tau => EmbeddedClusterMatrixBand(local_i)%MatrixBlock(ni)%kau_l
!
   write(6,'(/,a)')'************************************************'
   write(6,'( a )')'*    Output from printEmbeddedClusterMatrix    *'
   write(6,'(a,/)')'************************************************'
   write(6,'(80(''=''))')
!  
   if (present(e)) then
      write(6,*)'energy ::',e
   endif
   write(6,'(80(''=''))')
!
   write(6,'(/,a,i4,2x,i4,/)')"    Sites i :: ",EmbeddedClusterMatrixBand(local_i)%global_index
   call writeMatrix('Tau-matrix in Global Frame of Reference in Spin Space ::', &
                     tau,kmaxj,kmaxj,nSpinCant*nSpinCant)
!
   nullify(tau)
!
   end subroutine printEmbeddedClusterMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  subroutine setupHostMedium(t_host,t_inverse,getMediumTau)
   subroutine setupHostMedium(e,getSingleSiteTmat,configuration)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
   use CrystalMatrixModule, only : calCrystalMatrix, getMediumTau=>getTau
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: configuration(:)
!
   complex (kind=CmplxKind), intent(in) :: e
!
   complex (kind=CmplxKind), pointer :: pt(:)
   complex (kind=CmplxKind), pointer :: t_host(:,:)
   complex (kind=CmplxKind), pointer :: p_thinv(:,:)
!
   integer (kind=IntKind) :: i, j, ig, jg
!
   interface
      function getSingleSiteTmat(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleSiteTmat
   end interface
!
!  ===================================================================
!  Setup t_host_inv. This needs to be checked for the spin canting case
!  ===================================================================
   do i = 1, LocalNumSitesInCluster
      t_host => getSingleSiteTmat('TInv-Matrix',spin=1,&
                                        site=i,atom=configuration(i))
      host_kmax_kkr = size(t_host,1)
      pt => t_host_inv(:,i)
      p_thinv => aliasArray2_c(pt,host_kmax_kkr,host_kmax_kkr)
      p_thinv = t_host
   enddo
!
!  -------------------------------------------------------------------
   call calCrystalMatrix(e,getSingleSiteTmat,use_tmat=.true.,         &
                         tau_needed=.true.,configuration=configuration)
!  -------------------------------------------------------------------
!
   if (GlobalNumSitesInCluster == 1) then
      EmbeddedClusterMatrixBand(1)%MatrixBlock(1)%tau_l => getMediumTau()
      EmbeddedClusterMatrixBand(1)%MatrixBlock(1)%kau_l = CZERO
   else
      do j = 1, LocalNumSitesInCluster
         jg = EmbeddedClusterMatrixBand(j)%global_index
         do i = 1, GlobalNumSitesInCluster
            ig = EmbeddedClusterMatrixBand(j)%MatrixBlock(i)%global_index
            EmbeddedClusterMatrixBand(j)%MatrixBlock(i)%tau_l =>      &
                                         getMediumTau(global_i=ig,global_j=jg)
            EmbeddedClusterMatrixBand(j)%MatrixBlock(i)%kau_l = CZERO
         enddo
      enddo
   endif
!
   end subroutine setupHostMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine beginEmbedding()
!  ===================================================================
   implicit none
!
   NumEmbeddedAtoms = 0
!
   end subroutine beginEmbedding
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine embedScatterInHostMedium(id,tmat_a,isInverse)
!  ===================================================================
!
!  place a species described by tmat at site "id" of the medium
!
!  *******************************************************************
   use MatrixInverseModule, only : MtxInv_LU
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i, ni, n, ig, dsize
!
   logical, intent(in) :: isInverse
!
   complex (kind=CmplxKind), intent(in) :: tmat_a(:,:)
   complex (kind=CmplxKind), pointer :: p_thinv(:,:)
   complex (kind=CmplxKind), pointer :: p_tainv(:,:)
   complex (kind=CmplxKind), pointer :: tau_c(:,:)
   complex (kind=CmplxKind), pointer :: pt(:)
!
   dsize = size(tmat_a,1)
   if (dsize /= host_kmax_kkr) then
!     ----------------------------------------------------------------
      call ErrorHandler('embedScatterInHostMedium',                   &
           'embedded atom t-matrix has different size than the host', &
           dsize,host_kmax_kkr)
!     ----------------------------------------------------------------
   else if (LocalNumSitesInCluster == NumEmbeddedAtoms) then
!     ----------------------------------------------------------------
      call ErrorHandler('embedScatterInHostMedium',                   &
                        'Number of embedded atoms exceeds the cluster size')
!     ----------------------------------------------------------------
   endif
!
   if (NumEmbeddedAtoms == 0) then
      present_embed => EmbeddedCluster
   else if (.not.associated(present_embed)) then
!     ----------------------------------------------------------------
      call ErrorHandler('embedScatterInHostMedium',                   &
                        'Embedding procedure is broken.')
!     ----------------------------------------------------------------
   else
      present_embed => present_embed%next
   endif
!
   NumEmbeddedAtoms = NumEmbeddedAtoms + 1
!
   present_embed%kmax_kkr = dsize/nSpinCant
   present_embed%embed_index = NumEmbeddedAtoms
   present_embed%local_index = id
   ig = EmbeddedClusterMatrixBand(id)%global_index
   present_embed%global_index = ig
   present_embed%block_index = gid_array(ig)
!
!  -------------------------------------------------------------------
   pt => t_host_inv(:,id)
   p_thinv => aliasArray2_c(pt,dsize,dsize)
!  -------------------------------------------------------------------
   if (isInverse) then ! tmat_a is the inverse of t-matrix
      present_embed%diff_TmatInv = tmat_a - p_thinv
   else
      p_tainv => aliasArray2_c(WORK,dsize,dsize)
      p_tainv = tmat_a
!     ----------------------------------------------------------------
      call MtxInv_LU(p_tainv,dsize)
!     ----------------------------------------------------------------
      present_embed%diff_TmatInv = p_tainv - p_thinv
   endif
!
   nullify(p_tainv, p_thinv)
!
   end subroutine embedScatterInHostMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endEmbedding()
!  ===================================================================
   implicit none
!
   if (LocalNumSitesInCluster < NumEmbeddedAtoms) then
!     ----------------------------------------------------------------
      call ErrorHandler('endEmbedding','inconsistent linked data',    &
                        LocalNumSitesInCluster,NumEmbeddedAtoms)
!     ----------------------------------------------------------------
   endif
!
   nullify(present_embed)
!
   end subroutine endEmbedding
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calEmbeddedSiteMatrix(id,mat,compute_x)
!  ===================================================================
   use MatrixModule, only : computeAprojB
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i, dsize, jd, ni
!
   logical, intent(in), optional :: compute_x
   logical :: found, Xmat
!
   complex (kind=CmplxKind), intent(out) :: mat(:,:)
   complex (kind=CmplxKind), pointer :: m_diff(:,:), tau_c(:,:,:)
!
   type (ScatterInsertStruct), pointer :: head_embed
!
   dsize = size(mat,1)
   if (dsize /= host_kmax_kkr) then
!     ----------------------------------------------------------------
      call ErrorHandler('embedScatterInHostMedium',                   &
                        'mat has different size than the t matrix',   &
                        dsize,host_kmax_kkr)
!     ----------------------------------------------------------------
   endif
!
   head_embed => EmbeddedCluster
   found = .false.
   LOOP_i: do i = 1, NumEmbeddedAtoms
      if (head_embed%local_index == id) then
         found = .true.
         exit LOOP_i
      else
         head_embed => head_embed%next
      endif
   enddo LOOP_i
   if (.not.found) then
!     ----------------------------------------------------------------
      call ErrorHandler('calEmbeddedSiteMatrix','The atom is not found in database',id)
!     ----------------------------------------------------------------
   endif
!
   if (present(compute_x)) then
      Xmat = compute_x
   else
      Xmat = .false.
   endif
!
   jd = EmbeddedClusterMatrixBand(id)%global_index ! For printing the diagonal block
   ni = gid_array(jd)
   tau_c => EmbeddedClusterMatrixBand(id)%MatrixBlock(ni)%tau_l
   m_diff => head_embed%diff_TmatInv
!
   if (Xmat) then
!     ================================================================
!     mat = [1 + m_diff*tau_c]^{-1}*m_diff
!     ----------------------------------------------------------------
      call computeAprojB('L',dsize,m_diff,tau_c(:,:,1),mat)
!     ----------------------------------------------------------------
   else
!     ================================================================
!     mat = [1 + tau_c*m_diff]^{-1}*tau_c
!     ----------------------------------------------------------------
      call computeAprojB('L',dsize,tau_c(:,:,1),m_diff,mat)
!     ----------------------------------------------------------------
   endif
!
   nullify(head_embed, tau_c, m_diff)
!
   end subroutine calEmbeddedSiteMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calEmbeddedClusterMatrix()
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
!
   implicit none
!
   integer (kind=IntKind) :: ni, nj, ki, kj, i, j, nr, nc, nb
   integer (kind=IntKind) :: dsize, kmaxi_ns, kmaxj_ns
!
   complex (kind=CmplxKind), pointer :: tau_a(:,:,:)
   complex (kind=CmplxKind), pointer :: dmat(:,:)
   complex (kind=CmplxKind), pointer :: tau_c(:,:,:)
!
   type (ScatterInsertStruct), pointer :: pi_embed, pj_embed
!
   dsize = 0
   pi_embed => EmbeddedCluster
   do i = 1, NumEmbeddedAtoms
      dsize = dsize + pi_embed%kmax_kkr*nSpinCant
      pi_embed => pi_embed%next
   enddo
!
   TMP_MatrixBand = CZERO
   do i = 1, dsize
      TMP_MatrixBand((i-1)*dsize+i) = CONE
   enddo
!
   pj_embed => EmbeddedCluster
   nc = 0
   do j = 1, NumEmbeddedAtoms
      nj = pj_embed%local_index
      kj = pj_embed%block_index
      kmaxj_ns = EmbeddedClusterMatrixBand(nj)%MatrixBlock(kj)%kmax_kkr*nSpinCant
      pi_embed => EmbeddedCluster
      nr = 0
      do i = 1, NumEmbeddedAtoms
         ni = pi_embed%local_index
         ki = pi_embed%block_index
         nb = nc*dsize+nr
         tau_c => EmbeddedClusterMatrixBand(nj)%MatrixBlock(ki)%tau_l
         kmaxi_ns = EmbeddedClusterMatrixBand(nj)%MatrixBlock(ki)%kmax_kkr*nSpinCant
!        -------------------------------------------------------------
         call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxj_ns, CONE,     &
                    tau_c(1,1,1), kmaxi_ns, pj_embed%diff_TmatInv, kmaxj_ns, &
                    CONE, TMP_MatrixBand(nb+1), dsize)
!        -------------------------------------------------------------
         pi_embed => pi_embed%next
         nr = nr + kmaxi_ns
      enddo
      pj_embed => pj_embed%next
      nc = nc + kmaxj_ns
   enddo
!
   dmat => aliasArray2_c(TMP_MatrixBand,dsize,dsize)
!  -------------------------------------------------------------------
   call MtxInv_LU(dmat,dsize)
!  -------------------------------------------------------------------
!
   pj_embed => EmbeddedCluster
   nc = 0
   do j = 1, NumEmbeddedAtoms
      nj = pj_embed%local_index
      kj = pj_embed%block_index
      kmaxj_ns = EmbeddedClusterMatrixBand(nj)%MatrixBlock(kj)%kmax_kkr*nSpinCant
      pi_embed => EmbeddedCluster
      nr = 0
      do i = 1, NumEmbeddedAtoms
         ni = pi_embed%local_index
         ki = pi_embed%block_index
         nb = nc*dsize+nr
         tau_c => EmbeddedClusterMatrixBand(nj)%MatrixBlock(ki)%tau_l
         tau_a => EmbeddedClusterMatrixBand(nj)%MatrixBlock(ki)%kau_l
         kmaxi_ns = EmbeddedClusterMatrixBand(nj)%MatrixBlock(ki)%kmax_kkr*nSpinCant
!        -------------------------------------------------------------
         call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxj_ns, CONE,     &
                    TMP_MatrixBand(nb+1), dsize,                      &
                    tau_c(1,1,1), kmaxi_ns, CZERO, tau_a(1,1,1), kmaxi_ns)
!        -------------------------------------------------------------
         pi_embed => pi_embed%next
         nr = nr + kmaxi_ns
      enddo
      pj_embed => pj_embed%next
      nc = nc + kmaxj_ns
   enddo
!
   end subroutine calEmbeddedClusterMatrix
!  ===================================================================
end module EmbeddedClusterModule
