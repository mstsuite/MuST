module MediumMatrixModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1
   use MathParamModule, only : two, half, third
   use MathParamModule, only : ten2m8, ten2m6, ten2m14, ten2m9, ten2m10
   use MathParamModule, only : sqrtm1, pi2, pi
   use TimerModule, only: getTime
   use WriteMatrixModule, only  : writeMatrix
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
!
public :: initMediumMatrix,   &
          endMediumMatrix,    &
          calMediumMatrix,    &
          beginEmbedding,     &
          endEmbedding,       &
          embedScatterToMedium, &
          calEmbeddedmediumMatrix, &
          getTmatInvDiff,     &
          projectMatrix,      &
          getTauMedium
!
private
   type MatrixBlockStruct
      integer (kind=IntKind) :: lmax_kkr
      integer (kind=IntKind) :: kmax_kkr
      integer (kind=IntKind) :: row_index
      integer (kind=IntKind) :: global_index
      complex (kind=CmplxKind), pointer :: tau_c(:,:)
!     complex (kind=CmplxKind), pointer :: tau_a(:,:)
   end type MatrixBlockStruct
!
   type MatrixBandStruct
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: column_index  ! The column index in the band matrix, not the "big" KKR matrix
      integer (kind=IntKind) :: lmax_kkr
      integer (kind=IntKind) :: kmax_kkr
      type (MatrixBlockStruct), pointer :: MatrixBlock(:)
   end type MatrixBandStruct
!
   type ScmBlockStruct
      complex (kind=CmplxKind), pointer :: strcon_matrix(:,:)
   end type ScmBlockStruct
!
   type (ScmBlockStruct), allocatable :: sc_blocks(:,:)
!
   type (MatrixBandStruct), allocatable :: MatrixBand(:)  ! Matrix column band
!
   type ImpurtyInsertStruct
      integer (kind=IntKind) :: embed_index
      integer (kind=IntKind) :: local_index
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: block_index
      integer (kind=IntKind) :: kmax_kkr
      complex (kind=CmplxKind), pointer :: diff_TmatInv(:,:)
      type (ImpurtyInsertStruct), pointer :: next
      type (ImpurtyInsertStruct), pointer :: prev
   end type ImpurtyInsertStruct
!
   type (ImpurtyInsertStruct), pointer :: head_embed, tail_embed
!
   complex (kind=CmplxKind), allocatable, target :: Tau_MatrixBand(:,:,:)
   complex (kind=CmplxKind), allocatable, target :: TMP_MatrixBand(:)
   complex (kind=CmplxKind), allocatable, target :: KKR_MatrixBand(:)
   complex (kind=CmplxKind), allocatable, target :: SCM_MatrixBand(:)
   complex (kind=CmplxKind), allocatable, target :: strconrel(:,:)
!
   logical :: isRelativistic = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
   integer (kind=IntKind) :: nSpinCant
   integer (kind=IntKind) :: lmax_kkr_max, kmax_kkr_max, tsize
   integer (kind=IntKind) :: KKRMatrixSize
   integer (kind=IntKind) :: BandSize
   integer (kind=IntKind) :: KKRMatrixSizeCant
   integer (kind=IntKind) :: BandSizeCant
   integer (kind=IntKind) :: NumEmbeddedAtoms = 0
   integer (kind=IntKind) :: NumEmbeddedAtoms_save = 0
!
   integer (kind=IntKind), allocatable :: ip_array(:) ! Relates a matrix block row index to the mapped processor 
                                                      ! index (0, 1, ...)
   integer (kind=IntKind), allocatable :: il_array(:) ! Relates a matrix block row index to the local site index 
                                                      ! on the mapped processor
   integer (kind=IntKind), allocatable :: id_array(:) ! Relates a matrix block row index to the global index of 
                                                      ! the corresponding site
   integer (kind=IntKind), allocatable :: jd_array(:) ! Relates a local site index to the global index of the site
   integer (kind=IntKind), allocatable :: gid_array(:)! Relates a global index to the corresponding matrix block 
                                                      ! row index
   integer (kind=IntKind), allocatable :: lmaxi_array(:)
   integer (kind=IntKind), allocatable :: lmaxj_array(:)
!
   integer (kind=IntKind) :: LWORK, LIWORK
   integer (kind=IntKind), allocatable :: IWORK(:)
   complex (kind=CmplxKInd), allocatable, target :: WORK(:)
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
!
   complex (kind=CmplxKind) :: energy, kappa
   complex (kind=CmplxKind), allocatable, target :: tmat_g(:,:) ! t-matrix in global frame
   complex (kind=CmplxKind), allocatable, target :: tinv(:,:)   ! t-matrix inverse, in global spin frame, of the 
                                                                ! medium sites on local CPU
!
#ifdef USE_SCALAPACK
!  ===================================================================
!  *****      ScaLAPACK parameters
!  ===================================================================
   integer (kind=IntKind), parameter :: DLEN_ = 9
   integer (kind=IntKind) :: ICTXT
   integer (kind=IntKind) :: NPROW
   integer (kind=IntKind) :: NPCOL
   integer (kind=IntKind) :: MYROW
   integer (kind=IntKind) :: MYCOL
   integer (kind=IntKind) :: DESC_A( DLEN_ )
#endif
   integer (kind=IntKind) :: INFO
   integer (kind=IntKind), allocatable :: IPVT(:)
!  ===================================================================
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMediumMatrix(cant, rel, istop, iprint)
!  ===================================================================
   use MPPModule, only : MyPE
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalMaxInGroup, getGroupCommunicator
   use GroupCommModule, only : syncAllPEsInGroup
   use MediumHostModule, only  : getNumSites, getSitePosition,              &
                                 getLmaxKKR, getLmaxPhi, getBravaisLattice, &
                                 getGlobalIndex, getLocalNumSites
   use StrConstModule, only : initStrConst
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: cant, rel
   integer (kind=IntKind), intent(in) :: iprint
!
   character (len=14) :: sname = "initMediumMatrix"
!
   integer (kind=IntKind) :: i, j, ig, na, n, k, nk, nb
   integer (kind=IntKind) :: lmaxi, kmaxi, kmaxj, t0size, nsize, kmaxi_ns
   integer (kind=IntKind) :: status, mbuf(2)
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: global_posi(:,:)
!
   complex (kind=CmplxKind), pointer :: pmat(:)
!
   stop_routine = istop
!
   GlobalNumSites = getNumSites()
   LocalNumSites = getLocalNumSites()
   nSpinCant = cant
!
   GroupID = getGroupID('Medium Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   if (rel>1) then
      isRelativistic = .true.
   endif
!
!  ===================================================================
!  Initialize the structure constant module.
!  ===================================================================
   allocate( global_posi(1:3,1:GlobalNumSites) )
   do i = 1, GlobalNumSites
      global_posi(1:3,i) = getSitePosition(i)
   enddo
   bravais(1:3,1:3) = getBravaisLattice()
!  -------------------------------------------------------------------
   call initStrConst(getLmaxPhi(), GlobalNumSites, global_posi, bravais, &
                     istop, iprint)
!  -------------------------------------------------------------------
   deallocate(global_posi)
!  ===================================================================
!
   allocate ( MatrixBand(LocalNumSites), STAT=status )
!
   allocate( id_array(GlobalNumSites), jd_array(LocalNumSites) )
   allocate( ip_array(GlobalNumSites), il_array(GlobalNumSites) )
   allocate( lmaxi_array(GlobalNumSites), lmaxj_array(LocalNumSites) )
   allocate( gid_array(GlobalNumSites), sc_blocks(GlobalNumSites,LocalNumSites))
!
   print_level= iprint
   BandSize = 0
   n = 0
   do j = 1, LocalNumSites
      MatrixBand(j)%global_index = getGlobalIndex(j)
      MatrixBand(j)%lmax_kkr = getLmaxKKR(MatrixBand(j)%global_index)
      MatrixBand(j)%kmax_kkr = (MatrixBand(j)%lmax_kkr+1)**2
      BandSize = BandSize + MatrixBand(j)%kmax_kkr
      allocate( MatrixBand(j)%MatrixBlock(GlobalNumSites) )
      MatrixBand(j)%column_index = n + 1
      n = n + MatrixBand(j)%kmax_kkr*nSpinCant
      jd_array(j) = MatrixBand(j)%global_index
      lmaxj_array(j) = MatrixBand(j)%lmax_kkr
   enddo
!
   lmax_kkr_max = 0
   KKRMatrixSize = 0
   do i = 1, GlobalNumSites
      lmaxi = getLmaxKKR(i)
      kmaxi = (lmaxi+1)**2
      lmax_kkr_max = max( lmax_kkr_max,lmaxi )
      KKRMatrixSize = KKRMatrixSize + kmaxi
   enddo
   kmax_kkr_max = (lmax_kkr_max+1)**2
!
   KKRMatrixSizeCant = KKRMatrixSize*nSpinCant
   BandSizeCant = BandSize*nSpinCant
   tsize = kmax_kkr_max*kmax_kkr_max*nSpinCant*nSpinCant
!
   allocate ( Tau_MatrixBand(tsize,GlobalNumSites,LocalNumSites) )
   allocate ( TMP_MatrixBand(KKRMatrixSizeCant*BandSizeCant+tsize+2*kmax_kkr_max*nSpinCant) )
   allocate ( KKR_MatrixBand(KKRMatrixSizeCant*BandSizeCant+kmax_kkr_max*nSpinCant) )
   allocate ( SCM_MatrixBand(KKRMatrixSize*BandSize) )
   if (isRelativistic) then
      allocate( strconrel(kmax_kkr_max*nSpinCant,kmax_kkr_max*nSpinCant) )
   endif
!
   LWORK = KKRMatrixSizeCant*BandSizeCant
   LIWORK = 2*KKRMatrixSizeCant
   allocate( WORK(1:max(LWORK,2*BandSizeCant*BandSizeCant)), IWORK(1:LIWORK) )
   WORK = CZERO; IWORK = 0
!
   allocate( tinv(tsize,LocalNumSites) )
!
!  ===================================================================
!  Set up the matrix block row/column index in such way that each
!  a band of columns of the matrix is stored on a processor
!      n = the matrix block row/column index in the band matrix
!     ig = the global index of the site that corresponds to the matrix block
!  ===================================================================
   n = 0
   nk = 0
   nb = 0
   do k = 1, NumPEsInGroup
      na = getLocalNumSites(k-1)
      do i = 1, na
         n = n + 1
         ig = getGlobalIndex(i,k-1)
         lmaxi = getLmaxKKR(ig)
         kmaxi = (lmaxi+1)**2
         kmaxi_ns = kmaxi*nSpinCant
         t0size = kmaxi_ns*kmaxi_ns
         do j = 1, LocalNumSites
            MatrixBand(j)%MatrixBlock(n)%lmax_kkr = lmaxi
            MatrixBand(j)%MatrixBlock(n)%kmax_kkr = kmaxi
            MatrixBand(j)%MatrixBlock(n)%global_index = ig
            MatrixBand(j)%MatrixBlock(n)%row_index = nk + 1
            MatrixBand(j)%MatrixBlock(n)%tau_c => aliasArray2_c(Tau_MatrixBand(:,ig,j),  &
                                                                kmaxi_ns,kmaxi_ns)
!           pmat => KKR_MatrixBand(nb+1:nb+t0size) ! KKR_MatrixBand is to store tau_a
!           MatrixBand(j)%MatrixBlock(n)%tau_a => aliasArray2_c(pmat,kmaxi_ns,kmaxi_ns)
            nb = nb + t0size
         enddo
         nk = nk + kmaxi*nSpinCant
         id_array(n) = ig
         gid_array(ig) = n
         ip_array(n) = k-1
         il_array(n) = i
         lmaxi_array(n) = lmaxi
      enddo
   enddo
!
   nsize = 0
   do j = 1, LocalNumSites
      kmaxj = (MatrixBand(j)%lmax_kkr+1)**2
      do n = 1, GlobalNumSites
         kmaxi = MatrixBand(j)%MatrixBlock(n)%kmax_kkr
         sc_blocks(n,j)%strcon_matrix => aliasArray2_c(SCM_MatrixBand(nsize+1:nsize+kmaxi*kmaxj), &
                                                       kmaxi,kmaxj)
         nsize = nsize + kmaxi*kmaxj
      enddo
   enddo
!
!  ===================================================================
   if ( print_level >0 .or. MyPE==0 ) then
      write(6,*) "===================================================="
      write(6,*) "init MediumMatrix Module ::  "
      write(6,'(a,i5)') "   LocalNumSites   : ", LocalNumSites
      write(6,'(a,i5)') "   lmax            : ", lmax_kkr_max
      write(6,'(a,i5)') "   n_spin_cant     : ", nSpinCant
      write(6,'(a,i5)') "   KKR_row_dim     : ", KKRMatrixSize
      write(6,'(a,i5)') "   KKR_col_dim     : ", BandSize
      do j = 1,LocalNumSites
         write(6,'(a,i5)') "   Local Site Index : ", j
         write(6,'(a,i5)') "       LocalCol_Ind : ", MatrixBand(j)%column_index
         write(6,'(a,i5)') "       kmax         : ", MatrixBand(j)%kmax_kkr
         write(6,'(a,i5)') "       lmax         : ", MatrixBand(j)%lmax_kkr
         write(6,'(a,i5)') "       global index : ", MatrixBand(j)%global_index
         if (GlobalNumSites < 101) then
            do i = 1, GlobalNumSites
               write(6,'(a,4i5)')"   block index, starting row, kmax, global Index : ", &
                                 i, MatrixBand(j)%MatrixBlock(i)%row_index,             &
                                 MatrixBand(j)%MatrixBlock(i)%kmax_kkr,                 &
                                 MatrixBand(j)%MatrixBlock(i)%global_index
            enddo
         endif
      enddo
      write(6,*) "===================================================="
   endif
!
!  ===================================================================
!  check the distribution of the matrix blocks. If lmax_kkr is the same
!  for all sites and each process has the same number of sites, BandSize
!  will be the same on all the processes
!  ===================================================================
   mbuf(1) = BandSizeCant
   mbuf(2) = tsize+kmax_kkr_max
!  -------------------------------------------------------------------
   call GlobalMaxInGroup(GroupID,mbuf,2)
!  -------------------------------------------------------------------
   if (mbuf(1) /= BandSizeCant) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMediumMatrix',                                    &
                        'ScaLAPACK will fail as BandSize is unevenly distrubuted', &
                         BandSizeCant)
!     ----------------------------------------------------------------
   endif
   allocate( tmat_g(mbuf(2),GlobalNumSites) ) ! Note: an extra space of
                                              ! kmax_kkr_max size is used 
                                              ! in the spin canted case when 
                                              ! calling zgemm to avoid
                                              ! out of array bound
!
!  sould be dimension of the local matrix+blocksize
   allocate( IPVT(1:KKRMatrixSizeCant+BandSizeCant) )
!
   if (NumPEsInGroup == 1) then  ! ScaLapack will not be used for 1 process case
      return
   endif
!
#ifdef USE_SCALAPACK
!  ===================================================================
!  Initialize ScaLAPACK and set up matrix distribution
!  ===================================================================
   ICTXT = getGroupCommunicator(GroupID)
   call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEsInGroup )
!  ===================================================================
   call BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!  -------------------------------------------------------------------
   if (NPROW /= 1 .or. NPCOL /= NumPEsInGroup .or. MYROW /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMediumMatrix',                              &
              'Failed: NPROW /= 1 || NPCOL /= NumPEsInGroup || MYROW /= 0',  &
              NPROW, NPCOL, MYROW)
!     ----------------------------------------------------------------
   else if (MYCOL /= MyPEinGroup) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMediumtMatrix','MYCOL /= MyPEinGroup',MYCOL,MyPEinGroup)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call syncAllPEsInGroup(GroupID)
!  -------------------------------------------------------------------
   call DESCINIT( DESC_A, KKRMatrixSizeCant, KKRMatrixSizeCant,       &
                  BandSizeCant, BandSizeCant, 0, 0, ICTXT,            &
                  KKRMatrixSizeCant, INFO )
!  -------------------------------------------------------------------
#endif
!
   NumEmbeddedAtoms = 0
   NumEmbeddedAtoms_save = 0
!
   nullify(head_embed, tail_embed)
!
   end subroutine initMediumMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endMediumMatrix()
!  ===================================================================
   use StrConstModule, only :  endStrConst
!
   implicit none
!
   integer (kind=IntKind) :: i, j
!
   do j = 1, LocalNumSites
      do i = 1, GlobalNumSites
         nullify( MatrixBand(j)%MatrixBlock(i)%tau_c )
         nullify( sc_blocks(i,j)%strcon_matrix )
      enddo
      deallocate( MatrixBand(j)%MatrixBlock )
   enddo
!
   do j = 1, NumEmbeddedAtoms-1
      deallocate(tail_embed%diff_TmatInv)
      tail_embed => tail_embed%prev
      nullify(tail_embed%next%prev)
      deallocate(tail_embed%next)
   enddo
   nullify(tail_embed)
   deallocate(head_embed%diff_TmatInv)
   deallocate(head_embed)
!
   NumEmbeddedAtoms = 0
   NumEmbeddedAtoms_save = 0
!      
   deallocate( MatrixBand, sc_blocks )
   deallocate( Tau_MatrixBand )
   deallocate( TMP_MatrixBand )
   deallocate( KKR_MatrixBand )
   deallocate( SCM_MatrixBand )
   deallocate( id_array, jd_array, ip_array, il_array )
   deallocate( lmaxi_array, lmaxj_array, gid_array )
   deallocate( tmat_g )
!
   if (isRelativistic) then
      deallocate( strconrel )
   endif
!
   deallocate( WORK, IWORK )
   deallocate( tinv )
!
   isRelativistic = .false.
!
   deallocate( IPVT )
!
!  -------------------------------------------------------------------
   call endStrConst()
!  -------------------------------------------------------------------
!
   if (NumPEsInGroup == 1) then  ! ScaLapack will not be used for 1 process case
      return
   endif
!
#ifdef USE_SCALAPACK
!  -------------------------------------------------------------------
   call BLACS_GRIDEXIT(ICTXT)
!  -------------------------------------------------------------------
!  call BLACS_EXIT(1)
!
#endif
!
   end subroutine endMediumMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupSSMatrixBuf(getTMatrix)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use WriteMatrixModule,  only : writeMatrix
!
   use MatrixInverseModule, only : MtxInv_LU
!
   implicit none
!
   integer (kind=IntKind) :: t0size, i, ig, nsize
!
   complex (kind=CmplxKind), pointer :: tm(:,:)
   complex (kind=CmplxKind), pointer :: p_tinv(:,:)
!
   interface
      function getTMatrix(site,atom,dsize) result(tmat)
         use KindParamModule, only : IntKind, CmplxKind
         integer (kind=IntKind), intent(in) :: site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: tmat(:,:)
      end function getTMatrix
   end interface
!
   tmat_g = CZERO
   do i = 1, LocalNumSites
!     ================================================================
!     Obtain the t-matrix in Global frame
!     ================================================================
!     kmax_kkr = MatrixBand(i)%kmax_kkr
      tm => getTMatrix(site=i,atom=0,dsize=nsize) ! If "i" is a CPA medium site, it returns Tcpa
      t0size = nsize*nsize
!     ----------------------------------------------------------------
      p_tinv => aliasArray2_c(tinv(:,i),nsize,nsize)
      p_tinv = tm
      call MtxInv_LU(p_tinv,nsize)
!     ----------------------------------------------------------------
      ig = MatrixBand(i)%global_index ! "ig" is the global index of the corresponding site
!     ----------------------------------------------------------------
      call zcopy(t0size,tm,1,tmat_g(:,ig),1)
!     ----------------------------------------------------------------
   enddo
!
   nullify(tm)
!
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,tmat_g,size(tmat_g,1),GlobalNumSites)
!  -------------------------------------------------------------------
!
   end subroutine setupSSMatrixBuf
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calMediumMatrix(e,getTMatrix,diag_only)
!  ===================================================================
   use MPPModule, only : MyPE
   use MatrixModule, only : computeAStarT
   use BZoneModule, only : getNumKs, getAllWeights, getAllKPoints,    &
                           getNumRotations, getWeightSum
   use ProcMappingModule, only : isKPointOnMyProc, getNumKsOnMyProc,  &
                                 getKPointIndex, getNumRedundantKsOnMyProc
   use GroupCommModule, only : getGroupID, GlobalSumInGroup, getMyPEinGroup
   use StrConstModule, only : getStrConstMatrix
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   character (len=12) :: sname = "calMediumMatrix"
!
   integer (kind=IntKind) :: k_loc, k, i, row, col, MyPEinKGroup
   integer (kind=IntKind) :: NumKs, kGID, NumKsOnMyProc, NumRedunKs
!
   real (kind=RealKind), pointer :: kpts(:,:), weight(:)
   real (kind=RealKind) :: kfac, kaij
   real (kind=RealKind) :: weightSum, kvec(1:3), aij(3)
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: wfac, efac
   complex (kind=CmplxKind), pointer :: rotmat(:,:), scm(:,:)
!
   logical, intent(in), optional :: diag_only
   logical :: DiagonalOnly = .false.
!
   interface
      function getTMatrix(site,atom,dsize) result(tmat)
         use KindParamModule, only : IntKind, CmplxKind
         integer (kind=IntKind), intent(in) :: site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: tmat(:,:)
      end function getTMatrix
   end interface
!
   energy = e
   kappa = sqrt(e)
!
   if (GlobalNumSites == 1) then
      DiagonalOnly = .true.
   else if (present(diag_only)) then
      DiagonalOnly = diag_only
   else
      DiagonalOnly = .false.
   endif
!
!  -------------------------------------------------------------------
   call setupSSMatrixBuf(getTMatrix) ! Setup single site scattering 
                                     ! matrix data set
!  -------------------------------------------------------------------
!
   kGID = getGroupID('K-Mesh')
   NumKsOnMyProc = getNumKsOnMyProc()
   NumRedunKs = getNumRedundantKsOnMyProc()
   MyPEinKGroup = getMyPEinGroup(kGID)
!
   NumKs = getNumKs()
   kpts => getAllKPoints(kfac)
   weight => getAllWeights()
   weightSum = getWeightSum()
!
!  ===================================================================
!  Loop over k-points on mesh
!  ===================================================================
   KKR_MatrixBand = CZERO
   do k_loc = 1,NumKsOnMyProc
!     ================================================================
!     Normorlize BZ integration weights
!     ================================================================
      k = getKPointIndex(k_loc)
      kvec(1:3) = kpts(1:3,k)*kfac
!     ================================================================
!     get structure constant matrix for the k-point and energy
      do col = 1, LocalNumSites
         do row = 1, GlobalNumSites
!           ----------------------------------------------------------
            scm => getStrConstMatrix(kvec,kappa,id_array(row),jd_array(col), &
                                     lmaxi_array(row),lmaxj_array(col),aij)
!           ----------------------------------------------------------
!           kaij = kvec(1)*aij(1)+kvec(2)*aij(2)+kvec(3)*aij(3)
!           efac = exp(sqrtm1*kaij)
!           sc_blocks(row,col)%strcon_matrix = scm*efac
            sc_blocks(row,col)%strcon_matrix = scm
         enddo
      enddo
!
      wfac = weight(k)/weightSum
!     write(6,'(a,i4,a,3f12.5,a,2d12.5,a,2d12.5)')'k-ind = ',k,       &
!             ', kvec = ',kvec,', wfac = ',wfac,', weightsum = ',weightSum
!
!     ================================================================
!     Compute the the following KKR matrix, which is stored in TMP_MatrixBand
!         TMP_MatrixBand = [1 - (B+i*kappa) * Tmat]^{-1}
!     and is parallelly distributed.
!     ----------------------------------------------------------------
      call computeMatrixBand(TMP_MatrixBand)
!     ----------------------------------------------------------------
!     call checkMatrixBandRotation(kvec,TMP_MatrixBand)
!     ----------------------------------------------------------------
!
!     ================================================================
!     Sum TMP_MatrixBand over the k-points across the processors that 
!     take care different k-points. The loop will continue for the 
!     redundant k-points, if there are any. The redundant k-points are 
!     needed for the load balance purpose. The result is stored in 
!     KKR_MatrixBand
!     ================================================================
      if (k_loc <= NumKsOnMyProc - NumRedunKs .or. MyPEinKGroup == 0) then
!        -------------------------------------------------------------
         call zaxpy(KKRMatrixSizeCant*BandSizeCant,wfac,TMP_MatrixBand,1, &
                    KKR_MatrixBand,1)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call zaxpy(KKRMatrixSizeCant*BandSizeCant,CZERO,TMP_MatrixBand,1, &
                    KKR_MatrixBand,1)
!        -------------------------------------------------------------
      endif
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(kGID,KKR_MatrixBand,KKRMatrixSizeCant*BandSizeCant)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Compute the the Tau matrix as follows:
!     Tau = Tmat * KKR_MatrixBand
!  ===================================================================
   call computeTauMatrix()
!  -------------------------------------------------------------------
!
   if (getNumRotations() > 1) then
      if (DiagonalOnly) then
!        -------------------------------------------------------------
         call sumIBZRotation()
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call ErrorHandler('calMediumMatrix',                         &
                           'IBZ rotation is not applicable to non-diagonal Tau calculation')
!        -------------------------------------------------------------
      endif
   endif
!
   nullify( weight, kpts )
!
   if (trim(stop_routine) ==trim(sname)) then
      stop 'At calMediumMatrix'
   endif
!
   end subroutine calMediumMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeTauMatrix()
!  ===================================================================
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind) :: i, j, ig, ni, nr, n
   integer (kind=IntKind) :: kmaxj, kmaxj_ns, t0size
!
   complex (kind=CmplxKind), pointer :: tau_c(:,:)
!
!  ===================================================================
!  calculate Tau Matrix as follows:
!    Tau = [t(e)^{-1} - (B(k,e)+i*kappa)]^{-1}
!        = t(e) * [1 - (B(k,e)+i*kappa)*t(e)]^{-1}
!    !!  = [ kappa*Jost(e)*S(e)^{-1} - (B(k,e)+i*kappa) ]^{-1}
!    !!  = -[ kappa*C(e)*S(e)^{-1}) + B(k,e) ]^{-1}
!
!  Note:
!    KKR_MatrixBand = [ 1 - (B(k,e)+i*kappa) * t(e) ]^{-1}
!    The resulting Tauij matrices are in global frame reference in spin space
!  ===================================================================
   do i = 1, LocalNumSites
      ni = MatrixBand(i)%column_index
      do j = 1, GlobalNumSites
         kmaxj = MatrixBand(i)%MatrixBlock(j)%kmax_kkr
         kmaxj_ns = kmaxj*nSpinCant
         tau_c => MatrixBand(i)%MatrixBlock(j)%tau_c
         ig = MatrixBand(i)%MatrixBlock(j)%global_index
         nr = MatrixBand(i)%MatrixBlock(j)%row_index
         n = KKRMatrixSizeCant*(ni-1)+nr
!        -------------------------------------------------------------
         call zgemm('n', 'n', kmaxj_ns, kmaxj_ns, kmaxj_ns, CONE,     &
                    tmat_g(1,ig), kmaxj_ns,                           &
                    KKR_MatrixBand(n), KKRMatrixSizeCant,             &
                    CZERO, tau_c, kmaxj_ns)
!        -------------------------------------------------------------
      enddo
   enddo
!
   nullify( tau_c )
!
   end subroutine computeTauMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeMatrixBand(p_MatrixBand)
!  ===================================================================
   use MatrixModule, only : setupUnitMatrix
!
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind) :: kmaxj, kmaxj_ns, kmaxi, kmaxi_ns, t0size, n_BT
   integer (kind=IntKind) :: j, nj, ni, i, is, ig, nc, kl, klp, n, np, nr
!
   complex (kind=CmplxKind), pointer :: strcon(:,:)
   complex (kind=CmplxKind), pointer :: p_BT(:)
   complex (kind=CmplxKind), pointer :: p_Tmat(:)
   complex (kind=CmplxKind), pointer :: p_1mikt(:) ! 1 - i*kappa*Tmat
   complex (kind=CmplxKind), intent(out), target :: p_MatrixBand(:)
   complex (kind=CmplxKInd) :: cfac
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
!  ===================================================================
!  calculate the following KKR Matrix (or the M-matrix):
!    p_MatrixBand = [1 - (B(e,k)+i*kappa) * t(e)]^{-1}
!  ===================================================================
   p_MatrixBand = CZERO
   cfac = SQRTm1*kappa
   do j = 1, LocalNumSites
      kmaxj = MatrixBand(j)%kmax_kkr
      kmaxj_ns = kmaxj*nSpinCant
      t0size = kmaxj_ns*kmaxj_ns
      ig = MatrixBand(j)%global_index ! "ig" is the global index of the corresponding site
      p_Tmat => tmat_g(:,ig)
      p_1mikt => WORK(1:t0size)
!     ================================================================
!     setup p_1mikt = 1 - i*kappa*tmat
!     ================================================================
      p_1mikt = CZERO
      do kl = 1, kmaxj_ns
         p_1mikt((kl-1)*kmaxj_ns+kl) = CONE
      enddo
      call zaxpy(t0size,-cfac,p_Tmat,1,p_1mikt,1)
!     ================================================================
      nj = MatrixBand(j)%column_index
      do i = 1, GlobalNumSites        ! "i" is the row index of the matrix block
         nr = MatrixBand(j)%MatrixBlock(i)%row_index
         kmaxi = MatrixBand(j)%MatrixBlock(i)%kmax_kkr
         kmaxi_ns = kmaxi*nSpinCant
!        =============================================================
!        Obtain:
!           p_BT = - B(e,k) * t(e)
!        =============================================================
         p_BT => p_MatrixBand
         n_BT = KKRMatrixSizeCant*(nj-1)+nr-1
         strcon => sc_blocks(i,j)%strcon_matrix(:,:)
!        call writeMatrix('strcon',strcon,kmaxi,kmaxj)
         if (isRelativistic) then
!           ----------------------------------------------------------
            call convertGijToRel(strcon, strconrel, kmaxi, kmaxj, energy)
!           ----------------------------------------------------------
            strcon => strconrel
!           ----------------------------------------------------------
            call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxj_ns, -CONE, &
                       strcon, kmaxi_ns, p_Tmat, kmaxj_ns, CONE, p_BT, kmaxi_ns)
!           ----------------------------------------------------------
         else
!           ==========================================================
!           In spin canted case, for is = 1, p_BT = -B * (t11, t12)
!                                       = 2,      = -B * (t21, t22)
!           ==========================================================
            do is = 1, nSpinCant
!              Note: both p_Tmat and p_BT has an extra space of size at 
!              least kmaxj, which is necessary to avoid out of array 
!              bound in the spin canted case
!              -------------------------------------------------------
               call zgemm('n', 'n', kmaxi, kmaxj_ns, kmaxj, -CONE,    &
                          strcon, kmaxi,                              &
                          p_Tmat((is-1)*kmaxj+1), kmaxj_ns,           &
                          CONE, p_BT(n_BT+(is-1)*kmaxi+1), KKRMatrixSizeCant)
!              -------------------------------------------------------
            enddo
         endif
      enddo
!     ================================================================
!     Add 1 - i*kappa*tmat to -B * tmat, so that
!         p_MatrixBand = 1 - i*kappa*tmat - B*tmat
!     ================================================================
      nc = gid_array(ig) ! "nc" is the column index of the site block in the big matrix
      n = KKRMatrixSizeCant*(nj-1)+MatrixBand(j)%MatrixBlock(nc)%row_index-1
      np = 0
      do kl = 1, kmaxj_ns
         do klp = 1, kmaxj_ns
            p_BT(n+klp) = p_BT(n+klp) + p_1mikt(np+klp)
         enddo
         n = n + KKRMatrixSizeCant
         np = np + kmaxj_ns
      enddo
   enddo
!
   nullify(p_BT, p_Tmat, p_1mikt)
!
   if (NumPEsInGroup == 1) then  ! BandSizeCant = KKRMatrixSizeCant
!     ----------------------------------------------------------------
      call ZGETRF(KKRMatrixSizeCant, KKRMatrixSizeCant,               &
                  p_MatrixBand, KKRMatrixSizeCant, IPVT, INFO)
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('computeMatrixBand','Failed in ZGETRF',INFO)
!        -------------------------------------------------------------
      endif
   else
#ifdef USE_SCALAPACK
!     ----------------------------------------------------------------
      call PZGETRF(KKRMatrixSizeCant, KKRMatrixSizeCant,              &
!???  call PZGETRF(KKRMatrixSizeCant, BandSizeCant,                   &
                   p_MatrixBand, 1, 1, DESC_A, IPVT, INFO)
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('computeMatrixBand','Failed in PZGETRF',INFO)
!        -------------------------------------------------------------
      endif
#else
!     ----------------------------------------------------------------
      call ErrorHandler('computeMatrixBand','Compiling with -DUSE_SCALAPACK is needed!')
!     ----------------------------------------------------------------
      call ZGETRF(KKRMatrixSizeCant, KKRMatrixSizeCant,               &
                  p_MatrixBand, KKRMatrixSizeCant, IPVT, INFO)
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('computeMatrixBand','Failed in ZGETRF',INFO)
!        -------------------------------------------------------------
      endif
#endif
   endif
!
   if (NumPEsInGroup == 1) then  ! BandSizeCant = KKRMatrixSizeCant
!     ----------------------------------------------------------------
      call ZGETRI( KKRMatrixSizeCant, p_MatrixBand, KKRMatrixSizeCant, IPVT, WORK, LWORK, INFO )
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('computeMatrixBand','Failed in ZGETRS',INFO)
!        -------------------------------------------------------------
      endif
   else
#ifdef USE_SCALAPACK
!     ----------------------------------------------------------------
      call PZGETRI(KKRMatrixSizeCant, p_MatrixBand, 1, 1,             &
                   DESC_A, IPVT, WORK, LWORK, IWORK, LIWORK, INFO )
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('computeMatrixBand','Failed in PZGETRS',INFO)
!        -------------------------------------------------------------
      endif
#else
!     ----------------------------------------------------------------
      call ErrorHandler('computeMatrixBand','Compiling with -DUSE_SCALAPACK is needed!')
!     ----------------------------------------------------------------
      call ZGETRI( KKRMatrixSizeCant, p_MatrixBand, KKRMatrixSizeCant, IPVT, WORK, LWORK, INFO )
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('computeMatrixBand','Failed in ZGETRS',INFO)
!        -------------------------------------------------------------
      endif
#endif
   endif
!
   end subroutine computeMatrixBand
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sumIBZRotation()
!  ===================================================================
   use MatrixModule, only : computeUAUtc
   use IBZRotationModule, only : getNumIBZRotations, getIBZRotationMatrix
!
   implicit none
!
   integer (kind=IntKind) :: i, ig, id, irot, nrot, np, is, js
   integer (kind=IntKind) :: kmaxi, kmaxi_ns
!
   complex (kind=CmplxKind), pointer :: w0(:), tau(:), rotmat(:,:)
   complex (kind=CmplxKind), pointer :: sub_tau(:,:), sub_w0(:,:)
   complex (kind=CmplxKind) :: cfac
!
   nrot = getNumIBZRotations()
   cfac = CONE/real(nrot,RealKind)
!
   do i = 1, LocalNumSites
      ig = MatrixBand(i)%global_index
      id = gid_array(ig)  ! Here, we only consider&rotate the diagonal blocks of
                          ! the band matrix. For non-diagonal matrix blocks, there 
                          ! is factor of exp(i*(1-Rot(k_vector))*aij_vector) needs 
                          ! to be applied to each transformation and each k point
      kmaxi = MatrixBand(i)%kmax_kkr
      kmaxi_ns = kmaxi*nSpinCant
      if (kmaxi /= MatrixBand(i)%MatrixBlock(id)%kmax_kkr) then
!        -------------------------------------------------------------
         call ErrorHandler('sumIBZRotation','Inconsistent kkr size',  &
                           kmaxi,MatrixBand(i)%MatrixBlock(id)%kmax_kkr)
!        -------------------------------------------------------------
      endif
!     ================================================================
!     In the spin canted case, an extra space, at least of size kmaxi,
!     is needed to be part of tau and w0 for calling computeUAUtc...
!     ================================================================
      np = kmaxi_ns*kmaxi_ns+(nSpinCant-1)*kmaxi_ns
      w0 => TMP_MatrixBand(1:np)
      tau => TMP_MatrixBand(np+1:2*np)
!     ----------------------------------------------------------------
      call zcopy(kmaxi_ns*kmaxi_ns,MatrixBand(i)%MatrixBlock(id)%tau_c,1,&
                 tau,1)
!     ----------------------------------------------------------------
      w0 = CZERO
      do irot = 1, nrot
         rotmat => getIBZRotationMatrix('c',irot)
         do js = 1, nSpinCant
            np = (js-1)*kmaxi_ns*kmaxi
            do is = 1, nSpinCant
               sub_w0 => aliasArray2_c(w0(np+(is-1)*kmaxi+1:),kmaxi_ns,kmaxi)
               sub_tau => aliasArray2_c(tau(np+(is-1)*kmaxi+1:),kmaxi_ns,kmaxi)
!              -------------------------------------------------------
               call computeUAUtc(rotmat,kmaxi,kmaxi,rotmat,kmaxi,cfac, &
                                 sub_tau,kmaxi_ns,CONE,sub_w0,kmaxi_ns,WORK)
!              -------------------------------------------------------
            enddo
         enddo
      enddo
!     ----------------------------------------------------------------
      call zcopy(kmaxi_ns*kmaxi_ns,w0,1,MatrixBand(i)%MatrixBlock(id)%tau_c,1)
!     ----------------------------------------------------------------
   enddo
!
   nullify(w0, rotmat, tau, sub_w0, sub_tau)
!
   end subroutine sumIBZRotation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTauMedium(local_id,global_id,dsize) result(tau)
!  ===================================================================
   use Atom2ProcModule, only : getGlobalIndex
   implicit none
!
   integer (kind=IntKind), intent(in) :: local_id ! Local index
   integer (kind=IntKind), intent(in), optional :: global_id ! Global index
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind) :: ni, jd
!
   complex (kind=CmplxKind), pointer :: tau(:,:)
!
   if (local_id < 1 .or. local_id > LocalNumSites) then
      call ErrorHandler('getTau','invalid local index',local_id)
   else if (present(global_id)) then
      if (global_id < 1 .or. global_id > GlobalNumSites) then
         call ErrorHandler('getTau','invalid global index',global_id)
      else
         jd = global_id
      endif
   else
      jd = MatrixBand(local_id)%global_index ! For returning the diagonal block
   endif
!
   ni = gid_array(jd)
   tau => MatrixBand(local_id)%MatrixBlock(ni)%tau_c
!
   if (present(dsize)) then
      dsize = MatrixBand(local_id)%MatrixBlock(ni)%kmax_kkr*nSpinCant
   endif
!
   end function getTauMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printMediumMatrix(i,j)
!  ===================================================================
   use Atom2ProcModule, only : getAtom2ProcInGroup, getLocalIndex
   implicit none
!
   integer (kind=IntKind), intent(in) :: i ! Local index
   integer (kind=IntKind), intent(in), optional :: j ! Global index
   integer (kind=IntKind) :: ni, ig, id, jd, js
   integer (kind=IntKind) :: kmaxi, kmaxj, kmaxi_ns, kmaxj_ns
!
   complex (kind=CmplxKind), pointer :: tau(:,:)
!
   if (i < 1 .or. i > LocalNumSites) then
      call ErrorHandler('printMediumMatrix','invalid local index i',i)
   else if (present(j)) then
      if (j < 1 .or. j > GlobalNumSites) then
         call ErrorHandler('printMediumMatrix','invalid global index j',j)
      else
         jd = j
      endif
   else
      jd = MatrixBand(i)%global_index ! For printing the diagonal block
   endif
!
   ni = gid_array(jd)
   kmaxj = MatrixBand(i)%kmax_kkr
   kmaxj_ns = kmaxj*nSpinCant
   tau => MatrixBand(i)%MatrixBlock(ni)%tau_c
!
   write(6,'(/,a)')'***************************************'
   write(6,'( a )')'*    Output from printMediumMatrix    *'
   write(6,'(a,/)')'***************************************'
   write(6,'(80(''=''))')
!  
   write(6,*)'energy   ::',energy
   write(6,'(80(''=''))')
!
   write(6,'(/,a,i4,2x,i4,/)')"    Sites i :: ",MatrixBand(i)%global_index
   call writeMatrix('Tau-matrix in Global Frame of Reference in Spin Space ::', &
                     tau,kmaxj_ns,kmaxj_ns )
!
   nullify(tau)
!
   end subroutine printMediumMatrix
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
   if (.not.associated(head_embed) .and. NumEmbeddedAtoms_save > 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('beginEmbedding','inconsistent linked data',  &
                        NumEmbeddedAtoms_save)
!     ----------------------------------------------------------------
   endif
!
   end subroutine beginEmbedding
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endEmbedding()
!  ===================================================================
   implicit none
!
   if (NumEmbeddedAtoms_save < NumEmbeddedAtoms) then
!     ----------------------------------------------------------------
      call ErrorHandler('endEmbedding','inconsistent linked data',    &
                        NumEmbeddedAtoms_save,NumEmbeddedAtoms)
!     ----------------------------------------------------------------
   endif
!
   end subroutine endEmbedding
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine embedScatterToMedium(id,tmat_a,dsize)
!  ===================================================================
!
!  place a species described by tmat at site "id" of the medium
!
!  *******************************************************************
   use MatrixInverseModule, only : MtxInv_LU
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, dsize
   integer (kind=IntKind) :: ni, n, ig
!
   complex (kind=CmplxKind), target :: tmat_a(:,:)
   complex (kind=CmplxKind), pointer :: diff_tinv(:,:)
   complex (kind=CmplxKind), pointer :: p_tcinv(:,:)
   complex (kind=CmplxKind), pointer :: p_tainv(:,:)
   complex (kind=CmplxKind), pointer :: tau_c(:,:)
!
   if (.not.associated(head_embed)) then
      allocate(head_embed)
      nullify(head_embed%prev)
      tail_embed => head_embed
      NumEmbeddedAtoms_save = 1  ! number of linked data
      nullify(tail_embed%next)
      allocate(tail_embed%diff_TmatInv(dsize,dsize))
   else if (NumEmbeddedAtoms == 0) then
      tail_embed => head_embed
   else if (NumEmbeddedAtoms == NumEmbeddedAtoms_save) then
      allocate(tail_embed%next)
      tail_embed%next%prev => tail_embed
      tail_embed => tail_embed%next
      NumEmbeddedAtoms_save = NumEmbeddedAtoms_save + 1  ! number of linked data
      nullify(tail_embed%next)
      allocate(tail_embed%diff_TmatInv(dsize,dsize))
   else
      tail_embed => tail_embed%next
   endif
!
   NumEmbeddedAtoms = NumEmbeddedAtoms + 1
!
   tail_embed%kmax_kkr = dsize/nSpinCant
   tail_embed%embed_index = NumEmbeddedAtoms
   tail_embed%local_index = id
   ig = MatrixBand(id)%global_index
   tail_embed%global_index = ig
   tail_embed%block_index = gid_array(ig)
!
!  -------------------------------------------------------------------
   p_tcinv => aliasArray2_c(tinv(:,id),dsize,dsize)
!  -------------------------------------------------------------------
   p_tainv => aliasArray2_c(WORK,dsize,dsize)
   p_tainv = tmat_a
   call MtxInv_LU(p_tainv,dsize)
!  -------------------------------------------------------------------
!
   tail_embed%diff_TmatInv = p_tainv - p_tcinv
!
   end subroutine embedScatterToMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTmatInvDiff(id) result(diff)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), pointer :: diff(:,:)
   type (ImpurtyInsertStruct), pointer :: present_embed
!
   present_embed => tail_embed
   if (present_embed%local_index /= id) then
      present_embed => head_embed
      LOOP_i: do i = 1, NumEmbeddedAtoms
         if (present_embed%local_index == id) then
            exit LOOP_i
         else
            present_embed => present_embed%next
         endif
      enddo LOOP_i
   endif
   diff => present_embed%diff_TmatInv
   nullify(present_embed)
!
   end function getTmatInvDiff
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine projectMatrix(c,id,mat,proj,dsize)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
!
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind), intent(in) :: id, dsize
   integer (kind=IntKind) :: jd, ni, i
!
   complex (kind=CmplxKind), intent(in), target :: mat(:,:)
   complex (kind=CmplxKind), intent(in), target :: proj(:,:)
   complex (kind=CmplxKind), pointer :: tau_c(:,:), xm(:,:)
!
   if (dsize /= size(mat,1)) then
!     ----------------------------------------------------------------
      call ErrorHandler('projectMatrix','inconsistent mat size',dsize,size(mat,1))
!     ----------------------------------------------------------------
   else if (dsize /= size(proj,1)) then
!     ----------------------------------------------------------------
      call ErrorHandler('projectMatrix','inconsistent proj size',dsize,size(proj,1))
!     ----------------------------------------------------------------
   else if (dsize /= MatrixBand(id)%kmax_kkr*nSpinCant) then
!     ----------------------------------------------------------------
      call ErrorHandler('projectMatrix','inconsistent tau_c size',dsize,&
                        MatrixBand(id)%kmax_kkr*nSpinCant)
!     ----------------------------------------------------------------
   endif
!
   jd = MatrixBand(id)%global_index ! For printing the diagonal block
   ni = gid_array(jd)
   tau_c => MatrixBand(id)%MatrixBlock(ni)%tau_c
!
   xm => aliasArray2_c(WORK,dsize,dsize)
   xm = CZERO
   do i = 1, dsize
      xm(i,i) = CONE
   enddo
!  -------------------------------------------------------------------
   call zgemm('n','n',dsize,dsize,dsize,CONE,tau_c,dsize,mat,dsize,   &
              CONE,xm,dsize)
!  -------------------------------------------------------------------
   call MtxInv_LU(xm,dsize)
!  -------------------------------------------------------------------
!
   if (c == 'L' .or. c == 'l') then
!     ----------------------------------------------------------------
      call zgemm('n','n',dsize,dsize,dsize,CONE,mat,dsize,xm,dsize,   &
                 CZERO,proj,dsize)
!     ----------------------------------------------------------------
   else if (c == 'R' .or. c == 'r') then
!     ----------------------------------------------------------------
      call zgemm('n','n',dsize,dsize,dsize,CONE,xm,dsize,mat,dsize,   &
                 CZERO,proj,dsize)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('projectMatrix','invalid projection type',c)
!     ----------------------------------------------------------------
   endif
!
   end subroutine projectMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calEmbeddedMediumMatrix(tau_a)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
!
   implicit none
!
   integer (kind=IntKind) :: ni, nj, ki, kj, i, j, nr, nc, nb
   integer (kind=IntKind) :: dsize, kmaxi_ns, kmaxj_ns
!
   complex (kind=CmplxKind), intent(inout), target :: tau_a(:,:)
   complex (kind=CmplxKind), pointer :: dmat(:,:)
   complex (kind=CmplxKind), pointer :: tau_c(:,:)
!
   type (ImpurtyInsertStruct), pointer :: pi_embed, pj_embed
!
   dsize = 0
   pi_embed => head_embed
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
   pj_embed => head_embed
   nc = 0
   do j = 1, NumEmbeddedAtoms
      nj = pj_embed%local_index
      kj = pj_embed%block_index
      kmaxj_ns = MatrixBand(nj)%MatrixBlock(kj)%kmax_kkr*nSpinCant
      pi_embed => head_embed
      nr = 0
      do i = 1, NumEmbeddedAtoms
         ni = pi_embed%local_index
         ki = pi_embed%block_index
         nb = nc*dsize+nr
         tau_c => MatrixBand(nj)%MatrixBlock(ki)%tau_c
         kmaxi_ns = MatrixBand(nj)%MatrixBlock(ki)%kmax_kkr*nSpinCant
!        -------------------------------------------------------------
         call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxj_ns, CONE,     &
                    tau_c, kmaxi_ns, pj_embed%diff_TmatInv, kmaxj_ns, &
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
   pj_embed => head_embed
   nc = 0
   do j = 1, NumEmbeddedAtoms
      nj = pj_embed%local_index
      kj = pj_embed%block_index
      kmaxj_ns = MatrixBand(nj)%MatrixBlock(kj)%kmax_kkr*nSpinCant
      pi_embed => head_embed
      nr = 0
      do i = 1, NumEmbeddedAtoms
         ni = pi_embed%local_index
         ki = pi_embed%block_index
         nb = nc*dsize+nr
         tau_c => MatrixBand(nj)%MatrixBlock(ki)%tau_c
!        MatrixBand(nj)%MatrixBlock(ki)%tau_a => tau_a
         kmaxi_ns = MatrixBand(nj)%MatrixBlock(ki)%kmax_kkr*nSpinCant
!        -------------------------------------------------------------
         call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxj_ns, CONE,     &
                    TMP_MatrixBand(nb+1), dsize,                      &
                    tau_c, kmaxi_ns, CZERO, tau_a, kmaxi_ns)
!        -------------------------------------------------------------
         pi_embed => pi_embed%next
         nr = nr + kmaxi_ns
      enddo
      pj_embed => pj_embed%next
      nc = nc + kmaxj_ns
   enddo
!
   end subroutine calEmbeddedMediumMatrix
!  ===================================================================
end module MediumMatrixModule
