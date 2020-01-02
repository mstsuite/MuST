module SuperLUSolverModule
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : zero, one, czero, cone, ten2m10
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
   use MPPModule, only : MyPE
   use superlu_mod
!
   implicit none
!
   public ::  initSuperLUSolver,     &
              setRhs,                &
              setRhsToZero,          &
              getCompRMtrx,          &
              getSolution,           &
              delSLUstruct,          &
              insertSLUlist,         &
              SuperLUSolution,       &
              SLUCompRowMtx_Cmplx,   &
              writeSLU2Matrix,       &
              endSuperLUSolver
!
   interface getCompRMtrx
      module procedure getCompressedRowMatrix
   end interface getCompRMtrx
!
   type SLUCompRowMtx_Cmplx
      integer (kind=IntKind) :: ind_row
      integer (kind=IntKind) :: ind_col
      integer (kind=IntKind) :: nnz
      integer (kind=IntKind) :: n_row
      integer (kind=IntKind) :: n_col
      integer (kind=IntKind),   pointer :: colind(:)
      integer (kind=IntKind),   pointer :: rowptr(:)
      complex (kind=CmplxKind), pointer :: nzval(:)
   end type SLUCompRowMtx_Cmplx
!
private
   logical :: isInverse = .false.
   logical :: Initialized = .false.
   logical :: isDiag = .false.
!
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: NpRow
   integer (kind=IntKind) :: NpCol
   integer (kind=IntKind) :: FSTRow
   integer (kind=IntKind) :: NRow_g
   integer (kind=IntKind) :: NCol_g
   integer (kind=IntKind) :: init
   integer (kind=IntKind) :: Nrhs
   integer (kind=IntKind) :: col_dim      !  real dimensions of the
   integer (kind=IntKind) :: row_dim      !! inverted matrix
   integer (kind=IntKind) :: nnz          !  number of non-zero
                                          !! elemets of the big matrix
!
   integer (superlu_ptr) :: grid
   integer (superlu_ptr) :: options
   integer (superlu_ptr) :: SOLVEstruct
   integer (superlu_ptr) :: LUstruct
   integer (superlu_ptr) :: ScalePermstruct
   integer (superlu_ptr) :: A
   integer (superlu_ptr) :: stat
!
   integer (kind=IntKind) :: ldb          ! leading dimension of B_*
!
   logical :: isAllocWorkBlock = .false.
!
   integer (kind=IntKind) :: maxWorkSize = -1
   integer (kind=IntKind), allocatable, target :: workRow(:), workCol(:)
   complex (kind=CmplxKind), allocatable, target :: workNNZ_z(:)
   real    (kind=RealKind), allocatable :: berr(:)
!
   real    (kind=RealKind) :: tol
!
#ifdef TIMING
   integer (kind=IntKind)  :: count_t
   real    (kind=RealKind) :: time_pzgssvx
   real    (kind=RealKind) :: time_superLu
#endif
!
   complex (kind=CmplxKind),allocatable, target :: B_cmplx(:,:)
!
   integer (kind=IntKind) :: n_SLU_list = 0
   integer (kind=IntKind) :: n_max_list = 0
#ifdef DEBUG
   integer (kind=IntKind) :: reallocate_count
#endif
!
   type (SLUCompRowMtx_Cmplx), target, save :: SLU
   type (SLUCompRowMtx_Cmplx), allocatable, target :: SLU_l(:)
   type (SLUCompRowMtx_Cmplx), pointer :: pSLU, pSLU_l
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSuperLUSolver( numprow, numpcol, nrow, ncol,   &
                       fst_row, max_blcks, slutol, iprint, max_sz )
!  ===================================================================
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: numprow, numpcol
   integer(kind=IntKind), intent(in) :: nrow, ncol
   integer(kind=IntKind), intent(in) :: fst_row
   integer(kind=IntKind), intent(in) :: max_blcks
   integer(kind=IntKind), intent(in) :: iprint
   integer(kind=IntKind), optional   :: max_sz
!
   real (kind=RealKind), intent(in) :: slutol
!
   character (len=17) :: sname = "initSuperLUSolver"
!
   integer(kind=IntKind) :: iam, status, i
!
   include 'mpif.h'
!
   print_level = iprint
!
!  ===================================================================
!  Create Fortran handles for the C structures used in SuperLU_DIST
!  ===================================================================
!
   call f_create_gridinfo_handle(grid)
   call f_create_options_handle(options)
   call f_create_ScalePermstruct_handle(ScalePermstruct)
   call f_create_LUstruct_handle(LUstruct)
   call f_create_SOLVEstruct_handle(SOLVEstruct)
   call f_create_SuperMatrix_handle(A)
   call f_create_SuperLUStat_handle(stat)
!
!  ===================================================================
!  Initialize the SuperLU_DIST process grid
!  ===================================================================
!
   NRow_g  = nrow
   NCol_g  = ncol
   NpRow   = numprow
   NpCol   = numpcol
   FSTRow  = fst_row
   n_max_list = max_blcks
   if (present(max_sz)) isDiag = .true.
   if (isDiag) then
      Nrhs = max_sz
   else
      Nrhs = NCol_g
   endif
!
   allocate(berr(Nrhs), STAT=status)
   call check_allocate( status, sname )
!   if ( print_level>0 ) then
      write(6,*) "Initialized SLU Matrix Inversion Module, Node:", MyPE
      write(6,*) "Nrhs = ",Nrhs
!   endif
!
!  -------------------------------------------------------------------
   call f_superlu_gridinit(MPI_COMM_WORLD, NpRow, NpCol, Grid)
!  -------------------------------------------------------------------
   call get_GridInfo(grid, iam=iam)
   if ( iam /= MyPE ) then
      call ErrorHandler('initSuperLUSolver','wrong grid',iam,MyPE)
   endif
!
   allocate( SLU_l(n_max_list) )
   nullify( SLU%colind, SLU%rowptr, SLU%nzval )
   do i = 1,n_max_list
      nullify( SLU_l(i)%colind, SLU_l(i)%rowptr, SLU_l(i)%nzval )
   enddo
!
   if ( slutol>zero ) then
      tol = slutol
   else
      tol = zero
   endif
!
!  ===================================================================
!  Set the default input options
!  ===================================================================
   call f_set_default_options(options)
!  ===================================================================
!  Set one or more options
!  ===================================================================
   call set_superlu_options(options, Fact= DOFACT, Equil= NO,           &
                            ColPerm= NATURAL, RowPerm= NOROWPERM,       &
                            ReplaceTinyPivot= YES, Trans= NOTRANS,      &
                            IterRefine= NOREFINE, SolveInitialized= NO, &
                            RefineInitialized= NO, PrintStat= NO )
#ifdef TIMING
   count_t = 0
   time_pzgssvx = zero
   time_superLu = zero
#endif
!
#ifdef DEBUG
   reallocate_count = 0
#endif
   Initialized = .true.
!
   end subroutine initSuperLuSolver
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function SuperLUSolution()                          result(pb)
!  ===================================================================
#ifdef TIMING
   use TimerModule,only : getTime
#endif
!
   implicit none
!
   complex (kind=CmplxKInd), pointer :: pb(:,:)
!
#ifdef TIMING
   real (kind=RealKind) :: t0, t1
#endif
#ifdef DEBUG
   if ( print_level>0 ) then
      write(6,*) 'Enter SuperLUSolution_z '
   endif
#endif
!
#ifdef TIMING
   t0 = getTime()
#endif
!   call delSLUstruct(n_SLU_list)
!
   n_SLU_list = 0
   call distributeSparseMatrix()
#ifdef TIMING
   t1 = getTime()
   time_superLu = t1-t0 + time_SuperLu
#endif
!
   call calSolutionSLU()
!
#ifdef TIMING
   t0 = getTime()
#endif
   call destroyDistributedMatrix()
!
!   call delSLUstruct()
   pb => B_cmplx
#ifdef TIMING
   t1 = getTime()
   time_superLu = t1-t0 + time_SuperLu
#endif
!
!   end subroutine SuperLUSolution
   end function SuperLUSolution
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine distributeSparseMatrix()
!  ===================================================================
   use MPPModule , only : syncAllPEs
!
   implicit none
!
   integer (kind=IntKind), pointer :: rowptr(:), colind(:)
   integer (kind=IntKind) :: row_dim
   integer (kind=IntKind) :: nnz
!
   complex (kind=CmplxKind), pointer :: values(:)
!
!  ===================================================================
!  Bail out if MyNod does not belong in the grid.
!  ===================================================================
!
   pSLU => SLU
!
#ifdef DEBUG
   if ( MyPE>=(nprow*npcol) ) then
      write(6,*) " ERROR :: MyPE = ",MyPE
      write(6,*) "Stop program"
      return
   endif
   if ( print_level>0 .and. MyPE==0 ) then
      write(6,*) ' Process grid ', nprow, ' X ', npcol
   endif
#endif
!
!  ===================================================================
   nnz = pSLU%nnz
   row_dim = pSLU%n_row
   rowptr => pSLU%rowptr
   colind => pSLU%colind
   values => pSLU%nzval
!  ===================================================================
!  Distribute the matrix to the process gird
!  ===================================================================
   call  f_zcreate_dist_matrix( A, NRow_g, NCol_g, nnz, values,   &
                               colind, rowptr, grid, row_dim, FSTRow)
#ifdef DEBUG
   if (print_level>2)then
!      call f_zPrint_CompRowLoc_Matrix_dist(A)
   endif
#endif
!
!
   end subroutine distributeSparseMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine destroyDistributedMatrix()
!  ===================================================================
   implicit none
!
   call f_Destroy_CompRowLoc_Matrix_dist(A)
!
   end subroutine destroyDistributedMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSolutionSLU()
!  ===================================================================
   use WriteMatrixModule, only : writeMatrix
   use MPPModule, only : syncAllPEs
#ifdef TIMING
   use TimerModule, only : getTime
#endif
!
   implicit none
!
   integer (kind=IntKind) :: info, i
!
#ifdef TIMING
   real (kind=RealKind) :: time_0
#endif
!
!  ===================================================================
!  Setup the right hand side
!  ===================================================================
!  -------------------------------------------------------------------
   call  get_CompRowLoc_Matrix( A, nrow = row_dim, ncol = col_dim, &
                                nrow_loc = ldb)
!  -------------------------------------------------------------------
#ifdef DEBUG
   if ( print_level>0 ) then
      write(6,*) 'getInverse :: ldb, Nrhs',ldb,nrhs
   endif
#endif
!
!   if ( isInverse ) then
!     ----------------------------------------------------------------
!      call setRhstoUnit(ldb,Nrhs,cone)
!     ----------------------------------------------------------------
!   endif
!  ===================================================================
!  Initialize ScalePermstruct and LUstruct
!  ===================================================================
   call f_ScalePermstructInit(row_dim, col_dim, ScalePermstruct)
   call f_LUstructInit(row_dim, col_dim, LUstruct)
!  ===================================================================
!  Initialize the statistics variables
!  ===================================================================
   call f_PStatInit(stat)
!  ===================================================================
!  Call the linear equation solver
!  ===================================================================
#ifdef TIMING
   time_0 = getTime()
#endif
!  -------------------------------------------------------------------
   call f_pzgssvx(options, A, ScalePermstruct, B_cmplx, ldb, Nrhs,  &
                  grid, LUstruct, SOLVEstruct, berr, stat, info)
!  -------------------------------------------------------------------
#ifdef TIMING
   time_pzgssvx = getTime() - time_0 + time_pzgssvx
   count_t = count_t+1
#endif
!   call f_PStatPrint(options,stat,grid)
#ifdef DEBUG
   if ( info == 0 ) then
      if ( print_level>2 .and. MyPE==0 ) then
         write (*,*) 'Backward error: ', (berr(i), i = 1,Nrhs)
      endif
   else
      write(6,*) 'INFO from f_pzgssvx = ', info
   endif
#endif
!
!  ===================================================================
!  Deallocate SuperLU allocated storage
!  ===================================================================
!  -------------------------------------------------------------------
   call f_PStatFree(stat)
   call f_ScalePermstructFree(ScalePermstruct)
   call f_Destroy_LU(col_dim, grid, LUstruct)
   call f_LUstructFree(LUstruct)
   call f_zSolveFinalize(options, SOLVEstruct)
!
   end subroutine calSolutionSLU
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine setRhsToZero(dim1,dim2)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind),intent(in) :: dim1, dim2
!
   integer (kind=IntKind) :: i
!
   if ( .not.allocated(B_cmplx) ) then
      allocate(B_cmplx(dim1,dim2))
   else if ( dim1/=size(B_cmplx,dim=1)) then
      deallocate(B_cmplx)
      if (dim2/=Nrhs) then
         call ErrorHandler('setRhsToZero','Wrong dimensions',dim2,Nrhs)
      endif
      allocate( B_cmplx(dim1,dim2) )
   endif
!
   B_cmplx(1:dim1,1:dim2) = czero
!
   end subroutine setRhsToZero
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine setRhs( pb, prow, pcol, dim1, dim2, dim1_g, dim2_g)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind),intent(in) :: prow, pcol
   integer (kind=IntKind),intent(in) :: dim1, dim2
   integer (kind=IntKind),intent(in) :: dim1_g, dim2_g
!
   complex (kind=CmplxKind), target :: pb(:,:)
!
   character (len=6) :: sname = "setRhs"
!
   integer (kind=IntKind) :: i, j, status
!
#ifdef DEBUG
   if ( dim1>dim1_g .or. dim2>Nrhs ) then
      write(6,*) dim1_g, dim2_g, Nrhs
      call ErrorHandler( 'setRhs', 'Wrong block dimensions',    &
                         dim1, dim2 )
   endif
#endif
!
   if ( .not.allocated(B_cmplx) ) then
      allocate(B_cmplx(dim1_g,dim2_g), STAT=status)
#ifdef DEBUG
      call check_allocate( status, sname )
#endif
      do i = 1,dim2_g
         B_cmplx(:,i) = czero
      enddo
      isInverse = .false.
   else if ( dim1_g/=size(B_cmplx,dim=1)) then
      deallocate(B_cmplx)
#ifdef DEBUG
      if (dim2_g/=Nrhs) then
         call ErrorHandler('setRhs','Wrong dimensions',          &
                            dim2_g,Nrhs)
      endif
#endif
      allocate( B_cmplx(dim1_g,dim2_g), STAT=status )
#ifdef DEBUG
      call check_allocate( status, sname )
#endif
      do i = 1,dim2_g
         B_cmplx(:,i) = czero
      enddo
      isInverse = .false.
   endif
!
   do i = 1,dim2
#ifdef No_BLAS
      B_cmplx(prow+1:prow+dim1,pcol+i)=pb(1:dim1,i)
#else
      call zcopy(dim1,pb(1,i),1,B_cmplx(prow+1,pcol+i),1)
#endif
   enddo
!
   end subroutine setRhs
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getCompressedRowMatrix( size1, size2, max_size, pMtx,   &
                                      ind_row, ind_col, nnz_b )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: size1, size2, max_size
   integer (kind=IntKind), intent(in) :: ind_col, ind_row
   integer (kind=IntKind), intent(out):: nnz_b
!
   complex (kind=CmplxKind), target :: pMtx(:,:)
!
   character (len=22) :: sname = "getCompressedRowMatrix"
!
   integer (kind=IntKind) :: n_nz, r_nz, i, j, size, status
   integer (kind=IntKind), pointer :: pcol(:), prow(:)
!
   complex (kind=CmplxKind), pointer :: pnzval(:)
!
   size = max(size1,size2)
   if ( size>maxWorkSize .or. .not.isAllocWorkBlock ) then
      call setWorkBlock(size,1)
   endif
   prow => workRow
   pcol => workCol
   pnzval => workNNZ_z
!
   prow = 0
   pcol = 0
   pnzval = czero
!
   n_nz = 0
   r_nz = 0
   do i = 1,size1
      prow(i)= r_nz
      do j = 1,size2
         if ( abs( pMtx(i,j) ) > tol ) then
            n_nz = n_nz+1
            pnzval(n_nz) = pMtx(i,j)
            pcol(n_nz) = j-1 ! C compatible indexing
            r_nz = r_nz+1    ! in C array indecies start at 0
         endif
      enddo
   enddo
   prow(size1+1) = r_nz
!  ===================================================================
!  Save Matrix in the compressed row format
!  ===================================================================
!
   n_SLU_list = n_SLU_list + 1
!
   call resetSLUstruct(n_nz,size1,n_SLU_list)
!   allocate( pSLU_l%colind(n_nz), STAT=status )
!#ifdef DEBUG
!   call check_allocate( status, sname )
!#endif
!   allocate( pSLU_l%rowptr(size1+1), STAT=status )
!#ifdef DEBUG
!   call check_allocate( status, sname )
!#endif
!   allocate( pSLU_l%nzval(n_nz), STAT=status )
!#ifdef DEBUG
!   call check_allocate( status, sname )
!#endif
!
   pSLU_l => SLU_l(n_SLU_list)
#ifdef No_BLAS
   pSLU_l%nzval(1:n_nz) = pnzval(1:n_nz)
#else
   call zcopy(n_nz,pnzval(1),1,pSLU_l%nzval(1),1)
#endif
   pSLU_l%colind(1:n_nz) = pcol(1:n_nz)
   pSLU_l%rowptr(1:size1+1) = prow(1:size1+1)
!
   pSLU_l%ind_row = ind_row
   pSLU_l%ind_col = ind_col
   pSLU_l%nnz = n_nz
   pSLU_l%n_row = size1
   pSLU_l%n_col = size2
!
   nnz_b = n_nz
!
   nullify( pnzval, pcol, prow )
!
   end subroutine getCompressedRowMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setWorkBlock(size,type)
!  ===================================================================
   implicit none
!
   character (len=12) :: sname = "setWorkBlock"
!
   integer (kind=IntKind), intent(in) :: size, type
!
   integer (kind=IntKind) :: status
!
   if ( maxWorkSize< size .and. maxWorkSize>0) then
      call delWorkBlock()
   endif
   maxWorkSize = size
!
   allocate( workRow(size+1), workCol(size*size), STAT=status )
#ifdef DEBUG
   call check_allocate( status, sname )
#endif
   if ( type /= 0 ) then
      allocate( workNNZ_z(size*size), STAT=status )
#ifdef DEBUG
      call check_allocate( status, sname )
#endif
   endif
!
   isAllocWorkBlock = .true.
!
   end subroutine setWorkBlock
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delWorkBlock()
!  ===================================================================
   implicit none
!
   if ( .not.isAllocWorkBlock ) then
      return
   endif
   deallocate( workRow, workCol)
   if ( allocated( workNNZ_z) ) then
      deallocate( workNNZ_z )
   endif
   maxWorkSize = -1
   isAllocWorkBlock = .false.
!
   end subroutine delWorkBlock
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insertSLUlist(nnz, d_row, d_col, nlst )
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nnz  ! total number of nonzero
                                              ! elements in matrix
   integer (kind=IntKind), intent(in) :: d_row ! row dimension of matrix
   integer (kind=IntKind), intent(in) :: d_col ! column dimension of matrix
   integer (kind=IntKind), intent(in) :: nlst ! number of elements
                                              ! in the list
!
   integer (kind=IntKind) :: i, j, k, count_rowb, diff_rb, status
   integer (kind=IntKind) :: r_begin,r_end
   integer (kind=IntKind) :: rb_begin,rb_end
   integer (kind=IntKind) :: pos_b, dim_b
!
   type(SLUCompRowMtx_Cmplx), pointer :: ptmp
!
#ifdef DEBUG
   if ( nlst /= n_SLU_list ) then
      call ErrorHandler('insertSLUlist', 'Different number of blocks', &
                         nlst, n_SLU_list )
   endif
#endif
   n_max_list = max(n_SLU_list,n_max_list)
!
#ifdef DEBUG
   if ( print_level>1 ) then
      write(6,*)" Enter insertSLUlist_z"
   endif
#endif
!
   call resetSLUstruct(nnz,d_row)
!   allocate(SLU%colind(nnz), STAT=status)
!#ifdef DEBUG
!   call check_allocate(status)
!#endif
!   allocate(SLU%rowptr(d_row+1), STAT=status)
!#ifdef DEBUG
!   call check_allocate(status)
!#endif
!   allocate(SLU%nzval(nnz), STAT=status)
!#ifdef DEBUG
!   call check_allocate(status)
!#endif
!
   pSLU => SLU
   pSLU_l => SLU_l(1)
!
   pSLU%nnz = nnz
   pSLU%ind_row = 1
   pSLU%ind_col = 1
   pSLU%n_row = d_row
   pSLU%n_col = d_col
!
   r_begin = 0
   r_end = 0
   count_rowb = 0
   pSLU%rowptr = 0
   pSLU%colind = 0
   pSLU%nzval = zero
!
   if ( nlst == 1 ) then
      pSLU%colind(:) = pSLU_l%colind(:)+pSLU_l%ind_col-1
      pSLU%rowptr(:) = pSLU_l%rowptr(:)
      pSLU%nzval(:)  = pSLU_l%nzval(:)
      return
   endif
!
   do i = 1,d_row ! get the nnzval on the rows of the global matrix
!
      count_rowb = count_rowb+1
      pSLU%rowptr(i) = r_begin
      dim_b = 0
#ifdef DEBUG
      if ( print_level > 4 ) then
         write(6,*)'row: r_begin: count: ',   &
                   i,pSLU%rowptr(i),count_rowb
      endif
#endif
!
      ListLoop : do j = 1,nlst
!
         ptmp => SLU_l(j)
         pos_b = ptmp%ind_row
         if ( i==(pos_b+count_rowb-1) ) then
!
            dim_b = ptmp%n_row
            rb_begin = ptmp%rowptr(count_rowb)
            rb_end = ptmp%rowptr(count_rowb+1)
            diff_rb = rb_end-rb_begin
#ifdef DEBUG
            if ( print_level > 4 ) then
               write(6,*)"nlist: dim_b: rb_begin: rb_end",             &
                         j, dim_b, rb_begin, rb_end
            endif
#endif
            if ( diff_rb >= 1 ) then
               r_end = r_begin+diff_rb
               do k = 1,diff_rb
                  pSLU%colind(r_begin+k) =                             &
                        ptmp%colind(rb_begin+k)+ptmp%ind_col-1
                  pSLU%nzval(r_begin+k) = ptmp%nzval(rb_begin+k)
#ifdef DEBUG
                  if ( print_level > 4 ) then
                     write(6,*)"   cloind: nzval : ", &
                         pSLU%colind(r_begin+k),pSLU%nzval(r_begin+k)
                  endif
#endif
               enddo
               r_begin = r_end
            endif
!
         endif
!
      enddo ListLoop
!
      if (count_rowb==dim_b) then
         count_rowb=0
      endif
!
   enddo
!
   pSLU%rowptr(d_row+1) = r_begin
!
#ifdef DEBUG
   if ( pSLU%nnz /= pSLU%rowptr(d_row+1) ) then
      write(6,*) "Error insertSLUlist_c MyPe=",MyPE," :: nnz, rowptr", &
                 pSLU%nnz, pSLU%rowptr(d_row+1)
      call ErrorHandler( "insertSLUlist_z",                            &
                         "Wrong count number of non-zero elements",    &
                          pSLU%nnz, pSLU%rowptr(d_row+1) )
   endif
#endif
!
   end subroutine insertSLUlist
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine resetSLUstruct( n_nz, row_sz, i_list )
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n_nz, row_sz
   integer (kind=IntKind), optional :: i_list
!
   type (SLUCompRowMtx_Cmplx), pointer :: p_tmp
!
   if ( present(i_list) ) then
#ifdef DEBUG
      if ( i_list /= n_SLU_list ) then
         call ErrorHandler('resetSLUstruct', 'Different number of blocks', &
                         i_list, n_SLU_list )
      endif
#endif
      p_tmp => SLU_l(i_list)
      if ( .not.associated( p_tmp%colind ) ) then
         allocate( p_tmp%colind(n_nz) )
         allocate( p_tmp%nzval(n_nz) )
      else if ( n_nz > size( p_tmp%colind ) ) then
         deallocate( p_tmp%colind )
         deallocate( p_tmp%nzval )
         allocate( p_tmp%colind(n_nz) )
         allocate( p_tmp%nzval(n_nz) )
#ifdef DEBUG
         reallocate_count = reallocate_count+1
#endif
      endif
      if ( .not.associated( p_tmp%rowptr ) ) then
         allocate( p_tmp%rowptr(row_sz+1) )
      else if ( row_sz+1 > size( p_tmp%rowptr) ) then
         deallocate( p_tmp%rowptr )
         allocate( p_tmp%rowptr(row_sz+1) )
      endif
   else
      if ( .not.associated(SLU%colind) ) then
         allocate( SLU%colind(n_nz) )
         allocate( SLU%nzval(n_nz) )
      else if ( n_nz > size(SLU%colind) ) then
         deallocate( SLU%colind )
         deallocate( SLU%nzval )
         allocate( SLU%colind(n_nz) )
         allocate( SLU%nzval(n_nz) )
      endif
      if ( .not.associated( SLU%rowptr ) ) then
         allocate( SLU%rowptr(row_sz+1) )
      else if ( row_sz+1 > size( SLU%rowptr ) ) then
         deallocate( SLU%rowptr )
         allocate( SLU%rowptr(row_sz+1) )
      endif
   endif
!
   nullify(p_tmp)
!
   end subroutine resetSLUstruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delSLUstruct( n_list )
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), optional :: n_list
!
   integer (kind=IntKind) :: i
!
   type (SLUCompRowMtx_Cmplx), pointer :: p_tmp
!
   if ( present(n_list) ) then
#ifdef DEBUG
      if ( n_list /= n_SLU_list ) then
         call ErrorHandler('delSLUstruct', 'Different number of blocks', &
                         n_list, n_SLU_list )
      endif
#endif
      do i = 1,n_list
         p_tmp => SLU_l(i)
         deallocate(p_tmp%colind)
         deallocate(p_tmp%rowptr)
         deallocate(p_tmp%nzval)
      enddo
      n_SLU_list = 0
   else
      deallocate(SLU%colind)
      deallocate(SLU%rowptr)
      deallocate(SLU%nzval)
   endif
!
   nullify( p_tmp )
!
   end subroutine delSLUstruct
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeSLU2Matrix(fst_row, iprint)
!  ===================================================================
   use WriteMatrixModule, only : writeMatrix
   use MathParamModule, only : ten2m8, czero
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: fst_row, iprint
!
   type (SLUCompRowMtx_Cmplx), pointer :: pSLU
!
   integer (kind=IntKind)   :: i,j,n_ind,row_nnz,row_g,col_g,idvar
   integer (kind=IntKind)   :: n_row, n_col
!
   complex (kind=CmplxKind) :: nzval
   complex (kind=CmplxKind), allocatable :: matrix(:,:)
!
   pSLU => SLU
   n_row = pSLU%n_row
   n_col = pSlu%n_col
!
   allocate( matrix(n_row,n_col) )
   matrix = czero
!
   n_ind = 1
   i = 1
   do i=1,pSLU%n_row
      row_nnz = pSLU%rowptr(i+1)-pSLU%rowptr(i)
      row_g = i+fst_row               ! fortran indexing
      do j=1,row_nnz
         col_g = pSLU%colind(n_ind) +1   ! fortran indexing
         nzval = pSLU%nzval(n_ind)
         matrix(i,col_g) = nzval
         n_ind = n_ind+1
      enddo
   enddo
!
   if ( iprint >= 0 ) then
      call writeMatrix( 'SLU matrix', matrix, n_row, n_col )
      call writeMatrix( 'Rhds matrix', B_cmplx, n_row, Nrhs )
   endif
!
   if ( print_level > 4 ) then
      call writeBinaryMatrix( matrix, n_row, n_col, fst_row )
   endif
   if ( print_level >= 2 ) then
     print_level = 0
   endif
!
   deallocate( matrix )
!
   end subroutine writeSLU2Matrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSolution()                                result(pB)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), pointer :: pB(:,:)
!
   pB => B_Cmplx
!
   end function getSolution
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSuperLUSolver()
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind) :: i
!  ===================================================================
!  Release the SuperLU process grid
!  ===================================================================
!  -------------------------------------------------------------------
   call f_superlu_gridexit(grid)
!  -------------------------------------------------------------------
!  ===================================================================
!  Destroy C structures in superlu_matrix_type
!  ===================================================================
!  -------------------------------------------------------------------
   call f_destroy_gridinfo_handle(grid)
   call f_destroy_options_handle(options)
   call f_destroy_ScalePermstruct_handle(ScalePermstruct)
   call f_destroy_LUstruct_handle(LUstruct)
   call f_destroy_SOLVEstruct_handle(SOLVEstruct)
   call f_destroy_SuperMatrix_handle(A)
   call f_destroy_SuperLUStat_handle(stat)
!  -------------------------------------------------------------------
!
   if ( isAllocWorkBlock ) then
      call delWorkBlock()
   endif
!
   if ( allocated(B_cmplx) ) then
      deallocate(B_cmplx)
   endif
   if ( allocated(berr) ) then
      deallocate(berr)
   endif
!
   call delSLUstruct(n_SLU_list)
   deallocate(SLU_l)
   call delSLUstruct()
!
#ifdef DEBUG
   write(6,'(a,i16)') "SuperLU:: Total reallocations of SLU structures",&
                      reallocate_count
#endif
#ifdef TIMING
   write(6,'(a,d15.8)') "SuperLU:: Total time in pzgssvx: ", &
                        time_pzgssvx
   write(6,'(a,d15.8)') "          Total time in SuperLU: ", &
                        time_superLu
#endif
!
#ifdef DEBUG
   reallocate_count = 0
#endif
   isAllocWorkBlock = .false.
   isDiag = .false.
   Initialized = .false.
!
   end subroutine endSuperLuSolver
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine check_allocate( i, sname, size )
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
!
   character (len=*), optional :: sname
!
   integer (kind=IntKind), optional :: size
!
   if ( i==0 ) then
      return
   else
      if( present(sname) ) write(6,*) "Allocation in ", sname
      if( present(size) ) write(6,*) "  Size : ", size
      call ErrorHandler('check_allocate','Problems with allocation', i)
   endif
!
   end subroutine check_allocate
!  ===================================================================
!
end module SuperLUSolverModule
