module MatrixBlockInversionModule
!
use KindParamModule, only : IntKind, RealKind, CmplxKind
use MathParamModule, only : cone, czero, zero, one, &
                            ten2m8, ten2m6, ten2m10, ten2m20
use ErrorHandlerModule, only : ErrorHandler
!
public :: initMatrixBlockInv, &
          endMatrixBlockInv,  &
          invertMatrixBlock
!
private
   logical :: Initialized
   logical :: isLU = .true.
   logical :: isSymmetryOn = .false.
!
   integer (kind=IntKind), allocatable :: print_level(:)
!
   integer (kind=IntKind) :: Algorithm
   integer (kind=IntKind) :: NumMatrixes
   integer (kind=IntKind) :: NumSteps
   integer (kind=IntKind) :: MaxNumBlocks
   integer (kind=IntKind) :: MaxNumBlockSize
   integer (kind=IntKind) :: MinBlockSize = 16
!
   integer (kind=IntKind), allocatable :: AlgSwitch(:,:)
   integer (kind=IntKind), allocatable :: CountSteps(:)
   integer (kind=IntKind), allocatable :: NumKBlocks(:)
   integer (kind=IntKind), allocatable :: KBlocksSize(:,:)
   integer (kind=IntKind), allocatable :: NumBlocks(:)
   integer (kind=IntKind), allocatable :: MaxMatrixSize(:)     
   integer (kind=IntKind), allocatable :: Block_Size(:)
   integer (kind=IntKind), allocatable :: Block_Dir(:)
   integer (kind=IntKind), allocatable, target :: MatrixBlockSizes(:,:)
!
   integer (kind=IntKind), allocatable, target :: ipvt(:)
   integer (kind=IntKind), allocatable, target :: idcol(:,:)
!
   real (kind=RealKind), allocatable :: Time_Old(:),Time_Dir(:)
!
   integer (kind=IntKind), allocatable, target :: iwspace(:)
   real (kind=RealKind), allocatable, target :: rwspace(:)
   complex (kind=CmplxKind), allocatable, target :: cwspace(:) 
!
   complex (kind=CmplxKind), allocatable, target :: symop_wspace(:)
   complex (kind=CmplxKind), allocatable, target :: vecspace(:)  
!
contains 
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMatrixBlockInv( alg, nmtrx, nsteps, nblck, max_nblck, &
                                  k_blocks, iprint )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: alg, nmtrx, nsteps, max_nblck
   integer (kind=IntKind), intent(in) :: nblck(nmtrx)
   integer (kind=IntKind), intent(in) :: k_blocks(max_nblck,nmtrx)
   integer (kind=IntKind), intent(in) :: iprint(nmtrx)
!
   integer (kind=IntKind) :: i, j, nblk, max_nblk, min_sz
   integer (kind=IntKind) :: max_block, max_matrixsz, max_vsz
!
   if (Initialized) then
      return
   endif
!
   Algorithm = alg
   NumMatrixes = nmtrx
   NumSteps = nsteps+1
   MaxNumBlocks = max_nblck
!
   allocate( NumBlocks(NumMatrixes) )
   allocate( MaxMatrixSize(NumMatrixes) )
   allocate( AlgSwitch(NumSteps,NumMatrixes) )
   allocate( CountSteps(NumMatrixes) )
   allocate( NumKBlocks(NumMatrixes) ) 
   allocate( KBlocksSize(MaxNumBlocks,NumMatrixes) )
   allocate( Block_Size(NumMatrixes) )
   allocate( Block_Dir(NumMatrixes) )
   allocate( Time_Old(NumMatrixes) )
   allocate( Time_Dir(NumMatrixes) )
   allocate( MatrixBlockSizes(MaxNumBlocks,NumMatrixes) );  MatrixBlockSizes = 0
   allocate( print_level(NumMatrixes) )
!
   do i = 1,NumMatrixes
      NumKBlocks(i) = nblck(i) 
   enddo
!
#ifdef ACCEL
   do i = 1,NumMatrixes
      MaxMatrixSize(i) = maxval(k_blocks(:,i))
   enddo
   print *,'Calling initialize_matinv_accel...........................'
!  -------------------------------------------------------------------
   call initialize_matinv_accel(NumMatrixes,NumKBlocks,MaxMatrixSize)
!  -------------------------------------------------------------------
#endif
!
   CountSteps(1:NumMatrixes) = 0
   AlgSwitch(1:NumSteps,1:NumMatrixes) = alg
   Block_Size(1:NumMatrixes) = 3
   Block_Dir(1:NumMatrixes) = 1
   Time_Dir(1:NumMatrixes) = 1.d10
   Time_Old(1:NumMatrixes) = ZERO
   print_level(1:NumMatrixes) = iprint(1:NumMatrixes)
!
   MaxMatrixSize(1:NumMatrixes) = 0
   MaxNumBlockSize = 0
   do i = 1,NumMatrixes
      AlgSwitch(1,i) = 2
      do j = 1,nblck(i)
         MaxMatrixSize(i) = MaxMatrixSize(i) + k_blocks(j,i)
         KBlocksSize(j,i) = k_blocks(j,i)
      enddo
      MaxNumBlockSize = max(MaxNumBlockSize,MaxMatrixSize(i))
   enddo
!
!  ===================================================================
!  for a given algorithm setup the blocks of the big matrix
!  ===================================================================
!
   if ( Algorithm >= 10) then
      nblk = 2
      NumBlocks(1:NumMatrixes) = 2
!     allocate( MatrixBlockSizes(nblk,NumMatrixes) )
      do i = 1,NumMatrixes
         MatrixBlockSizes(1,i) = KBlocksSize(1,i)
!        MatrixBlockSizes(2,i) = 0 
         MatrixBlockSizes(2,i) = MaxMatrixSize(i)-MatrixBlockSizes(1,i)
      enddo
   else
!     ================================================================
!     The minimum block size is 16
!     ================================================================
      if ( Algorithm == 1) then
!         allocate( MatrixBlockSizes(MaxNumBlocks,NumMatrixes) )
         do i = 1,NumMatrixes
            MatrixBlockSizes(1,i) = KBlocksSize(1,i)
            do j = 2,nblck(i)
               MatrixBlockSizes(j,i) = KBlocksSize(j,i)
            enddo
!           NumBlocks(i) = NumKBlocks(i)
         enddo
         NumBlocks(1:NumMatrixes) = NumKBlocks(1:NumMatrixes)
      else
         do i = 1,NumMatrixes                   ! debugged
            if ( Block_Size(i) == 0 ) then
               Block_Size(i)=1
            endif
         enddo                                  ! debugged
         max_nblk = 0
         do i = 1,NumMatrixes
            nblk = MaxMatrixSize(i)-KBlocksSize(1,i)
!            NumBlocks(i) = min( Block_Size(i),nblk/MinBlockSize )+1
            NumBlocks(i) = min( Block_Size(i),nblk/KBlocksSize(1,i) )+1
         enddo
         max_nblk = maxval( NumBlocks )
!        allocate( MatrixBlockSizes(max_nblk,NumMatrixes) )
         do i = 1,NumMatrixes
            MatrixBlockSizes(1,i) = KBlocksSize(1,i)
            if ( NumBlocks(i) > 1 ) then
               min_sz = (MaxMatrixSize(i)-MatrixBlockSizes(1,i))/      & 
                                                (NumBlocks(i)-1)
!              do j = 2,NumBlocks(i)-1
!                 MatrixBlockSizes(j,i) = min_sz
!              enddo
               MatrixBlockSizes(2:NumBlocks(i)-1,i) = min_sz
               MatrixBlockSizes(NumBlocks(i),i) = MaxMatrixSize(i) -   &
                         MatrixBlockSizes(1,i)-(NumBlocks(i)-2)*min_sz
            endif
         enddo
      endif
   endif   
!
!  Allocate the working space
!
   max_block = 0
   do i = 1,NumMatrixes
      max_block = max(max_block,KBlocksSize(1,i) )
   enddo
   max_matrixsz = maxval( MaxMatrixSize )
!
   allocate( iwspace(max_block) )
   allocate( rwspace(max_block) )
   allocate( cwspace(max_block) )
!
   allocate( ipvt(max_matrixsz) )
   allocate( idcol(max_block,NumMatrixes) )
   allocate( symop_wspace(max_block*max_block*(max_block-1)/2) )
   do i = 1,NumMatrixes
      idcol(1,i) = 0
   enddo
!
   if ( Algorithm>=3 .and. Algorithm<=6 ) then
      max_vsz = 0
      do i = 1,NumMatrixes
         max_vsz = max(max_vsz,MaxMatrixSize(i)-KBlocksSize(1,i))
      enddo
      allocate( vecspace(max_vsz*(max_block*6+6)) )
   endif
!
   Initialized = .true.
!
   end subroutine initMatrixBlockInv
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endMatrixBlockInv()
!  ===================================================================
   implicit none
!
   deallocate( MatrixBlockSizes, NumBlocks, MaxMatrixSize, AlgSwitch )
   deallocate( NumKBlocks, KBlocksSize, CountSteps )
   deallocate( Block_Size, Block_Dir )
   deallocate( Time_Old, Time_Dir )
!
   deallocate( iwspace, ipvt, idcol )
   deallocate( rwspace )
   deallocate( symop_wspace, cwspace, print_level  )
!
   if ( Algorithm>=3 .and. Algorithm<=6 ) then
      deallocate(vecspace)
   endif
!
#ifdef ACCEL
   print *,'Calling finalize_matinv_accel...........................'
   call finalize_matinv_accel()
#endif
!
   Initialized = .false.
!
   end subroutine endMatrixBlockInv
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setBlockAlgorithm(id,alg)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: id, alg
!
   integer(kind=IntKind) :: dummy, blockSize
!
   if ( id<0 .or. id>NumMatrixes ) then
      call ErrorHandler("setBlockAlgorithm",'Invalid id',id)
   endif
   if ( .not.Initialized ) then
      call ErrorHandler('setBlockAlgorithm','initMatrixBlockInv have')
   endif
!
   if ( alg >= 10) then
      NumBlocks(id) = 2
      if ( NumBlocks(id)> MaxNumBlockSize ) then
         call ErrorHandler('setBlockAlgorithm','wrong algorithm switch')
      endif
      MatrixBlockSizes(1,id) = KBlocksSize(1,id)
!     MatrixBlockSizes(2,id) = 0 
      MatrixBlockSizes(2,id) = MaxMatrixSize(id)-MatrixBlockSizes(1,id)
   else
!     ================================================================
!     The minimum block size is 16
!     ================================================================
      if ( alg == 1) then
         NumBlocks(id) = NumKBlocks(id)
         if ( NumBlocks(id)> MaxNumBlockSize ) then
            call ErrorHandler('setBlockAlgorithm','wrong algorithm switch')
         endif
         MatrixBlockSizes(1:NumKBlocks(id),id) = KBlocksSize(1:NumKBlocks(id),id)
!        do j = 2,NumKBlocks(id)
!           MatrixBlockSizes(j,id) = KBlocksSize(j,id)
!        enddo
      else
         blockSize = Block_Size(id)
         if ( blockSize == 0 ) then
            blockSize = 1
         endif
         dummy = MaxMatrixSize(id)-KBlocksSize(1,id)
!         NumBlocks(id) = min( blockSize,dummy/MinBlockSize )+1
         NumBlocks(id) = min( blockSize,dummy/KBlocksSize(1,id) )+1
         if ( NumBlocks(id)>MaxNumBlockSize ) then
            call ErrorHandler('setBlockAlgorithm','wrong algorithm switch')
         endif
!
         MatrixBlockSizes(1,id) = KBlocksSize(1,id)
         if ( NumBlocks(id) > 1 ) then
            dummy = (MaxMatrixSize(id)-MatrixBlockSizes(1,id))/      & 
                                                  (NumBlocks(id)-1)
            MatrixBlockSizes(2:NumBlocks(id)-1,id) = dummy
!           do j = 2,NumBlocks(id)-1
!              MatrixBlockSizes(j,id) = dummy
!           enddo
            MatrixBlockSizes(NumBlocks(id),id) = MaxMatrixSize(id) -  &
                      MatrixBlockSizes(1,id)-(NumBlocks(id)-2)*dummy
         endif
      endif
   endif   
!
   AlgSwitch(CountSteps(id),id) = alg
!  
   end subroutine setBlockAlgorithm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine invertMatrixBlock( id_in, b, ldb, nb, a, lda, na )
!  ===================================================================
!
!  *******************************************************************
!  PURPOSE:   inverse the first block of a complex matrix: a
!             (a^{-1})_00=(a_00-b)^{-1}
!
!  INPUT:     a,      the complex matrix to be inverted
!                    [ a11, a12 ]                [ b11, b12 ]
!                a = [          ], b = a**{-1} = [          ]       
!                    [ a21, a22 ]                [ b21, b22 ]
!
!             alg,        the algorithm of the invertion
!                         = 1, UL method
!                         = 2, LU method
!                         = 3,4,5, W.A.S. method
!
!  OUTPUT:    a,   contains the first block of the inverted matrix
!             b = a12 * a22**{-1} * a21 = a11 - b11**{-1}
!             Therefore, b11 = [a11 - b]**{-1}
!  *******************************************************************
   use TimerModule, only : getTime
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id_in
   integer (kind=IntKind), intent(in) :: lda, na, ldb, nb
!
   complex (kind=CmplxKind), target :: a(lda,na)
   complex (kind=CmplxKind), target :: b(ldb,nb)
!
   character(len=3) :: alg_name
!
   integer (kind=IntKind) :: nlim, alg, id
   integer (kind=IntKind) :: i, j, k, ioff, iqmr
   integer (kind=Intkind) :: nblk, blk1, mp
   integer (kind=IntKind) :: blockDir, blockSize
   integer (kind=IntKind), pointer :: p_blksz(:)
   integer (kind=IntKind), pointer :: p_idcol(:)
!
   real (kind=RealKind) :: time, timeDirect, timeNew, timeOld, timeQmr
   real (kind=RealKind) :: tol
!
   complex (kind=CmplxKind), pointer :: vecs(:,:)
   complex (kind=CmplxKind), pointer :: pvecs(:,:)
   complex (kind=CmplxKind), pointer :: sym_ops(:,:,:)
!
   id = id_in
!
   if ( CountSteps(id)>=NumSteps ) then
!     for spin polarized case 
      CountSteps(id) = 0
   endif
!
   tol = ten2m6
   CountSteps(id) = CountSteps(id)+1
   timeDirect = Time_Dir(id)
   timeOld = Time_Old(id)
   timeNew    = timeDirect
   blockDir  = Block_Dir(id)
   blockSize = Block_Size(id)
!   
   time = getTime()
!  ===================================================================
!  Invert the KKR-Matrix using the method defined by 'alg'
!  if using QMR and either the method fails or it takes longer than
!  the direct method (LU or UL) then the routine uses LU thereafter
!  ===================================================================
!
   if ( id<1 .or. id>NumMatrixes ) then
      call ErrorHandler('invertBlockMatrix','Invalid id',id)
   endif
!
   nlim = 20
   iqmr = 0
!
   alg = AlgSwitch(CountSteps(id),id)
   blk1 = MatrixBlockSizes(1,id)
   p_blksz => MatrixBlockSizes(:,id)
   nblk = NumBlocks(id)
   p_idcol => idcol(1:blk1,id)
   p_idcol(1:blk1) = 0 
!
!  do i = 1,ldb
!     do j = 1,nb
!        a(j,i)=czero
!     enddo
!  enddo
   a(1:nb,1:ldb) = czero
   k = 0
!
   if ( alg<=2 .or. alg>=10 ) then
!     ================================================================
!     Use the  LU algorithm
!     ================================================================
      isLU = .true.
#ifdef ACCEL
      mp = maxval( MaxMatrixSize )
!     ----------------------------------------------------------------
      call zblock_lu_accel( a, lda, na, MatrixBlockSizes(1,id),       &
                            MaxNumBlocks, NumBlocks(id), p_idcol, blk1, ipvt, mp, k)
!     ----------------------------------------------------------------
#else
!     ----------------------------------------------------------------
      call zblock_lu1( id, a, lda, na, k)
!     ----------------------------------------------------------------
#endif
   else if ( alg>=3 .and. alg<=6 ) then
!      call ErrorHandler('block_inv','Error: not defined yet.',alg)
!      return
!     ================================================================
!     Use the QMR algorithm
!     ================================================================
      isLU = .false.
      call zblock_lu1( id, a, lda, na, k )
      vecs => gettarget_2c(na-blk1,blk1*6+6,vecspace)
      pvecs => vecs(1:na-blk1,blk1*6+1:blk1*6+6)
!     ----------------------------------------------------------------
      call wasinv( na-blk1,k,a(blk1+1,blk1+1),lda,   &
                   a(blk1+1,blk1-k+1),lda,vecs,            &
                   nlim,pvecs,                                         &
                   rwspace,cwspace,iwspace,tol,alg-2,iqmr)
!     ================================================================
!     if qmr algorithm fails then redo the current energy: Use LU..
!     ================================================================
      if ( iqmr>0 ) then
         isLU = .true.
         call zblock_lu1( id, a, lda, na, k )
      else
         call zgemm( 'n', 'n', blk1, k, na-blk1, -cone, a(1,blk1+1),   &
                      lda, a(blk1+1,blk1-k+1), lda, czero, a, lda )
      endif
   else
      call ErrorHandler('block_inv','incorect alg value',alg)
   endif
!
   b(1:blk1,1:blk1) = -a(1:blk1,1:blk1)
   sym_ops => gettarget_3c( blk1, blk1, (blk1-1)/2, symop_wspace )
!
   if ( p_idcol(1) == 0 ) then
!     do i = 1,blk1
!        do j = 1,blk1
!           b(j,i) = -a(j,i)
!        enddo
!     enddo
!      b(1:blk1,1:blk1) = -a(1:blk1,1:blk1)
!     ================================================================
!     Here look for symetry
!     ================================================================
!     if ( print_level(id) >= 0 ) then
!         write(6,'(''block_inv:: looking for sym'')')
!     endif
      call find_sym( id, b, blk1 )
   else
      k = 0
      ioff = 0
      do i = 1,blk1
        if ( p_idcol(i)==i ) then
            k = k+1
!           do j = 1,blk1
!              b(j,i) = -a(j,k)
!           enddo
            b(1:blk1,i) = -a(1:blk1,k)
         else
            j = p_idcol(i)
            ioff = ioff+1
            call zgemv( 'n', blk1, blk1, CONE, sym_ops(1,1,ioff),     &
                        blk1, b(1,j), 1, CZERO, b(1,i), 1 )
         endif
      enddo
   endif
!
   time = getTime() - time
   if ( alg<=2 .or. alg>=10 ) then
      alg_name = ' LU'
      if ( time>ZERO ) then
        timeNew = time
      endif
   else 
      timeQmr = ZERO
      if ( alg==6 ) then
         alg_name = 'MAR'
         if ( time>ZERO ) then
            timeQmr = time
         endif
      else
         alg_name = 'QMR'
         if ( time>ZERO ) then
            timeQmr = time
         endif
      endif
!     ================================================================
!     If qmr algorithm takes longer (or fails) than direct method.....
!     switch to LU for subsequent energies............................
!     ================================================================
      if ( timeQmr>=timeNew ) then
         do i = CountSteps(id),NumSteps
            AlgSwitch(i,id) = 2
         enddo
      endif
   endif

   if ( print_level(id) >= 0 .and. nblk > 1) then
      write( 6,'(/,a,a,a,i6,a7,f12.6,a4)') &
        'invertBlockMatrix     :: Using ',alg_name,': Block size = ', &
        MatrixBlockSizes(2,id),', Time=',time,' Sec' 
      call FlushFile(6)
   endif
!
   if ( timeDirect>=1.d10 ) then
      timeDirect = timeNew
      blockSize = blockSize + blockDir
   else if ( timeOld==ZERO) then
!     ================================================================
!     Make sure we are using LU at this energy and at the previous energy
!     ================================================================
      if ( CountSteps(id) >= 2 ) then
         if ( AlgSwitch(CountSteps(id)-1,id)==2 .and. &
              AlgSwitch(CountSteps(id),id)==2 ) then
            if ( timeNew<=timeDirect ) then
               timeOld = timeDirect
               timeDirect = timeNew
               blockSize = blockSize + blockDir
            else
               timeOld = timeNew
               blockDir = -blockDir
               blockSize = blockSize + 2*blockDir
            endif
         endif
      else
         if ( AlgSwitch(CountSteps(id),id)==2 ) then
            if ( timeNew<=timeDirect ) then
               timeOld = timeDirect
               timeDirect = timeNew
               blockSize = blockSize + blockDir
            else
               timeOld = timeNew
               blockDir = -blockDir
               blockSize = blockSize + 2*blockDir
            endif
         endif
      endif
   else if ( timeDirect<timeNew ) then
      blockSize = blockSize - blockDir
      blockDir = 0
   else if ( blockDir/=0 ) then
      timeOld = timeDirect
      timeDirect = timeNew
      blockSize = blockSize + blockDir
   endif
!
   if (blockSize <= 0) then      ! debugged
      blockSize = 1              ! debugged
   endif                         ! debugged
!
   Block_Dir(id) = blockDir
   Block_Size(id) = blockSize
   Time_Dir(id) = timeDirect
   Time_Old(id) = timeOld
!
   if (Algorithm /= 1) then
!     ----------------------------------------------------------------
      call setBlockAlgorithm( id,AlgSwitch(CountSteps(id),id) )
!     ----------------------------------------------------------------
   endif
!
   end subroutine invertMatrixBlock
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine zblock_lu1( id, a, lda, na, k )
!  ===================================================================
!  does a partitioning of a matrix and a partial inversion to
!  get the upper subblock of the inverse
!
!  a      -- complex*16 of (lda,size) the matrix to be inverted where
!            size is the sum of the all elements of blk_sz on return,
!            the upper subblock of a is untouched, and postprocessing
!            is needed to obtain the inverse
!
!  blk_sz -- integer of (nblk), the size of each subblock of matrix a
!
!  nblk   -- number of blocks (number of elements in blk_sz)
!
!  ipvt   -- integer of (mp), work array
!
!  idcol  -- integer of (blk_sz(1)), if idcol(i)=idcol(j), then
!            the two columns are equivalent by symmetry, and only
!            the min(i,j) th column of the inverse is calculated.
!
!  k      -- returns actual number of columns in the calculated inverse
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) ::  id, lda, na
   integer (kind=IntKind), intent(out) :: k
!
   complex (kind=CmplxKind), target :: a(lda,na)
!
   integer (kind=IntKind) :: i, m, n, info, blk1, nblk
   integer (kind=IntKind) :: ioff, joff, iblk
   integer (kind=IntKind), pointer :: p_ipvt(:)
   integer (kind=IntKind), pointer :: p_idcol(:)
!
!  print *, "start zblock_lu1"
!  ===================================================================
!  eliminate columns that are equiv due to symmetry
!  ===================================================================
!
   blk1 = MatrixBlockSizes(1,id)
   nblk = NumBlocks(id)
   p_idcol => idcol(1:blk1,id)
   k = blk1+1
!  print *,'k = ',k,', blk1 = ',blk1,', lda = ',lda,', na = ',na,', nblk = ',nblk
   do i = blk1,1,-1
      if ( p_idcol(1)==0 .or. p_idcol(i)==i ) then
         k = k-1
         if ( k/=i ) then
            call zcopy( na-blk1, a(blk1+1,i), 1, a(blk1+1,k), 1 )
         endif
      endif
   enddo
!
!  print *, 'nblk = ',nblk
   if ( isLU ) then
!     ================================================================
!     Do block LU
!     ================================================================
      n = MatrixBlockSizes(nblk,id)
      joff = na-n
      do iblk = nblk,2,-1
         m = n
         ioff = joff
         n = MatrixBlockSizes(iblk-1,id)
         joff = joff-n
!        =============================================================
!        invert the diagonal blk_sz(iblk) x blk_sz(iblk) block
!        =============================================================
         p_ipvt => ipvt(1:m)
!        -------------------------------------------------------------
         call zgetrf( m,m,a(ioff+1,ioff+1),lda,p_ipvt,info )
         if (info /= 0) then
            call ErrorHandler('BinvMatrix','zgetrf failed', info)
         endif
!        -------------------------------------------------------------
!        =============================================================
!        calculate the inverse of above multiplying the row block
!        blk_sz(iblk) x ioff
!        =============================================================
!        -------------------------------------------------------------
         call zgetrs( 'n',m,ioff,a(ioff+1,ioff+1),lda,p_ipvt,          &
                       a(ioff+1,1),lda,info )
         if (info /= 0) then
            call ErrorHandler('BinvMatrix','zgetrs failed', info)
         endif
         if (iblk.gt.2) then
!           ----------------------------------------------------------
            call zgemm( 'n','n',n,ioff-k+1,na-ioff,-cone,              &
                        a(joff+1,ioff+1),lda,a(ioff+1,k),lda,          &
                        cone,a(joff+1,k),lda )
!           ----------------------------------------------------------
            call zgemm( 'n','n',joff,n,na-ioff,-cone,a(1,ioff+1),lda,  &
                        a(ioff+1,joff+1),lda,cone,a(1,joff+1),lda )
!           ----------------------------------------------------------
         endif
      enddo
!     ----------------------------------------------------------------
      call zgemm( 'n', 'n', blk1, blk1-k+1, na-blk1, -cone,            &
                a(1,blk1+1), lda, a(blk1+1,k), lda, cone, a, lda )
!     ----------------------------------------------------------------
   endif
!
   k = blk1-k+1
!
   end subroutine zblock_lu1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine wasinv( NLEN,KKRSZ,AWAS,ndim,rhs,ldr,VECS,NLIM,         &
                      work,r0,rho,pointer_ind,TOL,ALG,IQMR)
!  ===================================================================
!
!  Purpose:
!
!  Parameters:
!  NDIM  = the declared size of A (input).
!  NLEN  = the size of the matrix (input).
!  KKRSZ =
!  A     = the matrix, NLEN by NLEN (input).
!  rhs   = the right hand sides, nlen by kkrsz
!  VECS  =
!  NLIM  = the maximum number of iterations (input).
!  work  = work space, nlen by 10
!  TOL   = the relative convergence tolerance (input).
!  ALG   = the choice of algorithm (input).
!  IQMR  = indicates an error and will switch to an LU based algorithm
!
!     External routines used:
!
!     Noel M. Nachtigal
!     March 4, 1994
!
!  *******************************************************************
!
   integer (kind=IntKind), intent(in) ::   ALG
   integer (kind=IntKind), intent(in) ::   KKRSZ
   integer (kind=IntKind), intent(in) ::   NDIM,ldr
   integer (kind=IntKind), intent(in) ::   NLEN
   integer (kind=IntKind), intent(in) ::   NLIM
   integer (kind=IntKind), intent(inout) ::   pointer_ind(kkrsz)
!
   complex (kind=CmplxKind), intent(inout), target :: AWAS(NDIM,nlen)
   complex (kind=CmplxKind), intent(inout), target :: rhs(ldr,kkrsz)
   complex (kind=CmplxKind), intent(inout), target :: VECS(nlen,kkrsz,6)
   complex (kind=CmplxKind), intent(inout), target :: work(nlen,6)
!
   real (kind=RealKind), intent(in) ::     TOL
!
!     Variables used by reverse communication.
!
   integer (kind=IntKind) ::  CB
   integer (kind=IntKind) ::  CX,cx0
   integer (kind=IntKind) ::  IERR
   integer (kind=IntKind) ::  REVCOM
   integer (kind=IntKind) ::  IQMR
   integer (kind=IntKind) ::  n
!
   real (kind=RealKind) ::   r0(kkrsz)
!
   complex (kind=CmplxKind) :: rho(kkrsz)
!
!     Local variables.
!
   integer (kind=IntKind) ::  I, J
   integer (kind=IntKind) ::  INFO(4)
   integer (kind=IntKind) ::  cols, cols_old
!
!  ===================================================================
   if ( kkrsz<=0 ) then
      call ErrorHandler("wasinv","kkrsz incorrect",kkrsz)
   endif
!  ===================================================================
!  if alg<=3 then there is no block-code.
!  ===================================================================
   if ( alg<=3 ) then 
!  -------------------------------------------------------------------
      call WASINV_p(NLEN,KKRSZ,AWAS,ndim,rhs,ldr,VECS,NLIM,           &
                      vecs(1,1,2),TOL,ALG,IQMR)
!  -------------------------------------------------------------------
      return
   endif
!  ===================================================================
!     Compute the modified right hand side.
!
!  ===================================================================
   do j = 1,kkrsz
      call zcopy(nlen,rhs(1,j),1,vecs(1,j,2),1)
   enddo
!  -------------------------------------------------------------------
   call zgemm('n','N',nlen,kkrsz,NLEN,-CONE,AWAS,NDIM,rhs,ldr,        &
              CONE,VECS(1,1,2),nlen)
!  -------------------------------------------------------------------
!  ===================================================================
!  precondition: note that diagonal is 1
!  ===================================================================
!  -------------------------------------------------------------------
   call ztrtrs('l','n','u',nlen,kkrsz,awas,ndim,vecs(1,1,2),nlen,     &
               info)
!  -------------------------------------------------------------------
!
   cols = kkrsz
   cx = 2
   do j = 1,kkrsz
      pointer_ind(j)=j
   enddo

   Loop_n: do n = 1,nlim
!
      cols_old = cols
      cols = 0
      cx0 = cx
!     ================================================================
!     Loop over all the columns.
!     ================================================================
      do J = 1,cols_old
!
!        Compute Jth column of the inverse.
!
         if ( n==1 ) then
!           ==========================================================
!           Set up call to linear systems solver.
!           Compute true residual norms, 
!           generate second starting vector.
!           ==========================================================
            INFO(1) = 0
            INFO(2) = 0
            call zcopy(nlen,vecs(1,j,2),1,work(1,2),1)
         else
            INFO(1) = 8
            INFO(2) = 1
            do i = 1,5
               call zcopy(nlen,vecs(1,j,i),1,work(1,i),1)
            enddo
!           ==========================================================
!           no preconditioner so we fake it
!           ==========================================================
           call zcopy(nlen,work(1,cx0),1,work(1,6),1)
         endif
!        =============================================================
!        call the solver.
!        =============================================================
         call zmar1( nlen,NLEN,work(1,1),TOL,INFO,                    &
                     n,rho(j),r0(j) )
!        =============================================================
!        Perform matrix-vector multiplications when needed.
!        =============================================================
         IERR   = INFO(1)
         REVCOM = INFO(2)

         if ( REVCOM==1 ) then
!           
            CX = INFO(3)
            CB = INFO(4)
            cols = cols + 1
            if ( j>cols ) then
               i = pointer_ind(cols)
               pointer_ind(cols) = pointer_ind(j)
               rho(cols) = rho(j)
               r0(cols) = r0(j)
               pointer_ind(j) = i
               call zcopy(nlen,vecs(1,cols,1),1,vecs(1,j,1),1)
            endif
            do i = 1,5
               call zcopy(nlen,work(1,i),1,vecs(1,cols,i),1)
            enddo
         else
            if ( REVCOM>1 ) then
               write(6,'(''REVCOM is greater than 1.'')')
               IQMR = 1
               return
            endif
!           ==========================================================
!           Check why the solver stopped (this could be more compact).
!           ==========================================================
            if ( IERR==0 ) then
!              write(6,'(A32)') 'The residual norm has converged.'
!              save the converged solution
               call zcopy (NLEN,work(1,1),1,vecs(1,j,1),1)
            else if ( IERR==1 ) then
               write(6,'(A35)') 'Invalid reverse communication call.'
               IQMR = IERR
               return
            else if ( IERR==2 ) then
               write(6,'(A27)') 'Invalid inputs encountered.'
               IQMR = IERR
               return
            else if ( IERR==4 ) then
               write(6,'(A31)') 'The algorithm did not converge.'
               IQMR = IERR
               return
            else if ( IERR==8 ) then
               write(6,'(A25)') 'The algorithm broke down.'
               IQMR = IERR
               return
            else 
               write(6,'(A19,I5)') 'Unknown INFO code: ', IERR
               IQMR = IERR
               return
            endif
         endif
!
!     Do the next column.
!
      enddo

      if ( cols/=0 ) then
!        ============================================================
!        Multiply work(1,CX) with the preconditioned matrix.
!
!        if ( cols>=1 ) then
!           call zgemm('N','n',nlen,cols,NLEN,CONE,AWAS,NDIM,         &
!                  vecs(1,1,cx),nlen,CZERO,vecs(1,1,CB),nlen)
!        else
!           call zgemv('N',nlen,NLEN,CONE,AWAS,NDIM,                  &
!                      vecs(1,1,cx),1,CZERO,vecs(1,1,CB),1)
!        endif
!        ============================================================
         call zcopy(nlen*cols,vecs(1,1,cx),1,vecs(1,1,6),1)
         call ztrtrs('u','n','u',nlen,cols,awas,ndim,vecs(1,1,6),     &
                     nlen,info)
         do i = 1,cols
            do j = 1,nlen
               vecs(j,i,cb)=vecs(j,i,cx)-vecs(j,i,6)
            enddo
         enddo
         call ztrtrs('l','n','u',nlen,cols,awas,ndim,vecs(1,1,cb),    &
                      nlen,info)
         do i = 1,cols
            do j = 1,nlen
               vecs(j,i,cb)=vecs(j,i,cb)+vecs(j,i,6)
            enddo
         enddo
         if ( n==nlim ) then
            iqmr = 1
            return
         endif
!
      else
         exit Loop_n
      endif
   enddo Loop_n
!
!  ================================================================
!  Compute the unpreconditioned solution.
!  ================================================================
   call ztrtrs('u','n','u',nlen,kkrsz,awas,ndim,vecs,nlen,info)
   do j = 1,kkrsz
      call zaxpy(nlen,CONE,vecs(1,j,1),1,rhs(1,pointer_ind(j)),1)
   enddo
!
   end subroutine wasinv
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine wasinv_p(NLEN,KKRSZ,AWAS,ndim,rhs,ldr,VECS,NLIM,         &
                       work,TOL,ALG,IQMR)
!  ===================================================================
!
!  Purpose:
!
!  Parameters:
!  NDIM  = the declared size of A (input).
!  NLEN  = the size of the matrix (input).
!  KKRSZ =
!  A     = the matrix, NLEN by NLEN (input).
!  rhs   = the right hand sides, nlen by kkrsz
!  VECS  =
!  NLIM  = the maximum number of iterations (input).
!  work  = work space, nlen by 10
!  TOL   = the relative convergence tolerance (input).
!  ALG   = the choice of algorithm (input).
!  IQMR  = indicates an error and will switch to an LU based algorithm
!
!     External routines used:
!
!     Noel M. Nachtigal
!     March 4, 1994
!
!  *******************************************************************
!
   implicit none
!
   integer (kind=IntKind), intent(in) ::   ALG
   integer (kind=IntKind), intent(in) ::   KKRSZ
   integer (kind=IntKind), intent(in) ::   NDIM,ldr
   integer (kind=IntKind), intent(in) ::   NLEN
   integer (kind=IntKind), intent(in) ::   NLIM
!
   complex (kind=CmplxKind), intent(inout), target :: AWAS(NDIM,nlen)
   complex (kind=CmplxKind), intent(inout), target :: rhs(ldr,kkrsz)
   complex (kind=CmplxKind), intent(inout), target :: VECS(nlen,kkrsz)
   complex (kind=CmplxKind), intent(inout), target :: work(nlen,9)
!
   real (kind=RealKind), intent(in) ::     TOL
!
!     Variables used by reverse communication.
!
   integer (kind=IntKind) ::  CB
   integer (kind=IntKind) ::  CX
   integer (kind=IntKind) ::  IERR
   integer (kind=IntKind) ::  REVCOM
   integer (kind=IntKind) ::  IQMR
!
!     Local variables.
!
   integer (kind=IntKind) ::  J
   integer (kind=IntKind) ::  INFO(4)
   integer (kind=IntKind) ::  INLIM
!
!  ==================================================================
!  Compute the modified right hand side.
!  ==================================================================
   do j = 1,kkrsz
      call zcopy(nlen,rhs(1,j),1,vecs(1,j),1)
   enddo
   call zgemm('n','N',nlen,kkrsz,NLEN,-CONE,AWAS,NDIM,rhs,ldr,       &
              CONE,VECS,nlen)
!
!  ==================================================================
!  Loop over all the columns.
!  ==================================================================
!
   do J = 1,KKRSZ
!     ===============================================================
!     Compute Jth column of the inverse.
!     ===============================================================
      call zcopy(nlen,vecs(1,j),1,work(1,2),1)
!     ===============================================================
!     Set up call to linear systems solver.
!     Compute true residual norms, generate second starting vector.
!     ===============================================================
      INFO(2) = 0
!        INFO(1) = 010006
!        INFO(1) = 000006
      INFO(1) = 000000
      INLIM   = NLIM
!     ===============================================================
!     call the solver.
!     ===============================================================
      REVCOM = 1
      do while ( REVCOM==1 .or. REVCOM==2 )
         if ( ALG==1 ) then
!           ---------------------------------------------------------
            call ZUCPX (nlen,NLEN,INLIM,work(1,1),TOL,INFO)
!           ---------------------------------------------------------
         else if ( ALG==2 ) then
!           ---------------------------------------------------------
            call ZUQMX (nlen,NLEN,INLIM,work(1,1),TOL,INFO)
!           ---------------------------------------------------------
         else if ( ALG==3 ) then
!           ---------------------------------------------------------
            call ZUTFX (nlen,NLEN,INLIM,work(1,1),TOL,INFO)
!           ---------------------------------------------------------
         endif
!        ============================================================
!        Perform matrix-vector multiplications when needed.
!        ============================================================
         IERR   = INFO(1)
         REVCOM = INFO(2)
         CX     = INFO(3)
         CB     = INFO(4)
         if (REVCOM==1) then
!           =========================================================
!           Multiply work(1,CX) with the preconditioned matrix.
!           =========================================================
!           ---------------------------------------------------------
            call zgemv('N',nlen,NLEN,CONE,AWAS,NDIM,work(1,CX),1,    &
                       CZERO,work(1,CB),1)
!           ---------------------------------------------------------
         else if (REVCOM==2) then
!           =========================================================
!           Multiply work(1,CX) with the preconditioned transpose.
!           =========================================================
!           ---------------------------------------------------------
            call zgemv('t',nlen,NLEN,CONE,AWAS,NDIM,work(1,CX),1,    &
                       CZERO,work(1,CB),1)
!           ---------------------------------------------------------
         endif
      enddo 
!     ===============================================================
!     Check why the solver stopped (this could be more compact).
!     ===============================================================
      if ( IERR==0 ) then
!        write(6,'(A32)') 'The residual norm has converged.'
!         call zcopy (NLEN,work(1,1),1,VECS(1,J),1)
      else if ( IERR==1 ) then
         write(6,'(A35)') 'Invalid reverse communication call.'
         IQMR = IERR
         return
      else if ( IERR==2 ) then
         write(6,'(A27)') 'Invalid inputs encountered.'
         IQMR = IERR
         return
      else if ( IERR==4 ) then
         write(6,'(A31)') 'The algorithm did not converge.'
         IQMR=IERR
         return
      else if ( IERR==8 ) then
         write(6,'(A25)') 'The algorithm broke down.'
         IQMR=IERR
         return
      else
         if ( (ALG==1 .OR. ALG==2) .and. IERR/=0 ) then
            if ( IERR==16 ) then
               write(6,'(A39)') 'An A-invariant subspace has been found.'
            else if ( IERR==32 ) then
               write(6,'(A41)') 'An A^T-invariant subspace has been found.'
            else if ( IERR==48 ) then
               write(6,'(A41)') 'Both invariant subspaces have been found.'
            endif
         else
            write(6,'(A19,I5)') 'Unknown INFO code: ', IERR
         endif
      endif
!     ===============================================================
!     Compute the unpreconditioned solution.
!     ===============================================================
!     call ZAXPY (NLEN,CONE,work(1,1),1,VECS(1,J),1)
      call zcopy (NLEN,work(1,1),1,VECS(1,J),1)
!     ===============================================================
!     Do the next column.
!     ===============================================================
   enddo
!
   do j = 1,kkrsz
        call zaxpy(nlen,CONE,vecs(1,j),1,rhs(1,j),1)
   enddo
!
   end subroutine wasinv_p
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine find_sym( id, delta, nlen )
!  ===================================================================
!  *******************************************************************
!     find the symmetry relations between the columns of delta
!     The symmetry matrices are stored in sym_ops
!  *******************************************************************
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, nlen
!
   complex (kind=CmplxKind), intent(in) :: delta(nlen,nlen)
!
   integer (kind=IntKind) :: i,j,k,ioff,k0,itmp
   integer (kind=IntKind) :: ks,ip(2)
!   integer (kind=IntKind), pointer :: p_ipvt(:), p_idcol(:)
!
   real (kind=RealKind) :: c, sum1, sum2, absq
!
   complex (kind=CmplxKind) :: a,b,v1,v2,s,x
   complex (kind=CmplxKind), pointer :: sym_ops(:,:,:)
!
!   interface
!      subroutine zeroout(xr,nx)
!      use KindParamModule, only : IntKind,RealKind
!      implicit none
!         integer (kind=IntKind), intent(in) :: nx
!         real (kind=RealKind), intent(out) :: x(nx)
!      end subroutine zeroout
!   end interface
!
   absq(x) = real(x*conjg(x),RealKind)
!  ===================================================================
!  For each column, check to see if the norm is the same as a previous
!  column
!  ===================================================================
   sym_ops => gettarget_3c( nlen, nlen, (nlen-1)/2, symop_wspace )
   ioff = 1
   idcol(1,id) = 1
   Loop_i: do i = 2,nlen
!     ================================================================
!     Symmetry is turned off by the following goto statement
!       goto 11
      if ( .not.isSymmetryOn ) then
         idcol(i,id)=i
         cycle Loop_i
      endif
!     if ioff exceeds the dimension of sym_op we forget about symmetry
!     ================================================================
      if ( ioff>(nlen-1)/2 ) then
         idcol(i,id)=i
         cycle Loop_i
      endif
      Loop_j: do j = 1,i-1
         sum1 = ZERO
         sum2 = ZERO
         do k = 1,nlen
            sum1 = sum1+absq(delta(k,j))
            sum2 = sum2+absq(delta(k,i))
         enddo
         if ( abs(sum1-sum2)<ten2m10*sum1 ) then
!           ==========================================================
!           If the norms are the same they are possibly equivelent
!           =======================================================
            do k = 1,nlen
               ipvt(k) = k
               sym_ops(k,1,ioff) = ZERO
            enddo
            ks = 0
            Loop_k: do k = 1,nlen
!              =======================================================
!              Do not do any zeros since they will mess up the matching
!              =======================================================
               if ( absq(delta(k,i))<ten2m20*sum1) then
                  cycle Loop_k
               endif
               do k0 = 1,nlen
                  if ( sym_ops(ipvt(k0),1,ioff)==ZERO ) then
                     if ( absq(delta(k,i)-delta(ipvt(k0),j)) &
                          <ten2m20*sum1 ) then
                        itmp     = ipvt(k)
                        ipvt(k)  = ipvt(k0)
                        ipvt(k0) =itmp
                        sym_ops(ipvt(k),1,ioff) = ONE
                        cycle Loop_k
                     endif
                  endif
               enddo  ! k0
!              =======================================================
!              We have not found any match
!              =======================================================
               ks = ks+1
!              =======================================================
!              If there are more than two elements not matched we
!              give up and quit
!              =======================================================
               if ( ks>2 ) then
                  idcol(i,id) = i
                  cycle Loop_i
               endif
               ip(ks) = k
            enddo  Loop_k
!           ==========================================================
!           Match the ks elements with the nonzero elements in jth col
!           ==========================================================
            do k = 1,ks
               if ( absq(delta(ipvt(ip(k)),j))<ten2m20*sum1 ) then
                  do k0 = 1,nlen
                     if ( sym_ops(ipvt(k0),1,ioff)==ZERO .and.     &
                          absq(delta(ipvt(k0),j))>ten2m20*sum1 ) then
                        itmp = ipvt(ip(k))
                        ipvt(ip(k)) = ipvt(k0)
                        ipvt(k0) = itmp
                     endif
                  enddo  ! k0
               endif
            enddo  ! ks
!           ==========================================================
!           Now match the zeros
!           ==========================================================
            do k = 1,nlen
               if ( absq(delta(k,i))<ten2m20*sum1 ) then
                  if ( absq(delta(ipvt(k),j))<ten2m20*sum1) then
!                    =================================================
!                    The zeros match
!                    =================================================
                     sym_ops(ipvt(k),1,ioff) = ONE
                  else
!                    =================================================
!                    The zero's don't match
!                    =================================================
                     ks = ks+1
                     if ( ks>2 ) then
                        idcol(i,id) = i
                        cycle Loop_i
                     endif
                     ip(ks) = k
                  endif
               endif
            enddo  ! k
!
!            call zeroout(sym_ops(1,1,ioff),2*nlen*nlen)
!
            sym_ops(1:nlen,1:nlen,ioff) = CZERO
            do k = 1,nlen
               sym_ops(ipvt(k),k,ioff) = ONE
            enddo
!
!           ==========================================================
!           now find the rotation matrices
!           ==========================================================
            if ( ks==2 ) then
               k  = ipvt(ip(1))
               k0 = ipvt(ip(2))
               v2 = delta(k,j)
               b  = delta(k0,j)
               call zrotg(v2,b,c,s)
               sym_ops(ip(1),k,ioff)  = c
               sym_ops(ip(2),k,ioff)  = -conjg(s)
               sym_ops(ip(1),k0,ioff) = s
               sym_ops(ip(2),k0,ioff) = c
!
               v1 = delta(ip(1),i)
               b  = delta(ip(2),i)
               call zrotg(v1,b,c,s)
               if ( absq(v1+v2)>ten2m20*sum1 ) then
                  c = -c
                  s = -s
               else if ( absq(v1-v2)>ten2m20*sum1 ) then
                  idcol(i,id) = i
                  cycle Loop_i
!                 ====================================================
!                 write(6,'(''find_sym:: A bug in find_sym!'')')
!                 call fstop('find_sym')
!                 ====================================================
               endif
               a = c*sym_ops(ip(1),k,ioff)-conjg(s)*sym_ops(ip(2),k,ioff)
               b = s*sym_ops(ip(1),k,ioff)+c*sym_ops(ip(2),k,ioff)
               sym_ops(ip(1),k,ioff) = a
               sym_ops(ip(2),k,ioff) = b
               a = c*sym_ops(ip(1),k0,ioff)-conjg(s)*sym_ops(ip(2),k0,ioff)
               b = s*sym_ops(ip(1),k0,ioff)+c*sym_ops(ip(2),k0,ioff)
               sym_ops(ip(1),k0,ioff) = a
               sym_ops(ip(2),k0,ioff) = b
            else if ( ks==1 ) then
                 idcol(i,id) = i
                 cycle Loop_i
!              =======================================================
!              ks cannot be 1
!              write(6,'(''FIND_SYM::Code doesnt work for the symmetry'')')
!              write(6,'(''FIND_SYM::ks is 1'')')
!              call fstop('find_sym')
!              =======================================================
            endif  ! ks.eq.2
!
            if ( print_level(id) > 0 ) then
               write(6,'(''col '',i3,'' is equiv to col '',i3)') i,j
!              =======================================================
!              do k = 1,nlen
!                 if ( absq(delta(k,j))>Ten2m20 .or.                     &
!                    absq(delta(k,i))>Ten2m20) then
!                    write(6,'(i3,4e15.6)') k,delta(k,j),delta(k,i)
!                 endif
!              enddo
!              call wrtmtx(sym_ops(1,1,ioff),nlen,' ')
!              =======================================================
            endif
            idcol(i,id) = j
            ioff = ioff+1
            cycle Loop_i
         endif
      enddo Loop_j
      idcol(i,id) = i
   enddo Loop_i
!
   end subroutine find_sym 
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function gettarget_2c(size1,size2,array)             result(parray)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: size1, size2
   complex(kind=CmplxKind), intent(in), target :: array(size1,size2)
!
   complex(kind=CmplxKind), pointer :: parray(:,:)
!
   parray=> array(1:size1,1:size2)
!
   end function gettarget_2c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function gettarget_3c(size1,size2,size3,array)       result(parray)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: size1, size2, size3
   complex(kind=CmplxKind), intent(in), target :: array(size1,size2,size3)
!
   complex(kind=CmplxKind), pointer :: parray(:,:,:)
!
   parray=> array(1:size1,1:size2,1:size3)
!
   end function gettarget_3c
!  ===================================================================
end module MatrixBlockInversionModule
