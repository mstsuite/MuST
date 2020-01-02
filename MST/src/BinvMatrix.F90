!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine BinvMatrix( nblk,blk_sz,b,ldb,nb,         &
                          a,lda,na,alg )
!  ================================================================
!
!  ****************************************************************
!  PURPOSE:   inverse the first block of a complex matrix: a
!             (a^{-1})_00=(a_00-b)^{-1}
!
!  INPUT:     a,      the complex matrix to be inverted
!             alg,        the algorithm of the invertion
!                         = 1, UL method
!                         = 2, LU method
!                         = 3,4,5, W.A.S. method
!
!  OUTPUT:    a,   contains the first block of the inverted
!                         matrix
!             b,
!  ****************************************************************
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : cone, czero, ten2m8
   use ErrorHandlerModule, only : ErrorHandler
!
   use TimerModule, only : getTime
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nblk
   integer (kind=IntKind), intent(in) :: lda, na, ldb, nb, alg
   integer (kind=IntKind), intent(in) :: blk_sz(nblk)
!
!  real (kind=RealKind), intent(out) :: time_direct
!
   complex (kind=CmplxKind), intent(out) :: a(lda,na)
   complex (kind=CmplxKind), intent(out) :: b(ldb,nb)
!
   integer (kind=IntKind) :: i, j, k, ioff, iqmr, blk1
   integer (kind=IntKind), allocatable :: idcol(:)
!
   complex (kind=CmplxKind), allocatable :: sym_ops(:,:,:)
!
!  write(6,*)"start BinvMatrix"
!  ===================================================================
!  Invert the KKR-Matrix using the method defined by 'alg'
!  If using QMR and either the method fails or it takes longer than
!  the direct method (LU or UL) then the routine uses LU thereafter
!  ===================================================================
!
!  t0 = getTime()
!
   blk1 = blk_sz(1)
   allocate( idcol(blk1))
   idcol(:) = 0
!
   do i = 1,ldb
      do j = 1,nb
         a(j,i)=czero
      enddo
   enddo
   if ( alg.le.2.or.alg.ge.10 ) then
!     ================================================================
!     Use the  LU algorithm
!     ================================================================
!     ----------------------------------------------------------------
      call zblock_lu(a,lda,na,blk_sz,nblk,idcol,k)
!     ----------------------------------------------------------------
   else if ( alg.ge.3 .and. alg.le.6 ) then
      call ErrorHandler('block_inv','Error: not defined yet.',alg)
      return
!     ================================================================
!     Use the QMR algorithm
!     ================================================================
      call zblock_lu(a,lda,na,blk_sz,nblk,idcol,k)
!     ----------------------------------------------------------------
!      call wasinv( na-blk_sz(1),k,a(blk_sz(1)+1,blk_sz(1)+1),lda,      &
!                   a(blk_sz(1)+1,blk_sz(1)-k+1),lda,vecs,              &
!                   nlim,vecs(1,blk_sz(1)*6+1),                         &
!                   rwork,work1,iwork,ten2m8,alg-2,iqmr)
!     ================================================================
!     If qmr algorithm fails then redo the current energy: Use LU..
!     ================================================================
      if ( iqmr .gt. 0 ) then
         call zblock_lu(a,lda,na,blk_sz,nblk,idcol,k)
      else
         call zgemm( 'n','n',blk1,k,na-blk1,-cone,a(1,blk1+1),lda,     &
                      a(blk1+1,blk1-k+1),lda,czero,a,lda)
      endif
   else
      call ErrorHandler('block_inv','incorect alg value',alg)
   endif
!
   if ( idcol(1) == 0 ) then
      do i = 1,blk1
         do j = 1,blk1
            b(j,i) = -a(j,i)
         enddo
      enddo
!     Here look for symetry
   else
      ioff = 0
      do i = 1,blk1
         if ( idcol(i) /= i ) then
            ioff = ioff+1
         endif
      enddo
      if (ioff > 0) then
         allocate( sym_ops(blk1,blk1,ioff) )
      endif
      k = 0
      ioff = 0
      do i = 1,blk1
         if ( idcol(i).eq.i ) then
            k = k+1
            do j = 1,blk1
               b(j,i) = -a(j,k)
            enddo
         else
            j = idcol(i)
            ioff = ioff+1
            call zgemv('n',blk1,blk1,cone,sym_ops(1,1,ioff),           &
                        blk1,b(1,j),1,czero,b(1,i),1)
         endif
      enddo
      if (ioff > 0) then
         deallocate( sym_ops )
      endif
   endif
!  print *, "end BinvMtx"
!  print *,'Inverse time = ',getTime()-t0,' sec'
!
   end subroutine BinvMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine zblock_lu( a, lda, na, blk_sz, nblk, idcol, k )
!  ===================================================================
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
   use KindParamModule, only : IntKind, CmplxKind
   use MathParamModule, only : cone
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   integer (kind=IntKind), intent(in) ::  lda, na, nblk
   integer (kind=IntKind), intent(in) :: blk_sz(nblk)
   integer (kind=IntKind), intent(out) :: k
   integer (kind=IntKind), intent(in) ::  idcol(blk_sz(1))
!
   integer (kind=IntKind) :: i, m, n, info
   integer (kind=IntKind) :: ioff, joff, iblk
   integer (kind=IntKind), allocatable :: ipvt(:)
   complex (kind=CmplxKind), intent(inout) :: a(lda,na)
!
!   print *, "start zblock_lu"
!  ===================================================================
!  eliminate columns that are equiv due to symmetry
!  ===================================================================
!
   k = blk_sz(1)+1
   do i = blk_sz(1),1,-1
      if ( idcol(1).eq.0 .or. idcol(i).eq.i ) then
         k = k-1
         if ( k.ne.i ) then
            call zcopy( na-blk_sz(1),a(blk_sz(1)+1,i),1,               &
                       a(blk_sz(1)+1,k),1 )
         endif
      endif
   enddo
!
!   print *, "idclo(1)",idcol(1)
   if ( nblk.gt.0 ) then
!     ================================================================
!     Do block LU
!     ================================================================
      n = blk_sz(nblk)
      joff = na-n
      do iblk = nblk,2,-1
         m = n
         ioff = joff
         n = blk_sz(iblk-1)
         joff = joff-n
!        =============================================================
!        invert the diagonal blk_sz(iblk) x blk_sz(iblk) block
!        =============================================================
         allocate( ipvt(m) )
!        -------------------------------------------------------------
         call zgetrf( m,m,a(ioff+1,ioff+1),lda,ipvt,info )
         if (info /= 0) then
            call ErrorHandler('BinvMatrix','zgetrf failed', info)
         endif
!        -------------------------------------------------------------
!        =============================================================
!        calculate the inverse of above multiplying the row block
!        blk_sz(iblk) x ioff
!        =============================================================
!        -------------------------------------------------------------
         call zgetrs( 'n',m,ioff,a(ioff+1,ioff+1),lda,ipvt,            &
                       a(ioff+1,1),lda,info )
         if (info /= 0) then
            call ErrorHandler('BinvMatrix','zgetrs failed', info)
         endif
         deallocate( ipvt )
!        -------------------------------------------------------------
         if ( iblk.gt.2 ) then
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
      call zgemm( 'n','n',blk_sz(1),blk_sz(1)-k+1,na-blk_sz(1),-cone,  &
                a(1,blk_sz(1)+1),lda,a(blk_sz(1)+1,k),lda,cone,a,lda )
!     ----------------------------------------------------------------
   endif
!
   k = blk_sz(1)-k+1
!
!  print *,"end zblock_lu"
!
   end subroutine zblock_lu
!  ===================================================================
