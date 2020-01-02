program testBlockInverse
   use KindParamModule, only : IntKind, realKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
   use TimerModule, only : initTimer
   use MatrixBlockInversionModule, only : initMatrixBlockInv, endMatrixBlockInv, invertMatrixBlock
!
   implicit none
!
   logical :: FileExist
!
   integer (kind=IntKind) :: isize, jsize, dsize, i, j, my_atom, max_nblck, n
   integer (kind=IntKind) :: Minv_alg, NumEs, NumAtoms, LocalNumAtoms
   integer (kind=IntKind), allocatable :: nblck(:),kblocks(:,:), iprint(:)
!
   complex (kind=CmplxKind), allocatable :: BigMatrix(:), BigMatrix_in(:)
   complex (kind=CmplxKind), allocatable :: BlockMatrix(:,:), BlockMatrix_in(:,:)
!
!  -------------------------------------------------------------------
   call initTimer()
!  -------------------------------------------------------------------
!
!  ===================================================================
   inquire(file='BigMatrix.dat',exist=FileExist)
   if (.not.FileExist) then
      call ErrorHandler('testBlockInverse','BigMatrix.dat does not exist')
   endif
   open(unit=11,file='BigMatrix.dat',form='formatted',status='old')
   read(11,'(i5)')dsize
   if (dsize < 1) then
      call ErrorHandler('testBlockInverse','Invalid dsize value',dsize)
   endif
!  ===================================================================
!
!  ===================================================================
   inquire(file='BlockMatrix.dat',exist=FileExist)
   if (.not.FileExist) then
      call ErrorHandler('testBlockInverse','BlockMatrix.dat does not exist')
   endif
   open(unit=12,file='BlockMatrix.dat',form='formatted',status='old')
   read(12,'(3i5)')my_atom,isize,jsize
   if (isize < 1 .or. jsize < 1) then
      call ErrorHandler('testBlockInverse','Invalid isize or jsize value',isize,jsize)
   endif
!  ===================================================================
!
!  ===================================================================
   allocate(BigMatrix(dsize*dsize), BigMatrix_in(dsize*dsize))
   read(11,'(4d20.12)')(BigMatrix_in(i),i=1,dsize*dsize)
   close(11)
!
   allocate(BlockMatrix(isize,jsize), BlockMatrix_in(isize,jsize))
   read(12,'(4d20.12)')((BlockMatrix_in(i,j),i=1,isize),j=1,jsize)
   close(12)
!  ===================================================================
!
!  ===================================================================
   inquire(file='InitMatrix.dat',exist=FileExist)
   if (.not.FileExist) then
      call ErrorHandler('testBlockInverse','InitMatrix.dat does not exist')
   endif
   open(unit=51,file='InitMatrix.dat',form='formatted',status='old')
   read(51,'(4i5)')Minv_alg,LocalNumAtoms,NumEs,max_nblck
   allocate( nblck(LocalNumAtoms), kblocks(max_nblck,LocalNumAtoms), iprint(LocalNumAtoms) )
   do i = 1,LocalNumAtoms
      read(51,'(2i5)')nblck(i),NumAtoms
      read(51,'(20i5)')(kblocks(j,i),j=1,NumAtoms+1)
   enddo
   close(51)
   iprint(:) = 0
!  -------------------------------------------------------------------
   call initMatrixBlockInv(Minv_alg, LocalNumAtoms, NumEs, nblck, max_nblck, kblocks, iprint )
!  -------------------------------------------------------------------
   deallocate(nblck, kblocks)
!  ===================================================================
!
   do n = 1, 30 ! Perform the matrix inverse 30 times
      BigMatrix = BigMatrix_in
      print *,' '
      print *,'Calling invertMatrixBlock #',n
!     ----------------------------------------------------------------
      call invertMatrixBlock( my_atom, BlockMatrix, isize, jsize,        &
                              BigMatrix, dsize, dsize )
!     ----------------------------------------------------------------
      do j = 1, jsize
         do i = 1, isize
            if (abs(BlockMatrix(i,j)-BlockMatrix_in(i,j)) > 1.0d-8) then
!              ----------------------------------------------------
               call ErrorHandler('testBlockInverse','Matrix difference > 1.0d-8', &
                                 BlockMatrix(i,j),BlockMatrix_in(i,j))
!              ----------------------------------------------------
            endif
         enddo
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call endMatrixBlockInv()
!  -------------------------------------------------------------------
!
   deallocate(BigMatrix, BigMatrix_in, BlockMatrix, BlockMatrix_in)
!
   stop 'Ok'
end program testBlockInverse
