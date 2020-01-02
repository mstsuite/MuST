program testMPP
   use KindParamModule, only : IntKind, RealKind
   use TimerModule, only : initTimer, getTime
   use MPPModule, only : initMPP, endMPP
   use MPPModule, only : NumPEs, MyPE
   use MPPModule, only : syncAllPEs
   use MPPModule, only : openLocalMemory, closeLocalMemory, readRemoteMemory
   implicit none
!
   integer (kind=IntKind) :: Win
   integer (kind=IntKind) :: i, j, rank1, rank2, rank3
   integer (kind=IntKind), parameter :: nsize = 81*81*2
!
   real (kind=RealKind), allocatable :: base(:), rslt1(:), rslt2(:), rslt3(:)
   real (kind=RealKind) :: t0
!
   call initMPP()
!
   call initTimer()
!
   allocate( base(nsize), rslt1(nsize), rslt2(nsize), rslt3(nsize) )
!
   do i=1,nsize
      base(i) = MyPE*nsize+i
   enddo
!
!  ===================================================================
!  open local data array, base, for remote processors to read.
!  Note: openLocalMemory is a global call, and all processors are synchronized.
!        Each call to openLocalMemoery returns a unique integer key, Win,
!        for later use to identify the openned data space.
!  -------------------------------------------------------------------
   call openLocalMemory(base, nsize, Win)
!  -------------------------------------------------------------------
!
   if (MyPE == 0) then
      print *,'One-sided communication test ========================='
   endif
!
   rank1 = mod(NumPEs+MyPE-1,NumPEs)
   rank2 = mod(MyPE+1,NumPEs)
   rank3 = mod(MyPE+2,NumPEs)
!
   call syncAllPEs()
   t0 = getTime()
!
   do j=1,100
!     ================================================================
!     call readRemoteMemory to read the data on read remote processors,
!     rank1, rank2, and rank3. The data is identified by the integer
!     key, Win.
!     ================================================================
      call readRemoteMemory(rslt1, nsize, rank1, Win, 1)
      call readRemoteMemory(rslt2, nsize, rank2, Win, 1)
      if (MyPE < NumPEs/2) then
         call readRemoteMemory(rslt3, nsize, rank3, Win, 1)
      endif
!
      do i=1,nsize
         if (abs(rslt1(i)-rank1*nsize-i) > 0.0000001d0) then
            write(6,'(a,3i5,5x,f5.0)')'MyPE, i, j, rslt1 = ',MyPE,i,j,rslt1(i)
         endif
         if (abs(rslt2(i)-rank2*nsize-i) > 0.0000001d0) then
            write(6,'(a,3i5,5x,f5.0)')'MyPE, i, j, rslt2 = ',MyPE,i,j,rslt2(i)
         endif
         if (MyPE < NumPEs/2) then
            if(abs(rslt3(i)-rank3*nsize-i) > 0.0000001d0) then
               write(6,'(a,3i5,5x,f5.0)')'MyPE, i, j, rslt3 = ',MyPE,i,j,rslt3(i)
            endif
         endif
      enddo
      rslt1(1:nsize) =0.0d0
      rslt2(1:nsize) =0.0d0
      rslt3(1:nsize) =0.0d0
   enddo
!
   call syncAllPEs()
!
   if (MyPE == 0) then
      print *,'CPU Time = ',getTime()-t0,' (second)'
   endif
!
!  ===================================================================
!  close the previously openned local data space from remote access.
!  The local data space is identified by integer key, Win.
!  ===================================================================
   call closeLocalMemory(Win)
!
!  ===================================================================
!
   deallocate( base, rslt1, rslt2, rslt3 )
!
   call endMPP()
!
end program testMPP
