program testPZGETRI
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
   use MathParamModule, only : CZERO, CONE, TEN2m8, SQRTm1
   use MPPModule, only : initMPP, syncAllPEs, endMPP, NumPEs, MyPE, GlobalSum, &
                         getCommunicator
   use GroupCommModule, only : initGroupComm, endGroupComm, &
                               createProcGrid, createBoxGroupFromGrid, &
                               getProcGridID, getGroupID, getGroupCommunicator, &
                               getNumPEsInGroup, getMyPEinGroup, GlobalSumInGroup
   implicit none
!! include 'mpif.h'
!
   integer (kind=IntKind), parameter :: ndim = 100
   integer (kind=IntKind), parameter :: DLEN_ = 9
   integer (kind=IntKind) :: ICTXT
   integer (kind=IntKind) :: NPROW
   integer (kind=IntKind) :: NPCOL
   integer (kind=IntKind) :: MYROW
   integer (kind=IntKind) :: MYCOL
   integer (kind=IntKind) :: DESC_A( DLEN_ )
   integer (kind=IntKind) :: INFO
   integer (kind=IntKind) :: LWORK
   integer (kind=IntKind) :: LIWORK
   integer (kind=IntKind) :: IPVT(ndim)
!
   integer (kind=IntKind), allocatable :: IWORK(:)
   complex (kind=CmplxKind), allocatable :: WORK(:)
   complex (kind=CmplxKind) :: tmp
   complex (kind=CmplxKind), parameter :: alpha = 0.01d0+SQRTm1
!
!  integer (kind=IntKind) :: MyPE, NumPEs
   integer (kind=IntKind) :: i, j, k, n, lk, ln, m
   integer (kind=IntKind) :: BlockSize
   integer (kind=IntKind) :: itmp, ictxt1, ictxt2, grid, group
   integer (kind=IntKind) :: dim(1), box(1), color, key, my_rank, group_size
   integer (kind=Intkind) :: new_comm, old_comm, old_group, new_group
!
   complex (kind=CmplxKind) :: am(ndim,ndim)
   complex (kind=CmplxKind) :: bm(ndim,ndim)
   complex (kind=CmplxKind) :: uij(ndim)
   complex (kind=CmplxKind) :: bam(ndim,ndim)
   complex (kind=CmplxKind), allocatable :: pam(:,:)
!
#ifdef USE_SCALAPACK
   interface 
      function BLACS_PNUM(ic,j,k) result(i)
         integer :: ic, j, k
         integer :: i
      end function BLACS_PNUM
   end interface
!
   interface 
      function NUMROC(ic,j,k,l,m) result(i)
         integer :: ic, j, k, l, m
         integer :: i
      end function NUMROC
   end interface
!
   interface 
      function indxg2p(ic,j,k,l,m) result(i)
         integer, intent(in) :: ic, j, k, l, m
         integer :: i
      end function indxg2p
   end interface
!
!  -------------------------------------------------------------------
   call initMPP()
!  call BLACS_PINFO(i, n)
!  -------------------------------------------------------------------
! MyPE = i
! NumPEs = n
!  call SL_INIT(ICTXT,1,NumPEs)
!  if (NumPEs /= n .or. MyPE /= i) then
!     print *,'Error in calling BLACS_PINFO:'
!     print *,'MyPE, i = ',MyPE, i
!     print *,'NumPEs, n = ',NumPEs, n
!  endif
!
!  -------------------------------------------------------------------
   call initMatrix(ndim, ndim, am)
   call zcopy(ndim*ndim,am,1,bm,1)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call ZGETRF(ndim,ndim,am,ndim,IPVT,INFO)
!  -------------------------------------------------------------------
   if (INFO /= 0) then
      call ErrorHandler('testPZGETRI','Return from ZGETRF: info <> 0',info)
   endif
!
!  -------------------------------------------------------------------
   call ZGETRI(ndim,am,ndim,IPVT,tmp,-1,INFO)
!  -------------------------------------------------------------------
   LWORK=int(real(tmp,kind=RealKind))
!  LWORK = ndim*ndim
!
   allocate( WORK(1:LWORK) )
!  -------------------------------------------------------------------
   call ZGETRI(ndim,am,ndim,IPVT,WORK,LWORK,INFO)
!  -------------------------------------------------------------------
   if (INFO /= 0) then
      call ErrorHandler('testPZGETRI','Return from ZGETRI: info <> 0',info)
   endif
!
   call zgemm('n','n',ndim,ndim,ndim,CONE,am,ndim,bm,ndim,CZERO,bam,ndim)
   do j=1,ndim
      do i=1,ndim
         if (i /= j .and. abs(bam(i,j)) > TEN2m8) then
            print *,'For i = ',i,', j = ',j,', bam = ',bam(i,j)
         else if (i == j .and. abs(bam(i,j)-CONE) > TEN2m8) then
            print *,'For i = ',i,', j = ',j,', bam = ',bam(i,j)
         endif
      enddo
   enddo
   deallocate( WORK )
!
   write(6,'(/,a)')'Start testing on ScaLapack ...'
!  -------------------------------------------------------------------
   call BLACS_GET(-1, 0, ICTXT) 
   call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEs )
   call BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!  -------------------------------------------------------------------
   write(6,'(a,2i5)')'MyPE, ICTXT = ',MyPE, ICTXT
   write(6,'(a,5i5)')'MyPE, NPROW, NPCOL, MYROW, MYCOL = ',MyPE, NPROW, NPCOL, MYROW, MYCOL
!  -------------------------------------------------------------------
   if (NPROW /= 1 .or. NPCOL /= NumPEs .or. MYROW /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('testPZGETRI',                                 &
              'Failed: NPROW /= 1 || NPCOL /= NumPEs || MYROW /= 0',  &
              NPROW, NPCOL, MYROW)
!     ----------------------------------------------------------------
   else if (MYCOL /= MyPE) then
!     ----------------------------------------------------------------
      call ErrorHandler('testPZGETRI','MYCOL /= MyPE',MYCOL,MyPE)
!     ----------------------------------------------------------------
   endif
!
   BlockSize = ndim/NumPEs
!
   if (MyPE==0) then
      do i=0,(NPROW*NPCOL)-1
         call blacs_pcoord(ictxt,i,j,k)
         print *,'Processor ',i,' is at grid position (p,q)=',j,k
      end do
      do j=0,NPROW-1
         do k=0,NPCOL-1
            print *,'At grid position (',j,',',k,') is processor', blacs_pnum(ictxt,j,k)
         end do
      end do
   end if
!
   call syncAllPEs()
!
   print * ,'MyPE: the number of rows in my local array is ',MyPE, NUMROC(ndim,BlockSize,MYROW,0,NPROW)
   print * ,'MyPE: the number of cols in my local array is ',MyPE, NUMROC(ndim,BlockSize,MYCOL,0,NPCOL)
!
   call syncAllPEs()
!  -------------------------------------------------------------------
   call DESCINIT( DESC_A, ndim, ndim, BlockSize, BlockSize, 0, 0, ICTXT, ndim, INFO )
!  -------------------------------------------------------------------
   if (INFO /= 0) then
      call ErrorHandler('testPZGETRI','Return from DESCINIT: info <> 0',info)
   endif
   print *,'DESC_A(6) = ',DESC_A(6)
!
!  -------------------------------------------------------------------
   call initMatrix(ndim, ndim, am)
!  -------------------------------------------------------------------
   allocate( pam(ndim,BlockSize) )
!
!  lk = 0
!  do k=1,ndim
!     j=INDXG2P(k,BlockSize,0,0,NPCOL)
!     if (j == MYCOL) then
!        lk = lk + 1
!        ln = 0
!        do n=1,ndim
!           i=INDXG2P(n,BlockSize,0,0,NPROW)
!           if (i == MYROW) then
!              ln = ln + 1
!              pam(ln,lk) = am(n,k)
!           endif
!           print *,'Element (',n,',',k,') is on proc (',i,',',j,')'
!        enddo
!     endif
!  end do
   pam(1:ndim,1:BlockSize) = am(1:ndim,MYCOL*BlockSize+1:(MYCOL+1)*BlockSize)
!
!  -------------------------------------------------------------------
   call PZGETRF(ndim, ndim, pam, 1, 1, DESC_A, IPVT, INFO)
!  -------------------------------------------------------------------
   if (INFO /= 0) then
      call ErrorHandler('testPZGETRI','Return from PZGETRF: info <> 0',info)
   endif
!
!  -------------------------------------------------------------------
   call PZGETRI(ndim, pam, 1, 1, DESC_A, IPVT, tmp, -1, itmp, -1, INFO )
!  -------------------------------------------------------------------
   LWORK = int(real(tmp, kind=RealKind))
   LIWORK = itmp
   print *,'MyPE, LWORK, LIWORK = ',MyPE, LWORK, LIWORK
!  LWORK = ndim*BlockSize
!  LIWORK = ndim*BlockSize
   allocate( WORK(1:LWORK), IWORK(1:LIWORK) )
!
!  -------------------------------------------------------------------
   call PZGETRI(ndim, pam, 1, 1, DESC_A, IPVT, WORK, LWORK, IWORK, LIWORK, INFO )
!  -------------------------------------------------------------------
   if (INFO /= 0) then
      call ErrorHandler('testPZGETRI','Return from PZGETRI: info <> 0',info)
   endif
   deallocate( WORK, IWORK )
!
!  do k = 1, NumPEs-1
!     n = mod(MyPE+k,NumPEs)
!     call sendMessage(BlockSize,BlockSize,pam,n)
!     m = mod(MyPE-k+NumPEs,NumPEs)
!     call recvMessage(BlockSize,BlockSize,qam,m)
!     call blacs_pcoord(ictxt,m,i,j)
!
   bm(1:ndim,1:ndim) = CZERO
!  lk = 0
!  do k=1,ndim
!     j=INDXG2P(k,BlockSize,0,0,NPCOL)
!     if (j == MYCOL) then
!        lk = lk + 1
!        ln = 0
!        do n=1,ndim
!           i=INDXG2P(n,BlockSize,0,0,NPROW)
!           if (i == MYROW) then
!              ln = ln + 1
!              bm(n,k) = pam(ln,lk)
!           endif
!        enddo
!     endif
!  end do
   bm(1:ndim,MYCOL*BlockSize+1:(MYCOL+1)*BlockSize) = pam(1:ndim,1:BlockSize)
   call GlobalSum(bm,ndim,ndim)
!
!  -------------------------------------------------------------------
   call zgemm('n','n',ndim,ndim,ndim,CONE,am,ndim,bm,ndim,CZERO,bam,ndim)
!  -------------------------------------------------------------------
   LOOP_j: do j=1,ndim
      do i=1,ndim
         if (i /= j .and. abs(bam(i,j)) > TEN2m8) then
            print *,'For i = ',i,', j = ',j,', bam = ',bam(i,j)
            exit LOOP_j
         else if (i == j .and. abs(bam(i,j)-CONE) > TEN2m8) then
            print *,'For i = ',i,', j = ',j,', bam = ',bam(i,j)
            exit LOOP_j
         endif
      enddo
   enddo LOOP_j
!
   deallocate( pam )
!
   call syncAllPEs()
   CALL BLACS_GRIDEXIT( ICTXT )
!
!  ===================================================================
!  the following code tests the the situation that two groups of the 
!  processors using different handlers
!  ===================================================================
   write(6,'(/,a)')'Start a new test ...'
   call initGroupComm()
   dim(1) = NumPEs
   call createProcGrid(1,'1-d Grid',dim)
   grid = getProcGridID('1-d Grid')
   box(1)=NumPEs/2
   call createBoxGroupFromGrid(grid,box,'Half')
   group = getGroupID('Half')
   new_comm = getGroupCommunicator(group)
   old_comm = getCommunicator()
   write(6,'(a,i5,2x,2i15)')'MyPE, new_comm, MPI_comm_world = ',MyPE, new_comm, old_comm
   write(6,'(a,3i5)')'MyPE, rank, size = ',MyPE, getMyPEinGroup(group), getNumPEsInGroup(group)
   call MPI_COMM_GROUP(old_comm,old_group,info)
   call MPI_COMM_GROUP(new_comm,new_group,info)
   write(6,'(a,i5,2x,2i15)')'MyPE, old_group, new_group = ',MyPE, old_group, new_group
   call syncAllPEs()
!
   ICTXT = new_comm
!
   call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEs/2 )
   call BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!  -------------------------------------------------------------------
   write(6,'(a,2i10)')'MyPE, ICTXT = ',MyPE, ICTXT
   write(6,'(a,5i5)')'MyPE, NPROW, NPCOL, MYROW, MYCOL = ',MyPE, NPROW, NPCOL, MYROW, MYCOL
!  -------------------------------------------------------------------
   if (NPROW /= 1 .or. NPCOL /= NumPEs/2 .or. MYROW /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('testPZGETRI',                                 &
              'Failed: NPROW /= 1 || NPCOL /= NumPEs/2 || MYROW /= 0', &
              NPROW, NPCOL, MYROW)
!     ----------------------------------------------------------------
   else if (MYCOL /= mod(MyPE,NumPEs/2)) then
!     ----------------------------------------------------------------
      call ErrorHandler('testPZGETRI','MYCOL /= mod(MyPE,NumPEs/2)',MYCOL,MyPE)
!     ----------------------------------------------------------------
   endif
!
   BlockSize = 2*ndim/NumPEs
!
   if (MyPE==0 .or. MyPE==NumPEs/2) then
      do i=0,(NPROW*NPCOL)-1
         call blacs_pcoord(ictxt,i,j,k)
         print *,'Processor ',i+MyPE,' is at grid position (p,q)=',j,k
      end do
      do j=0,NPROW-1
         do k=0,NPCOL-1
            print *,'At grid position (',j,',',k,') is processor', blacs_pnum(ictxt,j,k)+MyPE
         end do
      end do
   end if
!
   call syncAllPEs()
!
   print * ,'MyPE: the number of rows in my local array is ',MyPE, NUMROC(ndim,BlockSize,MYROW,0,NPROW)
   print * ,'MyPE: the number of cols in my local array is ',MyPE, NUMROC(ndim,BlockSize,MYCOL,0,NPCOL)
!
   call syncAllPEs()
!  -------------------------------------------------------------------
   call DESCINIT( DESC_A, ndim, ndim, BlockSize, BlockSize, 0, 0, ICTXT, ndim, INFO )
!  -------------------------------------------------------------------
   if (INFO /= 0) then
      call ErrorHandler('testPZGETRI','Return from DESCINIT: info <> 0',info)
   endif
   print *,'DESC_A(6) = ',DESC_A(6)
!
!  -------------------------------------------------------------------
   call initMatrix(ndim, ndim, am)
!  -------------------------------------------------------------------
   if (MyPE >= NumPEs/2) then
      am(1:ndim,1:ndim) = alpha*am(1:ndim,1:ndim)
   endif
!
   allocate( pam(ndim,BlockSize) )
   pam(1:ndim,1:BlockSize) = am(1:ndim,MYCOL*BlockSize+1:(MYCOL+1)*BlockSize)
!
!  -------------------------------------------------------------------
   call PZGETRF(ndim, ndim, pam, 1, 1, DESC_A, IPVT, INFO)
!  -------------------------------------------------------------------
   if (INFO /= 0) then
      call ErrorHandler('testPZGETRI','Return from PZGETRF: info <> 0',info)
   endif
!
!  -------------------------------------------------------------------
   call PZGETRI(ndim, pam, 1, 1, DESC_A, IPVT, tmp, -1, itmp, -1, INFO )
!  -------------------------------------------------------------------
   LWORK = int(real(tmp, kind=RealKind))
   LIWORK = itmp
   print *,'MyPE, LWORK, LIWORK = ',MyPE, LWORK, LIWORK
!  LWORK = ndim*BlockSize
!  LIWORK = ndim*BlockSize
   allocate( WORK(1:LWORK), IWORK(1:LIWORK) )
!
!  -------------------------------------------------------------------
   call PZGETRI(ndim, pam, 1, 1, DESC_A, IPVT, WORK, LWORK, IWORK, LIWORK, INFO )
!  -------------------------------------------------------------------
   if (INFO /= 0) then
      call ErrorHandler('testPZGETRI','Return from PZGETRI: info <> 0',info)
   endif
   deallocate( WORK, IWORK )
!
   bm(1:ndim,1:ndim) = CZERO
   bm(1:ndim,MYCOL*BlockSize+1:(MYCOL+1)*BlockSize) = pam(1:ndim,1:BlockSize)
   call GlobalSumInGroup(group,bm,ndim,ndim)
!  -------------------------------------------------------------------
   call zgemm('n','n',ndim,ndim,ndim,CONE,am,ndim,bm,ndim,CZERO,bam,ndim)
!  -------------------------------------------------------------------
!! bm(1:ndim,1:ndim) = CZERO
!! if (MyPE >= NumPEs/2) then
!!    bm(1:ndim,MYCOL*BlockSize+1:(MYCOL+1)*BlockSize) = pam(1:ndim,1:BlockSize)
!! endif
!! call GlobalSum(bm,ndim,ndim)
!! if (MyPE >= NumPEs/2) then
!     ----------------------------------------------------------------
!!    call zgemm('n','n',ndim,ndim,ndim,CONE,am,ndim,bm,ndim,CZERO,bam,ndim)
!     ----------------------------------------------------------------
!! endif
!
   LOOP_j1: do j=1,ndim
      do i=1,ndim
         if (i /= j .and. abs(bam(i,j)) > TEN2m8) then
            print *,'On MyPE = ',MyPE,', for i = ',i,', j = ',j,', bam = ',bam(i,j)
            exit LOOP_j1
         else if (i == j .and. abs(bam(i,j)-CONE) > TEN2m8) then
            print *,'On MyPE = ',MyPE,', for i = ',i,', j = ',j,', bam = ',bam(i,j)
            exit LOOP_j1
         endif
      enddo
   enddo LOOP_j1
!
   deallocate( pam )
!
   CALL BLACS_GRIDEXIT( ICTXT )
!
   call syncAllPEs()
   call endGroupComm()
   call endMPP()
#endif
   stop 'Ok'
!
end program testPZGETRI
!  ===================================================================
subroutine initMatrix(n,m,a)
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, m
   integer (kind=IntKind) :: i, j
!
   real (kind=RealKind) :: t, rr, ri
!
   complex (kind=CmplxKind), intent(out) :: a(n,m)
!
   do j = 1, m
      do i = 1, n
!        rr = (11.0d0*i-3.0d0*j)/real(i+j,kind=RealKind)
!        ri = (5.0d0*i-2.0d0*j)/real(i+j,kind=RealKind)
!        rr = i + 2.0d0*j
!        ri = 2.0d0*i-j
         t= abs(i-j)/real(n, kind=RealKind)
         ri = ONE-sqrt(t)
         rr = ONE+sqrt(t)
         a(i,j) = cmplx(rr,ri,kind=CmplxKind)
      enddo
   enddo
end subroutine initMatrix
