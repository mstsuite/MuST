#ifndef MaxOutProcs
#define MaxOutProcs MaxIOPEs
#endif
#ifndef MaxInProcs
#define MaxInProcs MaxIOPEs
#endif
!
module ParallelIOModule
   use KindParamModule, only : IntKind, RealKind, LongIntKind, CmplxKind
   use MathParamModule, only : ZERO
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initParallelIO,        &
          endParallelIO,         &
          getMyInputProc,        &  ! returns the proc that handles my input
          getNumInputClients,    &  ! returns the number of procs being my clients for input
          getInputClient,        &  ! returns a client proc for input
          getMyOutputProc,       &  ! returns the proc that handles my output
          getNumOutputClients,   &  ! returns the number of procs being my clients for output
          getOutputClient,       &  ! returns a client proc for output
          isInputProc,           &
          isInputHeadProc,       &
          isOutputProc,          &
          isOutputHeadProc,      &
          getIOCommunicator,     &
          getMyPEinIOGroup,      &
          getNumPEsinIOGroup,    &
!         openParallel,          &
!         closeParallel,         &
          readParallel,          &
          writeParallel,         &
          getNumOpenFiles,       &
          getFilePointer,        &
          setFilePointer,        &
          rewindFilePointer,     &
          forwardFilePointer
!
   interface initParallelIO
      module procedure initPIO1, initPIO2
   end interface
!
   interface readParallel
      module procedure readPll_char
      module procedure readPll_int
      module procedure readPll_real
      module procedure readPll_cmplx
   end interface
!
   interface writeParallel
      module procedure writePll_char
      module procedure writePll_int
      module procedure writePll_real
      module procedure writePll_cmplx
   end interface
!
private
   integer (kind=IntKind), parameter :: MaxIOPEs=4096
!
   integer (kind=IntKind) :: MaxInputPEs = min(MaxInProcs,MaxIOPEs)
   integer (kind=IntKind) :: MaxOutputPEs = min(MaxOutProcs,MaxIOPEs)
!
   integer (kind=IntKind) :: GroupID = -1
   integer (kind=IntKind) :: MyInputProc = -1
   integer (kind=IntKind) :: MyOutputProc = -1
   integer (kind=IntKind) :: NumInputClients = 0
   integer (kind=IntKind) :: NumOutputClients = 0
   integer (kind=IntKind), allocatable :: InputClient(:)
   integer (kind=IntKind), allocatable :: OutputClient(:)
!
   integer (kind=IntKind) :: MyProcInGroup
   integer (kind=IntKind) :: NumProcsInGroup
   integer (kind=IntKind) :: HeadProc_read, Headproc_write
!
   integer (kind=IntKind) :: IOCommunicator
!
   logical :: ImInputProc = .false.
   logical :: ImOutputProc = .false.
!
!  ===================================================================
!  the following variables are from the original code.................
!  ===================================================================
   integer (kind=IntKind) :: ParallelIOFlag = 0
   integer (kind=IntKind), parameter :: MaxOpenUnits = 10
   integer (kind=IntKind), parameter :: MessageType = 130005
   integer (kind=IntKind) :: UnitNumber(MaxOpenUnits)
   integer (kind=LongIntKind) :: CurrentRec(MaxOpenUnits)
   integer (kind=IntKind) :: NumOpenUnits
   integer (kind=IntKind) :: RECORD_LENGTH
   integer (kind=IntKind) :: NUMC_PER_REC
   integer (kind=IntKind) :: IOHostPE
   integer (kind=IntKind) :: NumClientPEs = 0 ! Eventually, it needs to
                                              ! be changed to NumInputClient
                                              ! or NumOutputClients
   integer (kind=IntKind), allocatable :: msgid(:)
!  ===================================================================
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initPIO1(ig,id,nin,nout)
!  ===================================================================
   use GroupCommModule, only : getMyPEinGroup, getNumPEsInGroup,      &
                               getMyClusterIndex, getGroupCommunicator
   implicit none
!
   integer (kind=IntKind), intent(in) :: ig
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), optional, intent(in) :: nin, nout
!
   if (present(nin)) then
      if (nin < 1) then
         call ErrorHandler('initParallelIO','Number of input processes < 1',nin)
      endif
      MaxInputPEs = nin
   else
      MaxInputPEs = min(MaxInProcs,MaxIOPEs)
   endif
   if (present(nout)) then
      if (nout < 1) then
         call ErrorHandler('initParallelIO','Number of output processes < 1',nout)
      endif
      MaxOutputPEs = nout
   else
      MaxOutputPEs = min(MaxOutProcs,MaxIOPEs)
   endif
!
   GroupID = ig
   MyProcInGroup = getMyPEinGroup(ig)
   NumProcsInGroup = getNumPEsInGroup(ig)
   IOCommunicator = getGroupCommunicator(ig)
!
   if (getMyClusterIndex(ig) == id) then
!     ----------------------------------------------------------------
      call setupIOClient(.true.)
!     ----------------------------------------------------------------
   else ! this group will not perform output
!     ----------------------------------------------------------------
      call setupIOClient(.false.)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   call determineHeadProc()
!  -------------------------------------------------------------------
!
   end subroutine initPIO1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initPIO2()
!  ===================================================================
   use MPPModule, only : MyPE, NumPEs
   implicit none
!
   MaxInputPEs = min(MaxInProcs,MaxIOPEs)
   MaxOutputPEs = min(MaxOutProcs,MaxIOPEs)
!
   MyProcInGroup = MyPE
   NumProcsInGroup = NumPEs
!
!  -------------------------------------------------------------------
   call setupIOClient(.true.)
!  -------------------------------------------------------------------
   call determineHeadProc()
!  -------------------------------------------------------------------
!
   end subroutine initPIO2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupIOClient(w)
!  ===================================================================
   implicit none
!
   logical, intent(in) :: w
!
   integer (kind=IntKind) :: n, i, j
!
!  ===================================================================
!  determine I/O processor and its clients.
!  Note: if the number of atoms on each processor is different, the IO
!        will hang.
!  ===================================================================
   MyInputProc = mod(MyProcInGroup,MaxInputPEs)
!
   if (MyInputProc == MyProcInGroup) then
      ImInputProc = .true.
      n = 0
      do i=1,NumProcsInGroup
         if (i-1 /= MyProcInGroup) then
            j=mod(i-1,MaxInputPEs)
            if (j == MyProcInGroup) then
               n = n + 1
            endif 
         endif
      enddo
      NumInputClients = n
      if (NumInputClients > 0) then
         allocate( InputClient(NumInputClients) )
         n = 0
         do i=1,NumProcsInGroup
            if (i-1 /= MyProcInGroup) then
               j=mod(i-1,MaxInputPEs)
               if (j == MyProcInGroup) then
                  n = n + 1
                  InputClient(n) = i-1
               endif
            endif
         enddo
      endif
   else
      ImInputProc = .false.
      NumInputClients = 0
   endif
!
   if ( w ) then
      MyOutputProc = mod(MyProcInGroup,MaxOutputPEs)
      if (MyProcInGroup == 0) then
         write(6,'(a,i6)')'The number of processes that write output data: ', &
                          min(NumProcsInGroup,MaxOutputPEs)
      endif
   else
      MyOutputProc = -1
   endif
!
   if (MyOutputProc == MyProcInGroup) then
      ImOutputProc = .true.
      n = 0
      do i=1,NumProcsInGroup
         if (i-1 /= MyProcInGroup) then
            j=mod(i-1,MaxOutputPEs)
            if (j == MyProcInGroup) then
               n = n + 1
            endif 
         endif
      enddo
      NumOutputClients = n
      if (NumOutputClients > 0) then
         allocate( OutputClient(NumOutputClients) )
         n = 0
         do i=1,NumProcsInGroup
            if (i-1 /= MyProcInGroup) then
               j=mod(i-1,MaxOutputPEs)
               if (j == MyProcInGroup) then
                  n = n + 1
                  OutputClient(n) = i-1
               endif
            endif
         enddo
      endif
   else
      ImOutputProc = .false.
      NumOutputClients = 0
   endif
!
   end subroutine setupIOClient
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine determineHeadProc()
!  ===================================================================
   use MPPModule, only : MyPE, GlobalMax
!
   implicit none
!
   integer (kind=IntKind) :: head_proc(2)
!
   if ( ImInputProc ) then
      head_proc(1) = MyPE
   else
      head_proc(1) = 0
   endif
!
   if ( ImOutputProc ) then
      head_proc(2) = MyPE
   else
      head_proc(2) = 0
   endif
!  -------------------------------------------------------------------
   call GlobalMax(head_proc,2)
!  -------------------------------------------------------------------
   HeadProc_read = head_proc(1)
   HeadProc_write = head_proc(2)
!
   end subroutine determineHeadProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupFileRecord()
!  ===================================================================
   implicit none
!
   character (len=128) :: c = ' '
!
   integer (kind=IntKind) :: length
!
   real (kind=RealKind) :: a = ZERO
   real (kind=RealKind) :: mp
!
!  ===================================================================
!  The following lines are from the original code. We do not use them
!  for now.......
!  ===================================================================
   if (ParallelIOFlag == 0) then
      ParallelIOFlag = 1
   else
      return
   endif
   inquire(iolength=RECORD_LENGTH) a
!
   NUMC_PER_REC=0
   length=1
   do while (length==1 .and. NUMC_PER_REC<=128)
      NUMC_PER_REC=NUMC_PER_REC+1
      inquire(iolength=length) c(1:NUMC_PER_REC)
   enddo
   NUMC_PER_REC=NUMC_PER_REC-1
!
   NumOpenUnits=0
   UnitNumber(1:MaxOpenUnits)=0
   CurrentRec(1:MaxOpenUnits)=0
!
   IOHostPE = mod(MyProcInGroup,MaxIOPEs)
   if (MyProcInGroup == IOHostPE) then
      mp=MaxIOPEs
      NumClientPEs=floor((NumProcsInGroup-IOHostPE-1.0)/mp)
   else
      NumClientPEs=0
   endif
   if (NumClientPEs > 0) then
      allocate(msgid(NumClientPEs))
   endif
!
   end subroutine setupFileRecord
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endParallelIO()
!  ===================================================================
   implicit none
!
   if (NumInputClients > 0) then
      deallocate(InputClient)
   endif
!
   if (NumOutputClients > 0) then
      deallocate(OutputClient)
   endif
!
   GroupID = -1
   MyInputProc = -1
   MyOutputProc = -1
   NumInputClients = 0
   NumOutputClients = 0
   ImInputProc = .false.
   ImOutputProc = .false.
!
   if (allocated(msgid)) then
      deallocate(msgid)
   endif
!
   end subroutine endParallelIO
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine openParallelIO(funit,fname,fstatus)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
   use GroupCommModule, only : SyncAllPEsInGroup
!
   implicit none
   character (len=*), intent(in) :: fname
   character (len=*), intent(in) :: fstatus
!
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: ios
!
   if (ParallelIOFlag == 0) then
      call ErrorHandler('openParallelIO','need to call initParallelIO first')
   else if (funit < 1) then
      call ErrorHandler('openParallelIO','funit < 1',funit)
   endif
!
   call SyncAllPEsInGroup(GroupID)
   if (MyProcInGroup < MaxIOPEs) then
      open(unit=funit,file=fname,status=fstatus,form='unformatted',      &
           access='direct',action='read',recl=RECORD_LENGTH,iostat=ios)
   else
      ios = 0
   endif
!
   if (ios > 0) then
      call ErrorHandler('openParallelIO','error opening file, iostat > 0',ios)
   endif
!
   NumOpenUnits=NumOpenUnits+1
   UnitNumber(NumOpenUnits)=funit
   CurrentRec(NumOpenUnits)=1
!
   end subroutine openParallelIO
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine closeParallelIO(funit)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
   use GroupCommModule, only : SyncAllPEsInGroup
!
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: ios, i
!
   if (funit < 1) then
      call ErrorHandler('closeParallelIO','funit < 1',funit)
   endif
!
   call SyncAllPEsInGroup(GroupID)
   if (MyProcInGroup < MaxIOPEs) then
      close(unit=funit,iostat=ios)
   else
      ios = 0
   endif
!
   if (ios > 0) then
      call ErrorHandler('closeParallelIO','error closing file, iostat > 0',ios)
   endif
!
   do i=1,NumOpenUnits
      if (UnitNumber(i) == funit) then
         UnitNumber(i)=UnitNumber(NumOpenUnits)
         CurrentRec(i)=CurrentRec(NumOpenUnits)
         UnitNumber(NumOpenUnits)=0
         CurrentRec(NumOpenUnits)=0
         NumOpenUnits=NumOpenUnits-1
         exit
      endif
   enddo
!
   end subroutine closeParallelIO
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readPll_char(funit,s,l)
!  ===================================================================
   use MPPModule, only : AnyPE, SyncAllPEs, waitMessage
   use MPPModule, only : packMessage, unpackMessage
   use MPPModule, only : sendPackage, recvPackage
   use MPPModule, only : nbsendMessage, nbrecvMessage, recvMessage
!
   implicit none
   integer (kind=IntKind), intent(in) :: l
   character (len=1), intent(out) :: s(l)
   character (len=NUMC_PER_REC) :: c
   character (len=1), allocatable :: st(:)
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: i,j,k,n,m,ncr,ncrt,lt,TargetPE
   integer (kind=LongIntKind) :: fp, fpt
   real (kind=RealKind) :: rtmp
!
   c(1:NUMC_PER_REC) = ' '
   rtmp=NUMC_PER_REC
   ncr = ceiling(l/rtmp)
   fp=CurrentRec(funit)-1
   fpt = 0
   if (MyProcInGroup == IOHostPE) then
      n=0
      do j=1,ncr
         read(funit,rec=fp+j)c
         m=min(l-n,NUMC_PER_REC)
         do i=1,m
            s(n+i)=c(i:i)
         enddo
         n=n+NUMC_PER_REC
      enddo
      do k=1,NumClientPEs
!        -------------------------------------------------------------
         call recvPackage(MessageType,AnyPE)
         call unpackMessage(TargetPE)
         call unpackMessage(fpt)
         call unpackMessage(lt)
!        -------------------------------------------------------------
         allocate(st(lt))
         ncrt = ceiling(lt/rtmp)
         n=0
         do j=1,ncrt
            read(funit,rec=fpt+j)c
            m=min(lt-n,NUMC_PER_REC)
            do i=1,m
               st(n+i)=c(i:i)
            enddo
            n=n+NUMC_PER_REC
         enddo
         if (k > 1) then
!           ----------------------------------------------------------
            call waitMessage(msgid(k-1))
!           ----------------------------------------------------------
         endif
!        -------------------------------------------------------------
         msgid(k)=nbsendMessage(st,lt,MessageType+1,TargetPE)
!        -------------------------------------------------------------
         deallocate(st)
      enddo
!     ----------------------------------------------------------------
      call waitMessage(msgid(NumClientPEs))
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call packMessage(MyProcInGroup)
      call packMessage(fp)
      call packMessage(l)
      call sendPackage(MessageType,IOHostPE)
      call recvMessage(s,l,MessageType+1,IOHostPE)
!     ----------------------------------------------------------------
   endif
   CurrentRec(funit)=CurrentRec(funit)+ncr
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
   end subroutine readPll_char
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readPll_int(funit,s,l)
!  ===================================================================
   use MPPModule, only : AnyPE, SyncAllPEs, waitMessage
   use MPPModule, only : packMessage, unpackMessage
   use MPPModule, only : sendPackage, recvPackage
   use MPPModule, only : nbsendMessage, nbrecvMessage, recvMessage
!
   implicit none
   integer (kind=IntKind), intent(in) :: l
   integer (kind=IntKind), intent(out) :: s(l)
   integer (kind=IntKind), allocatable :: st(:)
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: j,k,lt,TargetPE
   integer (kind=LongIntKind) :: fp,fpt
!
   fp=CurrentRec(funit)-1
   fpt = 0
   if (MyProcInGroup == IOHostPE) then
      do j=1,l
         read(funit,rec=fp+j)s(j)
      enddo
      do k=1,NumClientPEs
!        -------------------------------------------------------------
         call recvPackage(MessageType,AnyPE)
         call unpackMessage(TargetPE)
         call unpackMessage(fpt)
         call unpackMessage(lt)
!        -------------------------------------------------------------
         allocate(st(lt))
         do j=1,lt
            read(funit,rec=fpt+j)st(j)
         enddo
         if (k > 1) then
!           ----------------------------------------------------------
            call waitMessage(msgid(k-1))
!           ----------------------------------------------------------
         endif
!        -------------------------------------------------------------
         msgid(k)=nbsendMessage(st,lt,MessageType+1,TargetPE)
!        -------------------------------------------------------------
         deallocate(st)
      enddo
!     ----------------------------------------------------------------
      call waitMessage(msgid(NumClientPEs))
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call packMessage(MyProcInGroup)
      call packMessage(fp)
      call packMessage(l)
      call sendPackage(MessageType,IOHostPE)
      call recvMessage(s,l,MessageType+1,IOHostPE)
!     ----------------------------------------------------------------
   endif
   CurrentRec(funit)=CurrentRec(funit)+l
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
   end subroutine readPll_int
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readPll_real(funit,s,l)
!  ===================================================================
   use MPPModule, only : AnyPE, SyncAllPEs, waitMessage
   use MPPModule, only : packMessage, unpackMessage
   use MPPModule, only : sendPackage, recvPackage
   use MPPModule, only : nbsendMessage, nbrecvMessage, recvMessage
!
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: l
   real (kind=RealKind), intent(out) :: s(l)
   real (kind=RealKind), allocatable :: st(:)
   integer (kind=IntKind) :: j,k,lt,TargetPE
   integer (kind=LongIntKind) :: fp,fpt
!
   fp=CurrentRec(funit)-1
   fpt = 0
   if (MyProcInGroup == IOHostPE) then
      do j=1,l
         read(funit,rec=fp+j)s(j)
      enddo
      do k=1,NumClientPEs
!        -------------------------------------------------------------
         call recvPackage(MessageType,AnyPE)
         call unpackMessage(TargetPE)
         call unpackMessage(fpt)
         call unpackMessage(lt)
!        -------------------------------------------------------------
         allocate(st(lt))
         do j=1,lt
            read(funit,rec=fpt+j)st(j)
         enddo
         if (k > 1) then
!           ----------------------------------------------------------
            call waitMessage(msgid(k-1))
!           ----------------------------------------------------------
         endif
!        -------------------------------------------------------------
         msgid(k)=nbsendMessage(st,lt,MessageType+1,TargetPE)
!        -------------------------------------------------------------
         deallocate(st)
      enddo
!     ----------------------------------------------------------------
      call waitMessage(msgid(NumClientPEs))
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call packMessage(MyProcInGroup)
      call packMessage(fp)
      call packMessage(l)
      call sendPackage(MessageType,IOHostPE)
      call recvMessage(s,l,MessageType+1,IOHostPE)
!     ----------------------------------------------------------------
   endif
   CurrentRec(funit)=CurrentRec(funit)+l
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
   end subroutine readPll_real
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readPll_cmplx(funit,s,l)
!  ===================================================================
   use MPPModule, only : AnyPE, SyncAllPEs, waitMessage
   use MPPModule, only : packMessage, unpackMessage
   use MPPModule, only : sendPackage, recvPackage
   use MPPModule, only : nbsendMessage, nbrecvMessage, recvMessage
!
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: l
   complex (kind=CmplxKind), intent(out) :: s(l)
   complex (kind=CmplxKind), allocatable :: st(:)
   real (kind=RealKind) :: cr, ci
   integer (kind=IntKind) :: j,k,lt,n,TargetPE
   integer (kind=LongIntKind) :: fp,fpt
!
   fp=CurrentRec(funit)-1
   fpt = 0
   if (MyProcInGroup == IOHostPE) then
      n=2*l
      do j=1,n,2
         read(funit,rec=fp+j)cr
         read(funit,rec=fp+j+1)ci
         s(j)=cmplx(cr,ci,CmplxKind)
      enddo
      do k=1,NumClientPEs
!        -------------------------------------------------------------
         call recvPackage(MessageType,AnyPE)
         call unpackMessage(TargetPE)
         call unpackMessage(fpt)
         call unpackMessage(lt)
!        -------------------------------------------------------------
         allocate(st(lt))
         n=2*lt
         do j=1,n,2
            read(funit,rec=fpt+j)cr
            read(funit,rec=fpt+j+1)ci
            st(j)=cmplx(cr,ci,CmplxKind)
         enddo
         if (k > 1) then
!           ----------------------------------------------------------
            call waitMessage(msgid(k-1))
!           ----------------------------------------------------------
         endif
!        -------------------------------------------------------------
         msgid(k)=nbsendMessage(st,lt,MessageType+1,TargetPE)
!        -------------------------------------------------------------
         deallocate(st)
      enddo
!     ----------------------------------------------------------------
      call waitMessage(msgid(NumClientPEs))
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call packMessage(MyProcInGroup)
      call packMessage(fp)
      call packMessage(l)
      call sendPackage(MessageType,IOHostPE)
      call recvMessage(s,l,MessageType+1,IOHostPE)
!     ----------------------------------------------------------------
   endif
   CurrentRec(funit)=CurrentRec(funit)+2*l
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
   end subroutine readPll_cmplx
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writePll_char(funit,s,l)
!  ===================================================================
   use MPPModule, only : AnyPE, SyncAllPEs
   use MPPModule, only : packMessage, unpackMessage
   use MPPModule, only : sendPackage, recvPackage
!
   implicit none
   integer (kind=IntKind), intent(in) :: l
   character (len=1), intent(in) :: s(l)
   character (len=NUMC_PER_REC) :: c
   character (len=1), allocatable :: st(:)
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: i,j,k,n,m,ncr,ncrt,lt
   integer (kind=LongIntKind) :: fp,fpt
   real (kind=RealKind) :: rtmp
!
   rtmp=NUMC_PER_REC
   ncr = ceiling(l/rtmp)
   fp=CurrentRec(funit)-1
   fpt = 0
   if (MyProcInGroup == IOHostPE) then
      n=0
      do j=1,ncr
         m=min(l-n,NUMC_PER_REC)
         do i=1,m
            c(i:i)=s(n+i)
         enddo
         do i=m+1,NUMC_PER_REC
            c(i:i)=' '
         enddo
         n=n+NUMC_PER_REC
         write(funit,rec=fpt+j)c
      enddo
      do k=1,NumClientPEs
!        -------------------------------------------------------------
         call recvPackage(MessageType,AnyPE)
         call unpackMessage(fpt)
         call unpackMessage(lt)
         allocate(st(lt))
         call unpackMessage(st,lt)
!        -------------------------------------------------------------
         ncrt = ceiling(lt/rtmp)
         n=0
         do j=1,ncrt
            m=min(lt-n,NUMC_PER_REC)
            do i=1,m
               c(i:i)=st(n+i)
            enddo
            do i=m+1,NUMC_PER_REC
               c(i:i)=' '
            enddo
            n=n+NUMC_PER_REC
            write(funit,rec=fpt+j)c
         enddo
         deallocate(st)
      enddo
   else
!     ----------------------------------------------------------------
      call packMessage(fp)
      call packMessage(l)
      call packMessage(s,l)
      call sendPackage(MessageType,IOHostPE)
!     ----------------------------------------------------------------
   endif
   CurrentRec(funit)=CurrentRec(funit)+ncr
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
   end subroutine writePll_char
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writePll_int(funit,s,l)
!  ===================================================================
   use MPPModule, only : AnyPE, SyncAllPEs
   use MPPModule, only : packMessage, unpackMessage
   use MPPModule, only : sendPackage, recvPackage
!
   implicit none
   integer (kind=IntKind), intent(in) :: l
   integer (kind=IntKind), intent(in) :: s(l)
   integer (kind=IntKind), allocatable :: st(:)
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: j,k,lt
   integer (kind=LongIntKind) :: fp,fpt
!
   fp=CurrentRec(funit)-1
   fpt = 0
   if (MyProcInGroup == IOHostPE) then
      do j=1,l
         write(funit,rec=fp+j)s(j)
      enddo
      do k=1,NumClientPEs
!        -------------------------------------------------------------
         call recvPackage(MessageType,AnyPE)
         call unpackMessage(fpt)
         call unpackMessage(lt)
         allocate(st(lt))
         call unpackMessage(st,lt)
!        -------------------------------------------------------------
         do j=1,lt
            write(funit,rec=fpt+j)st(j)
         enddo
         deallocate(st)
      enddo
   else
!     ----------------------------------------------------------------
      call packMessage(fp)
      call packMessage(l)
      call packMessage(s,l)
      call sendPackage(MessageType,IOHostPE)
!     ----------------------------------------------------------------
   endif
   CurrentRec(funit)=CurrentRec(funit)+l
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
   end subroutine writePll_int
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writePll_real(funit,s,l)
!  ===================================================================
   use MPPModule, only : AnyPE, SyncAllPEs
   use MPPModule, only : packMessage, unpackMessage
   use MPPModule, only : sendPackage, recvPackage
!
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: l
   real (kind=RealKind), intent(in) :: s(l)
   real (kind=RealKind), allocatable :: st(:)
   integer (kind=IntKind) :: j,k,lt
   integer (kind=LongIntKind) :: fp,fpt
!
   fp=CurrentRec(funit)-1
   fpt = 0
   if (MyProcInGroup == IOHostPE) then
      do j=1,l
         write(funit,rec=fp+j)s(j)
      enddo
      do k=1,NumClientPEs
!        -------------------------------------------------------------
         call recvPackage(MessageType,AnyPE)
         call unpackMessage(fpt)
         call unpackMessage(lt)
         allocate(st(lt))
         call unpackMessage(st,lt)
!        -------------------------------------------------------------
         do j=1,lt
            write(funit,rec=fpt+j)st(j)
         enddo
         deallocate(st)
      enddo
   else
!     ----------------------------------------------------------------
      call packMessage(fp)
      call packMessage(l)
      call packMessage(s,l)
      call sendPackage(MessageType,IOHostPE)
!     ----------------------------------------------------------------
   endif
   CurrentRec(funit)=CurrentRec(funit)+l
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
   end subroutine writePll_real
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writePll_cmplx(funit,s,l)
!  ===================================================================
   use MPPModule, only : AnyPE, SyncAllPEs
   use MPPModule, only : packMessage, unpackMessage
   use MPPModule, only : sendPackage, recvPackage
!
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: l
   complex (kind=CmplxKind), intent(in) :: s(l)
   complex (kind=CmplxKind), allocatable :: st(:)
   integer (kind=IntKind) :: j,k,lt,n
   integer (kind=LongIntKind) :: fp,fpt
!
   fp=CurrentRec(funit)-1
   fpt = 0
   if (MyProcInGroup == IOHostPE) then
      n=2*l
      do j=1,n,2
         write(funit,rec=fp+j)real(s(j),RealKind)
         write(funit,rec=fp+j+1)aimag(s(j))
      enddo
      do k=1,NumClientPEs
!        -------------------------------------------------------------
         call recvPackage(MessageType,AnyPE)
         call unpackMessage(fpt)
         call unpackMessage(lt)
         allocate(st(lt))
         call unpackMessage(st,lt)
!        -------------------------------------------------------------
         n=2*lt
         do j=1,n,2
            write(funit,rec=fpt+j)real(st(j),RealKind)
            write(funit,rec=fpt+j+1)aimag(st(j))
         enddo
         deallocate(st)
      enddo
   else
!     ----------------------------------------------------------------
      call packMessage(fp)
      call packMessage(l)
      call packMessage(s,l)
      call sendPackage(MessageType,IOHostPE)
!     ----------------------------------------------------------------
   endif
   CurrentRec(funit)=CurrentRec(funit)+2*l
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
   end subroutine writePll_cmplx
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumOpenFiles() result(num)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: num
!
   num = NumOpenUnits
!
   end function getNumOpenFiles
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getFilePointer(funit) result(fp)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: fp
!
   fp = CurrentRec(funit)
!
   end function getFilePointer
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setFilePointer(funit,rec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: rec
!
   CurrentRec(funit) = rec
!
   end subroutine setFilePointer
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine forwardFilePointer(funit,rec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: rec
!
   CurrentRec(funit) = CurrentRec(funit) + rec
!
   end subroutine forwardFilePointer
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine rewindFilePointer(funit,rec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(in) :: rec
!
   CurrentRec(funit) = CurrentRec(funit) - rec
   if (CurrentRec(funit) < 1) then
      CurrentRec(funit) = 1
   endif
!
   end subroutine rewindFilePointer
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMyInputProc() result (pe)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: pe
!
   pe = MyInputProc
!
   end function getMyInputProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumInputClients() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumInputClients
!
   end function getNumInputClients
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInputClient(i) result(pe)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: pe
!
   if (i < 1 .or. i > NumInputClients) then
      call ErrorHandler('getInputClient','Invalid i',i)
   endif
!
   pe = InputClient(i)
!
   end function getInputClient
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMyOutputProc() result (pe)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: pe
!
   pe = MyOutputProc
!
   end function getMyOutputProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumOutputClients() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumOutputClients
!
   end function getNumOutputClients
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOutputClient(i) result(pe)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: pe
!
   if (i < 1 .or. i > NumOutputClients) then
      call ErrorHandler('getOutputClient','Invalid i',i)
   endif
!
   pe = OutputClient(i)
!
   end function getOutputClient
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isInputProc() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   y = ImInputProc
!
   end function isInputProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isInputHeadProc() result(y)
!  ===================================================================
   use MPPMOdule, only : MyPE
!
   implicit none
!
   logical :: y
!
   if ( ImInputProc .and. MyPE == HeadProc_read ) then
      y = .true.
   else
      y = .false.
   endif
!
   end function isInputHeadProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isOutputProc() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   y = ImOutputProc
!
   end function isOutputProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isOutputHeadProc() result(y)
!  ===================================================================
   use MPPMOdule, only : MyPE
!
   implicit none
!
   logical :: y
!
   if ( ImOutputProc .and. MyPE == HeadProc_write ) then
      y = .true.
   else
      y = .false.
   endif
!
   end function isOutputHeadProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIOCommunicator() result(comm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: comm
!
   comm = IOCommunicator
!
   end function getIOCommunicator
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMyPEinIOGroup() result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: p
!
   p = MyProcInGroup
!
   end function getMyPEinIOGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumPEsInIOGroup() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumProcsInGroup
!
   end function getNumPEsInIOGroup
!  ===================================================================
end module ParallelIOModule
