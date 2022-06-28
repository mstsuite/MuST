module CheckPointModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initCheckPoint,    &
          insertStopPoint,   &
          insertPrintPoint,  &
          isPrintPoint,      &
          isStopPoint,       &
          takePrintAction,   &
          takeStopAction,    &
          endCheckPoint
!
private
   integer (kind=IntKind) :: NumStopPoints = 0
   integer (kind=IntKind) :: NumPrintPoints = 0
!
   integer (kind=IntKind), parameter :: MaxNameSize = 50
!
   type CheckPointStruct
      integer (kind=IntKind) :: index
      character (len=MaxNameSize) :: rname
      type (CheckPointStruct), pointer :: next
   end type CheckPointStruct
!
   type(CheckPointStruct), pointer :: PrintPoint
   type(CheckPointStruct), pointer :: StopPoint
!
   interface
      function nocaseCompare(s1,s2) result(t)
         implicit none
         logical :: t
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
      end function nocaseCompare
   end interface
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initCheckPoint()
!  ===================================================================
   implicit none
!   
   NumStopPoints = 0
   NumPrintPoints = 0
!
   end subroutine initCheckPoint
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCheckPoint()
!  ===================================================================
   implicit none
!
   type(CheckPointStruct), pointer :: ptr, ptr_next
!
   integer (kind=IntKind) :: i
!
   ptr => StopPoint
   do i = 1, NumStopPoints
      ptr_next => ptr%next
      deallocate(ptr)
      ptr => ptr_next
   enddo
!
   ptr => PrintPoint
   do i = 1, NumPrintPoints
      ptr_next => ptr%next
      deallocate(ptr)
      ptr => ptr_next
   enddo
!   
   NumStopPoints = 0
   NumPrintPoints = 0
!
   end subroutine endCheckPoint
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insertStopPoint(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: p
!
   type(CheckPointStruct), pointer :: ptr
!
   integer (kind=IntKind) :: i
!
   if (NumStopPoints == 0) then
      allocate(StopPoint)
      ptr => StopPoint
   else
      do i = 1, NumStopPoints
         if (i == 1) then
            ptr => StopPoint
         else
            ptr => ptr%next
         endif
         if (nocaseCompare(p,ptr%rname)) then
            call WarningHandler('insertStopPoint','The print point has already inserted',trim(p))
            return
         endif
      enddo
      allocate(ptr%next)
      ptr => ptr%next
   endif
!
   NumStopPoints = NumStopPoints + 1
   ptr%index = NumStopPoints
   ptr%rname = trim(adjustl(p))
   nullify(ptr%next)
!
   end subroutine insertStopPoint
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insertPrintPoint(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: p
!
   type(CheckPointStruct), pointer :: ptr
!
   integer (kind=IntKind) :: i
!
   if (NumPrintPoints == 0) then
      allocate(PrintPoint)
      ptr => PrintPoint
   else
      do i = 1, NumPrintPoints
         if (i == 1) then
            ptr => PrintPoint
         else
            ptr => ptr%next
         endif
         if (nocaseCompare(p,ptr%rname)) then
            call WarningHandler('insertPrintPoint','The print point has already inserted',trim(p))
            return
         endif
      enddo
      allocate(ptr%next)
      ptr => ptr%next
   endif
!
   NumPrintPoints = NumPrintPoints + 1
   ptr%index = NumPrintPoints
   ptr%rname = trim(adjustl(p))
   nullify(ptr%next)
!
   end subroutine insertPrintPoint
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isPrintPoint(p,id) result(y)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: p
!
   logical :: y
!
   integer (kind=IntKind), intent(out) :: id
   integer (kind=IntKind) :: i
!
   type(CheckPointStruct), pointer :: ptr
!
   y = .false.; id = 0
   ptr => PrintPoint
   Loop_i: do i = 1, NumPrintPoints
      if (.not.associated(ptr)) then
!        -------------------------------------------------------------
         call ErrorHandler('isPrintPoint','End of data has reached at step',i)
!        -------------------------------------------------------------
      else if (i /= ptr%index) then
!        -------------------------------------------------------------
         call ErrorHandler('isPrintPoint','Invalid index value',i,ptr%index)
!        -------------------------------------------------------------
      endif
      if (nocaseCompare(p,ptr%rname)) then
         y = .true.
         id = ptr%index
         exit LOOP_i
      endif
      ptr => ptr%next
   enddo LOOP_i
!
   end function isPrintPoint
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isStopPoint(p,id) result(y)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: p
!
   logical :: y
!
   integer (kind=IntKind), intent(out) :: id
   integer (kind=IntKind) :: i
!
   type(CheckPointStruct), pointer :: ptr
!
   y = .false.; id = 0
   ptr => StopPoint
   Loop_i: do i = 1, NumStopPoints
      if (.not.associated(ptr)) then
!        -------------------------------------------------------------
         call ErrorHandler('isStopPoint','End of data has reached at step',i)
!        -------------------------------------------------------------
      else if (i /= ptr%index) then
!        -------------------------------------------------------------
         call ErrorHandler('isStopPoint','Invalid index value',i,ptr%index)
!        -------------------------------------------------------------
      endif
      if (nocaseCompare(p,ptr%rname)) then
         y = .true.
         id = ptr%index
         exit LOOP_i
      endif
      ptr => ptr%next
   enddo LOOP_i
!
   end function isStopPoint
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine takeStopAction(id,print_pe)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in), optional :: print_pe
   integer (kind=IntKind) :: i
!
   type(CheckPointStruct), pointer :: ptr
!
   if (id < 1 .or. id > NumStopPoints) then
      return
   endif
!
   ptr => StopPoint
   LOOP_i: do i = 1, NumStopPoints
      if (id == ptr%index) then
         exit LOOP_i
      endif
      ptr => ptr%next
   enddo LOOP_i
!
   if (present(print_pe)) then
      if (print_pe == MyPE .or. print_pe == -1) then
         write(6,'(2a)')'The process stops at check point: ',trim(ptr%rname)
      endif
   endif
!  -------------------------------------------------------------------
   call syncAllPEs()
!  -------------------------------------------------------------------
   stop
!
   end subroutine takeStopAction
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine takePrintAction(id,message,print_pe)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
!
   implicit none
!
   character (len=*), intent(in) :: message
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in), optional :: print_pe
   integer (kind=IntKind) :: i
!
   type(CheckPointStruct), pointer :: ptr
!
   if (id < 1 .or. id > NumPrintPoints) then
      return
   endif
!
   ptr => PrintPoint
   LOOP_i: do i = 1, NumPrintPoints
      if (id == ptr%index) then
         exit LOOP_i
      endif
      ptr => ptr%next
   enddo LOOP_i
!
   if (present(print_pe)) then
      if (print_pe == MyPE .or. print_pe == -1) then
         write(6,'(3a,i5,2a)')'Routine name: ',trim(ptr%rname),'; Process: ',MyPE, &
                              '; Message: ',trim(message)
      endif
   else
      write(6,'(4a)')'Routine name: ',trim(ptr%rname),'; Message: ',trim(message)
   endif
!  -------------------------------------------------------------------
   call syncAllPEs()
!  -------------------------------------------------------------------
!
   end subroutine takePrintAction
!  ===================================================================
end module CheckPointModule
