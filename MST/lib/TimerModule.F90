module TimerModule
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, TEN2m8
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initTimer, &
          endTimer,  &
          getTime,   &
          storeTime, &
          fetchStoredTime, &
          registerRoutine, &      ! register a routine to be recorded for its number of calls and timing.
          checkinTiming,   &      ! record the timing of a routine call. The timing will be accumulated.
          getRoutineCallTiming, & ! get the accumulated timing of a routine call.
          resetRoutineCallTiming  ! restart the timing of a routine call.
!
private
   character (len=8) :: date
   character (len=10) :: time
!
   integer (kind=IntKind) :: LastCount=0
   integer (kind=IntKind) :: OverFlow=0
   integer (kind=IntKind) :: hours, minutes
   integer (kind=IntKind) :: c, r, m
   integer (kind=IntKind) :: NumChecks
!
   integer (kind=IntKind) :: InitialYear
   integer (kind=IntKind) :: InitialMonth
   integer (kind=IntKind) :: InitialDay
   integer (kind=IntKind) :: InitialHour
   integer (kind=IntKind) :: InitialMinutes
   integer (kind=IntKind) :: InitialSeconds
   integer (kind=IntKind) :: InitialMilliSeconds
!
   real (kind=RealKind) :: InitialTime
   real (kind=RealKind) :: ClockTime
   real (kind=RealKind) :: MaxSec = -1.0
   real (kind=RealKind) :: seconds
   real (kind=RealKind) :: StoredTime
!
   type TimingStruct
      integer (kind=IntKind) :: index
      integer (kind=IntKind) :: num_calls
      real (kind=RealKind) :: accumulated_timing
      character (len=50) :: routine_name
      type (TimingStruct), pointer :: next
   end type TimingStruct
!
   type (TimingStruct), pointer :: check_list
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initTimer()
!  ===================================================================
   implicit none
   character (len=5) :: zone
!
   integer :: values(8)
   real (kind=RealKind) :: s
!
!  -------------------------------------------------------------------
   call date_and_time(date,time,zone,values)
   call system_clock(c, r, m)
!  -------------------------------------------------------------------
   if (m <= 0) then
      stop 'ERROR: invalid timer'
   endif
   read(time,'(i2,i2,f6.3)')hours,minutes,seconds
   ClockTime=hours*3.6d+3+minutes*6.0d+1+seconds
   LastCount=c
   s=r
   MaxSec=m/s
!  InitialTime=c/s
   InitialTime = ClockTime
   StoredTime = ZERO
!
   InitialYear = values(1)
   InitialMonth = values(2)
   InitialDay = values(3)
   InitialHour = values(5)
   InitialMinutes = values(6)
   InitialSeconds = values(7)
   InitialMilliSeconds = values(8)
!
   NumChecks = 0
   nullify(check_list)
!
   end subroutine initTimer
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endTimer()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i
!
   type (TimingStruct), pointer :: p
!
   do i = 1, NumChecks
      p => check_list
      check_list => p%next
      nullify(p%next)
      deallocate(p)
   enddo
   NumChecks = 0
   nullify(check_list)
!
   InitialTime = ZERO
   StoredTime = ZERO
!
   end subroutine endTimer
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine registerRoutine(func)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: func
!
   integer (kind=IntKind) :: i
!
   type (TimingStruct), pointer :: p
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (NumChecks == 0) then
      call createACheck(check_list,func)
   else
      p => check_list
      do i = 1, NumChecks
         if (nocaseCompare(func,p%routine_name)) then
            call WarningHandler('registerRoutine','The routine is already in check list',func)
            return
         endif
         if (i < NumChecks) then
            p => p%next
         endif
      enddo
      call createACheck(p%next,func)
   endif
!
   end subroutine registerRoutine
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createACheck(p,func)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: func
!
   type (TimingStruct), pointer :: p
!
   allocate(p)
!
   if (len_trim(adjustl(func)) > len(p%routine_name)) then
      call ErrorHandler('createACheck','The routine name exceeds the length limit',func)
   endif
!
   NumChecks = NumChecks + 1
!
   p%index = NumChecks
   p%routine_name = adjustl(func)
   p%num_calls = 0
   p%accumulated_timing = ZERO
   nullify(p%next)
!
   end subroutine createACheck
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkinTiming(func,t,nc)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: func
!
   real (kind=RealKind), intent(in) :: t
!
   integer (kind=IntKind), intent(in), optional :: nc
   integer (kind=IntKind) :: i, n
!
   type (TimingStruct), pointer :: p
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (present(nc)) then
      if (nc < 1) then
         call ErrorHandler('checkinTiming','Invalid number of routine calls',nc)
      else
         n = nc
      endif
   else
      n = 1
   endif
!
   if (t < TEN2m8) then
      call ErrorHandler('checkinTiming','The timing value < 10^(-8)',t)
   endif
!
   p => check_list
   LOOP_i: do i = 1, NumChecks
      if (nocaseCompare(func,p%routine_name)) then
         exit LOOP_i
      endif
      if (i < NumChecks) then
         p => p%next
      else
         call ErrorHandler('checkinTiming','The routine is not in the check list',func)
      endif
   enddo LOOP_i
!
   p%accumulated_timing = p%accumulated_timing + t
   p%num_calls = p%num_calls + n
!
   end subroutine checkinTiming
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRoutineCallTiming(func,nc) result(t)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: func
!
   integer (kind=IntKind), intent(out), optional :: nc
   integer (kind=IntKind) :: i
!
   real (kind=RealKind) :: t
!
   type (TimingStruct), pointer :: p
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   p => check_list
   LOOP_i: do i = 1, NumChecks
      if (nocaseCompare(func,p%routine_name)) then
         exit LOOP_i
      endif
      if (i < NumChecks) then
         p => p%next
      else
         call ErrorHandler('getRoutineCallTiming','The routine is not in the check list',func)
      endif
   enddo LOOP_i
!
   t = p%accumulated_timing
!
   if (present(nc)) then
      nc = p%num_calls
   endif
!
   end function getRoutineCallTiming
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine resetRoutineCallTiming(func,t,n)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: func
!
   integer (kind=IntKind), intent(out), optional :: n
   integer (kind=IntKind) :: i, n0
!
   real (kind=RealKind), intent(in), optional :: t
   real (kind=RealKind) :: t0
!
   type (TimingStruct), pointer :: p
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (present(t)) then
      t0 = t
   else
      t0 = ZERO
   endif
!
   if (present(n)) then
      n0 = n
   else
      n0 = 0
   endif
!
   p => check_list
   LOOP_i: do i = 1, NumChecks
      if (nocaseCompare(func,p%routine_name)) then
         exit LOOP_i
      endif
      if (i < NumChecks) then
         p => p%next
      else
         call ErrorHandler('resetRoutineCallTiming','The routine is not in the check list',func)
      endif
   enddo LOOP_i
!
   p%accumulated_timing = t0
   p%num_calls = n0
!
   end subroutine resetRoutineCallTiming
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTime() result(s)
!  ===================================================================
   implicit none
!  *******************************************************************
!  returns time in the unit of seconds
!  *******************************************************************
   character (len=5) :: zone
!
   integer :: values(8)
!
   real (kind=RealKind) :: newt
   real (kind=RealKind) :: s
!
!  -------------------------------------------------------------------
   call date_and_time(date,time,zone,values)
   call system_clock(c, r, m)
!  -------------------------------------------------------------------
!  print *,'time = ',time
!  print *,'values = ',values(1:8)
!  read(time,'(i2,i2,f6.3)')hours,minutes,seconds
!
   hours = values(5)
   minutes = values(6)
   seconds = values(7)+0.001d0*values(8)
   newt=hours*3.6d+3+minutes*6.0d+1+seconds
!      
   if (MaxSec < 0.0d0) then
      stop 'ERROR: need to call initTimer() first'
   else if (newt-ClockTime > MaxSec) then
      OverFlow=OverFlow+floor((newt-ClockTime)/MaxSec)
   endif
!
   LastCount=c
!! s=r
!! s=(c+m*OverFLow)/s
!
   ClockTime=newt
!
   if (values(3) - InitialDay == 0) then
      s = ClockTime - InitialTime
   else
      s = ClockTime - InitialTime + 24.0d0*3.6d+3
   endif
!
   end function getTime
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine storeTime(t)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: t
!
   StoredTime = t
!
   end subroutine storeTime
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function fetchStoredTime() result(t)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: t
!
   t = StoredTime
!
   end function fetchStoredTime
!  ===================================================================
end module TimerModule
