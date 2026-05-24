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
          storeTime, &            ! store the timing for a Timing Tag
          fetchStoredTime, &      ! fetch the timing associated with the stored Timing Tag
          registerTimingTag, &    ! register a Timing-Tag to be recorded for its number of calls and its timing.
          checkinTiming,   &      ! record the timing of a Timing-Tag. The timing will be accumulated.
          getAccumelatedTiming, & ! get the accumulated timing of a Timing-Tag
          resetTimingTag          ! restart the timing of a Timing Tag.
!
private
   character (len=8) :: date
   character (len=10) :: time
!
   integer (kind=IntKind) :: LastCount=0
   integer (kind=IntKind) :: OverFlow=0
   integer (kind=IntKind) :: hours, minutes
   integer (kind=IntKind) :: c, r, m
   integer (kind=IntKind) :: NumStores
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
!
   type StoredTimeStruct
      integer (kind=IntKind) :: index
      real (kind=RealKind) :: stored_time
      character (len=80) :: tag_name
      type (StoredTimeStruct), pointer :: next
   end type StoredTimeStruct
!
   type TimingStruct
      integer (kind=IntKind) :: index
      integer (kind=IntKind) :: num_calls
      real (kind=RealKind) :: accumulated_timing
      character (len=80) :: tag_name
      type (TimingStruct), pointer :: next
   end type TimingStruct
!
   type (StoredTimeStruct), pointer :: store_list
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
!
   InitialYear = values(1)
   InitialMonth = values(2)
   InitialDay = values(3)
   InitialHour = values(5)
   InitialMinutes = values(6)
   InitialSeconds = values(7)
   InitialMilliSeconds = values(8)
!
   NumStores = 0
   NumChecks = 0
   nullify(store_list)
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
   type (StoredTimeStruct), pointer :: ps
!
   type (TimingStruct), pointer :: p
!
   do i = 1, NumStores
      ps => store_list
      store_list => ps%next
      nullify(ps%next)
      deallocate(ps)
   enddo
   NumStores = 0
   nullify(store_list)
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
!
   end subroutine endTimer
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine registerTimingTag(tag)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tag
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
      call createACheck(check_list,tag)
   else
      p => check_list
      do i = 1, NumChecks
         if (nocaseCompare(tag,p%tag_name)) then
            call WarningHandler('registerTimingTag','The tag is already in check list',tag)
            return
         endif
         if (i < NumChecks) then
            p => p%next
         endif
      enddo
      call createACheck(p%next,tag)
   endif
!
   end subroutine registerTimingTag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createAStore(p,tag)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tag
!
   type (StoredTimeStruct), pointer :: p
!
   allocate(p)
!
   if (len_trim(adjustl(tag)) > len(p%tag_name)) then
      call ErrorHandler('createAStore','The tag name exceeds the length limit',tag)
   endif
!
   NumStores = NumStores + 1
!
   p%index = NumStores
   p%tag_name = adjustl(tag)
   p%stored_time = ZERO
   nullify(p%next)
!
   end subroutine createAStore
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createACheck(p,tag)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tag
!
   type (TimingStruct), pointer :: p
!
   allocate(p)
!
   if (len_trim(adjustl(tag)) > len(p%tag_name)) then
      call ErrorHandler('createACheck','The tag name exceeds the length limit',tag)
   endif
!
   NumChecks = NumChecks + 1
!
   p%index = NumChecks
   p%tag_name = adjustl(tag)
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
   subroutine checkinTiming(tag,t,nc)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tag
   character (len=100) :: text
!
   real (kind=RealKind), intent(in) :: t
   real (kind=RealKind) :: timing
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
         call ErrorHandler('checkinTiming','Invalid number of tag calls',nc)
      else
         n = nc
      endif
   else
      n = 1
   endif
!
   if (t < TEN2m8) then
      text = 'Timing Tag: '//tag//'. The input timing value < 10^(-8). It is reset to 0'
      call WarningHandler('checkinTiming',trim(text),t)
      timing = ZERO
   else 
      timing = t
   endif
!
   p => check_list
   LOOP_i: do i = 1, NumChecks
      if (nocaseCompare(tag,p%tag_name)) then
         exit LOOP_i
      endif
      if (i < NumChecks) then
         p => p%next
      else
         call ErrorHandler('checkinTiming','The tag is not in the check list',tag)
      endif
   enddo LOOP_i
!
   p%accumulated_timing = p%accumulated_timing + timing
   p%num_calls = p%num_calls + n
!
   end subroutine checkinTiming
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAccumelatedTiming(tag,aver,nc) result(t)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tag
!
   integer (kind=IntKind), intent(out), optional :: nc
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(out), optional :: aver
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
      if (nocaseCompare(tag,p%tag_name)) then
         exit LOOP_i
      endif
      if (i < NumChecks) then
         p => p%next
      else
         call ErrorHandler('getAccumelatedTiming','The tag is not in the check list',tag)
      endif
   enddo LOOP_i
!
   t = p%accumulated_timing
!
   if (present(aver)) then
      if (p%num_calls == 0) then
         call WarningHandler('getAccumelatedTiming','The tag has not been called',tag)
         aver = ZERO
      else
         aver = t/real(p%num_calls,kind=RealKind)
      endif
   endif
!
   if (present(nc)) then
      nc = p%num_calls
   endif
!
   end function getAccumelatedTiming
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine resetTimingTag(tag,t,n)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tag
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
      if (nocaseCompare(tag,p%tag_name)) then
         exit LOOP_i
      endif
      if (i < NumChecks) then
         p => p%next
      else
         call ErrorHandler('resetTimingTag','The tag is not in the check list',tag)
      endif
   enddo LOOP_i
!
   p%accumulated_timing = t0
   p%num_calls = n0
!
   end subroutine resetTimingTag
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
   subroutine storeTime(tag,t)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tag
!
   real (kind=RealKind), intent(in) :: t
!
   integer (kind=IntKind) :: i
!
   type (StoredTimeStruct), pointer :: p
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (NumStores == 0) then
      call createAStore(store_list,tag)
      p => store_list
   else
      p => store_list
      do i = 1, NumStores
         if (nocaseCompare(tag,p%tag_name)) then
            call WarningHandler('storeTime','The tag is already in store list',tag)
            return
         endif
         if (i < NumStores) then
            p => p%next
         endif
      enddo
      call createAStore(p%next,tag)
      p => p%next
   endif
!
   p%stored_time = t
!
   end subroutine storeTime
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function fetchStoredTime(tag) result(t)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tag
!
   real (kind=RealKind) :: t
!
   integer (kind=IntKind) :: i
!
   type (StoredTimeStruct), pointer :: p
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   p => store_list
   LOOP_i: do i = 1, NumStores
      if (nocaseCompare(tag,p%tag_name)) then
         exit LOOP_i
      endif
      if (i < NumStores) then
         p => p%next
      else
         call ErrorHandler('fetchStoredTime','The tag is not in the store list',tag)
      endif
   enddo LOOP_i
!
   t = p%stored_time
!
   end function fetchStoredTime
!  ===================================================================
end module TimerModule
