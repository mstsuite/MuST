module TimerModule
   use KindParamModule, only : IntKind, RealKind
!
public :: initTimer, &
          getTime
!
private
   character (len=8) :: date
   character (len=10) :: time
!
   integer (kind=IntKind) :: LastCount=0
   integer (kind=IntKind) :: OverFlow=0
   integer (kind=IntKind) :: hours, minutes
   integer (kind=IntKind) :: c, r, m
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
   end subroutine initTimer
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
end module TimerModule
