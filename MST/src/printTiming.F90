subroutine printTiming()
   use KindParamModule, only : IntKind, RealKind
   use TimerModule, only : getAccumelatedTiming
   use ScfDataModule, only : isRelativisticValence, isLSMS
!
   implicit none
!
   integer (kind=IntKind) :: num_calls
!
   real (kind=RealKind) :: time_acc, time_aver
!
   if (isRelativisticValence()) then
!     ----------------------------------------------------------------
      time_acc = getAccumelatedTiming(tag='SingleDiracScattering',aver=time_aver,nc=num_calls)
!     ----------------------------------------------------------------
      if (num_calls > 0) then
         write(6,'(a,t60,a,f12.3,a)')'Accumulated time in calling SingleDiracScattering',':',time_acc,' Sec'
         write(6,'(a,t60,a,f12.3,a,/)')'Averaged time in calling SingleDiracScattering',':',time_aver,' Sec'
      endif
   else
!     ----------------------------------------------------------------
      time_acc = getAccumelatedTiming(tag='solveSingleScattering',aver=time_aver,nc=num_calls)
!     ----------------------------------------------------------------
      if (num_calls > 0) then
         write(6,'(a,t60,a,f12.3,a)')'Accumulated time in calling solveSingleScattering',':',time_acc,' Sec'
         write(6,'(a,t60,a,f12.3,a,/)')'Averaged time in calling solveSingleScattering',':',time_aver,' Sec'
      endif
   endif
   if ( isLSMS() ) then
!     ----------------------------------------------------------------
      time_acc = getAccumelatedTiming(tag='Calculate the Gij matrix',aver=time_aver,nc=num_calls)
!     ----------------------------------------------------------------
      if (num_calls > 0) then
         write(6,'(a,t60,a,f12.3,a)')'Accumulated time for calculating Gij matrix',':',time_acc,' Sec'
         write(6,'(a,t60,a,f12.3,a,/)')'Averaged time for calculating Gij matrix',':',time_aver,' Sec'
      endif
!     ----------------------------------------------------------------
      time_acc = getAccumelatedTiming(tag='Construct the LIZ KKR matrix',aver=time_aver,nc=num_calls)
!     ----------------------------------------------------------------
      if (num_calls > 0) then
         write(6,'(a,t60,a,f12.3,a)')'Accumulated time for calculating LIZ KKR matrix',':',time_acc,' Sec'
         write(6,'(a,t60,a,f12.3,a,/)')'Averaged time for calculating LIZ KKR matrix',':',time_aver,' Sec'
      endif
!     ----------------------------------------------------------------
      time_acc = getAccumelatedTiming(tag='Copy the LIZ KKR matrix to GPU',aver=time_aver,nc=num_calls)
!     ----------------------------------------------------------------
      if (num_calls > 0) then
         write(6,'(a,t60,a,f12.3,a)')'Accumulated time for copying LIZ KKR matrix to GPU',':',time_acc,' Sec'
         write(6,'(a,t60,a,f12.3,a,/)')'Averaged time for copying LIZ KKR matrix to GPU',':',time_aver,' Sec'
      endif
!     ----------------------------------------------------------------
      time_acc = getAccumelatedTiming(tag='Invert the LIZ KKR matrix',aver=time_aver,nc=num_calls)
!     ----------------------------------------------------------------
      if (num_calls > 0) then
         write(6,'(a,t60,a,f12.3,a)')'Accumulated time for inverting LIZ KKR matrix',':',time_acc,' Sec'
         write(6,'(a,t60,a,f12.3,a,/)')'Averaged time for inverting LIZ KKR matrix',':',time_aver,' Sec'
      endif
   endif
!  -------------------------------------------------------------------
   time_acc = getAccumelatedTiming(tag='Calculate MS DOS',aver=time_aver,nc=num_calls)
!  -------------------------------------------------------------------
   if (num_calls > 0) then
      write(6,'(a,t60,a,f12.3,a)')'Accumulated time for calculating MS DOS',':',time_acc,' Sec'
      write(6,'(a,t60,a,f12.3,a)')'Averaged time for calculating MS DOS',':',time_aver,' Sec'
   endif
!
   write(6,'(80(''-''))')
!
end subroutine printTiming
