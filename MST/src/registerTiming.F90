subroutine registerTiming()
   use KindParamModule, only : IntKind, RealKind
   use TimerModule, only : registerTimingTag
   use ScfDataModule, only : isRelativisticValence, isLSMS
   use CmdLineOptionModule, only : getCmdLineOption
!
   implicit none
!
   if (isRelativisticValence()) then
      call registerTimingTag(tag='SingleDiracScattering')
   else
      call registerTimingTag(tag='solveSingleScattering')
   endif
   if ( isLSMS() ) then
!     ----------------------------------------------------------------
      call registerTimingTag(tag='Calculate the Gij matrix')
      call registerTimingTag(tag='Construct the LIZ KKR matrix')
      call registerTimingTag(tag='Copy the LIZ KKR matrix to GPU')
      call registerTimingTag(tag='Invert the LIZ KKR matrix')
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call registerTimingTag(tag='Calculate MS DOS')
!  -------------------------------------------------------------------

end subroutine registerTiming
