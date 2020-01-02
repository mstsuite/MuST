subroutine finishProcess()
   use KindParamModule, only : IntKind
!
   use MPPModule, only : endMPP, MyPE
!
   use GroupCommModule, only : endGroupComm
!
   use ParallelIOModule, only : endParallelIO
!
   use DataServiceCenterModule, only : endDataServiceCenter
!
   use InputModule, only : endInput
!
   use OutputModule, only : endOutput, getStandardOutputLevel
!
   use ScfDataModule, only : endScfData, isLSMS
!
   use PotentialTypeModule, only : endPotentialType
!
   use SystemModule, only : endSystem
!
   use SystemSymmetryModule, only : endSystemSymmetry
!
   use SystemVolumeModule, only : endSystemVolume
!
   use ProcMappingModule, only : endProcMapping
!
   use Atom2ProcModule, only : endAtom2Proc
!
   use AtomModule, only : endAtom
!
   use ContourModule, only : endContour
!
   use BZoneModule, only : endBZone
!
   implicit none
!
   character (len=8)  :: exec_date
   character (len=10) :: exec_time
!
   integer (kind=IntKind) :: my_id
!
   my_id = MyPE
!
   if (getStandardOutputLevel() >= 0) then
      write(6,*) " Start finishProcess! "
      call FlushFile(6)
   endif
!
   if (.not.isLSMS()) then
      call endBZone()
   endif
!
   call endAtom()
   call endOutput()
   call endSystemVolume()
   call endAtom2Proc()
   call endParallelIO()
   call endProcMapping()
   call endContour()
   call endSystem()
   call endPotentialType()
   call endScfData()
   call endInput()
   call endDataServiceCenter()
   call endGroupComm()
   call endMPP()
!
!  ===================================================================
   if (my_id == 0) then
      call date_and_time(exec_date,exec_time)
      write(6,'(/,12a)')'Execution ends at ',                         &
           exec_time(1:2),':',exec_time(3:4),':',exec_time(5:6),', ', &
           exec_date(5:6),'-',exec_date(7:8),'-',exec_date(1:4)
      write(6,'(80(''-''))')
      stop 'End of the program. OK!'
   else
      stop
   endif
!
end subroutine finishProcess
