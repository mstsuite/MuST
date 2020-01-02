subroutine FlushFile(fu)
   use KindParamModule, only : IntKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: fu
!
   logical :: FileExist
!
#ifdef NoFLUSH
   return
#else
   inquire(unit=fu,exist=FileExist)
   if (FileExist) then
      flush(fu)
!     call flush(fu)
   endif
#endif
end subroutine FlushFile
