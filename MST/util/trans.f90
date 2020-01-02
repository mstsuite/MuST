program trans
   implicit none
!
   integer :: zed, status
!
   real*8 :: x, y, z
!
   LOOP_DO: do
      read(5,*,iostat=status)zed,x,y,z
      if (status < 0) then
         exit LOOP_DO
      endif
      write(6,'(i5,3f12.5,a)')zed,x,y,z,' '
   enddo LOOP_DO
   stop 'Ok'
end program trans
