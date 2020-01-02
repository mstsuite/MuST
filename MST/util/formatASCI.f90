program formatASCI
   character (len=8000) :: text
   integer :: status, l
   do
      read(5,'(a)',iostat=status)text
      if (status /= 0) then
         exit
      endif
      l = len( trim(text) )
      write(6,'(a)')text(1:l)
   enddo
end program formatASCI
