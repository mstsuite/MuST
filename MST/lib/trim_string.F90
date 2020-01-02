function trim_string(s1) result(s2)
   implicit none
!
   character (len=*), intent(in) :: s1
   character (len=len(s1)) :: s2
   integer :: n, i, m
!
   n = len(trim(s1))
   m = 0
   LOOP_i : do i = n, 1, -1
      if (iachar(s1(i:i)) > 32) then
         m = i
         exit LOOP_i
      endif
   enddo LOOP_i
   if (m == 0) then
      s2 = ''
   else
      s2 = s1(1:m)
   endif
!
end function trim_string
