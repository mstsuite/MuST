function trim_string(s1,c) result(s2)
   implicit none
!
   character (len=*), intent(in) :: s1
   character (len=len(s1)) :: s2
   character (len=1), optional, intent(in) :: c
   integer :: n, i, m
!
   n = len(trim(s1))
!
   if (present(c)) then
      m = n
      LOOP_j : do i = 1, n
         if (s1(i:i) == c) then
            m = i-1
            exit LOOP_j
         endif
      enddo LOOP_j
   else
      m = 0
      LOOP_i : do i = n, 1, -1
         if (iachar(s1(i:i)) > 32) then
            m = i
            exit LOOP_i
         endif
      enddo LOOP_i
   endif
!
   if (m == 0) then
      s2 = ''
   else
      s2 = s1(1:m)
   endif
!
end function trim_string
