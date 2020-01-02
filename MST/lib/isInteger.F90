!  *******************************************************************
!
!  If string s is a integer number, return true, otherwise false.
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isInteger(s) result(t)
!  ===================================================================
   implicit none
!
   logical :: t
!
   character (len=*), intent(in) :: s
!
   character (len=len_trim(s)) :: st
   character (len=1) :: c
!
   integer :: i, j, n, a0, a9, L
!
   st = adjustl(s); L = len_trim(st)
!
   if (st(1:1) == '+' .or. st(1:1) == '-') then
      n = 2
   else
      n = 1
   endif
!
   if (n > L) then
      t = .false.
      return
   else
      t = .true.
   endif
!
   a0 = iachar('0'); a9 = iachar('9')
   do i = n, L
      j = iachar(st(i:i))
      if (j < a0 .or. j > a9) then
         t = .false.
         exit
      endif
   enddo
!
   end function isInteger
