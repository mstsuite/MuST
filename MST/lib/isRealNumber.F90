!  *******************************************************************
!
!  If string s is a real number, return true, otherwise false.
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isRealNumber(s) result(t)
!  ===================================================================
   implicit none
!
   logical :: t
   logical :: NoPoint
   logical :: isExp
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
   t = .true.
   NoPoint = .true.
   isExp = .false.
   a0 = iachar('0'); a9 = iachar('9')
   do i = n, L
      j = iachar(st(i:i))
      if (st(i:i) == '.' .and. NoPoint) then
         NoPoint = .false.
      else if (st(i:i) == 'd' .or. st(i:i) == 'D' .or. st(i:i) == 'e' .or.  &
          st(i:i) == 'E') then
         if (i > n .and. .not.isExp) then
            isExp = .true.
            NoPoint = .false.
         else
            t = .false.
            exit
         endif
      else if (st(i:i) == '+' .or. st(i:i) == '-') then
         c = st(i-1:i-1)
         if (c /= 'd' .and. c /= 'D' .and. c /= 'e' .and. c /= 'E') then
            t = .false.
            exit
         endif
      else if (j < a0 .or. j > a9) then
         t = .false.
         exit
      endif
   enddo
!
   if (NoPoint) then
      t = .false.
   endif
!
   end function isRealNumber
