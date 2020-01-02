!  *******************************************************************
!
!  Compare two character strings, s1 and s2, without considering their
!  case (upper or lower). If true, both strings are the same (case 
!  insensitive); if false, otherwise.
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function nocaseCompare(s1,s2) result(t)
!  ===================================================================
   implicit none
!
   logical :: t
!
   character (len=*), intent(in) :: s1
   character (len=*), intent(in) :: s2
!
   character (len=len_trim(s1)) :: s1t
   character (len=len_trim(s2)) :: s2t
!
   integer :: i, j, n, Ua, La
   integer :: l1, l2
!
   s1t = adjustl(s1); l1 = len_trim(s1t)
   s2t = adjustl(s2); l2 = len_trim(s2t)
!
   if (l1 /= l2) then
      t = .false.
      return
   else if (s1t == s2t) then
      t = .true.
      return
   endif
!
   n = iachar('a')-iachar('A')
   Ua = iachar('A'); La = iachar('a')
   do i=1,l1
      j = iachar(s1t(i:i))
      if (j >= Ua .and. j < La) then
         s1t(i:i) = achar(j+n)
      endif
      j = iachar(s2t(i:i))
      if (j >= Ua .and. j < La) then
         s2t(i:i) = achar(j+n)
      endif
   enddo
   if (s1t == s2t) then
      t = .true.
   else
      t = .false.
   endif
!
   end function nocaseCompare
