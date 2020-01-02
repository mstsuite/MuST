!  *******************************************************************
!
!  Given string s and the token index k, return token position in the
!  string. Tokens are separated by space, comma, semicolumn, or Tab.
!
!  if p = 0, the token position is not determined.
!
!  Updated on May 2nd, 2016: A Tab is added as a separator.
!     --- by Yang Wang
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTokenPosition(k,s,token_len) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: s
!
   integer, intent(in) :: k
   integer :: p, m
   integer, intent(out), optional :: token_len
!
   integer, parameter :: NumSeps = 4
!
   character (len=1), parameter :: sep(1:NumSeps) = (/' ', ',', ';',achar(9)/)
!
   logical :: isSep, wasSep
!
   integer :: i, j, n, NumTokens
!
   n = len_trim(s)
!
   p = 0
   m = n
   if (k < 1 .or. n < 1 .or. k > n-2) then
      if (present(token_len)) then
         token_len = m
      endif
      return
   endif
!
   NumTokens = 0
   wasSep = .true.
   LOOP_i: do i = 1, n
      isSep = .false.
      LOOP_j: do j=1,NumSeps
         if (s(i:i) == sep(j)) then
            isSep = .true.
            exit LOOP_j
         endif
      enddo LOOP_j
      if (.not.isSep .and. wasSep) then
         NumTokens = NumTokens + 1
         if (NumTokens == k) then
            p = i
!           return
            exit LOOP_i
         endif
      endif
      wasSep = isSep
   enddo LOOP_i
!
   if (present(token_len)) then
      m = 1
      LOOP_i1: do i = p+1, n
         isSep = .false.
         LOOP_j1: do j=1, NumSeps
            if (s(i:i) == sep(j)) then
               isSep = .true.
               exit LOOP_j1
            endif
         enddo LOOP_j1
         if (isSep) then
            m = i - p
            exit LOOP_i1
         else if (i == n) then
            m = n - p + 1
         endif
      enddo LOOP_i1
      token_len = m
   endif
!
   end function getTokenPosition
