!  *******************************************************************
!
!  Given string s and the token index k, return the token of the string.
!  Tokens are separated by space, comma, semicolumn, or Tab.
!
!  Created on Mar 18th, 2024
!     --- by Yang Wang
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getToken(k,s,token_len,err) result(t)
!  ===================================================================
   use KindParamModule, only : IntKind
!
   implicit none
!
   character (len=*), intent(in) :: s
   Character (len=len(s)) :: t
!
   integer (kind=IntKind), intent(in) :: k
   integer (kind=IntKind) :: p, m
   integer (kind=IntKind), intent(out), optional :: token_len, err
!
   integer (kind=IntKind), parameter :: NumSeps = 4
!
   character (len=1), parameter :: sep(1:NumSeps) = (/' ', ',', ';',achar(9)/)
!
   logical :: isSep, wasSep
!
   integer (kind=IntKind) :: i, j, n, NumTokens
!
   if (present(err)) then
      err = 0
   endif
!
   n = len_trim(s)
!
   p = 0
   m = n
   if (k < 1 .or. n < 1 .or. k > n) then
      if (present(token_len)) then
         token_len = m
      endif
      t = s
      if (present(err)) then
         err = 1
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
   t = s(p:p+m-1)
!
   if (present(token_len)) then
      token_len = m
   endif
!
   end function getToken
