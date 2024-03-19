!  *******************************************************************
!
!  Given string s, return the number of tokens in the string.
!  Tokens are separated by space, comma, semicolumn, or Tab.
!
!  Created on March 17nd, 2024
!     --- by Yang Wang
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumTokens(s) result(NumTokens)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: s
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
      endif
      wasSep = isSep
   enddo LOOP_i
!
   end function getNumTokens
