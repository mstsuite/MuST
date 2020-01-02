!  ********************************************************************
!  * MODULE NAME    : StringModule                                    *
!  *                                                                  *
!  * VERSION NUMBER : 1.0                                             *
!  * LAST MODIFIED  : MARCH 09, 2004                                  *
!  *                                                                  *
!  * DESCRIPTION: MODULE for string utility functions                 *
!  *                                                                  *
!  * EXTERNAL MODULE DEPENDENCE                                       *
!  * ================================================================ *
!  *    KindParamModule                                               *
!  *    ErrorHandlerModule                                            *
!  *                                                                  *
!  * PUBLIC FUNCTIONS                                                 *
!  * ================================================================ *
!  *    subroutine initString(str)                                    *
!  *    Purpose: initialize the module for the given string.          *
!  *    Input:   str   = character string                             *
!  *    Output:  none                                                 *
!  *    ============================================================= *
!  *    subroutine initString(str_len)                                *
!  *    Purpose: initialize the module for the given string length    *
!  *    Input:   str_len = character string length                    *
!  *    Output:  none                                                 *
!  *    ============================================================= *
!  *    subroutine endString()                                        *
!  *    Purpose: clean the memory allocated within the module.        *
!  *    Input:   none                                                 *
!  *    Output:  none                                                 *
!  *    ============================================================= *
!  *    subroutine setString(str)                                     *
!  *    Purpose: set the string to be str.                            *
!  *    Input:   str   = character string                             *
!  *    Output:  none                                                 *
!  *    ============================================================= *
!  *    function getNumTokens() result(n)                             *
!  *    Purpose: returns the number of tokens contained in the string *
!  *             Tokens are seperated by ' ', ',', ';'                *
!  *    Input:   none                                                 *
!  *    Output:  n     = integer, the number of tokens.               *
!  *    ============================================================= *
!  *    subroutine readToken(i,s,n)                                   *
!  *    Purpose: returns the i'th token and its length                *
!  *    Input:   i     = integer, i'th token                          *
!  *    Output:  s     = character string, token character string     *
!  *             n     = integer, the length of the token             *
!  *                                                                  *
!  ********************************************************************
module StringModule
   use KindParamModule, only : IntKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initString,    &
          setString,     &
          endString,     &
          getNumTokens,  &
          readToken
!
   interface initString
      module procedure initString_str, initString_len
   end interface
!
private
!
   logical :: Initialized = .false.
!
   character (len=1), allocatable :: c_str(:)
!
   integer (kind=IntKind) :: nstr
   integer (kind=IntKind) :: NumTokens
   integer (kind=IntKind), allocatable :: TokenIndex(:)
   integer (kind=IntKind), allocatable :: TokenLen(:)
!
   integer (kind=IntKind), parameter :: NumSeps = 4
!
   character (len=1), parameter :: sep(1:NumSeps) = (/' ',',',';',achar(9)/)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initString_str(str,separator,separators)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: str
   character (len=1), intent(in), optional :: separator
   character (len=1), intent(in), optional :: separators(:)
!
   if (Initialized) then
      call ErrorHandler('initString','String module is already initialized')
   endif
!
   nstr = len(str)
   allocate( c_str(1:nstr), TokenIndex(1:nstr), TokenLen(1:nstr) )
!
   if (present(separator)) then
      call processString(str,sp=separator)
   else if (present(separators)) then
      call processString(str,sps=separators)
   else
      call processString(str)
   endif
!
   Initialized = .true.
!
   end subroutine initString_str
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initString_len(l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: l
!
   if (Initialized) then
      call ErrorHandler('initString','String module is already initialized')
   endif
!
   nstr = l
   allocate( c_str(1:nstr), TokenIndex(1:nstr), TokenLen(1:nstr) )
!
   Initialized = .true.
!
   end subroutine initString_len
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setString(str,separator,separators)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: str
   character (len=1), intent(in), optional :: separator
   character (len=1), intent(in), optional :: separators(:)
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('setString','String module not initialized')
   endif
!
   n = len(str)
   if (n > nstr) then
      deallocate( c_str, TokenIndex, TokenLen )
      nstr = n
      allocate( c_str(1:nstr), TokenIndex(1:nstr), TokenLen(1:nstr) )
   else
      TokenIndex(1:nstr) = 0; TokenLen(1:nstr) = 0
   endif
!
   if (present(separator)) then
      call processString(str,sp=separator)
   else if (present(separators)) then
      call processString(str,sps=separators)
   else
      call processString(str)
   endif
!
   Initialized = .true.
!
   end subroutine setString
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine processString(str,sp,sps)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: str
   character (len=1), intent(in), optional :: sp
   character (len=1), intent(in), optional :: sps(:)
!
   logical :: isSep, wasSep
!
   integer (kind=IntKind) :: i, j, n, nsp
!
   if (present(sp)) then
      nsp = 0
   else if (present(sps)) then
      nsp = size(sps)
   else
      nsp = -1
   endif
!
   n = len_trim(str)
   do i=1,n
      c_str(i) = str(i:i)
   enddo
   do i=n+1,nstr
      c_str(i) = ' '
   enddo
!
   NumTokens = 0
   wasSep = .true.
   isSep = .false.
   LOOP_i: do i=1, n
      if (iachar(c_str(i)) < 32) then
         isSep = .true.
      else
         isSep = .false.
         LOOP_j: do j=1,NumSeps
            if (c_str(i) == sep(j)) then
               isSep = .true.
               exit LOOP_j
            endif
         enddo LOOP_j
         if (.not.isSep) then
            if (nsp == 0) then
               if (c_str(i) == sp) then
                  isSep = .true.
               endif
            else
               LOOP_jp: do j=1,nsp
                  if (c_str(i) == sps(j)) then
                     isSep = .true.
                     exit LOOP_jp
                  endif
               enddo LOOP_jp
            endif
         endif
      endif
      if (.not.isSep .and. wasSep) then
         NumTokens = NumTokens + 1
         TokenIndex(NumTokens) = i
      else if (isSep .and. .not.wasSep) then
         TokenLen(NumTokens) = i-TokenIndex(NumTokens)
      endif
      wasSep = isSep
   enddo LOOP_i
!
   if (.not.isSep .and. .not.wasSep) then
      TokenLen(NumTokens) = n-TokenIndex(NumTokens)+1
   endif
!
   end subroutine processString
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endString()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
      call ErrorHandler('endString','String module not initialized')
   endif
!
   deallocate( c_str, TokenIndex, TokenLen )
!
   nstr = 0
!
   Initialized = .false.
!
   end subroutine endString
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumTokens() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumTokens','String module not initialized')
   endif
!
   n = NumTokens
!
   end function getNumTokens
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readToken(i,s,nt)
!  ===================================================================
   implicit none
!
   character (len=*) :: s
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind), intent(out), optional :: nt
   integer (kind=IntKind) :: m, j, j0, n
!
   if (.not.Initialized) then
      call ErrorHandler('readToken','String module not initialized')
   else if (i < 1 .or. i > NumTokens) then
      call ErrorHandler('readToken','Token index is out of range',i)
   endif
!
   n = TokenLen(i)
   m = min(n, len(s))
   j0 = TokenIndex(i)-1
   do j=1,m
      s(j:j) = c_str(j0+j)
   enddo
   s(m+1:)=' '
   if (present(nt)) then
      nt = n
   endif
!
   end subroutine readToken
!  ===================================================================
end module StringModule
