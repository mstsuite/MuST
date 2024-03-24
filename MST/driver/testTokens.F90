program testTokens
   use kindParamModule, only : IntKind
!
   implicit none
!
   character (len=80) :: text
!
   integer (kind=IntKind) :: status, i, n
!
   interface 
      function getNumTokens(s) result (n)
         use kindParamModule, only : IntKind
         character (len=*), intent(in) :: s
         integer (kind=IntKind) :: n
      end function getNumtokens
   end interface
!
   interface 
      function getToken(k,s,n,e) result (t)
         use kindParamModule, only : IntKind
         character (len=*), intent(in) :: s
         character (len=len(s)) :: t
         integer (kind=IntKind), intent(in) :: k
         integer (kind=IntKind), intent(out), optional :: n, e
      end function getToken
   end interface
!
   status = 0
   do while (status == 0) 
      read(5,'(a)',iostat=status)text
      print *,'String: ',adjustl(trim(text))
      n = getNumTokens(text)
      print *,'        Number of token = ', n
      do i = 1, n
         print *,'        ',trim(getToken(i,text))
      enddo
   enddo
!
end program testTokens
