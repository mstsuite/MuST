program testTokens
   implicit none
!
   character (len=80) :: text
!
   integer :: status, i, n
!
   interface 
      function getNumTokens(s) result (n)
         character (len=*), intent(in) :: s
      end function getNumtokens
   end interface
!
   interface 
      function getToken(k,s,n,e) result (t)
         character (len=*), intent(in) :: s
         character (len=len(s)) :: t
         integer, intent(in) :: k
         integer, intent(out), optional :: n, e
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
