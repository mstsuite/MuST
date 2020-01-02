program tst_string
   use KindParamModule, only : IntKind
   use StringModule, only : initString, endString, getNumTokens, readToken
!
   implicit none
!
   character (len=101) :: s1 = ' I am a student, and he is a student too.'
   character (len=87) :: s2 = 'I am not a student,   but he is a student.  '
   character (len=10) :: t
   integer (kind=IntKind) :: n, i, m
!
   call initString(s1)
   n = getNumTokens()
   print *,'Num of tokens = ',n
   do i=1,n
      call readToken(i,t,m)
      print *,'token = ',t(1:m),', length = ',m
   enddo
   call endString()
!
   call initString(s2)
   n = getNumTokens()
   print *,'Num of tokens = ',n
   do i=1,n
      call readToken(i,t,m)
      print *,'token = ',t(1:m),', length = ',m
   enddo
   call endString()
!
end program tst_string
