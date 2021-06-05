program tst_string
   use KindParamModule, only : IntKind
   use StringModule, only : initString, endString, getNumTokens, readToken
!
   implicit none
!
   integer (kind=IntKind), parameter :: line_width = 40
!
   character (len=2000) :: s1 = &
' So long as there shall exist, by reason of law &
and custom, a social condemnation, which, in the face of civilization, &
artificially creates hells on earth, and complicates a destiny that is divine &
with human fatality; so long as the three problems of the age--the degradation &
of man by poverty, the ruin of women by starvation, and the dwarfing of &
childhood by physical and spiritual night--are not solved; so long as, in &
certain regions, social asphyxia shall be possible; in other words, and from a &
yet more extended point of view, so long as ignorance and misery remain on &
earth, books like this cannot be useless.'
   character (len=100) :: s2 = 'Quote from Les Miserables  '
   character (len=50) :: t
   character (len=500) :: text = &
'The novel Anna Karenina begins with one of its most oft-quoted lines: &
Happy families are all alike; every unhappy family is unhappy in its own way. '
   integer (kind=IntKind) :: n, i, m, nl
   character (len=line_width), pointer :: sml(:)
!
   interface
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine breakLine(s,lw,nl,sml)
!     ================================================================
      use KindParamModule, only : IntKind
!
      integer (kind=IntKind), intent(in) :: lw
      integer (kind=IntKind), intent(out) :: nl
!
      character (len=*), intent(in) :: s
      character (len=*), pointer, intent(out) :: sml(:)
!
      end subroutine breakLine
!     ================================================================
   end interface
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
   call breakLine(s1,line_width,nl,sml)
   write(6,'(/)')
   do i = 1, nl
      write(6,'(a)')sml(i)
   enddo
!
   call breakLine(s2,line_width,nl,sml)
   write(6,'(/)')
   do i = 1, nl
      write(6,'(a)')sml(i)
   enddo
!
   call breakLine(text,line_width,nl,sml)
   write(6,'(/)')
   do i = 1, nl
      write(6,'(a)')sml(i)
   enddo
!
end program tst_string
