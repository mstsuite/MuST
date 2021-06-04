!  ==================================================================
!
!  Purpose: Breaks a string into a text in multiple lines.
!           s  ( input)  : String to be broken into multiple lines
!           lw  (input)  : The width of each line
!           nl  (output) : The number of lines the string is broken into
!           sml (output) : The multi-line strings.
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine breakLine(s,lw,nl,sml)
!  ==================================================================
   use KindParamModule, only : IntKind
!
   integer (kind=IntKind), intent(in) :: lw
   integer (kind=IntKind), intent(out) :: nl
   integer (kind=IntKind) :: n, i, m, m0
!
   character (len=*), intent(in) :: s
   character (len=*), pointer, intent(out) :: sml(:)
!
   n = len(trim(s))
!
   if (lw < 1 .or. n < 1) then
      nl = 0
      return
   else
      nl = 1
      m0 = 0
      m = 0
      do while (n - m > lw)
         m0 = m
         m = m + lw
         do while (m > 1 .and. s(m:m) /= ' ')
            m = m - 1
         enddo
!        write(6,'(a)')s(m0+1:m)
         nl = nl + 1
      enddo
!     write(6,'(a)')s(m+1:n)
   endif
!
   if (associated(sml)) then
      if (size(sml) < nl) then
         deallocate(sml)
         allocate(sml(nl))
      else
         sml = ' '
      endif
   else
      allocate(sml(nl))
   endif
!
   nl = 1
   m0 = 0
   m = 0
   do while (n - m > lw)
      m0 = m
      m = m + lw
      do while (m > 1 .and. s(m:m) /= ' ')
         m = m - 1
      enddo
      sml(nl) = s(m0+1:m)
      nl = nl + 1
   enddo
   sml(nl) = s(m+1:n)
!
!  ==================================================================
   end subroutine breakLine
