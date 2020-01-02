!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   program slice
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, HALF
!
   use ChemElementModule, only : getZtot
!
   use StringModule, only : initString, endString, getNumTokens, readToken
!
   implicit none
!
   character (len=80) :: text
   character (len=40) :: posfile
   character (len=40) :: chgfile
   character (len=10) :: atom
   character (len=10) :: plane
   character (len=80) :: dummy
!
   logical :: ElementName, Jump1stLine
!
   integer (kind=IntKind) :: num_atoms
!
   integer (kind=IntKind), parameter :: funit = 10
   integer (kind=IntKind) :: ios, status, n, alen, i, j, k, m, n1, n2, n3
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind) :: DataForm
   integer (kind=IntKind) :: StartingLine, NB
!
   real (kind=RealKind) :: CenterX, CenterY, CenterZ
   real (kind=RealKind), allocatable :: AtomPositionX(:)
   real (kind=RealKind), allocatable :: AtomPositionY(:)
   real (kind=RealKind), allocatable :: AtomPositionZ(:)
   real (kind=RealKind), allocatable :: rdist(:)
   real (kind=RealKind), allocatable :: chg(:), mom(:)
   real (kind=RealKind) :: r, dq, moment, d
   real (kind=RealKind), parameter :: anstr2au = 1.88973d0
!
   real (kind=RealKind) :: BravaisLatticeVec(3,3)
!
   interface
      function isNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isNumber
   end interface
!
   interface
      function getTokenPosition(k,s) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer :: p
      end function getTokenPosition
   end interface
!
   write(6,'(a,$)')' Position Data File: '
   read(5,'(a)')text
   posfile = adjustl(text)
   print *,'File name: ',trim(posfile)
!
   LOOP_plane: do
      write(6,'(a,$)')' Slicing Plane (e.g. [1,1,1], [1,-1,0], etc): '
      read(5,'(a)')text
      n = len(trim(adjustl(text)))
      if (n >= 7 .and. n <= 10) then
         plane = adjustl(text)
         if (plane(1:1) == '[' .and. plane(n:n) == ']') then
            read(plane(2:n-1),*)n1,n2,n3
!           call initString(plane(2:n-1))
!           if (getNumTokens() == 3) then
!              call readToken(1,a1,alen)
!              call readToken(2,a2,alen)
!              call readToken(3,a3,alen)
!           endif 
!           call endString()
            print *,'Plane Index: ',n1,n2,n3
            exit LOOP_plane
         endif 
      endif
      print *,'Invalid plane: ',trim(adjustl(text))
   enddo LOOP_plane
!
   open(unit=funit,file=posfile,status='old')
!
!  ===================================================================
!  Determine the format of the position data
!  ===================================================================
   DataForm = -1
   Jump1stLine = .false.
   ElementName = .false.
   StartingLine = 0
   NB = 0
   num_atoms = 0
   do
      read(funit,'(a)',iostat=status)text
      if (status < 0) then
         exit
      else if (num_atoms == 0) then
         StartingLine = StartingLine + 1
      endif
      text=adjustl(text)
      if (text(1:1) == '#' .or. text(1:1) == '!'                      &
                           .or. len_trim(text) == 0) then
         cycle
      else if (DataForm == -1) then
!        -------------------------------------------------------------
         call initString(text)
!        -------------------------------------------------------------
         n = getNumTokens()
!        -------------------------------------------------------------
         if (n == 1) then   ! the first line may only have 1 number
            if (.not.Jump1stLine) then
               Jump1stLine = .true.
            else
               print *,'Invalid line: ',trim(text)
               stop 'Error'
            endif
         else if (n == 3) then   ! 3-column data
            NB = NB + 1
            if (NB > 3) then
               print *,'Invalid line: ',trim(text)
               stop 'Error'
            else
               read(text,*,iostat=status)BravaisLatticeVec(1:3,NB)
               if (status > 0) then
                  print *,'Invalid line: ',trim(text)
                  stop 'Error'
               endif
            endif
         else 
            if (n == 4) then   ! 4-column data
               DataForm = 0
!              -------------------------------------------------------
               call readToken(1,atom,alen)
!              -------------------------------------------------------
            else if (n == 5) then   ! 5-column data
               DataForm = 1
!              -------------------------------------------------------
               call readToken(2,atom,alen)
!              -------------------------------------------------------
            else if (n > 5) then
               DataForm = 2        ! info_table format
!              -------------------------------------------------------
               call readToken(1,atom,alen)
!              -------------------------------------------------------
            else
               print *,'Invalid data form: ',trim(text)
               stop 'Error'
            endif
            if (isNumber(atom)) then
               ElementName = .false.
            else
               ElementName = .true.
            endif
            StartingLine = StartingLine - 1
            num_atoms = 1
         endif
!        -------------------------------------------------------------
         call endString()
!        -------------------------------------------------------------
      else
         num_atoms = num_atoms + 1
      endif
   enddo
!
   write(6,'(a,i6)')' Total Number of Atoms: ',num_atoms
!
   if (NB < 3 .or. StartingLine < 3) then
      print *,'Invalid position data'
      stop 'Error'
   endif
!
   allocate( AtomicNumber(num_atoms) )
   allocate( AtomPositionX(num_atoms) )
   allocate( AtomPositionY(num_atoms) )
   allocate( AtomPositionZ(num_atoms) )
!
   rewind(funit)
!
!  ===================================================================
!  StartingLine+1 is the line number of the first position data.
!  ===================================================================
   do n = 1, StartingLine
      read(funit,'(a)')dummy
   enddo
!
   n=0
   do
      read(funit,'(a)',iostat=status)text
      if (status < 0 .or. n == num_atoms) then
         exit
      endif
      text=adjustl(text)
      if (text(1:1) == '#' .or. text(1:1) == '!'                      &
                           .or. len_trim(text) == 0) then
         cycle
      else
         n = n+1
         if (DataForm == 0) then
             if (ElementName) then
                j = getTokenPosition(1,text)
                k = getTokenPosition(2,text)
                atom = trim(text(j:k-1))
                AtomicNumber(n) = getZtot(atom(1:3))
                read(text(k:),*)                                      &
                     AtomPositionX(n),AtomPositionY(n),AtomPositionZ(n)
             else
                read(text,*)AtomicNumber(n),                          &
                     AtomPositionX(n),AtomPositionY(n),AtomPositionZ(n)
             endif
         else if (DataForm == 1) then
             if (ElementName) then
                j = getTokenPosition(2,text)
                k = getTokenPosition(3,text)
                atom = trim(text(j:k-1))
                AtomicNumber(n) = getZtot(atom(1:3))
                read(text(k:),*)                                      &
                         AtomPositionX(n),AtomPositionY(n),AtomPositionZ(n)
             else
                read(text,*)k,AtomicNumber(n),                        &
                         AtomPositionX(n),AtomPositionY(n),AtomPositionZ(n)
             endif
         else if (DataForm == 2) then
             if (ElementName) then
                j = getTokenPosition(1,text)
                k = getTokenPosition(2,text)
                atom = trim(text(j:k-1))
                AtomicNumber(n) = getZtot(atom(1:3))
                k = getTokenPosition(3,text)
                read(text(k:),*)                                      &
                         AtomPositionX(n),AtomPositionY(n),AtomPositionZ(n)
             else
                read(text,*)AtomicNumber(n),k,                        &
                         AtomPositionX(n),AtomPositionY(n),AtomPositionZ(n)
             endif
         else
             print *,'Unkown data format'
             stop 'Error'
         endif
      endif
   enddo
!
   if (n < num_atoms) then
      print *,'Position data is incomplete'
      print *,'n = ',n
      print *,'num_atoms = ',num_atoms
      stop 'Error'
   endif
!
   close(unit=funit)
!
   write(6,'(/,a,/)')' The Position Data is read successfully!'
!
   CenterX = ZERO
   CenterY = ZERO
   CenterZ = ZERO
   do i=1,num_atoms
      CenterX = CenterX + AtomPositionX(i)/real(num_atoms,kind=RealKind)
   enddo
   do i=1,num_atoms
      CenterY = CenterY + AtomPositionY(i)/real(num_atoms,kind=RealKind)
   enddo
   do i=1,num_atoms
      CenterZ = CenterZ + AtomPositionZ(i)/real(num_atoms,kind=RealKind)
   enddo
   write(6,'(a,f15.8)')' CenterX = ',CenterX
   write(6,'(a,f15.8)')' CenterY = ',CenterY
   write(6,'(a,f15.8)')' CenterZ = ',CenterZ
!
!
   write(6,'(a,$)')' Charge Data File: '
   read(5,'(a)')text
   chgfile = adjustl(text)
   print *,'File name: ',trim(chgfile)
   open(unit=11,file=chgfile,status='old')
!
   read(11,'(a)')dummy
   read(11,'(a)')dummy
   read(11,'(a)')dummy
   read(11,'(a)')dummy
!
   allocate( rdist(num_atoms) )
   allocate( chg(num_atoms), mom(num_atoms) )
!
   CenterX = ZERO; CenterY = ZERO; CenterZ = ZERO
   LOOP_i: do i=1,num_atoms
      rdist(i) = sqrt( (AtomPositionX(i) - CenterX)**2                &
                      +(AtomPositionY(i) - CenterY)**2                &
                      +(AtomPositionZ(i) - CenterZ)**2 )
      read(11,'(48x,f9.5,2x,f9.5)')chg(i),mom(i)
   enddo LOOP_i
!
   close(unit=11)
!
   open(unit=12,file='Sliced',status='unknown')
!
   write(12,'(a)')plane
   write(12,'(80("="))')
   write(12,'(a)')                                                    &
'Atom         x                 y                 z           Charge    Moment'
   write(12,'(80("-"))')
   do i=1,num_atoms
      d =  n1*(AtomPositionX(i) - CenterX)                            &
         + n2*(AtomPositionY(i) - CenterY) + n3*(AtomPositionZ(i) - CenterZ)
      if (d <= ZERO) then
         write(12,'(i3,1x,3f18.12,2f10.5)')AtomicNumber(i),           &
               AtomPositionX(i),AtomPositionY(i),AtomPositionZ(i),chg(i),mom(i)
      endif
   enddo
!
   close(unit=12)
!
   deallocate( rdist, chg, mom )
   deallocate( AtomicNumber, AtomPositionX, AtomPositionY, AtomPositionZ )
!
   end program slice
!  ===================================================================
