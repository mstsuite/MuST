!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   program chgaver
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
   character (len=10) :: atom
   character (len=1) :: dummy
!
   logical :: ElementName, Jump1stLine
!
   integer (kind=IntKind) :: num_atoms
   integer (kind=IntKind), parameter :: num_bins = 400
!
   integer (kind=IntKind), parameter :: funit = 5
   integer (kind=IntKind) :: ios, status, n, alen, i, j, k, m
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: dist_freq(:)
   integer (kind=IntKind) :: DataForm
   integer (kind=IntKind) :: StartingLine, NB
!
   real (kind=RealKind) :: CenterX, CenterY, CenterZ
   real (kind=RealKind), allocatable :: AtomPositionX(:)
   real (kind=RealKind), allocatable :: AtomPositionY(:)
   real (kind=RealKind), allocatable :: AtomPositionZ(:)
   real (kind=RealKind), allocatable :: rdist(:)
   real (kind=RealKind), allocatable :: chg(:), mom(:)
   real (kind=RealKind) :: bin_size, bin_half
   real (kind=RealKind) :: side_len, r, dq, moment
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
   side_len = ZERO
   do i=1,3
      r = sqrt( BravaisLatticeVec(1,i)*BravaisLatticeVec(1,i)         &
               +BravaisLatticeVec(2,i)*BravaisLatticeVec(2,i)         &
               +BravaisLatticeVec(3,i)*BravaisLatticeVec(3,i) )
      side_len = max(side_len,r)
   enddo
   bin_size = side_len/real(num_bins,kind=RealKind)
   bin_half = bin_size*HALF
!
   write(6,'(a,i6)')'Total Number of Atoms: ',num_atoms
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
print *,'CenterX = ',CenterX
print *,'CenterY = ',CenterY
print *,'CenterZ = ',CenterZ
!
   open(unit=10,file='ChargeTable',status='old')
   read(10,'(a)')dummy
   read(10,'(a)')dummy
   read(10,'(a)')dummy
!
   allocate( rdist(num_bins+1), dist_freq(num_bins+1) )
   allocate( chg(num_bins+1), mom(num_bins+1) )
   do i = 1, num_bins+1
      rdist(i) = bin_size*(i-1)
   enddo
   chg(1:num_bins+1) = ZERO
   mom(1:num_bins+1) = ZERO
   dist_freq(1:num_bins+1) = 0
!
   n = 0
   LOOP_i: do i=1,num_atoms
      r = sqrt( (AtomPositionX(i) - CenterX)**2                       &
               +(AtomPositionY(i) - CenterY)**2                       &
               +(AtomPositionZ(i) - CenterZ)**2 )
      read(10,'(48x,f9.5,2x,f9.5)')dq,moment
      LOOP_j: do j = 1, num_bins+1
         if (abs(r-rdist(j)) <= bin_half) then
            dist_freq(j) = dist_freq(j) + 1
            chg(j) = chg(j) + dq
            mom(j) = mom(j) + moment
            n = max(n,j)
            exit LOOP_j
         endif
      enddo LOOP_j
   enddo LOOP_i
!
   close(10)
!
   write(6,'(66(''=''))')
   write(6,'(a)')                  &
'  Bin Index       Distance(A)      Excess Electron         Moment'
   write(6,'(66(''-''))')
   do i = 1, n
      if (dist_freq(i) > 0) then
         write(6,'(3x,i5,3(9x,f10.5))')i,rdist(i)/anstr2au,           &
                         chg(i)/real(dist_freq(i),kind=RealKind),     &
                         mom(i)/real(dist_freq(i),kind=RealKind)
      endif
   enddo
!
   deallocate( rdist, dist_freq, chg, mom )
   deallocate( AtomicNumber, AtomPositionX, AtomPositionY, AtomPositionZ )
!
   end program chgaver
!  ===================================================================
