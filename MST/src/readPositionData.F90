!  *******************************************************************
!  The position data are usually in the following format:
!
!  ! comments
!  # comments
!  scaling_factor
!  bravais_lattice(1,1)  bravais_lattice(2,1)  bravais_lattice(3,1)
!  bravais_lattice(1,2)  bravais_lattice(2,2)  bravais_lattice(3,2)
!  bravais_lattice(1,3)  bravais_lattice(2,3)  bravais_lattice(3,3)
!  D(irect) or d(irect) or C(artition) or c(artition) or not specified (C is default)
!      atom_name      x  y  z
!  ! or
!      CPA            x  y  z  atom_name_1      content_1  atom_name_2      content_2
!  ! or
!      atomic_number  x  y  z
!  ! or
!      -1             x  y  z  atomic_number_1  content_1  atomic_number_2  content_2
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readPositionData(fname,NumSitesIn,NumSitesOut)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ZERO, ONE
   use MPPModule, only : MyPE, bcastMessage
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use DataServiceCenterModule, only : createDataStorage,             &
                                       getDataStorage,                &
                                       setDataStorageLDA,             &
                                       RealType, RealMark,            &
                                       IntegerType, IntegerMark,      &
                                       CharacterType, CharacterMark
!
   use ChemElementModule, only : MaxLenOfAtomName, getZtot, getName
!
   use StringModule, only : initString, setString, endString,         &
                            getNumTokens, readToken
!
   implicit none
!
   character (len=*), intent(in) :: fname
!
   character (len=160) :: text
   character (len=MaxLenOfAtomName) :: atom
   character (len=30) :: column1, column2, column3, column6
   character (len=1) :: dummy
   character (len=1), pointer :: AtomNameCharacter(:)
!
   logical :: ElementName, isDirect, ScaleBravais, isAlloy
!
   integer (kind=IntKind), intent(in), optional :: NumSitesIn
   integer (kind=IntKind), intent(out), optional :: NumSitesOut
   integer (kind=IntKind) :: num_sites
!
   integer (kind=IntKind), parameter :: funit = 103
   integer (kind=IntKind) :: ios, status, n, alen, i, j, k, tend
   integer (kind=IntKind), pointer :: AtomicNumber(:)
   integer (kind=IntKind) :: DataForm
   integer (kind=IntKind) :: StartingLine, NB
   integer (kind=IntKind) :: nsn, nst
   integer (kind=IntKind), parameter :: MaxSingleNumbers = 2
   integer (kind=IntKind), parameter :: MaxSingleStrings = 1
   integer (kind=IntKind), parameter :: IOPE = 0
!
   real (kind=RealKind), pointer :: a0
   real (kind=RealKind) :: x, y, z
   real (kind=RealKind), pointer :: AtomPosition(:,:), p_AtomPosition(:)
   real (kind=RealKind), pointer :: BravaisLatticeVec(:,:)
!
   integer (kind=IntKind), parameter :: MaxComponents = 8 ! This maximum should be more than enough.
   character (len=MaxLenOfAtomName) :: component
   type AlloyStruct
      integer (kind=IntKind) :: NumComponents
      integer (kind=IntKind) :: GlobalIndex
      integer (kind=IntKind) :: element(MaxComponents)  
      real (kind=RealKind) :: content(MaxComponents)
      type (AlloyStruct), pointer :: next
   end type AlloyStruct
!
   integer (kind=IntKind) :: NumAlloySubLatts = 0
   type (AlloyStruct), pointer :: p_header, p_curr
!
   integer (kind=IntKind), allocatable :: AlloyElement(:,:)
   real (kind=RealKind), allocatable :: ElementContent(:,:)
!
   integer (kind=IntKind), pointer :: p_AlloySublattIndex(:)
   integer (kind=IntKind), pointer :: p_NumComponents(:)
   integer (kind=IntKind), pointer :: p_AlloyElement(:,:)
   real (kind=RealKind), pointer :: p_AlloyContent(:,:)
!
   interface
      function isNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isNumber
   end interface
!
   interface
      function isRealNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isRealNumber
   end interface
!
   interface
      function isInteger(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isInteger
   end interface
!
   interface
      function getTokenPosition(k,s,n) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer, intent(out), optional :: n
         integer :: p
      end function getTokenPosition
   end interface
!
   interface
      function trim_string(s1) result(s2)
         character (len=*), intent(in) :: s1
         character (len=len(s1)) :: s2
      end function trim_string
   end interface
!
   interface
      subroutine copyString2CharArray(a,b,n)
         use KindParamModule, only : IntKind
         integer (kind=IntKind), intent(in) :: n
         character (len=*), intent(in) :: a
         character (len=1), intent(out) :: b(:)
      end subroutine copyString2CharArray
   end interface
!
   interface
      function nocaseCompare(s1,s2) result(t)
         implicit none
         logical :: t
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
      end function nocaseCompare
   end interface
!
   isAlloy = .false.
!
!  ===================================================================
!  Determine the number of sites in the position data
!  ===================================================================
   num_sites = 0
   if (MyPE == IOPE) then
      open(unit=funit,file=fname,form='formatted',status='old',iostat=ios, &
           action='read')
      if (ios > 0) then
         call ErrorHandler('readPositionData','iostat > 0',ios)
      endif
!
!     ----------------------------------------------------------------
      call initString(80)
!     ----------------------------------------------------------------
      LOOP_do: do
         read(funit,'(a)',iostat=status)text
         if (status < 0) then
            exit LOOP_do
         endif
         text=trim(adjustl(text))
         if (text(1:1) == '#' .or. text(1:1) == '!'                   &
                              .or. len_trim(text) == 0) then
            cycle LOOP_do
         else
!           ----------------------------------------------------------
            call setString(text)
!           ----------------------------------------------------------
            if (getNumTokens() >= 4) then   ! 4 or more columns of data
               num_sites = num_sites + 1
               if (.not.isAlloy) then       ! check if the system is an alloy
                  call readToken(1,column1,alen)
                  if (.not.isNumber(column1)) then
                     atom = column1(1:alen)
                     if (nocaseCompare(atom(1:3),'CPA')) then
                        isAlloy = .true.
                     endif
                  else
                     read(column1,*) i
                     if (i == -1) then
                        isAlloy = .true.
                     else 
                        call readToken(2,column2,alen)
                        if (.not.isNumber(column2)) then
                           atom = column2(1:alen)
                           if (nocaseCompare(atom(1:3),'CPA')) then
                              isAlloy = .true.
                           endif
                        else
                           read(column2,*) i
                           if (i == -1) then
                              isAlloy = .true.
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo LOOP_do
      close(funit)
!     ----------------------------------------------------------------
      call endString()
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call bcastMessage(num_sites,IOPE)
!  -------------------------------------------------------------------
!
   if (present(NumSitesIn)) then
      if (num_sites /= NumSitesIn) then
         call WarningHandler('readPositionData','num_sites <> NumSitesIn', &
                             num_sites,NumSitesIn)
      endif
   endif
!
   if (num_sites < 1) then
      call ErrorHandler('readPositionData','num_sites < 1',num_sites)
   else if (present(NumSitesOut)) then
      NumSitesOut = num_sites
   endif
!
!  -------------------------------------------------------------------
   call createDataStorage('Position Scaling Factor',RealType)
   call createDataStorage('Bravais Vector',9,RealType)
   call createDataStorage('Atomic Number',num_sites,IntegerType)
   call createDataStorage('Atomic Name',MaxLenOfAtomName*num_sites,CharacterType)
   call createDataStorage('Atomic Position',3*num_sites,RealType)
!  -------------------------------------------------------------------
   a0 => getDataStorage('Position Scaling Factor',RealMark)
   BravaisLatticeVec => getDataStorage('Bravais Vector',3,3,RealMark)
   AtomicNumber => getDataStorage('Atomic Number',num_sites,IntegerMark)
   AtomNameCharacter => getDataStorage('Atomic Name',MaxLenOfAtomName*num_sites,CharacterMark)
   AtomPosition => getDataStorage('Atomic Position',3,num_sites,RealMark)
   p_AtomPosition => getDataStorage('Atomic Position',3*num_sites,RealMark)
!  -------------------------------------------------------------------
!
   NumAlloySubLatts = 0
!
   if (MyPE == IOPE) then
      write(6,'(/,a,a)')' The Position Data File: ',trim(fname)
!
      open(unit=funit,file=fname,form='formatted',status='old',iostat=ios, &
           action='read')
      if (ios > 0) then
         call ErrorHandler('readPositionData','iostat > 0',ios)
      endif
!
!     ================================================================
!     Determine the format of the position data
!     ================================================================
      DataForm = -1
      isDirect = .false.
      ElementName = .false.
      ScaleBravais = .false.
      StartingLine = 0
      NB = 0
      nsn = 0
      nst = 0
      a0 = 1.0d0
!     ----------------------------------------------------------------
      call initString(80)
!     ----------------------------------------------------------------
      do
         read(funit,'(a)',iostat=status)text
         if (status < 0) then
            call ErrorHandler('readPositionData','Invalid file',fname)
         else
            StartingLine = StartingLine + 1
         endif
         text=trim(adjustl(text))
         if (text(1:1) == '#' .or. text(1:1) == '!'                   &
                              .or. len_trim(text) == 0) then
            cycle
         else
!           ----------------------------------------------------------
            call setString(text)
!           ----------------------------------------------------------
            n = getNumTokens()
!           ----------------------------------------------------------
            if (n == 1) then   ! the first line may only have 1 number
               if (isNumber(text)) then
                  nsn = nsn + 1
                  if (nsn > MaxSingleNumbers) then
                     call ErrorHandler('readPositionData','Invalid line #1',text)
                  endif
                  if (isRealNumber(text)) then
                     read(text,*,iostat=status)a0
                     if (status > 0) then
                        call ErrorHandler('readPositionData','Invalid line #2', &
                                          text)
                     else if (a0 <= 1.0d-6) then
                        call ErrorHandler('readPositionData','Invalid line #3', &
                                          text)
                     endif
                     ScaleBravais = .true.
                  endif
               else if (nst < MaxSingleStrings) then
                  nst = nst + 1
                  if (text(1:1) == 'D' .or. text(1:1) == 'd') then
                     isDirect = .true.
                  else if (text(1:1) == 'C' .or. text(1:1) == 'c') then
                     isDirect = .false.
                  else
                     call ErrorHandler('readPositionData','Invalid line #4',text)
                  endif
!                 if (.not.Jump1stLine) then
!                    Jump1stLine = .true.
               else
                  call ErrorHandler('readPositionData','Invalid line #5',text)
               endif
            else if (n == 3) then   ! 3-column data
               NB = NB + 1
               if (NB > 3) then
                  call ErrorHandler('readPositionData','Invalid line #6',text)
               else
                  read(text,*,iostat=status)BravaisLatticeVec(1:3,NB)
                  if (status > 0) then
                     call ErrorHandler('readPositionData','Invalid line #7',text)
                  endif
               endif
            else 
!              =======================================================
!              determine DataForm:
!                   = 0,  atomic_symbol  pos_x   pos_y   pos_z  ...
!                         atomic_number  pos_x   pos_y   pos_z  ...
!                   = 1,  index   atomic_symbol  pos_x   pos_y   pos_z  ...
!                         index   atomic_number  pos_x   pos_y   pos_z  ...
!                   = 2,  atomc_symbol node_number  pos_x   pos_y   pos_z  lmax  ...
!                         atomc_number node_number  pos_x   pos_y   pos_z  lmax  ...
!              determine ElementName:
!                   = true,  atomic_symbol is used
!                   = false, atomic_number is used
!              =======================================================
               if (n >= 4) then   ! 4 or more columns of data
!                 ----------------------------------------------------
                  call readToken(1,column1,alen)
                  call readToken(2,column2,alen)
                  call readToken(3,column3,alen)
!                 ----------------------------------------------------
                  if (n >= 6) then   ! 6 or more columns of data
!                    -------------------------------------------------
                     call readToken(6,column6,alen)
!                    -------------------------------------------------
                  else
                     column6 = 'xxxxx'
                  endif
               else
!                 ----------------------------------------------------
                  call ErrorHandler('readPositionData','Invalid data form #1',text)
!                 ----------------------------------------------------
               endif
               if ( .not.isNumber(column3) ) then
!                 ----------------------------------------------------
                  call ErrorHandler('readPositionData','Invalid data form #2',text,column3)
!                 ----------------------------------------------------
               else if ( isInteger(column2) .and. isInteger(column6) ) then
                  if ( isNumber(column1) ) then
                     ElementName = .false.
                  else
                     ElementName = .true.
                  endif
                  DataForm = 2
               else if ( .not.isNumber(column1) ) then
                  ElementName = .true.
                  DataForm = 0
               else if ( .not.isInteger(column1) ) then
!                 ----------------------------------------------------
                  call ErrorHandler('readPositionData','Invalid data form #3',text,column1)
!                 ----------------------------------------------------
               else if ( .not.isNumber(column2) ) then
                  ElementName = .true.
                  DataForm = 1
               else if (isInteger(column2) ) then
                  ElementName = .false.
                  DataForm = 1
               else
                  ElementName = .false.
                  DataForm = 0
               endif
               StartingLine = StartingLine - 1
               exit
            endif
         endif
      enddo
!
      if (NB < 3 .or. StartingLine < 3) then
!        -------------------------------------------------------------
         call ErrorHandler('readPositionData','Invalid position data',fname)
!        -------------------------------------------------------------
      else if (DataForm /= 0 .and. DataForm /= 1) then
!        -------------------------------------------------------------
         call ErrorHandler('readPositionData','Unknown data format',DataForm)
!        -------------------------------------------------------------
      endif
!
      rewind(funit)
!
!     ================================================================
!     StartingLine+1 is the line number of the first position data.
!     ================================================================
      do n = 1, StartingLine
         read(funit,'(a)')dummy
      enddo
!
      nullify(p_header, p_curr)
!
      n=0
      do
         read(funit,'(a)',iostat=status)text
         if (status < 0 .or. n == num_sites) then
            exit
         endif
         text=adjustl(text)
         if (text(1:1) == '#' .or. text(1:1) == '!'                   &
                              .or. len_trim(text) == 0) then
            cycle
         else
            n = n+1
            if (DataForm == 0) then
               if (ElementName) then
                  j = getTokenPosition(1,text,k)
!                 k = getTokenPosition(2,text)
!                 atom = trim_string(text(j:k-1))
                  atom = text(j:j+k-1)
!                 AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName) = atom(1:MaxLenOfAtomName)
!                 call copyString2CharArray(atom,AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName), &
!                                           MaxLenOfAtomName)
                  AtomicNumber(n) = getZtot(atom(1:MaxLenOfAtomName))
                  k = getTokenPosition(2,text)
                  read(text(k:),*) AtomPosition(1:3,n)
               else
                  read(text,*) AtomicNumber(n), AtomPosition(1:3,n)
                  if (AtomicNumber(n) == -1) then
                     atom = 'CPA'
                  else
                     atom = getName(AtomicNumber(n))
!                    AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName) = atom(1:MaxLenOfAtomName)
!                    call copyString2CharArray(atom,AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName), &
!                                              MaxLenOfAtomName)
                  endif
               endif
               if (atom(1:3) == 'CPA') then
                  tend = getTokenPosition(5,text)  ! The component information starts from the 5th column
               endif
            else if (DataForm == 1) then
               if (ElementName) then
                  j = getTokenPosition(2,text,k)
!                 k = getTokenPosition(3,text)
!                 atom = trim_string(text(j:k-1))
                  atom = text(j:j+k-1)
!                 AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName) = atom(1:MaxLenOfAtomName)
!                 call copyString2CharArray(atom,AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName), &
!                                           MaxLenOfAtomName)
                  AtomicNumber(n) = getZtot(atom(1:MaxLenOfAtomName))
                  k = getTokenPosition(3,text)
                  read(text(k:),*) AtomPosition(1:3,n) 
               else
                  read(text,*) k, AtomicNumber(n), AtomPosition(1:3,n)
                  if (AtomicNumber(n) == -1) then
                     atom = 'CPA'
                  else
                     atom = getName(AtomicNumber(n))
!                    AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName) = atom(1:MaxLenOfAtomName)
!                    call copyString2CharArray(atom,AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName), &
!                                              MaxLenOfAtomName)
                  endif
               endif
               if (atom(1:3) == 'CPA') then
                  tend = getTokenPosition(6,text)  ! The component information starts from the 6th column
               endif
            else if (DataForm == 2) then
               if (ElementName) then
                  j = getTokenPosition(1,text,k)
!                 k = getTokenPosition(2,text)
!                 atom = trim_string(text(j:k-1))
                  atom = text(j:j+k-1)
!                 AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName) = atom(1:MaxLenOfAtomName)
!                 call copyString2CharArray(atom,AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName), &
!                                           MaxLenOfAtomName)
                  AtomicNumber(n) = getZtot(atom(1:MaxLenOfAtomName))
                  k = getTokenPosition(3,text)
                  read(text(k:),*) AtomPosition(1:3,n)
               else
                  read(text,*) AtomicNumber(n),k, AtomPosition(1:3,n) 
                  if (AtomicNumber(n) == -1) then
                     atom = 'CPA'
                  else
                     atom = getName(AtomicNumber(n))
!                    AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName) = atom(1:MaxLenOfAtomName)
!                    call copyString2CharArray(atom,AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName), &
!                                              MaxLenOfAtomName)
                  endif
               endif
               if (atom(1:3) == 'CPA') then
                  tend = getTokenPosition(13,text) ! The component information starts from the 13th column
               endif
            endif
!           ----------------------------------------------------------
            call copyString2CharArray(atom,AtomNameCharacter((n-1)*MaxLenOfAtomName+1:n*MaxLenOfAtomName), &
                                      MaxLenOfAtomName)
!           ----------------------------------------------------------
            if (isAlloy) then
               if ( NumAlloySubLatts == 0 ) then
                  allocate( p_header )
                  p_curr => p_header
               else
                  allocate( p_curr%next )
                  p_curr => p_curr%next
               endif
               if (atom(1:3) == 'CPA') then
!                 ----------------------------------------------------
                  call setString(text(tend:))
!                 ----------------------------------------------------
                  k = getNumTokens()
!                 ----------------------------------------------------
                  if (mod(k,2) /= 0) then
!                 ----------------------------------------------------
                     call ErrorHandler('readPositionData','Invalid component data',text(tend:))
!                 ----------------------------------------------------
                  endif
                  p_curr%NumComponents = k/2
                  do i = 1, p_curr%NumComponents
                     k = tend + getTokenPosition(2*i-1,text(tend:)) - 1
                     read(text(k:),*) component,p_curr%content(i)
                     if (isInteger(component)) then
                        read(component,*) p_curr%element(i)
                     else
                        p_curr%element(i) = getZtot(component)
                     endif
                  enddo
               else
                  p_curr%NumComponents = 1
                  p_curr%content(1) = ONE
                  p_curr%element(1) = AtomicNumber(n)
               endif
               if (p_curr%NumComponents > MaxComponents) then
!                  ---------------------------------------------------
                   call ErrorHandler('readPositionData','Too many alloy components',text(tend:))
!                  ---------------------------------------------------
               endif
               NumAlloySubLatts = NumAlloySubLatts + 1
               p_curr%GlobalIndex = n
            endif
         endif
      enddo
!
      if (n < num_sites) then
!        -------------------------------------------------------------
         call ErrorHandler('readPositionData','Position data is incomplete')
!        -------------------------------------------------------------
      endif
!
!     ----------------------------------------------------------------
      call endString()
!     ----------------------------------------------------------------
!
      close(funit)
      write(6,'(/,a,/)')' The Position Data is read successfully!'
!
      if (ScaleBravais) then
         BravaisLatticeVec(1:3,1:3) = a0*BravaisLatticeVec(1:3,1:3)
      endif
!
      if (isDirect) then
         do n = 1, num_sites
            x = AtomPosition(1,n)
            y = AtomPosition(2,n)
            z = AtomPosition(3,n)
            AtomPosition(1,n) = x*BravaisLatticeVec(1,1) +            &
                                y*BravaisLatticeVec(1,2) +            &
                                z*BravaisLatticeVec(1,3)
            AtomPosition(2,n) = x*BravaisLatticeVec(2,1) +            &
                                y*BravaisLatticeVec(2,2) +            &
                                z*BravaisLatticeVec(2,3)
            AtomPosition(3,n) = x*BravaisLatticeVec(3,1) +            &
                                y*BravaisLatticeVec(3,2) +            &
                                z*BravaisLatticeVec(3,3)
         enddo
      else if (ScaleBravais) then
         do n = 1, num_sites
            AtomPosition(1,n) = a0*AtomPosition(1,n)
            AtomPosition(2,n) = a0*AtomPosition(2,n)
            AtomPosition(3,n) = a0*AtomPosition(3,n)
         enddo
      endif
   endif
!
!  -------------------------------------------------------------------
   call bcastMessage(BravaisLatticeVec,3,3,IOPE)
   call bcastMessage(AtomNameCharacter,MaxLenOfAtomName*num_sites,IOPE)
   call bcastMessage(AtomicNumber,num_sites,IOPE)
   call bcastMessage(p_AtomPosition,3*num_sites,IOPE)
   call bcastMessage(NumAlloySubLatts,IOPE)
!  -------------------------------------------------------------------
!
   if (NumAlloySubLatts > 0) then
      allocate( AlloyElement(MaxComponents+2,NumAlloySubLatts) )
      allocate( ElementContent(MaxComponents,NumAlloySubLatts) )
!
      if (MyPE == IOPE) then
         do i = 1, NumAlloySubLatts
            AlloyElement(1,i) = p_header%GlobalIndex
            AlloyElement(2,i) = p_header%NumComponents
            do j = 1, p_header%NumComponents
               AlloyElement(j+2,i) = p_header%element(j)
               ElementContent(j,i) = p_header%content(j)
            enddo
            p_curr => p_header%next
            nullify( p_header%next )
            deallocate( p_header )
            p_header => p_curr
         enddo
      endif
!
!     ----------------------------------------------------------------
      call bcastMessage(AlloyElement,MaxComponents+2,NumAlloySubLatts,IOPE)
      call bcastMessage(ElementContent,MaxComponents,NumAlloySubLatts,IOPE)
!     ----------------------------------------------------------------
      call createDataStorage('Alloy Sublattice Index',NumAlloySubLatts,IntegerType)
      call createDataStorage('Number of Components on Sublattice',NumAlloySubLatts,IntegerType)
      call createDataStorage('Alloy Element',MaxComponents*NumAlloySubLatts,IntegerType)
      call setDataStorageLDA('Alloy Element',MaxComponents)
      call createDataStorage('Alloy Content',MaxComponents*NumAlloySubLatts,RealType)
      call setDataStorageLDA('Alloy Content',MaxComponents)
!     ----------------------------------------------------------------
      p_AlloySublattIndex => getDataStorage('Alloy Sublattice Index',NumAlloySubLatts,IntegerMark)
      p_NumComponents => getDataStorage('Number of Components on Sublattice',NumAlloySubLatts,IntegerMark)
      p_AlloyElement => getDataStorage('Alloy Element',MaxComponents,NumAlloySubLatts,IntegerMark)
      p_AlloyContent => getDataStorage('Alloy Content',MaxComponents,NumAlloySubLatts,RealMark)
!     ----------------------------------------------------------------
      do i = 1, NumAlloySubLatts
         p_AlloySublattIndex(i) = AlloyElement(1,i)
         p_NumComponents(i) = AlloyElement(2,i)
         do j = 1, AlloyElement(2,i)
            p_AlloyElement(j,i) = AlloyElement(j+2,i)
            p_AlloyContent(j,i) = ElementContent(j,i)
         enddo
      enddo
!
      deallocate( AlloyElement, ElementContent )
      nullify( p_AlloySublattIndex, p_NumComponents, p_AlloyElement, p_AlloyContent )
   endif
!
   nullify( a0, BravaisLatticeVec, AtomNameCharacter, AtomicNumber, AtomPosition )
!
   end subroutine readPositionData
