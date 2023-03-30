!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   program generateAtomPosition
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, HALF, ONE, TWO, THREE, FOUR, TEN2m6, TEN2m10
!
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
!
   use SortModule, only : HeapSort
!
   use ChemElementModule, only : getZtot, getName
!
   use StringModule, only : initString, endString, setString, getNumTokens, readToken
!
   use MatrixInverseModule, only : MtxInv_GE
!
   use SampleModule, only : initSample, getUnitCell, placeAtoms,       &
                            getAtomPosition, getAtomName, getNumAtoms, &
                            endSample, MaxShells, MaxAtomTypes
!
   implicit none
!
   logical :: RandomMedium = .true.
   logical :: RandomCluster = .true.
!
!   integer (kind=IntKind), parameter :: MaxAtomTypes = 20
   integer (kind=IntKind), parameter :: MaxBounds = 50
   integer (kind=IntKind), parameter :: MaxBasis = 250
   integer (kind=IntKind), parameter :: MaxClusters = 10
   integer (kind=IntKind), parameter :: MaxECA = 20 ! max number of empty cell atoms in small cell
   integer (kind=IntKind), parameter :: VaZ = 200   ! Shift the vacancy Z by VaZ
!
   character (len=2) :: Cluster(MaxAtomTypes)
   character (len=2) :: Medium(MaxAtomTypes)
   character (len=2) :: MediumBasis(MaxBasis)
   character (len=2) :: ClusterBasis(MaxBasis)
   character (len=60) :: text, file_name, anm
   character (len=80) :: string80
!
   character (len=2), pointer :: AtomName_medium(:)
   character (len=2), pointer :: AtomName_cluster(:)
!
   integer (kind=IntKind) :: alen, ClusterShape, lattice, NumBounds, iconv
   integer (kind=IntKind) :: NumBasis(3), nbasis, na, nb, nc, NumAtoms
   integer (kind=IntKind) :: i, j, k, ib, n, ncl, status
   integer (kind=IntKind) :: NumMediumAtomTypes, NumClusterAtomTypes
   integer (kind=IntKind) :: NumClusters
   integer (kind=IntKind) :: NumMediumAtoms, NumClusterAtoms(MaxClusters)
   integer (kind=IntKind) :: ordered, embed, rloop
   integer (kind=IntKind) :: ftype
   integer (kind=IntKind) :: nshell
   integer (kind=IntKind) :: add_eca, num_eca, eca_index(MaxECA), basis2eca(MaxBasis)
!
   integer (kind=IntKind) :: NumMediumAtomsOfType(MaxAtomTypes)
   integer (kind=IntKind) :: NumClusterAtomsOfType(MaxAtomTypes)
   integer (kind=IntKind), allocatable :: ClusterFlag(:)
   integer (kind=IntKind), allocatable :: b_cluster(:)
   integer (kind=IntKind), allocatable :: AtomZ(:)
   integer (kind=IntKind), allocatable :: IndexN(:)
   integer (kind=IntKind), allocatable :: AtomZ_Cluster(:)
   integer (kind=IntKind), allocatable :: IndexN_Cluster(:)
!
   integer (kind=IntKind) :: str_type
! GDS
   integer (kind=IntKind) :: ipiv(3), lwork, info
!
   real (kind=RealKind) :: ClusterContent(MaxAtomTypes)
   real (kind=RealKind) :: MediumContent(MaxAtomTypes)
   real (kind=RealKind) :: weight(MaxShells)
   real (kind=RealKind) :: bins(MaxAtomTypes+1)
   real (kind=RealKind) :: BoundVec(3,MaxBounds), BoundV2(MaxBounds)
!
   real (kind=RealKind), pointer :: box(:,:)
   real (kind=RealKind) :: small_box(3,3), rpos(3), box_inv(3,3)
   real (kind=RealKind) :: BasisVec(3,4,3), bv(3,MaxBasis)
   real (kind=RealKind) :: srop(MaxAtomTypes*(MaxAtomTypes-1)/2,MaxShells), Tmax, Tstep
!
   real (kind=RealKind), pointer :: x_medium(:), y_medium(:), z_medium(:)
   real (kind=RealKind), allocatable :: x_cluster(:), y_cluster(:), z_cluster(:)
!
   real (kind=RealKind), parameter :: anstr2au = 1.88973000d0
   real (kind=RealKind), parameter :: au2anstr = 0.52917000d0
   real (kind=RealKind) :: uconv, fact, a0, b0, c0, x, y, z
! GDS   
   real (kind=RealKind) :: work(3)
!
   logical :: isDirect
!
!  data a0/3.8525, 3.8525, 3.7133/
!  data cut/19.3/
!
   data NumBasis(1:3)/4, 2, 1/
   data BasisVec(1:3,1:4,1) &
      & /ZERO, ZERO, ZERO, HALF, HALF, ZERO, HALF, ZERO, HALF, ZERO, HALF, HALF/
   data BasisVec(1:3,1:4,2) &
      & /ZERO, ZERO, ZERO, HALF, HALF, HALF, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
   data BasisVec(1:3,1:4,3) &
      & /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
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
   write(6,'(/)')
   write(6,'(10x,a)')'*******************************************************'
   write(6,'(10x,a)')'*                                                     *'
   write(6,'(10x,a)')'*    Utility Program for Generating Atomic Position   *'
   write(6,'(10x,a)')'*                                                     *'
   write(6,'(10x,a)')'*                    version 2.0                      *'
   write(6,'(10x,a)')'*                                                     *'
   write(6,'(10x,a)')'*******************************************************'
!
   MediumBasis(:)  = '   '
   ClusterBasis(:) = '   '
!
!  ==========================================================================
!  read input parameters for setting up the big box and the underlying lattice
!  ==========================================================================
   write(6,'(//,a)')  &
      'Please Enter the Underlying Lattice Information as Follows ---->'
   lattice = -1
   do while (lattice < 0 .or. lattice > 4)
      write(6,'(/,2x,a,$)')     &
  'Choose (1. Face Centered;  2. Body Centered;  3. Orthorhombic;  4. POSCAR;  0. MuST position file): '
      read(5,*)lattice
   enddo
   write(6,'(i3)')lattice
!
   if (lattice == 1) then
      str_type = -1
      do while (str_type < 0 .or. str_type > 2)
         write(6,'(/,2x,a,$)')'Choose crystal structure type (1. NaCl;  2. Diamond or ZnS;  0. Other): '
         read(5,*)str_type
      enddo
      write(6,'(a,i3)')'Structure type of the FCC lattice: ',str_type
   else
      str_type = 0
   endif
!
   iconv = -1
   write(6,'(//,a)')  &
      'Please Enter the Units for Input Lattice Constants ---->'
   do while (iconv /= 0 .and. iconv /= 1)
      write(6,'(/,2x,a,$)') 'Choose (0. Atomic Units;  1. Angstrom): '
      read(5,*)iconv
   enddo
   write(6,'(i3)')iconv
   if (iconv == 0) then
       uconv = ONE
   else
       uconv = anstr2au
   endif
!
   num_eca = 0
   eca_index = 0
   basis2eca = 0
!
   if (lattice == 0) then
      write(6,'(/,2x,a,$)') 'Name of the Underlying Lattice File: '
      read(5,'(a)')file_name
      write(6,'(a)')trim(adjustl(file_name))
      file_name = adjustl(file_name)
!     ================================================================
      open(unit=11,file=file_name,form='formatted',status='old')
!     ================================================================
!     read line-by-line of the position file and process the data
!     ----------------------------------------------------------------
      call initString(80)
!     ----------------------------------------------------------------
      isDirect = .false.  ! Using Cartisian coordinates is the default
      a0 = ONE
      bv = ZERO
      ib = 0
      nbasis = 0
      LOOP_do: do  
         read(11,'(a)',iostat=status) string80
         if (status < 0) then
            exit LOOP_do
         endif
         string80 = trim(adjustl(string80))
         if (string80(1:1) == '#' .or. string80(1:1) == '!' .or. len_trim(string80) == 0) then
            cycle LOOP_do
         endif
!        -------------------------------------------------------------
         call setString(string80)
!        -------------------------------------------------------------
         n = getNumTokens()
         if (n == 1) then
            if (isNumber(string80)) then
               read(string80,*,iostat=status)a0
               if (status > 0) then
                  call ErrorHandler('main','Invalid line',trim(string80))
               endif
            else
               if (string80(1:1) == 'D' .or. string80(1:1) == 'd') then
                  isDirect = .true.
               else if (string80(1:1) == 'C' .or. string80(1:1) == 'c') then
                  isDirect = .false.
               else
                  call ErrorHandler('main','Invalid line',trim(string80))
               endif
            endif
         else if (n == 3) then
            ib = ib + 1
            if (ib > 3) then
               call ErrorHandler('main','The number of Bravais lattice vectors > 3')
            endif
            read(string80,*,iostat=status)small_box(1:3,ib)
            if (status > 0) then
               call ErrorHandler('main','Invalid line',trim(string80))
            endif
         else if (n > 3) then
            nbasis = nbasis + 1
            if (nbasis > MaxBasis) then
               call ErrorHandler('main','nbasis > MaxBasis',nbasis,MaxBasis)
            endif
            k = getTokenPosition(2,string80)
            read(string80(1:k-1),'(a)')anm
            if (isNumber(anm)) then
               read(anm,*)n
               MediumBasis(nbasis) = getName(n)
            else if (len_trim(anm) > 2) then
               call ErrorHandler('main','Invalid element name',trim(anm))
            else
               MediumBasis(nbasis) = trim(anm)
            endif
            read(string80(k:),*)bv(1:3,nbasis)
         else
            call ErrorHandler('main','Invalid line',trim(string80))
         endif
      enddo LOOP_do
!     ================================================================
      close(11)
!     ================================================================
      if (ib /= 3) then
         call ErrorHandler('main','The number of Bravais lattice vectors /= 3')
      endif
      if (isDirect) then
         do i = 1, nbasis
            x = bv(1,i)
            y = bv(2,i)
            z = bv(3,i)
            bv(1,i) = x*small_box(1,1)+y*small_box(1,2)+z*small_box(1,3)
            bv(2,i) = x*small_box(2,1)+y*small_box(2,2)+z*small_box(2,3)
            bv(3,i) = x*small_box(3,1)+y*small_box(3,2)+z*small_box(3,3)
         enddo
      endif
!!    small_box = small_box*a0
!!    bv = bv*a0
   else if (lattice == 4 ) then
!     write(6,'(/,2x,a,$)') 'Name of the Underlying Lattice File (e.g., POSCAR): '
!     read(5,'(a)')file_name
!     file_name = adjustl(file_name)
      file_name = 'POSCAR'
      write(6,'(a)')trim(file_name)
      open(unit=11,file=file_name,form='formatted',status='old')
      read(11,'(a)')text
      read(11,*)fact
      read(11,*)small_box(1:3,1)
      read(11,*)small_box(1:3,2)
      read(11,*)small_box(1:3,3)
!!    uconv = uconv*fact
      a0 = fact*uconv
      small_box(1:3,1:3) = small_box(1:3,1:3)*uconv
      read(11,'(a)')text
      call initString(text)
      text = adjustl(text)
      na = getNumTokens()
      NumMediumAtomTypes = na
!     GDS added this line
      read(11,'(a)')text
      read(text,*) (NumMediumAtomsOfType(i),i=1,na)
!      print *, NumMediumAtomsOfType(1)
      write(6,'(/,2x,a,$)') 'The atom names (separated by spaces and/or comma) in the order defined in POSCAR (e.g., Cu Al Fe): '
      read(5,'(a)')text
!      print *, text
      text = adjustl(text)
      write(6,'(a)')trim(text)
      call setString(text)
      na = getNumTokens()
      if (na /= NumMediumAtomTypes) then
         call ErrorHandler('main','na <> NumMediumAtomTypes',na,NumMediumAtomTypes)
      endif
      j = 1
      do i = 1,na
         call readToken(i,anm,n)
         if (n >= 3) then
            Medium(i) = anm(1:2)
         else
            Medium(i) = anm
         endif
         MediumBasis(j:j+NumMediumAtomsOfType(i)-1) = Medium(i)
         j = j + NumMediumAtomsOfType(i)
      enddo
!     text = adjustl(text)
      call endString()
      nbasis = 0
      do i = 1,na
         nbasis = nbasis + NumMediumAtomsOfType(i)
      enddo
      if (nbasis < 1 .or. nbasis > MaxBasis) then
         call ErrorHandler('main','Number of basis is out of range',nbasis)
      endif
      read(11,'(a)')text
      text = adjustl(trim(text))
      if (text(1:1) == 'c' .or. text(1:1) == 'C' .or. text(1:1) == 'k' .or.  text(1:1) == 'K') then
         isDirect = .false.
      else
         isDirect = .true.
      endif
      do i = 1,nbasis
         read(11,'(a)')text
         call initString(text)
         text = adjustl(text)
         na = getNumTokens()
         if (na == 3) then
            read(text,*) rpos(1:3)
         else
            call ErrorHandler('main','Invalid input data',trim(text))
         endif
         call endString()  
         if (isDirect) then
            do j = 1,3
               bv(j,i) = small_box(j,1)*rpos(1)+small_box(j,2)*rpos(2)+small_box(j,3)*rpos(3)
            enddo
         else
            do j = 1,3
               bv(j,i) = rpos(j)*uconv
            enddo
         endif
      enddo
      close(11)
   else
      small_box(1:3,1:3) = ZERO
      if (lattice == 3) then
         if (iconv == 0) then
            write(6,'(/,2x,a,$)')     &
            'Lattice Constants a0, b0, c0 (in a.u. seperated by space or comma): '
         else
            write(6,'(/,2x,a,$)')     &
            'Lattice Constants a0, b0, c0 (in Angstrom seperated by space or comma): '
         endif
         read(5,*)a0, b0, c0
         small_box(1,1) = ONE
         small_box(2,2) = b0/a0
         small_box(3,3) = c0/a0
      else
         if (iconv == 0) then
            write(6,'(/,2x,a,$)')'Lattice Constants a0 (in a.u.): '
         else
            write(6,'(/,2x,a,$)')'Lattice Constants a0 (in Angstrom): '
         endif
         read(5,*)a0
         small_box(1,1) = ONE
         small_box(2,2) = ONE
         small_box(3,3) = ONE
      endif
      a0 = a0*uconv
      small_box(1,1) = small_box(1,1)*uconv
      small_box(2,2) = small_box(2,2)*uconv
      small_box(3,3) = small_box(3,3)*uconv
      if (str_type == 1) then ! NaCl
         nbasis = 0
         do i = 1, 4
            nbasis = nbasis + 1
            bv(1,nbasis) = BasisVec(1,i,1)
            bv(2,nbasis) = BasisVec(2,i,1)
            bv(3,nbasis) = BasisVec(3,i,1)
            nbasis = nbasis + 1
            bv(1,nbasis) = BasisVec(1,i,1) + HALF
            bv(2,nbasis) = BasisVec(2,i,1)
            bv(3,nbasis) = BasisVec(3,i,1)
         enddo
      else if (str_type == 2) then ! diamond or ZnS structure
         add_eca = -1
         do while (add_eca /= 0 .and. add_eca /= 1)
            write(6,'(/,2x,a,$)') 'Add empty cell atoms (0. No;  1. Yes): '
            read(5,*)add_eca
         enddo
         nbasis = 0
         do i = 1, 4
            nbasis = nbasis + 1
            bv(1,nbasis) = BasisVec(1,i,1)
            bv(2,nbasis) = BasisVec(2,i,1)
            bv(3,nbasis) = BasisVec(3,i,1)
            nbasis = nbasis + 1
            bv(1,nbasis) = BasisVec(1,i,1) + ONE/FOUR
            bv(2,nbasis) = BasisVec(2,i,1) + ONE/FOUR
            bv(3,nbasis) = BasisVec(3,i,1) + ONE/FOUR
            if (add_eca == 1) then
               nbasis = nbasis + 1
               bv(1,nbasis) = BasisVec(1,i,1) + HALF
               bv(2,nbasis) = BasisVec(2,i,1) + HALF
               bv(3,nbasis) = BasisVec(3,i,1) + HALF
               do j = 1, 3
                  if (bv(j,nbasis) > ONE) then
                     bv(j,nbasis) = bv(j,nbasis) - ONE
                  endif
               enddo
               num_eca = num_eca + 1
               if (num_eca > MaxECA) then
                  call ErrorHandler('main','num_eca > MaxECA',num_eca,MaxECA)
               endif
               eca_index(num_eca) = nbasis
               basis2eca(nbasis) = num_eca
               nbasis = nbasis + 1
               bv(1,nbasis) = BasisVec(1,i,1) + THREE/FOUR
               bv(2,nbasis) = BasisVec(2,i,1) + THREE/FOUR
               bv(3,nbasis) = BasisVec(3,i,1) + THREE/FOUR
               do j = 1, 3
                  if (bv(j,nbasis) > ONE) then
                     bv(j,nbasis) = bv(j,nbasis) - ONE
                  endif
               enddo
               num_eca = num_eca + 1
               if (num_eca > MaxECA) then
                  call ErrorHandler('main','num_eca > MaxECA',num_eca,MaxECA)
               endif
               eca_index(num_eca) = nbasis
               basis2eca(nbasis) = num_eca
            endif
         enddo
      else
         nbasis = NumBasis(lattice)
         do i=1, nbasis
            bv(1,i) = BasisVec(1,i,lattice)
            bv(2,i) = BasisVec(2,i,lattice)
            bv(3,i) = BasisVec(3,i,lattice)
         enddo
      endif
      do i = 1, nbasis
         bv(1,i) = bv(1,i)*small_box(1,1)
         bv(2,i) = bv(2,i)*small_box(2,2)
         bv(3,i) = bv(3,i)*small_box(3,3)
      enddo
   endif
!
   write(6,'(/,2x,a,$)')       &
           'Number of Repeats of Small Box Along A, B, C Directions: '
   read(5,*)na,nb,nc
   write(6,'(3i5)')na,nb,nc
!
   write(6,'(/,2x,a,$)')       &
           'Type of output format( 0. MuST; 1. Legacy LSMS; 2. VASP (Cartesian) ): '
   read(5,*) ftype
   write(6,'(i3)') ftype
!
!  --------------------------------------------------------------------------
   call initSample(nbasis,na,nb,nc,small_box,bv,basis2eca)
   box => getUnitCell()
   NumAtoms = getNumAtoms()
   x_medium => getAtomPosition(1)
   y_medium => getAtomPosition(2)
   z_medium => getAtomPosition(3)
!  --------------------------------------------------------------------------
!
!            print *, box(1,1), box(1,2), box(1,3)
!            print *, box(2,1), box(2,2), box(2,3)
!            print *, box(3,1), box(3,2), box(3,3)
!   call MtxInv_GE(3,box,box_inv)
! GDS
      do i=1, 3
         box_inv(1,i) = box(1,i)
         box_inv(2,i) = box(2,i)
         box_inv(3,i) = box(3,i)
      enddo
   call dgetrf(3,3,box_inv,3,ipiv,info)
   call dgetri(3,box_inv,3,ipiv,work,9,info)
!            print *, box_inv(1,1), box_inv(1,2), box_inv(1,3)
!            print *, box_inv(2,1), box_inv(2,2), box_inv(2,3)
!            print *, box_inv(3,1), box_inv(3,2), box_inv(3,3)
!
   if (lattice == 0) then
      RandomMedium = .false.
   else
!     =======================================================================
!     read input parameters for the solid solution occupying the lattice
!     =======================================================================
      ordered = -1; rloop = 0
      do while (ordered < 0 .or. ordered > 2)
         write(6,'(//,a)')        &
                 'Please Enter the Solid Solution Information as Follows ---->'
!
         write(6,'(2x,a)') 'Is the solid solution ordered or random?'
         write(6,'(4x,a,$)') 'Enter 0 for ORDERED; 1 for RANDOM; or 2 for RANDOM with SHORT-RANGE ORDER: '
         read(5,*)ordered
         write(6,'(i3)')ordered
         if (ordered == 0) then
            RandomMedium = .false.
         else if (ordered == 1 .or. ordered == 2) then
            RandomMedium = .true.
         else if (rloop <= 5) then
            write(6,'(a,i3)')'Undefined input: ',ordered
            rloop = rloop + 1
         else
            call ErrorHandler('main','Undefined input',ordered)
         endif
      enddo
   endif
!
   if (.not.RandomMedium) then
      if (MediumBasis(1) == '   ') then
         if (str_type == 0) then
            do i = 1, nbasis
               write(6,'(/,2x,a,3f10.5,a,$)')'Enter Atom Name Located at',   &
                                             bv(1:3,i)*uconv,' : '
               read(5,'(a)')MediumBasis(i)
               write(6,'(a)')trim(adjustl(MediumBasis(i)))
            enddo
         else ! This works for NaCl and ZnS or Diadmond structures
            do i = 1, 2
               write(6,'(/,2x,a,3f10.5,a,$)')'Enter Atom Name Located at',   &
                                             bv(1:3,i)*uconv,' : '
               read(5,'(a)')MediumBasis(i)
               write(6,'(a)')trim(adjustl(MediumBasis(i)))
            enddo
            do i = 3, nbasis
               if (basis2eca(i) == 0) then
!                 write(6,'(/,2x,a,3f10.5,a,$)')'Enter Atom Name Located at',   &
!                                             bv(1:3,i)*uconv,' : '
!                 read(5,'(a)')MediumBasis(i)
!                 write(6,'(a)')trim(adjustl(MediumBasis(i)))
                  j = mod(i-1,2) + 1
                  MediumBasis(i) = MediumBasis(j)
               else
                  MediumBasis(i) = 'Va'
               endif
            enddo
         endif
      endif
      k = 1
      NumMediumAtomTypes = 1
      Medium(1) = MediumBasis(1)
      do j = 2,nbasis
         NumMediumAtomTypes = NumMediumAtomTypes + 1
         LoopZ1: do i = j-1,1,-1
            if ( MediumBasis(j) == MediumBasis(i) ) then
               NumMediumAtomTypes = NumMediumAtomTypes -1
               exit LoopZ1
            endif
         enddo LoopZ1
         if ( i==0 ) then
            k = k+1
            Medium(k) = MediumBasis(j)
         endif
      enddo
!     -----------------------------------------------------------------------
      call placeAtoms(nbasis,MediumBasis)
!     -----------------------------------------------------------------------
   else if (RandomMedium) then
      if (str_type == 2 .and. add_eca == 1) then
         write(6,'(/,a)') '=================================================================='
         write(6,'(2x,a)')'A random alloy with empty cells is currently Under Construction!  '
         write(6,'(2x,a)')'                                                                  '
         write(6,'(2x,a)')'A temporary solution is to take the following steps to produce    '
         write(6,'(2x,a)')'the position data file:                                           '
         write(6,'(2x,a)')'  1. Run the random sample generation without empty cells, and    '
         write(6,'(2x,a)')'     rename the position data file (position.dat)                 '
         write(6,'(2x,a)')'  2. Run an ordered sample generation of the same size but        '
         write(6,'(2x,a)')'     with empty cells.                                            '
         write(6,'(2x,a)')'  3. Copy the Va position lines from this newly generated         '
         write(6,'(2x,a)')'     position.dat file and place them in the renamed position     '
         write(6,'(2x,a)')'     data file produced by the random sample generation.          '
         write(6,'(a)')   '=================================================================='
         call StopHandler('main','The requested calculation is Not Implemented!')
      endif
      write(6,'(/,2x,a,$)')     &
       'Constituents (atomic name separated by space or comma, e.g., Fe Ni): '
      read(5,'(a)')text
      write(6,'(a)')trim(adjustl(text))
      call initString(text)
      NumMediumAtomTypes = getNumTokens()
      if (NumMediumAtomTypes < 1 ) then
         call ErrorHandler('main','Number of atom types < 1',NumMediumAtomTypes)
      else if (NumMediumAtomTypes > MaxAtomTypes) then
         call ErrorHandler('main','Number of atom types > limit',NumMediumAtomTypes,MaxAtomTypes)
      endif
      do i = 1, NumMediumAtomTypes
         call readToken(i,Medium(i),alen)
         write(6,'(/,2x,3a,$)')'Content of ',Medium(i),' (>= 0 and =< 1): '
         read(5,*)MediumContent(i)
         write(6,'(f10.5)')MediumContent(i)
      enddo
      call endString()
      write(6,'(a)')' '
!
      bins(1) = ZERO
      do i = 2, NumMediumAtomTypes+1
         bins(i) = bins(i-1) + MediumContent(i-1)
      enddo
      if (abs(bins(NumMediumAtomTypes+1)-ONE) > TEN2m6) then
         call ErrorHandler('main','The summation of the contents is not 1',bins(NumMediumAtomTypes+1))
      endif

!
      nshell=-1
      if (ordered == 2) then   ! Random medium with short range order
         write(6,'(/,2x,a,$)') 'How many shells will be considered: '
         read(5,*)nshell
!        --------------------------------------------------------------------
         if (nshell < 1 ) then
            call ErrorHandler('main','Number of shells < 1',nshell)
         else if (nshell > MaxShells) then
            call ErrorHandler('main','Number of shells > MaxShells',nshell,MaxShells)
         endif

         write(6,'(2x,a,i3)') 'nshell is:', nshell

         write(6,'(/,2x,a,$)') 'For each shell please give the weight of SROs:'
!        --------------------------------------------------------------------
         ib=-1
         read(5,'(a)')text
         call initString(text)
         ib = getNumTokens()
!        -----------------------------------------------------------------
         if (ib /= nshell) then
            call ErrorHandler('main','Number of weight of SRO not correct, should be',nshell)
         endif

         do j = 1, nshell
           call readToken(j,anm,alen)
           read(anm,*)weight(j)
         enddo
         call endString()

         write(6,'(a)')' '
         write(6,'(2x,a,5f15.5)') 'The input weight is:',weight(1:nshell)

!        ------------------------------------------------------------------------
         write(6,'(/,2x,a,$)') 'For each shell please give N*(N-1)/2 SROs:'
         do i = 1, nshell
            ib=-1
            write(6,'(/,2x,a,i2,a,$)') 'shell',i,':  ' 
!        --------------------------------------------------------------------
            read(5,'(a)')text
            call initString(text)
            ib = getNumTokens()
!        --------------------------------------------------------------------
            if (ib /= NumMediumAtomTypes*(NumMediumAtomTypes-1)/2) then
               call ErrorHandler('main','Number of SRO not correct, should be',NumMediumAtomTypes*(NumMediumAtomTypes-1)/2)
            endif

            do j = 1, NumMediumAtomTypes*(NumMediumAtomTypes-1)/2
              call readToken(j,anm,alen)
              read(anm,*)srop(j,i)
            enddo 
            call endString()
         enddo


         write(6,'(/,2x,a)')'the input SROs are: '
         do i=1,nshell
           do j=1, NumMediumAtomTypes*(NumMediumAtomTypes-1)/2
              write(6,'(2x,f15.5,$)') srop(j,i)
           enddo
           write(6,'(a)')' '
         enddo 


         write(6,'(/,2x,a,$)')'Enter the temperature (K) to start annealing: '
         read(5,*)Tmax
         write(6,'(f12.5)')Tmax
         write(6,'(/,2x,a,$)')'Enter the temperature step (K) for cooling: '
         read(5,*)Tstep
         write(6,'(f12.5)')Tstep
         write(6,'(a)')' '
!        --------------------------------------------------------------------
         call placeAtoms(NumMediumAtomTypes,Medium,MediumContent,weight,nshell,srop,Tmax,Tstep)
!        --------------------------------------------------------------------
      else
!        --------------------------------------------------------------------
         call placeAtoms(NumMediumAtomTypes,Medium,MediumContent)
!        --------------------------------------------------------------------
      endif
   endif
!
   AtomName_medium => getAtomName()
!
!  ==========================================================================
!  The following code takes care the situation where there are embedded clusters
!  ==========================================================================
   embed = -1; rloop = 0
   do while (embed < 0 .or. embed > 1) 
      write(6,'(//,a)')'Do you need embed cluster(s) into the solid solution?'
      write(6,'(2x,a,$)')'Enter 0 for NO; or 1 for YES: '
      read(5,*)embed
      write(6,'(i3)')embed
      if (embed == 0) then
         NumClusters = 0
      else if (embed == 1) then
         write(6,'(/,2x,a,$)')'Number of Clusters: '
         read(5,*)NumClusters
         allocate( ClusterFlag(NumAtoms), b_cluster(NumAtoms) )
         allocate( x_cluster(NumAtoms), y_cluster(NumAtoms),            &
                   z_cluster(NumAtoms), AtomName_cluster(NumAtoms),     &
                   AtomZ_Cluster(NumAtoms), IndexN_Cluster(NumAtoms) )
      else if (rloop <= 5) then
         write(6,'(a,i3)')'Invalid input value: ',embed
         rloop = rloop + 1
      else
         call ErrorHandler('main','Invalid input value',embed)
      endif
   enddo
!
   if (NumClusters < 0 .or. NumClusters > MaxClusters) then
      call ErrorHandler('main','Number of clusters is out of range',NumClusters)
   endif
!
!  ===================================================================
!  Output the head lines of the position data
!  ===================================================================
   open(unit=14,file='position.dat',status='unknown',form='formatted')
   fact = ONE
   if (ftype==2) then
      write(14,'(a)') '#POSCAR file'
      write(14,*) a0*au2anstr
      fact = au2anstr
   else
      write(14,'(a)') '# Units:: atomic units'
      write(14,'(f12.8)')a0
      write(14,'(a)') '# If Angsrtroms is needed, comment out the line above and uncomment the following line'
      write(14,'(a,f12.8)')'# ', au2anstr*a0
      write(14,'(a)')      '# =============================================================='
   endif
   write(14,'(2x,3f19.11)')fact*box(1:3,1)
   write(14,'(2x,3f19.11)')fact*box(1:3,2)
   write(14,'(2x,3f19.11)')fact*box(1:3,3)
   if (ftype/=2) then
      write(14,'(a)')      '# =============================================================='
      write(14,'(a,i8)')   '# Number of clusters: ',NumClusters
   endif
!
   NumMediumAtoms = NumAtoms
   do ncl = 1, NumClusters
!     -----------------------------------------------------------------------
      stop 'insertClusters has not been implemented.'
!     call insertClusters(ncl)
!     -----------------------------------------------------------------------
      if (ftype/=2) then
         write(14,'(a,i6)')'#     Number of cluster atoms:  ',NumClusterAtoms(ncl)
         write(14,'(a,i8)')'#     Number of cluster atom type: ',NumClusterAtomTypes
      endif
      do ib = 1, NumClusterAtomTypes
         NumClusterAtomsOfType(ib) = 0
         do i = 1,NumClusterAtoms(ncl)
            if ( Cluster(ib) == AtomName_cluster(i) ) then
               NumClusterAtomsOfType(ib) = NumClusterAtomsOfType(ib)+1
            endif
         enddo
         if (ftype/=2) then
            write(14,'(a,a3,1x,'': '',i8)')'# Number of ',Cluster(ib),    &
                                               NumClusterAtomsOfType(ib)
         endif
      enddo
      if ( ftype==1 ) then
         do i = 1, NumClusterAtoms(ncl)
            k = IndexN_Cluster(i)
            n = AtomZ_Cluster(k)
            write(14,'(i3,1x,3f19.15,1x,f7.4,2x,i5,3x,a)') n,         &
                 x_cluster(k),y_cluster(k),z_cluster(k), 1.0d0, 0, "U"
         enddo
      else if ( ftype==2) then
         do i = 1, NumClusterAtoms(ncl)
            k = IndexN_Cluster(i)
            do j = 1,3
               rpos(j) = box_inv(j,1)*x_cluster(k) +               &
                         box_inv(j,2)*y_cluster(k) +               &
                         box_inv(j,3)*z_cluster(k)
            enddo
            write(14,'(2x,3f19.11)') rpos(1:3)
         enddo
      else
         do i = 1, NumClusterAtoms(ncl)
            k = IndexN_Cluster(i)
            write(14,'(2x,a3,2x,3f19.11)') AtomName_cluster(k),       &
                              x_cluster(k), y_cluster(k),z_cluster(k)
         enddo
      endif
   enddo
!
!  ===================================================================
!  Output the medium atom position data
!  ===================================================================
   if (ftype/=2) then
      write(14,'(a,i8)')'# Number of medium atoms:     ',NumMediumAtoms
      write(14,'(a,i8)')'# Number of medium atom type: ',NumMediumAtomTypes
   endif
   do ib = 1, NumMediumAtomTypes
      NumMediumAtomsOfType(ib) = 0
      do i = 1,NumAtoms
         if ( Medium(ib) == AtomName_medium(i) ) then
            if ( embed==1 ) then
               if ( ClusterFlag(i) == 0 ) then
                   NumMediumAtomsOfType(ib) = NumMediumAtomsOfType(ib)+1
               endif
            else 
               NumMediumAtomsOfType(ib) = NumMediumAtomsOfType(ib)+1
            endif
         endif
      enddo
      if (ftype/=2) then
         write(14,'(a,a3,1x,'': '',i8)')'# Number of ',Medium(ib),       &
                           NumMediumAtomsOfType(ib)
      endif
   enddo
   if ( ftype==2 ) then
      do ib = 1, NumMediumAtomTypes
         write(14,'(a6,$)')Medium(ib)
      enddo
      write(14,'(a)')' '
      do ib = 1, NumMediumAtomTypes
         write(14,'(i5,$)')NumMediumAtomsOfType(ib)
      enddo
      write(14,'(/,a)')'Direct'
   endif
!
   allocate( AtomZ(NumAtoms), IndexN(NumAtoms) )
   do i = 1,NumAtoms
      AtomZ(i) = getZtot(AtomName_medium(i))
      if (AtomZ(i) == 0) then
         AtomZ(i) = VaZ
      endif
   enddo
!  -------------------------------------------------------------------
   call HeapSort(NumAtoms, AtomZ, IndexN)
!  -------------------------------------------------------------------
!
   if (embed == 1) then
      do i = 1, NumMediumAtoms
         k = IndexN(i)
         if (ClusterFlag(i) == 0) then
            if ( ftype==1 ) then
               if (AtomZ(i) == VaZ) then
                  n = 0
               else
                  n = AtomZ(i)
               endif
               write(14,'(i3,1x,3f19.15,1x,f7.4,2x,i5,3x,a)') n,      &
                   x_medium(k),y_medium(k),z_medium(k), 1.0d0, 0, "U"
            else if ( ftype==2) then
               do j = 1,3 
                  rpos(j) = box_inv(j,1)*x_medium(k) +                &
                            box_inv(j,2)*y_medium(k) +                &
                            box_inv(j,3)*z_medium(k)
                  enddo
               write(14,'(2x,a3,2x,3f19.11)')AtomName_medium(k), rpos(1:3)
            else
               write(14,'(2x,a3,2x,3f19.11)') AtomName_medium(k),     &
                                   x_medium(k),y_medium(k),z_medium(k)
            endif
         endif
      enddo
      deallocate( ClusterFlag, b_cluster )
      deallocate( x_cluster, y_cluster, z_cluster, AtomName_cluster )
      deallocate( AtomZ_Cluster, IndexN_Cluster )
   else
      if ( ftype/=2 ) then
         write(14,'(a)')'# =============================================================='
      endif
      do i = 1, NumMediumAtoms
         k = IndexN(i)
         if ( ftype==1 ) then
            if (AtomZ(i) == VaZ) then
               n = 0
            else
               n = AtomZ(i)
            endif
            write(14,'(i3,1x,3f19.15,1x,f7.4,2x,i5,3x,a)') n,      &
                x_medium(k),y_medium(k),z_medium(k), 1.0d0, 0, "U"
         else if ( ftype==2) then
            do j = 1,3
               rpos(j) = box_inv(j,1)*x_medium(k) +                &
                         box_inv(j,2)*y_medium(k) +                &
                         box_inv(j,3)*z_medium(k)
            enddo
            write(14,'(2x,a3,2x,3f19.11)')AtomName_medium(k), rpos(1:3)
            !write(14,'(2x,3f19.11)') rpos(1:3)
!            print *, box_inv(1,1), box_inv(1,2), box_inv(1,3)
!            print *, box_inv(2,1), box_inv(2,2), box_inv(2,3)
!            print *, box_inv(3,1), box_inv(3,2), box_inv(3,3)
!            print *, x_medium(k), y_medium(k), z_medium(k)
!            print *, AtomName_medium(k), rpos(1:3)
         else
            write(14,'(2x,a3,2x,3f19.11)') AtomName_medium(k),     &
                                x_medium(k),y_medium(k),z_medium(k)
         endif
      enddo
   endif
   write(6,*) " "
!   write(14,*) " "
   close(14)
!
   nullify( x_medium, y_medium, z_medium, AtomName_medium )
   deallocate( AtomZ, IndexN )
!
!  -------------------------------------------------------------------
   call endSample()
!  -------------------------------------------------------------------
!
   end program generateAtomPosition
