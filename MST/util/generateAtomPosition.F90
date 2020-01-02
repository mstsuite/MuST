!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   program generateAtomPosition
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, HALF, ONE, TWO, TEN2m6, TEN2m10
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use SortModule, only : HeapSort
!
   use ChemElementModule, only : getZtot
!
   use StringModule, only : initString, endString, setString, getNumTokens, readToken
!
   use MatrixInverseModule, only : MtxInv_GE
!
   use SampleModule, only : initSample, getUnitCell, placeAtoms,       &
                            getAtomPosition, getAtomName, getNumAtoms, &
                            getNumSpecies, substituteAtoms,            &
                            isSampleAtom, endSample, MaxShells, MaxAtomTypes
!
   use DataServiceCenterModule, only : initDataServiceCenter,         &
                                       endDataServiceCenter,          &
                                       getDataStorage,                &
                                       RealType, RealMark,            &
                                       IntegerType, IntegerMark,      &
                                       CharacterType, CharacterMark
!
   implicit none
!
   logical :: RandomHost = .true.
   logical :: RandomCluster = .true.
   logical :: file_exist = .true.
   logical :: done = .true.
!
!   integer (kind=IntKind), parameter :: MaxAtomTypes = 20
   integer (kind=IntKind), parameter :: MaxBounds = 50
   integer (kind=IntKind), parameter :: MaxBasis = 1024
   integer (kind=IntKind), parameter :: MaxClusters = 10
   integer (kind=IntKind), parameter :: NumLattices = 4
!
   character (len=2) :: Impurity(MaxAtomTypes)
   character (len=2) :: ReplacedHost(MaxAtomTypes)
   character (len=2) :: Cluster(MaxAtomTypes)
   character (len=2) :: Host(MaxAtomTypes)
   character (len=2) :: HostBasis(MaxBasis)
   character (len=2) :: ClusterBasis(MaxBasis)
   character (len=60) :: text, file_name, anm
!
   character (len=2), allocatable :: AtomName_cluster(:)
!
   integer (kind=IntKind) :: alen, ClusterShape, lattice, NumBounds, iconv
   integer (kind=IntKind) :: NumBasis(NumLattices), nbasis, na, nb, nc, NumAtoms
   integer (kind=IntKind) :: i, j, k, ib, n, ncl
   integer (kind=IntKind) :: NumHostAtomTypes, NumClusterAtomTypes
   integer (kind=IntKind) :: NumClusters, NumImpuritySpecies, NumAtomSpecies
   integer (kind=IntKind) :: NumHostAtoms, NumClusterAtoms(MaxClusters)
   integer (kind=IntKind) :: ordered, substitution, rloop
   integer (kind=IntKind) :: ftype
   integer (kind=IntKind) :: nshell
   integer (kind=IntKind) :: ios
   integer (kind=IntKind) :: iseed
!
   integer (kind=IntKind) :: NumHostAtomsOfType(MaxAtomTypes)
   integer (kind=IntKind) :: NumClusterAtomsOfType(MaxAtomTypes)
   integer (kind=IntKind), allocatable :: ClusterFlag(:)
   integer (kind=IntKind), allocatable :: b_cluster(:)
   integer (kind=IntKind), allocatable :: AtomZ(:)
   integer (kind=IntKind), allocatable :: IndexN(:)
   integer (kind=IntKind), allocatable :: AtomZ_Cluster(:)
   integer (kind=IntKind), allocatable :: IndexN_Cluster(:)
!
   real (kind=RealKind) :: ImpurityContent(MaxAtomTypes)
   real (kind=RealKind) :: ClusterContent(MaxAtomTypes)
   real (kind=RealKind) :: HostContent(MaxAtomTypes)
   real (kind=RealKind) :: weight(MaxShells)
   real (kind=RealKind) :: bins(MaxAtomTypes+1)
   real (kind=RealKind) :: BoundVec(3,MaxBounds), BoundV2(MaxBounds)
!
   real (kind=RealKind), pointer :: box(:,:)
   real (kind=RealKind) :: small_box(3,3), rpos(3), box_inv(3,3)
   real (kind=RealKind) :: BasisVec(3,16,NumLattices), bv(3,MaxBasis)
   real (kind=RealKind) :: srop(MaxAtomTypes*(MaxAtomTypes-1)/2,MaxShells), Tmax, Tstep
!
   real (kind=RealKind), pointer :: Bravais(:,:), AtomPosition(:,:)
   real (kind=RealKind), pointer :: x_host(:), y_host(:), z_host(:)
   real (kind=RealKind), allocatable :: x_cluster(:), y_cluster(:), z_cluster(:)
!
   real (kind=RealKind), parameter :: anstr2au = 1.8897518760
   real (kind=RealKind), parameter :: au2anstr = 0.52917000d0
   real (kind=RealKind) :: uconv, fact, a0
!
!  data a0/3.8525, 3.8525, 3.7133/
!  data cut/19.3/
!
   data NumBasis(1:NumLattices)/4, 2, 1, 16/
   data BasisVec(1:3,1:4,1) &
      & /ZERO, ZERO, ZERO, HALF, HALF, ZERO, HALF, ZERO, HALF, ZERO, HALF, HALF/
   data BasisVec(1:3,1:4,2) &
      & /ZERO, ZERO, ZERO, HALF, HALF, HALF, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
   data BasisVec(1:3,1:4,3) &
      & /ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO/
   data BasisVec(1:3,1:16,4) &
      & /ZERO, ZERO, ZERO,   0.25d0, 0.25d0, 0.25d0, &
      &  HALF, HALF, ZERO,   0.75d0, 0.75d0, 0.25d0, &
      &  HALF, ZERO, HALF,   0.75d0, 0.25d0, 0.75d0, &
      &  ZERO, HALF, HALF,   0.25d0, 0.75d0, 0.75d0, &
      &  HALF, HALF, HALF,   0.75d0, 0.75d0, 0.75d0, &
      &  ZERO, ZERO, HALF,   0.25d0, 0.25d0, 0.75d0, &
      &  ZERO, HALF, ZERO,   0.25d0, 0.75d0, 0.25d0, &
      &  HALF, ZERO, ZERO,   0.75d0, 0.25d0, 0.25d0/
!
   interface
      subroutine readPositionData(fname,NumAtomsIn,NumAtomsOut)
         use KindParamModule, only : IntKind
         character (len=*), intent(in) :: fname
         integer (kind=IntKind), intent(in), optional :: NumAtomsIn
         integer (kind=IntKind), intent(out), optional :: NumAtomsOut
      end subroutine readPositionData
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
   HostBasis(:)  = '   '
   ClusterBasis(:) = '   '
   a0 = ONE
   call initString('Arbitrary')
!
!  ==========================================================================
!  read input parameters for setting up the big box and the underlying lattice
!  ==========================================================================
   write(6,'(//,a)')  &
      'Please Enter the Underlying Lattice Information as Follows ---->'
   lattice = -1
   do while (lattice < 0 .or. lattice > 5)
      write(6,'(2x,a)')'For the following options:'
      write(6,'(2x,a)')'    1. Face Centered;'
      write(6,'(2x,a)')'    2. Body Centered;'
      write(6,'(2x,a)')'    3. Orthorhombic;'
      write(6,'(2x,a)')'    4. Diamond;'
      write(6,'(2x,a)')'    5. Read Underlying Lattice Data from VASP POSCAR File'
      write(6,'(2x,a)')'    0. Read Underlying Lattice Data from MST Position File'
      write(6,'(2x,a,$)')'Enter your your choice: '
      read(5,*)lattice
   enddo
!
   iconv = -1
   write(6,'(//,a)')  &
      'Please Enter the Units for Input Lattice Constants ---->'
   do while (iconv /= 0 .and. iconv /= 1)
      write(6,'(2x,a,$)') 'Choose (0. Atomic Units;  1. Angstrom): '
      read(5,*)iconv
   enddo
   if (iconv == 0) then
       uconv = ONE
   else
       uconv = anstr2au
   endif
!
   if (lattice == 0) then
      file_exist = .false.
      do while (.not.file_exist)
         write(6,'(/,2x,a,$)') 'Name of the MST position data file: '
         read(5,'(a)')file_name
         file_name = adjustl(file_name)
         inquire(file=file_name,exist=file_exist)
         if (.not.file_exist) then
            write(6,'(/,2x,a,$)') 'The file just entered does not exist. Try again...'
         endif
      enddo
!     ----------------------------------------------------------------
      call initDataServiceCenter()
      call readPositionData(file_name,NumAtomsOut=nbasis)
      Bravais => getDataStorage('Bravais Vector',3,3,RealMark)  
      a0 = getDataStorage('Position Scaling Factor',RealMark)
      AtomPosition => getDataStorage('Atomic Position',3,nbasis,RealMark)
!     ----------------------------------------------------------------
      if (nbasis < 1 .or. nbasis > MaxBasis) then
         call ErrorHandler('main','Number of basis is out of range',nbasis)
      endif
      small_box = Bravais/a0
      do i=1, nbasis
         bv(1,i) = AtomPosition(1,i)/a0
         bv(2,i) = AtomPosition(2,i)/a0
         bv(3,i) = AtomPosition(3,i)/a0
      enddo
      nullify(Bravais,AtomPosition)
!     ----------------------------------------------------------------
      call endDataServiceCenter()
!     ----------------------------------------------------------------
   else if (lattice == 5 ) then
      write(6,'(/,2x,a,$)') 'Name of the Underlying Lattice File (e.g., POSCAR): '
      read(5,'(a)')file_name
      file_name = adjustl(file_name)
      open(unit=11,file=file_name,form='formatted',status='old')
      read(11,'(a)')text
      read(11,*)fact
      read(11,*)small_box(1:3,1)
      read(11,*)small_box(1:3,2)
      read(11,*)small_box(1:3,3)
      uconv = uconv*fact
      small_box(1:3,1:3) = small_box(1:3,1:3)*uconv
      read(11,'(a)')text
      call setString(text)
      text = adjustl(text)
      na = getNumTokens()
      NumHostAtomTypes = na
      read(text,*) (NumHostAtomsOfType(i),i=1,na)
      write(6,'(/,2x,a,$)') 'The atom names (separated by spaces and/or comma) in the order defined in POSCAR (e.g., Cu Al Fe): '
      read(5,'(a)')text
      text = adjustl(text)
      call setString(text)
      na = getNumTokens()
      if (na /= NumHostAtomTypes) then
         call ErrorHandler('main','na <> NumHostAtomTypes',na,NumHostAtomTypes)
      endif
      j = 1
      do i = 1,na
         call readToken(i,anm,n)
         if (n >= 3) then
            Host(i) = anm(1:2)
         else
            Host(i) = anm
         endif
         HostBasis(j:j+NumHostAtomsOfType(i)-1) = Host(i)
         j = j + NumHostAtomsOfType(i)
      enddo
!     text = adjustl(text)
      nbasis = 0
      do i = 1,na
         nbasis = nbasis + NumHostAtomsOfType(i)
      enddo
      if (nbasis < 1 .or. nbasis > MaxBasis) then
         call ErrorHandler('main','Number of basis is out of range',nbasis)
      endif
      read(11,'(a)')text
      do i = 1,nbasis
         read(11,'(a)')text
         call setString(text)
         text = adjustl(text)
         na = getNumTokens()
         if (na == 3) then
            read(text,*) rpos(1:3)
         else
            call ErrorHandler('main','Invalid input data',trim(text))
         endif
         do j = 1,3
            bv(j,i) = small_box(j,1)*rpos(1)+small_box(j,2)*rpos(2)+small_box(j,3)*rpos(3)
         enddo
      enddo
      close(11)
   else
      small_box(1:3,1:3) = ZERO
      ios = 1
      do while (ios /= 0)
         if (iconv == 0) then
            if (lattice == 3) then
               write(6,'(/,2x,a,$)')     &
               'Lattice Constants a, b, c (in a.u. seperated by space or comma): '
            else
               write(6,'(/,2x,a,$)')     &
               'Lattice Constants a (in a.u.): '
            endif
         else
            if (lattice == 3) then
               write(6,'(/,2x,a,$)')     &
               'Lattice Constants a, b, c (in Angstrom seperated by space or comma): '
            else
               write(6,'(/,2x,a,$)')     &
               'Lattice Constants a (in Angstrom): '
            endif
         endif
         if (lattice == 3) then
            read(5,*,iostat=ios)small_box(1,1),small_box(2,2),small_box(3,3)
         else
            read(5,*,iostat=ios) a0
            small_box(1,1) = a0
            small_box(2,2) = a0
            small_box(3,3) = a0
         endif
         if (small_box(1,1) < TEN2m6 .or. small_box(2,2) < TEN2m6 .or. &
             small_box(3,3) < TEN2m6) then
            ios = 1
         endif
         if (ios /= 0) then
            write(6,'(a)')'Invalid input! Try again...'
         endif
      enddo
!     ======================================================================
!     small_box is an orthorombic box
!     ======================================================================
      a0 = small_box(1,1)
      small_box(1,1) = small_box(1,1)*uconv/a0
      small_box(2,2) = small_box(2,2)*uconv/a0
      small_box(3,3) = small_box(3,3)*uconv/a0
      nbasis = NumBasis(lattice)
      do i=1, nbasis
         bv(1,i) = BasisVec(1,i,lattice)*small_box(1,1)
         bv(2,i) = BasisVec(2,i,lattice)*small_box(2,2)
         bv(3,i) = BasisVec(3,i,lattice)*small_box(3,3)
      enddo
   endif
!
   write(6,'(/,2x,a,$)')       &
           'Number of Repeats of Small Box Along A, B, C Directions: '
   read(5,*)na,nb,nc
!
   write(6,'(/,2x,a,$)')       &
           'Type of output format( 0. generic(x,y.z); 1. i_bigcell; 2. VASP (Direct) ): '
   read(5,*) ftype
!
!  --------------------------------------------------------------------------
   call initSample(nbasis,na,nb,nc,small_box,bv)
   box => getUnitCell()
   NumAtoms = getNumAtoms()
!  --------------------------------------------------------------------------
!
   call MtxInv_GE(3,box,box_inv)
!
!  ==========================================================================
!  read input parameters for the solid solution occupying the lattice
!  ==========================================================================
   ordered = -1; rloop = 0
   do while (ordered < 0 .or. ordered > 2)
      write(6,'(//,a)')        &
              'Please Enter the Solid Solution Information as Follows ---->'
!
      write(6,'(2x,a)') 'Is the solid solution ordered or random?'
      write(6,'(4x,a,$)') 'Enter 0 for ORDERED; 1 for RANDOM; or 2 for RANDOM with SHORT-RANGE ORDER: '
      read(5,*)ordered
      if (ordered == 0) then
         RandomHost = .false.
      else if (ordered == 1 .or. ordered == 2) then
         RandomHost = .true.
      else if (rloop <= 5) then
         write(6,'(a,i3)')'Undefined input: ',ordered
         rloop = rloop + 1
      else
         call ErrorHandler('main','Undefined input','Too many tries')
      endif
   enddo
!
   if (.not.RandomHost .and. HostBasis(1) == '   ') then
      if (lattice /= 4) then
         do i = 1, nbasis
            write(6,'(/,2x,a,3f10.5,a,$)')'Enter Atom Name Located at',   &
                                        bv(1:3,i)*uconv,' : '
            read(5,'(a)')HostBasis(i)
         enddo
      else
         do i = 1, 2
            write(6,'(/,2x,a,3f10.5,a,$)')'Enter Atom Name Located at',   &
                                           bv(1:3,i)*uconv,' : '
            read(5,'(a)')HostBasis(i)
         enddo
         do i = 3, 8, 2
            HostBasis(i) = HostBasis(1)
            HostBasis(i+1) = HostBasis(2)
         enddo
         do i = 9, 16
            HostBasis(i) = 'Va'
         enddo
      endif
      k = 1
      NumHostAtomTypes = 1
      Host(1) = HostBasis(1)
      do j = 2,nbasis
         NumHostAtomTypes = NumHostAtomTypes + 1
         LoopZ1: do i = j-1,1,-1
            if ( HostBasis(j) == HostBasis(i) ) then
               NumHostAtomTypes = NumHostAtomTypes -1
               exit LoopZ1
            endif
         enddo LoopZ1
         if ( i==0 ) then
            k = k+1
            Host(k) = HostBasis(j)
         endif
      enddo
!     -----------------------------------------------------------------------
      call placeAtoms(nbasis,HostBasis)
!     -----------------------------------------------------------------------
   else if (RandomHost) then
!     -----------------------------------------------------------------------
      NumHostAtomTypes = 0
      do while (NumHostAtomTypes < 1 .or. NumHostAtomTypes > MaxAtomTypes)
         write(6,'(/,2x,a,$)')     &
          'Constituents (atomic name separated by space or comma, e.g., Fe Ni): '
         read(5,'(a)')text
         call setString(text)
         NumHostAtomTypes = getNumTokens()
         if (NumHostAtomTypes < 1 ) then
            call WarningHandler('main','Number of atom types < 1',NumHostAtomTypes)
            write(6,'(a)')'Try again ...'
         else if (NumHostAtomTypes > MaxAtomTypes) then
            call WarningHandler('main','Number of atom types > limit',NumHostAtomTypes,MaxAtomTypes)
            write(6,'(a)')'Try again ...'
         endif
      enddo
!     -----------------------------------------------------------------------
      done = .false.
      do while (.not.done)
         do i = 1, NumHostAtomTypes
            call readToken(i,Host(i),alen)
            write(6,'(/,2x,3a,$)')'Content of ',Host(i),' (>= 0 and =< 1): '
            read(5,*)HostContent(i)
         enddo
         write(6,'(a)')' '
!
         bins(1) = ZERO
         do i = 2, NumHostAtomTypes+1
            bins(i) = bins(i-1) + HostContent(i-1)
         enddo
         if (abs(bins(NumHostAtomTypes+1)-ONE) > TEN2m6) then
            call WarningHandler('main','The summation of the contents is not 1',bins(NumHostAtomTypes+1))
            write(6,'(a)')'Try again ...'
         else
            done = .true.
         endif
      enddo
!     -----------------------------------------------------------------------
      iseed = 0
      write(6,'(/,2x,a,$)') 'Initialize random seed by entering an integer (> 0), other wise using the default seed: '
      read(5,'(i10)',iostat=ios) iseed
      if (ios == 0) then
         if (iseed > 0) then
            write(6,'(/,2x,a,i10)')'Initial seed for random number generator is: ',iseed
         endif
      else
         iseed = 0
      endif
      if (iseed < 1) then
         write(6,'(/,2x,a)')'Initial seed for random number generator will be the default.'
      endif
!     -----------------------------------------------------------------------
!
      nshell=-1
      if (ordered == 2) then   ! Random host with short range order
!        --------------------------------------------------------------------
         done = .false.
         do while (.not.done)
            write(6,'(/,2x,a,$)') 'How many shells will be considered (e.g. 4): '
            read(5,*)nshell
            if (nshell < 1 ) then
               call WarningHandler('main','Number of shells < 1',nshell)
            else if (nshell > MaxShells) then
               call WarningHandler('main','Number of shells > MaxShells',nshell,MaxShells)
            else
               done = .true.
            endif
         enddo
         write(6,'(2x,a,i3)') 'nshell is:', nshell
!        --------------------------------------------------------------------
         done = .false.
         do while (.not.done)
            write(6,'(/,2x,a,$)') 'For each shell please give the weight of SROs (e.g. 1.0, 1.0, 1.0, 1.0): '
            ib=-1
            read(5,'(a)')text
            call setString(text)
            ib = getNumTokens()
!
            if (ib /= nshell) then
               call WarningHandler('main','Number of weight of SRO not correct, should be',nshell)
            else
               done = .true.
            endif
         enddo
         do j = 1, nshell
           call readToken(j,anm,alen)
           read(anm,*)weight(j)
         enddo
         write(6,'(a)')' '
         write(6,'(2x,a,5f15.5)') 'The input weight is:',weight(1:nshell)

!        ------------------------------------------------------------------------
         write(6,'(/,2x,a,$)') 'For each shell please give N*(N-1)/2 SROs, where N is the number of species on the shell:'
         do i = 1, nshell
            ib=-1
            do while (ib /= NumHostAtomTypes*(NumHostAtomTypes-1)/2)
               if (NumHostAtomTypes == 2) then
                  write(6,'(/,2x,a,i2,a,$)') 'shell',i,' (1 number is expected):  ' 
               else
                  write(6,'(/,2x,a,i2,a,i2,a,$)') 'shell',i,' (',NumHostAtomTypes,' numbers are expected):  ' 
               endif
!              --------------------------------------------------------------
               read(5,'(a)')text
               call setString(text)
               ib = getNumTokens()
!              --------------------------------------------------------------
               if (ib /= NumHostAtomTypes*(NumHostAtomTypes-1)/2) then
                  call WarningHandler('main','Number of SRO not correct, should be',NumHostAtomTypes*(NumHostAtomTypes-1)/2)
               endif
            enddo

            do j = 1, NumHostAtomTypes*(NumHostAtomTypes-1)/2
              call readToken(j,anm,alen)
              read(anm,*)srop(j,i)
            enddo 
         enddo

         write(6,'(/,2x,a)')'the input SROs are: '
         do i=1,nshell
           do j=1, NumHostAtomTypes*(NumHostAtomTypes-1)/2
              write(6,'(2x,f15.5,$)') srop(j,i)
           enddo
           write(6,'(a)')' '
         enddo 

         write(6,'(/,2x,a,$)')'Enter the temperature (K) to start annealing: '
         read(5,*)Tmax
         write(6,'(/,2x,a,$)')'Enter the temperature step (K) for cooling: '
         read(5,*)Tstep
         write(6,'(a)')' '
         if (iseed > 0) then
!           -----------------------------------------------------------------
            call placeAtoms(NumHostAtomTypes,Host,HostContent,weight,nshell,srop,Tmax,Tstep,iseed)
!           -----------------------------------------------------------------
         else
!           -----------------------------------------------------------------
            call placeAtoms(NumHostAtomTypes,Host,HostContent,weight,nshell,srop,Tmax,Tstep)
!           -----------------------------------------------------------------
         endif
      else
         if (iseed > 0) then
!           -----------------------------------------------------------------
            call placeAtoms(NumHostAtomTypes,Host,HostContent,iseed)
!           -----------------------------------------------------------------
         else
!           -----------------------------------------------------------------
            call placeAtoms(NumHostAtomTypes,Host,HostContent)
!           -----------------------------------------------------------------
         endif
      endif
   endif
!
!  ==========================================================================
!  The following code takes care the situation where there are impurities or
!  clusters to substitute the atoms in the sample.
!  ==========================================================================
   substitution = -1; rloop = 0
   NumImpuritySpecies = 0
   NumClusters = 0
   LOOP_dowhile: do while (substitution < 0 .or. substitution > 1) 
      write(6,'(//,a)')'Do you need to modify the sameple?'
      write(6,'(2x,a)')'Choose  0 No modification;'
      write(6,'(2x,a)')'        1 Subsituting atoms by impurities;'
      write(6,'(2x,a)')'        2 Embedding clusters into the sample (not yet available).'
      write(6,'(2x,a,$)')'Enter your choice: '
      read(5,*)substitution
      if (substitution == 0) then
         exit LOOP_dowhile
      else if (substitution == 1) then
         write(6,'(/,2x,a)')   'The host sample consists of:'
         do ib = 1, getNumSpecies()
            write(6,'(8x,i8,2x,a,a)')getNumAtoms(Species=ib), &
                                     getAtomName(ib,SpeciesIndex=.true.),' atoms'
         enddo
         write(6,'(/,2x,a,$)')'Impurity species (atomic name separated by space or comma, e.g., Fe Ni): '
         read(5,'(a)')text
!        ------------------------------------------------------------------------------
         call setString(text)
         NumImpuritySpecies = getNumTokens()
!        ------------------------------------------------------------------------------
         do ib = 1, NumImpuritySpecies
!           ---------------------------------------------------------------------------
            call readToken(ib,Impurity(ib),alen)
!           ---------------------------------------------------------------------------
            write(6,'(/,2x,3a,$)')'Content of ',Impurity(ib),' (>= 0 and =< 1): '
            read(5,*)ImpurityContent(ib)
            do
               write(6,'(2x,a,$)')'Which host atom it replaces? '
               read(5,'(a)')anm
               anm = adjustl(anm)
               if (isSampleAtom(anm)) then
                  exit
               else
                  write(6,'(/,2x,2a)')trim(anm),' is not a host atom species. Try again...'
               endif
            enddo
            ReplacedHost(ib) = anm(1:2)
         enddo
!        ------------------------------------------------------------------------------
         call substituteAtoms(NumImpuritySpecies,Impurity,ImpurityContent,ReplacedHost)
!        ------------------------------------------------------------------------------
      else if (substitution == 2) then
         write(6,'(/,2x,a)')'Not yet implemented...'
         rloop = rloop + 1
         cycle
         write(6,'(/,2x,a,$)')'Number of Clusters: '
         read(5,*)NumClusters
         allocate( ClusterFlag(NumAtoms), b_cluster(NumAtoms) )
         allocate( x_cluster(NumAtoms), y_cluster(NumAtoms),            &
                   z_cluster(NumAtoms), AtomName_cluster(NumAtoms),     &
                   AtomZ_Cluster(NumAtoms), IndexN_Cluster(NumAtoms) )
      else if (rloop <= 5) then
         write(6,'(a,i3)')'Invalid input value: ',substitution
         rloop = rloop + 1
      else
         call ErrorHandler('main','Invalid input value',substitution)
      endif
   enddo LOOP_dowhile
!
   if (NumClusters < 0 .or. NumClusters > MaxClusters) then
      call ErrorHandler('main','Number of clusters is out of range',NumClusters)
   endif
!
   NumAtomSpecies = getNumSpecies()
!
!  ===================================================================
!  Output the head lines of the position data
!  ===================================================================
   file_exist = .true.
   do while (file_exist)
      write(6,'(/,2x,a,$)')'Name of output position data file: '
      read(5,'(a)')file_name
      inquire(file=file_name,exist=file_exist)
      if (file_exist) then
         write(6,'(/,2x,a)')'A file with the name just entered already exists. Try a different name...'
      endif
   enddo
   open(unit=14,file=file_name,status='new',form='formatted')
   fact = ONE
   if (ftype==2) then
      write(14,'(a)') '#POSCAR file'
      write(14,*) a0
      fact = au2anstr
   else
      write(14,'(a)') '# Default units:: atomic units'
      write(14,'(a)') '# '
      write(14,'(a)') '# Uncomment following number if using Angsrtroms units'
      write(14,'(a)') '# =================================================='
      write(14,'(a,f12.8)')"# ", au2anstr*a0
      write(14,'(a)') '# '
      write(14,'(a)') '# If using Angstroms units, comment out the following number'
      write(14,'(a)') '# =================================================='
      write(14,'(f12.8)')a0
      write(14,'(a)') '# '
      write(14,'(a)') '# The following three lines define the unit cell box'
      write(14,'(a)') '# =================================================='
   endif
   write(14,'(2x,3f19.11)')fact*box(1:3,1)
   write(14,'(2x,3f19.11)')fact*box(1:3,2)
   write(14,'(2x,3f19.11)')fact*box(1:3,3)
   if (ftype/=2) then
      write(14,'(a)') '# '
      write(14,'(a,i8)')    '# Number of atoms in unit cell:  ',NumAtoms
      do ib = 1, NumAtomSpecies
         write(14,'(3a,i8)')'# Atom type: ',getAtomName(ib,SpeciesIndex=.true.), &
                            ' Number of Atoms: ',getNumAtoms(ib)
      enddo
      write(14,'(a)') '# '
      write(14,'(a)') '# The following lines are the position data of the atoms'
      write(14,'(a)') '# ======================================================'
   endif
!
if (.false.) then
   NumHostAtoms = NumAtoms
!  =============
!  The following codes are used for inserting a cluster. We will work on it in
!  the future.
!  =============
   if (ftype/=2) then
      write(14,'(a,i8)')'# Number of clusters: ',NumClusters
   endif
!
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
!  Output the host atom position data
!  ===================================================================
   do ib = 1, NumHostAtomTypes
      NumHostAtomsOfType(ib) = 0
      do i = 1,NumAtoms
         if ( Host(ib) == getAtomName(i) ) then
            if ( substitution==2 ) then
               if ( ClusterFlag(i) == 0 ) then
                   NumHostAtomsOfType(ib) = NumHostAtomsOfType(ib)+1
               endif
            else 
               NumHostAtomsOfType(ib) = NumHostAtomsOfType(ib)+1
            endif
         endif
      enddo
      if (ftype/=2) then
         write(14,'(a,a3,1x,'': '',i8)')'# Number of ',Host(ib),       &
                           NumHostAtomsOfType(ib)
      endif
   enddo
endif
!
   if ( ftype==2 ) then
      do ib = 1, NumAtomSpecies
         write(14,'(a6,$)')getAtomName(ib,SpeciesIndex=.true.)
      enddo
      write(14,'(a)')' '
      do ib = 1, NumAtomSpecies
         write(14,'(i5,$)')getNumAtoms(ib)
      enddo
      write(14,'(/,a)')'Direct'
   endif
!
   allocate( AtomZ(NumAtoms), IndexN(NumAtoms) )
   do i = 1,NumAtoms
      AtomZ(i) = getZtot(getAtomName(i))
   enddo
!  -------------------------------------------------------------------
   call HeapSort(NumAtoms, AtomZ, IndexN)
!  -------------------------------------------------------------------
!
   x_host => getAtomPosition(1)
   y_host => getAtomPosition(2)
   z_host => getAtomPosition(3)
!
   if (substitution < 2) then
      do i = 1, NumAtoms
         k = IndexN(i)
         if ( ftype==1 ) then
            n = AtomZ(i)
            write(14,'(i3,1x,3f19.15,1x,f7.4,2x,i5,3x,a)') n,      &
                x_host(k),y_host(k),z_host(k), 1.0d0, 0, "U"
         else if ( ftype==2) then
            do j = 1,3
               rpos(j) = box_inv(j,1)*x_host(k) +                &
                         box_inv(j,2)*y_host(k) +                &
                         box_inv(j,3)*z_host(k)
            enddo
            write(14,'(2x,a3,2x,3f19.11)')getAtomName(k), rpos(1:3)
            !write(14,'(2x,3f19.11)') rpos(1:3)
         else
            write(14,'(2x,a3,2x,3f19.11)') getAtomName(k),     &
                                x_host(k),y_host(k),z_host(k)
         endif
      enddo
   else
      do i = 1, NumHostAtoms
         k = IndexN(i)
         if (ClusterFlag(i) == 0) then
            if ( ftype==1 ) then
               n = AtomZ(i)
               write(14,'(i3,1x,3f19.15,1x,f7.4,2x,i5,3x,a)') n,      &
                   x_host(k),y_host(k),z_host(k), 1.0d0, 0, "U"
            else if ( ftype==2) then
               do j = 1,3 
                  rpos(j) = box_inv(j,1)*x_host(k) +                &
                            box_inv(j,2)*y_host(k) +                &
                            box_inv(j,3)*z_host(k)
                  enddo
               write(14,'(2x,a3,2x,3f19.11)')getAtomName(k), rpos(1:3)
            else
               write(14,'(2x,a3,2x,3f19.11)') getAtomName(k),     &
                                   x_host(k),y_host(k),z_host(k)
            endif
         endif
      enddo
      deallocate( ClusterFlag, b_cluster )
      deallocate( x_cluster, y_cluster, z_cluster, AtomName_cluster )
      deallocate( AtomZ_Cluster, IndexN_Cluster )
   endif
!
   close(14)
!
   write(6,*) " "
!
   nullify( x_host, y_host, z_host )
   deallocate( AtomZ, IndexN )
!
   call endString()
!
!  -------------------------------------------------------------------
   call endSample()
!  -------------------------------------------------------------------
!
   end program generateAtomPosition
