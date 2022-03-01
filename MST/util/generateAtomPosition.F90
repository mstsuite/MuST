!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   program generateAtomPosition
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, HALF, ONE, TWO, TEN2m6, TEN2m10
!
   use ErrorHandlerModule, only : ErrorHandler
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
!
   character (len=2) :: Cluster(MaxAtomTypes)
   character (len=2) :: Medium(MaxAtomTypes)
   character (len=2) :: MediumBasis(MaxBasis)
   character (len=2) :: ClusterBasis(MaxBasis)
   character (len=60) :: text, file_name, anm
!
   character (len=2), pointer :: AtomName_medium(:)
   character (len=2), pointer :: AtomName_cluster(:)
!
   integer (kind=IntKind) :: alen, ClusterShape, lattice, NumBounds, iconv
   integer (kind=IntKind) :: NumBasis(3), nbasis, na, nb, nc, NumAtoms
   integer (kind=IntKind) :: i, j, k, ib, n, ncl
   integer (kind=IntKind) :: NumMediumAtomTypes, NumClusterAtomTypes
   integer (kind=IntKind) :: NumClusters
   integer (kind=IntKind) :: NumMediumAtoms, NumClusterAtoms(MaxClusters)
   integer (kind=IntKind) :: ordered, embed, rloop
   integer (kind=IntKind) :: ftype
   integer (kind=IntKind) :: nshell
!
   integer (kind=IntKind) :: NumMediumAtomsOfType(MaxAtomTypes)
   integer (kind=IntKind) :: NumClusterAtomsOfType(MaxAtomTypes)
   integer (kind=IntKind), allocatable :: ClusterFlag(:)
   integer (kind=IntKind), allocatable :: b_cluster(:)
   integer (kind=IntKind), allocatable :: AtomZ(:)
   integer (kind=IntKind), allocatable :: IndexN(:)
   integer (kind=IntKind), allocatable :: AtomZ_Cluster(:)
   integer (kind=IntKind), allocatable :: IndexN_Cluster(:)
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
   real (kind=RealKind) :: uconv, fact
! GDS   
   real (kind=RealKind) :: work(3)
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
  'Choose (1. Face Centered;  2. Body Centered;  3. Orthorhombic;  4. VASP positions:  0. Other): '
      read(5,*)lattice
   enddo
   write(6,'(i3)')lattice
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
   if (lattice == 0) then
      write(6,'(/,2x,a,$)') 'Name of the Underlying Lattice File: '
      read(5,'(a)')file_name
      write(6,'(a)')trim(adjustl(file_name))
      file_name = adjustl(file_name)
      open(unit=11,file=file_name,form='formatted',status='old')
      read(11,*)small_box(1:3,1)
      read(11,*)small_box(1:3,2)
      read(11,*)small_box(1:3,3)
      small_box(1:3,1:3) = small_box(1:3,1:3)*uconv
      read(11,*)nbasis
      if (nbasis < 1 .or. nbasis > MaxBasis) then
         call ErrorHandler('main','Number of basis is out of range',nbasis)
      endif
      do i=1, nbasis
         read(11,'(a)')text
         call initString(text)
         text = adjustl(text)
         na = getNumTokens()
         if (na == 4) then
            read(text(1:3),'(a)')MediumBasis(i)
            read(text(4:),*)bv(1:3,i)
         else if (na == 3) then
            read(text,*)bv(1:3,i)
         else
            call ErrorHandler('main','Invalid input data',trim(text))
         endif
         bv(1:3,i) = bv(1:3,i)*uconv
         call endString()
      enddo
      close(11)
   else if (lattice == 4 ) then
      write(6,'(/,2x,a,$)') 'Name of the Underlying Lattice File (e.g., POSCAR): '
      read(5,'(a)')file_name
      file_name = adjustl(file_name)
      write(6,'(a)')trim(file_name)
      open(unit=11,file=file_name,form='formatted',status='old')
      read(11,'(a)')text
      read(11,*)fact
      read(11,*)small_box(1:3,1)
      read(11,*)small_box(1:3,2)
      read(11,*)small_box(1:3,3)
      uconv = uconv*fact
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
         do j = 1,3
            bv(j,i) = small_box(j,1)*rpos(1)+small_box(j,2)*rpos(2)+small_box(j,3)*rpos(3)
         enddo
      enddo
      close(11)
   else
      small_box(1:3,1:3) = ZERO
      if (iconv == 0) then
         write(6,'(/,2x,a,$)')     &
         'Lattice Constants a, b, c (in a.u. seperated by space or comma): '
      else
         write(6,'(/,2x,a,$)')     &
         'Lattice Constants a, b, c (in Angstrom seperated by space or comma): '
      endif
      read(5,*)small_box(1,1),small_box(2,2),small_box(3,3)
      small_box(1,1) = small_box(1,1)*uconv
      small_box(2,2) = small_box(2,2)*uconv
      small_box(3,3) = small_box(3,3)*uconv
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
   write(6,'(3i5)')na,nb,nc
!
   write(6,'(/,2x,a,$)')       &
           'Type of output format( 0. generic(x,y.z); 1. i_bigcell; 2. VASP (Direct) ): '
   read(5,*) ftype
   write(6,'(i3)') ftype
!
!  --------------------------------------------------------------------------
   call initSample(nbasis,na,nb,nc,small_box,bv)
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
!
   if (.not.RandomMedium .and. MediumBasis(1) == '   ') then
      do i = 1, nbasis
         write(6,'(/,2x,a,3f10.5,a,$)')'Enter Atom Name Located at',   &
                                     bv(1:3,i)*uconv,' : '
         read(5,'(a)')MediumBasis(i)
         write(6,'(a)')trim(adjustl(MediumBasis(i)))
      enddo
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
!     --------------------------------------------------------------------
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
      write(6,'(2x,a,$)')'Enter 0 for YES; or 1 for NO: '
      read(5,*)embed
      write(6,'(i3)')embed
      if (embed == 0) then
         write(6,'(/,2x,a,$)')'Number of Clusters: '
         read(5,*)NumClusters
         allocate( ClusterFlag(NumAtoms), b_cluster(NumAtoms) )
         allocate( x_cluster(NumAtoms), y_cluster(NumAtoms),            &
                   z_cluster(NumAtoms), AtomName_cluster(NumAtoms),     &
                   AtomZ_Cluster(NumAtoms), IndexN_Cluster(NumAtoms) )
      else if (embed == 1) then
         NumClusters = 0
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
      write(14,*) 1.0d0
      fact = au2anstr
   else
      write(14,'(a)') '# Units:: atomic units( uncomment next line for Angsrtroms )'
      write(14,'(a,f12.8)')"# ", au2anstr
   endif
   write(14,'(2x,3f19.11)')fact*box(1:3,1)
   write(14,'(2x,3f19.11)')fact*box(1:3,2)
   write(14,'(2x,3f19.11)')fact*box(1:3,3)
   if (ftype/=2) then
      write(14,'(a,i8)')'# Number of clusters: ',NumClusters
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
            if ( embed==0 ) then
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
   enddo
!  -------------------------------------------------------------------
   call HeapSort(NumAtoms, AtomZ, IndexN)
!  -------------------------------------------------------------------
!
   if (embed == 0) then
      do i = 1, NumMediumAtoms
         k = IndexN(i)
         if (ClusterFlag(i) == 0) then
            if ( ftype==1 ) then
               n = AtomZ(i)
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
      do i = 1, NumMediumAtoms
         k = IndexN(i)
         if ( ftype==1 ) then
            n = AtomZ(i)
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
