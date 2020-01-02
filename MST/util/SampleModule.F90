module SampleModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initSample,        &
          endSample,         &
          placeAtoms,        &
          substituteAtoms,   &
          getAtomPosition,   &
          getAtomName,       &
          getUnitCell,       &
          getNumSpecies,     &
          getNumAtoms,       &
          isSampleAtom
!
   interface placeAtoms
      module procedure placeAtomsOrderly, placeAtomsRandomly, placeAtomsRandomlyWithSRO
   end interface
!
   interface getAtomPosition
      module procedure getAtomPosition0, getAtomPosition1
   end interface
!
   interface getAtomName
      module procedure getAtomName0, getAtomName1
   end interface
!
   integer (kind=IntKind), parameter, public :: MaxAtomTypes = 20
   integer (kind=IntKind), parameter, public :: MaxShells = 20
!
private
   character (len=2), target, allocatable :: AtomName(:)
   character (len=2), target, allocatable :: TypeName(:)
!
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: NumTypes
   integer (kind=IntKind) :: nshell      ! input shell number
   integer (kind=IntKind) :: NumRepeatsA, NumRepeatsB, NumRepeatsC
   integer (kind=IntKind) :: NumShell = 0   ! calculated max shell number for the supercell
   integer (kind=IntKind) :: iseed = -30
   integer (kind=IntKind), allocatable :: AtomType(:)
   integer (kind=IntKind), allocatable :: NumAtomsOfType(:)
!
   real (kind=RealKind), target :: large_cell(3,3)
   real (kind=RealKind), target, allocatable :: AtomPosition(:,:)
   real (kind=RealKind), pointer :: atom_position_x(:)
   real (kind=RealKind), pointer :: atom_position_y(:)
   real (kind=RealKind), pointer :: atom_position_z(:)
!
!  ===================================================================
!  dr_shell: the atoms within (r,r+dr_shell) in shell radius are considered
!            on the same shell
!  ===================================================================
   real (kind=RealKind), parameter :: dr_shell = 0.001
!
   type PairStruct
      integer (kind=IntKind) :: npairs
      real (kind=RealKind), pointer :: srop(:,:)   ! array of SRO for each shell, MaxShells x Ntype*(Ntype-1)/2
      real (kind=RealKind)   :: rpair2
      logical :: isAllocated
      integer (kind=IntKind), pointer :: i_of_pair(:)  ! stores the index of AtomName
      integer (kind=IntKind), pointer :: j_of_pair(:)  ! stores the index of AtomName
      integer (kind=IntKind), pointer :: npairs_matrix(:,:)
   end type PairStruct
   type (PairStruct) :: sr_pairs(MaxShells)
!
   logical :: Initialized = .false.
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSample(nbasis,na,nb,nc,small_cell,bv,ct)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nbasis, na, nb, nc
!
   real (kind=RealKind), target, intent(in) :: small_cell(3,3)
   real (kind=RealKind), target, intent(in) :: bv(3,nbasis)
!
   logical, optional, intent(in) :: ct
   logical :: centered
!
   integer (kind=IntKind) :: i
!
   if (na < 1) then
      call ErrorHandler('initSample','Invalid parameter: na',na)
   else if (nb < 1) then
      call ErrorHandler('initSample','Invalid parameter: nb',nb)
   else if (nc < 1) then
      call ErrorHandler('initSample','Invalid parameter: nc',nc)
   else if (nbasis < 1) then
      call ErrorHandler('initSample','Invalid parameter: nbasis',nbasis)
   endif
!
   if ( present(ct) ) then
      centered = ct
   else
      centered = .false.
   endif
!
   NumAtoms = nbasis*na*nb*nc
   NumRepeatsA = na
   NumRepeatsB = nb
   NumRepeatsC = nc
!
   allocate( AtomPosition(NumAtoms,3), AtomName(1:NumAtoms), AtomType(1:NumAtoms) )
!
   large_cell(1:3,1) = small_cell(1:3,1)*na
   large_cell(1:3,2) = small_cell(1:3,2)*nb
   large_cell(1:3,3) = small_cell(1:3,3)*nc
!
   atom_position_x => AtomPosition(1:NumAtoms,1)
   atom_position_y => AtomPosition(1:NumAtoms,2)
   atom_position_z => AtomPosition(1:NumAtoms,3)
!
!  -------------------------------------------------------------------
   call genLatticePoints(nbasis,na,nb,nc,small_cell,bv,centered)
!  -------------------------------------------------------------------
!
   do i = 1, MaxShells
      sr_pairs(i)%npairs = 0
      sr_pairs(i)%isAllocated = .false.
   enddo
!
   Initialized = .true.
!
   end subroutine initSample
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSample()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i
!
   nullify( atom_position_x,atom_position_y,atom_position_z )
!
   deallocate( AtomName, AtomPosition, AtomType )
!
   if (allocated(NumAtomsOfType)) then
      deallocate(NumAtomsOfType)
   endif
   if (allocated(TypeName)) then
      deallocate(TypeName)
   endif
!
   do i = 1, MaxShells
      sr_pairs(i)%npairs = 0
      if (sr_pairs(i)%isAllocated) then
         deallocate( sr_pairs(i)%npairs_matrix, sr_pairs(i)%i_of_pair, &
                     sr_pairs(i)%j_of_pair )
         nullify(  sr_pairs(i)%npairs_matrix, sr_pairs(i)%i_of_pair,   &
                   sr_pairs(i)%j_of_pair )
         sr_pairs(i)%isAllocated = .false.
      endif
   enddo
!
   Initialized = .false.
!
   end subroutine endSample
!  ===================================================================
!
!  *******************************************************************
!              
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genLatticePoints(nbasis,na,nb,nc,small_box,bv,centered)
!  ===================================================================
   implicit none
!                          
   integer (kind=IntKind), intent(in) :: nbasis,na,nb,nc
!
   real (kind=RealKind), intent(in) :: small_box(3,3)
   real (kind=RealKind), intent(in) :: bv(3,nbasis)
!
   logical, intent(in) :: centered
!
   integer (kind=IntKind) :: i, j, k, n, isub
   integer (kind=IntKind) :: imax, imin
   integer (kind=IntKind) :: jmax, jmin
   integer (kind=IntKind) :: kmax, kmin
!        
   real (kind=RealKind) :: x0, y0, z0
!  
   if (centered) then
      imax = int(na/2)
      imin = imax-na+1
      jmax = int(nb/2)
      jmin = jmax-nb+1
      kmax = int(nc/2)
      kmin = kmax-nc+1
   else            
      imax = na-1
      imin = 0 
      jmax = nb-1 
      jmin = 0              
      kmax = nc-1           
      kmin = 0    
   endif       
!           
   n = 0       
   do k = kmin, kmax               
      do j = jmin, jmax
         do i = imin, imax
            x0 = small_box(1,1)*i + small_box(1,2)*j + small_box(1,3)*k
            y0 = small_box(2,1)*i + small_box(2,2)*j + small_box(2,3)*k
            z0 = small_box(3,1)*i + small_box(3,2)*j + small_box(3,3)*k
            do isub = 1, nbasis
               n = n + 1
               atom_position_x(n) = bv(1,isub) + x0
               atom_position_y(n) = bv(2,isub) + y0
               atom_position_z(n) = bv(3,isub) + z0
            enddo
         enddo
      enddo     
   enddo 
!           
   end subroutine genLatticePoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine placeAtomsOrderly(nbasis,Basis)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nbasis
   integer (kind=IntKind) :: i, j, k, na, isub
!
   character (len=*), intent(in) :: Basis(nbasis)
   character (len=2) :: atn(150)      ! Number of atom types < 150
!
   integer (kind=IntKind) :: nt(150)  ! Number of atom types < 150
   integer (kind=IntKind), allocatable :: ati(:)
!
   logical :: same_atom
!
   allocate( ati(nbasis) )
!
   nt(:) = 0
   atn(1) = Basis(1)
   ati(1) = 1
   nt(1) = 1
   NumTypes = 1
   do i = 2, nbasis
      same_atom = .false.
      do j = 1, NumTypes
         k = j
         if(Basis(i) == atn(j)) then
            same_atom = .true.
            exit
         endif
      enddo
      if (.not.same_atom) then
         NumTypes = NumTypes + 1
         atn(NumTypes) = Basis(i)
         ati(i) = NumTypes
         nt(NumTypes) = 1
      else
         ati(i) = k
         nt(k) = nt(k) + 1
      endif
   enddo
!
   allocate( TypeName(NumTypes), NumAtomsOfType(NumTypes) )
   do i = 1, NumTypes
      TypeName(i) = atn(i)
      NumAtomsOfType(i) = nt(i)*NumRepeatsA*NumRepeatsB*NumRepeatsC
   enddo
!
   AtomType(1) = 1
!  NumAtomsOfType(1:NumTypes) = 0
   do na = 1, NumAtoms, nbasis
      do isub = 1, nbasis
         AtomName(na+isub-1) = Basis(isub)
         AtomType(na+isub-1) = ati(isub)
      enddo
   enddo
!
   deallocate( ati )
!
   end subroutine placeAtomsOrderly
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine placeAtomsRandomly(nt,atn,content,seed_in)
!  ===================================================================
   use MathParamModule, only : ZERO, ONE, TEN2m6
   implicit none
!
   integer (kind=IntKind), intent(in) :: nt
   integer (kind=IntKind), intent(in), optional :: seed_in
!
   character (len=*), intent(in) :: atn(nt)
!
   real (kind=RealKind), intent(inout) :: content(nt)
!
   integer (kind=IntKind) :: i, na, isub, ib, jb, ka, n
!
   real (kind=RealKind) :: r, ran2
   real (kind=RealKind), allocatable :: bins(:)
!
   if (present(seed_in)) then
      iseed = seed_in
   endif
!
   NumTypes = nt
   allocate( TypeName(NumTypes), NumAtomsOfType(NumTypes) )
   do i = 1, NumTypes
      TypeName(i) = atn(i)
      NumAtomsOfType(i) = 0
   enddo
!
   allocate( bins(NumTypes+1) )
   bins(1) = ZERO
   do i = 2, NumTypes+1
      bins(i) = bins(i-1) + content(i-1)
   enddo
   if (abs(bins(NumTypes+1)-ONE) > TEN2m6) then
      call ErrorHandler('placeAtoms','The summation of the content <> 1',bins(NumTypes+1))
   endif
!
   r = ran2(iseed)
   do na = 1, NumAtoms
      r = ran2(iseed)
      LOOP_i: do i = 1, NumTypes
         if (bins(i) <= r .and. r <= bins(i+1)) then
            AtomName(na) = atn(i)
            NumAtomsOfType(i) = NumAtomsOfType(i) + 1
            AtomType(na) = i
            exit LOOP_i
         endif
      enddo LOOP_i
   enddo

!
   do ib = 1, NumTypes-1
      n = NumAtomsOfType(ib) - nint(NumAtoms*content(ib))
!  write(6,'(11x,a,t30,a,5i4)')'lillll111','=',ib,n, NumAtomsOfType(ib) ,floor(NumAtoms*content(ib)),nint(NumAtoms*content(ib))
      if (n /= 0) then
         jb = ib+1
         LOOP_dowhile: do while (jb < NumTypes)
            if (NumAtomsOfType(jb) /= nint(NumAtoms*content(jb))) then
               exit LOOP_dowhile
            else
               jb = jb + 1
            endif
         enddo LOOP_dowhile
!  write(6,'(11x,a,t30,a,5i4)')'2222222','=',ib,jb,n
         if (n < 0) then
            na = -n
         else
            na = n
         endif
         do i = 1, na
            if (n < 0) then
               ka = -1
               do while (ka == -1)
                  r = ran2(iseed)
                  ka = floor(r*NumAtoms)+1
                  if (AtomName(ka) == atn(jb)) then
                     AtomName(ka) = atn(ib)
                     AtomType(ka) = ib
                     if (NumAtomsOfType(jb) > 1) then
                        NumAtomsOfType(ib) = NumAtomsOfType(ib) + 1
                        NumAtomsOfType(jb) = NumAtomsOfType(jb) - 1
                     endif
                  else
                     ka = -1
                  endif
               enddo
            else
               ka = -1
               do while (ka == -1)
                  r = ran2(iseed)
                  ka = floor(r*NumAtoms)+1
                  if (AtomName(ka) == atn(ib)) then
                     if (NumAtomsOfType(ib) > 1) then
                        NumAtomsOfType(ib) = NumAtomsOfType(ib) - 1
                        NumAtomsOfType(jb) = NumAtomsOfType(jb) + 1
                        AtomName(ka) = atn(jb)
                        AtomType(ka) = jb
                     endif
                  else
                     ka = -1
                  endif
               enddo
            endif
         enddo
      endif
   enddo

   write(6,'(/,11x,a,t30,a,5i8)')'Atom number','=',NumAtomsOfType(1:NumTypes)
!
   deallocate( bins)
!
   do i=1,NumTypes
      content(i) = NumAtomsOfType(i)/real(NumAtoms,kind=RealKind)
   enddo
!
   write(6,'(11x,a,t30,a,$)')'Actual content','='
   do i = 1, NumTypes
      write(6,'(f10.5,$)')content(i)
   enddo
   write(6,'(a)')' '
!
!  ===================================================================
!  calculate the the short range order parameters of the sample
!  ===================================================================
   if (NumAtoms < 10000) then
!     ----------------------------------------------------------------
      call getsrop(content)
!     ----------------------------------------------------------------
   else
      write(6,'(a)')'Since number of atoms > 10000, the SROP for the random sample is not calculated.'
   endif
!
   end subroutine placeAtomsRandomly
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine placeAtomsRandomlyWithSRO(nt,atn,content,weight,nshell,srop_in, &
                                        Tmax,Tstep,seed_in)
!  ===================================================================
   use MathParamModule, only : ZERO, ONE, TEN2m6
   implicit none
!
   integer (kind=IntKind), intent(in) :: nt
   integer (kind=IntKind), intent(in) :: nshell
   integer (kind=IntKind), intent(in), optional :: seed_in
!
   character (len=*), intent(in) :: atn(nt)
!
   real (kind=RealKind), intent(inout) :: content(nt)
   real (kind=RealKind), intent(inout) :: weight(nshell)
   real (kind=RealKind), intent(in) :: srop_in(MaxAtomTypes*(MaxAtomTypes-1)/2,MaxShells)
   real (kind=RealKind), intent(in) :: Tmax
   real (kind=RealKind), intent(in) :: Tstep
 
   integer (kind=IntKind) :: i,j

!
   if (nshell > MaxShells) then
      call ErrorHandler('placeAtoms','nshell > MaxShells',nshell,MaxShells)
   endif
!
   if (present(seed_in)) then
!     ----------------------------------------------------------------
      call placeAtomsRandomly(nt,atn,content,seed_in)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call placeAtomsRandomly(nt,atn,content)
!     ----------------------------------------------------------------
   endif
! 
!  ===================================================================
!  re-arrange atom positions according to the desired short range
!  order parameters taken from the input
!  -------------------------------------------------------------------
   call sroplace(content,weight,nshell,srop_in,Tmax,Tstep)
!  -------------------------------------------------------------------
!
   end subroutine placeAtomsRandomlyWithSRO
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getsrop(content)
!  ===================================================================
!
!  *******************************************************************
!  calculate short range order parameters for various shells.
!
!  npairs_matrix[i,j] : the number of i-j pairs in the shell, where i, j are the atom types
!  npairs             : the number of pairs of any type in the shell.
!  srop               : contains a short range order parameter.
!  rpair2             : the distance square of the atom pair.
!  NumShell           : the number of short range order parameters.
!  
!  dr_shell: the atoms within (r,r+dr_shell) in shell radius are considered
!            on the same shell
!  *******************************************************************
   use MathParamModule, only : ZERO, HALF, ONE, FOUR, TEN2m8
!
   implicit   none
!
   integer (kind=IntKind), parameter :: max_shifts = 1
!
   integer (kind=IntKind) :: ia,ja
   integer (kind=IntKind) :: jp
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: id
   integer (kind=IntKind) :: n1
   integer (kind=IntKind) :: n2
   integer (kind=IntKind) :: n3
   integer (kind=IntKind) :: nmax
!
   real (kind=RealKind), intent(in) :: content(NumTypes)
   real (kind=RealKind) :: atom_position_sq, atom_position_r
   real (kind=RealKind) :: rij2
   real (kind=RealKind) :: r2max
   real (kind=RealKind) :: p1
   real (kind=RealKind) :: p2
   real (kind=RealKind) :: p3
   real (kind=RealKind) :: pmin
   real (kind=RealKind) :: pmax
   real (kind=RealKind) :: qmin
   real (kind=RealKind) :: xshift((2*max_shifts+1)**3)
   real (kind=RealKind) :: yshift((2*max_shifts+1)**3)
   real (kind=RealKind) :: zshift((2*max_shifts+1)**3)
   real (kind=RealKind) :: r2s
   real (kind=RealKind) :: r2t
   real (kind=RealKind), parameter :: rtol = TEN2m8
!
   logical :: same_atom
!
!  ===================================================================
!  determine the smallest (pmin) and the largest (pmax) distance squares
!  from the origin to the 7 other corners of the unit cell.
!  ===================================================================
   p1=large_cell(1,1)**2+large_cell(2,1)**2+large_cell(3,1)**2
   p2=large_cell(1,2)**2+large_cell(2,2)**2+large_cell(3,2)**2
   p3=large_cell(1,3)**2+large_cell(2,3)**2+large_cell(3,3)**2
   pmin=min(p1,p2,p3)
   pmax=max(p1,p2,p3)
   p1= (large_cell(1,1)+large_cell(1,2))**2+   &
       (large_cell(2,1)+large_cell(2,2))**2+   &
       (large_cell(3,1)+large_cell(3,2))**2
   pmin=min(pmin,p1)
   pmax=max(pmax,p1)
   p1= (large_cell(1,2)+large_cell(1,3))**2+   &
       (large_cell(2,2)+large_cell(2,3))**2+   &
       (large_cell(3,2)+large_cell(3,3))**2
   pmin=min(pmin,p1)
   pmax=max(pmax,p1)
   p1= (large_cell(1,3)+large_cell(1,1))**2+   &
       (large_cell(2,3)+large_cell(2,1))**2+   &
       (large_cell(3,3)+large_cell(3,1))**2
   pmin=min(pmin,p1)
   pmax=max(pmax,p1)
   p1= (large_cell(1,1)+large_cell(1,2)+large_cell(1,3))**2+   &
       (large_cell(2,1)+large_cell(2,2)+large_cell(2,3))**2+   &
       (large_cell(3,1)+large_cell(3,2)+large_cell(3,3))**2
   pmin=min(pmin,p1)
   pmax=max(pmax,p1)
!
!  ===================================================================
!  Determine qmin: the square of the distance within which the atom pairs
!  are counted for the SRO parameters.
!  Due to the periodic boundary condition effect, we will not counts
!  the pairs whoe separation is larger than the square root of qmin.
!  ===================================================================
   qmin=pmin/four-rtol
!
!  ===================================================================
!  determine the number of fuzzy shells (NumShell) and the shell distance
!  square (rpair2)
!  ===================================================================
   NumShell=0
   do ja = 1, NumAtoms-1
      LOOP_ns: do ia=ja+1,NumAtoms
         atom_position_sq = (atom_position_x(ia)-atom_position_x(ja))**2 + &
                            (atom_position_y(ia)-atom_position_y(ja))**2 + &
                            (atom_position_z(ia)-atom_position_z(ja))**2
         atom_position_r = sqrt(atom_position_sq)
         if(atom_position_sq <= qmin) then
            r2s = atom_position_sq
!           ==========================================================
!           setup an ordered rpair2 list:
!                  rpair2(1) < rpair2(2) < ... < rpair2(NumShell) 
!           ==========================================================
            LOOP_i: do i=1,NumShell      
               if(abs(sqrt(sr_pairs(i)%rpair2)-atom_position_r) <= dr_shell) then
                   cycle LOOP_ns
               else if(sr_pairs(i)%rpair2 > atom_position_sq) then
                   do j=i,NumShell
                      r2t = sr_pairs(j)%rpair2
                      sr_pairs(j)%rpair2=r2s
                      r2s=r2t
                   enddo
                   exit LOOP_i
               endif
            enddo LOOP_i
            if (NumShell < MaxShells) then
               NumShell=NumShell+1
               sr_pairs(NumShell)%rpair2=r2s
            endif
         endif
      enddo LOOP_ns
   enddo
!
   if (NumShell < 1) then
      write(6,'(/,5x,a)') "Given the small size of the unit cell sample, the SROP is not calculated ..."
      return
   endif
!
   r2max = ZERO
   do i=1,NumShell      
      r2max = max(r2max,sr_pairs(NumShell)%rpair2)
   enddo
!
   nmax=0
   do n3=-max_shifts,max_shifts
      do n2=-max_shifts,max_shifts
         do n1=-max_shifts,max_shifts
            nmax=nmax+1
            xshift(nmax)=n1*large_cell(1,1)+n2*large_cell(1,2)+n3*large_cell(1,3)
            yshift(nmax)=n1*large_cell(2,1)+n2*large_cell(2,2)+n3*large_cell(2,3)
            zshift(nmax)=n1*large_cell(3,1)+n2*large_cell(3,2)+n3*large_cell(3,3)
         enddo
      enddo
   enddo
!
   do jp=1,NumShell
      sr_pairs(jp)%npairs=0
   enddo
   do j=1,NumAtoms-1
      do i=j+1,NumAtoms
         LOOP_k: do k=1,nmax
            rij2= (xshift(k)+atom_position_x(i)-atom_position_x(j))**2+    &
                  (yshift(k)+atom_position_y(i)-atom_position_y(j))**2+    &
                  (zshift(k)+atom_position_z(i)-atom_position_z(j))**2
            if (rij2 > rtol .and. rij2 < r2max+rtol) then
               do jp=1,NumShell
                  if (abs(sqrt(rij2)-sqrt(sr_pairs(jp)%rpair2)) <= dr_shell) then
                     sr_pairs(jp)%npairs=sr_pairs(jp)%npairs+1
                     cycle LOOP_k
                  endif
               enddo
               write(6,'(/,a,3f10.5,1f12.5)')'           Rij,rij2 =',      &
                     xshift(k)+atom_position_x(i)-atom_position_x(j),      &
                     yshift(k)+atom_position_y(i)-atom_position_y(j),      &
                     zshift(k)+atom_position_z(i)-atom_position_z(j),rij2
!              -------------------------------------------------------
               call ErrorHandler('getsrop','The Rij is not in the table.')
!              -------------------------------------------------------
            endif
         enddo LOOP_k
      enddo
   enddo
!
   do jp=1,NumShell
      if (sr_pairs(jp)%isAllocated) then
         deallocate( sr_pairs(jp)%i_of_pair, sr_pairs(jp)%j_of_pair,       &
                     sr_pairs(jp)%srop, sr_pairs(jp)%npairs_matrix )
         nullify( sr_pairs(jp)%i_of_pair, sr_pairs(jp)%j_of_pair,          &
                  sr_pairs(jp)%srop, sr_pairs(jp)%npairs_matrix )
      endif
      allocate( sr_pairs(jp)%i_of_pair(sr_pairs(jp)%npairs) )
      allocate( sr_pairs(jp)%j_of_pair(sr_pairs(jp)%npairs) )
      allocate( sr_pairs(jp)%srop(NumTypes*(NumTypes-1)/2,MaxShells) )
      allocate( sr_pairs(jp)%npairs_matrix(NumTypes,NumTypes) )
      sr_pairs(jp)%isAllocated = .true.
   enddo
!
   do jp=1,NumShell
      id = 0
      do ja=1,NumAtoms-1
         do ia=ja+1,NumAtoms
            do k=1,nmax
               rij2= (xshift(k)+atom_position_x(ia)-atom_position_x(ja))**2+    &
                     (yshift(k)+atom_position_y(ia)-atom_position_y(ja))**2+    &
                     (zshift(k)+atom_position_z(ia)-atom_position_z(ja))**2
               if (rij2 > rtol .and. rij2 < r2max+rtol) then
                  if (abs(sqrt(rij2)-sqrt(sr_pairs(jp)%rpair2)) <= dr_shell) then
                     id = id + 1
                     sr_pairs(jp)%i_of_pair(id)=ia
                     sr_pairs(jp)%j_of_pair(id)=ja
                  endif
               endif
            enddo
         enddo
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call calSROP(content,.true.)
!  -------------------------------------------------------------------
!
   end subroutine getsrop
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSROP(content,lprint)
!  ===================================================================
   use MathParamModule, only : HALF, ONE
   implicit none
!
   real (kind=RealKind), intent(in) :: content(NumTypes)
!
   logical, intent(in), optional :: lprint
!
   integer (kind=IntKind) :: i, j, ia, ja, id, ip, isro
!
   real (kind=RealKind) :: pab
!
   do ip=1,NumShell
      sr_pairs(ip)%npairs_matrix(1:NumTypes,1:NumTypes)=0
      do id=1,sr_pairs(ip)%npairs
         i=sr_pairs(ip)%i_of_pair(id)
         j=sr_pairs(ip)%j_of_pair(id)
         ia = AtomType(i)
         ja = AtomType(j)
         if (ia == ja) then
            sr_pairs(ip)%npairs_matrix(ia,ia)=sr_pairs(ip)%npairs_matrix(ia,ia)+1
         else
            sr_pairs(ip)%npairs_matrix(ia,ja)=sr_pairs(ip)%npairs_matrix(ia,ja)+1
            sr_pairs(ip)%npairs_matrix(ja,ia)=sr_pairs(ip)%npairs_matrix(ja,ia)+1
         endif
      enddo
      if (NumTypes > 1) then
         isro=1
         do j = 1, NumTypes-1
            do i = j+1, NumTypes
              pab=sr_pairs(ip)%npairs_matrix(i,j)/real(sr_pairs(ip)%npairs,kind=RealKind)
              sr_pairs(ip)%srop(isro,ip)=one-half*pab/(content(i)*content(j)) 
              isro=isro+1
            enddo
         enddo
      endif
   enddo
   if ( present(lprint) ) then
      if (lprint) then
         write(6,'(/,a,t70,a,2i8)')' calSROP:: the (max/required) no. of shells and SRO for each shell', &
               '=',NumShell, NumTypes*(NumTypes-1)/2
         write(6,'(/,10x,60(''=''))')
         if (NumTypes == 1) then
            write(6,'(10x,a)')'ishell   npairs    pair    naa'
         else
            write(6,'(10x,a)')'ishell   npairs    pair      naa     nab     nbb      srop'
         endif
         write(6,'(10x,60(''-''))')
         do ip=1,NumShell
           if (NumTypes == 1) then
             write(6,'(11x,1i5,1i8,$)')ip,sr_pairs(ip)%npairs
             write(6,'(i10)')sr_pairs(ip)%npairs_matrix(1,1)
           else
             isro=1
             do j = 1, NumTypes-1
               do i = j+1, NumTypes
                 write(6,'(11x,1i5,1i8,$)')ip,sr_pairs(ip)%npairs
                 write(6,'(A5,3x,A,3i8,3x,1f9.5)')                                                              &
                     TypeName(i), TypeName(j),sr_pairs(ip)%npairs_matrix(i,i),sr_pairs(ip)%npairs_matrix(i,j),  &
                     sr_pairs(ip)%npairs_matrix(j,j),sr_pairs(ip)%srop(isro,ip)
                 isro=isro+1
                 if (j < NumTypes-1 .or. i < NumTypes) then
                   write(6,'(24x,$)')
                 endif
               enddo
             enddo
           endif
        enddo
!       if (NumTypes > 2) then
!         write(6,'(a)')' '
!       endif
        write(6,'(10x,60(''=''))')
!        do i=1,NumAtoms
!          write(6,'(a5,3f10.5)')AtomName(i),AtomPosition(i,1:3)*0.529177208
!        enddo
      endif
   endif
!
   end subroutine calSROP
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sroplace(content,weight,nshell_in,srop_in,Tmax,Tstep)
!  ===================================================================
   use MathParamModule, only : ZERO, HALF, ONE, TWO
   use PhysParamModule, only : Boltzmann
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: nshell_in
!
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: ip,id
   integer (kind=IntKind) :: i, myloopflag,myi,myj,flag,mysro
   integer (kind=IntKind) :: min_step
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: it
   integer (kind=IntKind) :: ia, ja, isro
   integer (kind=IntKind) :: idum = -1
   integer (kind=IntKind) :: num_temp_steps                ! No. of steps to reduce T before -> 0
   integer (kind=IntKind) :: num_picks, NumSwaps
   integer (kind=IntKind) :: pos_i
   integer (kind=IntKind) :: pos_j
   integer (kind=IntKind) :: change = 0
   integer (kind=IntKind), parameter :: max_picks = 50000 ! Maximum no. of picks at each step
!
   real (kind=RealKind), intent(in) :: content(NumTypes)
   real (kind=RealKind), intent(inout) :: weight(nshell_in)
   real (kind=RealKind), intent(in) :: srop_in(MaxAtomTypes*(MaxAtomTypes-1)/2,MaxShells)
   real (kind=RealKind), intent(in) :: Tmax                ! The initial temperature (T)
   real (kind=RealKind), intent(in) :: Tstep               ! The temperature descreasing step
!
   real (kind=RealKind) :: ran2
!   real (kind=RealKind), allocatable :: weight(:)
   real (kind=RealKind) :: weight_sum
   real (kind=RealKind) :: aver_diff_single, aver_diff_total, aver_diff_total2
   real (kind=RealKind) :: aver_diff_new
   real (kind=RealKind) :: fac
   real (kind=RealKind) :: tfac, r
   real (kind=RealKind) :: temperature
   real (kind=RealKind) :: diff_check
   real (kind=RealKind) :: min_diff
   real (kind=RealKind) :: pab=0
   real (kind=RealKind), parameter :: tol = 0.1d0
   real (kind=RealKind), parameter :: perror = 0.25  ! on average, no more than perror of total pairs of A-B 
                                                     ! atoms are allowed off from the desired value.
                                                     ! To speed up, you may increase this value.
   real (kind=RealKind), parameter :: kbf =1.0d3     ! This is the factor in the front of Boltzmann
!
   integer (kind=IntKind), allocatable :: nab_in(:,:)
   real (kind=RealKind), allocatable :: sro_best(:,:)
   real (kind=RealKind), allocatable :: weight_norm(:)

!
   if(nshell_in > NumShell) then
      write(6,'(a,i3)')'The max shell number is only:  ',NumShell
   endif

   NumShell=min(NumShell,nshell_in)    ! number of shells considered

   fac=ran2(idum)
!
   allocate( nab_in(NumTypes*(NumTypes-1)/2,NumShell))
   allocate( sro_best(NumTypes*(NumTypes-1)/2,NumShell))
   allocate( weight_norm(NumShell))
!
!  ===================================================================
!  determine weight for each shell....................................
!  ===================================================================
!   weight(1)=one
!   do ip=2,NumShell   ! let weight decreases linearly so that weight(1)=1.0
!      weight(ip)=(half*ip+half-NumShell)/(one-NumShell) ! while weight(NumShell)=0.5
!   enddo
!   write(6,'(2x,a,5f15.5)') 'mmm weight is:',weight(1:NumShell)
   weight_sum=zero
    do ip=1,NumShell
      weight_sum=weight_sum+weight(ip)   ! weight_sum is the normalization
   enddo                                 ! factor.
!   write(6,'(2x,a,5f15.5)') 'weight sum is:',weight_sum

   do ip=1,NumShell
      weight_norm(ip)=weight(ip)/weight_sum
   enddo                   
   write(6,'(a,8f12.6)')'      The weight of each shell is:  ',weight_norm(1:NumShell)

!  ===================================================================
!  get nab_in and initial aver_diff...................................
!  ===================================================================

   aver_diff_total=zero
   do ip=1,NumShell
     isro=1
     do j = 1, NumTypes-1
       do i = j+1, NumTypes
          aver_diff_total= aver_diff_total +                          &
             weight_norm(ip)*2*content(i)*content(j)*abs(sr_pairs(ip)%srop(isro,ip)- &
                                                         srop_in(isro,ip))*sr_pairs(ip)%npairs
          nab_in(isro,ip)=int(2*content(i)*content(j)*(one-srop_in(isro,ip))*sr_pairs(ip)%npairs+half)

          pab=nab_in(isro,ip)/real(sr_pairs(ip)%npairs,kind=RealKind)
          sro_best(isro,ip)=one-half*pab/(content(i)*content(j))
          isro=isro+1
       enddo
     enddo
   enddo

!   write(6,'(/)')
   write(6,'('' SROPLACE:: nab corresponding to the input SRO parameters:'')')
   write(6,'(''            input sro_in, nab_in, sro_best          ='')')
   do ip=1,NumShell
     do isro=1, NumTypes*(NumTypes-1)/2
       write(6,'(9x,f10.5,i10,f12.5,$)')srop_in(isro,ip),nab_in(isro,ip),sro_best(isro,ip)
     enddo
     write(6,'(a)')' '
   enddo


   num_temp_steps = ceiling(Tmax/Tstep)
   num_temp_steps=max(int(log(aver_diff_total/perror)/log(2.0d0))+1,num_temp_steps)

   write(6,'(''            initial temperature ='',1f12.5)')Tmax
   write(6,'(''            initial aver_diff   ='',1f12.5)')aver_diff_total
   write(6,'(''            no. of steps to 0 K ='',1i5)')num_temp_steps+1
   write(6,'(/,7x,70(''=''))')
   write(6,'(10x,a,$)')'T       aver_diff    num_swaps'
   do ip = 1, NumShell
     do isro=1, NumTypes*(NumTypes-1)/2
      write(6,'(4x,a,i1,i2,a,$)')'srop(',isro,ip,')'
     enddo
   enddo
   write(6,'(a)')' '
   write(6,'(7x,70(''-''))')

!  diff_check=0.500d0
   diff_check=ONE


      
   do it = 0, num_temp_steps
      temperature=Tmax*(ONE-it/real(num_temp_steps,kind=RealKind))
      tfac=kbf*Boltzmann*temperature  
      num_picks=0
      NumSwaps = 0
      do while( aver_diff_total > 0 .and. (NumSwaps < NumAtoms .and. num_picks < max_picks) )
!        =============================================================
!        determine the position of the first randomly-selected site.
!        =============================================================
         pos_i=floor(NumAtoms*ran2(idum))+1
!        =============================================================
!        determine the position of the second randomly-selected site.
!        =============================================================
         pos_j=floor(NumAtoms*ran2(idum))+1
         ia = AtomType(pos_i)
         ja = AtomType(pos_j)
         if (ja /= ia) then  ! if the two atoms are different
!         if ((ia ==myi .and. ja == myj ) .or. (ia==myj .and. ja==myi)) then  ! if the two atoms are different
            AtomType(pos_i) = ja
            AtomType(pos_j) = ia
!           ----------------------------------------------------------
!           here to calculate srop in regardless exchange or not

            call calSROP(content)
!           ----------------------------------------------------------
!           ==========================================================
!           determine if atoms on pos_i and pos_j need to be exchanged
!           ==========================================================
            num_picks=num_picks+1
            aver_diff_new=zero

            do ip=1,NumShell
              isro=1
              do j = 1, NumTypes-1
                do i = j+1, NumTypes
                   aver_diff_new = aver_diff_new + &
                                   weight_norm(ip)*2*content(i)*content(j)*abs(sr_pairs(ip)%srop(isro,ip)- &
                                                                               srop_in(isro,ip))*sr_pairs(ip)%npairs
                   isro=isro+1
                enddo
              enddo
            enddo

!           ==========================================================
!           if exp(-diff/tfac) > ran2 accept the exchange
!           ==========================================================
            r = tfac*log(ran2(idum))
!            write(6,'(6x,f8.2)')r
!            write(6,'(6x,f8.2)')aver_diff-aver_diff_new
            if (aver_diff_total-aver_diff_new > r) then  ! accept the pair swap
               NumSwaps = NumSwaps + 1
               AtomName(pos_i) = TypeName(ja)
               AtomName(pos_j) = TypeName(ia)
               aver_diff_total=aver_diff_new
            else                                   ! undo the pair swap..
               AtomType(pos_i) = ia
               AtomType(pos_j) = ja
               call calSROP(content)  ! no exchange, no change of srop
            endif
         endif                         
         
!         if (aver_diff_total==0) then
!           exit
!         endif

      enddo
      write(6,'(6x,f8.2,2x,f11.5,5x,1i6,$)')temperature,aver_diff_total,NumSwaps
      do ip = 1, NumShell
         do isro=1, NumTypes*(NumTypes-1)/2
           write(6,'(5x,f8.5,$)')sr_pairs(ip)%srop(isro,ip)
         enddo
      enddo
      write(6,'(a)')' '
   enddo
   write(6,'(7x,70(''=''))')
   write(6,'(/,11x,''aver_diff_single ='',f15.8)')aver_diff_total

   call calSROP(content,.true.)

   aver_diff_total=zero
   aver_diff_total2=zero
   do ip=1,NumShell
     isro=1
     do j = 1, NumTypes-1
       do i = j+1, NumTypes
          aver_diff_total= aver_diff_total + &
                           weight_norm(ip)*2*content(i)*content(j)*abs(sr_pairs(ip)%srop(isro,ip)- &
                                                                       srop_in(isro,ip))*sr_pairs(ip)%npairs
          aver_diff_total2= aver_diff_total2 + &
                           weight_norm(ip)*2*content(i)*content(j)*abs(sr_pairs(ip)%srop(isro,ip)- &
                                                                       sro_best(isro,ip))*sr_pairs(ip)%npairs
          isro=isro+1
       enddo
     enddo
   enddo
   write(6,'(/,11x,''aver_diff_total compared with srop_in and srop_best ='',2f15.8)')aver_diff_total,aver_diff_total2
   write(6,'(7x,70(''=''))')
      
!
!  =====================================================================
!  recalculate the SRO parameters....................................
!  ---------------------------------------------------------------------
!   call calSROP(content,.true.)
!  ---------------------------------------------------------------------
!   deallocate( nab_in, weight, sro_best)
   deallocate( nab_in, sro_best)
   deallocate( weight_norm)
!
   end subroutine sroplace
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomPosition0() result(ap)
!  =====================================================================
   implicit none
!
   real (kind=RealKind), pointer :: ap(:,:)
!
   ap => AtomPosition(1:NumAtoms,1:3)
!
   end function getAtomPosition0
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomPosition1(i) result(ap)
!  =====================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
!
   real (kind=RealKind), pointer :: ap(:)
!
   if (i < 1 .or. i > 3) then
      call ErrorHandler('getAtomPosition','Invalid component index',i)
   endif
!
   ap => AtomPosition(1:NumAtoms,i)
!
   end function getAtomPosition1
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomName0() result(an)
!  =====================================================================
   implicit none
!
   character (len=2), pointer :: an(:)
!
   an => AtomName(1:NumAtoms)
!
   end function getAtomName0
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomName1(i,SpeciesIndex) result(an)
!  =====================================================================
   implicit none
!
   character (len=2) :: an
!
   integer (kind=IntKind), intent(in) :: i
!
   logical, intent(in), optional :: SpeciesIndex
!
   if (present(SpeciesIndex)) then
      if (SpeciesIndex) then
         if (i < 1 .or. i > NumTypes) then
            call ErrorHandler('getAtomName','Species type index is out of range',i)
         endif
         an = TypeName(i)
         return
      endif
   endif
!
   if (i < 1 .or. i > NumAtoms) then
      call ErrorHandler('getAtomName','Atom index is out of range',i)
   endif
!
   an = AtomName(i)
!
   end function getAtomName1
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getUnitCell() result(uc)
!  =====================================================================
   implicit none
!
   real (kind=RealKind), pointer :: uc(:,:)
!
   uc => large_cell(1:3,1:3)
!
   end function getUnitCell
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumAtoms(Species) result(n)
!  =====================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: Species
   integer (kind=IntKind) :: n
!
   if (present(Species)) then
      if (Species < 1 .or. Species > NumTypes) then
         call ErrorHandler('getNumAtoms','Species type index is out of range',Species)
      else
         n = NumAtomsOfType(Species)
      endif
   else
      n = NumAtoms
   endif
!
   end function getNumAtoms
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumSpecies() result(n)
!  =====================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumTypes
!
   end function getNumSpecies
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSampleAtom(an) result(t)
!  =====================================================================
   implicit none
!
   character (len=2), intent(in) :: an
!
   integer (kind=IntKind) :: i
   logical :: t
!
   t = .false.
   do i = 1, NumTypes
     if (an == TypeName(i)) then
        t = .true.
        exit
     endif
   enddo
!
   end function isSampleAtom
!  =====================================================================
!
!  *********************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine substituteAtoms(nimp,imp,impc,repl)
!  =====================================================================
   implicit none
!
!  *********************************************************************
!  Purpose: Given nimp impurity types, perform:
!           do i = 1, nimp
!              replacing host atom repl(i) by 
!              impurity imp(i) of content impc(i)
!           enddo
!  *********************************************************************
!
   logical :: replaced
!
   integer (kind=IntKind), intent(in) :: nimp
   integer (kind=IntKind) :: i, j, k, m, n, add_types
   integer (kind=IntKind) :: np(nimp), nh(NumTypes), impurity_type(nimp)
   integer (kind=IntKind) :: num_host_species(NumTypes)
!
   character (len=2), intent(in) :: imp(nimp), repl(nimp)
   character (len=2) :: host_species(NumTypes)
!
   real (kind=RealKind), intent(in) :: impc(nimp)
   real (kind=RealKind) :: r, ran2
!
   add_types = nimp
   do j = 1, nimp
      if (repl(j) == imp(j)) then
!        ---------------------------------------------------------------
         call ErrorHandler('substituteAtoms','The host atom is the same as the impurity', &
                           imp(j))
!        ---------------------------------------------------------------
      endif
      impurity_type(j) = 0
      do i = 1, NumTypes
         if (imp(j) == TypeName(i)) then
            add_types = add_types - 1
            impurity_type(j) = i
         endif
      enddo
   enddo
!
   do j = 1, nimp
      np(j) = 0
      nh(j) = 0
      LOOP_i: do i = 1, NumTypes
         if (repl(j) == TypeName(i)) then
            np(j) = max(int(impc(j)*NumAtomsOfType(i)),1)
            nh(j) = i
            exit LOOP_i
         endif
      enddo LOOP_i
   enddo
!
   host_species = TypeName
   num_host_species = NumAtomsOfType
!
   deallocate(TypeName, NumAtomsOfType)
   allocate(TypeName(NumTypes+add_types))
   allocate(NumAtomsOfType(NumTypes+add_types))
!
   do i = 1, NumTypes
      TypeName(i) = host_species(i)
      NumAtomsOfType(i) = num_host_species(i)
   enddo
   i = NumTypes
   do j = 1, nimp
      if (impurity_type(j) == 0) then
         i = i + 1
         TypeName(i) = imp(j)
         impurity_type(j) = i
         NumAtomsOfType(i) = 0
      endif
   enddo
!
   do j = 1, nimp
      do i = 1, np(j)
         r = ran2(iseed)
      enddo
   enddo
!
   do j = 1, nimp
      m = nh(j) ! The type index of the host atom to be replaced
      n = impurity_type(j) ! the type index of the impurity atom
      do i = 1, np(j)
         LOOP_do: do
            r = ran2(iseed)
            k = floor(r*NumAtoms)+1
            if (AtomName(k) == repl(j)) then
               AtomName(k) = imp(j)
               AtomType(k) = n
               NumAtomsOfType(m) = NumAtomsOfType(m) - 1
               NumAtomsOfType(n) = NumAtomsOfType(n) + 1
               exit LOOP_do
            endif
         enddo LOOP_do
      enddo
   enddo
!
   NumTypes = NumTypes + add_types
!
   end subroutine substituteAtoms
!  =====================================================================
end module SampleModule
