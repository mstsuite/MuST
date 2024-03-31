!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readPOSCAR(fpath,NumAtomsIn,NumAtomsOut)
!  ===================================================================
!
!  Read POSCAR file used for VASP
!
!  *******************************************************************
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ZERO, ONE, TEN2m8
   use PhysParamModule, only : Angstrom2Bohr
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
   use ChemElementModule, only : MaxLenOfAtomName, getZtot, isAtomName, getName
!
   use StringModule, only : initString, setString, endString,         &
                            getNumTokens, readToken
!
   implicit none
!
   character (len=*), intent(in) :: fpath
!
   integer (kind=IntKind), parameter :: text_len = 80
!
   character (len=text_len) :: text, tmp
   character (len=MaxLenOfAtomName) :: atom
   character (len=20), allocatable :: atom_name_of_type(:)
   character (len=1), pointer :: AtomNameCharacter(:)
!
   logical :: FileExist
   logical :: cart = .false.
   logical :: seldyn = .false.
!
   integer (kind=IntKind), intent(in), optional :: NumAtomsIn
   integer (kind=IntKind), intent(out), optional :: NumAtomsOut
   integer (kind=IntKind) :: num_atoms
!
   integer (kind=IntKind), parameter :: IOPE = 0
   integer (kind=IntKind), parameter :: funit = 103
   integer (kind=IntKind) :: ios, status, n, i, j, k
   integer (kind=IntKind) :: line_number, num_tokens, pre_pos_line
   integer (kind=IntKind) :: num_atom_types
   integer (kind=IntKind), allocatable :: num_atoms_of_type(:)
   integer (kind=IntKind), pointer :: AtomicNumber(:)
!
   real (kind=RealKind), pointer :: pa
   real (kind=RealKind), pointer :: AtomPosition(:,:)
   real (kind=RealKind), pointer :: BravaisLatticeVec(:,:)
   real (kind=RealKind) :: x, y, z
   real (kind=RealKind) :: fac1, fac2, fac3
   real (kind=RealKind) :: bra(3,3), vec(3), v
   real (kind=RealKind), parameter :: tol = TEN2m8
!
   interface
      function isNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isNumber
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
      function trim_string(s1,c) result(s2)
         character (len=*), intent(in) :: s1
         character (len=len(s1)) :: s2
         character (len=1), optional, intent(in) :: c
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
!  ===================================================================
!  Check the status of the POSCAR file
!  ===================================================================
   if (MyPE == IOPE) then
      inquire(file=trim(fpath)//'POSCAR',exist=FileExist)
      if (FileExist) then
         open(unit=funit,file=trim(fpath)//'POSCAR',form='formatted',status='old', &
              iostat=ios,action='read')
         if (ios > 0) then
            call ErrorHandler('realPOSCAR','iostat > 0',ios)
         endif
      else
         call ErrorHandler('readPOSCAR','POSCAR does not exist in the current directory')
      endif
   endif
!
!  ===================================================================
!  Determine the number of atoms in the POSCAR data
!  ===================================================================
   num_atoms = 0
   if (MyPE == IOPE) then
!     ----------------------------------------------------------------
      call initString(text_len)
!     ----------------------------------------------------------------
      line_number = 0
      LOOP_dowhile: do while (num_atoms == 0)
         read(funit,'(a)',iostat=status) text
         if (status /= 0) then
            call ErrorHandler('readPOSCAR','Read POSCAR data',trim(text))
         endif
!
         tmp = adjustl(text)
         tmp = trim_string(tmp,'#')
         text = trim_string(tmp,'!')
         if (len_trim(text) == 0) then
            cycle LOOP_dowhile
         endif
!
         line_number = line_number + 1
!
         if (line_number == 1) then
!           ==========================================================
!           Comment (mandatory)
!           The first line is reserved for a free user comment, e.g. a system description. The maximum line length is 40 characters,
!           extra characters are truncated
!           ==========================================================
            cycle LOOP_dowhile
         endif
!
!
!        -------------------------------------------------------------
         call setString(text)
!        -------------------------------------------------------------
         num_tokens = getNumTokens()
!
         if (line_number == 2) then
            fac1 = ONE; fac2 = ONE; fac3 = ONE
!           ==========================================================
!           Scaling factor(s) (mandatory)
!           This line may contain one or three numbers. If one number is provided it specifies a universal lattice scaling factor s. It is
!           multiplied with the three vectors in the following section to obtain the lattice vectors of the unit cell. Also, the ion positions
!           are scaled with this factor if the "Cartesian" mode is selected (see section "Ion positions"). If the number is negative, it is
!           interpreted as the desired cell volume. Then, the scaling factor s is computed automatically to obtain the desired volume. If
!           three numbers are provided in this line they act as individual scaling factors for the x-,y- and z-Cartesian components for the
!           lattice vectors (and "Cartesian" mode ion positions). In this case all three numbers must be positive. 
!           ==========================================================
            if (num_tokens == 1) then
               read(text,*,iostat=status)fac1
               if (status /= 0) then
                  call ErrorHandler('readPOSCAR','Invalid scaling factor line',trim(text))
               else if (abs(fac1) < tol) then
                  call ErrorHandler('readPOSCAR','Invalid scaling factor line',trim(text))
               else
                  fac2 = fac1
                  fac3 = fac1
               endif
            else if (num_tokens == 3) then
               read(text,*,iostat=status)fac1, fac2, fac3
               if (status /= 0) then
                  call ErrorHandler('readPOSCAR','Invalid scaling factor line',trim(text))
               else if (fac1 < tol .or. fac2 < tol .or. fac3 < tol) then
                  call ErrorHandler('readPOSCAR','Invalid scaling factor line',trim(text))
               endif
            else
               call ErrorHandler('readPOSCAR','invalid scaling factor line',trim(text))
            endif
         else if (line_number <= 5) then
!           ==========================================================
!           Lattice (mandatory)
!           This sections contains three lines defining the lattice vectors. Each line holds the unscaled Cartesian components of
!           one lattice vector. The actual lattice vectors a1, a2, and a3 (in Angstrom ) are the product of the given numbers with 
!           the lattice scaling factor s. Set the universal scaling factor to 1 if you want to enter the lattice vectors directly 
!           and avoid any additional scaling.
!           ==========================================================
            if (num_tokens /= 3) then
               call ErrorHandler('readPOSCAR','Invalid Bravais lattice vector line',trim(text))
            endif
            read(text,*,iostat=status)bra(1:3,line_number-2)
            if (status /= 0) then
               call ErrorHandler('readPOSCAR','read Bravais lattice vector',trim(text))
            endif
         else if (line_number <= 7) then
!           ==========================================================
!           Species names (optional)
!           This optional line lists the species of the present ions. The given order should match the order of 
!           species appearing in the POTCAR file. This line is optional, if omitted the species names are taken from the POTCAR file. 
!
!           Ions per species (mandatory)
!           This mandatory line lists how many ions of each species are present. The given order should match the
!           order of species appearing in the POTCAR file. 
!           ==========================================================
            if (num_tokens < 1) then
               call ErrorHandler('readPOSCAR','Invalid line format',trim(text))
            endif
!           ==========================================================
!           check the format of the line: whether it contains species names
!           ==========================================================
            i = getTokenPosition(1,text,k)
            if (line_number == 6) then
               num_atom_types = num_tokens
               allocate(atom_name_of_type(num_atom_types))
               allocate(num_atoms_of_type(num_atom_types))
               if (isAtomName(text(i:i+k-1))) then ! The atom name is assumed to the element name with 1 or 2 letters
                                                   ! In the future, we should allow for full element name case.
                  read(text,*,iostat=status) (atom_name_of_type(n), n = 1, num_atom_types)
                  if (status /= 0) then
                     call ErrorHandler('readPOSCAR','Invalid atom name list line',trim(text))
                  endif
                  cycle LOOP_dowhile
               else
                  inquire(file=trim(fpath)//'POTCAR',exist=FileExist)
                  if (.not.FileExist) then
                     write(6,'(a)')'The last item read:'
                     write(6,'(a)')text(i:k-1)
                     call ErrorHandler('readPOSCAR','POTCAR is needed but it does not exist in the current directory')
                  endif
                  call system('grep VRHFIN '//trim(fpath)//'POTCAR'//' |cut -d"=" -f2|cut -d":" -f1 > GREP_POTCAR_OUT')
                  open(unit=1010,file='GREP_POTCAR_OUT',iostat=ios,status='old',form='formatted')
                  if (ios > 0) then
                     call ErrorHandler('readPOSCAR','Open GREP_POTCAR_OUT',ios)
                  endif
                  do n = 1, num_atom_types
                     read(1010,'(a)',iostat=status)atom_name_of_type(n)
                     if (status /= 0) then
                        call ErrorHandler('readPOSCAR','Invalid atom name line',atom_name_of_type(n))
                     endif
                  enddo
                  close(unit=1010,status='delete')
               endif
            endif
            read(text,*,iostat=status) (num_atoms_of_type(n), n = 1, num_atom_types)
            if (status /= 0) then
               call ErrorHandler('readPOSCAR','Invalid atom number list line',trim(text))
            endif
            do n = 1, num_atom_types
               num_atoms = num_atoms + num_atoms_of_type(n)
            enddo
         else
            call ErrorHandler('readPOSCAR','Invalid POSCAR data format -> num_atoms = 0')
         endif
      enddo LOOP_dowhile
   endif
!  -------------------------------------------------------------------
   call bcastMessage(num_atoms,IOPE)
!  -------------------------------------------------------------------
   if (num_atoms < 1) then
      call ErrorHandler('readPOSCAR','num_atoms < 1',num_atoms)
   endif 
!
   if (present(NumAtomsIn)) then
      if (num_atoms /= NumAtomsIn) then
         call ErrorHandler('readPOSCAR','num_atoms <> NumAtomsIn',num_atoms,NumAtomsIn)
      endif
   endif
!
   if (present(NumAtomsOut)) then
      NumAtomsOut = num_atoms
   endif
!
!  -------------------------------------------------------------------
   call createDataStorage('Position Scaling Factor',RealType)
   call createDataStorage('Bravais Vector',9,RealType)
   call createDataStorage('Atomic Number',num_atoms,IntegerType)
   call createDataStorage('Atomic Alt Name',MaxLenOfAtomName*num_atoms,CharacterType)
   call createDataStorage('Atomic Position',3*num_atoms,RealType)
!  -------------------------------------------------------------------
   pa => getDataStorage('Position Scaling Factor',RealMark); pa = ONE
   BravaisLatticeVec => getDataStorage('Bravais Vector',3,3,RealMark)
   AtomicNumber => getDataStorage('Atomic Number',num_atoms,IntegerMark)
   AtomNameCharacter => getDataStorage('Atomic Alt Name',MaxLenOfAtomName*num_atoms,CharacterMark)
   AtomPosition => getDataStorage('Atomic Position',3,num_atoms,RealMark)
!  -------------------------------------------------------------------
!
   if (MyPE == IOPE) then
      if (fac1 < ZERO) then
!        =============================================================
!        If the scaling number is negative, it is interpreted as the desired cell volume. 
!        Then, the scaling factor is computed to obtain the desired volume.
!        -------------------------------------------------------------
         call vcross(bra(:,1),bra(:,2),vec)
!        -------------------------------------------------------------
         v = vec(1)*bra(1,3)+vec(2)*bra(2,3)+vec(3)*bra(3,3)
         fac1 = abs(fac1/v)
         fac2 = fac1; fac3 = fac1
      endif
!
      do i = 1, 3
         BravaisLatticeVec(1,i) = fac1*bra(1,i)*Angstrom2Bohr
         BravaisLatticeVec(2,i) = fac2*bra(2,i)*Angstrom2Bohr
         BravaisLatticeVec(3,i) = fac3*bra(3,i)*Angstrom2Bohr
      enddo
!
      j = 0
      do n = 1, num_atom_types
         atom = atom_name_of_type(n)
         if (nocaseCompare(atom,'Va')) then
            atom = 'V'
         else if (.not.isAtomName(atom)) then
            call ErrorHandler('readPOSCAR','species name is invalid',atom)
         endif
         k = getZtot(atom)
         do i = 1, num_atoms_of_type(n)
            j = j + 1
            AtomicNumber(j) = k
            call copyString2CharArray(atom,AtomNameCharacter((j-1)*MaxLenOfAtomName+1:j*MaxLenOfAtomName), &
                                      MaxLenOfAtomName)
         enddo
      enddo
   endif
!  -------------------------------------------------------------------
   call bcastMessage(BravaisLatticeVec,3,3,IOPE)
   call bcastMessage(AtomNameCharacter,MaxLenOfAtomName*num_atoms,IOPE)
   call bcastMessage(AtomicNumber,num_atoms,IOPE)
!  -------------------------------------------------------------------
!
   if (MyPE == IOPE) then
      seldyn = .false.
      cart = .false.
      pre_pos_line = line_number
      n = 0
      LOOP_do: do
         read(funit,'(a)',iostat=status) text
         if (status < 0) then
            exit LOOP_do
         else if (status > 0) then
            call ErrorHandler('readPOSCAR','Read POSCAR data',trim(text))
         endif
!
         tmp = adjustl(text)
         tmp = trim_string(tmp,'#')
         text = trim_string(tmp,'!')
         if (len_trim(text) == 0) then
            cycle LOOP_do
         endif
!
         line_number = line_number + 1
!
         if (line_number-pre_pos_line == 1) then
!           ==========================================================
!           Selective dynamics (optional)
!           If the line after the "Ions per species" section contains Selective dynamics it enables the "selective dynamics" feature
!           (actually only the first character is relevant and must be S or s). This allows to provide extra flags for each atom
!           signaling whether the respective coordinate(s) of this atom will be allowed to change during the ionic relaxation. This
!           setting is useful if only certain shells around a defect or layers near a surface should relax. See also the IBRION tag. 
!           ==========================================================
            if (nocaseCompare(text(1:1),'S')) then
               if (.not.seldyn) then
                  seldyn = .true.
                  pre_pos_line = line_number
                  cycle LOOP_do
               else
                  call ErrorHandler('readPOSCAR','Invalid line',trim(text))
               endif
            endif
!
!           ==========================================================
!           Ion position mode (mandatory)
!           This line selects one of the two possible modes how the coordinates x1, x2, and x3 given in the following lines are 
!           interpreted:
!           * "Direct" means the positions are provided in direct (fractional) coordinates:
!             R = x1*a1 + x2*a2 + x3*a3
!             where R is the ion position vector
!           * "Cartesian" specifies that positions are provided in a Cartesian coordinate system. However, the actual ion positions
!             are also multiplied with the universal scaling factor, i.e.
!             R = s*(x1, x2, x3)
!           Actually, only the first character on the line is significant and the only key characters recognized are C, c, K or k for
!           switching to the "Cartesian" mode. Everything else will be interpreted as "Direct" mode.
!           ==========================================================
            if (nocaseCompare(text(1:1),'C') .or. nocaseCompare(text(1:1),'K')) then
               cart = .true.
            else
               cart = .false.
            endif 
         else
!           ==========================================================
!           Ion positions (mandatory)
!           Here, the ion positions (in Angstrom) are listed. 
!           The total number of lines with positions must match the total number of ions given in the "Ions per species" section. The ion
!           species are also derived from there, e.g. if the "Ions per species" section lists 5 8, then there must be five ion position
!           lines for the first species, followed by eight ions of the second species. If your are not sure whether you have a correct 
!           input please check the OUTCAR file, which contains both the final Cartesian components of the vector R and the positions 
!           in direct (fractional) coordinates.
!           If the selective dynamics feature is enabled on each coordinate triplet is followed by three additional logical flags, i.e.
!           each is either T or F for true and false, respectively. This determines whether to allow changes of the coordinates or not.
!           If the line selective dynamics is removed from the POSCAR file this flag will be ignored (and internally set to T).
!           Mind: The flags refer to the positions of the ions in direct coordinates, no matter whether the positions are entered in 
!                 "Cartesian" or "Direct" coordinate modes.
!           For example, consider the following ion specification:
!           ...
!           Selective dynamics
!           Cartesian
!           0.00 0.00 0.00 T F T
!           1.27 0.98 0.32 F T F
!           ...
!           Here, the first atom is allowed to move into the direction of the first and third direct lattice vector. The second atom may only
!           move in the second lattice vector direction.
!           If no initial velocities are provided, the file may end here.
!           ==========================================================
            read(text,*,iostat=status) x, y, z
            if (status /= 0) then
               call ErrorHandler('readPOSCAR','Read x, y, z coordinates',trim(text))
            endif
            n = n + 1
            if (cart) then
!              =======================================================
!              Cartesian coordinates
!              =======================================================
               AtomPosition(1,n) = fac1*x*Angstrom2Bohr
               AtomPosition(2,n) = fac2*y*Angstrom2Bohr
               AtomPosition(3,n) = fac3*z*Angstrom2Bohr
            else
!              =======================================================
!              Direct coordinates
!              =======================================================
               AtomPosition(1:3,n) = BravaisLatticeVec(1:3,1)*x+BravaisLatticeVec(1:3,2)*y+BravaisLatticeVec(1:3,3)*z
            endif
         endif
         if (n == num_atoms) then
            exit LOOP_do
         endif
      enddo LOOP_do
!
      if (n < num_atoms) then
!        -------------------------------------------------------------
         call ErrorHandler('readPOSCAR','Position data POSCAR is incomplete')
!        -------------------------------------------------------------
      endif
!
      deallocate(atom_name_of_type)
      deallocate(num_atoms_of_type)
!     ----------------------------------------------------------------
      call endString()
!     ----------------------------------------------------------------
!
      close(funit)
!
      write(6,'(/,a,/)')' The Position data POSCAR is read successfully!'
   endif
!
!  -------------------------------------------------------------------
   call bcastMessage(AtomPosition,3,num_atoms,IOPE)
!  -------------------------------------------------------------------
!
   nullify( pa, BravaisLatticeVec, AtomNameCharacter, AtomicNumber, AtomPosition )
!
   end subroutine readPOSCAR
