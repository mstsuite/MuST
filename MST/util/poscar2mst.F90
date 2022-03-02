!  convert POTCAR to the position data format acceptable by MST
program poscar2mst
   use KindParamModule, only : IntKind, RealKind
   use Matrix3dModule, only : invm3
   use StringModule, only : initString, endString, getNumTokens, setString
   use PhysParamModule, only : Angstrom2Bohr
   use MathParamModule, only : ONE, TEN2m6
   implicit none
!
   character (len=80) :: text, tmp
   character (len=3), allocatable :: atom_name(:)
! 
   logical :: cart = .false.
!
   real (kind=RealKind) :: b0(3), b1(3), fac, a0
   real (kind=RealKind) :: bra(3,3), brainv(3,3)
!
   integer (kind=IntKind) :: i, j, k, n, m, status
   integer (kind=IntKind) :: num_atom_types, tot_num_atoms
   integer (kind=IntKind), allocatable :: num_atoms(:)
!
   interface
      function getTokenPosition(k,s,m) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer, intent(out), optional :: m
         integer :: p
      end function getTokenPosition
   end interface
!
   interface
      function isNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isNumber
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
   interface
      function trim_string(s1,c) result(s2)
         character (len=*), intent(in) :: s1
         character (len=len(s1)) :: s2
         character (len=1), optional, intent(in) :: c
      end function trim_string
   end interface
!
   call initString(80)
!
   a0 = ONE
   tot_num_atoms = 0
   j = 0
   do 
      read(5,'(a)',iostat=status) text
      if (status /= 0) then
         exit
      else
         tmp = adjustl(text)
         if (tmp(1:1) == '#' .or. tmp(1:1) == '!') then
            write(6,'(a)')text
            cycle
         endif
      endif
      text = trim_string(tmp,'#')
      if (j <= 3) then
         if (j == 0) then
            if (isNumber(text)) then
               read(text,*,iostat=status)fac
               if (status /= 0) then
                  write(6,'(a)')text
                  stop 'Error'
               else if (fac < TEN2m6) then
                  write(6,'(a,d15.8)')'# WARNING: the scaling factor < 1.0^-6. Check the results!'
                  fac = ONE
               endif
               write(6,'(a)')'#****************************************************'
               write(6,'(a)')'# Insert the scaling factor in the following line:  *'
               write(6,'(a)')'#     The length will be in atomic units            *'
               write(6,'(a)')'#****************************************************'
            else
               write(6,'(a,a)')'# ',trim(text)
               cycle
            endif
         else
            read(text,*)bra(1:3,j)
            if (j == 1) then
               if (abs(fac*bra(1,1)) > ONE) then
                  a0 = abs(bra(1,1))
               else if (abs(fac*bra(2,1)) > ONE) then
                  a0 = abs(bra(2,1))
               else if (abs(fac*bra(3,1)) > ONE) then
                  a0 = abs(bra(3,1))
               endif
               write(6,'(1d15.8)')a0*fac*Angstrom2Bohr
               write(6,'(a)')' '
               write(6,'(a)')'# Bravais Lattice..........'
            endif
            write(6,'(6x,3d24.16)')bra(1:3,j)/a0
            if (j == 3) then
               write(6,'(a)')' '
               write(6,'(a)')'CartCoordinates'
 !             call invm3(bra,brainv,status)
            endif
         endif
      else if (j == 4) then
         call setString(text)
         num_atom_types = getNumTokens()
         allocate(atom_name(num_atom_types))
         allocate(num_atoms(num_atom_types))
!        =============================================================
!        check the format of the line: whether it contains atoms names
!                                      or atoms nnumber
!        =============================================================
         i = getTokenPosition(1,text)
         k = getTokenPosition(2,text)
         if (isNumber(text(i:k-1))) then
            do n = 1, num_atom_types-1
               i = getTokenPosition(n,text)
               k = getTokenPosition(n+1,text)
               read(text(i:k-1),*)num_atoms(n)
               if (n < 10) then
                  write(atom_name(n),'(a,i2)')'A',10+n
                  atom_name(n)(2:2)='0'
               else
                  write(atom_name(n),'(a,i2)')'A',n
               endif
            enddo
            read(text(k:),*)num_atoms(num_atom_types)
            if (num_atom_types < 10) then
               write(atom_name(num_atom_types),'(a,i2)')'A',10+num_atom_types
               atom_name(num_atom_types)(2:2)='0'
            else
               write(atom_name(num_atom_types),'(a,i2)')'A',num_atom_types
            endif
            j = j + 1
         else
            do n = 1, num_atom_types-1
               i = getTokenPosition(n,text)
               k = getTokenPosition(n+1,text,m)
               if (m > 2) then
                  write(6,'(a,a)')'Invalid element name: ',trim(text(i:k-1))
                  stop 'Error'
               else
                  atom_name(n) = trim(text(i:k-1))
               endif
            enddo
            atom_name(num_atom_types) = trim(text(k:))
         endif
      else if (j == 5) then
         call setString(text)
         do n = 1, num_atom_types-1
            i = getTokenPosition(n,text)
            k = getTokenPosition(n+1,text)
            read(text(i:k-1),*,iostat=status)num_atoms(n)
            if (status /= 0) then
               write(6,'(a,2i3,2x,a)')'string range: ',i,k-1,text
               stop 'Error!'
            endif
         enddo
         read(text(k:),*)num_atoms(num_atom_types)
      else if (nocaseCompare(text(1:3),'Sel')) then
         cycle
      else if (nocaseCompare(text(1:1),'C')) then
         cart = .true.
         i = 0
         n = 1
      else if (nocaseCompare(text(1:1),'D')) then
         cart = .false.
         i = 0
         n = 1
      else if (n <= num_atom_types) then
         tot_num_atoms = tot_num_atoms + 1
         i = i + 1
         read(text,*) b0(1:3)
         if (cart) then
            b1 = b0
         else
!           b1(1:3) = brainv(1:3,1)*b0(1)+brainv(1:3,2)*b0(2)+brainv(1:3,3)*b0(3)
            b1(1:3) = bra(1:3,1)*b0(1)+bra(1:3,2)*b0(2)+bra(1:3,3)*b0(3)
         endif
         write(6,'(2x,a,2x,3d24.16)')atom_name(n),b1(1:3)/a0
!        b0(1:3) = bra(1:3,1)*b1(1)+bra(1:3,2)*b1(2)+bra(1:3,3)*b1(3)
!        write(6,'(a,2x,3d24.16)')'check: ',b0(1:3)*fac
         if (i == num_atoms(n)) then
            i = 0
            n = n + 1
         endif
      else
         exit
      endif
      j = j + 1
   enddo
!
   call endString()
!
   stop 'Ok'
end program poscar2mst
