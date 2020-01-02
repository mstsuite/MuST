!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   program charge
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, HALF, ONE, TEN2m6
!
   use SortModule, only : HeapSort
!
   use ChemElementModule, only : getName
!
   use StringModule, only : initString, endString, getNumTokens, readToken
!
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, &
                                       isDataStorageExisting, &
                                       getDataStorage, IntegerMark, RealMark
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use MPPModule, only : initMPP, endMPP
!
   implicit none
!
   logical :: outside
!
   character (len=80) :: text
   character (len=60) :: natom_dat = ' '
   character (len=60) :: position_dat = ' '
   character (len=60) :: charge_dat   = ' '
   character (len=60) :: sorted_dat   = ' '
   character (len=60) :: boundary_dat = ' '
   character (len=60) :: madelung_dat = ' '
   character (len=60) :: complete_dat = ' '
   character (len=60) :: corner_dat = ' '
   character (len=60) :: direction_dat = ' '
   character (len=60) :: shift_x = ' '
   character (len=60) :: shift_y = ' '
   character (len=60) :: shift_z = ' '
   character (len=60) :: neighbor_cut = ' '
   character (len=60) :: number_of_cuts = ' '
   character (len=60) :: skey
   character (len=10) :: atom
   character (len=1) :: dummy
!
   logical :: ElementName, Jump1stLine
!
   integer (kind=IntKind), parameter :: sunit = 5
   integer (kind=IntKind), parameter :: funit = 11
   integer (kind=IntKind), parameter :: MaxCorners = 20
!
   integer (kind=IntKind) :: num_atoms
   integer (kind=IntKind) :: ios, status, n, alen, i, j, k, m, jdx, sign_flag
   integer (kind=IntKind) :: nc, k0, k1, k2, ns, num_shells
   integer (kind=IntKind) :: DataForm
   integer (kind=IntKind) :: StartingLine, NB
!
   integer (kind=IntKind), pointer :: AtomicNumber(:)
   integer (kind=IntKind), allocatable :: idx(:), idx_inv(:)
   integer (kind=IntKind), allocatable :: unlike(:,:)
!
   real (kind=RealKind), pointer :: AtomPosition(:,:)
   real (kind=RealKind), pointer :: AtomPositionX(:)
   real (kind=RealKind), pointer :: AtomPositionY(:)
   real (kind=RealKind), pointer :: AtomPositionZ(:)
   real (kind=RealKind), allocatable :: rdist(:)
   real (kind=RealKind), allocatable :: chg(:), mom(:), mad(:), rcut(:)
!
   real (kind=RealKind) :: CenterX, CenterY, CenterZ
   real (kind=RealKind) :: r, dq, moment, mmt
   real (kind=RealKind) :: a, b, c, d, epsi, rd, x0, y0, z0
   real (kind=RealKind) :: vx(MaxCorners), vy(MaxCorners), vz(MaxCorners)
   real (kind=RealKind) :: cc1(3), cc2(3), ccp(3), orig_sign, vdotv_sign
   real (kind=RealKind), pointer :: BravaisLatticeVec(:,:)
   real (kind=RealKind) :: x, y, z, r1, r2
   real (kind=RealKind) :: bshiftx(27)
   real (kind=RealKind) :: bshifty(27)
   real (kind=RealKind) :: bshiftz(27)
!
   real (kind=RealKind), parameter :: anstr2au = 1.88973d0
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
   interface
      subroutine readPositionData(fname,NumAtomsIn,NumAtomsOut)
         use KindParamModule, only : IntKind
         character (len=*), intent(in) :: fname
         integer (kind=IntKind), intent(in), optional :: NumAtomsIn
         integer (kind=IntKind), intent(out), optional :: NumAtomsOut
      end subroutine readPositionData
   end interface
!
   call initMPP()
   call initDataServiceCenter()
!
!  ===================================================================
!  get the input data from standard input (unit = 5)
!  ===================================================================
   do 
      read(sunit,'(a)',iostat=status)text
      if (status < 0) then
         exit
      endif
!
      text=adjustl(text)
      if (text(1:1) == '#' .or. text(1:1) == '!'                      &
                           .or. len_trim(text) == 0) then
         cycle
      endif
!
      jdx = index(text,'::')
      if (jdx > 1) then
         skey = text(:jdx-1)
      else
         write(6,'(a)')trim(text)
         stop 'Error: invalid input format'
      endif
!
      if (skey == 'Number of Atoms') then
         natom_dat = adjustl(text(jdx+2:))
         read(natom_dat,*)num_atoms
      else if (skey == 'Position Data') then
         position_dat = adjustl(text(jdx+2:))
      else if (skey == 'Charge Table') then
         charge_dat = adjustl(text(jdx+2:))
      else if (skey == 'Boundary Data') then
         boundary_dat = adjustl(text(jdx+2:))
      else if (skey == 'Corner Data') then
         corner_dat = adjustl(text(jdx+2:))
      else if (skey == 'Direction Data') then
         direction_dat = adjustl(text(jdx+2:))
      else if (skey == 'Sorted Output') then
         sorted_dat = adjustl(text(jdx+2:))
      else if (skey == 'Shift Center X') then
         shift_x = adjustl(text(jdx+2:))
      else if (skey == 'Shift Center Y') then
         shift_y = adjustl(text(jdx+2:))
      else if (skey == 'Shift Center Z') then
         shift_z = adjustl(text(jdx+2:))
      else if (skey == 'Number of Cuts') then
         number_of_cuts = adjustl(text(jdx+2:))
      else if (skey == 'Neighbor Cutoff') then
         neighbor_cut = adjustl(text(jdx+2:))
      else if (skey == 'Madelung Data') then
         madelung_dat = adjustl(text(jdx+2:))
      else if (skey == 'Complete Output') then
         complete_dat = adjustl(text(jdx+2:))
      endif
   enddo
!
   write(6,'(a,i6)')'Total Number of Atoms: ',num_atoms
!
   if (position_dat == ' ') then
      stop 'Error: Position data is unknown!'
   else if (charge_dat == ' ') then
      stop 'Error: Charge distribution data is unknown'
   endif
!
!  ===================================================================
!  Determine the format of the position data
!  ===================================================================
   call readPositionData(position_dat,NumAtomsIn=num_atoms)
!
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      BravaisLatticeVec => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('charge','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
   if (isDataStorageExisting('Atomic Position')) then
      AtomPosition => getDataStorage('Atomic Position',num_atoms,3,RealMark)
      AtomPositionX => AtomPosition(1:num_atoms,1)
      AtomPositionY => AtomPosition(1:num_atoms,2)
      AtomPositionZ => AtomPosition(1:num_atoms,3)
   else
!     ----------------------------------------------------------------
      call ErrorHandler('charge','Atom position vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
   if (isDataStorageExisting('Atomic Number')) then
      AtomicNumber => getDataStorage('Atomic Number',num_atoms,IntegerMark)
   else
!     ----------------------------------------------------------------
      call ErrorHandler('charge','Atomic number data does not exist')
!     ----------------------------------------------------------------
   endif
!
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
   write(6,'(a,f15.8)')'CenterX = ',CenterX
   write(6,'(a,f15.8)')'CenterY = ',CenterY
   write(6,'(a,f15.8)')'CenterZ = ',CenterZ
!
   open(unit=10,file=charge_dat,status='old')
   read(10,'(a)')dummy
   read(10,'(a)')dummy
   read(10,'(a)')dummy
   read(10,'(a)')dummy
!
   allocate( rdist(num_atoms) )
   allocate( chg(num_atoms), mom(num_atoms), mad(num_atoms) )
!
   if (shift_x /= ' ') then
      read(shift_x,*)CenterX
   else
      CenterX = ZERO
   endif
!
   if (shift_y /= ' ') then
      read(shift_y,*)CenterY
   else
      CenterY = ZERO
   endif
!
   if (shift_z /= ' ') then
      read(shift_z,*)CenterZ
   else
      CenterZ = ZERO
   endif
!
   LOOP_i: do i=1,num_atoms
      rdist(i) = sqrt( (AtomPositionX(i) - CenterX)**2                &
                      +(AtomPositionY(i) - CenterY)**2                &
                      +(AtomPositionZ(i) - CenterZ)**2 )
      read(10,'(46x,3(2x,f9.5))')chg(i),mmt,mom(i)
   enddo LOOP_i
!
   close(10)
!
   allocate( idx(num_atoms), idx_inv(num_atoms) )
   call HeapSort(num_atoms,rdist,idx)
   do i=1,num_atoms
      j=idx(i)
      idx_inv(j)=i
   enddo
!
   if (number_of_cuts /= ' ') then
      read(number_of_cuts,*)num_shells
      if (num_shells < 1) then
         stop 'Error: num_shells < 1'
      endif
   else
      num_shells = 1
   endif
   allocate( unlike(num_atoms,num_shells), rcut(0:num_shells) )
   rcut(0) = 0.0d0
!
   if (neighbor_cut /= ' ') then
      read(neighbor_cut,*)rcut(1:num_shells)
      write(6,'(a,10f10.5)')'Neighbor Cut Radius = ',rcut(1:num_shells)
      n = 0
      do k = -1, 1
         do j = -1, 1
            do i = -1, 1
               n = n + 1
               bshiftx(n) = i*BravaisLatticeVec(1,1)+                 &
                            j*BravaisLatticeVec(1,2)+k*BravaisLatticeVec(1,3)
               bshifty(n) = i*BravaisLatticeVec(2,1)+                 &
                            j*BravaisLatticeVec(2,2)+k*BravaisLatticeVec(2,3)
               bshiftz(n) = i*BravaisLatticeVec(3,1)+                 &
                            j*BravaisLatticeVec(3,2)+k*BravaisLatticeVec(3,3)
            enddo
         enddo
      enddo
      do ns = 1, num_shells
         do i = 1, num_atoms
            unlike(i,ns) = 0
            do j = 1, num_atoms
               if (AtomicNumber(i) == AtomicNumber(j)) then
                  cycle
               else
                  x0 = AtomPositionX(j) - AtomPositionX(i)
                  y0 = AtomPositionY(j) - AtomPositionY(i)
                  z0 = AtomPositionZ(j) - AtomPositionZ(i)
                  do k = 1, 27
                     x = x0 + bshiftx(k)
                     y = y0 + bshifty(k)
                     z = z0 + bshiftz(k)
                     r = sqrt(x*x+y*y+z*z)
                     if (r <= rcut(ns) .and. r > rcut(ns-1)) then
                        unlike(i,ns) = unlike(i,ns) + 1
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo
   endif
!
   if (sorted_dat /= ' ' .and. neighbor_cut /= ' ') then
      open(unit=21,file=sorted_dat,form='formatted',status='unknown')
      write(21,'(75(''=''))')
      write(21,'(a)')                  &
! '  Z     Distance(A)    Unlike Neighbor     Excess Electron         Moment'
  ' Atom   Distance(A)    Unlike Neighbor     Excess Electron         Moment'
      write(21,'(75(''-''))')
      do i = 1, num_atoms
         j = idx(i)
!        write(21,'(i3,4x,f10.5,11x,i3,13x,f10.5,9x,f10.5)')       &
!            AtomicNumber(j), rdist(i)/anstr2au, unlike(j,1), chg(j), mom(j)
         write(21,'(1x,a2,4x,f10.5,11x,i3,13x,f10.5,9x,f10.5)')    &
             getName(AtomicNumber(j)), rdist(i)/anstr2au, unlike(j,1), chg(j), mom(j)
      enddo
      close(21)
   else if (sorted_dat /= ' ') then
      open(unit=21,file=sorted_dat,form='formatted',status='unknown')
      write(21,'(66(''=''))')
      write(21,'(a)')                  &
!     '      Z          Distance(A)      Excess Electron         Moment'
      '    Atom         Distance(A)      Excess Electron         Moment'
      write(21,'(66(''-''))')
      do i = 1, num_atoms
         j = idx(i)
!        write(21,'(5x,i3,3(9x,f10.5))')AtomicNumber(j), rdist(i)/anstr2au, &
         write(21,'(6x,a2,3(9x,f10.5))')getName(AtomicNumber(j)), rdist(i)/anstr2au, &
                                        chg(j), mom(j)
      enddo
      close(21)
   else
      open(unit=21,file='qmvsr.dat',form='formatted',status='unknown')
      write(21,'(66(''=''))')
      write(21,'(a)')                  &
!     '      Z          Distance(A)      Excess Electron         Moment'
      '    Atom         Distance(A)      Excess Electron         Moment'
      write(21,'(66(''-''))')
      do i = 1, num_atoms
         j = idx_inv(i)
         write(21,'(6x,a2,3(9x,f10.5))')getName(AtomicNumber(i)), rdist(j)/anstr2au, &
                                        chg(i), mom(i)
      enddo
      close(21)
   endif
!
   if (madelung_dat /= ' ' .and. complete_dat /= ' ' .and.       &
       neighbor_cut /= ' ') then
      open(unit=19,file=trim(madelung_dat),form='formatted',status='old')
      read(19,'(a)')dummy
      read(19,'(a)')dummy
      read(19,'(a)')dummy
      read(19,'(a)')dummy
      do i = 1, num_atoms
         read(19,'(48x,f9.5)')mad(i)
      enddo
      close(unit=19)
!
      open(unit=22,file=trim(complete_dat),form='formatted',status='unknown')
      write(22,'(75(''=''))')
      if (num_shells > 1) then
         write(22,'(a)')                  &
!        '  Z     Dist(A)    1st-UL    2nd-UL      dQ        Moment       Vmad'
         ' Atom   Dist(A)    1st-UL    2nd-UL      dQ        Moment       Vmad'
         write(22,'(75(''-''))')
         do j = 1, num_atoms
            i = idx_inv(j)
!           write(22,'(i3,3x,f10.5,2(4x,i3,2x),3(2x,f10.5))')       &
!              AtomicNumber(j),rdist(i)/anstr2au,unlike(j,1),unlike(j,2), &
            write(22,'(1x,a2,3x,f10.5,2(4x,i3,2x),3(2x,f10.5))')    &
               getName(AtomicNumber(j)),rdist(i)/anstr2au,unlike(j,1),unlike(j,2), &
               chg(j),mom(j),mad(j)
         enddo
      else 
         write(22,'(a)')                  &
!        '  Z     Dist(A)    1st-UL     dQ        Moment       Vmad'
         ' Atom   Dist(A)    1st-UL     dQ        Moment       Vmad'
         write(22,'(75(''-''))')
         do j = 1, num_atoms
            i = idx_inv(j)
!           write(22,'(i3,3x,f10.5,4x,i3,2x,3(2x,f10.5))')          &
!              AtomicNumber(j),rdist(i)/anstr2au,unlike(j,1),       &
            write(22,'(1x,a2,3x,f10.5,4x,i3,2x,3(2x,f10.5))')          &
               getName(AtomicNumber(j)),rdist(i)/anstr2au,unlike(j,1), &
               chg(j),mom(j),mad(j)
         enddo
      endif
      close(22)
   endif
!
!  ===================================================================
!  analyze the charge and moment distribution on the surfaces of a given 
!  volume
!  ===================================================================
   if (boundary_dat /= ' ') then
      open(unit = 12,file=boundary_dat,form='formatted',status='old')
!
      do 
         read(12,'(a)',iostat=status)text
         if (status < 0) then
            exit
         endif
         read(text,*)a,b,c,d,epsi
         write(6,'(//,a)')'Charge and Moment distribution on the surface'
         write(6,'(4f10.5)')a,b,c,d
         write(6,'(75(''=''))')
         write(6,'(a)')                  &
  ' Z         x           y           z      Distance(A)     -dQ        Moment'
         write(6,'(75(''-''))')
         read(12,*)nc
         if (nc > MaxCorners) then
            write(6,'(a,i5)')'Too many boundary corners for the surface: ',nc
            stop 'Error!'
         endif
         r1 = d - epsi
         r2 = r1 + 2.0d0*epsi
         x0 = ZERO; y0 = ZERO; z0 = ZERO
         do k=1,nc
            read(12,*)i,vx(k),vy(k),vz(k)
            r=sqrt(vx(k)**2+vy(k)**2+vz(k)**2) + epsi
            if (r > r2) then
               r2 = r
            endif
            x0 = x0 + vx(k)/real(nc,kind=RealKind)
            y0 = y0 + vy(k)/real(nc,kind=RealKind)
            z0 = z0 + vz(k)/real(nc,kind=RealKind)
         enddo
!        -------------------------------------------------------------
         k1 = 1
         call hunt(num_atoms,rdist,r1,k1)
         if (k1 == 0) then
            k1 = 1
         endif
         k2 = k1 + 1
         call hunt(num_atoms,rdist,r2,k2)
         if (k2 < num_atoms) then
            k2 = k2 + 1
         endif
!        -------------------------------------------------------------
         do i = k1, k2
            j = idx(i)
            rd = (AtomPositionX(j)-CenterX)*a                         &
                +(AtomPositionY(j)-CenterY)*b                         &
                +(AtomPositionZ(j)-CenterZ)*c
            if (abs(rd-d) > epsi) then
               cycle
            endif
!
!           ==========================================================
!           check if the point is inside the boundaries defined by the
!           corners.            
!           ==========================================================
            outside = .false.
            sign_flag = 0
            orig_sign = ZERO
            LOOP_k: do k=1,nc
               if (k == nc) then
                  k0 = 1
               else
                  k0 = k + 1
               endif
!
!              =======================================================
!              if the point is near the edge within a tolenrent range,
!              we consider the point is inside the boundaries.
!              =======================================================
               cc1(1) = AtomPositionX(j)-CenterX-vx(k0)
               cc1(2) = AtomPositionY(j)-CenterY-vy(k0)
               cc1(3) = AtomPositionZ(j)-CenterZ-vz(k0)
               cc2(1) = vx(k)-vx(k0)
               cc2(2) = vy(k)-vy(k0)
               cc2(3) = vz(k)-vz(k0)
               r = sqrt(cc2(1)*cc2(1)+cc2(2)*cc2(2)+cc2(3)*cc2(3))
               if (r < TEN2m6) then
                  write(6,'(a,d15.8)')'Distance between vertices is near 0: ',r
                  stop 'Error!'
               endif
!              -------------------------------------------------------
               call vcross(cc2,cc1,ccp)
!              -------------------------------------------------------
               if ((ccp(1)**2+ccp(2)**2+ccp(3)**2)/r < epsi) then
                  outside = .false.
                  exit LOOP_k
               endif
!
!              =======================================================
!              For the point far away from the edge, we need to check
!              if the point is inside the boundaries.
!              =======================================================
               cc1(1) = vx(k)-x0
               cc1(2) = vy(k)-y0
               cc1(3) = vz(k)-z0
               cc2(1) = vx(k0)-x0
               cc2(2) = vy(k0)-y0
               cc2(3) = vz(k0)-z0
!              -------------------------------------------------------
               call vcross(cc2,cc1,ccp)
!              -------------------------------------------------------
               if (ccp(1)*ccp(1)+ccp(2)*ccp(2)+ccp(3)*ccp(3) >= TEN2m6) then
                  vdotv_sign = sign(ONE,ccp(1)*a+ccp(2)*b+ccp(3)*c)
                  if (sign_flag == 0) then
                     orig_sign = vdotv_sign
                     sign_flag = 1
                  else if (abs(orig_sign-vdotv_sign) > TEN2m6) then
                     outside = .true.
                  endif
               endif
            enddo LOOP_k
            if (.not.outside) then
!              write(6,'(i3,6(2x,f10.5))')AtomicNumber(j),            &
               write(6,'(1x,a2,6(2x,f10.5))')getName(AtomicNumber(j)),&
                                          AtomPositionX(j)/anstr2au,  &
                                          AtomPositionY(j)/anstr2au,  &
                                          AtomPositionZ(j)/anstr2au,  &
                                          rdist(i)/anstr2au, chg(j), mom(j)
            endif
         enddo
      enddo
!
      close (unit=12)
   endif
!
!  ===================================================================
!  analyze the charge and moment distribution around a point
!  ===================================================================
   if (corner_dat /= ' ') then
      open(unit = 13,file=corner_dat,form='formatted',status='old')
!
      do 
         read(13,'(a)',iostat=status)text
         if (status < 0) then
            exit
         endif
         read(text,*)x,y,z,epsi
!
         write(6,'(//,a)')'Charge and Moment distribution around the corner:'
         write(6,'(3f10.5)')x,y,z
         write(6,'(75(''=''))')
         write(6,'(a)')                                                   &
  ' Z         x           y           z      Distance(A)     -dQ        Moment'
         write(6,'(75(''-''))')
         r2= sqrt(x*x+y*y+z*z)+epsi
         r1= r2-2.0d0*epsi
!        -------------------------------------------------------------
         k1 = 1
         call hunt(num_atoms,rdist,r1,k1)
         if (k1 == 0) then
            k1 = 1
         endif
         k2 = k1 + 1
         call hunt(num_atoms,rdist,r2,k2)
         if (k2 < num_atoms) then
            k2 = k2 + 1
         endif
!        -------------------------------------------------------------
         do i = k1, k2
            j = idx(i)
            r = sqrt((AtomPositionX(j)-CenterX-x)**2 +                &
                     (AtomPositionY(j)-CenterY-y)**2 +                &
                     (AtomPositionZ(j)-CenterZ-z)**2)
            if (r <= epsi) then
!              write(6,'(i3,6(2x,f10.5))')AtomicNumber(j),            &
               write(6,'(1x,a2,6(2x,f10.5))')getName(AtomicNumber(j)),&
                                          AtomPositionX(j)/anstr2au,  &
                                          AtomPositionY(j)/anstr2au,  &
                                          AtomPositionZ(j)/anstr2au,  &
                                          rdist(i)/anstr2au, chg(j), mom(j)
            endif
         enddo
      enddo
!
      close (unit=13)
   endif
!
!  ===================================================================
!  analyze the charge and moment distribution along certain directions
!  ===================================================================
   if (direction_dat /= ' ') then
      open(unit = 14,file=direction_dat,form='formatted',status='old')
!
      do 
         read(14,'(a)',iostat=status)text
         if (status < 0) then
            exit
         endif
         read(text,*)x,y,z,r1,r2,epsi
         r = sqrt(x*x+y*y+z*z)
         if (r < TEN2m6) then
            write(6,'(a,3d15.5)')'Invalid direction parameters: ',x,y,z
            stop 'Error!'
         endif
         x = x/r; y = y/r; z = z/r
!
         write(6,'(//,a)')'Charge and Moment distribution along the direction:'
         write(6,'(3f10.5)')x,y,z
         write(6,'(75(''=''))')
         write(6,'(a)')                                               &
  ' Z         x           y           z      Distance(A)     -dQ        Moment'
         write(6,'(75(''-''))')
!        -------------------------------------------------------------
         k1 = 1
         call hunt(num_atoms,rdist,r1,k1)
         if (k1 == 0) then
            k1 = 1
         endif
         k2 = k1 + 1
         call hunt(num_atoms,rdist,r2,k2)
         if (k2 < num_atoms) then
            k2 = k2 + 1
         endif
!        -------------------------------------------------------------
         do i = k1, k2
            j = idx(i)
            r = sqrt((AtomPositionX(j)-CenterX)**2 +                  &
                     (AtomPositionY(j)-CenterY)**2 +                  &
                     (AtomPositionZ(j)-CenterZ)**2)
            rd = (AtomPositionX(j)-CenterX)*x                         &
                +(AtomPositionY(j)-CenterY)*y                         &
                +(AtomPositionZ(j)-CenterZ)*z
            if (abs(r-rd) <= epsi) then
               write(6,'(1x,a2,6(2x,f10.5))')getName(AtomicNumber(j)),&
                                          AtomPositionX(j)/anstr2au,  &
                                          AtomPositionY(j)/anstr2au,  &
                                          AtomPositionZ(j)/anstr2au,  &
                                          rdist(i)/anstr2au, chg(j), mom(j)
            endif
         enddo
      enddo
!
      close (unit=14)
   endif
!
   deallocate( rdist, idx, idx_inv, chg, mom )
   deallocate( unlike, rcut )
!
   call endMPP()
!
   end program charge
!  ===================================================================
