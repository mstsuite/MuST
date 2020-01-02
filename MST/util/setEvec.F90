!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   program setEvec
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, HALF, ONE, TEN2m6, TWO, PI, PI2
!
   use MPPModule, only : initMPP, endMPP
!
   use StringModule, only : initString, endString, getNumTokens, readToken
!
   use DataServiceCenterModule, only : initDataServiceCenter,         &
                                       createDataStorage,             &
                                       endDataServiceCenter
!
   use DataServiceCenterModule, only : getDataStorage,                &
                                       RealType, RealMark,            &
                                       IntegerType, IntegerMark,      &
                                       CharacterType, CharacterMark
!
   implicit none
!
   logical :: outside
!
   character (len=80) :: text
   character (len=60) :: position_dat = ' '
   character (len=60) :: boundary_dat = ' '
   character (len=60) :: corner_dat = ' '
   character (len=60) :: direction_dat = ' '
   character (len=60) :: shift_x = ' '
   character (len=60) :: shift_y = ' '
   character (len=60) :: shift_z = ' '
   character (len=60) :: evec_in_dat = ' '
   character (len=60) :: evec_out_dat = ' '
   character (len=60) :: skey
   character (len=60) :: sthe_max = ' '
   character (len=60) :: sthe_min = ' '
   character (len=60) :: sphi_max = ' '
   character (len=60) :: sphi_min = ' '
   character (len=10) :: atom
   character (len=1) :: dummy
   character (len=3) :: AtomName
   character (len=1), pointer :: AtomicName(:)
!
   logical :: ElementName, Jump1stLine
!
   integer (kind=IntKind), parameter :: sunit = 5
   integer (kind=IntKind), parameter :: funit = 11
   integer (kind=IntKind), parameter :: MaxCorners = 20
!
   integer (kind=IntKind) :: num_atoms
   integer (kind=IntKind) :: ios, status, n, i, jdx
   integer (kind=IntKind) :: ran_in, ran_out, iseed
!
   integer (kind=IntKind), pointer :: AtomicNumber(:)
!
   real (kind=RealKind), pointer :: AtomPosition(:,:)
   real (kind=RealKind), pointer :: AtomEvec(:,:)
   real (kind=RealKind), pointer :: ConstrField(:,:)
!
   real (kind=RealKind) :: radius
   real (kind=RealKind) :: cos_theta_min, cos_theta_max, cos_theta_range
   real (kind=RealKind) :: theta_min, theta_max
   real (kind=RealKind) :: phi_min, phi_max, phi_range
   real (kind=RealKind) :: alp_in, alp_out, alpev
   real (kind=RealKind) :: CenterX, CenterY, CenterZ
   real (kind=RealKind) :: r, f, ran2, cos_theta, phi
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
!  -------------------------------------------------------------------
   call initMPP()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  call initDataServiceCenter to startup the Data Storage Service
!  -------------------------------------------------------------------
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
         read(text(jdx+2:),*)num_atoms
         write(6,'(/,a,i6)')'Number of atoms: ',num_atoms
      else if (skey == 'Transition Radius') then
         read(text(jdx+2:),*)radius
         write(6,'(/,a,f10.5,a)')'Transition radius = ',radius,' Angstrom'
         write(6,'(  a,f10.5,a)')'                  = ',radius*anstr2au,' a.u.'
         radius = radius*anstr2au
      else if (skey == 'Inside Random Setting') then
         read(text(jdx+2:),*)ran_in
      else if (skey == 'Outside Random Setting') then
         read(text(jdx+2:),*)ran_out
      else if (skey == 'Random Function Seed') then
         read(text(jdx+2:),*)iseed
      else if (skey == 'Inside Mixing Param') then
         read(text(jdx+2:),*)alp_in
      else if (skey == 'Outside Mixing Param') then
         read(text(jdx+2:),*)alp_out
      else if (skey == 'Position Data') then
         position_dat = adjustl(text(jdx+2:))
      else if (skey == 'Boundary Data') then
         boundary_dat = adjustl(text(jdx+2:))
      else if (skey == 'Corner Data') then
         corner_dat = adjustl(text(jdx+2:))
      else if (skey == 'Direction Data') then
         direction_dat = adjustl(text(jdx+2:))
      else if (skey == 'Output Evec Data') then
         evec_out_dat = adjustl(text(jdx+2:))
      else if (skey == 'Input Evec Data') then
         evec_in_dat = adjustl(text(jdx+2:))
      else if (skey == 'Shift Center X') then
         shift_x = adjustl(text(jdx+2:))
      else if (skey == 'Shift Center Y') then
         shift_y = adjustl(text(jdx+2:))
      else if (skey == 'Shift Center Z') then
         shift_z = adjustl(text(jdx+2:))
      else if (skey == 'Minimum Theta (in the units of PI)') then
         sthe_min = adjustl(text(jdx+2:))
      else if (skey == 'Maximum Theta (in the units of PI)') then
         sthe_max = adjustl(text(jdx+2:))
      else if (skey == 'Minimum Phi (in the units of PI)') then
         sphi_min = adjustl(text(jdx+2:))
      else if (skey == 'Maximum Phi (in the units of PI)') then
         sphi_max = adjustl(text(jdx+2:))
      endif
   enddo
!
   if (position_dat == ' ') then
      stop 'Error: Position data is unknown!'
   endif
!
   call readPositionData(position_dat,NumAtomsIn=num_atoms)
!
   AtomicName => getDataStorage('Atomic Name',3*num_atoms,CharacterMark)
   AtomPosition => getDataStorage('Atomic Position',3,num_atoms,RealMark)
!
   CenterX = ZERO
   CenterY = ZERO
   CenterZ = ZERO
   do i=1,num_atoms
      CenterX = CenterX + AtomPosition(1,i)/real(num_atoms,kind=RealKind)
      CenterY = CenterY + AtomPosition(2,i)/real(num_atoms,kind=RealKind)
      CenterZ = CenterZ + AtomPosition(3,i)/real(num_atoms,kind=RealKind)
   enddo
   write(6,'(a,f15.8)')'CenterX = ',CenterX
   write(6,'(a,f15.8)')'CenterY = ',CenterY
   write(6,'(a,f15.8)')'CenterZ = ',CenterZ
!
   if (evec_in_dat /= ' ') then
      call readMomentDirectionData(evec_in_dat,num_atoms)
      AtomEvec => getDataStorage('Moment Direction',3,num_atoms,RealMark)
      ConstrField => getDataStorage('Constrain Field',3,num_atoms,RealMark)
   else
      call createDataStorage('Moment Direction',3*num_atoms,RealType)
      call createDataStorage('Constrain Field',3*num_atoms,RealType)
      AtomEvec => getDataStorage('Moment Direction',3,num_atoms,RealMark)
      ConstrField => getDataStorage('Constrain Field',3,num_atoms,RealMark)
      AtomEvec = ZERO
      ConstrField = ZERO
   endif
!
!
   if (shift_x /= ' ') then
      read(shift_x,*)CenterX
      CenterX = CenterX*anstr2au
   else
      CenterX = ZERO
   endif
!
   if (shift_y /= ' ') then
      read(shift_y,*)CenterY
      CenterZ = CenterX*anstr2au
   else
      CenterY = ZERO
   endif
!
   if (shift_z /= ' ') then
      read(shift_z,*)CenterZ
      CenterZ = CenterX*anstr2au
   else
      CenterZ = ZERO
   endif
!
   if (sthe_min /= ' ') then
      read(sthe_min,*)theta_min
   else
      theta_min = ZERO
   endif
   if (sthe_max /= ' ') then
      read(sthe_max,*)theta_max
   else
      theta_min = ONE
   endif
   if (sphi_min /= ' ') then
      read(sphi_min,*)phi_min
   else
      phi_min = ZERO
   endif
   if (sphi_max /= ' ') then
      read(sphi_max,*)phi_max
   else
      phi_min = TWO
   endif
!
   phi_range = (phi_max - phi_min)*PI
   phi_min = phi_min*PI
   phi_max = phi_max*PI
   cos_theta_min = cos(theta_max*PI)
   cos_theta_max = cos(theta_min*PI)
   cos_theta_range = cos_theta_max - cos_theta_min
!
   f = ran2(iseed)
!
   open(unit=funit,file=evec_out_dat,form='formatted',status='unknown')
   write(funit,'(a)')'! ****************************************************'
   write(funit,'(a)')'! *                                                  *'
   write(funit,'(a)')'! *  This data file is created by setEvec            *'
   write(funit,'(a)')'! *                                                  *'
   write(funit,'(a)')'! ****************************************************'
   LOOP_i: do i=1,num_atoms
      r = sqrt( (AtomPosition(1,i) - CenterX)**2                      &
               +(AtomPosition(2,i) - CenterY)**2                      &
               +(AtomPosition(3,i) - CenterZ)**2 )
      if ((r > radius .and. ran_out == 1) .or.                        &
          (r <= radius .and. ran_in == 1)) then
!        cos_theta = cos_theta_min
!        do while (cos_theta < HALF)
         cos_theta = cos_theta_range*ran2(iseed)+cos_theta_min
!        enddo
         phi=phi_range*ran2(iseed)+phi_min
         AtomEvec(1,i)=cos(phi)*sqrt(ONE-cos_theta*cos_theta)
         AtomEvec(2,i)=sin(phi)*sqrt(ONE-cos_theta*cos_theta)
         AtomEvec(3,i)=cos_theta
      endif
!
      if (r > radius) then
         alpev = alp_out
      else
         alpev = alp_in
      endif
!
!     ----------------------------------------------------------------
      call copyCharArray2String(AtomicName((i-1)*3+1:i*3),AtomName,3)
!     ----------------------------------------------------------------
      write(funit,'(a3,2x,6f22.13,1f10.5)')AtomName,                  &
                    AtomEvec(1:3,i), ConstrField(1:3,i),alpev
   enddo LOOP_i
   close(funit)
!
   nullify( AtomicName )
   nullify( AtomPosition )
   nullify( AtomEvec, ConstrField )
!
   call endDataServiceCenter()
   call endMPP()
!
   end program setEvec
!  ===================================================================
