!  *******************************************************************
!  The evec data are usually in the following format:
!
!  ! comments
!  # comments
!  ! optional 3 lines of rotation matrix to be applied to the evec vectors:
!  rot(1,1)  rot(2,1)  rot(3,1)
!  rot(1,2)  rot(2,2)  rot(3,2)
!  rot(1,3)  rot(2,3)  rot(3,3)
!      atom_name      evec_x  evec_y  evec_z
!  ! or with 3 columns of constrained field (cf)
!      atom_name      evec_x  evec_y  evec_z  cf_x  cf_y  cf_z
!  ! or with 1 mixing parameter
!      atom_name      evec_x  evec_y  evec_z  cf_x  cf_y  cf_z  e_mix
!  ! or
!      atomic_number  evec_x  evec_y  evec_z
!  ! or with 3 columns of constrained field (cf)
!      atomic_number  evec_x  evec_y  evec_z  cf_x  cf_y  cf_z
!  ! or with 1 mixing parameter
!      atomic_number  evec_x  evec_y  evec_z  cf_x  cf_y  cf_z  e_mix
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readMomentDirectionData(fname,num_atoms)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ZERO, ONE, TEN2m8
   use MPPModule, only : MyPE, bcastMessage
   use ErrorHandlerModule, only : ErrorHandler
   use Matrix3dModule, only : determ3
!
   use DataServiceCenterModule, only : createDataStorage,             &
                                       getDataStorage,                &
                                       RealType, RealMark,            &
                                       IntegerType, IntegerMark,      &
                                       CharacterType, CharacterMark
!
   use ChemElementModule, only : getZtot, getName
!
   use StringModule, only : initString, setString, endString,         &
                            getNumTokens, readToken
!
   implicit none
!
   character (len=*), intent(in) :: fname
!
   character (len=240) :: text
   character (len=10) :: atom
   character (len=1) :: dummy
!
   logical :: ElementName, isDirect, RotateEvec
!
   integer (kind=IntKind), intent(in) :: num_atoms
!
   integer (kind=IntKind), parameter :: funit = 103
   integer (kind=IntKind) :: ios, status, n, alen, j, k, anum
   integer (kind=IntKind), pointer :: AtomicNumber(:)
   integer (kind=IntKind) :: DataForm
   integer (kind=IntKind) :: StartingLine, NB
   integer (kind=IntKind) :: nsn, nst
   integer (kind=IntKind), parameter :: MaxSingleNumbers = 2
   integer (kind=IntKind), parameter :: MaxSingleStrings = 1
!
   real (kind=RealKind) :: x, y, z
   real (kind=RealKind), pointer :: AtomEvec(:,:), p_AtomEvec(:)
   real (kind=RealKind), pointer :: ConstrField(:,:), p_ConstrField(:)
   real (kind=RealKind), pointer :: emix(:)
   real (kind=RealKind) :: RotationMatrix(3,3)
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
      function getTokenPosition(k,s,n) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer, intent(out), optional :: n
         integer :: p
      end function getTokenPosition
   end interface
!
   if (num_atoms < 1) then
      call ErrorHandler('readMomentDirectionData','num_atoms < 1',num_atoms)
   endif
!
!  -------------------------------------------------------------------
   call createDataStorage('Moment Direction',3*num_atoms,RealType)
   call createDataStorage('Constrain Field',3*num_atoms,RealType)
   call createDataStorage('Evec Mixing Parameters',num_atoms,RealType)
!  -------------------------------------------------------------------
   AtomicNumber => getDataStorage('Atomic Number',num_atoms,IntegerMark)
   AtomEvec => getDataStorage('Moment Direction',3,num_atoms,RealMark)
   ConstrField => getDataStorage('Constrain Field',3,num_atoms,RealMark)
   p_AtomEvec => getDataStorage('Moment Direction',3*num_atoms,RealMark)
   p_ConstrField => getDataStorage('Constrain Field',3*num_atoms,RealMark)
   emix => getDataStorage('Evec Mixing Parameters',num_atoms,RealMark)
!  -------------------------------------------------------------------
!
   if (MyPE == 0) then
      write(6,'(/,a,a)')' The Evec Data File: ',trim(fname)
!
      open(unit=funit,file=fname,form='formatted',status='old',iostat=ios, &
           action='read')
      if (ios > 0) then
         call ErrorHandler('readMomentDirectionData','iostat > 0',ios)
      endif
!
!     ================================================================
!     Determine the format of the Evec data
!     ================================================================
      DataForm = -1
      isDirect = .false.
      ElementName = .false.
      RotateEvec = .false.
      StartingLine = 0
      NB = 0
      nsn = 0
      nst = 0
      RotationMatrix(1:3,1:3) = ZERO
      RotationMatrix(1,1) = ONE
      RotationMatrix(2,2) = ONE
      RotationMatrix(3,3) = ONE
      ConstrField(:,:) = ZERO
      emix(:) = -ONE
!     ----------------------------------------------------------------
      call initString(80)
!     ----------------------------------------------------------------
      do
         read(funit,'(a)',iostat=status)text
         if (status < 0) then
            call ErrorHandler('readMomentDirectionData','Invalid file',fname)
         else
            StartingLine = StartingLine + 1
         endif
         text=adjustl(text)
         if (text(1:1) == '#' .or. text(1:1) == '!'                   &
                              .or. len_trim(text) == 0) then
            cycle
         else
!           ----------------------------------------------------------
            call setString(text)
!           ----------------------------------------------------------
            n = getNumTokens()
!           ----------------------------------------------------------
            if (n == 3) then   ! 3-column data
               NB = NB + 1
               if (NB > 3) then
                  call ErrorHandler('readMomentDirectionData',        &
                                    'Invalid line',text)
               else
                  read(text,*,iostat=status)RotationMatrix(1:3,NB)
                  if (status > 0) then
                     call ErrorHandler('readMomentDirectionData',     &
                                       'Invalid line',text)
                  endif
               endif
               RotateEvec = .true.
            else
               if (n == 4) then   ! 4-column data
                  DataForm = 0
!                 ----------------------------------------------------
                  call readToken(1,atom,alen)
!                 ----------------------------------------------------
               else if (n == 7) then   ! 7-column data
                  DataForm = 1
!                 ----------------------------------------------------
                  call readToken(1,atom,alen)
!                 ----------------------------------------------------
               else if (n == 8) then
                  DataForm = 2        ! with mixing parameter added to each line
!                 ----------------------------------------------------
                  call readToken(1,atom,alen)
!                 ----------------------------------------------------
               else if (n > 8) then
                  DataForm = 3        ! info_evec format
!                 ----------------------------------------------------
                  call readToken(1,atom,alen)
!                 ----------------------------------------------------
                  RotateEvec = .false.
               else
!                 ----------------------------------------------------
                  call ErrorHandler('readMomentDirectionData',        &
                                    'Invalid data form',text)
!                 ----------------------------------------------------
               endif
               if (isNumber(atom)) then
                  ElementName = .false.
               else
                  ElementName = .true.
               endif
               StartingLine = StartingLine - 1
               exit
            endif
         endif
      enddo
!     ----------------------------------------------------------------
      call endString()
!     ----------------------------------------------------------------
!
      if (RotateEvec) then
         if (NB < 3) then
!           ----------------------------------------------------------
            call ErrorHandler('readMomentDirectionData',              &
                              'Rotation Matrix is ill defined',fname)
!           ----------------------------------------------------------
         else if ( abs(abs(determ3(RotationMatrix))-ONE) > TEN2m8 ) then
!           ----------------------------------------------------------
            call ErrorHandler('readMomentDirectionData',              &
                              'Rotation Matrix is not unitary')
!           ----------------------------------------------------------
         endif
      else if (StartingLine < NB) then
!        -------------------------------------------------------------
         call ErrorHandler('readMomentDirectionData',                 &
                           'Invalid Evec Data',fname)
!        -------------------------------------------------------------
      endif
!
      rewind(funit)
!
!     ================================================================
!     StartingLine+1 is the line number of the first Evec data.
!     ================================================================
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
         if (text(1:1) == '#' .or. text(1:1) == '!'                   &
                              .or. len_trim(text) == 0) then
            cycle
         else
            n = n+1
            if (DataForm == 0) then
                if (ElementName) then
                   j = getTokenPosition(1,text,k)
!                  k = getTokenPosition(2,text)
!                  atom = trim(text(j:k-1))
                   atom = text(j:j+k-1)
                   if ( AtomicNumber(n) /= getZtot(atom(1:3)) ) then
!                     ------------------------------------------------
                      call ErrorHandler('readMomentDirectionData',    &
                         'Incompatible Evec Data and Position Data for index',n)
!                     ------------------------------------------------
                   endif
                   read(text(j+k:),*) AtomEvec(1:3,n)
                else
                   read(text,*) anum, AtomEvec(1:3,n)
                   if ( AtomicNumber(n) /= anum ) then
!                     ------------------------------------------------
                      call ErrorHandler('readMomentDirectionData',    &
                         'Incompatible Evec Data and Position Data for index',n)
!                     ------------------------------------------------
                   endif
                endif
            else if (DataForm == 1) then
                if (ElementName) then
                   j = getTokenPosition(1,text,k)
!                  k = getTokenPosition(2,text)
!                  atom = trim(text(j:k-1))
                   atom = text(j:j+k-1)
                   if ( AtomicNumber(n) /= getZtot(atom(1:3)) ) then
!                     ------------------------------------------------
                      call ErrorHandler('readMomentDirectionData',    &
                         'Incompatible Evec Data and Position Data for index',n)
!                     ------------------------------------------------
                   endif
                   read(text(j+k:),*) AtomEvec(1:3,n), ConstrField(1:3,n)
                else
                   read(text,*)anum, AtomEvec(1:3,n), ConstrField(1:3,n)
                   if ( AtomicNumber(n) /= anum ) then
!                     ------------------------------------------------
                      call ErrorHandler('readMomentDirectionData',    &
                         'Incompatible Evec Data and Position Data for index',n)
!                     ------------------------------------------------
                   endif
                endif
            else if (DataForm == 2) then
                if (ElementName) then
                   j = getTokenPosition(1,text,k)
!                  k = getTokenPosition(2,text)
!                  atom = trim(text(j:k-1))
                   atom = text(j:j+k-1)
                   if ( AtomicNumber(n) /= getZtot(atom(1:3)) ) then
!                     ------------------------------------------------
                      call ErrorHandler('readMomentDirectionData',    &
                         'Incompatible Evec Data and Position Data for index',n)
!                     ------------------------------------------------
                   endif
                   read(text(j+k:),*) AtomEvec(1:3,n), ConstrField(1:3,n), emix(n)
                else
                   read(text,*) anum, AtomEvec(1:3,n), &
                                ConstrField(1:3,n), emix(n)
                   if ( AtomicNumber(n) /= anum ) then
!                     ------------------------------------------------
                      call ErrorHandler('readMomentDirectionData',    &
                         'Incompatible Evec Data and Position Data for index',n)
!                     ------------------------------------------------
                   endif
                endif
            else if (DataForm == 3) then
                if (ElementName) then
                   j = getTokenPosition(1,text,k)
!                  k = getTokenPosition(2,text)
!                  atom = trim(text(j:k-1))
                   atom = text(j:j+k-1)
                   if ( AtomicNumber(n) /= getZtot(atom(1:3)) ) then
!                     ------------------------------------------------
                      call ErrorHandler('readMomentDirectionData',    &
                         'Incompatible Evec Data and Position Data for index',n)
!                     ------------------------------------------------
                   endif
                   k = getTokenPosition(3,text)
                   read(text(j+k:),*)x,y,z,                             &
                          AtomEvec(1:3,n), emix(n), ConstrField(1:2,n)
                   ConstrField(3,n) = ZERO
                else
                   read(text,*)anum,k,x,y,z,                          &
                          AtomEvec(1:3,n), emix(n), ConstrField(1:2,n)
                   ConstrField(3,n) = ZERO
                   if ( AtomicNumber(n) /= anum ) then
!                     ------------------------------------------------
                      call ErrorHandler('readMomentDirectionData',            &
                        'Incompatible Evec Data and Position Data for index',n)
!                     ------------------------------------------------
                   endif
                endif
            else
!               ------------------------------------------------------
                call ErrorHandler('readMomentDirectionData',          &
                                  'Unkown data format')
!               ------------------------------------------------------
            endif
         endif
      enddo
!
      if (n < num_atoms) then
!        -------------------------------------------------------------
         call ErrorHandler('readMomentDirectionData','Evec data is incomplete')
!        -------------------------------------------------------------
      endif
!
      close(funit)
      write(6,'(/,a,/)')' The Evec Data is read successfully!'
!
      if (RotateEvec) then
         do n = 1, num_atoms
            x = AtomEvec(1,n)
            y = AtomEvec(2,n)
            z = AtomEvec(3,n)
            AtomEvec(1,n) = RotationMatrix(1,1)*x +                   &
                            RotationMatrix(1,2)*y +                   &
                            RotationMatrix(1,3)*z
            AtomEvec(2,n) = RotationMatrix(2,1)*x +                   &
                            RotationMatrix(2,2)*y +                   &
                            RotationMatrix(2,3)*z
            AtomEvec(3,n) = RotationMatrix(3,1)*x +                   &
                            RotationMatrix(3,2)*y +                   &
                            RotationMatrix(3,3)*z
!
            x = ConstrField(1,n)
            y = ConstrField(2,n)
            z = ConstrField(3,n)
            ConstrField(1,n) = RotationMatrix(1,1)*x +                &
                               RotationMatrix(1,2)*y +                &
                               RotationMatrix(1,3)*z
            ConstrField(2,n) = RotationMatrix(2,1)*x +                &
                               RotationMatrix(2,2)*y +                &
                               RotationMatrix(2,3)*z
            ConstrField(3,n) = RotationMatrix(3,1)*x +                &
                               RotationMatrix(3,2)*y +                &
                               RotationMatrix(3,3)*z
         enddo
      endif
   endif
!
!  -------------------------------------------------------------------
   call bcastMessage(p_AtomEvec,3*num_atoms,0)
   call bcastMessage(p_ConstrField,3*num_atoms,0)
   call bcastMessage(emix,num_atoms,0)
!  -------------------------------------------------------------------
!
   nullify( AtomicNumber )
   nullify( AtomEvec, ConstrField, emix )
   nullify( p_AtomEvec, p_ConstrField )
!
   end subroutine readMomentDirectionData
