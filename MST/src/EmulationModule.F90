module EmulationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
!
   use MathParamModule, only : ZERO
!
public :: initEmulation,          &
          endEmulation,           &
          isEmulationEnabled,     &
          setEmulationPrec,       &
          getPrecEnvVar,          &
          getPrecEnvVarValue,     &
          getExtraEnvVar,         &
          getExtraEnvVarValue
!
   integer (kind=IntKind), parameter, public :: VarStrLen = 32
!
private
   integer (kind=IntKind), parameter :: MaxSchemes = 8 ! Maximum number of emulation schemes
   integer (kind=IntKind), parameter :: MaxEnvVars = 8 ! Maximum number of environment
                                                       ! variables for an emulation scheme
!
   integer (kind=IntKind) :: NumSchemes                ! The number of emulation schemes defined
   integer (kind=IntKind) :: NumEnvVars(MaxSchemes)    ! The number of defined environment variables
                                                       ! for each given emulation scheme
!
!  ===================================================================
!  Options defined in CmdLineOptions.inc
!  ===================================================================
   character (len=VarStrLen), target :: SchemeName(MaxSchemes)
   character (len=VarStrLen), target :: PrecEnvVar(MaxSchemes)
   character (len=VarStrLen), target :: PrecEnvVarValue(MaxSchemes)
   character (len=VarStrLen), target :: ExtraEnvVar(MaxEnvVars,MaxSchemes)
   character (len=VarStrLen), target :: ExtraEnvVarValue(MaxEnvVars,MaxSchemes)
!
   integer (kind=IntKind) :: SchemeIndex = 0  ! This equals the actual emulation scheme applied
   integer (kind=IntKind) :: print_instr = -1
!
   integer (kind=IntKind), parameter :: MaxEnergyZones = 8 ! Maximum number of contour zones
   integer (kind=IntKind) :: num_energy_zones              ! Actual number of contour zones
   character (len=8) ::      prec_ez_param(MaxEnergyZones) ! The emulation precision parameter for
                                                           ! each contour zone, starting from the
                                                           ! most left to the most right
!  ===================================================================
!  energy_zone_param are defined as follows:
!      energy_zone_param(1): energy contour zone 1, for which Re[energy] starts from the bottom
!                            of the valence band up to energy_zone_param(1)
!      energy_zone_param(2): energy contour zone 2, which is right next to zone 1 and ends where
!                            Im[energy] = energy_zone_param(2)
!      energy_zone_param(3): energy contour zone 3, which is right next to zone 2 and ends where
!                            Im[energy] = energy_zone_param(3)
!      ...,...
!      energy_zone_param(N): energy contour zone N, which is right next to zone N-1 and ends where
!                            Im[energy] = energy_zone_param(N) = 0.0. Here N = num_energy_zones.
!  The constrain: energy_zone_param(2) > energy_zone_param(3) > ... > energy_zone_param(N) = 0.0
!  ===================================================================
   real (kind=RealKind) ::   energy_zone_param(MaxEnergyZones)
!
   logical :: emul_enabled = .false.
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initEmulation(iprint)
!  ===================================================================
   use CmdLineOptionModule, only : getCmdLineOption
!
   implicit none
!
   character (len=VarStrLen) :: emul_scheme, eval
!
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: n, i, j, nn
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   print_instr = iprint
   SchemeIndex = 0
!
   NumEnvVars = 0
   PrecEnvVar = ' '
   PrecEnvVarValue = ' '
   ExtraEnvVar = ' '
   ExtraEnvVarValue = ' '
!
   n = 0
!
   include '../src/EmulationEnvVar.inc'
!
   NumSchemes = n
!
   do i = 1, n
      NumEnvVars(i) = 0
      do j = 1, MaxEnvVars
         if (len_trim(ExtraEnvVar(j,i)) > 0) then
            if (j == NumEnvVars(i) + 1) then
               NumEnvVars(i) = NumEnvVars(i) + 1
            else
               call ErrorHandler('initEmulation','An error is found in EmulationEnvVar.inc for the scheme',i)
            endif
         endif
      enddo
   enddo
!
!  ===================================================================
!  Determine if emulation is enabled in the job environment setup
!  ===================================================================
   emul_enabled = .false.
   LOOP_i: do i = 1, NumSchemes
!     ----------------------------------------------------------------
      call get_environment_variable(PrecEnvVar(i),PrecEnvVarValue(i), &
                                    status=n)
!     ----------------------------------------------------------------
      if (n == 0) then
         SchemeIndex = i
         if (print_instr >= 0) then
            write(6,'(/,a)') '========================================'
            write(6,'(2a)')  'Emulation Scheme : ',SchemeName(i)
            write(6,'(a)')   'Initial Emulation Envrionment Variables:'
            write(6,'(a)')   '----------------------------------------'
            write(6,'(4a)')'*  ',trim(PrecEnvVar(i)),' = ',trim(PrecEnvVarValue(i))
         endif
         do j = 1, NumEnvVars(i)
!           ----------------------------------------------------------
            call get_environment_variable(ExtraEnvVar(j,i),           &
                                          ExtraEnvVarValue(j,i),status=nn)
!           ----------------------------------------------------------
            if (print_instr >= 0) then
               if (nn == 0) then
                  write(6,'(4a)')'*  ',trim(ExtraEnvVar(j,i)),' = ',trim(ExtraEnvVarValue(j,i))
               else
                  write(6,'(3a)')'*  ',trim(ExtraEnvVar(j,i)),' is undefined!'
               endif
            endif
         enddo
         if (print_instr >= 0) then
            write(6,'(a)')   '========================================'
         endif
         num_energy_zones = 1
         emul_enabled = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
!
!  ===================================================================
!  Obtain the emulation instructions from the command line.
!  -------------------------------------------------------------------
   call read_cmdline_emul_param()
!  -------------------------------------------------------------------
!
   end subroutine initEmulation
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endEmulation()
!  ===================================================================
   implicit none
!
   NumEnvVars = 0
   PrecEnvVar = ' '
   PrecEnvVarValue = ' '
   ExtraEnvVar = ' '
   ExtraEnvVarValue = ' '
   emul_enabled = .false.
   SchemeIndex = 0
!
   end subroutine endEmulation
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isEmulationEnabled() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   y = emul_enabled
!
   end function isEmulationEnabled
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setEmulationPrec(en,eval)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), intent(in), optional :: en
!
   character (len=*), intent(in), optional :: eval
!
   real (kind=RealKind) :: er, ei
!
   integer (kind=IntKind) :: i, zone_id, set_status
!
   if (.not.emul_enabled .or. num_energy_zones < 1) then
      if (print_instr >= 0) then
         write(6,'(a,2d15.8)')'Emulation is not performed at this energy',en
      endif
   else if (num_energy_zones == 1) then
      if (print_instr >= 0) then
         write(6,'(a)')'Emulation is performed uniformly for all energies'
         write(6,'(a,a)')'---- Emulation precision is set to: ',trim(PrecEnvVarValue(SchemeIndex))
      endif
   else if (present(eval)) then
      PrecEnvVarValue(SchemeIndex) = eval
!     ----------------------------------------------------------------
      call setEnvVar(PrecEnvVar(SchemeIndex),eval,                       &
                     set_status)
!     ----------------------------------------------------------------
      if (set_status /= 0) then
         call WarningHandler('setEmulationPrec','Setting precision failed', &
                             PrecEnvVar(SchemeIndex),eval)
      else if (print_instr >= 0) then
         write(6,'(a,a)')'In setEmulationPrec, emulation precision is set to: ',trim(eval) 
      endif
   else if (present(en)) then
      er = real(en,kind=RealKind)
      ei = aimag(en)
      zone_id = 0
      if (er <= energy_zone_param(1)) then
         zone_id = 1
      else
         do i = 2, num_energy_zones
            if (ei > energy_zone_param(i)) then
               zone_id = i
               exit
            endif
         enddo
      endif
      if (zone_id > 0) then
         PrecEnvVarValue(SchemeIndex) = prec_ez_param(zone_id)
      endif
!     ----------------------------------------------------------------
      call setEnvVar(PrecEnvVar(SchemeIndex),                            &
                     PrecEnvVarValue(SchemeIndex),                       &
                     set_status)
!     ----------------------------------------------------------------
      if (set_status /= 0) then
         write(6,'(a,2d15.8)')'At energy = ',en
         call WarningHandler('setEmulationPrec','Setting precision failed', &
                             PrecEnvVar(SchemeIndex))
      else if (print_instr >= 0) then
         write(6,'(a,2d15.8)')'In setEmulationPrec, at energy = ',en
         write(6,'(a,a)')     '---- Emulation precision is set to: ',trim(PrecEnvVarValue(SchemeIndex))
      endif
   endif
!
   end subroutine setEmulationPrec
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPrecEnvVar() result(evar)
!  ===================================================================
   implicit none
!
   character (len=VarStrLen) :: evar
!
   if (emul_enabled) then
      evar = PrecEnvVar(SchemeIndex)
   else
      evar = ' '
   endif
!
   end function getPrecEnvVar
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPrecEnvVarValue() result(eval)
!  ===================================================================
   implicit none
!
   character (len=VarStrLen) :: eval
!
   if (emul_enabled) then
      eval = PrecEnvVarValue(SchemeIndex)
   else
      eval = ' '
   endif
!
   end function getPrecEnvVarValue
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExtraEnvVar(num_evs) result(evar)
!  ===================================================================
   implicit none
!
   character (len=VarStrLen), pointer :: evar(:)
!
   integer (kind=IntKind), intent(out) :: num_evs
!
   if (emul_enabled) then
      evar => ExtraEnvVar(:,SchemeIndex)
      num_evs = NumEnvVars(SchemeIndex)
   else
      num_evs = 0
      nullify(evar)
   endif
!
   end function getExtraEnvVar
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExtraEnvVarValue(id_var) result(eval)
!  ===================================================================
   implicit none
!
   character (len=VarStrLen) :: eval
!
   integer (kind=IntKind), intent(in) :: id_var
   integer (kind=IntKind) :: i, j
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (emul_enabled) then
      if (id_var < 1 .or. id_var >  NumEnvVars(SchemeIndex)) then
         call ErrorHandler('getExtraEnvVarValue',               &
                           'The environment variable index is out of range',id_var)
      else
         eval = ExtraEnvVarValue(id_var)
      endif
   else
      eval = ' '
   endif
!
   end function getExtraEnvVarValue
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine read_cmdline_emul_param()
!  ===================================================================
   use CmdLineOptionModule, only : initCmdLineOption,  &
                                   endCmdLineOption,   &
                                   getCmdLineOption,   &
                                   printCmdLineOption, &
                                   isOptionDefined,    &
                                   getCmdLineOptionValue
!
   use StringModule, only : initString, endString, setString, getNumTokens, readToken
!
   implicit none
!
   integer (kind=intKind), parameter :: string_len = 32
!
   character (len=string_len) :: sval, zval, s
   character (len=VarStrLen) :: sname
!
   integer (kind=IntKind) :: i, num_zone_params, rstatus, set_status
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
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
      function isNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isNumber
   end interface
!
   interface
      subroutine setEnvVar(env_name,env_value,set_status)
         use KindParamModule, only : IntKind
         implicit none
!
         character (len=*), intent(in) :: env_name
         character (len=*), intent(in) :: env_value
!
         integer (kind=intKind), intent(out) :: set_status
      end subroutine setEnvVar
   end interface
!
   num_zone_params = 0
   prec_ez_param = ' '
   energy_zone_param = ZERO
!
!  ===================================================================
!  Check if the scheme defined on the command line is consistent with
!  the job environment setup.
!  ===================================================================
   if (isOptionDefined('Emulation Scheme')) then
      rstatus = getCmdLineOptionValue('Emulation Scheme',sname)
      if (emul_enabled) then
         if (nocaseCompare(sname,SchemeName(SchemeIndex)) /= 0) then
            call ErrorHandler('read_cmdline_emul_param','Inconsistent schemes are given', &
                              sname,SchemeName(SchemeIndex))
         endif
      else
!        =============================================================
!        In this case, the environment variables for emulation are not
!        set before the job starts. They will be set here.
!        =============================================================
         SchemeIndex = 0
         do i = 1, NumSchemes
            if (nocaseCompare(sname,SchemeName(i)) == 0) then
               SchemeIndex = i
               exit
            endif
         enddo
         if (SchemeIndex == 0) then
!           ----------------------------------------------------------
            call WarningHandler('read_cmdline_emul_param',            &
                                'The given emulation scheme is invalid',sname)
!           ----------------------------------------------------------
            emul_enabled = .false.
         else
!           ==========================================================
!           set the emulation environment variables here .............
!           ----------------------------------------------------------
            call setEnvVar(PrecEnvVar(SchemeIndex),                   &
                           PrecEnvVarValue(SchemeIndex),              &
                           set_status)
!           ----------------------------------------------------------
            if (set_status == 0) then
               emul_enabled = .true.
               do i = 1, NumEnvVars(SchemeIndex)
!                 ----------------------------------------------------
                  call setEnvVar(ExtraEnvVar(i,SchemeIndex),             &
                                 ExtraEnvVarValue(i,SchemeIndex),        &
                                 set_status)
!                 ----------------------------------------------------
                  if (set_status /= 0) then
!                    -------------------------------------------------
                     call WarningHandler('read_cmdline_emul_param',      &
                                   'Setting envrionment variable failed',&
                                   ExtraEnvVar(i,SchemeIndex))
!                    -------------------------------------------------
                     emul_enabled = .false.
                     exit
                  endif
               enddo
            else
!              -------------------------------------------------------
               call WarningHandler('read_cmdline_emul_param',            &
                                   'Setting envrionment variable failed',&
                                   PrecEnvVar(SchemeIndex))
!              -------------------------------------------------------
               emul_enabled = .false.
            endif
         endif
      endif
   endif
!
   if (.not.emul_enabled) then
      SchemeIndex = 0
      return
   endif
!
   if (.not.isOptionDefined('Emulation Precision Parameters')) then
      return
   endif
!
!  -------------------------------------------------------------------
   call initString(string_len)
!  -------------------------------------------------------------------
!
   if (getCmdLineOptionValue('Emulation Precision Parameters',sval) == 0) then
!     ----------------------------------------------------------------
      call setString(sval)
      num_energy_zones = getNumTokens()
!     ----------------------------------------------------------------
      if (num_energy_zones > MaxEnergyZones) then
         call ErrorHandler('read_cmdline_emul_param',                 &
                           'The number of zones is beyond the internal size',num_energy_zones)
      endif
      do i = 1, num_energy_zones
!        -------------------------------------------------------------
         call readToken(i,s)
!        -------------------------------------------------------------
         prec_ez_param(i) = s
         if (print_instr >= 0) then
            write(6,'(a,i2,2a)')'Scheme ',i,': ',prec_ez_param(i)
         endif
         if (isInteger(s) .and. print_instr >= 0) then
            write(6,'(a)')'This precision parameter is an integer'
         endif
      enddo
   else
      num_energy_zones = 1
   endif
!
   if (num_energy_zones > 1) then
      energy_zone_param = ZERO
      if (getCmdLineOptionValue('Energy Contour Zone Parameters',zval) == 0) then
!        -------------------------------------------------------------
         call setString(zval)
         num_zone_params = getNumTokens()
!        -------------------------------------------------------------
         if (num_zone_params+1 == num_energy_zones) then
            do i = 1, num_zone_params
!              -------------------------------------------------------
               call readToken(i,s)
!              -------------------------------------------------------
               if (isNumber(s)) then
                  read(s,*)energy_zone_param(i)
                  if (print_instr >= 0) then
                     write(6,'(a,i2,a,f8.5)')'Zone ',i,': ',energy_zone_param(i)
                  endif
                  if (i > 1 .and. energy_zone_param(i) < ZERO) then
                     if (print_instr >= 0) then
                        write(6,'(a,a)')'ERROR: The contour zone parameter < 0.0: ',zval
                     endif
                     emul_enabled = .false.
                     exit
                  endif
               else
                  if (print_instr >= 0) then
                     write(6,'(a,a)')'Invalid format of the contour zone parameters: ',zval
                  endif
                  emul_enabled = .false.
                  exit
               endif
            enddo
            do i = 2, num_zone_params
               if (energy_zone_param(i) < energy_zone_param(i+1)) then
                  if (print_instr >= 0) then
                     write(6,'(a,f8.5,a,f8.5)')'Invalid order of contour zone parameter: ', &
                                               energy_zone_param(i),' <',energy_zone_param(i+1)
                     write(6,'(a,a)')'Invalid format of the contour zone parameters: ',zval
                  endif
                  emul_enabled = .false.
                  exit
               endif
            enddo
         else
            if (print_instr >= 0) then
               write(6,'(2(a,i4))')'The number of contour zone parameters, ',num_zone_params, &
                                   ', is inconsistent with the expected number of zones,',num_energy_zones
            endif
            emul_enabled = .false.
         endif
      endif
   endif
!  -------------------------------------------------------------------
   call endString()
!  -------------------------------------------------------------------
!
   if (.not.emul_enabled) then
!     ================================================================
!     if the contour zone parameters are not given properly, the 
!     emulation is disabled
!     ================================================================
      SchemeIndex = 0
      num_energy_zones = 0
      prec_ez_param = ' '
      energy_zone_param = ZERO
   endif
!
end subroutine read_cmdline_emul_param
!=====================================================================
!
end module EmulationModule
