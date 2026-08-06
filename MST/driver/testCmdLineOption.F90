program testCmdlineOption
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use CmdLineOptionModule, only : initCmdLineOption,  &
                                   endCmdLineOption,   &
                                   getCmdLineOption,   &
                                   printCmdLineOption, &
                                   getCmdLineOptionValue
   use StringModule, only : initString, endString, setString, getNumTokens, readToken
   implicit none
!
   integer (kind=IntKind) :: num_args, i, n, rstatus, num_params, num_zones
   integer (kind=IntKind) :: scheme_params(10)
!
   integer (kind=intKind), parameter :: string_len = 32
!
   character (len=string_len) :: sval, zval, evar, s
!
   real (kind=RealKind) :: zone_params(10)
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
   call initCmdLineOption(num_args)
!
   write(6,'(a,i5)')'Number of command line options: ',num_args
!
   do i = 1, num_args
      call printCmdLineOption(i)
   enddo
!
   rstatus = getCmdLineOptionValue('Emulation Scheme',n) 
   if (rstatus == 0 .and. n > 0) then
      call initString(string_len)
!
      scheme_params = 0
      if (getCmdLineOptionValue('Emulation Scheme Parameters',sval) == 0) then
         call setString(sval)
         num_params = getNumTokens()
         do i = 1, num_params
            call readToken(i,s)
            if (isInteger(s)) then
               read(s,*)scheme_params(i)
               write(6,'(a,i2,a,i5)')'Scheme ',i,': ',scheme_params(i)
            else
               write(6,'(a,a)')'Invalid format of the scheme parameters: ',sval
            endif
         enddo
      endif
!
      zone_params = 0.0d0
      if (getCmdLineOptionValue('Contour Zone Parameters',zval) == 0) then
         call setString(zval)
         num_zones = getNumTokens() + 1
         do i = 1, num_zones-1
            call readToken(i,s)
            if (isNumber(s)) then
               read(s,*)zone_params(i)
               write(6,'(a,i2,a,f8.5)')'Zone ',i,': ',zone_params(i)
            else
               write(6,'(a,a)')'Invalid format of the contour zone parameters: ',zval
            endif
         enddo
      else
         num_zones = 1
      endif
!
!     ----------------------------------------------------------------
      call setEnvVar('EMULATION','Ozaki-1',n)
!     ----------------------------------------------------------------
      if (n == 0) then
!        -------------------------------------------------------------
         call get_environment_variable('EMULATION', evar)
!        -------------------------------------------------------------
         write(6,'(a,a)')'The environment variable is set: EMULATION = ',trim(evar)
      else
         write(6,'(a)')'Setting the environment variable has failed.'
      endif
!
      call endString()
   endif
!
   call endCmdLineOption()
!
end program testCmdlineOption
