subroutine setEnvVar(env_name,env_value,set_status)
    use iso_c_binding
    use KindParamModule, only : IntKind
    implicit none
!
    character (len=*), intent(in) :: env_name
    character (len=*), intent(in) :: env_value
!
    integer (kind=intKind), intent(out) :: set_status
!
    interface
        function c_setenv(name, value, overwrite) bind(C, name="setenv")
            import :: c_char, c_int
            character(kind=c_char), dimension(*) :: name
            character(kind=c_char), dimension(*) :: value
            integer(c_int), value :: overwrite
            integer(c_int) :: c_setenv
        end function c_setenv
    end interface
!
    integer(c_int) :: ierr
!
!   ------------------------------------------------------------------
    ierr = c_setenv(trim(env_name)//C_NULL_CHAR,                      &
                    trim(env_value)//C_NULL_CHAR,                     &
                    1_c_int)
!   ------------------------------------------------------------------
!   if (ierr == 0_c_int) then
!       print *, "MY_VAR set successfully."
!   else
!       print *, "Failed to set MY_VAR."
!   end if
!
!
    set_status = int(ierr,kind=IntKind)
!
end subroutine setEnvVar
