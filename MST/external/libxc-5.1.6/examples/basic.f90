program lxctest
  use xc_f90_lib_m

  implicit none

  TYPE(xc_f90_func_t) :: xc_func
  TYPE(xc_f90_func_info_t) :: xc_info
  real(8) :: rho(5) = (/0.1, 0.2, 0.3, 0.4, 0.5/)
  real(8) :: sigma(5) = (/0.2, 0.3, 0.4, 0.5, 0.6/)
  real(8) :: exc(5)
  integer :: i, vmajor, vminor, vmicro, func_id = 1

  call xc_f90_version(vmajor, vminor, vmicro)
  write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro

  call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
  xc_info = xc_f90_func_get_info(xc_func)

  select case (xc_f90_func_info_get_family(xc_info))
  case(XC_FAMILY_LDA)
    call xc_f90_lda_exc(xc_func, 5_8, rho(1), exc(1))
  case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
    call xc_f90_gga_exc(xc_func, 5_8, rho(1), sigma(1), exc(1))
  end select

  do i = 1, 5
    write(*,"(F8.6,1X,F9.6)") rho(i), exc(i)
  end do

  call xc_f90_func_end(xc_func)

end program lxctest
