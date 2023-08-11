!! Copyright (C) 2016 Micael Oliveira
!!               2020 Susi Lehtola
!! All rights reserved.
!!
!! This Source Code Form is subject to the terms of the Mozilla Public
!! License, v. 2.0. If a copy of the MPL was not distributed with this
!! file, You can obtain one at mozilla.org/MPL/2.0/.

module xc_f03_lib_m
  use, intrinsic :: iso_c_binding
  implicit none

  private
  public :: &
    ! version
    xc_f03_version, &
    xc_f03_version_string, &
    ! literature reference
    xc_f03_reference, &
    xc_f03_reference_doi, &
    ! func_info
    xc_f03_func_info_t, &
    xc_f03_func_info_get_number, &
    xc_f03_func_info_get_kind, &
    xc_f03_func_info_get_name, &
    xc_f03_func_info_get_family, &
    xc_f03_func_info_get_references, &
    xc_f03_func_info_get_flags, &
    xc_f03_func_info_get_n_ext_params, &
    xc_f03_func_info_get_ext_params_name, &
    xc_f03_func_info_get_ext_params_description, &
    xc_f03_func_info_get_ext_params_default_value, &
    ! func_reference
    xc_f03_func_reference_t, &
    xc_f03_func_reference_get_ref, &
    xc_f03_func_reference_get_doi, &
    xc_f03_func_reference_get_bibtex, &
    ! func
    xc_f03_func_t, &
    xc_f03_func_init, &
    xc_f03_func_end, &
    xc_f03_func_get_info, &
    xc_f03_functional_get_name, &
    xc_f03_functional_get_number, &
    xc_f03_family_from_id, &
    xc_f03_number_of_functionals, &
    xc_f03_maximum_name_length, &
    xc_f03_available_functional_numbers, &
    xc_f03_available_functional_names, &
    xc_f03_func_set_dens_threshold, &
    xc_f03_func_set_zeta_threshold, &
    xc_f03_func_set_sigma_threshold, &
    xc_f03_func_set_tau_threshold, &
    xc_f03_func_set_ext_params, &
    xc_f03_func_set_ext_params_name, &
    ! lda
    xc_f03_lda, &
    xc_f03_lda_exc, &
    xc_f03_lda_exc_vxc, &
    xc_f03_lda_exc_vxc_fxc, &
    xc_f03_lda_exc_vxc_fxc_kxc, &
    xc_f03_lda_vxc, &
    xc_f03_lda_vxc_fxc, &
    xc_f03_lda_vxc_fxc_kxc, &
    xc_f03_lda_fxc, &
    xc_f03_lda_kxc, &
    xc_f03_lda_lxc, &
    ! gga
    xc_f03_gga, &
    xc_f03_gga_exc, &
    xc_f03_gga_exc_vxc, &
    xc_f03_gga_exc_vxc_fxc, &
    xc_f03_gga_exc_vxc_fxc_kxc, &
    xc_f03_gga_vxc, &
    xc_f03_gga_vxc_fxc, &
    xc_f03_gga_vxc_fxc_kxc, &
    xc_f03_gga_fxc, &
    xc_f03_gga_kxc, &
    xc_f03_gga_lxc, &
    xc_f03_gga_ak13_get_asymptotic, &
    xc_f03_hyb_exx_coef, &
    xc_f03_hyb_cam_coef, &
    xc_f03_nlc_coef, &
    ! mgga
    xc_f03_mgga, &
    xc_f03_mgga_exc, &
    xc_f03_mgga_exc_vxc, &
    xc_f03_mgga_exc_vxc_fxc, &
    xc_f03_mgga_exc_vxc_fxc_kxc, &
    xc_f03_mgga_vxc, &
    xc_f03_mgga_vxc_fxc, &
    xc_f03_mgga_vxc_fxc_kxc, &
    xc_f03_mgga_fxc, &
    xc_f03_mgga_kxc, &
    xc_f03_mgga_lxc

  integer(c_int), parameter, public :: &
    XC_UNPOLARIZED          =   1,     &  ! Spin unpolarized
    XC_POLARIZED            =   2         ! Spin polarized

    integer(c_int), parameter, public :: &
    XC_NON_RELATIVISTIC     =   0,     &  ! Functional includes or not relativistic
    XC_RELATIVISTIC         =   1         ! corrections. Only available in some functionals.

  ! Kinds
  integer(c_int), parameter, public :: &
    XC_EXCHANGE             =   0,     &
    XC_CORRELATION          =   1,     &
    XC_EXCHANGE_CORRELATION =   2,     &
    XC_KINETIC              =   3

  ! Families of xc functionals
  integer(c_int), parameter, public :: &
    XC_FAMILY_UNKNOWN       =  -1,     &
    XC_FAMILY_NONE          =   0,     &
    XC_FAMILY_LDA           =   1,     &
    XC_FAMILY_GGA           =   2,     &
    XC_FAMILY_MGGA          =   4,     &
    XC_FAMILY_LCA           =   8,     &
    XC_FAMILY_OEP           =  16,     &
    XC_FAMILY_HYB_GGA       =  32,     &
    XC_FAMILY_HYB_MGGA      =  64,     &
    XC_FAMILY_HYB_LDA       = 128

  integer(c_int), parameter, public :: &
    XC_FLAGS_HAVE_EXC        =     1,   &
    XC_FLAGS_HAVE_VXC        =     2,   &
    XC_FLAGS_HAVE_FXC        =     4,   &
    XC_FLAGS_HAVE_KXC        =     8,   &
    XC_FLAGS_HAVE_LXC        =    16,   &
    XC_FLAGS_HAVE_ALL        =    31,   & ! The most common case
    XC_FLAGS_1D              =    32,   &
    XC_FLAGS_2D              =    64,   &
    XC_FLAGS_3D              =   128,   &
    XC_FLAGS_HYB_CAM         =   256,   &
    XC_FLAGS_HYB_CAMY        =   512,   &
    XC_FLAGS_VV10            =  1024,   &
    XC_FLAGS_HYB_LC          =  2048,   &
    XC_FLAGS_HYB_LCY         =  4096,   &
    XC_FLAGS_STABLE          =  8192,   &
    XC_FLAGS_DEVELOPMENT     = 16384,   &
    XC_FLAGS_NEEDS_LAPLACIAN = 32768

  integer(c_int), parameter, public :: &
    XC_TAU_EXPLICIT         =     0,   &
    XC_TAU_EXPANSION        =     1

  integer(c_int), parameter, public :: &
    XC_MAX_REFERENCES       =     5

  ! List of functionals
#include "libxc_inc.f90"

  ! These are old names kept for compatibility
  integer(c_int), parameter, public :: &
    XC_LDA_X_1D             =  21,     &
    XC_GGA_X_BGCP           =  38,     &
    XC_GGA_C_BGCP           =  39,     &
    XC_GGA_C_BCGP           =  39,     &
    XC_GGA_C_VPBE           =  83,     &
    XC_GGA_XC_LB            = 160,     &
    XC_MGGA_C_CC06          = 229,     &
    XC_GGA_K_ABSR1          = 506,     &
    XC_GGA_K_ABSR2          = 507,     &
    XC_LDA_C_LP_A           = 547,     &
    XC_LDA_C_LP_B           = 548,     &
    XC_MGGA_C_LP90          = 564

  !----------------------------------------------------------------
  interface
    subroutine xc_version(major, minor, micro) bind(c)
      import
      integer(c_int), intent(out) :: major, minor, micro
    end subroutine xc_version

    type(c_ptr) function xc_version_string() bind(c)
      import
    end function xc_version_string

    type(c_ptr) function xc_reference() bind(c)
      import
    end function xc_reference

    type(c_ptr) function xc_reference_doi() bind(c)
      import
    end function xc_reference_doi
  end interface


  !----------------------------------------------------------------
  type :: xc_f03_func_info_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_info_t

  interface
    integer(c_int) function xc_func_info_get_number(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_number

    integer(c_int) function xc_func_info_get_kind(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_kind

    type(c_ptr) function xc_func_info_get_name(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_name

    integer(c_int) function xc_func_info_get_family(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_family

    integer(c_int) function xc_func_info_get_flags(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_flags

    type(c_ptr) function xc_func_info_get_references(info, number) bind(c)
      import
      type(c_ptr),    value :: info
      integer(c_int), value :: number
    end function xc_func_info_get_references

    integer(c_int) function xc_func_info_get_n_ext_params(info) bind(c)
      import
      type(c_ptr), value :: info
    end function xc_func_info_get_n_ext_params

    type(c_ptr) function xc_func_info_get_ext_params_name(info, number) bind(c)
      import
      type(c_ptr),    value :: info
      integer(c_int), value :: number
    end function xc_func_info_get_ext_params_name

    type(c_ptr) function xc_func_info_get_ext_params_description(info, number) bind(c)
      import
      type(c_ptr),    value :: info
      integer(c_int), value :: number
    end function xc_func_info_get_ext_params_description

    real(c_double) function xc_func_info_get_ext_params_default_value(info, number) bind(c)
      import
      type(c_ptr),    value :: info
      integer(c_int), value :: number
    end function xc_func_info_get_ext_params_default_value

  end interface

  !----------------------------------------------------------------
  type :: xc_f03_func_reference_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_reference_t

  interface
    type(c_ptr) function xc_func_reference_get_ref(reference) bind(c)
      import
      type(c_ptr), value :: reference
    end function xc_func_reference_get_ref

    type(c_ptr) function xc_func_reference_get_doi(reference) bind(c)
      import
      type(c_ptr), value :: reference
    end function xc_func_reference_get_doi

    type(c_ptr) function xc_func_reference_get_bibtex(reference) bind(c)
      import
      type(c_ptr), value :: reference
    end function xc_func_reference_get_bibtex
  end interface

  !----------------------------------------------------------------
  type :: xc_f03_func_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_t

  interface
    type(c_ptr) function xc_func_alloc() bind(c)
      import
    end function xc_func_alloc

    integer(c_int) function xc_func_init(p, functional, nspin) bind(c)
      import
      type(c_ptr),    value :: p
      integer(c_int), value :: functional, nspin
    end function xc_func_init

    subroutine xc_func_end(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine xc_func_end

    subroutine xc_func_free(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine xc_func_free

    subroutine libxc_free(p) bind(c)
      import
      type(c_ptr), value :: p
    end subroutine libxc_free

    type(c_ptr) function xc_func_get_info(p) bind(c)
      import
      type(c_ptr), value :: p
    end function xc_func_get_info

    type(c_ptr) function xc_functional_get_name(number) bind(c)
      import
      integer(c_int), value :: number
    end function xc_functional_get_name

    integer(c_int) function xc_functional_get_number(func_string) bind(c)
      import
      character(kind=c_char), intent(in) :: func_string(*)
    end function xc_functional_get_number

    integer(c_int) function xc_family_from_id(id, family, number) bind(c)
      import
      integer(c_int), value :: id
      type(c_ptr),    value :: family, number
    end function xc_family_from_id

    integer(c_int) function xc_f03_number_of_functionals() bind(c, name="xc_number_of_functionals")
      import
    end function xc_f03_number_of_functionals

    integer(c_int) function xc_f03_maximum_name_length() bind(c, name="xc_maximum_name_length")
      import
    end function xc_f03_maximum_name_length

    subroutine xc_f03_available_functional_numbers(list) bind(c, name="xc_available_functional_numbers")
      import
      integer(c_int), intent(out) :: list(*)
    end subroutine xc_f03_available_functional_numbers

    subroutine xc_available_functional_names(list) bind(c)
      import
      type(c_ptr) :: list(*)
    end subroutine xc_available_functional_names

    subroutine xc_func_set_dens_threshold(p, dens_threshold) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: dens_threshold
    end subroutine xc_func_set_dens_threshold

    subroutine xc_func_set_zeta_threshold(p, zeta_threshold) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: zeta_threshold
    end subroutine xc_func_set_zeta_threshold

    subroutine xc_func_set_sigma_threshold(p, sigma_threshold) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: sigma_threshold
    end subroutine xc_func_set_sigma_threshold

    subroutine xc_func_set_tau_threshold(p, tau_threshold) bind(c)
      import
      type(c_ptr), value :: p
      real(c_double), value :: tau_threshold
    end subroutine xc_func_set_tau_threshold

    subroutine xc_func_set_ext_params(p, ext_params) bind(c)
      import
      type(c_ptr), value      :: p
      real(c_double), intent(in) :: ext_params(*)
    end subroutine xc_func_set_ext_params

    subroutine xc_func_set_ext_params_name(p, name, par) bind(c)
      import
      type(c_ptr), value      :: p
      character(kind=c_char), intent(in) :: name(*)
      real(c_double), value   :: par
    end subroutine xc_func_set_ext_params_name
end interface
    
  ! LDAs
  !----------------------------------------------------------------
  interface
    subroutine xc_lda(p, np, rho, zk, vrho, v2rho2, v3rho3, v4rho4) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*), v4rho4(*)
    end subroutine xc_lda

    subroutine xc_lda_exc(p, np, rho, zk) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: zk(*)
    end subroutine xc_lda_exc

    subroutine xc_lda_exc_vxc(p, np, rho, zk, vrho) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: zk(*), vrho(*)
    end subroutine xc_lda_exc_vxc

    subroutine xc_lda_exc_vxc_fxc(p, np, rho, zk, vrho, v2rho2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), v2rho2(*)
    end subroutine xc_lda_exc_vxc_fxc

    subroutine xc_lda_exc_vxc_fxc_kxc(p, np, rho, zk, vrho, v2rho2, v3rho3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)
    end subroutine xc_lda_exc_vxc_fxc_kxc

    subroutine xc_lda_vxc(p, np, rho, vrho) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: vrho(*)
    end subroutine xc_lda_vxc

    subroutine xc_lda_vxc_fxc(p, np, rho, vrho, v2rho2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: vrho(*), v2rho2(*)
    end subroutine xc_lda_vxc_fxc

    subroutine xc_lda_vxc_fxc_kxc(p, np, rho, vrho, v2rho2, v3rho3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: vrho(*), v2rho2(*), v3rho3(*)
    end subroutine xc_lda_vxc_fxc_kxc

    subroutine xc_lda_fxc(p, np, rho, v2rho2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: v2rho2(*)
    end subroutine xc_lda_fxc

    subroutine xc_lda_kxc(p, np, rho, v3rho3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: v3rho3(*)
    end subroutine xc_lda_kxc

    subroutine xc_lda_lxc(p, np, rho, v4rho4) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*)
      real(c_double),    intent(out) :: v4rho4(*)
    end subroutine xc_lda_lxc
  end interface


  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_gga(p, np, rho, sigma, zk, vrho, vsigma,        &
         v2rho2, v2rhosigma, v2sigma2,                            &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,              &
         v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 &
         ) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
      real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
      real(c_double),    intent(out) :: v4rho4(*), v4rho3sigma(*), v4rho2sigma2(*), v4rhosigma3(*), v4sigma4(*)
    end subroutine xc_gga

    subroutine xc_gga_exc(p, np, rho, sigma, zk) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: zk(*)
    end subroutine xc_gga_exc

    subroutine xc_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*)
    end subroutine xc_gga_exc_vxc

    subroutine xc_gga_exc_vxc_fxc(p, np, rho, sigma, zk, vrho, vsigma,        &
         v2rho2, v2rhosigma, v2sigma2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    end subroutine xc_gga_exc_vxc_fxc

    subroutine xc_gga_exc_vxc_fxc_kxc(p, np, rho, sigma, zk, vrho, vsigma,        &
         v2rho2, v2rhosigma, v2sigma2,                            &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
      real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    end subroutine xc_gga_exc_vxc_fxc_kxc

    subroutine xc_gga_vxc(p, np, rho, sigma, vrho, vsigma) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: vrho(*), vsigma(*)
    end subroutine xc_gga_vxc

    subroutine xc_gga_vxc_fxc(p, np, rho, sigma, vrho, vsigma,    &
         v2rho2, v2rhosigma, v2sigma2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: vrho(*), vsigma(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    end subroutine xc_gga_vxc_fxc

    subroutine xc_gga_vxc_fxc_kxc(p, np, rho, sigma, vrho, vsigma, &
         v2rho2, v2rhosigma, v2sigma2,                             &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: vrho(*), vsigma(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
      real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    end subroutine xc_gga_vxc_fxc_kxc
    
    subroutine xc_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    end subroutine xc_gga_fxc

    subroutine xc_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    end subroutine xc_gga_kxc

    subroutine xc_gga_lxc(p, np, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*)
      real(c_double),    intent(out) :: v4rho4(*), v4rho3sigma(*), v4rho2sigma2(*), v4rhosigma3(*), v4sigma4(*)
    end subroutine xc_gga_lxc

  end interface


  interface
    real(c_double) function xc_gga_ak13_get_asymptotic(homo) bind(c)
      import
      real(c_double), value :: homo
    end function xc_gga_ak13_get_asymptotic
  end interface


  interface
    real(c_double) function xc_hyb_exx_coef(p) bind(c)
      import
      type(c_ptr), value :: p
    end function xc_hyb_exx_coef

    subroutine xc_hyb_cam_coef(p, omega, alpha, beta) bind(c)
      import
      type(c_ptr), value       :: p
      real(c_double), intent(out) :: omega, alpha, beta
    end subroutine xc_hyb_cam_coef

    subroutine xc_nlc_coef(p, nlc_b, nlc_c) bind(c)
      import
      type(c_ptr), value       :: p
      real(c_double), intent(out) :: nlc_b, nlc_c
    end subroutine xc_nlc_coef
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,     &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3,                                                            &
         v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
         v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
         v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
         v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
         v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
         v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
         v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
         ) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
           v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
      real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
           v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),                  &
           v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),         &
           v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),       &
           v3lapltau2(*), v3tau3(*)
      real(c_double),    intent(out) :: &
           v4rho4(*), v4rho3sigma(*), v4rho3lapl(*), v4rho3tau(*), v4rho2sigma2(*), v4rho2sigmalapl(*),     &
           v4rho2sigmatau(*), v4rho2lapl2(*), v4rho2lapltau(*), v4rho2tau2(*), v4rhosigma3(*),              &
           v4rhosigma2lapl(*), v4rhosigma2tau(*), v4rhosigmalapl2(*), v4rhosigmalapltau(*),                 &
           v4rhosigmatau2(*), v4rholapl3(*), v4rholapl2tau(*), v4rholapltau2(*), v4rhotau3(*), v4sigma4(*), &
           v4sigma3lapl(*), v4sigma3tau(*), v4sigma2lapl2(*), v4sigma2lapltau(*), v4sigma2tau2(*),          &
           v4sigmalapl3(*), v4sigmalapl2tau(*), v4sigmalapltau2(*), v4sigmatau3(*), v4lapl4(*),             &
           v4lapl3tau(*), v4lapl2tau2(*), v4lapltau3(*), v4tau4(*)
    end subroutine xc_mgga

    subroutine xc_mgga_exc(p, np, rho, sigma, lapl, tau, zk) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: zk(*)
    end subroutine xc_mgga_exc

    subroutine xc_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    end subroutine xc_mgga_exc_vxc

    subroutine xc_mgga_exc_vxc_fxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
           v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    end subroutine xc_mgga_exc_vxc_fxc

    subroutine xc_mgga_exc_vxc_fxc_kxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,     &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
           v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
      real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
           v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),                  &
           v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),         &
           v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),       &
           v3lapltau2(*), v3tau3(*)
    end subroutine xc_mgga_exc_vxc_fxc_kxc

    subroutine xc_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    end subroutine xc_mgga_vxc

    subroutine xc_mgga_vxc_fxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau, &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
           v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    end subroutine xc_mgga_vxc_fxc

    subroutine xc_mgga_vxc_fxc_kxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,     &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
           v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
      real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
           v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),                  &
           v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),         &
           v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),       &
           v3lapltau2(*), v3tau3(*)
    end subroutine xc_mgga_vxc_fxc_kxc

    subroutine xc_mgga_fxc(p, np, rho, sigma, lapl, tau, &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
           v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    end subroutine xc_mgga_fxc

    subroutine xc_mgga_kxc(p, np, rho, sigma, lapl, tau, &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
           v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
           v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
           v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
           v3lapltau2(*), v3tau3(*)
    end subroutine xc_mgga_kxc

    subroutine xc_mgga_lxc(p, np, rho, sigma, lapl, tau, &
         v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
         v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
         v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
         v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
         v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
         v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
         v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
         ) bind(c)
      import
      type(c_ptr),       value       :: p
      integer(c_size_t), value       :: np
      real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
      real(c_double),    intent(out) :: &
           v4rho4(*), v4rho3sigma(*), v4rho3lapl(*), v4rho3tau(*), v4rho2sigma2(*), v4rho2sigmalapl(*),     &
           v4rho2sigmatau(*), v4rho2lapl2(*), v4rho2lapltau(*), v4rho2tau2(*), v4rhosigma3(*),              &
           v4rhosigma2lapl(*), v4rhosigma2tau(*), v4rhosigmalapl2(*), v4rhosigmalapltau(*),                 &
           v4rhosigmatau2(*), v4rholapl3(*), v4rholapl2tau(*), v4rholapltau2(*), v4rhotau3(*), v4sigma4(*), &
           v4sigma3lapl(*), v4sigma3tau(*), v4sigma2lapl2(*), v4sigma2lapltau(*), v4sigma2tau2(*),          &
           v4sigmalapl3(*), v4sigmalapl2tau(*), v4sigmalapltau2(*), v4sigmatau3(*), v4lapl4(*),             &
           v4lapl3tau(*), v4lapl2tau2(*), v4lapltau3(*), v4tau4(*)
    end subroutine xc_mgga_lxc
  end interface
    
  contains

  !----------------------------------------------------------------
  subroutine xc_f03_version(major, minor, micro)
    integer(c_int), intent(out) :: major, minor, micro

    call xc_version(major, minor, micro)

  end subroutine xc_f03_version

  subroutine xc_f03_version_string(version)
    character(len=*), intent(out) :: version

    type(c_ptr) :: c_version

    c_version = xc_version_string()
    call c_to_f_string_ptr(c_version, version)

  end subroutine xc_f03_version_string

  subroutine xc_f03_reference(ref)
    character(len=*), intent(out) :: ref

    type(c_ptr) :: c_ref

    c_ref = xc_reference()
    call c_to_f_string_ptr(c_ref, ref)

  end subroutine xc_f03_reference

  subroutine xc_f03_reference_doi(doi)
    character(len=*), intent(out) :: doi

    type(c_ptr) :: c_doi

    c_doi = xc_reference_doi()
    call c_to_f_string_ptr(c_doi, doi)

  end subroutine xc_f03_reference_doi

  !----------------------------------------------------------------
  integer(c_int) function xc_f03_func_info_get_number(info) result(number)
    type(xc_f03_func_info_t), intent(in) :: info

    number = xc_func_info_get_number(info%ptr)

  end function xc_f03_func_info_get_number

  integer(c_int) function xc_f03_func_info_get_kind(info) result(kind)
    type(xc_f03_func_info_t), intent(in) :: info

    kind = xc_func_info_get_kind(info%ptr)

  end function xc_f03_func_info_get_kind

  character(len=128) function xc_f03_func_info_get_name(info) result(name)
    type(xc_f03_func_info_t), intent(in) :: info

    call c_to_f_string_ptr(xc_func_info_get_name(info%ptr), name)

  end function xc_f03_func_info_get_name

  integer(c_int) function xc_f03_func_info_get_family(info) result(family)
    type(xc_f03_func_info_t), intent(in) :: info

    family = xc_func_info_get_family(info%ptr)

  end function xc_f03_func_info_get_family

  integer(c_int) function xc_f03_func_info_get_flags(info) result(flags)
    type(xc_f03_func_info_t), intent(in) :: info

    flags = xc_func_info_get_flags(info%ptr)

  end function xc_f03_func_info_get_flags

  type(xc_f03_func_reference_t) function xc_f03_func_info_get_references(info, number) result(reference)
    type(xc_f03_func_info_t), intent(in)    :: info
    integer(c_int),           intent(inout) :: number ! number of the reference. Must be 0 in the first call

    type(c_ptr) :: next_ref

    reference%ptr = xc_func_info_get_references(info%ptr, number)
    next_ref = xc_func_info_get_references(info%ptr, INT(number + 1, c_int))
    if (c_associated(next_ref)) then
      number = number + 1
    else
      number = -1
    end if

  end function xc_f03_func_info_get_references

  integer(c_int) function xc_f03_func_info_get_n_ext_params(info) result(n_ext_params)
    type(xc_f03_func_info_t), intent(in) :: info

    n_ext_params = xc_func_info_get_n_ext_params(info%ptr)

  end function xc_f03_func_info_get_n_ext_params

  character(len=128) function xc_f03_func_info_get_ext_params_name(info, number) result(name)
    type(xc_f03_func_info_t), intent(in) :: info
    integer(c_int),           intent(in) :: number

    call c_to_f_string_ptr(xc_func_info_get_ext_params_name(info%ptr, number), name)

  end function xc_f03_func_info_get_ext_params_name

  character(len=128) function xc_f03_func_info_get_ext_params_description(info, number) result(description)
    type(xc_f03_func_info_t), intent(in) :: info
    integer(c_int),           intent(in) :: number

    call c_to_f_string_ptr(xc_func_info_get_ext_params_description(info%ptr, number), description)

  end function xc_f03_func_info_get_ext_params_description

  real(c_double) function xc_f03_func_info_get_ext_params_default_value(info, number) result(val)
    type(xc_f03_func_info_t), intent(in) :: info
    integer(c_int),           intent(in) :: number

    val = xc_func_info_get_ext_params_default_value(info%ptr, number)

  end function xc_f03_func_info_get_ext_params_default_value

  !----------------------------------------------------------------
  character(len=1024) function xc_f03_func_reference_get_ref(reference) result(ref)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_ref(reference%ptr), ref)

  end function xc_f03_func_reference_get_ref

  character(len=1024) function xc_f03_func_reference_get_doi(reference) result(doi)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_doi(reference%ptr), doi)

  end function xc_f03_func_reference_get_doi

  character(len=1024) function xc_f03_func_reference_get_bibtex(reference) result(bibtex)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_bibtex(reference%ptr), bibtex)

  end function xc_f03_func_reference_get_bibtex


  !----------------------------------------------------------------
  subroutine xc_f03_func_init(p, functional, nspin, err)
    type(xc_f03_func_t),      intent(inout) :: p
    integer(c_int),           intent(in)    :: functional
    integer(c_int),           intent(in)    :: nspin
    integer(c_int), optional, intent(out)   :: err

    integer(c_int) :: ierr

    p%ptr = xc_func_alloc()
    ierr = xc_func_init(p%ptr, functional, nspin)

    if(present(err)) err = ierr
  end subroutine xc_f03_func_init

  subroutine xc_f03_func_end(p)
    type(xc_f03_func_t), intent(inout) :: p

    call xc_func_end(p%ptr)
    call xc_func_free(p%ptr)

  end subroutine xc_f03_func_end

  type(xc_f03_func_info_t) function xc_f03_func_get_info(p) result(info)
    type(xc_f03_func_t), intent(in) :: p

    info%ptr = xc_func_get_info(p%ptr)

  end function xc_f03_func_get_info

  character(len=128) function xc_f03_functional_get_name(number) result(name)
    integer(c_int), intent(in) :: number
    type(c_ptr) :: cstr

    cstr = xc_functional_get_name(number)
    call c_to_f_string_ptr(cstr, name)
    call libxc_free(cstr)

  end function xc_f03_functional_get_name

  integer(c_int) function xc_f03_functional_get_number(func_string) result(number)
    character(len=*), intent(in) :: func_string

    number = xc_functional_get_number(f_to_c_string(func_string))

  end function xc_f03_functional_get_number

  integer(c_int) function xc_f03_family_from_id(id, family, number)
    integer(c_int), intent(in)                    :: id
    integer(c_int), intent(out), optional, target :: family, number

    type(c_ptr) c_family, c_number
    integer(c_int), pointer :: f_family, f_number

    if (present(family)) then
      f_family => family
      call c_f_pointer(c_family, f_family)
    else
      c_family = C_NULL_PTR
    end if
    if (present(number)) then
      f_number => number
      call c_f_pointer(c_number, f_number)
    else
      c_number = C_NULL_PTR
    end if

    xc_f03_family_from_id = xc_family_from_id(id, c_family, c_number)

  end function xc_f03_family_from_id

  subroutine xc_f03_available_functional_names(list)
    character(len=*), intent(out) :: list(*)

    integer(c_int) :: n, i, maxlen
    character(kind=c_char), allocatable, target :: names(:,:)
    type(c_ptr), allocatable :: c_list(:)

    n = xc_f03_number_of_functionals()
    maxlen = xc_f03_maximum_name_length()

    allocate(names(maxlen, n))
    allocate(c_list(n))
    do i = 1, n
      c_list(i) = c_loc(names(1,i))
    end do

    call xc_available_functional_names(c_list)

    do i = 1, n
      call c_to_f_string_ptr(c_list(i), list(i))
    end do

    deallocate(c_list)
    deallocate(names)

  end subroutine xc_f03_available_functional_names


  subroutine xc_f03_func_set_dens_threshold(p, dens_threshold)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double),      intent(in) :: dens_threshold

    call xc_func_set_dens_threshold(p%ptr, dens_threshold)

  end subroutine xc_f03_func_set_dens_threshold

  subroutine xc_f03_func_set_zeta_threshold(p, zeta_threshold)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double),      intent(in) :: zeta_threshold

    call xc_func_set_zeta_threshold(p%ptr, zeta_threshold)

  end subroutine xc_f03_func_set_zeta_threshold

  subroutine xc_f03_func_set_sigma_threshold(p, sigma_threshold)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double),      intent(in) :: sigma_threshold

    call xc_func_set_sigma_threshold(p%ptr, sigma_threshold)

  end subroutine xc_f03_func_set_sigma_threshold

  subroutine xc_f03_func_set_tau_threshold(p, tau_threshold)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double),      intent(in) :: tau_threshold

    call xc_func_set_tau_threshold(p%ptr, tau_threshold)

  end subroutine xc_f03_func_set_tau_threshold

  subroutine xc_f03_func_set_ext_params(p, ext_params)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double),      intent(in) :: ext_params(*)

    call xc_func_set_ext_params(p%ptr, ext_params)

  end subroutine xc_f03_func_set_ext_params

  subroutine xc_f03_func_set_ext_params_name(p, name, par)
    type(xc_f03_func_t), intent(in) :: p
    character(len=*), intent(in) :: name
    real(c_double),      intent(in) :: par

    call xc_func_set_ext_params_name(p%ptr, f_to_c_string(name), par)

  end subroutine xc_f03_func_set_ext_params_name

  ! LDAs
  !----------------------------------------------------------------
  subroutine xc_f03_lda(p, np, rho, zk, vrho, v2rho2, v3rho3, v4rho4)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*), v4rho4(*)

    call xc_lda(p%ptr, np, rho, zk, vrho, v2rho2, v3rho3, v4rho4)

  end subroutine xc_f03_lda

  subroutine xc_f03_lda_exc(p, np, rho, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*)

    call xc_lda_exc(p%ptr, np, rho, zk)

  end subroutine xc_f03_lda_exc

  subroutine xc_f03_lda_exc_vxc(p, np, rho, zk, vrho)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*), vrho(*)

    call xc_lda_exc_vxc(p%ptr, np, rho, zk, vrho)

  end subroutine xc_f03_lda_exc_vxc

  subroutine xc_f03_lda_exc_vxc_fxc(p, np, rho, zk, vrho, v2rho2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), v2rho2(*)

    call xc_lda_exc_vxc_fxc(p%ptr, np, rho, zk, vrho, v2rho2)

  end subroutine xc_f03_lda_exc_vxc_fxc

  subroutine xc_f03_lda_exc_vxc_fxc_kxc(p, np, rho, zk, vrho, v2rho2, v3rho3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)

    call xc_lda_exc_vxc_fxc_kxc(p%ptr, np, rho, zk, vrho, v2rho2, v3rho3)

  end subroutine xc_f03_lda_exc_vxc_fxc_kxc

  subroutine xc_f03_lda_vxc(p, np, rho, vrho)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: vrho(*)

    call xc_lda_vxc(p%ptr, np, rho, vrho)

  end subroutine xc_f03_lda_vxc

  subroutine xc_f03_lda_vxc_fxc(p, np, rho, vrho, v2rho2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: vrho(*), v2rho2(*)

    call xc_lda_vxc_fxc(p%ptr, np, rho, vrho, v2rho2)

  end subroutine xc_f03_lda_vxc_fxc

  subroutine xc_f03_lda_vxc_fxc_kxc(p, np, rho, vrho, v2rho2, v3rho3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: vrho(*), v2rho2(*), v3rho3(*)

    call xc_lda_vxc_fxc_kxc(p%ptr, np, rho, vrho, v2rho2, v3rho3)

  end subroutine xc_f03_lda_vxc_fxc_kxc

  subroutine xc_f03_lda_fxc(p, np, rho, v2rho2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: v2rho2(*)

    call xc_lda_fxc(p%ptr, np, rho, v2rho2)

  end subroutine xc_f03_lda_fxc

  subroutine xc_f03_lda_kxc(p, np, rho, v3rho3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: v3rho3(*)

    call xc_lda_kxc(p%ptr, np, rho, v3rho3)

  end subroutine xc_f03_lda_kxc

  subroutine xc_f03_lda_lxc(p, np, rho, v4rho4)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: v4rho4(*)

    call xc_lda_lxc(p%ptr, np, rho, v4rho4)

  end subroutine xc_f03_lda_lxc

  ! GGAs
  !----------------------------------------------------------------
  subroutine xc_f03_gga(p, np, rho, sigma, zk, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2,                            &
       v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,              &
       v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 &
    )
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    real(c_double),      intent(out) :: v4rho4(*), v4rho3sigma(*), v4rho2sigma2(*), v4rhosigma3(*), v4sigma4(*)
    
    call xc_gga(p%ptr, np, rho, sigma, zk, vrho, vsigma,          &
         v2rho2, v2rhosigma, v2sigma2,                            &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,              &
         v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 &
         )

  end subroutine xc_f03_gga

  subroutine xc_f03_gga_exc(p, np, rho, sigma, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*)

    call xc_gga_exc(p%ptr, np, rho, sigma, zk)

  end subroutine xc_f03_gga_exc

  subroutine xc_f03_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*)

    call xc_gga_exc_vxc(p%ptr, np, rho, sigma, zk, vrho, vsigma)

  end subroutine xc_f03_gga_exc_vxc

  subroutine xc_f03_gga_exc_vxc_fxc(p, np, rho, sigma, zk, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    
    call xc_gga_exc_vxc_fxc(p%ptr, np, rho, sigma, zk, vrho, vsigma,          &
         v2rho2, v2rhosigma, v2sigma2)

  end subroutine xc_f03_gga_exc_vxc_fxc

  subroutine xc_f03_gga_exc_vxc_fxc_kxc(p, np, rho, sigma, zk, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2,                            &
       v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    
    call xc_gga_exc_vxc_fxc_kxc(p%ptr, np, rho, sigma, zk, vrho, vsigma,          &
         v2rho2, v2rhosigma, v2sigma2,                            &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)

  end subroutine xc_f03_gga_exc_vxc_fxc_kxc
  
  subroutine xc_f03_gga_vxc(p, np, rho, sigma, vrho, vsigma)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*)

    call xc_gga_vxc(p%ptr, np, rho, sigma, vrho, vsigma)

  end subroutine xc_f03_gga_vxc

  subroutine xc_f03_gga_vxc_fxc(p, np, rho, sigma, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    
    call xc_gga_vxc_fxc(p%ptr, np, rho, sigma, vrho, vsigma, &
         v2rho2, v2rhosigma, v2sigma2)

  end subroutine xc_f03_gga_vxc_fxc

  subroutine xc_f03_gga_vxc_fxc_kxc(p, np, rho, sigma, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2,                            &
       v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    
    call xc_gga_vxc_fxc_kxc(p%ptr, np, rho, sigma, vrho, vsigma,  &
         v2rho2, v2rhosigma, v2sigma2,                            &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)

  end subroutine xc_f03_gga_vxc_fxc_kxc

  subroutine xc_f03_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)

    call xc_gga_fxc(p%ptr, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)

  end subroutine xc_f03_gga_fxc

  subroutine xc_f03_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)

    call xc_gga_kxc(p%ptr, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)

  end subroutine xc_f03_gga_kxc

  subroutine xc_f03_gga_lxc(p, np, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: v4rho4(*), v4rho3sigma(*), v4rho2sigma2(*), v4rhosigma3(*), v4sigma4(*)

    call xc_gga_lxc(p%ptr, np, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4)
  end subroutine xc_f03_gga_lxc
  
  real(c_double) function xc_f03_gga_ak13_get_asymptotic(homo) result(asymptotic)
    real(c_double), intent(in) :: homo

    asymptotic = xc_gga_ak13_get_asymptotic(homo)

  end function xc_f03_gga_ak13_get_asymptotic

  real(c_double) function xc_f03_hyb_exx_coef(p) result(coef)
    type(xc_f03_func_t), intent(in) :: p

    coef = xc_hyb_exx_coef(p%ptr)

  end function xc_f03_hyb_exx_coef

  subroutine xc_f03_hyb_cam_coef(p, omega, alpha, beta)
    type(xc_f03_func_t), intent(in)  :: p
    real(c_double),       intent(out) :: omega, alpha, beta

    call xc_hyb_cam_coef(p%ptr, omega, alpha, beta)

  end subroutine xc_f03_hyb_cam_coef

  subroutine xc_f03_nlc_coef(p, nlc_b, nlc_c)
    type(xc_f03_func_t), intent(in)  :: p
    real(c_double),       intent(out) :: nlc_b, nlc_c

    call xc_nlc_coef(p%ptr, nlc_b, nlc_c)

  end subroutine xc_f03_nlc_coef


  ! the meta-GGAs
  !----------------------------------------------------------------
  subroutine xc_f03_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,   &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3,                                                            &
         v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
         v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
         v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
         v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
         v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
         v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
         v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
         )
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
         v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
         v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
         v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
         v3lapltau2(*), v3tau3(*)
    real(c_double),    intent(out) :: &
           v4rho4(*), v4rho3sigma(*), v4rho3lapl(*), v4rho3tau(*), v4rho2sigma2(*), v4rho2sigmalapl(*),     &
           v4rho2sigmatau(*), v4rho2lapl2(*), v4rho2lapltau(*), v4rho2tau2(*), v4rhosigma3(*),              &
           v4rhosigma2lapl(*), v4rhosigma2tau(*), v4rhosigmalapl2(*), v4rhosigmalapltau(*),                 &
           v4rhosigmatau2(*), v4rholapl3(*), v4rholapl2tau(*), v4rholapltau2(*), v4rhotau3(*), v4sigma4(*), &
           v4sigma3lapl(*), v4sigma3tau(*), v4sigma2lapl2(*), v4sigma2lapltau(*), v4sigma2tau2(*),          &
           v4sigmalapl3(*), v4sigmalapl2tau(*), v4sigmalapltau2(*), v4sigmatau3(*), v4lapl4(*),             &
           v4lapl3tau(*), v4lapl2tau2(*), v4lapltau3(*), v4tau4(*)

    call xc_mgga(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3,                                                            &
         v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
         v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
         v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
         v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
         v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
         v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
         v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
         )
  end subroutine xc_f03_mgga

  subroutine xc_f03_mgga_exc(p, np, rho, sigma, lapl, tau, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*)

    call xc_mgga_exc(p%ptr, np, rho, sigma, lapl, tau, zk)

  end subroutine xc_f03_mgga_exc

  subroutine xc_f03_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)

    call xc_mgga_exc_vxc(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)

  end subroutine xc_f03_mgga_exc_vxc

  subroutine xc_f03_mgga_exc_vxc_fxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,   &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)

    call xc_mgga_exc_vxc_fxc(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2)
  end subroutine xc_f03_mgga_exc_vxc_fxc

  subroutine xc_f03_mgga_exc_vxc_fxc_kxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,   &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
         v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
         v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
         v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
         v3lapltau2(*), v3tau3(*)

    call xc_mgga_exc_vxc_fxc_kxc(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3)
  end subroutine xc_f03_mgga_exc_vxc_fxc_kxc

  subroutine xc_f03_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)

    call xc_mgga_vxc(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)

  end subroutine xc_f03_mgga_vxc

  subroutine xc_f03_mgga_vxc_fxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,   &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)

    call xc_mgga_vxc_fxc(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2)
  end subroutine xc_f03_mgga_vxc_fxc

  subroutine xc_f03_mgga_vxc_fxc_kxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,   &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
         v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
         v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
         v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
         v3lapltau2(*), v3tau3(*)

    call xc_mgga_vxc_fxc_kxc(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3)
  end subroutine xc_f03_mgga_vxc_fxc_kxc

  subroutine xc_f03_mgga_fxc(p, np, rho, sigma, lapl, tau, &
       v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
       v2lapl2, v2lapltau, v2tau2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)

    call xc_mgga_fxc(p%ptr, np, rho, sigma, lapl, tau,   &
      v2rho2, v2rhosigma, v2rholapl, v2rhotau,           &
      v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2)

  end subroutine xc_f03_mgga_fxc

  subroutine xc_f03_mgga_kxc(p, np, rho, sigma, lapl, tau, &
       v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
       v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
       v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
       v3lapltau2, v3tau3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
           v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
           v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
           v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
           v3lapltau2(*), v3tau3(*)

    call xc_mgga_kxc(p%ptr, np, rho, sigma, lapl, tau,  &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3)
  end subroutine xc_f03_mgga_kxc

  subroutine xc_f03_mgga_lxc(p, np, rho, sigma, lapl, tau, &
       v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
       v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
       v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
       v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
       v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
       v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
       v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
       )
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: &
         v4rho4(*), v4rho3sigma(*), v4rho3lapl(*), v4rho3tau(*), v4rho2sigma2(*), v4rho2sigmalapl(*),     &
         v4rho2sigmatau(*), v4rho2lapl2(*), v4rho2lapltau(*), v4rho2tau2(*), v4rhosigma3(*),              &
         v4rhosigma2lapl(*), v4rhosigma2tau(*), v4rhosigmalapl2(*), v4rhosigmalapltau(*),                 &
         v4rhosigmatau2(*), v4rholapl3(*), v4rholapl2tau(*), v4rholapltau2(*), v4rhotau3(*), v4sigma4(*), &
         v4sigma3lapl(*), v4sigma3tau(*), v4sigma2lapl2(*), v4sigma2lapltau(*), v4sigma2tau2(*),          &
         v4sigmalapl3(*), v4sigmalapl2tau(*), v4sigmalapltau2(*), v4sigmatau3(*), v4lapl4(*),             &
         v4lapl3tau(*), v4lapl2tau2(*), v4lapltau3(*), v4tau4(*)

    call xc_mgga_lxc(p%ptr, np, rho, sigma, lapl, tau,  &
         v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
         v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
         v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
         v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
         v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
         v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
         v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
         )
  end subroutine xc_f03_mgga_lxc


  ! Helper functions to convert between C and Fortran strings
  ! Based on the routines by Joseph M. Krahn
  function f_to_c_string(f_string) result(c_string)
    character(len=*), intent(in) :: f_string
    character(kind=c_char,len=1) :: c_string(len_trim(f_string)+1)

    integer :: i, strlen

    strlen = len_trim(f_string)

    forall (i=1:strlen)
      c_string(i) = f_string(i:i)
    end forall
    c_string(strlen+1) = C_NULL_CHAR

  end function f_to_c_string

  subroutine c_to_f_string(c_string, f_string)
    character(kind=c_char,len=1), intent(in)  :: c_string(*)
    character(len=*),             intent(out) :: f_string

    integer :: i

    i = 1
    do while(c_string(i) /= C_NULL_CHAR .and. i <= len(f_string))
      f_string(i:i) = c_string(i)
      i = i + 1
    end do
    if (i < len(f_string)) f_string(i:) = ' '

  end subroutine c_to_f_string

  subroutine c_to_f_string_ptr(c_string, f_string)
    type(c_ptr),      intent(in)  :: c_string
    character(len=*), intent(out) :: f_string

    character(len=1, kind=c_char), pointer :: p_chars(:)
    integer :: i

    if (.not. c_associated(c_string)) then
      f_string = ' '
    else
      call c_f_pointer(c_string, p_chars, [huge(0)])
      i = 1
      do while(p_chars(i) /= C_NULL_CHAR .and. i <= len(f_string))
        f_string(i:i) = p_chars(i)
        i = i + 1
      end do
      if (i < len(f_string)) f_string(i:) = ' '
    end if

  end subroutine c_to_f_string_ptr

end module xc_f03_lib_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
