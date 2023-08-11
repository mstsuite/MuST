cmake_policy(SET CMP0048 NEW)  # project_VERSION* variables populated from project(... VERSION x.x.x) string
project(Libxc
  LANGUAGES C)
set(Libxc_AUTHORS      "Miguel A.L. Marques and others")
set(Libxc_DESCRIPTION  "Exchange-correlation functionals for density-functional theory")
set(Libxc_EMAIL        "libxc@tddft.org")
set(Libxc_URL          "http://www.tddft.org/programs/Libxc")
set(Libxc_LICENSE      "Mozilla Public License, version 2.0 (MPL-2.0)")

cmake_minimum_required(VERSION 3.1)
set(CMAKE_C_STANDARD 99)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

################################### Options ####################################
include(psi4OptionsTools)
option_with_default(CMAKE_BUILD_TYPE "Build type" Release)
option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static. For Windows, also add CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=ON" OFF)
option_with_print(BUILD_TESTING "Compile the testing infrastructure" ON)
option_with_default(BUILD_FPIC "Libraries will be compiled with position independent code" ON)
if((${BUILD_SHARED_LIBS}) AND NOT ${BUILD_FPIC})
  message(FATAL_ERROR "BUILD_SHARED_LIBS ON and BUILD_FPIC OFF are incompatible, as shared library requires position independent code")
endif()
option_with_default(NAMESPACE_INSTALL_INCLUDEDIR "Location within CMAKE_INSTALL_INCLUDEDIR to which headers are installed. Psi4 wants namespace (e.g., /libxc)" /)
option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)
option_with_print(ENABLE_FORTRAN "Build Fortran 2003 interface" OFF)
option_with_print(ENABLE_PYTHON "Install Python API interface" OFF)
if(${ENABLE_PYTHON} AND NOT ${BUILD_SHARED_LIBS})
  message(FATAL_ERROR "ENABLE_PYTHON ON requires BUILD_SHARED_LIBS ON because only shared libraries can be dynamically loaded")
endif()

# GCC appears to be catching -xHost as "-x Host" which breaks the detection
if(CMAKE_C_COMPILER_ID MATCHES Intel)
  option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON "-xHost" "-march=native" "/arch:AVX2")
else()
  option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON  "-march=native" "-xHost" "/arch:AVX2")
endif()

option_with_print(DISABLE_VXC "Don't compile first derivative code" OFF)
option_with_print(DISABLE_FXC "Don't compile second derivative code" OFF)
option_with_print(DISABLE_KXC "Don't compile third derivative code" ON)
option_with_print(DISABLE_LXC "Don't compile fourth derivative code" ON)

######################### Process & Validate Options ###########################
include(autocmake_safeguards)
include(custom_static_library)

################################# Main Project #################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})

# Need BASH to run test suite
find_program(BASH_EXECUTABLE bash)

# link -lm only if necessary
find_package(StandardMathLibraryC)
# check if cbrt exists and declare HAVE_CBRT if it does
check_c_source_compiles (
  "#include <math.h>
int main() {
  return (int)(cbrt(0.8));
}" HAVE_CBRT)
if (HAVE_CBRT)
  add_definitions (-DHAVE_CBRT)
endif (HAVE_CBRT)

# Turn off derivative code
if(DISABLE_VXC)
  add_definitions (-DXC_DONT_COMPILE_VXC)
  add_definitions (-DXC_DONT_COMPILE_FXC)
  add_definitions (-DXC_DONT_COMPILE_KXC)
  add_definitions (-DXC_DONT_COMPILE_LXC)
endif(DISABLE_VXC)
if(DISABLE_FXC)
  add_definitions (-DXC_DONT_COMPILE_FXC)
  add_definitions (-DXC_DONT_COMPILE_KXC)
  add_definitions (-DXC_DONT_COMPILE_LXC)
endif(DISABLE_FXC)
if(DISABLE_KXC)
  add_definitions (-DXC_DONT_COMPILE_KXC)
  add_definitions (-DXC_DONT_COMPILE_LXC)
endif(DISABLE_KXC)
if(DISABLE_LXC)
  add_definitions (-DXC_DONT_COMPILE_LXC)
endif(DISABLE_LXC)

# <<<  Build  >>>

# extract project version from source
file(STRINGS "configure.ac" _libxc_configure_ac
     REGEX "AC_INIT")
foreach(ver ${_libxc_configure_ac})
    if (ver MATCHES "^AC_INIT..libxc...([0-9]+).([0-9]+).([0-9]+).*$")
        set(PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})
        set(PROJECT_VERSION_MINOR ${CMAKE_MATCH_2})
        set(PROJECT_VERSION_MICRO ${CMAKE_MATCH_3})
    endif()
endforeach()

set(PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_MICRO})
message(STATUS "Version: Full ${PROJECT_VERSION}")

# repurpose xc_version from Make for CMake
set(VERSION ${PROJECT_VERSION})
set(XC_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set(XC_MINOR_VERSION ${PROJECT_VERSION_MINOR})
set(XC_MICRO_VERSION ${PROJECT_VERSION_MICRO})
configure_file(xc_version.h.in xc_version.h @ONLY)

# create dummy config.h
configure_file(config.h.cmake.in config.h @ONLY)
configure_file(cmake/libxc.pc.in libxc.pc @ONLY)

# extract project soversion from source
file(STRINGS "configure.ac" _libxc_configure_ac
  REGEX "XC_(CURRENT|REVISION|AGE)=")
foreach(ver ${_libxc_configure_ac})
  if (ver MATCHES "XC_(CURRENT|REVISION|AGE)=+([^ ]+)$")
    set(XC_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
  endif()
endforeach()
set(${PROJECT_NAME}_SOVERSION ${XC_CURRENT}:${XC_REVISION}:${XC_AGE})
math(EXPR ${PROJECT_NAME}_SOMAJOR "${XC_CURRENT} - ${XC_AGE}")
message(STATUS "SO Version: Full ${${PROJECT_NAME}_SOVERSION} Major ${${PROJECT_NAME}_SOMAJOR}")

set(raw_sources_list
  bessel.c
  expint_e1.c
  faddeeva.c
  func_info.c
  func_reference.c
  functionals.c
  hybrids.c
  gga.c
  gga_c_acgga.c
  gga_c_acggap.c
  gga_c_am05.c
  gga_c_bmk.c
  gga_c_chachiyo.c
  gga_c_cs1.c
  gga_c_ft97.c
  gga_c_gapc.c
  gga_c_gaploc.c
  gga_c_hcth_a.c
  gga_c_lm.c
  gga_c_lyp.c
  gga_c_lypr.c
  gga_c_op_b88.c
  gga_c_op_g96.c
  gga_c_op_pbe.c
  gga_c_op_pw91.c
  gga_c_op_xalpha.c
  gga_c_optc.c
  gga_c_p86.c
  gga_c_p86vwn.c
  gga_c_pbe.c
  gga_c_pbe_vwn.c
  gga_c_pbeloc.c
  gga_c_pw91.c
  gga_c_q2d.c
  gga_c_regtpss.c
  gga_c_revtca.c
  gga_c_scan_e0.c
  gga_c_sg4.c
  gga_c_sogga11.c
  gga_c_tca.c
  gga_c_w94.c
  gga_c_wi.c
  gga_c_wl.c
  gga_c_zpbeint.c
  gga_c_zvpbeint.c
  gga_c_zvpbeloc.c
  gga_k_dk.c
  gga_k_pg.c
  gga_k_rational_p.c
  gga_k_meyer.c
  gga_k_ol1.c
  gga_k_ol2.c
  gga_x_ol2.c
  gga_k_pearson.c
  gga_k_tflw.c
  gga_k_thakkar.c
  gga_k_gds08.c
  gga_k_exp4.c
  gga_x_2d_b86.c
  gga_x_2d_b86_mgc.c
  gga_x_2d_b88.c
  gga_x_2d_pbe.c
  gga_x_airy.c
  gga_x_ak13.c
  gga_x_am05.c
  gga_x_b86.c
  gga_x_b88.c
  gga_k_llp.c
  gga_k_pw86.c
  gga_k_mpbe.c
  gga_k_apbe.c
  gga_k_apbeint.c
  gga_k_lc94.c
  gga_k_lkt.c
  gga_k_vt84f.c
  gga_x_bayesian.c
  gga_x_beefvdw.c
  gga_x_bpccac.c
  gga_x_c09x.c
  gga_x_cap.c
  gga_x_chachiyo.c
  gga_x_dk87.c
  gga_x_ev93.c
  gga_x_ft97.c
  gga_x_g96.c
  gga_x_gg99.c
  gga_x_hcth_a.c
  gga_x_herman.c
  gga_x_hjs.c
  gga_x_hjs_b88_v2.c
  gga_x_htbs.c
  gga_x_ityh.c
  gga_x_ityh_optx.c
  gga_x_ityh_pbe.c
  gga_x_kt.c
  gga_x_lag.c
  gga_x_lb.c
  gga_x_lg93.c
  gga_x_lspbe.c
  gga_x_lsrpbe.c
  gga_x_lv_rpw86.c
  gga_x_mpbe.c
  gga_x_n12.c
  gga_x_ncap.c
  gga_x_optx.c
  gga_x_pbe.c
  gga_x_pbea.c
  gga_x_pbeint.c
  gga_x_pbepow.c
  gga_x_pbetrans.c
  gga_x_pw86.c
  gga_x_pw91.c
  gga_x_q1d.c
  gga_x_q2d.c
  gga_x_rge2.c
  gga_x_rpbe.c
  gga_x_s12.c
  gga_x_sfat.c
  gga_x_sfat_pbe.c
  gga_x_sg4.c
  gga_x_sogga11.c
  gga_x_ssb_sw.c
  gga_x_vmt.c
  gga_x_vmt84.c
  gga_x_wc.c
  gga_x_wpbeh.c
  gga_x_fd_lb94.c
  gga_xc_1w.c
  gga_xc_b97.c
  gga_xc_edf1.c
  gga_xc_oblyp_d.c
  gga_xc_th1.c
  gga_xc_th2.c
  gga_xc_th3.c
  gga_xc_vv10.c
  gga_k_lgap.c
  gga_k_lgap_ge.c
  gga_c_ccdf.c
  hyb_gga_xc_b1wc.c
  hyb_gga_xc_b3lyp.c
  hyb_gga_xc_cam_b3lyp.c
  hyb_gga_xc_camy_b3lyp.c
  hyb_gga_xc_camy_blyp.c
  hyb_gga_xc_edf2.c
  hyb_gga_xc_hse.c
  hyb_gga_xc_lc_blyp.c
  hyb_gga_xc_lcy_blyp.c
  hyb_gga_xc_lcy_pbe.c
  hyb_gga_xc_o3lyp.c
  hyb_gga_xc_pbeh.c
  hyb_gga_xc_wb97.c
  hyb_lda_xc_cam_lda0.c
  hyb_lda_xc_bn05.c
  hyb_mgga_x_dldf.c
  hyb_mgga_x_m05.c
  hyb_mgga_x_mvsh.c
  hyb_mgga_xc_b88b95.c
  hyb_mgga_xc_kcis.c
  hyb_mgga_xc_tpssh.c
  hyb_mgga_xc_wb97mv.c
  hyb_mgga_x_js18.c
  hyb_mgga_x_pjs18.c
  hyb_gga_xc_cam_o3lyp.c
  integrate.c
  lda.c
  lda_c_1d_csc.c
  lda_c_1d_loos.c
  lda_c_2d_amgb.c
  lda_c_2d_prm.c
  lda_c_chachiyo.c
  lda_c_chachiyo_mod.c
  lda_c_gombas.c
  lda_c_hl.c
  lda_c_lp96.c
  lda_c_ml1.c
  lda_c_pk09.c
  lda_c_pw.c
  lda_c_pz.c
  lda_c_rc04.c
  lda_c_rpa.c
  lda_c_vwn.c
  lda_c_vwn_1.c
  lda_c_vwn_2.c
  lda_c_vwn_3.c
  lda_c_vwn_4.c
  lda_c_vwn_rpa.c
  lda_c_wigner.c
  lda_c_gk72.c
  lda_c_w20.c
  lda_k_tf.c
  lda_k_zlp.c
  lda_k_gds08_worker.c
  lda_x.c
  lda_x_1d_soft.c
  lda_x_1d_exponential.c
  lda_x_2d.c
  lda_x_erf.c
  lda_x_rel.c
  lda_xc_1d_ehwlrg.c
  lda_xc_ksdt.c
  lda_xc_teter93.c
  lda_xc_zlp.c
  lda_c_pmgb06.c
  lda_xc_tih.c
  lda_x_sloc.c
  lda_x_yukawa.c
  mgga.c
  mgga_c_b88.c
  mgga_c_b94.c
  mgga_c_bc95.c
  mgga_c_cs.c
  mgga_c_kcis.c
  mgga_c_kcisk.c
  mgga_xc_lp90.c
  mgga_c_m05.c
  mgga_c_m06l.c
  mgga_c_m08.c
  mgga_c_pkzb.c
  mgga_c_revscan.c
  mgga_c_revtpss.c
  mgga_c_rscan.c
  mgga_c_scan.c
  mgga_c_scanl.c
  mgga_c_tpss.c
  mgga_c_tpssloc.c
  mgga_c_vsxc.c
  mgga_k_pc07.c
  mgga_k_csk.c
  mgga_k_csk_loc.c
  mgga_k_pgslb.c
  mgga_x_2d_prhg07.c
  mgga_x_2d_prp10.c
  mgga_x_br89.c
  mgga_x_br89_explicit.c
  mgga_x_tb09.c
  mgga_x_gvt4.c
  mgga_x_gx.c
  mgga_x_jk.c
  mgga_x_lta.c
  mgga_x_m06l.c
  mgga_x_m08.c
  mgga_x_m11.c
  mgga_x_m11_l.c
  mgga_x_mbeef.c
  mgga_x_mbeefvdw.c
  mgga_x_mn12.c
  mgga_x_ms.c
  mgga_x_msb.c
  mgga_x_mvs.c
  mgga_x_mvsb.c
  mgga_x_pbe_gx.c
  mgga_x_pkzb.c
  mgga_x_rscan.c
  mgga_x_sa_tpss.c
  mgga_x_scan.c
  mgga_x_scanl.c
  mgga_x_tau_hcth.c
  mgga_x_tm.c
  mgga_x_regtm.c
  mgga_x_revtm.c
  mgga_x_tpss.c
  mgga_x_regtpss.c
  mgga_x_vt84.c
  mgga_x_rtpss.c
  mgga_xc_b97mv.c
  mgga_xc_cc06.c
  mgga_xc_hle17.c
  mgga_xc_otpss_d.c
  mgga_xc_zlp.c
  mgga_xc_b98.c
  mgga_x_2d_js17.c
  mgga_x_edmgga.c
  mgga_x_gdme.c
  mgga_x_rlda.c
  mgga_x_scanl.c
  mgga_c_scanl.c
  mgga_x_mbrxh_bg.c
  mgga_x_mbrxc_bg.c
  mgga_x_task.c
  mgga_x_mggac.c
  mgga_x_th.c
  mgga_x_mbr.c
  mgga_c_ltapw.c
  mgga_x_r2scan.c
  mgga_c_r2scan.c
  mgga_k_lk.c
  mgga_k_rda.c
  mgga_k_gea2.c
  mgga_k_gea4.c
  mgga_x_r2scanl.c
  mgga_c_r2scanl.c
  mgga_c_rregtm.c
  mgga_x_mcml.c
  mgga_c_rmggac.c
  mix_func.c
  deorbitalize_func.c
  references.c
  special_functions.c
  math_brent.c
  util.c
  version.c
  )

set(src_prefix "src/")
string(REGEX REPLACE "([^;]+)" "${src_prefix}\\1" sources_list "${raw_sources_list}")

set(raw_sources_list_f03
  src/libxc_master.F90
  )
set(raw_sources_list_f90
  src/libxc_master_f90.F90
  )

# when headers namespaced, xc_version include in xc.h needs to be local, not
#   system to be found
file(READ ${src_prefix}xc.h _src_contents)
string(REPLACE "<xc_version.h>" "\"xc_version.h\"" _quoted_src "${_src_contents}")
file(WRITE ${PROJECT_BINARY_DIR}/${src_prefix}xc.h "${_quoted_src}")

# STATIC/SHARED on below governed by BUILD_SHARED_LIBS
add_library(xc ${sources_list})
target_link_libraries(xc INTERFACE ${STANDARD_MATH_LIBRARY})
set_target_properties(xc PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
  SOVERSION ${${PROJECT_NAME}_SOMAJOR})
if(${BUILD_SHARED_LIBS})
  target_link_libraries(xc PRIVATE ${LIBC_INTERJECT})
  if(APPLE)
    set_target_properties(xc PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif()
endif()

if(ENABLE_FORTRAN)
  enable_language(Fortran)
  configure_file(cmake/libxcf03.pc.in libxcf03.pc @ONLY)
  configure_file(cmake/libxcf90.pc.in libxcf90.pc @ONLY)
  add_library(xcf03 ${raw_sources_list_f03})
  target_link_libraries(xcf03 xc)
  set_target_properties(xcf03 PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
    SOVERSION ${${PROJECT_NAME}_SOMAJOR})
  add_library(xcf90 ${raw_sources_list_f90})
  target_link_libraries(xcf90 xc)
  set_target_properties(xcf90 PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
    SOVERSION ${${PROJECT_NAME}_SOMAJOR})
endif()

# run get_funcs.pl to generate headers for xc inclusion
find_package(Perl REQUIRED)
set(gen_dir "${PROJECT_BINARY_DIR}/gen_funcidx")
execute_process(
  COMMAND ${CMAKE_COMMAND} -E make_directory ${gen_dir}
  TIMEOUT 5
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
execute_process(
  COMMAND ${PERL_EXECUTABLE} scripts/get_funcs.pl ${src_prefix} ${gen_dir}
  TIMEOUT 60
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(xc-info "${src_prefix}/xc-info.c")
target_link_libraries(xc-info xc)

add_executable(xc-threshold "${src_prefix}/xc-threshold.c")
target_link_libraries(xc-threshold xc)


include_directories(${PROJECT_SOURCE_DIR}/${src_prefix}  # for util.h
  ${PROJECT_BINARY_DIR}/${src_prefix}  # for xc.h
  ${PROJECT_BINARY_DIR}  # for xc_version.h, config.h
  ${gen_dir})  # for xc_funcs.h

if(BUILD_TESTING)
  find_program(BZip2_EXECUTABLE
    NAMES bzip2
    DOC "Path to zipping utility")
  if(BZip2_EXECUTABLE)
    enable_testing ()
    add_subdirectory(testsuite)  
  else()
    message(FATAL_ERROR "Install `bzip2` command to enable tests")
  endif()
endif()

# <<<  Install  >>>

# by default, headers NOT namespace protected
install(FILES ${PROJECT_BINARY_DIR}/${src_prefix}/xc.h
  ${gen_dir}/xc_funcs.h
  ${gen_dir}/xc_funcs_worker.h
  ${PROJECT_BINARY_DIR}/xc_version.h
  ${PROJECT_SOURCE_DIR}/${src_prefix}/xc_funcs_removed.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}${NAMESPACE_INSTALL_INCLUDEDIR})
install(TARGETS xc-info
  OPTIONAL
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS xc
  EXPORT c_interface
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ENABLE_FORTRAN)
  install(TARGETS xc xcf03 xcf90
    EXPORT f_interface
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(FILES "${PROJECT_BINARY_DIR}/xc_f03_lib_m.mod"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}${NAMESPACE_INSTALL_INCLUDEDIR})
  install(FILES "${PROJECT_BINARY_DIR}/xc_f90_lib_m.mod"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}${NAMESPACE_INSTALL_INCLUDEDIR})
  if(NOT MSVC)
    install(FILES ${PROJECT_BINARY_DIR}/libxcf03.pc
                  ${PROJECT_BINARY_DIR}/libxcf90.pc
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/)
  endif()
endif()

# <<<  Export Interface  >>>

target_compile_definitions(xc INTERFACE USING_${PN})
target_include_directories(xc INTERFACE
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# <<<  Export Config  >>>

# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
configure_package_config_file(cmake/${PN}Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
  INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
  VERSION ${${PN}_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
  DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT c_interface
  NAMESPACE "${PN}::"
  FILE "${PN}Targets-C.cmake"
  DESTINATION ${CMAKECONFIG_INSTALL_DIR})

if(ENABLE_FORTRAN)
    install(EXPORT f_interface
      NAMESPACE "${PN}::"
      FILE "${PN}Targets-Fortran.cmake"
      DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()

if(ENABLE_PYTHON)
  set(SOURCE_PYTHON_API
    pylibxc/__init__.py
    pylibxc/core.py
    pylibxc/flags.py
    pylibxc/functional.py
    pylibxc/structs.py
    pylibxc/util.py
    pylibxc/version.py
  )
  install(
    FILES ${SOURCE_PYTHON_API}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pylibxc
  )
  if(UNIX)
    install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink \
      ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libxc.so.${${PROJECT_NAME}_SOMAJOR} \
      ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pylibxc/libxc.so)")
  endif()
endif()

export(EXPORT c_interface
    NAMESPACE "${PN}::"
    FILE "${PROJECT_BINARY_DIR}/${PN}Targets.cmake")

if(NOT MSVC)
  install(FILES ${PROJECT_BINARY_DIR}/libxc.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/)
endif()
