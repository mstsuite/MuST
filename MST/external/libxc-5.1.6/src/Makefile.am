## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

SUBDIRS = maple2c

if HAVE_CUDA
#this is a hack for nvcc, this should be improved.
CC+="-dc"
endif

bin_PROGRAMS = xc-info xc-threshold
noinst_PROGRAMS = genwiki xc-sanity

xc_info_SOURCES = xc-info.c
xc_info_LDADD = libxc.la -lm
xc_info_LDFLAGS =

xc_threshold_SOURCES = xc-threshold.c
xc_threshold_LDADD = libxc.la -lm
xc_threshold_LDFLAGS =

genwiki_SOURCES = genwiki.c
genwiki_LDADD = libxc.la -lm
genwiki_LDFLAGS =

xc_sanity_SOURCES = xc-sanity.c
xc_sanity_LDADD = libxc.la -lm
xc_sanity_LDFLAGS =

lib_LTLIBRARIES = libxc.la

libxc_la_FUNC_SOURCES = \
	special_functions.c bessel.c expint_e1.c faddeeva.c integrate.c		\
	math_brent.c util.c mix_func.c deorbitalize_func.c func_info.c functionals.c hybrids.c \
	version.c func_reference.c \
	lda.c \
	lda_x.c lda_x_2d.c lda_c_wigner.c lda_c_gombas.c lda_c_rpa.c	\
	lda_c_hl.c lda_c_vwn.c lda_c_vwn_1.c lda_c_vwn_2.c		\
	lda_c_vwn_3.c lda_c_vwn_4.c lda_c_vwn_rpa.c lda_c_pz.c		\
	lda_c_pw.c lda_c_pk09.c lda_c_lp96.c lda_c_ml1.c lda_c_rc04.c	\
	lda_xc_teter93.c lda_c_1d_csc.c lda_c_1d_loos.c			\
	lda_c_2d_amgb.c lda_c_2d_prm.c lda_k_tf.c lda_xc_zlp.c		\
	lda_xc_ksdt.c lda_c_chachiyo.c lda_c_chachiyo_mod.c		\
	lda_xc_1d_ehwlrg.c lda_x_rel.c lda_x_erf.c lda_k_zlp.c		\
	lda_c_gk72.c hyb_lda_xc_bn05.c lda_c_pmgb06.c			\
	lda_k_gds08_worker.c lda_xc_tih.c lda_x_1d_exponential.c	\
	lda_x_1d_soft.c hyb_lda_xc_cam_lda0.c lda_x_sloc.c lda_c_w20.c  \
	lda_x_yukawa.c 	                                                \
	gga.c \
	gga_x_lg93.c gga_x_pbe.c gga_x_rpbe.c gga_x_pbea.c		\
	gga_x_mpbe.c gga_x_herman.c gga_x_b86.c gga_x_b88.c		\
	gga_x_g96.c gga_x_pw86.c gga_x_pw91.c gga_x_optx.c		\
	gga_x_airy.c gga_x_c09x.c gga_x_dk87.c gga_x_ft97.c		\
	gga_x_wc.c gga_x_am05.c gga_x_bayesian.c gga_x_kt.c		\
	gga_x_htbs.c gga_c_pbe.c gga_c_lyp.c gga_c_p86.c gga_c_pw91.c	\
	gga_c_am05.c gga_c_lm.c gga_c_wl.c gga_c_wi.c gga_x_lb.c	\
	gga_xc_b97.c gga_xc_edf1.c gga_xc_1w.c gga_c_optc.c		\
	gga_c_tca.c gga_x_bpccac.c gga_x_sogga11.c gga_c_sogga11.c	\
	gga_xc_th1.c gga_x_wpbeh.c gga_x_hjs.c gga_x_ityh.c		\
	gga_x_sfat.c gga_x_sfat_pbe.c gga_x_ev93.c gga_x_ak13.c		\
	gga_x_q1d.c gga_x_q2d.c gga_c_q2d.c gga_x_ssb_sw.c gga_c_ft97.c		\
	gga_x_n12.c gga_x_lv_rpw86.c gga_x_2d_b86.c gga_k_llp.c		\
	gga_x_2d_b86_mgc.c gga_x_2d_b88.c gga_x_2d_pbe.c		\
	gga_x_lspbe.c gga_x_lsrpbe.c gga_k_tflw.c gga_k_pearson.c	\
	gga_x_ol2.c gga_k_ol1.c gga_k_ol2.c gga_k_dk.c gga_k_meyer.c	\
	gga_x_vmt.c gga_x_hcth_a.c gga_x_cap.c gga_x_beefvdw.c		\
	gga_x_sg4.c gga_x_gg99.c gga_x_pbepow.c gga_x_lag.c		\
	gga_k_thakkar.c gga_xc_oblyp_d.c hyb_gga_xc_b3lyp.c		\
	hyb_gga_xc_o3lyp.c gga_x_pbeint.c gga_x_rge2.c gga_x_vmt84.c	\
	gga_x_hjs_b88_v2.c gga_c_revtca.c gga_c_acgga.c			\
	gga_c_acggap.c gga_c_zpbeint.c gga_c_pbeloc.c gga_c_sg4.c	\
	gga_c_regtpss.c gga_c_scan_e0.c gga_xc_th2.c gga_xc_th3.c	\
	gga_c_hcth_a.c hyb_gga_xc_wb97.c gga_c_bmk.c gga_c_gapc.c	\
	gga_c_gaploc.c gga_c_zvpbeint.c gga_c_w94.c gga_c_cs1.c		\
	gga_c_op_xalpha.c gga_c_op_b88.c gga_c_op_pbe.c			\
	gga_c_op_g96.c gga_c_op_pw91.c gga_x_pbetrans.c			\
	gga_k_apbeint.c gga_x_chachiyo.c gga_k_gds08.c gga_k_exp4.c	\
	hyb_gga_xc_pbeh.c hyb_gga_xc_b1wc.c hyb_gga_xc_hse.c		\
	hyb_gga_xc_cam_b3lyp.c  gga_k_lc94.c	                        \
	gga_c_zvpbeloc.c hyb_gga_xc_camy_blyp.c				\
	hyb_gga_xc_camy_b3lyp.c hyb_gga_xc_lcy_pbe.c			\
	hyb_gga_xc_lcy_blyp.c gga_xc_vv10.c hyb_gga_xc_edf2.c		\
	gga_k_pw86.c gga_k_mpbe.c gga_k_apbe.c gga_x_fd_lb94.c		\
	gga_x_ncap.c gga_c_chachiyo.c hyb_gga_xc_lc_blyp.c		\
	gga_c_pbe_vwn.c gga_x_s12.c   \
	gga_k_lkt.c gga_k_lgap.c gga_k_vt84f.c   \
	gga_k_lgap_ge.c	gga_c_ccdf.c gga_k_rational_p.c gga_k_pg.c      \
	gga_c_p86vwn.c gga_x_ityh_optx.c hyb_gga_xc_cam_o3lyp.c \
	gga_x_ityh_pbe.c gga_c_lypr.c \
	mgga.c \
	mgga_x_lta.c mgga_x_tpss.c mgga_x_br89.c mgga_x_gvt4.c		\
	mgga_x_m06l.c mgga_x_tau_hcth.c mgga_x_2d_prhg07.c		\
	mgga_x_2d_prp10.c mgga_x_pkzb.c hyb_mgga_x_m05.c mgga_x_m08.c	\
	mgga_x_ms.c mgga_x_mn12.c mgga_c_vsxc.c mgga_c_pkzb.c		\
	mgga_c_bc95.c mgga_c_m08.c mgga_xc_cc06.c mgga_c_cs.c		\
	mgga_xc_otpss_d.c mgga_xc_zlp.c hyb_mgga_xc_b88b95.c		\
	mgga_x_mbeef.c mgga_xc_b97mv.c mgga_x_scan.c mgga_c_scan.c	\
	hyb_mgga_xc_tpssh.c mgga_x_mvs.c hyb_mgga_x_mvsh.c		\
	mgga_x_tm.c mgga_x_vt84.c mgga_x_sa_tpss.c mgga_k_pc07.c	\
	mgga_x_mbeefvdw.c mgga_x_m11.c mgga_x_m11_l.c			\
	hyb_mgga_xc_wb97mv.c mgga_c_m05.c mgga_c_m06l.c			\
	hyb_mgga_x_dldf.c mgga_c_tpss.c mgga_c_revtpss.c		\
	mgga_c_tpssloc.c mgga_xc_hle17.c mgga_c_kcis.c mgga_xc_lp90.c	\
	hyb_mgga_xc_kcis.c mgga_c_b88.c mgga_x_gx.c mgga_x_pbe_gx.c	\
	mgga_c_revscan.c mgga_x_br89_explicit.c mgga_xc_b98.c		\
	mgga_x_rtpss.c mgga_x_msb.c mgga_x_mvsb.c mgga_x_regtpss.c	\
	mgga_x_rscan.c mgga_c_rscan.c mgga_x_2d_js17.c mgga_x_tb09.c	\
	mgga_x_edmgga.c mgga_x_gdme.c mgga_x_rlda.c mgga_x_revtm.c	\
	mgga_x_scanl.c mgga_c_scanl.c hyb_mgga_x_js18.c			\
	hyb_mgga_x_pjs18.c mgga_x_mbrxh_bg.c mgga_x_mbrxc_bg.c		\
	mgga_x_task.c mgga_x_mggac.c mgga_x_th.c mgga_x_mbr.c		\
	mgga_c_ltapw.c mgga_x_jk.c mgga_x_r2scan.c mgga_c_r2scan.c      \
	mgga_c_b94.c mgga_k_lk.c mgga_k_rda.c mgga_k_gea2.c             \
	mgga_k_gea4.c mgga_k_csk.c mgga_k_csk_loc.c mgga_x_r2scanl.c    \
        mgga_c_r2scanl.c mgga_k_pgslb.c mgga_c_kcisk.c mgga_x_regtm.c   \
	mgga_c_rregtm.c mgga_x_mcml.c mgga_c_rmggac.c

libxc_la_SOURCES = references.c $(libxc_la_FUNC_SOURCES)
libxc_la_LDFLAGS = -version-info $(XC_LT_VERSION)

if COMPILE_FORTRAN
  lib_LTLIBRARIES += libxcf03.la libxcf90.la

  nodist_libxcf03_la_SOURCES = libxcf03.f90
  nodist_libxcf90_la_SOURCES = libxcf90.f90

  # libtool stuff
  libxcf03_la_LDFLAGS = -version-info $(XC_LT_VERSION)
  libxcf03_la_LIBADD = libxc.la
  # libtool stuff
  libxcf90_la_LDFLAGS = -version-info $(XC_LT_VERSION)
  libxcf90_la_LIBADD = libxc.la
endif

noinst_HEADERS = \
	string_f.h references.h util.h \
	decl_lda.h decl_gga.h decl_mgga.h \
  work_lda.c work_gga.c work_mgga.c \
	libxc_master.F90 libxc_master_f90.F90
include_HEADERS = xc.h xc_funcs_removed.h
nodist_include_HEADERS = xc_funcs.h xc_funcs_worker.h

BUILT_SOURCES = xc_funcs.h

## declare fortran modules so they are installed
fmoddir = @FMODDIR@
fmod_DATA =

if COMPILE_FORTRAN

# Double precision modules are always compiled
if F90_MOD_UPPERCASE
  XCF03LIBMODS = XC_F03_LIB_M.@ax_cv_f90_modext@
  XCF90LIBMODS = XC_F90_LIB_M.@ax_cv_f90_modext@
else
  XCF03LIBMODS = xc_f03_lib_m.@ax_cv_f90_modext@
  XCF90LIBMODS = xc_f90_lib_m.@ax_cv_f90_modext@
endif
fmod_DATA += $(XCF03LIBMODS) $(XCF90LIBMODS)

endif # COMPILE_FORTRAN


CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* funcs_*.c funcs.h \
  libxc_inc.f90 libxcf03.f90 libxcf90.f90 xc_funcs.h xc_funcs_worker.h \
  xc_f_inc.c libxc_docs.txt

xc_f.o : xc_f_inc.c

xc_funcs.h xc_funcs_worker.h libxc_inc.f90 \
funcs_lda.c funcs_gga.c funcs_mgga.c funcs_key.c: $(srcdir)/../scripts/get_funcs.pl
	$(srcdir)/../scripts/get_funcs.pl $(srcdir) $(top_builddir)/src

references.h references.c: $(srcdir)/../scripts/get_references.pl $(srcdir)/../libxc.bib
	(cd $(top_builddir)/src && $(srcdir)/../scripts/get_references.pl $(srcdir)/../libxc.bib)

LTPREF = $(subst xc_f.lo,,$(firstword $(am_libxcf03_la_OBJECTS)))
# Surprisingly, libtool/automake do not seem to provide a macro or variable for what the object names will be,
# so we have to extract it ourselves. LTPREF=libxcf03_la-, or blank for older libtool/automake versions.
# Produces these warnings on some systems. Does not seem to matter.
#src/Makefile.am:140: subst xc_f.lo,,$(firstword $(am_libxcf03_la_OBJECTS: non-POSIX variable name
#src/Makefile.am:140: (probably a GNU make extension)

LTPREF03 = $(subst libxcf03.lo,,$(firstword $(nodist_libxcf03_la_OBJECTS)))
libxcf03.f90 : libxc_master.F90 libxc_inc.f90
	@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I$(top_builddir)/src $(srcdir)/libxc_master.F90 > $(top_builddir)/src/libxcf03.f90
	@if [ "@F90_ACCEPTS_LINE_NUMBERS@" = "no" ]; then \
		grep -v "^#" $(top_builddir)/src/libxcf03.f90 > $(top_builddir)/src/libxcf03.f91; \
		mv -f $(top_builddir)/src/libxcf03.f91 $(top_builddir)/src/libxcf03.f90; \
	fi

$(XCF03LIBMODS) : $(LTPREF03)libxcf03.lo

LTPREF90 = $(subst libxcf90.lo,,$(firstword $(nodist_libxcf90_la_OBJECTS)))
libxcf90.f90 : libxc_master_f90.F90 libxc_inc.f90
	@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I$(top_builddir)/src $(srcdir)/libxc_master_f90.F90 > $(top_builddir)/src/libxcf90.f90
	@if [ "@F90_ACCEPTS_LINE_NUMBERS@" = "no" ]; then \
		grep -v "^#" $(top_builddir)/src/libxcf90.f90 > $(top_builddir)/src/libxcf90.f91; \
		mv -f $(top_builddir)/src/libxcf90.f91 $(top_builddir)/src/libxcf90.f90; \
	fi

$(XCF90LIBMODS) : $(LTPREF90)libxcf90.lo
