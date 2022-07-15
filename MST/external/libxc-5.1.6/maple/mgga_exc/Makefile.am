## Copyright (C) 2017 M. Marques
##
## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

mgga_order_4 = mgga_x_lta.mpl mgga_c_ltapw.mpl mgga_x_th.mpl		\
  mgga_x_mbeefvdw.mpl hyb_mgga_x_dldf.mpl mgga_k_pc07.mpl		\
  mgga_x_2d_prhg07.mpl hyb_mgga_x_m05.mpl mgga_x_br89_explicit.mpl	\
  mgga_x_mn12.mpl mgga_x_msb.mpl mgga_x_ms.mpl mgga_x_mvsb.mpl		\
  mgga_x_mvs.mpl mgga_x_pbe_gx.mpl mgga_x_gvt4.mpl mgga_x_pkzb.mpl	\
  mgga_x_gx.mpl mgga_x_rtpss.mpl mgga_x_sa_tpss.mpl mgga_x_m06l.mpl	\
  mgga_x_scan.mpl mgga_x_m08.mpl mgga_x_tau_hcth.mpl mgga_x_m11_l.mpl	\
  mgga_x_tm.mpl mgga_x_m11.mpl mgga_x_tpss.mpl mgga_x_mbeef.mpl		\
  mgga_x_vt84.mpl mgga_x_regtpss.mpl mgga_x_rscan.mpl			\
  mgga_x_edmgga.mpl mgga_x_gdme.mpl mgga_x_revtm.mpl mgga_x_rlda.mpl	\
  mgga_x_mbrxh_bg.mpl mgga_x_mbrxc_bg.mpl mgga_x_task.mpl		\
  hyb_mgga_x_js18.mpl hyb_mgga_x_pjs18.mpl hyb_mgga_xc_wb97mv.mpl	\
  mgga_c_vsxc.mpl mgga_c_b88.mpl mgga_xc_b97mv.mpl mgga_c_bc95.mpl	\
  mgga_xc_b98.mpl mgga_c_cs.mpl mgga_xc_cc06.mpl mgga_c_kcis.mpl	\
  mgga_xc_lp90.mpl mgga_c_m05.mpl mgga_xc_zlp.mpl mgga_c_m06l.mpl	\
  mgga_c_m08.mpl mgga_c_pkzb.mpl mgga_c_revscan.mpl			\
  mgga_c_revtpss.mpl mgga_c_scan.mpl mgga_c_tpssloc.mpl			\
  mgga_c_tpss.mpl mgga_c_rscan.mpl mgga_x_br89.mpl mgga_x_2d_js17.mpl	\
  mgga_x_mggac.mpl mgga_x_mbr.mpl mgga_x_jk.mpl mgga_x_r2scan.mpl	\
  mgga_c_r2scan.mpl mgga_k_lk.mpl mgga_k_rda.mpl mgga_k_gea2.mpl        \
  mgga_k_gea4.mpl mgga_k_csk.mpl mgga_k_csk_loc.mpl mgga_c_b94.mpl      \
  mgga_k_pgslb.mpl mgga_c_kcisk.mpl mgga_x_regtm.mpl mgga_c_rregtm.mpl  \
  mgga_x_mcml.mpl mgga_c_rmggac.mpl

mgga_order_3 = mgga_c_b94.mpl

dist_noinst_HEADERS = $(mgga_order_4) $(mgga_order_3)

maple2c:
	@for funct in $(mgga_order_4); do \
	  echo Processing: $${funct%.*}; \
	  $(top_srcdir)/scripts/maple2c.pl $(top_srcdir) $${funct%.*} 4; \
	done
	@for funct in $(mgga_order_3); do \
	  echo Processing: $${funct%.*}; \
	  $(top_srcdir)/scripts/maple2c.pl $(top_srcdir) $${funct%.*} 3; \
	done
