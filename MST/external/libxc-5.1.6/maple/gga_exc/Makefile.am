## Copyright (C) 2017 M. Marques
##
## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

gga_order_4 = gga_k_exp4.mpl gga_x_hjs_b88_v2.mpl gga_k_lc94.mpl	\
  gga_x_hjs.mpl gga_k_llp.mpl gga_x_htbs.mpl gga_k_meyer.mpl		\
  gga_x_ityh.mpl gga_k_mpbe.mpl gga_x_kt.mpl gga_k_ol1.mpl		\
  gga_x_lag.mpl gga_k_ol2.mpl gga_x_lg93.mpl gga_k_pearson.mpl		\
  gga_x_lspbe.mpl gga_k_pw86.mpl gga_x_lsrpbe.mpl gga_k_tflw.mpl	\
  gga_x_lv_rpw86.mpl gga_k_thakkar.mpl gga_x_mpbe.mpl			\
  gga_x_2d_b86_mgc.mpl gga_x_n12.mpl gga_x_2d_b86.mpl gga_x_ol2.mpl	\
  gga_x_2d_b88.mpl gga_x_optx.mpl gga_x_2d_pbe.mpl gga_x_pbea.mpl	\
  gga_x_airy.mpl gga_x_pbeint.mpl gga_x_ak13.mpl gga_x_pbe.mpl		\
  gga_x_am05.mpl gga_x_pbepow.mpl gga_x_b86.mpl gga_x_pbetrans.mpl	\
  gga_x_b88.mpl gga_x_pw86.mpl gga_x_bayesian.mpl gga_x_pw91.mpl	\
  gga_x_beefvdw.mpl gga_x_q2d.mpl gga_x_bpccac.mpl gga_x_rge2.mpl	\
  gga_x_c09x.mpl gga_x_rpbe.mpl gga_x_cap.mpl gga_x_sfat.mpl		\
  gga_x_sfat_pbe.mpl gga_x_chachiyo.mpl gga_x_sg4.mpl			\
  gga_x_sogga11.mpl gga_x_ssb_sw.mpl gga_x_vmt84.mpl gga_x_dk87.mpl	\
  gga_x_vmt.mpl gga_x_ev93.mpl gga_x_wc.mpl gga_x_ft97.mpl		\
  gga_k_apbeint.mpl gga_x_g96.mpl gga_k_apbe.mpl gga_x_hcth_a.mpl	\
  gga_k_dk.mpl gga_x_herman.mpl gga_x_gg99.mpl gga_x_fd_lb94.mpl	\
  gga_x_ncap.mpl gga_x_s12.mpl gga_k_lkt.mpl gga_c_acgga.mpl		\
  gga_c_acggap.mpl gga_c_am05.mpl gga_c_bmk.mpl gga_c_cs1.mpl		\
  gga_c_ft97.mpl gga_c_gapc.mpl gga_c_gaploc.mpl gga_c_hcth_a.mpl	\
  gga_c_lm.mpl gga_c_lyp.mpl gga_c_op_b88.mpl gga_c_op_g96.mpl		\
  gga_c_op_pbe.mpl gga_c_op_pw91.mpl gga_c_optc.mpl			\
  gga_c_op_xalpha.mpl gga_c_p86.mpl gga_c_pbeloc.mpl gga_c_pbe.mpl	\
  gga_c_pw91.mpl gga_c_q2d.mpl gga_c_regtpss.mpl gga_c_revtca.mpl	\
  gga_c_scan_e0.mpl gga_c_sg4.mpl gga_xc_b97.mpl gga_c_sogga11.mpl	\
  gga_c_tca.mpl gga_xc_th1.mpl gga_c_w94.mpl gga_xc_th2.mpl		\
  gga_c_wi.mpl gga_xc_th3.mpl gga_c_wl.mpl gga_c_zpbeint.mpl		\
  gga_c_zvpbeint.mpl hyb_gga_xc_wb97.mpl gga_c_zvpbeloc.mpl		\
  gga_c_chachiyo.mpl gga_c_pbe_vwn.mpl gga_x_wpbeh.mpl gga_k_lgap.mpl   \
  gga_k_vt84f.mpl gga_k_lgap_ge.mpl gga_c_ccdf.mpl gga_c_p86vwn.mpl     \
  gga_x_ityh_optx.mpl gga_x_ityh_pbe.mpl gga_c_lypr.mpl gga_x_q1d.mpl


dist_noinst_SCRIPTS = $(gga_order_4)

maple2c:
	@for funct in $(gga_order_4); do \
	  echo Processing: $${funct%.*}; \
    $(top_srcdir)/scripts/maple2c.pl $(top_srcdir) $${funct%.*} 4; \
	done
