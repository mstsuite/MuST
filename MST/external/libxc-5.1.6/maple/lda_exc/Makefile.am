## Copyright (C) 2017 M. Marques
##
## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
##


lda_order_4 = lda_c_1d_csc.mpl lda_c_pz.mpl lda_x_1d_exponential.mpl	\
  lda_c_1d_loos.mpl lda_c_rc04.mpl lda_x_1d_soft.mpl			\
  lda_c_2d_amgb.mpl lda_c_rpa.mpl lda_x_2d.mpl lda_c_2d_prm.mpl		\
  lda_c_vwn_1.mpl lda_xc_1d_ehwlrg.mpl lda_c_chachiyo.mpl		\
  lda_c_vwn_3.mpl hyb_lda_xc_bn05.mpl lda_c_gk72.mpl lda_c_vwn_3.mpl	\
  lda_xc_ksdt.mpl lda_c_gombas.mpl lda_c_vwn_4.mpl lda_xc_teter93.mpl	\
  lda_c_hl.mpl lda_c_vwn.mpl lda_xc_zlp.mpl lda_c_lp96.mpl		\
  lda_c_vwn_rpa.mpl lda_x_erf.mpl lda_c_ml1.mpl lda_c_wigner.mpl	\
  lda_x.mpl lda_k_gds08_worker.mpl lda_x_rel.mpl lda_c_pmgb06.mpl	\
  lda_k_tf.mpl lda_c_pw.mpl lda_k_zlp.mpl lda_x_sloc.mpl		\
  lda_c_chachiyo_mod.mpl lda_c_w20.mpl lda_x_yukawa.mpl

lda_order_3 = lda_c_pk09.mpl

dist_noinst_SCRIPTS = $(lda_order_4) $(lda_order_3)

maple2c:
	@for funct in $(lda_order_4); do \
	  echo Processing: $${funct%.*}; \
	  $(top_srcdir)/scripts/maple2c.pl $(top_srcdir) $${funct%.*} 4; \
	done
	@for funct in $(lda_order_3); do \
	  echo Processing: $${funct%.*}; \
	  $(top_srcdir)/scripts/maple2c.pl $(top_srcdir) $${funct%.*} 3; \
	done

