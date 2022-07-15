/*
 Copyright (C) 2006-2007 M.A.L. Marques and
                    2015 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_EDF2        476 /* Empirical functional from Lin, George and Gill */

static void
hyb_gga_xc_edf2_init(xc_func_type *p)
{
  static int    funcs_id  [6] = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_X_B88, XC_LDA_C_VWN, XC_GGA_C_LYP, XC_GGA_C_LYP};
  static double funcs_coef[6] = {0.2811, 0.6227, -0.0551, 0.3029, 0.5998, -0.0053};

  static double par_x_b88[] = {0.0035, 6.0};
  static double par_c_lyp[] = {0.055, 0.158, 0.25, 0.3505};

  xc_mix_init(p, 6, funcs_id, funcs_coef);
  xc_func_set_ext_params(p->func_aux[2], par_x_b88);
  xc_func_set_ext_params(p->func_aux[5], par_c_lyp);

  xc_hyb_init_hybrid(p, 0.1695);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_edf2 = {
  XC_HYB_GGA_XC_EDF2,
  XC_EXCHANGE_CORRELATION,
  "EDF2",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Lin2004_365, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  hyb_gga_xc_edf2_init,
  NULL, NULL, NULL, NULL
};
