/*
 Copyright (C) 2013 Rolf Wuerdemann, M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define  XC_HYB_GGA_XC_LCY_PBE 467  /* PBE with yukawa screening */

void
xc_hyb_gga_xc_lcy_pbe_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_SFAT_PBE, XC_GGA_C_PBE};
  static double funcs_coef[2];

  double gamma = 0.75;

  funcs_coef[0] = 1.0;
  funcs_coef[1] = 1.0;

  xc_mix_init(p, 2, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[0], &gamma);

  xc_hyb_init_camy(p, 1.0, -1.0, gamma);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lcy_pbe = {
  XC_HYB_GGA_XC_LCY_PBE,
  XC_EXCHANGE_CORRELATION,
  "LCY version of PBE",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Seth2012_901, &xc_ref_Seth2013_2286, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAMY | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_lcy_pbe_init, NULL,
  NULL, NULL, NULL
};
