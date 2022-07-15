/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_LDA_XC_CAM_LDA0    178 /* CAM version of LDA0 */

void
xc_hyb_lda_xc_cam_lda0_init(xc_func_type *p)
{
  static int funcs_id[3] = {XC_LDA_X, XC_LDA_X_ERF, XC_LDA_C_PW_MOD};
  double funcs_coef[3];
  double omega, exxglobal, exxlr;
  double alpha, beta;

  /* omega = 1/3 */
  omega = 1.0/3.0;
  /* 1/4 global exact exchange */
  exxglobal = 1.0/4.0;
  /* 1/2 long-range exact exchange */
  exxlr = 1.0/2.0;

  /* Convert to libxc parameters */
  alpha = exxlr;
  beta = exxglobal-exxlr;

  /* Calculate mixing coefficients: full LDA */
  funcs_coef[0] = 1.0 - alpha;
  funcs_coef[1] = -beta;
  funcs_coef[2] = 1.0;

  /* Initialize functional */
  xc_mix_init(p, 3, funcs_id, funcs_coef);

  /* Set parameters */
  xc_func_set_ext_params(p->func_aux[1], &omega);

  xc_hyb_init_cam(p, alpha, beta, omega);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_lda_xc_cam_lda0 = {
  XC_HYB_LDA_XC_CAM_LDA0,
  XC_EXCHANGE_CORRELATION,
  "CAM version of LDA0",
  XC_FAMILY_HYB_LDA,
  {&xc_ref_Mosquera2016_1605, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_lda_xc_cam_lda0_init, NULL,
  NULL, NULL, NULL
};
