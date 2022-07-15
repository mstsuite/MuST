/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_CAM_O3LYP   395 /* range separated hybrid based on the optx functional */

#define N_PAR 6
static const char *names[N_PAR] = {"_csr", "_b", "_c", "_clyp", "_clr", "_omega"};
static const char *desc[N_PAR] = {
  "fraction of short-range HF exchange",
  "fraction of LDA exchage",
  "fraction of OPTX gradient correction",
  "fraction of LYP correlation",
  "fraction of long-range HF exchange",
  "range separation parameter"
};
/* Functional is defined in main text in section 5.3 as 11%
   short-range as in O3LYP, and 80% long-range exact exchange as in
   one variant of CAM-B3LYP. However, as O3LYP actually has 11.61% of
   exact exchange, we use the more accurate value to minimize the
   difference to O3LYP.
 */
static const double par_cam_o3lyp[N_PAR] = {0.1161, 0.9262, 0.8133, 0.81, 0.80, 0.33};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  /* This is the fraction of LDA exchange in OPTX */
  const double a1 = 1.05151;
  double csr, b, c, clyp, clr, omega;

  assert(p != NULL);
  csr   = get_ext_param(p, ext_params, 0);
  b     = get_ext_param(p, ext_params, 1);
  c     = get_ext_param(p, ext_params, 2);
  clyp  = get_ext_param(p, ext_params, 3);
  clr   = get_ext_param(p, ext_params, 4);
  omega = get_ext_param(p, ext_params, 5);

  /* Remove double counting of LDA exchange */
  p->mix_coef[0] = b - a1*c;
  p->mix_coef[1] = c;
  p->mix_coef[2] = 1.0 - clyp;
  p->mix_coef[3] = clyp;

  /* Set range separation */
  xc_func_set_ext_params_name(p->func_aux[0], "_omega", omega);
  xc_func_set_ext_params_name(p->func_aux[1], "_omega", omega);

  /* Set hybrid terms */
  p->cam_beta = csr-clr;
  p->cam_omega = omega;
  p->cam_alpha = clr;
}

static void
hyb_gga_xc_cam_o3lyp_init(xc_func_type *p)
{
  static int funcs_id[4] = {XC_LDA_X_ERF, XC_GGA_X_ITYH_OPTX, XC_LDA_C_VWN, XC_GGA_C_LYP};
  double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0};
  xc_mix_init(p, 4, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_o3lyp = {
  XC_HYB_GGA_XC_CAM_O3LYP,
  XC_EXCHANGE_CORRELATION,
  "CAM-O3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Bircher2018_3184, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, par_cam_o3lyp, set_ext_params},
  hyb_gga_xc_cam_o3lyp_init,
  NULL, NULL, NULL, NULL
};

