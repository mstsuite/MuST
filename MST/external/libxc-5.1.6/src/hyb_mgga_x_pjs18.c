/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_MGGA_X_PJS18       706 /* a screened version of TM */
#define XC_HYB_MGGA_XC_LC_TMLYP  720 /* long-range corrected TM-LYP */

static void
hyb_mgga_x_pjs18_init(xc_func_type *p)
{
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

#define PJS18_N_PAR 1
static const char  *pjs18_names[PJS18_N_PAR]  = {"_omega"};
static const char  *pjs18_desc[PJS18_N_PAR]   = {
  "Range separation parameter"
};
static const double par_pjs18[PJS18_N_PAR] = {0.33};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/hyb_mgga_x_pjs18.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_pjs18 = {
  XC_HYB_MGGA_X_PJS18,
  XC_EXCHANGE,
  "Patra, Jana and Samal 2018, screened range-separated TM exchange",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Patra2018_8991, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-14,
  {PJS18_N_PAR, pjs18_names, pjs18_desc, par_pjs18, set_ext_params_cpy_lc},
  hyb_mgga_x_pjs18_init, NULL,
  NULL, NULL, work_mgga
};


#define LC_TMLYP_N_PAR 1
static const char  *lc_tmlyp_names[LC_TMLYP_N_PAR]  = {"_omega"};
static const char  *lc_tmlyp_desc[LC_TMLYP_N_PAR]   = {
  "Range separation parameter"
};
static const double par_lc_tmlyp[LC_TMLYP_N_PAR] = {0.28};

static void
hyb_mgga_xc_lc_tmlyp_init(xc_func_type *p)
{
  static int    funcs_id  [2] = {XC_HYB_MGGA_X_PJS18, XC_GGA_C_LYP};
  static double funcs_coef[2] = {1.0, 1.0};
  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

static void lc_tmlyp_set_ext_params(xc_func_type *p, const double *ext_params) {
  double omega;

  omega = get_ext_param(p, ext_params, 0);

  /* 100% long-range exact exchange */
  p->cam_alpha = 1.0;
  p->cam_beta = -1.0;
  p->cam_omega = omega;

  /* Set the parameters for js18 */
  xc_func_set_ext_params(p->func_aux[0], &omega);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_lc_tmlyp = {
  XC_HYB_MGGA_XC_LC_TMLYP,
  XC_EXCHANGE,
  "Long-range corrected TM-LYP by Jana et al",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Jana2018_1, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-14,
  {LC_TMLYP_N_PAR, lc_tmlyp_names, lc_tmlyp_desc, par_lc_tmlyp, lc_tmlyp_set_ext_params},
  hyb_mgga_xc_lc_tmlyp_init, NULL,
  NULL, NULL, NULL
};
