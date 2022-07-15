/*
 Copyright (C) 2006-2007 M.A.L. Marques
                    2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_O3LYP   404 /* hybrid using the optx functional */
#define XC_HYB_GGA_XC_X3LYP   411 /* hybrid by Xu and Goddard */


/*************************************************************/
#define N_PAR_O3LYP 4
static const char *names_o3lyp[N_PAR_O3LYP] = {"_a", "_b", "_c", "_clyp"};
static const char *desc_o3lyp[N_PAR_O3LYP] = {
  "fraction of HF exchange",
  "fraction of LDA exchage",
  "fraction of OPTX gradient correction",
  "fraction of LYP correlation"
};
static const double par_o3lyp[N_PAR_O3LYP] = {0.1161, 0.9262, 0.8133, 0.81};

static void
o3lyp_set_ext_params(xc_func_type *p, const double *ext_params)
{
  /* This is the fraction of LDA exchange in OPTX */
  const double a1 = 1.05151;
  double a, b, c, clyp;

  assert(p != NULL);
  a    = get_ext_param(p, ext_params, 0);
  b    = get_ext_param(p, ext_params, 1);
  c    = get_ext_param(p, ext_params, 2);
  clyp = get_ext_param(p, ext_params, 3);

  /* Remove double counting of LDA exchange */
  p->mix_coef[0] = b - a1*c;
  p->mix_coef[1] = c;
  p->mix_coef[2] = 1.0 - clyp;
  p->mix_coef[3] = clyp;

  p->cam_alpha = a;
}

static void
hyb_gga_xc_o3lyp_init(xc_func_type *p)
{
  static int funcs_id[4] = {XC_LDA_X, XC_GGA_X_OPTX, XC_LDA_C_VWN, XC_GGA_C_LYP};
  double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0};
  xc_mix_init(p, 4, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_o3lyp = {
  XC_HYB_GGA_XC_O3LYP,
  XC_EXCHANGE_CORRELATION,
  "O3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Cohen2001_607, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR_O3LYP, names_o3lyp, desc_o3lyp, par_o3lyp, o3lyp_set_ext_params},
  hyb_gga_xc_o3lyp_init,
  NULL, NULL, NULL, NULL
};


/*************************************************************/
#define N_PAR_X3LYP 5
static const char *names_x3lyp[N_PAR_X3LYP] = {"_a0", "_ax", "_ac", "_ax1", "_ax2"};
static const char *desc_x3lyp[N_PAR_X3LYP] = {
  "fraction of HF exchange",
  "fraction of XLYP gradient correction",
  "fraction of VWN correction",
  "weight of B88 enhancement in XLYP exchange",
  "weight of PW91 enhancement in XLYP exchange"
};
static const double par_x3lyp[N_PAR_X3LYP] = {0.218, 0.709, 0.871, 0.765, 0.235};

static void
x3lyp_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double a0, ax, ac, ax1, ax2;

  assert(p != NULL);
  a0    = get_ext_param(p, ext_params, 0);
  ax    = get_ext_param(p, ext_params, 1);
  ac    = get_ext_param(p, ext_params, 2);
  ax1   = get_ext_param(p, ext_params, 3);
  ax2   = get_ext_param(p, ext_params, 4);

  p->mix_coef[0] = 1.0 - a0 - ax*(ax1 + ax2);
  p->mix_coef[1] = ax*ax1;
  p->mix_coef[2] = ax*ax2;
  p->mix_coef[3] = 1.0 - ac;
  p->mix_coef[4] = ac;

  p->cam_alpha = a0;
}

static void
hyb_gga_xc_x3lyp_init(xc_func_type *p)
{
  static int funcs_id[5] = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_X_PW91, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  double funcs_coef[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

  xc_mix_init(p, 5, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_x3lyp = {
  XC_HYB_GGA_XC_X3LYP,
  XC_EXCHANGE_CORRELATION,
  "X3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Xu2004_2673, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR_X3LYP, names_x3lyp, desc_x3lyp, par_x3lyp, x3lyp_set_ext_params},
  hyb_gga_xc_x3lyp_init,
  NULL, NULL, NULL, NULL
};
