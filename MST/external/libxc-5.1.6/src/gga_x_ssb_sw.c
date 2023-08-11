/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_SSB_SW       90  /* Swart, Sola and Bickelhaupt correction to PBE  */
#define XC_GGA_X_SSB          91  /* Swart, Sola and Bickelhaupt  */
#define XC_GGA_X_SSB_D        92  /* Swart, Sola and Bickelhaupt dispersion  */
#define XC_GGA_X_REVSSB_D    312  /* Revised Swart, Sola and Bickelhaupt dispersion  */

typedef struct{
  double A, B, C, D, E;
} gga_x_ssb_sw_params;

static void
gga_x_ssb_sw_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ssb_sw_params));
}

#define N_PAR 5
static const char  *names[N_PAR]  = {"_A", "_B", "_C", "_D", "_E"};
static const char  *desc[N_PAR]   = {
  "A, Constant s limit",
  "B in B s^2/(1 + C s^2)",
  "C in B s^2/(1 + C s^2)",
  "D in D s^2/(1 + E s^4)",
  "E in D s^2/(1 + E s^4)"
};
static const double par_ssb_sw[N_PAR] =
  {1.0515, 0.191458, 0.254443, 0.180708, 4.036674};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ssb_sw.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ssb_sw = {
  XC_GGA_X_SSB_SW,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt correction to PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_69, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_ssb_sw, set_ext_params_cpy},
  gga_x_ssb_sw_init, NULL,
  NULL, work_gga, NULL
};

#define SSB_N_PAR 8
static const char  *ssb_names[SSB_N_PAR]  = {"_A", "_B", "_C", "_D", "_E", "_F", "_u", "_delta"};
static const char  *ssb_desc[SSB_N_PAR]   = {
  "A, Constant s limit",
  "B in B s^2/(1 + C s^2)",
  "C in B s^2/(1 + C s^2)",
  "D in D s^2/(1 + E s^4)",
  "E in D s^2/(1 + E s^4)",
  "F, prefactor for KT term",
  "u, reweighting of KT and SSB terms",
  "delta, KT parameter"
};

static const double par_ssb[SSB_N_PAR] =
  {1.071769, 0.137574, 0.187883, 0.137574, 6.635315, 0.995010, -1.205643, 0.1};
static const double par_ssb_d[SSB_N_PAR] =
  {1.079966, 0.197465, 0.272729, 0.197465, 5.873645, 0.949488, -0.749940, 0.1};
static const double par_revssb_d[SSB_N_PAR] =
  {1.082138, 0.177998, 0.246582, 0.177998, 6.284673, 1.0, -0.618168, 0.1};

static void
gga_x_ssb_init(xc_func_type *p) {
  static int    funcs_id  [3] = {XC_LDA_X, XC_GGA_X_SSB_SW, XC_GGA_X_KT1};
  static double funcs_coef[3] = {-1.0, 1.0, 1.0};
  xc_mix_init(p, 3, funcs_id, funcs_coef);
}

static void
ssb_set_ext_params(xc_func_type *p, const double *ext_params) {
  double A, B, C, D, E, F, u, delta;
  double tmp_ssb[5];
  double tmp_kt[2];

  /* Get the parameters */
  A = get_ext_param(p, ext_params, 0);
  B = get_ext_param(p, ext_params, 1);
  C = get_ext_param(p, ext_params, 2);
  D = get_ext_param(p, ext_params, 3);
  E = get_ext_param(p, ext_params, 4);
  F = get_ext_param(p, ext_params, 5);
  u = get_ext_param(p, ext_params, 6);
  delta = get_ext_param(p, ext_params, 7);

  /* Calculate the parameters for SSB-SW and KT1 */
  tmp_ssb[0] = A;
  tmp_ssb[1] = B;
  tmp_ssb[2] = C;
  tmp_ssb[3] = D * (1.0-u);
  tmp_ssb[4] = E;
  /* Swart et al define KT in terms of s^2 not x^2  */
  tmp_kt[0] = u * F * B * (X2S*X2S) * X_FACTOR_C;
  tmp_kt[1] = delta;

  /* Set the SSB and KT parameters */
  xc_func_set_ext_params(p->func_aux[1], tmp_ssb);
  xc_func_set_ext_params(p->func_aux[2], tmp_kt);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ssb = {
  XC_GGA_X_SSB,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_094103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {SSB_N_PAR, ssb_names, ssb_desc, par_ssb, ssb_set_ext_params},
  gga_x_ssb_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ssb_d = {
  XC_GGA_X_SSB_D,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt dispersion",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_094103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {SSB_N_PAR, ssb_names, ssb_desc, par_ssb_d, ssb_set_ext_params},
  gga_x_ssb_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_revssb_d = {
  XC_GGA_X_REVSSB_D,
  XC_EXCHANGE,
  "Revised Swart, Sola and Bickelhaupt dispersion",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2011_1117, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {SSB_N_PAR, ssb_names, ssb_desc, par_revssb_d, ssb_set_ext_params},
  gga_x_ssb_init, NULL,
  NULL, NULL, NULL
};
