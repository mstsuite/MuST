/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MS0          221 /* MS exchange of Sun, Xiao, and Ruzsinszky */
#define XC_MGGA_X_MS1          222 /* MS1 exchange of Sun, et al */
#define XC_MGGA_X_MS2          223 /* MS2 exchange of Sun, et al */
#define XC_HYB_MGGA_X_MS2H     224 /* MS2 hybrid exchange of Sun, et al */
#define XC_MGGA_X_MS2_REV      228 /* MS2 exchange of Sun, et al with a revised value for c  */

typedef struct{
  double kappa, c, b;
} mgga_x_ms_params;

static void
mgga_x_ms_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_ms_params));
}

#define MS0_N_PAR 3
static const char  *ms0_names[MS0_N_PAR]  = {"_kappa", "_c", "_b"};
static const char  *ms0_desc[MS0_N_PAR]   = {
  "kappa parameter",
  "c parameter",
  "exponent b"
};
static const double ms0_values[MS0_N_PAR]  = {0.29, 0.28771, 1.0};
static const double ms1_values[MS0_N_PAR] = {0.404, 0.18150, 1.0};
static const double ms2_values[MS0_N_PAR] = {0.504, 0.14601, 4.0};
static const double ms2rev_values[MS0_N_PAR] = {0.504, 0.14607, 4.0};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_ms.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms0 = {
  XC_MGGA_X_MS0,
  XC_EXCHANGE,
  "MS exchange of Sun, Xiao, and Ruzsinszky",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2012_051101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {MS0_N_PAR, ms0_names, ms0_desc, ms0_values, set_ext_params_cpy},
  mgga_x_ms_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms1 = {
  XC_MGGA_X_MS1,
  XC_EXCHANGE,
  "MS1 exchange of Sun, et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {MS0_N_PAR, ms0_names, ms0_desc, ms1_values, set_ext_params_cpy},
  mgga_x_ms_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms2 = {
  XC_MGGA_X_MS2,
  XC_EXCHANGE,
  "MS2 exchange of Sun, et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {MS0_N_PAR, ms0_names, ms0_desc, ms2_values, set_ext_params_cpy},
  mgga_x_ms_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms2_rev = {
  XC_MGGA_X_MS2_REV,
  XC_EXCHANGE,
  "MS2 exchange of Sun, et al with revised value for c",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2013_044113, &xc_ref_Furness2019_041119, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {MS0_N_PAR, ms0_names, ms0_desc, ms2rev_values, set_ext_params_cpy},
  mgga_x_ms_init, NULL,
  NULL, NULL, work_mgga,
};

static void
hyb_mgga_x_ms2h_init(xc_func_type *p)
{
  static int   funcs_id  [1] = {XC_MGGA_X_MS2};
  static double funcs_coef[1] = {0.91};

  xc_mix_init(p, 1, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.09);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_ms2h = {
  XC_HYB_MGGA_X_MS2H,
  XC_EXCHANGE,
  "MS2 hybrid exchange of Sun, et al",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  hyb_mgga_x_ms2h_init, NULL,
  NULL, NULL, NULL
};
