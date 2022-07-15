/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_KT1          145 /* Exchange part of Keal and Tozer version 1 */
#define XC_GGA_XC_KT1         167 /* Keal and Tozer version 1                  */
#define XC_GGA_XC_KT2         146 /* Keal and Tozer version 2                  */
#define XC_GGA_XC_KT3         587 /* Keal and Tozer version 3                  */

typedef struct{
  double gamma, delta;
} gga_x_kt_params;

static void
gga_x_kt_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_kt_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_kt.c"
#include "work_gga.c"

#define KT_N_PAR 2
static const char  *kt_names[KT_N_PAR]  = {"_gamma", "_delta"};
static const char  *kt_desc[KT_N_PAR]   = {
  "gamma",
  "delta"};
static const double kt_values[KT_N_PAR] =
  {-0.006, 0.1};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_kt1 = {
  XC_GGA_X_KT1,
  XC_EXCHANGE,
  "Exchange part of Keal and Tozer version 1",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {KT_N_PAR, kt_names, kt_desc, kt_values, set_ext_params_cpy},
  gga_x_kt_init, NULL,
  NULL, work_gga, NULL
};


static void
gga_xc_kt1_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_KT1, XC_LDA_C_VWN};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_kt1 = {
  XC_GGA_XC_KT1,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 1",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_xc_kt1_init, NULL,
  NULL, NULL, NULL
};


static void
gga_xc_kt2_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_LDA_X, XC_GGA_X_KT1, XC_LDA_C_VWN};
  static double funcs_coef[3] = {1.07173 - 1.0, 1.0, 0.576727};

  xc_mix_init(p, 3, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_kt2 = {
  XC_GGA_XC_KT2,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 2",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_xc_kt2_init, NULL,
  NULL, NULL, NULL
};


static void
gga_xc_kt3_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_C_LYP, XC_GGA_X_KT1, XC_GGA_X_OPTX};
  double funcs_coef[4];

  /* Equation (2) */
  static const double alpha = 1.092; /* full E_x^Dirac */
  static const double beta  = 0.864409; /* full E_c^LYP */
  static const double par_kt[2] = {-0.004, 0.1}; /* parameters gamma and delta of KT1 gradient term */
  static const double eps = +0.925452; /* gradient part of E_x^OPTX; changed to include the minus sign! */

  /* The OPTX part is given in terms of the OPTX gradient part */
  static const double a_optx = 1.05151; /* Coefficient of LDA term in OPTX */
  static const double b_optx = 1.43169; ///X_FACTOR_C; /* Coefficient of GGA term in OPTX */

  /* Total LDA weight: substract doubly counted part in KT1 and OPTX */
  funcs_coef[0] = alpha - eps*a_optx/b_optx - 1.0;
  funcs_coef[1] = beta; /* LYP */
  funcs_coef[2] = 1.0; /* KT1; gamma is already included by the parameter! */
  funcs_coef[3] = eps/b_optx; /* OPTX */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
  set_ext_params_cpy(p->func_aux[2], par_kt); /* kt parameters */
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_kt3 = {
  XC_GGA_XC_KT3,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 3",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2004_5654, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_xc_kt3_init, NULL,
  NULL, NULL, NULL
};
