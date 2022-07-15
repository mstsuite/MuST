/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/************************************************************************
 Implements Tao, Perdew, Staroverov & Scuseria
   meta-Generalized Gradient Approximation.

  Exchange part
************************************************************************/

#define XC_MGGA_X_TPSS          202 /* Tao, Perdew, Staroverov & Scuseria exchange */
#define XC_MGGA_X_MODTPSS       245 /* Modified Tao, Perdew, Staroverov & Scuseria exchange */
#define XC_MGGA_X_REVTPSS       212 /* revised Tao, Perdew, Staroverov & Scuseria exchange */
#define XC_MGGA_X_BLOC          244 /* functional with balanced localization */

typedef struct{
  double b, c, e, kappa, mu;
  double BLOC_a, BLOC_b;
} mgga_x_tpss_params;


static void
mgga_x_tpss_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_tpss_params));
}

#define TPSS_N_PAR 7
static const char  *tpss_names[TPSS_N_PAR]  = {"_b", "_c", "_e", "_kappa", "_mu", "_BLOC_a", "_BLOC_b"};
static const char  *tpss_desc[TPSS_N_PAR]   = {
  "b", "c", "e",
  "Asymptotic value of the enhancement function",
  "Coefficient of the 2nd order expansion",
  "BLOC_a", "BLOC_b"
};
static const double tpss_values[TPSS_N_PAR] = {
  0.40, 1.59096, 1.537, 0.8040, 0.21951, 2.0, 0.0
};
static const double modtpss_values[TPSS_N_PAR] = {
  0.40, 1.38496, 1.37, 0.804, 0.252, 2.0, 0.0
};
static const double revtpss_values[TPSS_N_PAR] = {
  0.40, 2.35203946, 2.16769874, 0.804, 0.14, 3.0, 0.0
};
static const double bloc_values[TPSS_N_PAR]    = {
  0.40, 1.59096, 1.537, 0.804, 0.21951, 4.0, -3.3
};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_tpss.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_tpss = {
  XC_MGGA_X_TPSS,
  XC_EXCHANGE,
  "Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Tao2003_146401, &xc_ref_Perdew2004_6898, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TPSS_N_PAR, tpss_names, tpss_desc, tpss_values, set_ext_params_cpy},
  mgga_x_tpss_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_modtpss = {
  XC_MGGA_X_MODTPSS,
  XC_EXCHANGE,
  "Modified Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2007_042506, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TPSS_N_PAR, tpss_names, tpss_desc, modtpss_values, set_ext_params_cpy},
  mgga_x_tpss_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_revtpss = {
  XC_MGGA_X_REVTPSS,
  XC_EXCHANGE,
  "revised Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2009_026403, &xc_ref_Perdew2009_026403_err, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TPSS_N_PAR, tpss_names, tpss_desc, revtpss_values, set_ext_params_cpy},
  mgga_x_tpss_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_bloc = {
  XC_MGGA_X_BLOC,
  XC_EXCHANGE,
  "functional with balanced localization",
  XC_FAMILY_MGGA,
  {&xc_ref_Constantin2013_2256, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TPSS_N_PAR, tpss_names, tpss_desc, bloc_values, set_ext_params_cpy},
  mgga_x_tpss_init, NULL,
  NULL, NULL, work_mgga,
};
