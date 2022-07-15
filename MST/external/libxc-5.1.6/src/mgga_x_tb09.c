/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_BJ06         207 /* Becke & Johnson correction to Becke-Roussel 89  */
#define XC_MGGA_X_TB09         208 /* Tran & Blaha correction to Becke & Johnson  */
#define XC_MGGA_X_RPP09        209 /* Rasanen, Pittalis, and Proetto correction to Becke & Johnson  */

typedef struct{
  double c;
  double alpha;
} mgga_x_tb09_params;

static void
mgga_x_tb09_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_tb09_params));
}

#define XC_NO_EXC
#include "decl_mgga.h"
#include "maple2c/mgga_vxc/mgga_x_tb09.c"
#include "work_mgga.c"

#define BJ06_N_PAR 2
static const char  *bj06_names[BJ06_N_PAR]  = {"c", "_alpha"};
static const char  *bj06_desc[BJ06_N_PAR]   = {
  "This parameter involves an average over the unit cell and must be calculated by the calling program.",
  "alpha = 0 for BJ06 and 1 for RPP"
};
static const double bj06_values[BJ06_N_PAR]  = {1.0, 0.0};
static const double rpp09_values[BJ06_N_PAR] = {1.0, 1.0};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_bj06 = {
  XC_MGGA_X_BJ06,
  XC_EXCHANGE,
  "Becke & Johnson 06",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke2006_221101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-15,
  {BJ06_N_PAR, bj06_names, bj06_desc, bj06_values, set_ext_params_cpy},
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_tb09 = {
  XC_MGGA_X_TB09,
  XC_EXCHANGE,
  "Tran & Blaha 09",
  XC_FAMILY_MGGA,
  {&xc_ref_Tran2009_226401, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-23,
  {BJ06_N_PAR, bj06_names, bj06_desc, bj06_values, set_ext_params_cpy},
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_rpp09 = {
  XC_MGGA_X_RPP09,
  XC_EXCHANGE,
  "Rasanen, Pittalis & Proetto 09",
  XC_FAMILY_MGGA,
  {&xc_ref_Rasanen2010_044112, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-15,
  {BJ06_N_PAR, bj06_names, bj06_desc, rpp09_values, set_ext_params_cpy},
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga,
};


