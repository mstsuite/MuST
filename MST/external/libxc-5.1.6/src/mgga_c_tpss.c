/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_TPSS          231 /* Tao, Perdew, Staroverov & Scuseria correlation */
#define XC_MGGA_C_TM            251 /* Tao and Mo 2016 correlation */

typedef struct{
  double beta, d;
  double C0_c[4];
} mgga_c_tpss_params;

static void
mgga_c_tpss_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_tpss_params));
}

#define TPSS_N_PAR 6
static const char  *tpss_names[TPSS_N_PAR]  = {"_beta", "_d", "_C0_c0", "_C0_c1", "_C0_c2", "_C0_c3"};
static const char  *tpss_desc[TPSS_N_PAR]   = {"beta", "d", "C0_c0", "C0_c1", "C0_c2", "C0_c3"};
static const double tpss_values[TPSS_N_PAR] = {
  0.06672455060314922, 2.8, 0.53, 0.87, 0.50, 2.26
};
static const double tm_values[TPSS_N_PAR]   = {
  0.06672455060314922, 2.8, 0.0, 0.1, 0.32, 0.0
};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_tpss.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_tpss = {
  XC_MGGA_C_TPSS,
  XC_CORRELATION,
  "Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Tao2003_146401, &xc_ref_Perdew2004_6898, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15, /* densities smaller than 1e-26 give NaNs */
  {TPSS_N_PAR, tpss_names, tpss_desc, tpss_values, set_ext_params_cpy},
  mgga_c_tpss_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_tm = {
  XC_MGGA_C_TM,
  XC_CORRELATION,
  "Tao and Mo 2016 correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Tao2016_073001, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15, /* densities smaller than 1e-26 give NaNs */
  {TPSS_N_PAR, tpss_names, tpss_desc, tm_values, set_ext_params_cpy},
  mgga_c_tpss_init, NULL,
  NULL, NULL, work_mgga,
};
