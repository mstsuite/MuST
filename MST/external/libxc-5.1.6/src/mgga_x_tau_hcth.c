/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_TAU_HCTH        205 /* tau-HCTH from Boese and Handy */
#define XC_HYB_MGGA_X_BMK         279 /* Boese-Martin for kinetics     */
#define XC_HYB_MGGA_X_TAU_HCTH    282 /* Hybrid version of tau-HCTH    */

#define N_PAR_PURE 8
static const char  *pure_names[N_PAR_PURE]  = {"_cxl0", "_cxl1", "_cxl2", "_cxl3", "_cxnl0", "_cxnl1", "_cxnl2", "_cxnl3"};
static const char  *pure_desc[N_PAR_PURE]   = {"Local exchange, u^0 coefficient", "Local exchange, u^1 coefficient", "Local exchange, u^2 coefficient", "Local exchange, u^3 coefficient", "Non-local exchange, u^0 coefficient", "Non-local exchange, u^1 coefficient", "Non-local exchange, u^2 coefficient", "Non-local exchange, u^3 coefficient"};

#define N_PAR_HYB 9
static const char  *hyb_names[N_PAR_HYB]  = {"_cxl0", "_cxl1", "_cxl2", "_cxl3", "_cxnl0", "_cxnl1", "_cxnl2", "_cxnl3", "_ax"};
static const char  *hyb_desc[N_PAR_HYB]   = {"Local exchange, u^0 coefficient", "Local exchange, u^1 coefficient", "Local exchange, u^2 coefficient", "Local exchange, u^3 coefficient", "Non-local exchange, u^0 coefficient", "Non-local exchange, u^1 coefficient", "Non-local exchange, u^2 coefficient", "Non-local exchange, u^3 coefficient", "Fraction of exact exchange"};

const double tHCTH_val [N_PAR_PURE] = {1.10734, -1.0534, 6.3491, -2.5531, 0.00110, -0.3041, 6.9543, -0.7235};
const double BMK_val [N_PAR_HYB] = { 0.474302, 2.77701, -11.4230, 11.7167, -0.192212, 4.73936, -26.6188, 22.4891, 0.42};
const double hyb_tHCTH_val [N_PAR_HYB] = { 0.86735,  0.3008, 1.2208,   0.1574, -0.00230, -0.2849, 5.4146, -10.909, 0.15};

typedef struct{
  double cx_local[4];
  double cx_nlocal[4];
} mgga_x_tau_hcth_params;

static void
mgga_x_tau_hcth_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_tau_hcth_params));
  if(p->info->number == XC_HYB_MGGA_X_BMK || p->info->number == XC_HYB_MGGA_X_TAU_HCTH)
    xc_hyb_init_hybrid(p, 0.0);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_tau_hcth.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_tau_hcth = {
  XC_MGGA_X_TAU_HCTH,
  XC_EXCHANGE,
  "tau-HCTH from Boese and Handy",
  XC_FAMILY_MGGA,
  {&xc_ref_Boese2002_9559, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_PURE, pure_names, pure_desc, tHCTH_val, set_ext_params_cpy},
  mgga_x_tau_hcth_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_bmk = {
  XC_HYB_MGGA_X_BMK,
  XC_EXCHANGE,
  "Boese-Martin for kinetics",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Boese2004_3405, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_HYB, hyb_names, hyb_desc, BMK_val, set_ext_params_cpy_exx},
  mgga_x_tau_hcth_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_tau_hcth = {
  XC_HYB_MGGA_X_TAU_HCTH,
  XC_EXCHANGE,
  "Hybrid version of tau-HCTH",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Boese2002_9559, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_HYB, hyb_names, hyb_desc, hyb_tHCTH_val, set_ext_params_cpy_exx},
  mgga_x_tau_hcth_init,  NULL,
  NULL, NULL, work_mgga,
};
