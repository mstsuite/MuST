/*
 Copyright (C) 2016 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_SCAN          263 /* SCAN exchange of Sun, Ruzsinszky, and Perdew  */
#define XC_HYB_MGGA_X_SCAN0     264 /* SCAN hybrid exchange */
#define XC_MGGA_X_REVSCAN       581 /* revised SCAN */
#define XC_HYB_MGGA_X_REVSCAN0  583 /* revised SCAN hybrid exchange */

typedef struct{
  double c1, c2, d, k1;
} mgga_x_scan_params;

#define N_PAR_SCAN 4
static const char *scan_names[N_PAR_SCAN] = {"_c1", "_c2", "_d", "_k1"};
static const char *scan_desc[N_PAR_SCAN] = {"c1 parameter", "c2 parameter", "d parameter",
                                            "k1 parameter"};

static const double par_scan[N_PAR_SCAN] = {0.667, 0.8, 1.24, 0.065};
static const double par_revscan[N_PAR_SCAN] = {0.607, 0.7, 1.37, 0.065};

static void
mgga_x_scan_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_scan_params));
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_scan.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_scan = {
  XC_MGGA_X_SCAN,
  XC_EXCHANGE,
  "SCAN exchange of Sun, Ruzsinszky, and Perdew",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2015_036402, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_SCAN, scan_names, scan_desc, par_scan, set_ext_params_cpy},
  mgga_x_scan_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_revscan = {
  XC_MGGA_X_REVSCAN,
  XC_EXCHANGE,
  "revised SCAN",
  XC_FAMILY_MGGA,
  {&xc_ref_Mezei2018_2469, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_SCAN, scan_names, scan_desc, par_revscan, set_ext_params_cpy},
  mgga_x_scan_init, NULL,
  NULL, NULL, work_mgga
};

typedef struct{
  double exx;
} hyb_mgga_x_scan0_params;

#define N_PAR_SCAN0 1
static const char *scan0_names[N_PAR_SCAN0] = {"_exx"};
static const char *scan0_desc[N_PAR_SCAN0] = {"fraction of exact exchange"};

static const double scan0_pars[N_PAR_SCAN0] = {0.25};

static void
scan0_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double a0;
  assert(p != NULL);
  a0 = get_ext_param(p, ext_params, 0);

  p->mix_coef[0] = 1.0 - a0;
  p->cam_alpha = a0;
}

static void
hyb_mgga_x_scan0_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_scan_params));

  static int   funcs_id  [1] = {XC_MGGA_X_SCAN};
  static double funcs_coef[1] = {0.0}; /* set by ext_params */

  xc_mix_init(p, 1, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0); /* set by ext_params */
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_scan0 = {
  XC_HYB_MGGA_X_SCAN0,
  XC_EXCHANGE,
  "SCAN hybrid exchange (SCAN0)",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Hui2016_044114, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_SCAN0, scan0_names, scan0_desc, scan0_pars, scan0_set_ext_params},
  hyb_mgga_x_scan0_init, NULL,
  NULL, NULL, NULL /* this is taken care of by the generic routine */
};

static void
hyb_mgga_x_revscan0_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_scan_params));

  static int   funcs_id  [1] = {XC_MGGA_X_REVSCAN};
  static double funcs_coef[1] = {0.0}; /* set by ext_params */

  xc_mix_init(p, 1, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0); /* set by ext_params */
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_revscan0 = {
  XC_HYB_MGGA_X_REVSCAN0,
  XC_EXCHANGE,
  "revised SCAN hybrid exchange (SCAN0)",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Mezei2018_2469, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_SCAN0, scan0_names, scan0_desc, scan0_pars, scan0_set_ext_params},
  hyb_mgga_x_revscan0_init, NULL,
  NULL, NULL, NULL /* this is taken care of by the generic routine */
};
