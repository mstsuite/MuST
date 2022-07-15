/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_R2SCAN         497 /* Re-regularized SCAN exchange */

typedef struct{
  double c1, c2, d, k1;
  double taur, alphar;
  double eta, dp2;
} mgga_x_r2scan_params;

#define N_PAR 8
static const char *names[N_PAR] = {"_c1", "_c2", "_d", "_k1", "_taur", "_alphar", "_eta", "_dp2"};
static const char *desc[N_PAR] = {"c1 parameter", "c2 parameter", "d parameter",
  "k1 parameter", "taur parameter", "alphar parameter", "eta parameter", "dp2 parameter"};

static const double par_r2scan[N_PAR] = {0.667, 0.8, 1.24, 0.065, 1.0e-4, 1.0e-3, 0.001, 0.361};

static void
mgga_x_r2scan_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_r2scan_params));
}


#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_r2scan.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_r2scan = {
  XC_MGGA_X_R2SCAN,
  XC_EXCHANGE,
  "Re-regularized SCAN exchange by Furness et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2020_8208, &xc_ref_Furness2020_9248, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-11,
  {N_PAR, names, desc, par_r2scan, set_ext_params_cpy},
  mgga_x_r2scan_init, NULL,
  NULL, NULL, work_mgga
};
