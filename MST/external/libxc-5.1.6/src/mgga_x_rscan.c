/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_RSCAN         493 /* Regularized SCAN exchange */

typedef struct{
  double c1, c2, d, k1;
  double taur, alphar;
} mgga_x_rscan_params;

#define N_PAR 6
static const char *names[N_PAR] = {"_c1", "_c2", "_d", "_k1", "_taur", "_alphar"};
static const char *desc[N_PAR] = {"c1 parameter", "c2 parameter", "d parameter",
                                  "k1 parameter", "taur parameter", "alphar parameter"};

static const double par_rscan[N_PAR] = {0.667, 0.8, 1.24, 0.065, 1.0e-4, 1.0e-3};

static void
mgga_x_rscan_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_rscan_params));
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_rscan.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_rscan = {
  XC_MGGA_X_RSCAN,
  XC_EXCHANGE,
  "Regularized SCAN exchange by Bartok and Yates",
  XC_FAMILY_MGGA,
  {&xc_ref_Bartok2019_161101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-11,
  {N_PAR, names, desc, par_rscan, set_ext_params_cpy},
  mgga_x_rscan_init, NULL,
  NULL, NULL, work_mgga,
};
