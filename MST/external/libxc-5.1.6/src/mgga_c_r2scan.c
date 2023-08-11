/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_R2SCAN         498 /* Re-regularized SCAN correlation */
#define XC_MGGA_C_R2SCAN01       642 /* Re-regularized SCAN correlation with larger value for eta */

typedef struct{
  double eta;   /* regularization parameter */
} mgga_c_r2scan_params;

#define N_PAR 1
static const char  *names[N_PAR]  = {"_eta"};
static const char  *desc[N_PAR]   = {
  "Regularization parameter"};

static const double r2scan_values[N_PAR] = {0.001};
static const double r2scan01_values[N_PAR] = {0.01};
  
#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_r2scan.c"
#include "work_mgga.c"

static void
mgga_c_r2scan_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_r2scan_params));
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_r2scan = {
  XC_MGGA_C_R2SCAN,
  XC_CORRELATION,
  "Re-regularized SCAN correlation by Furness et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2020_8208, &xc_ref_Furness2020_9248, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, r2scan_values, set_ext_params_cpy},
  mgga_c_r2scan_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_r2scan01 = {
  XC_MGGA_C_R2SCAN01,
  XC_CORRELATION,
  "Re-regularized SCAN correlation with larger value for eta",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2020_8208, &xc_ref_Furness2020_9248, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, r2scan01_values, set_ext_params_cpy},
  mgga_c_r2scan_init, NULL,
  NULL, NULL, work_mgga,
};
