/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_K_CSK1     629 /* mGGA-rev functional by Cancio, Stewart, and Kuna (a=1) */
#define XC_MGGA_K_CSK4     630 /* mGGA-rev functional by Cancio, Stewart, and Kuna (a=4) */

typedef struct{
  double csk_a;
} mgga_k_csk_params;

#define CSK_N_PAR 1
static const char  *csk_names  [CSK_N_PAR] = {"_a"};
static const char  *csk_desc   [CSK_N_PAR] = {"exponent used in the interpolation"};
static const double csk1_values[CSK_N_PAR] = {1.0};
static const double csk4_values[CSK_N_PAR] = {4.0};

static void
mgga_k_csk_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_k_csk_params));
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_k_csk.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_csk1 = {
  XC_MGGA_K_CSK1,
  XC_KINETIC,
  "mGGA-rev functional by Cancio, Stewart, and Kuna (a=1)",
  XC_FAMILY_MGGA,
  {&xc_ref_Cancio2016_084107, NULL, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {CSK_N_PAR, csk_names, csk_desc, csk1_values, set_ext_params_cpy},
  mgga_k_csk_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_csk4 = {
  XC_MGGA_K_CSK4,
  XC_KINETIC,
  "mGGA-rev functional by Cancio, Stewart, and Kuna (a=4)",
  XC_FAMILY_MGGA,
  {&xc_ref_Cancio2016_084107, NULL, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {CSK_N_PAR, csk_names, csk_desc, csk4_values, set_ext_params_cpy},
  mgga_k_csk_init, NULL,
  NULL, NULL, work_mgga
};
