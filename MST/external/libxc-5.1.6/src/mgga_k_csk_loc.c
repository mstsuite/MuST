/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_K_CSK_LOC1   631 /* mGGAloc-rev functional by Cancio, Stewart, and Kuna (a=1) */
#define XC_MGGA_K_CSK_LOC4   632 /* mGGAloc-rev functional by Cancio, Stewart, and Kuna (a=4) */

typedef struct{
  double csk_a, csk_cp, csk_cq;
} mgga_k_csk_loc_params;

#define CSK_N_PAR 3
static const char  *csk_names  [CSK_N_PAR] = {"_a", "_cp", "_cq"};
static const char  *csk_desc   [CSK_N_PAR] =
  {"exponent used in the interpolation", "coefficient of the p term", "coefficient of the q term"};
static const double csk1_values[CSK_N_PAR] = {1.0, -0.275, 2.895};
static const double csk4_values[CSK_N_PAR] = {4.0, -0.275, 2.895};

static void
mgga_k_csk_loc_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_k_csk_loc_params));
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_k_csk_loc.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_csk_loc1 = {
  XC_MGGA_K_CSK_LOC1,
  XC_KINETIC,
  "mGGAloc-rev functional by Cancio, Stewart, and Kuna (a=1)",
  XC_FAMILY_MGGA,
  {&xc_ref_Cancio2016_084107, NULL, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {CSK_N_PAR, csk_names, csk_desc, csk1_values, set_ext_params_cpy},
  mgga_k_csk_loc_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_csk_loc4 = {
  XC_MGGA_K_CSK_LOC4,
  XC_KINETIC,
  "mGGAloc-rev functional by Cancio, Stewart, and Kuna (a=4)",
  XC_FAMILY_MGGA,
  {&xc_ref_Cancio2016_084107, NULL, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {CSK_N_PAR, csk_names, csk_desc, csk4_values, set_ext_params_cpy},
  mgga_k_csk_loc_init, NULL,
  NULL, NULL, work_mgga
};
