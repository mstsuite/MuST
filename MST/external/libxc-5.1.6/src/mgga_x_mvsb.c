/*
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MVSB          302 /* MVSBeta exchange of Furness and Sun */
#define XC_MGGA_X_MVSBS         303 /* MVSBeta* exchange of Furness and Sun */

typedef struct {
  double e1, c1, k0, b;
} mgga_x_mvsb_params;

static void
mgga_x_mvsb_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_mvsb_params));
}

#define MVSB_N_PAR 4
static const char  *mvsb_names[MVSB_N_PAR]  = {"_e1", "_c1", "_k0", "_b",};
static const char  *mvsb_desc[MVSB_N_PAR]   = {
  "e1 parameter",
  "c1 parameter",
  "k0 parameter",
  "b parameter"
};
static const double mvsb_values[MVSB_N_PAR]  = {-1.6665, 7.8393, 0.174, 0.0233};
static const double mvsbs_values[MVSB_N_PAR] = {-2.3800, 6.3783, 0.174, 0.0233};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_mvsb.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mvsb = {
  XC_MGGA_X_MVSB,
  XC_EXCHANGE,
  "MVSbeta exchange by Furness and Sun",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2018, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {MVSB_N_PAR, mvsb_names, mvsb_desc, mvsb_values, set_ext_params_cpy},
  mgga_x_mvsb_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mvsbs = {
  XC_MGGA_X_MVSBS,
  XC_EXCHANGE,
  "MVSbeta* exchange by Furness and Sun",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2018, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {MVSB_N_PAR, mvsb_names, mvsb_desc, mvsbs_values, set_ext_params_cpy},
  mgga_x_mvsb_init, NULL,
  NULL, NULL, work_mgga,
};
