/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_X_MBR 716 /* modified Becke-Roussel by Patra et al */

typedef struct{
  double gamma, beta, lambda;
} mgga_x_mbr_params;

static void
mgga_x_mbr_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_mbr_params));
}

#define MBR_N_PAR 3
static const char  *mbr_names[MBR_N_PAR]  = {"_gamma", "_beta", "_lambda"};
static const char  *mbr_desc[MBR_N_PAR]   = { "gamma",  "beta",  "lambda"};
static const double mbr_values[MBR_N_PAR] = { 1.0, 20.0, 0.877};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_mbr.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mbr = {
  XC_MGGA_X_MBR,
  XC_EXCHANGE,
  "modified Becke-Roussel by Patra et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Patra2019_19639, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {MBR_N_PAR, mbr_names, mbr_desc, mbr_values, set_ext_params_cpy},
  mgga_x_mbr_init, NULL,
  NULL, NULL, work_mgga,
};
