/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_PBEINT        60 /* PBE for hybrid interfaces                      */

typedef struct{
  double kappa, alpha, muPBE, muGE;
} gga_x_pbeint_params;


static void
gga_x_pbe_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pbeint_params));
}

#define PBEINT_N_PAR 4
static const char  *pbeint_names[PBEINT_N_PAR]  = {"_kappa", "_alpha", "_muPBE", "_muGE"};
static const char  *pbeint_desc[PBEINT_N_PAR]   = {
  "Asymptotic value of the enhancement function",
  "defines the width of the interpolation",
  "Limiting value for large s",
  "Limiting value for small s"
};
static const double pbeint_values[PBEINT_N_PAR] = {0.8040, 0.197, MU_PBE, MU_GE};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_pbeint.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbeint = {
  XC_GGA_X_PBEINT,
  XC_EXCHANGE,
  "PBE for hybrid interfaces",
  XC_FAMILY_GGA,
  {&xc_ref_Fabiano2010_113104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {PBEINT_N_PAR, pbeint_names, pbeint_desc, pbeint_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};
