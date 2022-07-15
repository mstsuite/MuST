/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_OPTX         110 /* Handy & Cohen OPTX 01                          */

typedef struct{
  double a, b, gamma;
} gga_x_optx_params;


static void
gga_x_optx_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_optx_params));
}

#define OPTX_N_PAR 3
static const char  *optx_names[OPTX_N_PAR]  = {"_a", "_b", "_gamma"};
static const char  *optx_desc[OPTX_N_PAR]   = {
  "a",
  "b",
  "gamma"};
static const double optx_values[OPTX_N_PAR] =
  {1.05151, 1.43169/X_FACTOR_C, 0.006};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_optx.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_optx = {
  XC_GGA_X_OPTX,
  XC_EXCHANGE,
  "Handy & Cohen OPTX 01",
  XC_FAMILY_GGA,
  {&xc_ref_Handy2001_403, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {OPTX_N_PAR, optx_names, optx_desc, optx_values, set_ext_params_cpy},
  gga_x_optx_init, NULL,
  NULL, work_gga, NULL
};
