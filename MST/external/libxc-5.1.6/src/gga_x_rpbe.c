/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_RPBE  117 /* Hammer, Hansen & Norskov (PBE-like) */


typedef struct{
  double rpbe_kappa, rpbe_mu;
} gga_x_rpbe_params;


static void
gga_x_rpbe_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_rpbe_params));
}

#define RPBE_N_PAR 2
static const char  *rpbe_names[RPBE_N_PAR]  = {"_kappa", "_mu"};
static const char  *rpbe_desc[RPBE_N_PAR]   = {
  "Asymptotic value of the enhancement function",
  "Coefficient of the 2nd order expansion"};
static const double rpbe_values[RPBE_N_PAR] =
  {0.8040, MU_PBE};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_rpbe.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_rpbe = {
  XC_GGA_X_RPBE,
  XC_EXCHANGE,
  "Hammer, Hansen, and Norskov",
  XC_FAMILY_GGA,
  {&xc_ref_Hammer1999_7413, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {RPBE_N_PAR, rpbe_names, rpbe_desc, rpbe_values, set_ext_params_cpy},
  gga_x_rpbe_init, NULL,
  NULL, work_gga, NULL
};
