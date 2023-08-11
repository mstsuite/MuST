/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_LC94         521 /* Lembarki & Chermette */

typedef struct{
  double a, b, c, d, f, alpha, expo;
} gga_k_lc94_params;

#define N_PAR 7
static const char *names[N_PAR]  = {"_a", "_b", "_c", "_d", "_f", "_alpha", "_expo"};
static const char *desc[N_PAR]   = {"a parameter", "b parameter", "c parameter", "d parameter", "f parameter", "alpha parameter", "exponent"};

static const double par_k_lc94[N_PAR] =
  {0.093907, 76.320, 0.26608, -0.0809615, 0.000057767, 100.0, 4.0};

static void
gga_k_lc94_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_lc94_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_lc94.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_lc94 = {
  XC_GGA_K_LC94,
  XC_KINETIC,
  "Lembarki & Chermette",
  XC_FAMILY_GGA,
  {&xc_ref_Lembarki1994_5328, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_k_lc94, set_ext_params_cpy},
  gga_k_lc94_init, NULL,
  NULL, work_gga, NULL
};
