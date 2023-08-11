/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_LKT  613 /* Luo-Karasiev-Trickey kinetic GGA */

typedef struct{
  double a;
} gga_k_lkt_params;

#define N_PAR 1
static const char  *names[N_PAR]  = {"_a"};
static const char  *desc[N_PAR]   = {"a"};
static const double lkt_val[N_PAR] = {1.3};

static void
gga_k_lkt_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_lkt_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_lkt.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_lkt = {
  XC_GGA_K_LKT,
  XC_KINETIC,
  "Luo-Karasiev-Trickey GGA kinetic",
  XC_FAMILY_GGA,
  {&xc_ref_Luo2018_041111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, lkt_val, set_ext_params_cpy},
  gga_k_lkt_init, NULL,
  NULL, work_gga, NULL
};
