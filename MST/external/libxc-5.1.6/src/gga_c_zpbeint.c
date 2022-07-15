/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_ZPBEINT       61 /* spin-dependent gradient correction to PBEint       */
#define XC_GGA_C_ZPBESOL       63 /* spin-dependent gradient correction to PBEsol       */

typedef struct{
  double beta, alpha;
} gga_c_zpbeint_params;

static void
gga_c_zpbeint_init(xc_func_type *p)
{
  gga_c_zpbeint_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_zpbeint_params));
  params = (gga_c_zpbeint_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_C_ZPBEINT:
    params->beta  = 0.052;
    params->alpha = 2.4;
    break;
  case XC_GGA_C_ZPBESOL:
    params->beta  = 0.046;
    params->alpha = 4.8;
    break;
  default:
    fprintf(stderr, "Internal error in gga_c_zpbeint\n");
    exit(1);
  }
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_zpbeint.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_zpbeint = {
  XC_GGA_C_ZPBEINT,
  XC_CORRELATION,
  "spin-dependent gradient correction to PBEint",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_233103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {0, NULL, NULL, NULL, NULL},
  gga_c_zpbeint_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_zpbesol = {
  XC_GGA_C_ZPBESOL,
  XC_CORRELATION,
  "spin-dependent gradient correction to PBEsol",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_233103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {0, NULL, NULL, NULL, NULL},
  gga_c_zpbeint_init, NULL,
  NULL, work_gga, NULL
};
