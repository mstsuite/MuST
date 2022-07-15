/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_OL2          513 /* Ou-Yang and Levy v.2 */

typedef struct{
  double aa, bb, cc;
} gga_k_ol2_params;

static void
gga_k_ol2_init(xc_func_type *p)
{
  gga_k_ol2_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_ol2_params));
  params = (gga_k_ol2_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_K_OL2:
    params->aa = 1.0;
    params->bb = 1.0/K_FACTOR_C;
    params->cc = 0.0887*M_CBRT4/K_FACTOR_C;
    break;
  }
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_ol2.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_ol2 = {
  XC_GGA_K_OL2,
  XC_KINETIC,
  "Ou-Yang and Levy v.2",
  XC_FAMILY_GGA,
  {&xc_ref_OuYang1991_379, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_ol2_init, NULL,
  NULL, work_gga, NULL
};
