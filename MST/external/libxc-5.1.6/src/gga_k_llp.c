/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_LLP          522 /* Lee, Lee & Parr */
#define XC_GGA_K_FR_B88       514 /* Fuentealba & Reyes (B88 version) */

typedef struct{
  double beta, gamma;
} gga_k_llp_params;


static void
gga_k_llp_init(xc_func_type *p)
{
  gga_k_llp_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_llp_params));
  params = (gga_k_llp_params *) (p->params);

  /* value of beta in standard Becke 88 functional */
  switch(p->info->number){
  case XC_GGA_K_LLP:
    params->beta  = X_FACTOR_C*0.0044188;
    params->gamma = 0.0253/(X_FACTOR_C*0.0044188);
    break;
  case XC_GGA_K_FR_B88:
    params->beta  = X_FACTOR_C*0.004596;
    params->gamma = 0.02774/(X_FACTOR_C*0.004596);
    break;
  default:
    fprintf(stderr, "Internal error in gga_k_llp\n");
    exit(1);
  }
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_llp.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_llp = {
  XC_GGA_K_LLP,
  XC_KINETIC,
  "Lee, Lee & Parr",
  XC_FAMILY_GGA,
  {&xc_ref_Lee1991_768, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_llp_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_fr_b88 = {
  XC_GGA_K_FR_B88,
  XC_KINETIC,
  "Fuentealba & Reyes (B88 version)",
  XC_FAMILY_GGA,
  {&xc_ref_Fuentealba1995_31, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_llp_init, NULL,
  NULL, work_gga, NULL
};
