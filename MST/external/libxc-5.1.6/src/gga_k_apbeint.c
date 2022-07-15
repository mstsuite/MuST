/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_APBEINT       54 /* interpolated version of APBE                   */
#define XC_GGA_K_REVAPBEINT    53 /* interpolated version of REVAPBE                */


typedef struct{
  double kappa, alpha, muPBE, muGE;
} gga_k_apbeint_params;


static void
gga_k_apbe_init(xc_func_type *p)
{
  gga_k_apbeint_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_apbeint_params));
  params = (gga_k_apbeint_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_K_APBEINT:
    params->kappa = 0.8040;
    params->alpha = 5.0/3.0;
    params->muPBE = 0.23899;
    params->muGE  = 5.0/27.0;
    break;
  case XC_GGA_K_REVAPBEINT:
    params->kappa = 1.245;
    params->alpha = 5.0/3.0;
    params->muPBE = 0.23899;
    params->muGE  = 5.0/27.0;
    break;
  default:
    fprintf(stderr, "Internal error in gga_k_apbeint\n");
    exit(1);
  }
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_apbeint.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_apbeint = {
  XC_GGA_K_APBEINT,
  XC_KINETIC,
  "interpolated version of APBE",
  XC_FAMILY_GGA,
  {&xc_ref_Laricchia2011_2439, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_apbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_revapbeint = {
  XC_GGA_K_REVAPBEINT,
  XC_KINETIC,
  "interpolated version of revAPBE",
  XC_FAMILY_GGA,
  {&xc_ref_Laricchia2011_2439, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_apbe_init, NULL,
  NULL, work_gga, NULL
};
