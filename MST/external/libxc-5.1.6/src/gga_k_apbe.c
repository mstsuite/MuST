/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_APBE         185 /* mu fixed from the semiclassical neutral atom   */
#define XC_GGA_K_TW1          187 /* Tran and Wesolowski set 1 (Table II)           */
#define XC_GGA_K_TW2          188 /* Tran and Wesolowski set 2 (Table II)           */
#define XC_GGA_K_TW3          189 /* Tran and Wesolowski set 3 (Table II)           */
#define XC_GGA_K_TW4          190 /* Tran and Wesolowski set 4 (Table II)           */
#define XC_GGA_K_REVAPBE       55 /* revised APBE                                   */


typedef struct{
  double kappa, mu;
  double lambda;   /* parameter used in the Odashima & Capelle versions */
} gga_k_apbe_params;


static void
gga_k_apbe_init(xc_func_type *p)
{
  gga_k_apbe_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_apbe_params));
  params = (gga_k_apbe_params *) (p->params);

  params->lambda = 0.0;

  switch(p->info->number){
  case XC_GGA_K_APBE:
    params->kappa = 0.8040;
    params->mu    = 0.23889;
    break;
  case XC_GGA_K_TW1:
    params->kappa = 0.8209;
    params->mu    = 0.2335;
    break;
  case XC_GGA_K_TW2:
    params->kappa = 0.6774;
    params->mu    = 0.2371;
    break;
  case XC_GGA_K_TW3:
    params->kappa = 0.8438;
    params->mu    = 0.2319;
    break;
  case XC_GGA_K_TW4:
    params->kappa = 0.8589;
    params->mu    = 0.2309;
    break;
  case XC_GGA_K_REVAPBE:
    params->kappa = 1.245;
    params->mu    = 0.23889;
    break;
  default:
    fprintf(stderr, "Internal error in gga_k_apbe\n");
    exit(1);
  }
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_apbe.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_apbe = {
  XC_GGA_K_APBE,
  XC_KINETIC,
  "mu fixed from the semiclassical neutral atom",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_186406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_apbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_revapbe = {
  XC_GGA_K_REVAPBE,
  XC_KINETIC,
  "revised APBE",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_186406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_apbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_tw1 = {
  XC_GGA_K_TW1,
  XC_KINETIC,
  "Tran and Wesolowski set 1 (Table II)",
  XC_FAMILY_GGA,
  {&xc_ref_Tran2002_441, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_apbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_tw2 = {
  XC_GGA_K_TW2,
  XC_KINETIC,
  "Tran and Wesolowski set 2 (Table II)",
  XC_FAMILY_GGA,
  {&xc_ref_Tran2002_441, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_apbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_tw3 = {
  XC_GGA_K_TW3,
  XC_KINETIC,
  "Tran and Wesolowski set 3 (Table II)",
  XC_FAMILY_GGA,
  {&xc_ref_Tran2002_441, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_apbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_tw4 = {
  XC_GGA_K_TW4,
  XC_KINETIC,
  "Tran and Wesolowski set 4 (Table II)",
  XC_FAMILY_GGA,
  {&xc_ref_Tran2002_441, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_apbe_init, NULL,
  NULL, work_gga, NULL
};
