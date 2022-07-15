/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MS2B          300  /* MS2beta exchange by Furness and Sun */
#define XC_MGGA_X_MS2BS         301  /* MS2beta* exchange by Furness and Sun */

typedef struct{
  double kappa, c, b;
} mgga_x_msb_params;

static void
mgga_x_msb_init(xc_func_type *p)
{
  mgga_x_msb_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_msb_params));
  params = (mgga_x_msb_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_MS2B:
    params->kappa = 0.504;
    params->b     = (27.0*4.0 - 9.0)/64.0;
    params->c     = 0.14607;
    break;
  case XC_MGGA_X_MS2BS:
    params->kappa = 0.6263;
    params->b     = 4.3011;
    params->c     = 0.12268;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_msb\n");
    exit(1);
  }
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_msb.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms2b = {
  XC_MGGA_X_MS2B,
  XC_EXCHANGE,
  "MS2beta exchange of Furness and Sun",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2019_041119, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_msb_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms2bs = {
  XC_MGGA_X_MS2BS,
  XC_EXCHANGE,
  "MS2beta* exchange of Furness and Sun",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2018, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_msb_init, NULL,
  NULL, NULL, work_mgga,
};
