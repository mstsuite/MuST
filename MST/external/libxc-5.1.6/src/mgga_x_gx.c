/*
 Copyright (C) 2017 Miguel Marques, Mario Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_GX          575 /* GX functional of Loos */

typedef struct{
  double c0, c1, alphainf;
} mgga_x_gx_params;

static void
mgga_x_gx_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_gx_params));
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_gx.c"
#include "work_mgga.c"

#define GX_N_PAR 3
static const char  *gx_names[GX_N_PAR]  = {"_c0", "_c1", "_alphainf"};
static const char  *gx_desc[GX_N_PAR]   = {"c0", "c1", "alphainf"};
static const double gx_values[GX_N_PAR] = {0.827411L, -0.643560L, 0.852};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gx = {
  XC_MGGA_X_GX,
  XC_EXCHANGE,
  "GX functional of Loos",
  XC_FAMILY_MGGA,
  {&xc_ref_Loos2017_114108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {GX_N_PAR, gx_names, gx_desc, gx_values, set_ext_params_cpy},
  mgga_x_gx_init, NULL,
  NULL, NULL, work_mgga,
};
