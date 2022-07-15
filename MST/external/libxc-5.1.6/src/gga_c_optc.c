/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_OPTC       200 /* Optimized correlation functional of Cohen and Handy */

typedef struct{
  double c1, c2;
} gga_c_optc_params;

static void
gga_c_optc_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_optc_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_optc.c"
#include "work_gga.c"

#define OPTC_N_PAR 2
static const char  *optc_names[OPTC_N_PAR]  = {"_c1", "_c2"};
static const char  *optc_desc[OPTC_N_PAR]   = {"c1", "c2"};
static const double optc_values[OPTC_N_PAR] =
  {1.1015L, 0.6625L};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_optc = {
  XC_GGA_C_OPTC,
  XC_CORRELATION,
  "Optimized correlation functional of Cohen and Handy",
  XC_FAMILY_GGA,
  {&xc_ref_Cohen2001_607, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {OPTC_N_PAR, optc_names, optc_desc, optc_values, set_ext_params_cpy},
  gga_c_optc_init, NULL,
  NULL, work_gga, NULL
};
