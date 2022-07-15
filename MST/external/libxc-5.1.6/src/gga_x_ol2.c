/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_OL2          183 /* Exchange form based on Ou-Yang and Levy v.2 */

typedef struct{
  double aa, bb, cc;
} gga_x_ol2_params;

static void
gga_x_ol2_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ol2_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ol2.c"
#include "work_gga.c"

#define OL2_N_PAR 3
static const char  *ol2_names[OL2_N_PAR]  = {"_aa", "_bb", "_cc"};
static const char  *ol2_desc[OL2_N_PAR]   = {
  "aa", "bb", "cc"
};
static const double ol2_values[OL2_N_PAR] = {
  0.09564574034649151285038696952714226444963L, /* M_CBRT2*0.07064/X_FACTOR_C */
  0.09564574034649151285038696952714226444963L, /* M_CBRT2*0.07064/X_FACTOR_C */
  4.098833606342553442039881031486386917472L    /* M_CBRT2*M_CBRT2*0.07064*34.0135/X_FACTOR_C */
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ol2 = {
  XC_GGA_X_OL2,
  XC_EXCHANGE,
  "Exchange form based on Ou-Yang and Levy v.2",
  XC_FAMILY_GGA,
  {&xc_ref_Fuentealba1995_31, &xc_ref_OuYang1991_379, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {OL2_N_PAR, ol2_names, ol2_desc, ol2_values, set_ext_params_cpy},
  gga_x_ol2_init, NULL,
  NULL, work_gga, NULL
};
