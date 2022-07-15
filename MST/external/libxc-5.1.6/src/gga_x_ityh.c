/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_ITYH 529 /* short-range recipe B88 functionals - erf */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ityh.c"
#include "work_gga.c"

static void
xc_gga_x_ityh_init(xc_func_type *p)
{
  xc_hyb_init_hybrid(p, 0.0);
}

static const char  *omega_names[]  = {"_omega"};
static const char  *omega_desc[]   = {"screening parameter"};
static const double omega_values[] = {0.2};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ityh = {
  XC_GGA_X_ITYH,
  XC_EXCHANGE,
  "Short-range recipe for B88 functional - erf",
  XC_FAMILY_GGA,
  {&xc_ref_Iikura2001_3540, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14, /* functional does not seem to be stable below this value */
  {1, omega_names, omega_desc, omega_values, set_ext_params_cpy_omega},
  xc_gga_x_ityh_init, NULL,
  NULL, work_gga, NULL
};
