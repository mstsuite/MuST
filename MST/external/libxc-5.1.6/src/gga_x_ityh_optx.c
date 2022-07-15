/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_X_ITYH_OPTX 622 /* short-range OPTX functional */

typedef struct{
  double a, b, gamma;
} gga_x_ityh_optx_params;

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ityh_optx.c"
#include "work_gga.c"

static void
xc_gga_x_ityh_optx_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ityh_optx_params));

  xc_hyb_init_hybrid(p, 0.0);
}

#define N_PAR 4
static const char  *names[N_PAR]  = {"_a", "_b", "_gamma", "_omega"};
static const char  *desc[N_PAR]   = {
  "a",
  "b",
  "gamma",
  "omega"};
static const double values[N_PAR] =
  {1.05151, 1.43169/X_FACTOR_C, 0.006, 0.2};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ityh_optx = {
  XC_GGA_X_ITYH_OPTX,
  XC_EXCHANGE,
  "Short-range recipe for OPTX functional",
  XC_FAMILY_GGA,
  {&xc_ref_Handy2001_403, &xc_ref_Iikura2001_3540, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, values, set_ext_params_cpy_omega},
  xc_gga_x_ityh_optx_init, NULL,
  NULL, work_gga, NULL
};
