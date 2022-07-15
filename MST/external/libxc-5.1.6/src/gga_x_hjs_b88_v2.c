/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_HJS_B88_V2   46 /* HJS screened exchange corrected B88 version */

typedef struct{
  double a[6], b[9]; /* pointers to the a and b parameters */
} gga_x_hjs_params;

#define N_PARS 16
static const char  *names[N_PARS]  = {"_a0", "_a1", "_a2", "_a3", "_a4", "_a5", "_b0", "_b1", "_b2", "_b3", "_b4", "_b5", "_b6", "_b7", "_b8", "_omega"};
static const char  *desc[N_PARS]   = {"a0", "a1", "a2", "a3", "a4", "a5", "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "omega"};

static const double pars_B88_V2[N_PARS] =
  {0.0253933, -0.0673075, 0.0891476, -0.0454168, -0.00765813, 0.0142506,
   -2.6506, 3.91108, -3.31509, 1.54485, -0.198386, -0.136112, 0.0647862, 0.0159586, -0.000245066, 0.11};

static void
gga_x_hjs_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_hjs_params));

  xc_hyb_init_hybrid(p, 0.0);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_hjs_b88_v2.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_b88_v2 = {
  XC_GGA_X_HJS_B88_V2,
  XC_EXCHANGE,
  "HJS screened exchange B88 corrected version",
  XC_FAMILY_GGA,
  {&xc_ref_Weintraub2009_754, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-6, /* densities smaller than 1e-6 yield NaNs */
  {N_PARS, names, desc, pars_B88_V2, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, work_gga, NULL
};
