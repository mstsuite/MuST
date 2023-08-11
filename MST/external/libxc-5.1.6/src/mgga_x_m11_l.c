/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_M11_L        226 /* M11-L exchange functional from Minnesota  */

typedef struct{
  const double a[12], b[12], c[12], d[12];
} mgga_x_m11_l_params;

#define N_PAR 49
static const char *names[N_PAR] = {
  "_a0", "_a1", "_a2", "_a3", "_a4", "_a5",
  "_a6", "_a7", "_a8", "_a9", "_a10", "_a11",
  "_b0", "_b1", "_b2", "_b3", "_b4", "_b5",
  "_b6", "_b7", "_b8", "_b9", "_b10", "_b11",
  "_c0", "_c1", "_c2", "_c3", "_c4", "_c5",
  "_c6", "_c7", "_c8", "_c9", "_c10", "_c11",
  "_d0", "_d1", "_d2", "_d3", "_d4", "_d5",
  "_d6", "_d7", "_d8", "_d9", "_d10", "_d11",
  "_omega"
};
static const char *desc[N_PAR] = {
  "a0 parameter", "a1 parameter", "a2 parameter", "a3 parameter",
  "a4 parameter", "a5 parameter", "a6 parameter", "a7 parameter",
  "a8 parameter", "a9 parameter", "a10 parameter", "a11 parameter",
  "b0 parameter", "b1 parameter", "b2 parameter", "b3 parameter",
  "b4 parameter", "b5 parameter", "b6 parameter", "b7 parameter",
  "b8 parameter", "b9 parameter", "b10 parameter", "b11 parameter",
  "c0 parameter", "c1 parameter", "c2 parameter", "c3 parameter",
  "c4 parameter", "c5 parameter", "c6 parameter", "c7 parameter",
  "c8 parameter", "c9 parameter", "c10 parameter", "c11 parameter",
  "d0 parameter", "d1 parameter", "d2 parameter", "d3 parameter",
  "d4 parameter", "d5 parameter", "d6 parameter", "d7 parameter",
  "d8 parameter", "d9 parameter", "d10 parameter", "d11 parameter",
  "range separation"
};

static const double par_m11_l[N_PAR] = {
  8.121131e-01,  1.738124e+01,  1.154007e+00,  6.869556e+01,  1.016864e+02, -5.887467e+00,
  4.517409e+01, -2.773149e+00, -2.617211e+01,  0.000000e+00,  0.000000e+00,  0.000000e+00,
  1.878869e-01, -1.653877e+01,  6.755753e-01, -7.567572e+01, -1.040272e+02,  1.831853e+01,
  -5.573352e+01, -3.520210e+00,  3.724276e+01,  0.000000e+00,  0.000000e+00,  0.000000e+00,
  -4.386615e-01, -1.214016e+02, -1.393573e+02, -2.046649e+00,  2.804098e+01, -1.312258e+01,
  -6.361819e+00, -8.055758e-01,  3.736551e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,
  1.438662e+00,  1.209465e+02,  1.328252e+02,  1.296355e+01,  5.854866e+00, -3.378162e+00,
  -4.423393e+01,  6.844475e+00,  1.949541e+01,  0.000000e+00,  0.000000e+00,  0.000000e+00,
  0.25
};

static void
mgga_x_m11_l_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_m11_l_params));

  xc_hyb_init_hybrid(p, 0.0);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_m11_l.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_m11_l = {
  XC_MGGA_X_M11_L,
  XC_EXCHANGE,
  "Minnesota M11-L exchange functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Peverati2012_117, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-13,
  {N_PAR, names, desc, par_m11_l, set_ext_params_cpy_omega},
  mgga_x_m11_l_init, NULL,
  NULL, NULL, work_mgga,
};
