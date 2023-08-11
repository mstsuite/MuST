/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_WI0 153 /* Wilson & Ivanov initial version */
#define XC_GGA_C_WI  148 /* Wilson & Ivanov */

typedef struct {
  double a, b, c, d, k;
} gga_c_wi_params;

#define N_PAR 5
static const char *names[N_PAR] = {"_a", "_b", "_c", "_d", "_k"};
static const char *desc[N_PAR] = {"a parameter", "b parameter", "c parameter",
                                  "d parameter", "k parameter"};

static const double wi0_par[N_PAR] = {-0.44, 0.0032407, 7.8, 0.0073,
                                         0.000311};
static const double wi_par[N_PAR] = {-0.00652, 0.0007, 0.21, 0.002, 0.001};

static void gga_c_wi_init(xc_func_type *p) {
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_wi_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_wi.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_wi0 = {
  XC_GGA_C_WI0,
  XC_CORRELATION,
  "Wilson & Ivanov initial version",
  XC_FAMILY_GGA,
  {&xc_ref_Wilson1998_523, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, wi0_par, set_ext_params_cpy},
  gga_c_wi_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_wi = {
  XC_GGA_C_WI,
  XC_CORRELATION,
  "Wilson & Ivanov",
  XC_FAMILY_GGA,
  {&xc_ref_Wilson1998_523, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-10,
  {N_PAR, names, desc, wi_par, set_ext_params_cpy},
  gga_c_wi_init, NULL,
  NULL, work_gga, NULL
};
