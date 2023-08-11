/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_DK87_R1 111 /* dePristo & Kress 87 (version R1) */
#define XC_GGA_X_DK87_R2 112 /* dePristo & Kress 87 (version R2) */

typedef struct {
  double a1, b1, alpha;
} gga_x_dk87_params;

#define N_PAR 3
static const char *names[N_PAR] = {"_a1", "_b1", "_alpha"};
static const char *desc[N_PAR] = {"a1 parameter", "b1 parameter",
                                  "alpha parameter"};

static const double par_dk87_r1[N_PAR] = {0.861504, 0.044286, 1.0};

static const double par_dk87_r2[N_PAR] = {0.861213, 0.042076, 0.98};

static void gga_x_dk87_init(xc_func_type *p) {
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_dk87_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_dk87.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_dk87_r1 = {
  XC_GGA_X_DK87_R1,
  XC_EXCHANGE,
  "dePristo & Kress 87 version R1",
  XC_FAMILY_GGA,
  {&xc_ref_DePristo1987_1425, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_dk87_r1, set_ext_params_cpy},
  gga_x_dk87_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_dk87_r2 = {
  XC_GGA_X_DK87_R2,
  XC_EXCHANGE,
  "dePristo & Kress 87 version R2",
  XC_FAMILY_GGA,
  {&xc_ref_DePristo1987_1425, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_dk87_r2, set_ext_params_cpy},
  gga_x_dk87_init, NULL,
  NULL, work_gga, NULL
};
