/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_JK          256 /* Jemmer-Knowles meta-GGA exchange */

typedef struct {
  double beta;
  double gamma;
} mgga_x_jk_params;

static void
mgga_x_jk_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_jk_params));
}

#define N_PAR 2
static const char  *names[N_PAR]  = {"_beta", "_gamma"};
static const char  *desc[N_PAR]   = {
  "beta/X_FACTOR_C is the coefficient of the gradient expansion",
  "gamma should be 6 to get the right asymptotics of Ex"};
static const double jk_values[N_PAR] =
  {0.0586, 6.0};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_jk.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_jk = {
  XC_MGGA_X_JK,
  XC_EXCHANGE,
  "Jemmer-Knowles meta-GGA exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Jemmer1995_3571, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, jk_values, set_ext_params_cpy},
  mgga_x_jk_init, NULL,
  NULL, NULL, work_mgga,
};
