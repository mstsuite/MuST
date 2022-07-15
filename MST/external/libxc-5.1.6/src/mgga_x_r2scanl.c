/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_X_R2SCANL       718 /* Deorbitalized r^2SCAN exchange */

#define N_PAR 10
static const char *names[N_PAR] = {"_c1", "_c2", "_d", "_k1", "_taur", "_alphar", "_eta", "_dp2", "_a", "_b"};
static const char *desc[N_PAR] = {"c1 parameter", "c2 parameter", "d parameter",
  "k1 parameter", "taur parameter", "alphar parameter", "eta parameter", "dp2 parameter", "a parameter", "b parameter"};

static const double par_r2scanl[N_PAR] = {0.667, 0.8, 1.24, 0.065, 1.0e-4, 1.0e-3, 0.001, 0.361, 1.784720, 0.258304};

static void
r2scanl_set_ext_params(xc_func_type *p, const double *ext_params)
{
  const double *par_r2scan = NULL, *par_pc07 = NULL;
  if(ext_params != NULL) {
    par_r2scan = ext_params;
    par_pc07 = ext_params+8;
  }
  assert(p != NULL && p->func_aux != NULL);
  xc_func_set_ext_params(p->func_aux[0], par_r2scan);
  xc_func_set_ext_params(p->func_aux[1], par_pc07);
}

static void
mgga_x_r2scanl_init(xc_func_type *p)
{
  xc_deorbitalize_init(p, XC_MGGA_X_R2SCAN, XC_MGGA_K_PC07_OPT);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_r2scanl = {
  XC_MGGA_X_R2SCANL,
  XC_EXCHANGE,
  "Deorbitalized re-regularized SCAN (r2SCAN-L) exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2020_121109, &xc_ref_Furness2020_8208, &xc_ref_Furness2020_9248, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, par_r2scanl, r2scanl_set_ext_params},
  mgga_x_r2scanl_init, NULL,
  NULL, NULL, xc_deorbitalize_func
};

