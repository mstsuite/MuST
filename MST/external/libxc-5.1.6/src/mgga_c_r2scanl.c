/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_R2SCANL        719 /* Deorbitalized r^2SCAN correlation */

#define N_PAR 2
static const char *names[N_PAR] = {
  "_a", "_b"  /* parameters of pc07 */
};
static const char *desc[N_PAR] = {
  "pc07 a", "pc07 b"
};
static const double par_r2scanl[N_PAR] = {
  1.784720, 0.258304
};

static void
r2scanl_set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL && p->func_aux != NULL);
  xc_func_set_ext_params(p->func_aux[1], &ext_params[0]);
}

static void
mgga_c_r2scanl_init(xc_func_type *p)
{
  xc_deorbitalize_init(p, XC_MGGA_C_R2SCAN, XC_MGGA_K_PC07_OPT);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_r2scanl = {
  XC_MGGA_C_R2SCANL,
  XC_CORRELATION,
  "Deorbitalized re-regularized SCAN (r2SCAN-L) correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2020_121109, &xc_ref_Furness2020_8208, &xc_ref_Furness2020_9248, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, par_r2scanl, r2scanl_set_ext_params},
  mgga_c_r2scanl_init, NULL,
  NULL, NULL, xc_deorbitalize_func
};
