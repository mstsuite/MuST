/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_LGAP_GE          633  /* LGAP_GE by Constantin et al */

typedef struct{
  double mu[3];
} gga_k_lgap_ge_params;

static void
gga_k_lgap_ge_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_lgap_ge_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_lgap_ge.c"
#include "work_gga.c"

#define N_PAR 3
static const char  *names[N_PAR]    = {"_mu1", "_mu2", "_mu3"};
static const char  *desc[N_PAR]     = {"mu1 parameter", "mu2 parameter", "mu3 parameter"};
static const double par_lgap_ge[N_PAR] = {0.0131, 0.18528, 0.0262};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_lgap_ge = {
  XC_GGA_K_LGAP_GE,
  XC_KINETIC,
  "LGAP-GE by Constantin et al",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2017_115153, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_lgap_ge, set_ext_params_cpy},
  gga_k_lgap_ge_init, NULL,
  NULL, work_gga, NULL
};
