/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_LGAP          620  /* LGAP by Constantin et al */

typedef struct{
  double kappa;
  double mu[3];
} gga_k_lgap_params;


static void
gga_k_lgap_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_lgap_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_lgap.c"
#include "work_gga.c"

#define N_PAR 4
static const char  *names[N_PAR]    = {"_kappa", "_mu1", "_mu2", "_mu3"};
static const char  *desc[N_PAR]     = {"kappa parameter", "mu1 parameter", "mu2 parameter", "mu3 parameter"};

/* Annoyingly the parameters are not given in closed form; I solved these with Maple with 20 digits. */
static const double par_lgap[N_PAR] = {0.8, 0.016375000000000000000, 0.23173407031250000000, 0.036536516531615524505};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_lgap = {
  XC_GGA_K_LGAP,
  XC_KINETIC,
  "LGAP by Constantin et al",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2017_115153, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_lgap, set_ext_params_cpy},
  gga_k_lgap_init, NULL,
  NULL, work_gga, NULL
};
