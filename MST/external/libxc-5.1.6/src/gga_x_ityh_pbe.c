/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_X_ITYH_PBE 623 /* short-range PBE functional */

typedef struct{
  double kappa, mu;
  double lambda; /* parameter used in the Odashima & Capelle versions */
} gga_x_ityh_pbe_params;

#define N_PAR 3
static const char  *names[N_PAR]  = {"_kappa", "_mu", "_omega"};
static const char  *desc[N_PAR]   = {
  "Asymptotic value of the enhancement function",
  "Coefficient of the 2nd order expansion",
  "Range-separation parameter"};

static const double par_pbe[N_PAR] = {0.8040, MU_PBE, 0.33};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ityh_pbe.c"
#include "work_gga.c"

static void
gga_x_ityh_pbe_init(xc_func_type *p)
{
  gga_x_ityh_pbe_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ityh_pbe_params));
  params = (gga_x_ityh_pbe_params *) (p->params);

  /* This has to be explicitly initialized here */
  params->lambda = 0.0;

  xc_hyb_init_hybrid(p, 0.0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ityh_pbe = {
  XC_GGA_X_ITYH_PBE,
  XC_EXCHANGE,
  "Short-range recipe for PBE functional",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1996_3865, &xc_ref_Perdew1996_3865_err, &xc_ref_Iikura2001_3540, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_pbe, set_ext_params_cpy_omega},
  gga_x_ityh_pbe_init, NULL,
  NULL, work_gga, NULL
};
