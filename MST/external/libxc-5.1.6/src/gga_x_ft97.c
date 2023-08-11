/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_FT97_A       114 /* Filatov & Thiel 97 (version A) */
#define XC_GGA_X_FT97_B       115 /* Filatov & Thiel 97 (version B) */

typedef struct{
  double beta0, beta1, beta2;
} gga_x_ft97_params;

static void
gga_x_ft97_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ft97_params));
}

#define FT97_N_PAR 3
static const char  *ft97_names[FT97_N_PAR]  = {"_beta0", "_beta1", "_beta2"};
static const char  *ft97_desc[FT97_N_PAR]   = {
  "beta0", "beta1", "beta2"
};
static const double ft97a_values[FT97_N_PAR] =
  {0.00293, 0.0, 0.0};
/* These parameters are what Filatov and Thiel actually used, not
   the ones they published in the paper... the differences being that
   beta1 has one more digit, and beta2 is squared: 2501.149^2 */
static const double ft97b_values[FT97_N_PAR] =
  {0.002913644, 0.0009474169, 6255746.320201};


#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ft97.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ft97_a = {
  XC_GGA_X_FT97_A,
  XC_EXCHANGE,
  "Filatov & Thiel 97 (version A)",
  XC_FAMILY_GGA,
  {&xc_ref_Filatov1997_847, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {FT97_N_PAR, ft97_names, ft97_desc, ft97a_values, set_ext_params_cpy},
  gga_x_ft97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ft97_b = {
  XC_GGA_X_FT97_B,
  XC_EXCHANGE,
  "Filatov & Thiel 97 (version B)",
  XC_FAMILY_GGA,
  {&xc_ref_Filatov1997_847, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {FT97_N_PAR, ft97_names, ft97_desc, ft97b_values, set_ext_params_cpy},
  gga_x_ft97_init, NULL,
  NULL, work_gga, NULL
};
