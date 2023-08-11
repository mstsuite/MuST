/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_SOGGA11        151 /* Second-order generalized gradient approximation 2011 */
#define XC_HYB_GGA_X_SOGGA11_X  426 /* Hybrid based on SOGGA11 form */

typedef struct{
  double kappa, mu, a[6], b[6];
} gga_x_sogga11_params;

#define N_PAR_PURE 14
static const char  *pure_names[N_PAR_PURE]  = {
  "_kappa", "_mu", "_a0", "_a1", "_a2", "_a3", "_a4",
  "_a5", "_b0", "_b1", "_b2", "_b3", "_b4", "_b5"
};
static const char  *pure_desc[N_PAR_PURE]   = {
  "kappa", "mu", "a0", "a1", "a2", "a3", "a4",
  "a5", "b0", "b1", "b2", "b3", "b4", "b5"
};

#define N_PAR_HYB 15
static const char  *hyb_names[N_PAR_HYB]  = {
  "_kappa", "_mu", "_a0", "_a1", "_a2", "_a3", "_a4",
  "_a5", "_b0", "_b1", "_b2", "_b3", "_b4", "_b5", "_cx"
};
static const char  *hyb_desc[N_PAR_HYB]   = {
  "kappa", "mu", "a0", "a1", "a2", "a3", "a4",
  "a5", "b0", "b1", "b2", "b3", "b4", "b5",
  "Fraction of exact exchange"
};

static const double par_sogga11[N_PAR_PURE] = {
  0.552, MU_GE,
  0.50000, -2.95535,  15.7974, -91.1804,  96.2030,  0.18683,
  0.50000,  3.50743, -12.9523,  49.7870, -33.2545, -11.1396
};

/* These coefficients include the factor (1-X) in the functional definition. */
static const double par_sogga11_x[N_PAR_HYB] = {
  0.552, MU_GE,
  0.29925,  3.21638, -3.55605,  7.65852, -11.2830, 5.25813,
  0.29925, -2.88595,  3.23617, -2.45393, -3.75495,  3.96613,
  0.4015
};

static void
gga_x_sogga11_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_sogga11_params));

  if(p->info->number == XC_HYB_GGA_X_SOGGA11_X)
    xc_hyb_init_hybrid(p, 0.0);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_sogga11.c"
#include "work_gga.c"


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_sogga11 = {
  XC_GGA_X_SOGGA11,
  XC_EXCHANGE,
  "Second-order generalized gradient approximation 2011",
  XC_FAMILY_GGA,
  {&xc_ref_Peverati2011_1991, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_PURE, pure_names, pure_desc, par_sogga11, set_ext_params_cpy},
  gga_x_sogga11_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_sogga11_x = {
  XC_HYB_GGA_X_SOGGA11_X,
  XC_EXCHANGE,
  "Hybrid based on SOGGA11 form",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Peverati2011_191102, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_HYB, hyb_names, hyb_desc, par_sogga11_x, set_ext_params_cpy_exx},
  gga_x_sogga11_init, NULL,
  NULL, work_gga, NULL
};
