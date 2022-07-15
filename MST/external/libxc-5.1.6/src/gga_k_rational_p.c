/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_RATIONAL_P      218 /* Lehtomaki and Lopez-Acevedo */

typedef struct{
  double C2; /* prefactor for s^2 term */
  double p; /* exponent */
} gga_k_rational_p_params;

static void
gga_k_rational_p_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_rational_p_params));
}

#define N_PAR 2
static const char  *names[N_PAR]  = {"_C2", "_p"};
static const char  *desc[N_PAR]   = {
  "Coefficient for s^2",
  "Exponent"};

static const double par_p32[N_PAR] =
  {0.7687, 3.0/2.0};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_rational_p.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_rational_p = {
  XC_GGA_K_RATIONAL_P,
  XC_KINETIC,
  "RATIONAL^p by Lehtomaki and Lopez-Acevedo (by default p=3/2, C2=0.7687)",
  XC_FAMILY_GGA,
  {&xc_ref_Lehtomaki2019_165111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_p32, set_ext_params_cpy},
  gga_k_rational_p_init, NULL,
  NULL, work_gga, NULL
};
