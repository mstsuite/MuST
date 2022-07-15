/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_PG1      219 /* PG1 functional by Constantin, Fabiano, and Della Sala */

typedef struct{
  double pg_mu;
} gga_k_pg_params;

static void
gga_k_pg_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_pg_params));
}

#define N_PAR 1
static const char  *names[N_PAR]  = {"_mu"};
static const char  *desc[N_PAR]   = {"Prefactor in exponent"};
static const double par_pg1[N_PAR] = {1.0};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_pg.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_pg1 = {
  XC_GGA_K_PG1,
  XC_KINETIC,
  "PG1 (Pauli-Gaussian) functional by Constantin, Fabiano, and Della Sala",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2018_4385, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_pg1, set_ext_params_cpy},
  gga_k_pg_init, NULL,
  NULL, work_gga, NULL
};
