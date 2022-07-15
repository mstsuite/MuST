/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_K_PGSL025  220 /* PGSL0.25 functional by Constantin, Fabiano, and Della Sala */

typedef struct{
  double pgslb_mu, pgslb_beta;
} mgga_k_pgslb_params;

static void
mgga_k_pgslb_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_k_pgslb_params));
}

#define N_PAR 2
static const char  *names[N_PAR]  = {"_mu", "_beta"};
static const char  *desc[N_PAR]   = {"Prefactor in exponent", "Coefficient of Laplacian term"};
static const double par_pgsl025[N_PAR] = {40.0/27.0, 0.25};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_k_pgslb.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_pgsl025 = {
  XC_MGGA_K_PGSL025,
  XC_KINETIC,
  "PGSL025 (Pauli-Gaussian) functional by Constantin, Fabiano, and Della Sala",
  XC_FAMILY_MGGA,
  {&xc_ref_Constantin2018_4385, NULL, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_pgsl025, set_ext_params_cpy},
  mgga_k_pgslb_init, NULL,
  NULL, NULL, work_mgga
};
