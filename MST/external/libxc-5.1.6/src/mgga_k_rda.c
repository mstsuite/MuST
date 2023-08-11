/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_K_RDA          621 /* RDA by Karasiev et al */

typedef struct{
  double A0, A1, A2, A3;
  double beta1, beta2, beta3;
  double a, b, c;
} mgga_k_rda_params;

static void
mgga_k_rda_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_k_rda_params));
}

#define N_PAR 10
static const char  *names[N_PAR]  = {"_A0", "_A1", "_A2", "_A3", "_beta1", "_beta2", "_beta3", "_a", "_b", "_c"};
static const char  *desc[N_PAR]   = {"A0", "A1", "A2", "A3", "beta1", "beta2", "beta3", "a", "b", "c"};

static const double par_rda[N_PAR] = {0.50616,  3.04121, -0.34567, -1.89738,
  1.29691,  0.56184, 0.21944, 46.47662, 18.80658, -0.90346};
  
#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_k_rda.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_rda = {
  XC_MGGA_K_RDA,
  XC_KINETIC,
  "Reduced derivative approximation by Karasiev et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Karasiev2009_245120, NULL, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_rda, set_ext_params_cpy},
  mgga_k_rda_init, NULL,
  NULL, NULL, work_mgga
};
