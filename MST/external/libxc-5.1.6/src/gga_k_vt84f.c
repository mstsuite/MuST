/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_VT84F       619 /* VT84F by Karasiev et al */

typedef struct{
  double mu;
  double alpha;
} gga_k_vt84f_params;

static void
gga_k_vt84f_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_vt84f_params));
}

#define N_PAR 2
static const char  *names[N_PAR]  = {"_mu", "_alpha"};
static const char  *desc[N_PAR]   = { "mu parameter", "alpha parameter"};
static const double par_vt84f[N_PAR] = {2.778, 1.2965};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_vt84f.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_vt84f = {
  XC_GGA_K_VT84F,
  XC_KINETIC,
  "VT84F by Karasiev et al",
  XC_FAMILY_GGA,
  {&xc_ref_Karasiev2013_161108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_vt84f, set_ext_params_cpy},
  gga_k_vt84f_init, NULL,
  NULL, work_gga, NULL
};
