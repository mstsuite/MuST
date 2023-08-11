/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_K_L04      617 /* L0.4 by Laricchia et al */
#define XC_MGGA_K_L06      618 /* L0.6 by Laricchia et al */

typedef struct{
  double kappa;
} mgga_k_lk_params;

static void
mgga_k_lk_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_k_lk_params));
}

#define N_PAR 1
static const char  *names[N_PAR]  = {"_kappa"};
static const char  *desc[N_PAR]   = {"kappa parameter"};
static const double par_l04[N_PAR] = {0.402};
static const double par_l06[N_PAR] = {0.623};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_k_lk.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_l04 = {
  XC_MGGA_K_L04,
  XC_KINETIC,
  "L0.4 by Laricchia et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Laricchia2014_164, NULL, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_l04, set_ext_params_cpy},
  mgga_k_lk_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_l06 = {
  XC_MGGA_K_L06,
  XC_KINETIC,
  "L0.6 by Laricchia et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Laricchia2014_164, NULL, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_l06, set_ext_params_cpy},
  mgga_k_lk_init, NULL,
  NULL, NULL, work_mgga
};
