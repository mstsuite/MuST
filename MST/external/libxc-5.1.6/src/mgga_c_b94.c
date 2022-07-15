/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_B94           397 /* Meta-GGA correlation by Becke */
#define XC_HYB_MGGA_XC_B94_HYB  398 /* Hybrid meta-GGA by Becke      */

typedef struct{
  double gamma; /* gamma parameter in mgga_x_br89 */
  double css; /* same-spin constant */
  double cab; /* opposite-spin constant */
} mgga_c_b94_params;

#define N_PAR 3
static const char  *names[N_PAR]      = {"_gamma", "_css", "_cab"};
static const char  *desc[N_PAR]       = {"gamma", "css", "cab"};

static const double par_b94[N_PAR]    = {1.0, 0.88, 0.63};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_b94.c"
#include "work_mgga.c"

static void
mgga_c_b94_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_b94_params));
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_b94 = {
  XC_MGGA_C_B94,
  XC_CORRELATION,
  "Becke 1994 meta-GGA correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1994_625, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, par_b94, set_ext_params_cpy},
  mgga_c_b94_init, NULL,
  NULL, NULL, work_mgga
};

#define N_PAR_HYB 4
static const char  *names_hyb[N_PAR_HYB]      = {"_gamma", "_css", "_cab", "_cx"};
static const char  *desc_hyb[N_PAR_HYB]       = {"gamma", "css", "cab", "cx"};
static const double par_b94_hyb[N_PAR_HYB]    = {1.0, 0.88, 0.66, 0.154};

static void
hyb_mgga_xc_b94_hyb_init(xc_func_type *p)
{
  static int   funcs_id   [2] = {XC_MGGA_X_BR89, XC_MGGA_C_B94};
  static double funcs_coef[2] = {0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0);
}

static void
hyb_mgga_xc_b94_hyb_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double gamma, css, cab, cx;

  assert(p != NULL);

  gamma = get_ext_param(p, ext_params, 0);
  css = get_ext_param(p, ext_params, 1);
  cab = get_ext_param(p, ext_params, 2);
  cx = get_ext_param(p, ext_params, 3);

  xc_func_set_ext_params_name(p->func_aux[0], "_at", 0.0);
  xc_func_set_ext_params_name(p->func_aux[0], "_gamma", gamma);
  xc_func_set_ext_params_name(p->func_aux[1], "_gamma", gamma);
  xc_func_set_ext_params_name(p->func_aux[1], "_css", css);
  xc_func_set_ext_params_name(p->func_aux[1], "_cab", cab);

  p->mix_coef[0] = 1.0 - cx;
  p->mix_coef[1] = 1.0;

  p->cam_alpha = cx;
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_b94_hyb = {
  XC_HYB_MGGA_XC_B94_HYB,
  XC_EXCHANGE_CORRELATION,
  "Becke 1994 hybrid meta-GGA",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Becke1994_625, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR_HYB, names_hyb, desc_hyb, par_b94_hyb, hyb_mgga_xc_b94_hyb_set_ext_params},
  hyb_mgga_xc_b94_hyb_init, NULL,
  NULL, NULL, NULL
};
