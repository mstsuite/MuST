/*
 Copyright (C) 2006-2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* Local tau approximation */

#define XC_MGGA_X_RLDA          688 /* Reparametrized local-density approximation */
#define XC_MGGA_X_MK00          230 /* Exchange for accurate virtual orbital energies */
#define XC_MGGA_X_MK00B         243 /* Exchange for accurate virtual orbital energies (v. B) */

typedef struct{
  double prefactor;
} mgga_x_rlda_params;

static void
mgga_x_rlda_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_rlda_params));
}

#define RLDA_N_PAR 1
static const char  *rlda_names[RLDA_N_PAR]  = {"_prefactor"};
static const char  *rlda_desc[RLDA_N_PAR]   = {
  "Prefactor that multiplies functional"
};
static const double rlda_values[RLDA_N_PAR]  = {1.0};
static const double mk00_values[RLDA_N_PAR]  = {4.0/5.0};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_rlda.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_rlda = {
  XC_MGGA_X_RLDA,
  XC_EXCHANGE,
  "Reparametrized local-density approximation",
  XC_FAMILY_MGGA,
  {&xc_ref_Campi1978_263, &xc_ref_Koehl1996_835, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS | XC_FLAGS_DEVELOPMENT,
  1e-15,
  {RLDA_N_PAR, rlda_names, rlda_desc, rlda_values, set_ext_params_cpy},
  mgga_x_rlda_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mk00 = {
  XC_MGGA_X_MK00,
  XC_EXCHANGE,
  "Exchange for accurate virtual orbital energies",
  XC_FAMILY_MGGA,
  {&xc_ref_Manby2000_7002, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS | XC_FLAGS_DEVELOPMENT,
  1e-15,
  {RLDA_N_PAR, rlda_names, rlda_desc, mk00_values, set_ext_params_cpy},
  mgga_x_rlda_init, NULL,
  NULL, NULL, work_mgga,
};


static void
mgga_x_mk00b_init(xc_func_type *p)
{
  static int    funcs_id  [3] = {XC_LDA_X, XC_GGA_X_B88, XC_MGGA_X_MK00};
  static double funcs_coef[3] = {-1.0, 1.0, 1.0};

  static double par_x_b88[] = {0.0016, 6.0};

  xc_mix_init(p, 3, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[1], par_x_b88);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mk00b = {
  XC_MGGA_X_MK00B,
  XC_EXCHANGE,
  "Exchange for accurate virtual orbital energies (v. B)",
  XC_FAMILY_MGGA,
  {&xc_ref_Manby2000_7002, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS | XC_FLAGS_DEVELOPMENT,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_mk00b_init, NULL,
  NULL, NULL, NULL,
};

