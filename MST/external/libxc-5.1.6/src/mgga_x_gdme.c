/*
 Copyright (C) 2006-2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* Local tau approximation */

#define XC_MGGA_X_GDME_NV   687 /* Generalized density-matrix with a=1/2      */
#define XC_MGGA_X_GDME_0    689 /* Generalized density-matrix with a=0        */
#define XC_MGGA_X_GDME_KOS  690 /* Generalized density-matrix with a=0.00638  */
#define XC_MGGA_X_GDME_VT   691 /*  Varied-terms (VT) mGGA of Koehl, Odom, and Scuseria */

typedef struct{
  double a, AA, BB;
} mgga_x_gdme_params;


static void
mgga_x_gdme_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_gdme_params));
}

#define GDME_N_PAR 3
static const char  *gdme_names[GDME_N_PAR]  = {"_a", "_AA", "_BB"};
static const char  *gdme_desc[GDME_N_PAR]   = {
  "center of the s expansion of density-matrix",
  "parameter of the first (LDA) term",
  "parameter of the correction term"
};
static const double gdme_nv_values[GDME_N_PAR]  = {0.5, 9.0*M_PI/4.0, 35.0*M_PI/12.0};
static const double gdme_0_values[GDME_N_PAR]   = {0.0, 9.0*M_PI/4.0, 35.0*M_PI/12.0};
static const double gdme_kos_values[GDME_N_PAR] = {0.00638, 9.0*M_PI/4.0, 35.0*M_PI/12.0};
static const double gdme_vt_values[GDME_N_PAR]  = {0.0, 7.31275, 5.43182};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_gdme.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gdme_nv = {
  XC_MGGA_X_GDME_NV,
  XC_EXCHANGE,
  "Generalized density-matrix with a=1/2",
  XC_FAMILY_MGGA,
  {&xc_ref_Negele1972_1472, &xc_ref_Koehl1996_835, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-15,
  {GDME_N_PAR, gdme_names, gdme_desc, gdme_nv_values, set_ext_params_cpy},
  mgga_x_gdme_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gdme_0 = {
  XC_MGGA_X_GDME_0,
  XC_EXCHANGE,
  "Generalized density-matrix with a=0",
  XC_FAMILY_MGGA,
  {&xc_ref_Koehl1996_835, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-15,
  {GDME_N_PAR, gdme_names, gdme_desc, gdme_0_values, set_ext_params_cpy},
  mgga_x_gdme_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gdme_kos = {
  XC_MGGA_X_GDME_KOS,
  XC_EXCHANGE,
  "Generalized density-matrix with a=0.00638",
  XC_FAMILY_MGGA,
  {&xc_ref_Koehl1996_835, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-15,
  {GDME_N_PAR, gdme_names, gdme_desc, gdme_kos_values, set_ext_params_cpy},
  mgga_x_gdme_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gdme_vt = {
  XC_MGGA_X_GDME_VT,
  XC_EXCHANGE,
  "Varied-terms (VT) mGGA of Koehl, Odom, and Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Koehl1996_835, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-15,
  {GDME_N_PAR, gdme_names, gdme_desc, gdme_vt_values, set_ext_params_cpy},
  mgga_x_gdme_init, NULL,
  NULL, NULL, work_mgga,
};

