/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_BR89_EXPLICIT    586 /* Becke-Roussel 89 with an explicit inversion of x(y), gamma = 0.8 */
#define XC_MGGA_X_BR89_EXPLICIT_1  602 /* Becke-Roussel 89 with an explicit inversion of x(y), gamma = 1.0 */

typedef struct{
  double gamma;
} mgga_x_br89_params;

#define BR89_N_PAR 1
static const char  *br89_names[BR89_N_PAR]    = {"_gamma"};
static const char  *br89_desc[BR89_N_PAR]     = {"gamma"};
static const double br89_values[BR89_N_PAR]   = {0.8};
static const double br89_1_values[BR89_N_PAR] = {1.0};

static void
mgga_x_br89_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_br89_params));
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_br89_explicit.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_br89_explicit = {
  XC_MGGA_X_BR89_EXPLICIT,
  XC_EXCHANGE,
  "Becke-Roussel 89 with an explicit inversion of x(y), gamma = 0.8",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1989_3761, &xc_ref_Proynov2008_103, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-12,
  {BR89_N_PAR, br89_names, br89_desc, br89_values, set_ext_params_cpy},
  mgga_x_br89_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_br89_explicit_1 = {
  XC_MGGA_X_BR89_EXPLICIT_1,
  XC_EXCHANGE,
  "Becke-Roussel 89 with an explicit inversion of x(y), gamma = 1.0",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1989_3761, &xc_ref_Proynov2008_103, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-12,
  {BR89_N_PAR, br89_names, br89_desc, br89_1_values, set_ext_params_cpy},
  mgga_x_br89_init, NULL,
  NULL, NULL, work_mgga
};
