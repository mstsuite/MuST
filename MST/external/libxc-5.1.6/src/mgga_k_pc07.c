/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_K_PC07          543 /* Perdew and Constantin 2007 */
#define XC_MGGA_K_PC07_OPT      634 /* Reoptimized version by Mejia-Rodriguez and Trickey */

typedef struct{
  double a, b;
} mgga_k_pc07_params;

static void
mgga_k_pc07_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_k_pc07_params));
}

#define PC07_N_PAR 2
static const char  *pc07_names[PC07_N_PAR]  = {"_a", "_b"};
static const char  *pc07_desc[PC07_N_PAR]   = { "a",  "b"};
static const double pc07_values[PC07_N_PAR] = {0.5389, 3};
static const double pc07opt_values[PC07_N_PAR] = {1.784720, 0.258304};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_k_pc07.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_pc07 = {
  XC_MGGA_K_PC07,
  XC_KINETIC,
  "Perdew and Constantin 2007",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2007_155109, NULL, NULL, NULL, NULL},
  XC_FLAGS_DEVELOPMENT | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PC07_N_PAR, pc07_names, pc07_desc, pc07_values, set_ext_params_cpy},
  mgga_k_pc07_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_pc07_opt = {
  XC_MGGA_K_PC07_OPT,
  XC_KINETIC,
  "Reoptimized PC07 by Mejia-Rodriguez and Trickey",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Perdew2007_155109, NULL, NULL, NULL},
  XC_FLAGS_DEVELOPMENT | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PC07_N_PAR, pc07_names, pc07_desc, pc07opt_values, set_ext_params_cpy},
  mgga_k_pc07_init, NULL,
  NULL, NULL, work_mgga,
};
