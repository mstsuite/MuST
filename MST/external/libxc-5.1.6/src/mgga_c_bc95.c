/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_BC95          240 /* Becke correlation 95 */

typedef struct{
  double css, copp;
} mgga_c_bc95_params;


static void
mgga_c_bc95_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_bc95_params));
}

#define BC95_N_PAR 2
static const char  *bc95_names[BC95_N_PAR]  = {"_css", "_copp"};
static const char  *bc95_desc[BC95_N_PAR]   = {
  "Parallel spin",
  "Opposite spin"
};
static const double bc95_values[BC95_N_PAR] = {0.038, 0.0031};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_bc95.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_bc95 = {
  XC_MGGA_C_BC95,
  XC_CORRELATION,
  "Becke correlation 95",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1996_1040, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {BC95_N_PAR, bc95_names, bc95_desc, bc95_values, set_ext_params_cpy},
  mgga_c_bc95_init, NULL,
  NULL, NULL, work_mgga,
};

