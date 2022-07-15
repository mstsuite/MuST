/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_CCDF         313 /* ccDF, coupled-cluster based density functional */

typedef struct{
  double c1;
  double c2;
  double c3;
  double c4;
  double c5;
} gga_c_ccdf_params;

static void gga_c_ccdf_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_ccdf_params));
}

#define N_PAR 5
static const char  *names[N_PAR]  = {"_c1", "_c2", "_c3", "_c4", "_c5"};
static const char  *desc[N_PAR]   = {
  "c1 parameter",
  "c2 parameter",
  "c3 parameter",
  "c4 parameter",
  "c5 parameter"};

static const double par_ccdf[N_PAR] = {-0.0468, 0.023, 0.544, 23.401, 0.479};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_ccdf.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_ccdf = {
  XC_GGA_C_CCDF,
  XC_CORRELATION,
  "ccDF: coupled-cluster motivated density functional",
  XC_FAMILY_GGA,
  {&xc_ref_Margraf2019_244116, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {N_PAR, names, desc, par_ccdf, set_ext_params_cpy},
  gga_c_ccdf_init, NULL,
  NULL, work_gga, NULL
};
