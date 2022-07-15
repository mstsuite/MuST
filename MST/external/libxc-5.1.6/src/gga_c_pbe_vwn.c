/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_PBE_VWN      216 /* Perdew, Burke & Ernzerhof correlation based on VWN LDA */

typedef struct{
  double beta, gamma, BB;
} gga_c_pbe_vwn_params;

static void gga_c_pbe_vwn_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_pbe_vwn_params));
}

#define PBEVWN_N_PAR 3
static const char  *pbevwn_names[PBEVWN_N_PAR]  = {"_beta", "_gamma", "_B"};
static const char  *pbevwn_desc[PBEVWN_N_PAR]   = {
  "beta constant",
  "(1 - ln(2))/Pi^2 in the PBE",
  "Multiplies the A t^2 term. Used in the SPBE functional"
};
static const double pbevwn_values[PBEVWN_N_PAR] = {
  0.06672455060314922, 0.031090690869654895034, 1.0
};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_pbe_vwn.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_pbe_vwn = {
  XC_GGA_C_PBE_VWN,
  XC_CORRELATION,
  "Perdew, Burke & Ernzerhof based on VWN correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Kraisler2010_042516, &xc_ref_Perdew1996_3865, &xc_ref_Perdew1996_3865_err, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {PBEVWN_N_PAR, pbevwn_names, pbevwn_desc, pbevwn_values, set_ext_params_cpy},
  gga_c_pbe_vwn_init, NULL,
  NULL, work_gga, NULL
};
