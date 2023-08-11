/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_REVTPSS       241 /* revised TPSS correlation */
#define XC_MGGA_C_REVTM         694 /* revised Tao and Mo 2016 correlation */

typedef struct{
  double d;
  double C0_c[4];
} mgga_c_revtpss_params;

static void
mgga_c_revtpss_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_revtpss_params));
}

#define REVTPSS_N_PAR 5
static const char  *revtpss_names[REVTPSS_N_PAR]  = {"_d", "_C0_c0", "_C0_c1", "_C0_c2", "_C0_c3"};
static const char  *revtpss_desc[REVTPSS_N_PAR]   = {"d", "C0_c0", "C0_c1", "C0_c2", "C0_c3"};
static const double revtpss_values[REVTPSS_N_PAR] = {
  2.8, 0.59, 0.9269, 0.6225, 2.1540
};
static const double revtm_values[REVTPSS_N_PAR] = {
  2.8, 0.0, 0.1, 0.32, 0.0
};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_revtpss.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_revtpss = {
  XC_MGGA_C_REVTPSS,
  XC_CORRELATION,
  "revised TPSS correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2009_026403, &xc_ref_Perdew2009_026403_err, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-13, /* densities smaller than 1e-26 give NaNs */
  {REVTPSS_N_PAR, revtpss_names, revtpss_desc, revtpss_values, set_ext_params_cpy},
  mgga_c_revtpss_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_revtm = {
  XC_MGGA_C_REVTM,
  XC_CORRELATION,
  "revised Tao and Mo 2016 exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Jana2019_6356, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-13, /* densities smaller than 1e-26 give NaNs */
  {REVTPSS_N_PAR, revtpss_names, revtpss_desc, revtm_values, set_ext_params_cpy},
  mgga_c_revtpss_init, NULL,
  NULL, NULL, work_mgga
};

