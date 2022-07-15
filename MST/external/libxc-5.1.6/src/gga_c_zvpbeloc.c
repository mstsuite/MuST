/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_ZVPBELOC       606 /* PBEloc variation with enhanced compatibility with exact exchange */
#define XC_HYB_GGA_XC_APBE0     607 /* Hybrid based on APBE */
#define XC_HYB_GGA_XC_HAPBE     608 /* Hybrid based in APBE and zvPBEloc */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_zvpbeloc.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_zvpbeloc = {
  XC_GGA_C_ZVPBELOC,
  XC_CORRELATION,
  "PBEloc variation with enhanced compatibility with exact exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Fabiano2015_122, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-10,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};

void
xc_hyb_gga_xc_apbe0_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_APBE, XC_GGA_C_APBE};
  static double funcs_coef[2] = {0.75, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.25);
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_apbe0 = {
  XC_HYB_GGA_XC_APBE0,
  XC_EXCHANGE_CORRELATION,
  "Hybrid based on APBE",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Fabiano2015_122, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-10,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_apbe0_init, NULL,
  NULL, NULL, NULL
};

void
xc_hyb_gga_xc_hapbe_init(xc_func_type *p)
{
  static int    funcs_id  [3] = {XC_GGA_X_APBE, XC_GGA_C_APBE, XC_GGA_C_ZVPBELOC};
  static double funcs_coef[3] = {0.80, 0.80, 0.20};

  xc_mix_init(p, 3, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.20);
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_hapbe = {
  XC_HYB_GGA_XC_HAPBE,
  XC_EXCHANGE_CORRELATION,
  "Hybrid based in APBE and zvPBEloc",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Fabiano2015_122, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-10,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_hapbe_init, NULL,
  NULL, NULL, NULL
};
