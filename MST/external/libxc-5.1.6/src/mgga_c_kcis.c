/*
 Copyright (C) 2008 Georg Madsen

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_C_KCIS         562 /* Krieger, Chen, Iafrate, and Savin */
#define XC_HYB_MGGA_XC_B0KCIS  563 /* Hybrid based on KCIS */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_kcis.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_kcis = {
  XC_MGGA_C_KCIS,
  XC_CORRELATION,
  "Krieger, Chen, Iafrate, and Savin",
  XC_FAMILY_MGGA,
  {&xc_ref_Rey1998_581, &xc_ref_Krieger1999_463, &xc_ref_Krieger2001_48, &xc_ref_Kurth1999_889, &xc_ref_Toulouse2002_10465},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga
};

/*************************************************************/
void
xc_hyb_mgga_xc_b0kcis_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_B88, XC_MGGA_C_KCIS};
  static double funcs_coef[2] = {1.0 - 0.25, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.25);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_b0kcis = {
  XC_HYB_MGGA_XC_B0KCIS,
  XC_EXCHANGE_CORRELATION,
  "Hybrid based on KCIS",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Toulouse2002_10465, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_mgga_xc_b0kcis_init, NULL,
  NULL, NULL, work_mgga
};
