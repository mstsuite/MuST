/*
 Copyright (C) 2014 Jess Wellendorff, M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_X_BEEFVDW          285 /* BEEF-vdW exchange */
#define XC_GGA_XC_BEEFVDW         286 /* BEEF-vdW exchange-correlation */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_beefvdw.c"
#include "work_gga.c"


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_beefvdw = {
  XC_GGA_X_BEEFVDW,
  XC_EXCHANGE,
  "BEEF-vdW exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Wellendorff2012_235149, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL,
};


void
gga_xc_beefvdw_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_GGA_X_BEEFVDW, XC_LDA_C_PW_MOD, XC_GGA_C_PBE};
  static double funcs_coef[3] = {1.0, 0.6001664769, 1.0 - 0.6001664769};

  xc_mix_init(p, 3, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_beefvdw = {
  XC_GGA_XC_BEEFVDW,
  XC_EXCHANGE_CORRELATION,
  "BEEF-vdW exchange-correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Wellendorff2012_235149, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {0, NULL, NULL, NULL, NULL},
  gga_xc_beefvdw_init, NULL,
  NULL, NULL, NULL,
};
