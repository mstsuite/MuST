/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MBRXH_BG  697 /* Modified Becke-Roussel for band gaps - hydrogen hole */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_mbrxh_bg.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mbrxh_bg = {
  XC_MGGA_X_MBRXH_BG,
  XC_EXCHANGE,
  "Modified Becke-Roussel for band gaps - hydrogen hole",
  XC_FAMILY_MGGA,
  {&xc_ref_Patra2019_045147, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS | XC_FLAGS_DEVELOPMENT,
  1.0e-12,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
