/*
 Copyright (C) 2017 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_PBETRANS  291 /* Gradient-based interpolation between PBE and revPBE */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_pbetrans.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbetrans = {
  XC_GGA_X_PBETRANS,
  XC_EXCHANGE,
  "Gradient-regulated connection-based correction for the PBE exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Bremond2016_1059, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};
