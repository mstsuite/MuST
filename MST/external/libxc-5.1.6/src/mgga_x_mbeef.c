/*
 Copyright (C) 2014 Jess Wellendorff, M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MBEEF          249 /* mBEEF exchange */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_mbeef.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mbeef = {
  XC_MGGA_X_MBEEF,
  XC_EXCHANGE,
  "mBEEF exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Wellendorff2014_144107, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
