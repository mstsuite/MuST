/*
 Copyright (C) 2008 Georg Madsen

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_PBEA  121 /* Madsen (PBE-like) */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_pbea.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbea = {
  XC_GGA_X_PBEA,
  XC_EXCHANGE,
  "Madsen 07",
  XC_FAMILY_GGA,
  {&xc_ref_Madsen2007_195108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};
