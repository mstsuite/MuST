/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_ACGGA           39 /* acGGA, asymptotically corrected GGA */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_acgga.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_acgga = {
  XC_GGA_C_ACGGA,
  XC_CORRELATION,
  "acGGA, asymptotically corrected GGA correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Cancio2018_084116, &xc_ref_Burke2014_4834, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};
