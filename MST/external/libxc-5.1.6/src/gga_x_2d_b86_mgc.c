/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_2D_B86_MGC      124 /* Becke 86 MGC for 2D systems */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_2d_b86_mgc.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_2d_b86_mgc = {
  XC_GGA_X_2D_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction for 2D",
  XC_FAMILY_GGA,
  {&xc_ref_Pittalis2009_012503, NULL, NULL, NULL, NULL},
  XC_FLAGS_2D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};
