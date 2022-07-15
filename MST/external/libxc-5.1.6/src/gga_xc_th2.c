/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_XC_TH2          155 /* Tozer and Handy v. 2 */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_xc_th2.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_th2 = {
  XC_GGA_XC_TH2,
  XC_EXCHANGE_CORRELATION,
  "Tozer and Handy v. 2",
  XC_FAMILY_GGA,
  {&xc_ref_Tozer1998_3162, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};
