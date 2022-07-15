/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_TCA          100 /* Tognetti, Cortona, Adamo */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_tca.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_tca = {
  XC_GGA_C_TCA,
  XC_CORRELATION,
  "Tognetti, Cortona, Adamo",
  XC_FAMILY_GGA,
  {&xc_ref_Tognetti2008_034101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};
