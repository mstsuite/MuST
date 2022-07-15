/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_CS1          565 /* A dynamical correlation functional */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_cs1.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_cs1 = {
  XC_GGA_C_CS1,
  XC_CORRELATION,
  "A dynamical correlation functional",
  XC_FAMILY_GGA,
  {&xc_ref_Handy2002_5411, &xc_ref_Proynov2006_436, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-20,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};
