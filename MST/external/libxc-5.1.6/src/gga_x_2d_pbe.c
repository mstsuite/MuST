/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_2D_PBE          129 /* Perdew, Burke & Ernzerhof exchange in 2D          */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_2d_pbe.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_2d_pbe = {
  XC_GGA_X_2D_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof in 2D",
  XC_FAMILY_GGA,
  {&xc_ref_Vilhena2014, NULL, NULL, NULL, NULL},
  XC_FLAGS_2D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};
