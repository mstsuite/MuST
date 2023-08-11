/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MGGAC  711 /* MGGAC exchange of Patra et al */

GPU_FUNCTION double xc_mgga_x_mbrxc_get_x(double Q);

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_mggac.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mggac = {
  XC_MGGA_X_MGGAC,
  XC_EXCHANGE,
  "MGGAC exchange of Patra et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Patra2019_155140, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-12,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
