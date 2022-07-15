/*
 Copyright (C) 2008 Georg Madsen

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_C_KCISK         638 /* Krieger, Chen, and Kurth */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_kcisk.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_kcisk = {
  XC_MGGA_C_KCISK,
  XC_CORRELATION,
  "Krieger, Chen, and Kurth",
  XC_FAMILY_MGGA,
  {&xc_ref_Rey1998_581, &xc_ref_Krieger1999_463, &xc_ref_Krieger2001_48, &xc_ref_Kurth1999_889, &xc_ref_Toulouse2002_10465},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga
};
