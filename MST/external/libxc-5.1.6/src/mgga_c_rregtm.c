/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_RREGTM         391 /* revised regTM by Jana et al */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_rregtm.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_rregtm = {
  XC_MGGA_C_RREGTM,
  XC_CORRELATION,
  "Revised regTM correlation by Jana et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Jana2021_024103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
