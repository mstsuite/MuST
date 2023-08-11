/*
 Copyright (C) 2016 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_REVSCAN       582 /* revised SCAN correlation */
#define XC_MGGA_C_REVSCAN_VV10  585 /* revised SCAN correlation */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_revscan.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_revscan = {
  XC_MGGA_C_REVSCAN,
  XC_CORRELATION,
  "revised SCAN",
  XC_FAMILY_MGGA,
  {&xc_ref_Mezei2018_2469, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};

static void
mgga_c_revscan_vv10_init(xc_func_type *p)
{
  static int   funcs_id  [1] = {XC_MGGA_C_REVSCAN};
  static double funcs_coef[1] = {1.0};

  xc_mix_init(p, 1, funcs_id, funcs_coef);

  p->nlc_b = 9.8;
  p->nlc_C = 0.0093;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_revscan_vv10 = {
  XC_MGGA_C_REVSCAN_VV10,
  XC_CORRELATION,
  "REVSCAN + VV10 correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mezei2018_2469, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_VV10 | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_c_revscan_vv10_init, NULL,
  NULL, NULL, NULL
};
