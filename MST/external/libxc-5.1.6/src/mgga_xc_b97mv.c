/*
 Copyright (C) 2015 Narbe Mardirossian and Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_XC_B97M_V        254 /* Mardirossian and Head-Gordon */

static void
mgga_xc_b97mv_init(xc_func_type *p)
{
  /* Non-local correlation parameters */
  p->nlc_b = 6.0;
  p->nlc_C = 0.01;
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_xc_b97mv.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_xc_b97m_v = {
  XC_MGGA_XC_B97M_V,
  XC_EXCHANGE_CORRELATION,
  "B97M-V exchange-correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Mardirossian2015_074111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_VV10 | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_xc_b97mv_init, NULL,
  NULL, NULL, work_mgga,
};
