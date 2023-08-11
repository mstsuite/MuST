/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HYB_MGGA_XC_B98         598 /* Becke 98 */

static void
hyb_mgga_xc_b98_init(xc_func_type *p)
{
  xc_hyb_init_hybrid(p, 0.1985);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_xc_b98.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_b98 = {
  XC_HYB_MGGA_XC_B98,
  XC_EXCHANGE_CORRELATION,
  "Becke 98",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Becke1998_2092, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  hyb_mgga_xc_b98_init, NULL,
  NULL, NULL, work_mgga,
};

