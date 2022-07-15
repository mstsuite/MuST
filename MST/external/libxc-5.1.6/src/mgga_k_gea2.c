/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_K_GEA2      627 /* Second-order gradient expansion */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_k_gea2.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_gea2 = {
  XC_MGGA_K_GEA2,
  XC_KINETIC,
  "Second-order gradient expansion",
  XC_FAMILY_MGGA,
  {&xc_ref_Kompaneets1956_427, &xc_ref_Kirznits1957_115, NULL, NULL, NULL},
  XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga
};
