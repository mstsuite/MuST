/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_X_TH          225 /* Tsuneda and Hirao */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_th.c"
#include "work_mgga.c"

/*
Although it reproduces published values, the functional does not
appear to be numerically stable; self-consistent calculations with it
diverge; this is why the functional is marked development.
*/

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_th = {
  XC_MGGA_X_TH,
  XC_EXCHANGE,
  "Tsuneda and Hirao",
  XC_FAMILY_MGGA,
  {&xc_ref_Tsuneda2000_15527, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS | XC_FLAGS_DEVELOPMENT,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
