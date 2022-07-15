/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_C_OP_PW91    262 /* one-parameter progressive functional (PW91 version)  */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_op_pw91.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_op_pw91 = {
  XC_GGA_C_OP_PW91,
  XC_CORRELATION,
  "one-parameter progressive functional (PW91 version)",
  XC_FAMILY_GGA,
  {&xc_ref_Tsuneda1999_10664, &xc_ref_Tsuneda1999_5656, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  /* This functional is ridiculously sensitive to the density
     threshold and will give ridiculously low or high energies for too
     low thresholds. Due to the sensitivity, it's not clear what the
     cutoff should even be, since the energy doesn't converge at any
     point. 1e-10 is a decent lower limit, since anything below that
     blows up atomic correlation energies, even though Li is already
     wrong with this value.
   */
  1e-10,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};

