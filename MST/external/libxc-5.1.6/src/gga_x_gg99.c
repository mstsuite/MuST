/*
 Copyright (C) 2015 M.A.L. Marques, Markus Patzold

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_GG99   535 /* Gilbert and Gill 1999 */
#define XC_GGA_X_KGG99  544 /* Gilbert and Gill 1999 (mixed) */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_gg99.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_gg99 = {
  XC_GGA_X_GG99,
  XC_EXCHANGE,
  "Gilbert and Gill 1999",
  XC_FAMILY_GGA,
  {&xc_ref_Gilbert1999_511, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  5e-7,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};


/*************************************************************/
static void
gga_x_kgg_init(xc_func_type *p)
{
  /* defined in Eq. (25) of the paper */
  static int   funcs_id  [2] = {XC_LDA_X, XC_GGA_X_GG99};
  static double funcs_coef[2] = {-0.047/X_FACTOR_C, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_kgg99 = {
  XC_GGA_X_KGG99,
  XC_EXCHANGE,
  "Gilbert and Gill 1999 (mixed)",
  XC_FAMILY_GGA,
  {&xc_ref_Gilbert1999_511, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  5e-7,
  {0, NULL, NULL, NULL, NULL},
  gga_x_kgg_init, NULL,
  NULL, NULL, NULL
};

