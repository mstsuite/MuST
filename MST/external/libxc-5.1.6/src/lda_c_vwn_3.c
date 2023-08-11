/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_VWN_3   30   /* Vosko, Wilk, & Nusair (3)   */

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_vwn_3.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_vwn_3 = {
  XC_LDA_C_VWN_3,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair (VWN3)",
  XC_FAMILY_LDA,
  {&xc_ref_Vosko1980_1200, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL, work_lda,
  NULL, NULL
};
