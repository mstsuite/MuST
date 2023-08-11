/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_XC_ZLP     43   /* Zhao, Levy & Parr, Eq. (20)  */

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_xc_zlp.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_zlp = {
  XC_LDA_XC_ZLP,
  XC_EXCHANGE_CORRELATION,
  "Zhao, Levy & Parr, Eq. (20)",
  XC_FAMILY_LDA,
  {&xc_ref_Zhao1993_918, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  work_lda, NULL, NULL
};
