/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_XC_TIH   599   /* Neural network LDA from Tozer et al */

#define XC_NO_EXC
#include "decl_lda.h"
#include "maple2c/lda_vxc/lda_xc_tih.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_tih = {
  XC_LDA_XC_TIH,
  XC_EXCHANGE_CORRELATION,
  "Neural network LDA from Tozer et al",
  XC_FAMILY_LDA,
  {&xc_ref_Tozer1996_9200, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  5e-24,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  work_lda, NULL, NULL
};
