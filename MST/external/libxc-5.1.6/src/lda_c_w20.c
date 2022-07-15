/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_C_W20  317   /* Xie, Wu, and Zhao correlation  */

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_w20.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_w20 = {
  XC_LDA_C_W20,
  XC_CORRELATION,
  "Xie, Wu, and Zhao interpolation ansatz without fitting parameters",
  XC_FAMILY_LDA,
  {&xc_ref_Xie2021_045130, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  work_lda, NULL, NULL
};
