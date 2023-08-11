/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_X_REL   532   /* Relativistic exchange        */

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_x_rel.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_x_rel = {
  XC_LDA_X_REL,
  XC_EXCHANGE,
  "Slater exchange with relativistic corrections",
  XC_FAMILY_LDA,
  {&xc_ref_Rajagopal1978_L943, &xc_ref_MacDonald1979_2977, &xc_ref_Engel1995_2750, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  work_lda, NULL, NULL
};
